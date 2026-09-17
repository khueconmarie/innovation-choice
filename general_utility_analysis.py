"""Independent finite-horizon checks of the general-utility theorem.

The marginal-utility equation is solved directly in consumption units. These
checks exercise the theorem's finite analogue; the infinite-horizon claim is
proved in the manuscript, not inferred from a finite simulation.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

P = Path(__file__).resolve().parent


class Utility:
    def __init__(self, gammas, weights):
        self.gammas = np.asarray(gammas, dtype=float)
        self.weights = np.asarray(weights, dtype=float)

    def value(self, c):
        c = np.asarray(c)
        return sum(w * (np.log(c) if g == 1 else c ** (1-g)/(1-g))
                   for g, w in zip(self.gammas, self.weights))

    def marginal(self, c):
        return sum(w * np.asarray(c) ** (-g)
                   for g, w in zip(self.gammas, self.weights))

    def demand(self, log_marginal):
        x = -np.asarray(log_marginal) / np.dot(self.gammas, self.weights)
        for _ in range(80):
            a = np.log(self.weights)[:, None] - self.gammas[:, None]*x
            lm = logsumexp(a, axis=0)
            elasticity = np.sum(np.exp(a-lm)*self.gammas[:, None], axis=0)
            step = (lm-log_marginal)/elasticity
            x += step
            if np.max(np.abs(step)/(1+np.abs(x))) < 2e-14:
                return np.exp(x)
        raise AssertionError('Marginal-utility inversion did not converge')


def optimum(utility, returns, beta, resources):
    log_g = np.r_[0., np.cumsum(np.log(returns))]
    t = np.arange(len(log_g))
    disc = beta**t
    def allocation(log_lam):
        return utility.demand(log_lam-t*np.log(beta)-log_g)
    def budget(log_lam):
        return np.sum(allocation(log_lam)*np.exp(-log_g))-resources
    ll = brentq(budget, -100, 100, xtol=5e-14)
    c = allocation(ll)
    return dict(c=c, lam=np.exp(ll), value=np.dot(disc, utility.value(c)),
                budget=np.sum(c*np.exp(-log_g)), log_g=log_g)


def reservation(utility, r0, ri, beta, x):
    outside = optimum(utility, r0, beta, x)
    cap = -np.expm1(-np.sum(np.log(ri/r0)))
    def gap(p):
        return optimum(utility, ri, beta, x*(1-p))['value']-outside['value']
    p = brentq(gap, 0, cap, xtol=5e-14)
    access = optimum(utility, ri, beta, x*(1-p))
    return p*x, outside, access


def main():
    rng = np.random.default_rng(410917)
    specifications = [([.3, .7], [.4, .6], 1, 'curvature_below_one'),
                      ([1.3, 1.7], [.6, .4], -1, 'curvature_above_one'),
                      ([1.], [1.], 0, 'logarithmic'),
                      ([.5, 1.5], [.5, .5], None, 'crossing_curvature')]
    records = []
    for gs, ws, expected_sign, name in specifications:
        u = Utility(gs, ws)
        for _ in range(20):
            horizon = int(rng.integers(22, 42))
            T = int(rng.integers(2, 9))
            k = int(rng.integers(T, horizon-1))
            beta = float(rng.uniform(.90, .97))
            x = float(rng.uniform(.5, 3.))
            r0 = rng.uniform(.99, 1.055, horizon)
            ri = r0.copy()
            ri[:T] *= np.exp(rng.uniform(.003, .035, T))
            fee, b, a = reservation(u, r0, ri, beta, x)
            ell = np.exp(a['log_g']-b['log_g'])
            ratio = a['lam']/b['lam']
            assert 1 < ratio < ell[-1]
            assert a['c'][0] < b['c'][0]
            assert np.all(a['c'][T:] > b['c'][T:])
            assert np.all(np.sign(a['c']-b['c']) == np.sign(ell-ratio))
            t = np.arange(horizon+1)
            derivative = np.sum((beta**t)[t>k]*(a['c']*u.marginal(a['c'])-
                              b['c']*u.marginal(b['c']))[t>k])/a['lam']
            step = 2e-5
            fees = []
            for z in [-step, step]:
                rb, ra = r0.copy(), ri.copy()
                rb[k] *= np.exp(z); ra[k] *= np.exp(z)
                fees.append(reservation(u, rb, ra, beta, x)[0])
            finite_difference = (fees[1]-fees[0])/(2*step)
            # Integrating u'(c)*(1-R_u(c)) gives c*u'(c) at the endpoints.
            primitive_difference = sum(w*(a['c']**(1-g)-b['c']**(1-g))
                                       for g, w in zip(gs, ws))
            integral = np.sum((beta**t)[t>k]*primitive_difference[t>k])/a['lam']
            if expected_sign == 0:
                assert abs(derivative) < 1e-12
            elif expected_sign is not None:
                assert derivative*expected_sign > 0
            # A finite common improvement also obeys the curvature comparison.
            rb, ra = r0.copy(), ri.copy()
            rb[k:] *= 1.004; ra[k:] *= 1.004
            improved_fee = reservation(u, rb, ra, beta, x)[0]
            if expected_sign:
                assert expected_sign*(improved_fee-fee) > 0
            elif expected_sign == 0:
                assert abs(improved_fee-fee) < 2e-12
            records.append(dict(utility=name, horizon=horizon, T=T, k=k,
                beta=beta, X=x, fee=fee, fee_fraction=fee/x,
                multiplier_ratio=ratio, cumulative_advantage=float(ell[-1]),
                resource_cap=float(1-1/ell[-1]), derivative=float(derivative),
                finite_difference=float(finite_difference), improved_fee=improved_fee,
                derivative_error=abs(derivative-finite_difference),
                integral_error=abs(derivative-integral),
                budget_error=max(abs(b['budget']-x),abs(a['budget']-(x-fee))),
                compensation_error=abs(a['value']-b['value'])))
    summary = dict(seed=410917, cases=len(records),
        scope='Independent finite-horizon dual allocations, not an infinite-horizon numerical proof.',
        max_derivative_error=max(r['derivative_error'] for r in records),
        max_integral_error=max(r['integral_error'] for r in records),
        max_budget_error=max(r['budget_error'] for r in records),
        max_compensation_error=max(r['compensation_error'] for r in records),
        min_resource_cap_slack=min(r['resource_cap']-r['fee_fraction'] for r in records),
        records=records)
    assert summary['max_derivative_error'] < 3e-8
    assert summary['max_budget_error'] < 2e-12
    assert summary['max_integral_error'] < 2e-12
    assert summary['min_resource_cap_slack'] > 0
    (P/'checks').mkdir(exist_ok=True)
    (P/'checks/general_utility_analysis.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k != 'records'}, indent=2))


if __name__ == '__main__':
    main()
