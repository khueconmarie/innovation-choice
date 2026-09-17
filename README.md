# Substitution and Complementarity in Innovation Choice

Reproduction materials for manuscript **v0.41**, prepared for *Economic Theory*.
This identifies the submission target, not an acceptance or publication status.

**Author:** Hyunkyu Lee, Department of Economics, Kyung Hee University, Seoul,
Republic of Korea. **Contact:** richardhk2@khu.ac.kr.

A better replacement technology can raise the reservation payment for the
innovation it makes obsolete. These materials reproduce the exact examples,
resource bounds, arrival comparisons, and one-state concave-economy results,
and check the finite-horizon analogue of the general-utility theorem.

## Download and run

Repository: https://github.com/khueconmarie/innovation-choice

Download the repository ZIP using GitHub's **Code → Download ZIP**, or clone
the repository:

```sh
git clone https://github.com/khueconmarie/innovation-choice.git
cd innovation-choice
```
 Python 3.12 is the tested runtime. No GPU, credentials, network
API, proprietary dataset, or paid service is required after dependencies are installed.

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
```

On Windows, activate with `.venv\Scripts\activate` instead. The complete run
reoptimizes all nonlinear allocations, including the longer-horizon checks;
it does not merely read the saved table. Timings and package versions are
recorded in `checks/reproduction_run.json`. Individual script output is saved
in `logs/`.

To inspect the bundled numerical results without solving again:

```sh
python verify_outputs.py
```

That command verifies saved numbers against the frozen reference outputs.
It is a consistency check, not a fresh solution of the model.

## Result-to-program map

| Manuscript result | Program | Data / checks |
|---|---|---|
| Theorem 1, general utility | `general_utility_analysis.py` | `checks/general_utility_analysis.json` (80 finite-horizon comparisons) |
| Corollary 1; Table 1; Fig. 1 | `dynamic_analysis.py` and `render_figure.py` | `generated/dynamic_examples.csv`, `checks/dynamic_analysis.json`, `figures/Fig1.pdf` |
| Table 2; resource ceiling | `scale_analysis.py` | `generated/scale_examples.csv`, `checks/scale_analysis.json` |
| Table 3; curvature × survival grid | `nonlinear_analysis.py` | `generated/nonlinear_examples.csv`, `generated/nonlinear_paths.npz`, `checks/nonlinear_analysis.json` |
| Appendix Table 4; joint arrival and quality | `joint_analysis.py` | `generated/joint_examples.csv`, `checks/joint_analysis.json` |

`generated/*_table.tex` are the exact table fragments included in the paper.
The scripts run sequentially; `scale_analysis.py` and `joint_analysis.py`
import formulas from `dynamic_analysis.py`.

## Data dictionary and scope

All numerical inputs are **chosen theoretical examples**. The package contains
no observational data or estimated parameters. Fee fractions are dimensionless:
linear examples use initial resources X; nonlinear examples use initial output
Y0. A percentage-point effect is 100 times the difference of two fee fractions.
The `fund` column uses the strict criterion F < F-star. EIS equals 1/gamma.

In `nonlinear_paths.npz`, keys have the form
`eta0.9_sigma0.2_gamma0.5_q1.041_consumption` (and `_investment`, `_path`).
`path` contains capital at dates 0 through N; consumption and investment
contain dates 0 through N-1. These are the access allocations at their
compensating fees. Outside-economy paths are recomputed within the solver.

The 2×2 experiment changes eta and sigma, with A=(1/beta-sigma)/(eta*qbase)
to maintain the initial steady capital stock K0=1. A is therefore not fixed
across cells. Numerical value intervals use a feasible continuation, an affine
upper bound, a global gradient correction, and an explicit rounding allowance;
they are not formal directed-rounding interval certificates. The infinite-horizon
general-utility theorem is proved analytically; its 80 numerical checks are
finite-horizon tests of consumption ordering and the fee derivative.

## Reproduction evidence

`reference/` preserves the release's CSV/JSON/NPZ/TeX output. A full run
regenerates `generated/` and `checks/`, then compares table entries, saved paths,
and the declared numerical bounds against the reference. Numerical tolerances
allow minor floating-point differences. PDF binary hashes need not match because
PDF creation metadata can change; figure data are regenerated from the exact formulas.

The code and generated results are provided for scholarly verification of this
paper. No third-party source articles, private working notes, credentials,
or unrelated research materials are included.

The journal figure uses distinct line styles as well as color. `render_figure.py`
regenerates `Fig1.pdf` from the same closed-form data with embedded TrueType fonts.
