# Temporary Innovations and Future Investment Opportunities

Reproduction materials for manuscript **v0.43**, prepared for *Economic Theory*.
This identifies the submission target, not an acceptance or publication status.

**Author:** Hyunkyu Lee, Department of Economics, Kyung Hee University, Seoul,
Republic of Korea. **Contact:** richardhk2@khu.ac.kr.
**ORCID:** https://orcid.org/0000-0001-9134-1748.

A better replacement technology can raise the reservation payment for the
innovation it makes obsolete. These materials reproduce the exact examples,
resource bounds, arrival comparisons, and one-state concave-economy results,
and the compensated ordering and universal-curvature results. The central
comparison lets both access and no-access economies exploit the same later
technology. Theorem proofs are in the manuscript; the computations check
consequences and provide the concave-economy examples.

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
| Theorem 2; constructive necessity | `characterization_analysis.py` | `generated/characterization_witnesses.csv` (48 exact geometric-stream comparisons) |
| Proposition 3; logarithmic neutrality boundaries | `neutrality_analysis.py` | `generated/neutrality_boundaries.csv` (24 policy/value comparisons) |
| Section 3.2; one preference, two resource levels | `resource_example.py` | `generated/resource_level_example.csv`, `checks/resource_example_analysis.json` (70-decimal check) |
| Corollary 1; Table 1; Fig. 1 | `dynamic_analysis.py` and `render_figure.py` | `generated/dynamic_examples.csv`, `checks/dynamic_analysis.json`, `figures/Fig1.pdf` |
| Proposition 2; Table 2; resource ceiling | `scale_analysis.py` | `generated/scale_examples.csv`, `checks/scale_analysis.json` |
| Table 3; curvature × survival grid | `nonlinear_analysis.py` | `generated/nonlinear_examples.csv`, `generated/nonlinear_paths.npz`, `checks/nonlinear_analysis.json` |
| Appendix Table 4; joint arrival and quality | `joint_analysis.py` | `generated/joint_examples.csv`, `checks/joint_analysis.json` |

Table 2 uses `generated/scale_table_main.tex` (T=5,10). The full
`scale_examples.csv` and `scale_table.tex` additionally retain T=20,40
as different, longer-lived innovations. Other table fragments are included
without subsetting. Figure 1 data are in `reservation_curves.csv` and
`compensated_consumption_paths.csv`. Its right panel uses each preference's
own reservation fee at common return 1.06, not the fixed adoption fee.
For the resource-level example, the derivative column reports the local response
of F-star/X to log R1, holding X and the preference fixed.
The scripts run sequentially; `scale_analysis.py` and `joint_analysis.py`
import formulas from `dynamic_analysis.py`.
The resource-level script imports the independent floating-point implementation
from `characterization_analysis.py`.

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
finite-horizon tests of consumption, capital and investment ordering and
the fee derivative. `characterization_witnesses.csv` checks the necessity
construction with an exact three-group reduction of the infinite stream.
`neutrality_boundaries.csv` compares the two logarithmic formulas with direct
policy simulation plus exact infinite continuation. Bounded quality paths
and interior investments are enforced in these examples.

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

## Revision v0.43

Physical initial capital is now explicitly held fixed across access and outside
allocations; the fee is charged to the date-zero resource budget. The programs
check that budget separately from all later accumulation equations. No published
table entry or nonlinear path changes as a result of the clarification.
The logarithmic program additionally checks the initial-fee interiority bound.

`resource_example.py` adds two exact infinite-stream comparisons using the same
utility and technology but different initial resources. It solves the budget
and compensation equations at 70-decimal precision, compares the analytic fee
derivative with a central finite difference, and cross-checks the earlier
floating-point solver. These are local theoretical examples, not an empirical
wealth threshold or a global wealth monotonicity result. The figure now displays
positive ticks above zero in the compensated-consumption panel.

The full v0.43 run uses Python 3.12.14 on macOS ARM64 and the pinned versions in
`requirements.txt`. It takes approximately ten seconds in the author's
environment; execution time depends on hardware. The manuscript and submission
guide identify an immutable Git commit, so later repository updates need
not change the downloadable submission version.
