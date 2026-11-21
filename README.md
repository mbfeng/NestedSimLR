# Efficient Nested Simulation Experiment Design via the Likelihood Ratio Method — code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24409919.svg)](https://doi.org/10.5281/zenodo.24409919)

Code that reproduces every table and figure in Section 7 of

> Feng, B. M., & Song, E. (2025). Efficient nested simulation experiment design via the likelihood
> ratio method. *INFORMS Journal on Computing*, 37(3), 723–742.
> https://doi.org/10.1287/ijoc.2022.0392

This is the **maintained version** of the code. The snapshot that accompanied the paper is
https://github.com/INFORMSJoC/2022.0392 (code DOI 10.1287/ijoc.2022.0392.cd). If you use this code,
please cite the paper. The snapshot is frozen; this version is easier to read and extend, adds
missing seeds, runs about 10–70 times faster, and checks itself against the published tables.

## Quick start

```r
# R >= 4.3; from this folder:
install.packages(c("lpSolve", "lpSolveAPI", "mvtnorm", "randcorr", "matrixStats",
                   "future", "future.apply", "ggplot2", "latex2exp"))
source("reproduce_paper.R")   # or: Rscript reproduce_paper.R
```

By default (`RERUN_SIMULATIONS <- FALSE`) this redraws all figures (`figures/`) and tables
(`tables/`) from the summaries of the simulation runs reported in the paper
(`results/published_run/`, 13 MB) in about 30 seconds, and then runs `verify_against_paper.R`.

Set `RERUN_SIMULATIONS <- TRUE` in `reproduce_paper.R` to rerun every simulation from scratch (see
[Run times](#run-times)).

## Paper result → code

| Paper | Produced by | Simulation code | Summary file(s) in `results/…/` |
|---|---|---|---|
| Figure 1(a), 1(b) | `figure1()` in `R/figures.R` | closed form | — |
| Figure 2(a), 2(b) | `figure2()` | `scripts/sec71_simulate.R fig2` | `sec71_fig2_bands.csv` |
| Table 1 | `table1()` in `R/tables.R` | `scripts/sec71_simulate.R table1`, `truth` | `sec71_table1_*` |
| Table 2 | `table2()` | `scripts/sec72_simulate.R` | `sec72_budgets_per_macro.csv`, `sec72_cpu_times_*` |
| Table 3 | `table3()` | `scripts/sec72_simulate.R` | `sec72_per_macro.csv`, `sec72_truth.csv` |
| Table 4 | `table4()` | `scripts/sec73_simulate.R table4` | `sec73_table4_*` |
| Figure 3(a) | `figure3()` | `scripts/sec73_simulate.R fig3a` | `sec73_fig3a_var_ratio.csv` |
| Figure 3(b) | `figure3()` | `scripts/sec73_simulate.R fig3b` | `sec73_fig3b_budgets.csv` |
| Numbers in the text of §7 | `text73()`, `verify_against_paper.R` | as above | as above |

## Folder structure

```
reproduce_paper.R          entry point: figures + tables + verification
verify_against_paper.R     checks every published number (see Verification)
R/                         library (sourced by everything)
  lr_design.R              the method: Psi matrices, LP (11), pooling weights, pooled
                           self-normalized LR estimator, log likelihood ratios
  option_pricing.R         Black-Scholes, geometric Asian and barrier prices; barrier survival
  risk_measures.R          nested statistics (quantile, indicator, hockey stick, square), MSE
  sec71_single_asset.R     Section 7.1 model and designs (functions single_asset_*)
  sec72_multi_asset.R      Section 7.2 model and designs (functions multi_asset_*)
  sec73_newsvendor.R       Section 7.3 model and designs (functions newsvendor_*)
  utilities.R              loading the library, reproducible serial/parallel macro runs, timing
  tables.R, figures.R      Tables 1-4 and Figures 1-3 from the summaries
scripts/                   simulation drivers, one per section (write summaries)
results/published_run/     summaries of the runs reported in the paper (see Provenance)
tests/                     run_tests.R (does the library still give the reference results?),
                           test_barrier_survival.R (fast vs original barrier method)
tools/                     extract_published_runs.R (how the published-run summaries were made),
                           run_all_reruns.sh, assemble_rerun.R (the November 2026 rerun)
verify/                    the paper's printed values and the verification reports
figures/, tables/          output
```

## How the simulations work

All three examples use the same building blocks in `R/lr_design.R`:

1. `optimal_design()` solves LP (11) for the normalized budgets `c_j`, keeps the scenarios with
   `c_j > 0` ("sampled scenarios", `design$sampled`), sets `n_j = ceiling(N c_j)`
   (`design$n_inner`) and the pooling weights `gamma_ij ∝ n_j Psi_ij` (`design$weights`, Eq. (10)).
2. Inner replications are run only at the sampled scenarios.
3. `pooled_estimate()` computes the self-normalized LR estimate from each sampled scenario to every
   target scenario (`self_normalized_estimate()`) and pools them with the weights.

Each example file follows the same pattern, e.g. for §7.1: `single_asset_params()` (model),
`single_asset_scenarios()` (outer scenarios), `single_asset_true_value()` (closed-form μ(θ)),
`single_asset_inner()` (inner replications at one scenario), `single_asset_design()`,
`single_asset_optimal()`, `single_asset_standard()` (SNS, SNS+), `single_asset_regression()`.

| | §7.1 single-asset ERM | §7.2 multi-asset ERM | §7.3 newsvendor |
|---|---|---|---|
| Outer scenarios θ | M quantiles of S_τ (fixed) | M = 1000 draws of 10 correlated S_τ per macro run | M = 1000 posterior draws of λ per macro run |
| Inner output g(X) | straddle payoff | payoff of 30 options per stock | one-day profit |
| LR input | terminal price S_T | first-step price of each stock | 10 Poisson demands |
| Psi | `psi_lognormal()` | `psi_lognormal()` per stock | `psi_poisson()` |
| LP solver | lpSolve | lpSolveAPI (bounds, ≥10 reps, Psi < 1e-6 → 0) | lpSolve |
| Macro runs | 10,000 | 1,000 | 1,000 |
| SNS | ⌈Γ^{2/3}⌉ outer × ⌈Γ^{1/3}⌉ inner | same (Γ = average per-stock budget) | same |
| Regression | Laguerre order 0–2, ⌈Γ⌉ points | Laguerre order 0–2 per stock (31 coef.) | quadratic, no cross terms |

The optimal design in §7.2 is solved separately for each stock (the portfolio is a sum over stocks).

## Code notation → paper notation

| Code | Paper |
|---|---|
| `theta`, `lambda` | outer scenarios θ_i (S_τ in §7.1–7.2, λ in §7.3) |
| `n_outer`, `N_OUTER` | M |
| `target_ess`, `TARGET_ESS` | N (target ESS; also the SNS+ inner sample size) |
| `psi[i, j]` | Ψ_ij = 1 / E_j[W_ij²] |
| `design$sampled`, `design$n_inner`, `design$budget` | the scenarios with c_j > 0, n_j, and Γ = Σ_j n_j |
| `design$weights` | γ_ij |
| `log_lr` | log W_ij of every inner replication |
| `optimal`, `SNS`, `SNS+`, `regression`, `oracle` | the designs of the same names |
| `alpha`, `threshold` | VaR level and the threshold u in I(μ>u), (μ−u)⁺, ((μ−u)⁺)² |

## Reproducibility

* Every macro run k has its own L'Ecuyer-CMRG random-number stream derived from one seed per
  experiment (`run_macros()`, via `future.apply`), so results do not depend on the number of cores.
* The quantities that the 2024 code drew with fixed seeds are reproduced **exactly**: the
  correlation matrix in §7.2 (`set.seed(2023)`), the newsvendor posterior, the 10⁶ posterior
  draws used for CrI coverage, and the fixed scenarios of Figure 3(a) (`set.seed(2020)`). All
  optimal designs (budgets, sampled scenarios, weights) are deterministic given the scenarios and
  match the paper exactly (e.g., budget 2,148 for M = N = 1,000 in §7.1).
* Everything else uses different random numbers than in 2024, so reruns agree with the paper
  **within Monte Carlo error**, not digit for digit.

## Verification

`verify_against_paper.R <results_dir> [exact|mc]` compares the published values of Tables 1–4 (the
estimates; standard errors are not compared) and of the design budgets quoted in §7 with the
values computed from a results folder.

* **`results/published_run/` (exact mode):** the summaries come from the 2024 runs; Tables 1 and 4
  agree to the printed digits. For Tables 2–3 see Provenance.
* **Reruns (mc mode):** z = (x − v)/√(se_x² + se_v² + h²/3), where h is half a unit in the last
  printed digit; |z| > 2.5 is flagged.

**Full rerun (November 2026, `results/rerun_full/`, new random numbers, same sample sizes as the
paper):** every entry agrees within Monte Carlo error, except two Table 1 entries for the optimal
design at M = 512 (hockey-stick and square losses). Their squared errors are very heavy-tailed: a
handful of the 10,000 macro runs carry most of the MSE, so the MSE and its standard error vary
much more between runs than the standard error suggests. The per-macro estimates of these entries
have the same distribution as in the paper's run (Kolmogorov–Smirnov test);
`verify_against_paper.R` labels such entries "heavy-tailed" or "ok (KS)". The average SNS budget
in Table 2 is 0.4% below the published value (z = −3.1); the run behind Tables 2–3 was not kept
(see Provenance).

The figures redrawn from the rerun agree with the paper's. The deterministic quantities (optimal
design budgets, true values, correlation matrix, posterior) agree exactly.
`verify/verification_rerun_full.csv` has every entry. `tools/assemble_rerun.R` stacks the two
parts of the 10,000-run Table 1 rerun into `rerun_full`.

## Provenance of `results/published_run/`

The 2024 code saved complete R workspaces (9.9 GB in total). `tools/extract_published_runs.R`
reduced them to per-macro summaries (13 MB). The workspaces are kept offline by the authors.

* **Figure 2, Table 1:** the 10,000-macro runs reported in the paper.
* **Tables 2–3:** the run reported in the paper was **not kept**. The summaries come from a run of
  the same code with seeds 6001–6025 (25 × 40 = 1,000 macro runs, May 2024). It agrees with all
  28 tested entries of Tables 2–3 within Monte Carlo error (largest |z| = 1.8), but not digit for
  digit. A third saved run (seeds 2001–2025) used an earlier version of the code and does **not**
  agree with the paper; it was not used.
* **Table 4, Figure 3:** the runs reported in the paper.

## Run times

Full rerun, November 2026, with the paper's sample sizes. "Server": Linux, 2 × Xeon Gold 6126
(48 cores), 24 workers. "PC": Windows 11, Intel 14-core CPU, 12 workers.

| Part | Script | Wall-clock |
|---|---|---|
| Figure 2 (10,000 macro runs) | `sec71_simulate.R fig2` | 2 min (server) |
| Table 1, M = 128–1,024 (10,000 runs each) | `sec71_simulate.R table1` | 25 min (PC) |
| Table 1, M = 2,048 and 4,096 (10,000 runs each) | `sec71_simulate.R table1` | 50 min (server) |
| True risk measures (10⁸ scenarios) | `sec71_simulate.R truth` | 2 min (PC) |
| Tables 2–3 (1,000 runs) | `sec72_simulate.R` | 2.4 h (server); ≈ 3.4 min per macro run per core |
| Table 4 (1,000 runs), Figure 3(a) | `sec73_simulate.R table4`, `fig3a` | 11 min (PC) |
| Figure 3(b) (M up to 10,000, 100 runs each) | `sec73_simulate.R fig3b` | 49 min (server) |

Everything together: about 4.5 hours on 24 cores. Figure 3(b) at M = 10,000 needs about 3 GB of
memory per worker; set `MEM_BUDGET_GB` (see `tools/run_all_reruns.sh`) to cap the workers.

On Linux, R is often linked to a multithreaded BLAS. Set `OPENBLAS_NUM_THREADS=1` (the driver
script does), or every parallel worker starts one BLAS thread per core.

## Change log

* **2024-05** Code published with the paper (INFORMSJoC/2022.0392).
* **2026-11** Maintained version (this folder):
  * one library for all three examples; the LP design and pooling were in four slightly different
    copies;
  * compact per-macro summaries instead of `save.image()` (9.9 GB → 13 MB);
  * `reproduce_paper.R` and `verify_against_paper.R`;
  * per-macro RNG streams; seeds added where they were missing (Figure 2, Figure 3(b));
  * speed-ups (results unchanged in distribution): vectorized option pricing and a vectorized
    bivariate normal CDF (≈140× for barrier prices); matrix form of the LR pooling; one uniform
    per path for barrier survival instead of one per time step; vectorized newsvendor simulator
    (≈67× for the optimal design);
  * **fixed (no effect on the results):** a loop variable that overwrote the strike (`str_k`); the upper-barrier argument
    overwritten by uniforms in `BarrierOptionSim` (no effect, not used); geometric average via
    `prod()` (overflow risk); section numbers and a formula in the code comments; absolute `setwd()` paths;
  * **robustness:** in §7.2 the tuned lpSolveAPI settings fail numerically (status 5) on about 1 in
    5,000 per-stock LPs, when two outer scenarios nearly coincide; the 2024 code then stopped the
    whole batch (`quit()`). `optimal_design()` now retries with lpSolveAPI's default settings and
    then with lpSolve (geometric scaling), which gives the same LP optimum;
  * `gamma_params` hard-coded in the Figure 3(b) script ("QUESTION: where do these come from?") are
    now derived: they are the posterior from `set.seed(2020)`, as in the Table 4 script;
  * **readability:** descriptive names (`single_asset_*`, `multi_asset_*`, `newsvendor_*`,
    `target_ess`, `design$weights`, ...), plain loops instead of `apply`/`vapply`/`sweep`/
    `do.call` constructions, no nested helper functions or `<<-`, comments in the paper's
    notation; `tests/run_tests.R` checks that the library still gives the reference results. The
    rewrite was verified to give identical results to the previous version (library functions
    and small runs of all three simulation scripts with the same seeds);
  * **vectorized newsvendor regression design:** the replications at all design points are
    simulated at once (about 100× faster for this step, which is however only a small part of
    Table 4's run time; it draws its random numbers in a different order than the 2024 code, with
    the same distribution). Table 4 in `results/rerun_full/` was
    rerun with it and re-verified. Vectorizing SNS/SNS+ the same way was tested and gave no
    speed-up (1.0×), so they keep the 2024 order of random draws.

## License

MIT (see [LICENSE](LICENSE)), as the IJOC snapshot. Data: none (all inputs are simulated).
