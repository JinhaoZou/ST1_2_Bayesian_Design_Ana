# BayesianDesignAna

**Bayesian design and analysis for multi-arm trials with binary endpoints**

Version 0.6.0 · R package repository: **ST1_1_Bayesian_Design_Ana**

`BayesianDesignAna` provides Bayesian operating-characteristic simulation, sample-size planning with exact allocation ratios, optional historical-control borrowing, MAP prior effective sample size (ESS) assessment, and Bayesian analysis of binary endpoints. The common-control design functions support one or more active doses. Functions for patient-level Bayesian logistic regression support zero or more stratification factors, including categorical factors with multiple levels.

**Scope:** Development software. Its simulation-based power engine uses a Beta-binomial *posterior analysis* with fixed binomial response rates as the default data generator. The logistic analysis uses a different, covariate-adjusted model. Neither familywise error control nor numerical equivalence with other software should be assumed without independent validation. Outputs are intended for research and method development, not unvalidated confirmatory decisions.

## Installation

**Install the package from GitHub (after the repository is public):**

```r
install.packages("remotes")
remotes::install_github("JinhaoZou/ST1_1_Bayesian_Design_Ana", upgrade="never")
library(BayesianDesignAna)
packageVersion("BayesianDesignAna")
```

**Install locally:** Download this repository as ZIP from GitHub and extract it. In RStudio, replace the path with the extracted folder containing `DESCRIPTION` and `R/`:

```r
install.packages("remotes")
remotes::install_local("~/Downloads/BayesianDesignAna", upgrade="never", force=TRUE)
library(BayesianDesignAna)
packageVersion("BayesianDesignAna")
```

The GitHub **repository** name and the installed R **package** name are different intentionally. Use `JinhaoZou/ST1_1_Bayesian_Design_Ana` in `install_github()` and `BayesianDesignAna` in `library()`. A browser-only publishing walkthrough is available in [GITHUB_UPLOAD_GUIDE.md](GITHUB_UPLOAD_GUIDE.md).

To test a local source folder:

```r
install.packages("devtools")
devtools::test("~/Downloads/BayesianDesignAna")
devtools::check("~/Downloads/BayesianDesignAna")
```

`historical_map_ess()` requires `RBesT`, which is an optional dependency:

```r
install.packages("RBesT")
```

## 1. Exact allocation ratios

Each arm gets an **integer** sample size with an exact allocation ratio. All ratios must be named positive integer weights, with placebo/control first.

```r
r1 <- c(Placebo=1, Dose1=1)
allocation_unit(r1)                         # 2
candidate_totals(r1, 80, 100)              # 80, 82, ..., 100
allocate_arms(84, r1)                      # Placebo=42, Dose1=42

r4 <- c(Placebo=1, Low=1, Mid=1, High=1)
allocation_unit(r4)                         # 4
candidate_totals(r4, 160, 176)            # 160, 164, 168, 172, 176
allocate_arms(164, r4)                     # 41, 41, 41, 41
# allocate_arms(165, r4)                   # Error: off-ratio total

r_unequal <- c(Placebo=1, Low=2, Mid=2, High=1)
candidate_totals(r_unequal, 120, 144)      # 120, 126, ..., 144
allocate_arms(126, r_unequal)             # 21, 42, 42, 21
```

Use `step_blocks=5` in `candidate_totals()` to screen every fifth *valid allocation block* (rather than every fifth person):

```r
candidate_totals(r4, 160, 240, step_blocks=5) # 160, 180, 200, 220, 240
```

## 2. Power from specified arm sizes

```r
n <- c(Placebo=41, Dose1=41, Dose2=41)
rates <- c(Placebo=.30, Dose1=.50, Dose2=.50)

power <- design_power(
  n=n, response_rates=rates,
  borrowed_n=41, historical_rate=.30,
  overall="any", multiplicity="hochberg",
  alpha=.025, margin=0, seed=2026,
  n_sim=1000, n_post=2000,
  base_control_prior=c(.5,.5), treatment_prior=c(.5,.5)
)
power$per_arm       # Individual dose power and Monte Carlo SE
power$overall_power # Probability >= 1 dose succeeds, for overall='any'
power$scenario
```

`overall="any"` (default) requires **at least one** dose to succeed; `overall="all"` requires **every** dose. The default multiplicity convention is `"hochberg"`; alternatives are `"none"`, `"bonferroni"`, `"holm"`, and `"BH"`. These methods adjust Bayesian posterior-tail decision scores, not frequentist p-values, and do not automatically establish conventional FWER or FDR control.

Two ways to specify borrowing in **design** functions (choose only one):

- `borrowed_n=20`: fixed borrowed ESS.
- `borrowing_prop=.40`: borrowed ESS / (borrowed ESS + concurrent placebo N). With 30 concurrent controls, ESS is 20.

Historical ESS calculations are **independent**; they are not called by design simulation functions.

## 3. Sample-size planning

Full grid: every candidate is evaluated at the same simulation precision.

```r
ratios <- c(Placebo=1, Low=1, Mid=1, High=1)
rates <- c(Placebo=.30, Low=.40, Mid=.50, High=.60)

ss <- design_sample_size(
  total_n=candidate_totals(ratios, 120, 240),
  ratios=ratios, response_rates=rates,
  borrowing_prop=.50, target_power=.80,
  search="grid", overall="any", multiplicity="hochberg",
  seed=2026, n_sim=1000, n_post=2000, plot=TRUE
)
ss$results
ss$minimum_candidate
plot(ss)  # Plot existing results; no extra simulation
```

**Two-stage search:** First screen a coarser set of feasible totals using fewer simulations; then run more simulations at neighboring **valid totals** around the preliminary crossing. For 1:1:1:1, only totals such as 164 and 168 are possible; never 165, 166 or 167.

```r
ss2 <- design_sample_size(
  total_n=candidate_totals(ratios, 120, 240),
  ratios=ratios, response_rates=rates,
  borrowing_prop=.5, target_power=.80,
  search="two_stage",
  initial_step=5, n_sim_initial=500,
  refine_span=5, n_sim=10000,
  n_post=2000, overall="any", multiplicity="hochberg",
  seed=2026, plot=TRUE
)
ss2$initial_results
ss2$final_results
ss2$results
ss2$minimum_candidate
ss2$note
```

The selected sample size is the first **final-stage** candidate whose *estimated* overall power reaches the target. It is not a proof of the true minimum; consider confirmatory simulations and Monte Carlo uncertainty. A `historical_ess` value, if supplied, is only a comparison threshold: it does **not** change the sample-size calculation, but the output notes when proposed borrowing exceeds it.

## 4. Power curves for multiple dose-effect scenarios

Each row of `delta_matrix` is one combination of absolute risk differences relative to placebo; each row produces one curve. The plot has total randomized N on the x-axis and a table beneath with one row per arm, borrowed ESS, and scenario power. It calls `design_sample_size()` to calculate results; no further simulations are needed to draw the plot.

```r
ratios <- c(Placebo=1, Low=2, Mid=2, High=1)
deltas <- rbind(
  c(.05, .10, .20),
  c(.10, .20, .30),
  c(.15, .25, .35)
)
colnames(deltas) <- c("Low", "Mid", "High")

curves <- plot_power(
  delta_matrix=deltas,
  total_n=candidate_totals(ratios, 120, 240),
  ratios=ratios, placebo_rate=.30,
  borrowed_n=15, target_power=.80,
  search="two_stage", initial_step=3, refine_span=3,
  n_sim_initial=300, n_sim=3000,
  overall="any", multiplicity="hochberg", seed=2026,
  table_max_cols=10, plot=TRUE
)
curves$results
```

## 5. Borrowing and independent historical MAP ESS

```r
borrowing_ess(n_control=30, borrowing_prop=.40) # 20
borrowing_prior(historical_rate=.30, ess=20, base_prior=c(.5,.5))
borrowing_limits(n_control=30, historical_ess=50,
  max_fraction=.40, max_absolute=25, max_ratio=2)
```

Use `historical_map_ess()` separately with historical controls and `RBesT` to estimate ELIR, moment and Morita ESS:

```r
historical <- data.frame(
  study=paste0("Study",1:6),
  n=c(154,163,183,180,72,74),
  r=c(40,45,57,61,21,24)
)
map <- historical_map_ess(historical, seed=2026)
map$ess
reference_ess <- min(map$ess)

# Optional reference only: does not change trial simulations
ss_ref <- design_sample_size(
  total_n=candidate_totals(c(Placebo=1,Dose1=1), 60, 120),
  ratios=c(Placebo=1,Dose1=1),
  response_rates=c(Placebo=.30,Dose1=.50),
  borrowed_n=20, historical_ess=reference_ess,
  n_sim=200, n_post=500, plot=FALSE
)
ss_ref$note
```

The minimum across the three ESS estimates is a **user-selected conservative reference**, not a validated borrowing cap. MAP prior effective sample size and MCMC effective sample size are different concepts.

## 6. Analyze binary responder counts

```r
fit_bin <- analyze_binary(
  responses=c(Placebo=12, Low=18, High=25),
  n=c(Placebo=40, Low=40, High=40),
  control_prior=c(.5,.5), treatment_prior=c(.5,.5),
  margin=0, alpha=.025, multiplicity="hochberg",
  n_post=4000, seed=2026
)
fit_bin
```

For advanced use of the core simulator and legacy search helper:

```r
simulate_power(c(Placebo=40,Dose1=40),c(Placebo=.30,Dose1=.50),
  borrow_ess=20,n_sim=200,n_post=500,seed=2026)

find_sample_size(candidate_totals(c(Placebo=1,Dose1=1),40,100),
  ratios=c(Placebo=1,Dose1=1),
  response_rates=c(Placebo=.30,Dose1=.50),target="any",
  n_sim=100,n_post=500)
```

For new sample-size studies, prefer `design_sample_size()`.

## 7. Simulated patient-level data and Bayesian logistic regression

The following example uses two categorical stratification factors, each with three levels. Each factor with K levels contributes K−1 indicator coefficients.

```r
n <- c(Placebo=40, Low=40, Mid=60, High=60)
rates <- c(Placebo=.30, Low=.40, Mid=.50, High=.60)

dat <- simulate_logistic_data(
  n, rates,
  strata=list(region=c(North=.4,South=.35,West=.25),
              severity=c(Mild=.65,Moderate=.25,Severe=.10)),
  stratum_log_or=list(region=c(.15,-.20),severity=c(.30,.55)),
  seed=2026
)
table(dat$arm, dat$region)

fit <- analyze_logistic(
  data=dat, outcome="y", arm="arm", control="Placebo",
  covariates=c("region","severity"), center="overall",
  control_prior=c(.5,.5), intercept_sd=100, covariate_sd=100,
  interval_method="both", credible_level=.95,
  sas_compatible_summary=TRUE,
  margin=0, alpha=.025, multiplicity="hochberg",
  n_iter=12000, burnin=2000, thin=5, n_chains=2,
  proposal_sd=.15, seed=2026
)
fit$summary              # RD posterior means and both interval types
fit$reference_means      # Arm response probabilities at centered reference profile
fit$sas_summary          # One row, D1/D2/etc. named columns (HPD bounds)
fit$mcmc                 # Approximate MCMC R-hat and ESS
fit$acceptance           # Proposal acceptance fraction
print(fit)
mcmc_diagnostics(fit$chains)
```

**Credible interval reporting:**

- `interval_method="equal_tail"` computes the lower and upper quantiles with equal posterior tail probabilities.
- `interval_method="HPD"` computes the shortest *contiguous* interval of retained posterior draws containing approximately the requested probability; for unimodal posteriors it approximates the HPD interval.
- `interval_method="both"` (default) returns both named interval pairs; `rd_lower`/`rd_upper` remain equal-tailed for compatibility.
- `credible_level=.95` determines interval mass. Interval choice does **not** change posterior means or superiority probabilities.
- `sas_compatible_summary=TRUE` produces one-row SAS-style output. It uses HPD-style limits (`RD_D1_LCL`, `RD_D1_UCL`, etc.) and `Pr_Sup_D1 = Pr(RD_D1 > 0 | data)` regardless of `margin`. The D1/D2 labels refer to the order of active doses in `fit$arms`; verify this order when comparing external results.

Use an explicit prior when historical data justify one. Example:

```r
fit_historical <- analyze_logistic(
  dat, covariates=c("region","severity"), center="overall",
  control_prior=c(5.52,14.48),
  interval_method="HPD", credible_level=.95,
  sas_compatible_summary=TRUE, multiplicity="none",
  n_iter=20000,burnin=5000,n_chains=4,seed=2026
)
fit_historical$sas_summary
```

For one two-level factor, use `strata=list(severity=c(Mild=.6,Severe=.4))`, `stratum_log_or=list(severity=.5)`, and `covariates="severity"`. With no covariates, omit `strata` and set `covariates=character()`. `center="overall"` evaluates arm differences at a common centered covariate profile; `center="arm"` uses different arm-specific reference profiles and should not be interpreted as a common-profile contrast. Centered reference effects are not necessarily population-averaged marginal effects.

## 8. Function options and defaults

The tables summarize the exported functions. Arguments marked **required** have no default. Use `?function_name` in R for complete details.

| Function | Inputs and options |
|---|---|
| `allocation_unit(ratios)` | `ratios` **required**: named positive integer allocation weights, control first; returns allocation block size. |
| `candidate_totals(ratios, min_total, max_total, step_blocks=1)` | **Required:** weights and inclusive total range; `step_blocks` advances by multiples of the allocation block. |
| `allocate_arms(total_n, ratios=c(Control=1,Dose1=1))` | Exact randomized total and arm weights; off-ratio totals are rejected. |
| `borrowing_ess(n_control, borrowing_prop=0)` | Concurrent control N and fraction of combined nominal control information borrowed, in `[0,1)`. |
| `borrowing_prior(historical_rate, ess=0, base_prior=c(.5,.5))` | Historical control rate, borrowed pseudo-subject ESS, baseline Beta shapes. |
| `borrowing_limits(n_control, historical_ess, max_fraction=.5, max_absolute=Inf, max_ratio=Inf)` | Available historical ESS and optional caps; returns reference limits, not an automatically validated borrowing rule. |
| `historical_map_ess(data, seed=123, ...)` | **Required:** data with `study`, `n`, `r`; `...` forwarded to `RBesT::gMAP()`; returns MAP model, mixture and three ESS estimates. |
| `analyze_binary(responses, n, control_prior=c(.5,.5), treatment_prior=c(.5,.5), n_post=4000, margin=0, alpha=.025, multiplicity="hochberg", seed=NULL)` | Named responder counts and arm sizes; separate control/treatment priors, superiority margin, posterior-tail decision convention, optional seed. |
| `simulate_power(n, response_rates, borrow_ess=0, historical_rate=NULL, n_sim=500, n_post=2000, multiplicity="hochberg", alpha=.025, margin=0, seed=123, base_control_prior=c(.5,.5), treatment_prior=c(.5,.5))` | Low-level simulation; `historical_rate=NULL` uses assumed control rate for borrowing; trial replicates and posterior draws set Monte Carlo accuracy. |
| `design_power(n, response_rates, borrowed_n=NULL, borrowing_prop=NULL, historical_rate=NULL, overall=c("any","all"), multiplicity="hochberg", seed=123, base_control_prior=c(.5,.5), treatment_prior=c(.5,.5), ...)` | Power for a fixed exact named allocation; specify at most one borrowing option; `...`: `n_sim`, `n_post`, `alpha`, `margin`. |
| `design_sample_size(total_n, ratios, response_rates, target_power=.8, borrowed_n=NULL, borrowing_prop=NULL, historical_ess=NULL, plot=FALSE, overall=c("any","all"), multiplicity="hochberg", seed=123, base_control_prior=c(.5,.5), treatment_prior=c(.5,.5), search="grid", n_sim=1000, n_sim_initial=500, initial_step=5, refine_span=NULL, table_max_cols=12, ...)` | Evaluates feasible sample-size candidates. `search="two_stage"` screens and refines; `refine_span=NULL` uses `initial_step`. `historical_ess` is warning/reference only; `...`: `n_post`, `historical_rate`, `alpha`, `margin`. |
| `plot.bda_sample_size(x, ...)` / `plot(ss)` | Displays cached power estimates and allocation table, no new simulation; x is a result from `design_sample_size()`. |
| `plot_power(delta_matrix, total_n, ratios, placebo_rate, borrowed_n=NULL, borrowing_prop=NULL, historical_ess=NULL, plot=TRUE, overall=c("any","all"), multiplicity="hochberg", seed=123, base_control_prior=c(.5,.5), treatment_prior=c(.5,.5), table_max_cols=12, ...)` | One line per delta-vector row. `...` forwards sample-size search controls (`target_power`, `search`, `initial_step`, `refine_span`, `n_sim_initial`, `n_sim`, `n_post`) and simulation settings. |
| `find_sample_size(total_n, ratios, response_rates, target_power=.8, target="any", ...)` | Simple/legacy search; `target` can be `any`, `all` or a dose name; `...` is forwarded to `simulate_power()`. Prefer `design_sample_size()`. |
| `simulate_logistic_data(n, response_rates, strata=list(), stratum_log_or=list(), seed=123)` | Named arm sizes/rates, optional list of category probabilities and log-odds effects for nonreference levels, seed. |
| `analyze_logistic(data, outcome="y", arm="arm", control="Placebo", covariates=character(), center="overall", control_prior=c(.5,.5), intercept_sd=100, covariate_sd=100, n_iter=12000, burnin=2000, thin=5, n_chains=2, proposal_sd=.15, margin=0, alpha=.025, multiplicity="hochberg", seed=123, interval_method="both", credible_level=.95, sas_compatible_summary=FALSE)` | Patient-level logistic analysis. `center` is `overall`, `arm`, `none`. Credible-interval options: `both`, `HPD`, `equal_tail`. `sas_compatible_summary=TRUE` returns a one-row wide results table. |
| `mcmc_diagnostics(chains)` | List of retained MCMC parameter matrices; approximate split R-hat and ESS. |
| `print.bbd_logistic(x, ...)` / `print(fit)` | Pretty-prints the logistic analysis results; pass an analysis object. |

## 9. Interpretation and software status

- **Priors:** Default binary response priors are Beta(.5,.5); logistic control reference prior is Beta(.5,.5) and active intercept/covariate coefficient priors are diffuse proper Normal(0,100²). Historical borrowing is not automatic.
- **Multiplicity:** Standard adjustment names applied to Bayesian posterior-tail scores do **not** imply corresponding frequentist error control; evaluate global and partial-null operating characteristics.
- **Borrowing:** Borrowed ESS is nominal pseudo-patient information. `historical_map_ess()` is separate from the design and analysis functions. Prior-data conflict can reduce the appropriateness of borrowing.
- **MCMC:** `analyze_logistic()` uses a custom component-wise Metropolis sampler; assess convergence, sensitivity to proposals and priors, and numerical agreement on shared datasets. HPD intervals from MCMC are empirical approximations, not exact posterior regions for arbitrary multimodal densities.
- **Power:** Monte Carlo uncertainty may move a power curve up or down. A two-stage search is a computational heuristic and should be confirmed near the selected candidate.
- **Evidence:** The R code and unit tests require local execution and review. A successful software package check does not establish statistical validation.

## 10. References and acknowledgement

The independent historical MAP ESS functionality uses the **R Bayesian Evidence Synthesis Tools (`RBesT`)** package by Sebastian Weber and collaborators, including `gMAP()`, `automixfit()`, and `ess()`:

Weber S, Li Y, Seaman JW, Kakizume T, Schmidli H (2021). *Applying Meta-Analytic-Predictive Priors with the R Bayesian Evidence Synthesis Tools.* **Journal of Statistical Software**, 100(19), 1–32. https://doi.org/10.18637/jss.v100.i19

- RBesT project and reference: https://opensource.nibr.com/RBesT/
- RBesT binary tutorial: https://opensource.nibr.com/RBesT/articles/introduction.html

RBesT is an optional dependency for `historical_map_ess()` only. Results from the MAP function are not passed automatically into power calculations or logistic analysis.

## License and reporting issues

MIT license; see `LICENSE.md`. For questions and software issues, use the repository's [GitHub Issues page](https://github.com/JinhaoZou/ST1_1_Bayesian_Design_Ana/issues). Before a public release, the maintainer must replace the placeholder email in `DESCRIPTION`, confirm ownership of any incorporated work, run tests, and validate the statistical outputs.
