# jmsequential — sequential designs and monitoring for jamovi

A [jamovi](https://www.jamovi.org) module for studies analysed repeatedly as their data
accumulate — clinical trials with interim analyses, sequential experiments, safety
surveillance — with control of the type I and type II errors.

It needs no packages beyond those bundled with jamovi: `jmvcore`, `R6` and `ggplot2`.

Every analysis works in the same order: the design options come first, and the data go in the
**Analysis of observed data** box at the end of the options panel (variables, or patients and
responses for phase II). Without data, each analysis is a design tool; with data, it also
analyses them against that design.

## Analyses (menu: Sequential)

### Group Sequential Design

*Sequential → Group Sequential Design*.

**Efficacy boundaries**

- Classical: O'Brien–Fleming, Pocock, Wang–Tsiatis (Δ), Haybittle–Peto,
  Pampallona–Tsiatis (with its own futility boundary).
- Error spending: Lan–DeMets O'Brien–Fleming and Pocock types, Kim–DeMets (ρ),
  Hwang–Shih–DeCani (γ), exponential (ν), or user-defined cumulative proportions.
- Whitehead's triangular test, and the double triangular test for two-sided designs.

**Futility boundaries**

- Beta spending, fixed z-values or a conditional-power threshold.
- Binding or non-binding.
- For two-sided designs, an inner wedge.

**Design settings**

- One- or two-sided tests, equal or unequal timing of the analyses, any allocation ratio.

**Endpoints**

- Means, proportions, or a generic effect.
- Time to event: events (Schoenfeld or Freedman), then sample size, accrual duration or
  minimum follow-up with uniform accrual and dropout, and the expected calendar time of each
  analysis under H1 and H0.

**Interim monitoring**

A dataset with one row per analysis holds:

- the cumulative N or events;
- the estimate and its standard error (log hazard ratio for survival), or a Z statistic;
- optionally, the information fraction.

From it the analysis gives:

- the boundaries recomputed at the observed information, and the decision;
- repeated confidence intervals;
- conditional power (under the design effect and under the current trend) and predictive power;
- once a boundary is crossed or the final analysis is reached, the stagewise-ordering p-value,
  median-unbiased estimate and interval.

**Output**

- A design summary, the boundaries (Z, nominal p, effect scale) and the stopping
  probabilities.
- Operating characteristics.
- Plots of the boundaries, the error spending, and power and expected sample size.

### Phase II Multistage Design

*Sequential → Phase II Multistage Design*: single-arm designs for a binary response.

**Designs**

- Simon two-stage: optimal, minimax, or admissible for a chosen weight (Jung et al. 2004).
- Exact single-stage design.
- Fleming's multistage test.
- Bayesian monitoring with Lee & Liu's predictive probability or with posterior probabilities
  (Beta prior; futility, efficacy and final thresholds).

**Output**

- Operating characteristics, computed exactly from the binomial distribution: type I error,
  power, expected sample size, and probability of early stopping.
- For observed data: the decision, the Koyama & Chen p-value and exact interval, Jung & Kim's
  UMVUE after a Simon design, and the posterior and predictive probabilities for the Bayesian
  designs.
- For observed data, also the probability that the trial ends promising given the data so far
  (exact, under p0, p1 and the observed rate), at planned or unplanned looks.
- Stopping boundaries plot: shaded decision regions (not promising, continue, promising) at each
  look and the responses expected under p0 and p1. With observed data it also shows the current
  status, the triangle of trajectories still possible, and the projection at the observed rate
  with its interval (a fan of slopes). Each layer can be switched on or off.
- Outcome map for two-stage designs: every combination of stage-1 and stage-2 responses, coloured
  by the decision it leads to (promising in green). With observed data, the current status is
  outlined in red and the outcomes that can no longer be reached are faded.

### Biomarker-Stratified Phase II

*Sequential → Biomarker-Stratified Phase II*: phase II trials stratified by a biomarker (BM+ /
BM−) that decide both whether to go to phase III and in which population: **BM+ only**
(enrichment), **biomarker-stratified**, **all-comers**, or **stop**.

**Randomized designs** (experimental vs control)

- The decision algorithm of Freidlin et al. (2012): a one-sided test in BM+; if negative, an
  overall stratified test; if positive, the confidence interval of the effect in BM− decides
  between BM+ only, stratified and all-comers. Every constant is editable.
- Time to event (hazard ratio; events per stratum from the power target in BM+, patients from the
  control median per stratum, accrual, follow-up and dropout) or binary response (difference in
  rates).
- Stratum sizes fixed (a BM− : BM+ ratio) or at the natural prevalence of BM+.
- Optional futility look in each stratum at a fraction of the information (an extension: the
  original authors recommend it without specifying it).

**Single-arm designs** (two stages, both strata at stage 1)

- Optimal or minimax stratified design (Parashar et al. 2016): searched to minimise the expected
  sample size under H0, subject to the probability of any go under H0, the power to recommend
  all-comers when both strata respond, and the power to recommend BM+ only when only BM+ responds.
- A user-defined design with the rules of Jones & Holmgren (2007) or the isotonic rules of Li &
  Ivanova (2024), in which each stratum is judged against the pooled rate when the observed rates
  contradict their expected order.
- Parallel Simon designs, one per stratum at α/2.
- Uninteresting and target rates per stratum; for time to event, medians converted to the rate of
  patients event-free at a landmark (exponential survival).
- Stratum sizes fixed, or all-comers enrolment with the cut-offs re-derived at the realised
  stratum sizes (an extension).

**Output**

- Design summary, sample size per stratum, decision rules.
- Operating characteristics: the probability of each recommendation and the expected size under
  no effect, effect in BM+ only, equal and BM+-larger effects (randomized) or both strata active
  (single arm), plus a user-defined scenario. All exact: normal approximation for the randomized
  designs, binomial for the single-arm designs.
- **Decision map:** for randomized designs, the recommendation given the estimated effects in BM+
  and BM−; for single-arm designs, the route after every stage-1 outcome.
- **Recommendation by true effects:** the most likely recommendation over a grid of true effects,
  shaded by its probability.
- **Patient-level data:** biomarker, arm, and response or time and event (in arrival order for
  single-arm trials). Randomized: HR per stratum (Cox) with log-rank tests, the stratified
  overall test, the interaction and Gail & Simon's test for a qualitative interaction, and the
  decision path. Single arm: responses per stage and stratum, the route and the recommendation.

### Sequential Probability Ratio Test

*Sequential → Sequential Probability Ratio Test*: Wald's test after every observation.

**Tests**

- Proportion (one sample); mean with known SD (one sample); difference in means with known SD
  (two samples).
- Sequential one- and two-sample t-tests (Rushton 1950; Hajnal 1961), one- or two-sided.

**Boundaries**

- Wald SPRT, or Whitehead's continuous triangular test.
- Optional truncation at a maximum sample size, with the decision there taken from the sign of
  the LLR or only at the upper threshold.

**Output**

- Thresholds, attained α and power, expected n under H0 and H1, and the fixed-sample n.
- Operating characteristics: exact values (dynamic programming for proportions, the recursive
  integration engine for known-SD means) or seeded simulation (t-tests), next to Wald's
  approximations.
- For data in arrival order (with a group variable for two samples): the decision, the n at the
  decision and the estimate, a path table, and a plot of the path against the boundaries.

### Safety Surveillance (MaxSPRT)

*Sequential → Safety Surveillance (MaxSPRT)*: the maximized SPRT of Kulldorff et al. (2011) for
post-marketing safety surveillance.

- **Models:** Poisson (observed against expected counts), or binomial (cases in the risk window
  against the control window, with matching ratio z).
- **Design:** α, the surveillance length, a minimum number of cases for a signal, and continuous
  or group monitoring (equally spaced or custom looks).
- **Output:** the exact critical value and attained α, the cases needed for a signal, and
  power, expected time to signal and expected length of surveillance for chosen relative risks.
- **Data:** cumulative cases and expected counts (or all cases), one row per look, giving the
  LLR, the relative risk and the signal decision, and a plot of the counts against the
  signalling boundary.

### Bayesian Sequential Monitoring

*Sequential → Bayesian Sequential Monitoring*: Bayesian stopping rules with their frequentist
operating characteristics.

- **Models:** normal approximation (estimate and standard error per look; flat, sceptical,
  enthusiastic or custom normal prior); two arms with a binary outcome (Beta priors); JZS Bayes
  factor t-test on raw data (one or two samples).
- **Rules:** success when P(effect beyond δ | data) reaches a threshold, at every look or only at
  the end; futility when the predictive probability of final success is low; for the Bayes factor,
  stop when BF₁₀ ≥ k or ≤ 1/k between a minimum and a maximum n.
- **Operating characteristics:** exact for the normal model (the rules are thresholds on the
  estimate, handled by the group sequential engine) and for two arms (dynamic programming over the
  responses in both arms); seeded Bayes factor design analysis for the Bayes factor.
- **Monitoring:** posterior mean and interval, P(benefit), P(beyond δ), predictive probability and
  decision at each look, or the Bayes factor path; plots of the thresholds, the probabilities, the
  prior and posterior, and the operating characteristics.

### Anytime-Valid Monitoring

*Sequential → Anytime-Valid Monitoring*: e-values and confidence sequences, valid at every sample
size and under optional stopping.

- **Tests:** mean or difference in means. With a known SD the test is the normal-mixture mSPRT
  (Robbins 1970; Johari et al. 2022); with an unknown SD it is the safe t-test (Grünwald, de Heide
  & Koolen 2024). Both are two-sided, with a Gaussian mixing variance (default δ²).
- **Design:** the N for a target power when testing after every observation (exact for a known
  SD, simulated otherwise), the expected and median stopping N, the type I error by the horizon,
  and the N for one analysis with the e-value and for the fixed-sample test.
- **Data:** outcome in arrival order, with a group variable for two samples. The analysis gives
  the e-value, the always-valid p-value, the confidence sequence and the decision, with plots of
  the e-process and of the confidence sequence.

## Example data

*Data library → jmsequential*: interim results of a trial (difference in means, survival, two arms
binary), patients in arrival order (binary and continuous outcomes, two arms), weekly
surveillance counts, and a randomized and a single-arm biomarker-stratified phase II trial.

## Tutorial

A [website on sequential methods](https://victor-moreno.github.io/jamovi-sequential/) (in
`docs/`): boundaries, spending, futility, triangular tests, sample size, monitoring, inference,
SPRT, MaxSPRT, phase II, biomarker-stratified phase II, Bayesian and anytime-valid designs, with
examples computed by this
module. Each method page ends with its example worked in jamovi, with a screenshot of the options
and the results. The site is rendered with Quarto from a source kept outside this repository.

## Install

Prebuilt `.jmo` files are attached to the [Releases](../../releases) page — four per version, one
for each combination of jamovi series and CPU:

| your jamovi | Apple silicon | Intel / AMD |
| --- | --- | --- |
| **current** (bundles R 4.6.0) | `jmsequential_<version>_current_R4.6.0_arm64.jmo` | `jmsequential_<version>_current_R4.6.0_x64.jmo` |
| **solid** (bundles R 4.5.0) | `jmsequential_<version>_solid_R4.5.0_arm64.jmo` | `jmsequential_<version>_solid_R4.5.0_x64.jmo` |

The same file works on macOS, Windows and Linux: jamovi checks the R version and the CPU, not the
operating system (see Help → About). In jamovi, go to Modules (⊕) → Sideload and choose the
`.jmo`.

## Build and install (developers)

```sh
bash tools/install.sh                 # desktop and Docker (container "jamovi")
bash tools/install.sh desktop         # jamovi.app, via jmvtools::install()
bash tools/install.sh docker [name]   # jmc --install in a running container
bash tools/build-jmo.sh current       # release: dist/*_current_R4.6.0_{x64,arm64}.jmo
bash tools/build-jmo.sh solid         # release: dist/*_solid_R4.5.0_{x64,arm64}.jmo
bash tools/release.sh                 # build all four and publish one GitHub release
```

After each install, `tools/smoke.R` checks the installed module (a survival sample size and a
Simon design). The tests run with `Rscript -e 'testthat::test_local(".")'`; their reference
values in `tests/testthat/fixtures/` were generated from the R packages named below, which are
not needed to run them.

## Notes on methods

- The numerical engine is pure R: the recursive integration of Armitage, McPherson & Rowe (1969)
  as described by Jennison & Turnbull (2000, ch. 19), used for the group sequential boundaries,
  the known-SD SPRT, the normal Bayesian rules and the known-SD anytime-valid design.
- Phase II designs, the binomial SPRT and the two-arm Bayesian rules are computed exactly from the
  binomial distribution.
- Biomarker-stratified designs: randomized probabilities from independent normal statistics per
  stratum, with one Gauss–Legendre integral for the overall test and the bivariate normal
  distribution (Genz 2004) for the interim look; the single-arm search enumerates the BM− and BM+
  parts of the design separately and combines them by branch and bound.

## Acknowledgements

This module was designed and built together with Claude Code, over several rounds of design,
implementation and testing in jamovi.

Methods follow the publications cited in the tutorial, in particular:

- Jennison & Turnbull, *Group Sequential Methods with Applications to Clinical Trials* (2000);
- Whitehead, *The Design and Analysis of Sequential Clinical Trials* (1997);
- Lan & DeMets (1983);
- Lachin & Foulkes (1986);
- Wald, *Sequential Analysis* (1947);
- Kulldorff, Davis, Kolczak, Lewis, Lieu & Platt (2011);
- Spiegelhalter, Freedman & Parmar (1994); Rouder, Speckman, Sun, Morey & Iverson (2009);
- Robbins (1970); Grünwald, de Heide & Koolen (2024);
- Freidlin, McShane, Polley & Korn (2012); Jones & Holmgren (2007); Parashar, Bowden, Starr,
  Wernisch & Mander (2016); Li & Ivanova (2024).

Ideas and reference values for the tests come from these R packages:
[rpact](https://www.rpact.org) (G. Wassmer, F. Pahlke),
[gsDesign](https://keaven.github.io/gsDesign/) (K. Anderson), ldbounds (C. Casper, T. Cook,
O. A. Perez), seqtest (T. Yanagida), clinfun (V. E. Seshan), ph2bayes (K. Nagashima), SPRT
(H. Budihal), sprtt (M. Snijder-Steinhilber, M. Schnuerch, A.-L. Schubert), Sequential (I. R. Silva,
M. Kulldorff), BayesFactor (R. D. Morey, J. N. Rouder) and
[safestats](https://github.com/AlexanderLyNL/safestats) (A. Ly, R. Turner and co-authors), ph2hetero
(P. Sfumato, T. Filleron, B. Cabarrou), mvtnorm (A. Genz, F. Bretz, T. Hothorn and co-authors) and
survival (T. Therneau), and from the NCI program RP2BM of Freidlin et al. They are
not dependencies of the module, and no code was copied from packages whose licences are
incompatible with GPL-3.

This module has been developed with support of the Instituto de Salud Carlos III (ISCIII),
“Programa FORTALECE del Ministerio de Ciencia e Innovación”, through the project number
FORT23/00032 and the Consortium for Biomedical Research in Epidemiology and Public Health
(CIBERESP), action Genrisk.

## Licence

GPL (≥ 3).
