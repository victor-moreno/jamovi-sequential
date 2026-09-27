# jmsequential — sequential designs and monitoring for jamovi

A [jamovi](https://www.jamovi.org) module for studies analysed repeatedly as their data
accumulate — clinical trials with interim analyses, sequential experiments, safety
surveillance — with control of the type I and type II errors.

## Version 0.5.10

Every analysis works in the same order: the design options come first, and the data go in the
**Analysis of observed data** box at the end of the options panel (variables, or patients and
responses for phase II). Without data, each analysis is a design tool; with data, it also
analyses them against that design.

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

For one-sided designs these results agree with rpact's `getAnalysisResults` to 10⁻⁶.

**Output**

- A design summary, the boundaries (Z, nominal p, effect scale) and the stopping
  probabilities.
- Operating characteristics.
- Plots of the boundaries, the error spending, and power and expected sample size.

The numerical engine is pure R: the recursive integration of Armitage, McPherson & Rowe
(1969) as described by Jennison & Turnbull (2000, ch. 19). Its results agree with
**rpact 4.4.0** and **gsDesign 3.11.0** to within 10⁻⁵ over the 144 designs in the test suite.
Sample sizes, events, accrual, follow-up and analysis times agree with rpact to within 10⁻⁶
(relative).

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

**Validation**

- The Simon optimal and minimax designs agree with clinfun in 56 scenarios.
- The UMVUE and p-values agree with clinfun's `twostage.inference`.
- The predictive probabilities and boundaries agree with ph2bayes.
- Fleming's boundaries reproduce two published examples.

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

**Validation:** the LLR paths agree with the SPRT package and with sprtt; the exact binomial
operating characteristics agree with gsDesign's `binomialSPRT` to 10⁻⁷.

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

**Validation:** critical values, attained α, power, time to signal and surveillance length agree
with the Sequential package (Silva & Kulldorff) to 10⁻⁶ in 10 designs. Sequential reports the
lower end of the interval of critical values that give the same discrete design; so does
jmsequential.

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

**Validation:** the Bayes factor agrees with BayesFactor's `ttest.tstat` (to 2·10⁻⁵ on the log
scale, BayesFactor's own accuracy; to 10⁻⁹ against adaptive integration). The two-arm posterior
probability agrees with the exact formula to 10⁻¹². The exact operating characteristics agree with
simulation.

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

**Validation:** the e-values equal safestats' `eGauss` e-values; the single-analysis sample sizes
are identical to safestats' designs, and the continuous-monitoring sample sizes agree within
simulation error.

### Example data

*Data library → jmsequential*: interim results of a trial (difference in means, survival, two arms
binary), patients in arrival order (binary and continuous outcomes, two arms), and weekly
surveillance counts.

Planned next: individual-level data for interim analyses; a classical fixed-sample size module.

## Installation

Download the `.jmo` that matches your jamovi (Help → About shows its R version) from the
Releases page, then *Modules → jamovi library → Sideload*:

| your jamovi | Apple silicon | Intel / AMD |
| --- | --- | --- |
| current (R 4.6.0) | `jmsequential_<version>_current_R4.6.0_arm64.jmo` | `…_current_R4.6.0_x64.jmo` |
| solid (R 4.5.0) | `jmsequential_<version>_solid_R4.5.0_arm64.jmo` | `…_solid_R4.5.0_x64.jmo` |

## Tutorial

`tutorial/` is a Quarto website that introduces sequential methods — boundaries, spending,
futility, triangular tests, sample size, monitoring, inference, SPRT, MaxSPRT, phase II,
Bayesian and anytime-valid designs — with examples computed by this module. Each method page
ends with its example worked in jamovi, with a screenshot of the options and the results
(captured by `tutorial/screenshots/capture.py`). See `tutorial/README.md` to build it.

## Development

- **Tests:** `Rscript -e 'testthat::test_local(".")'`. The reference values under
  `tests/testthat/fixtures/` were generated from rpact, gsDesign, clinfun, ph2bayes, SPRT, sprtt,
  Sequential, BayesFactor and safestats.
- **Build:** `bash tools/build-jmo.sh {current|solid}` writes `dist/`.
- **Install:** `bash tools/install.sh [desktop|docker|both]`, which runs `tools/smoke.R` against
  the installed module.

## Acknowledgements

Methods follow the publications cited in the tutorial, in particular:

- Jennison & Turnbull, *Group Sequential Methods with Applications to Clinical Trials* (2000);
- Whitehead, *The Design and Analysis of Sequential Clinical Trials* (1997);
- Lan & DeMets (1983);
- Lachin & Foulkes (1986);
- Wald, *Sequential Analysis* (1947);
- Kulldorff, Davis, Kolczak, Lewis, Lieu & Platt (2011);
- Spiegelhalter, Freedman & Parmar (1994); Rouder, Speckman, Sun, Morey & Iverson (2009);
- Robbins (1970); Grünwald, de Heide & Koolen (2024).

The R packages [rpact](https://www.rpact.org) (G. Wassmer, F. Pahlke),
[gsDesign](https://keaven.github.io/gsDesign/) (K. Anderson), ldbounds, seqtest, clinfun,
ph2bayes, SPRT, sprtt, Sequential, BayesFactor and
[safestats](https://github.com/AlexanderLyNL/safestats) (A. Ly, R. Turner, J. ter Schure) served as
references and test oracles. No code was copied from packages whose licences are incompatible
with GPL-3.

## Licence

GPL (≥ 3).
