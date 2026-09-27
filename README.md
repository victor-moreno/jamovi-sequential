# jmsequential — sequential designs and monitoring for jamovi

A [jamovi](https://www.jamovi.org) module for studies analysed repeatedly as their data
accumulate — clinical trials with interim analyses, sequential experiments, safety
surveillance — with control of the type I and type II errors.

## Version 0.1.0

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

**Validation**

- The Simon optimal and minimax designs agree with clinfun in 56 scenarios.
- The UMVUE and p-values agree with clinfun's `twostage.inference`.
- The predictive probabilities and boundaries agree with ph2bayes.
- Fleming's boundaries reproduce two published examples.

Planned next: interim monitoring and inference after stopping; SPRT and MaxSPRT; Bayesian and
anytime-valid monitoring; a classical fixed-sample size module.

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
Bayesian and anytime-valid designs — with examples computed by this module.
See `tutorial/README.md` to build it.

## Development

- **Tests:** `Rscript -e 'testthat::test_local(".")'`. The reference values under
  `tests/testthat/fixtures/` were generated from rpact and gsDesign.
- **Build:** `bash tools/build-jmo.sh {current|solid}` writes `dist/`.
- **Install:** `bash tools/install.sh [desktop|docker|both]`, which runs `tools/smoke.R` against
  the installed module.

## Acknowledgements

Methods follow the publications cited in the tutorial, in particular:

- Jennison & Turnbull, *Group Sequential Methods with Applications to Clinical Trials* (2000);
- Whitehead, *The Design and Analysis of Sequential Clinical Trials* (1997);
- Lan & DeMets (1983);
- Lachin & Foulkes (1986).

The R packages [rpact](https://www.rpact.org) (G. Wassmer, F. Pahlke),
[gsDesign](https://keaven.github.io/gsDesign/) (K. Anderson), ldbounds and seqtest served as
references and test oracles. No code was copied from packages whose licences are incompatible
with GPL-3.

## Licence

GPL (≥ 3).
