fx <- readRDS(test_path("fixtures", "bayes_oracle.rds"))

test_that("Gaussian-mixture e-values match safestats (eGauss, two-sided)", {
    ev <- fx$evalue
    two <- !is.na(ev$n2)
    nEff <- ifelse(two, ev$n1 * ev$n2 / (ev$n1 + ev$n2), ev$n1)
    nu <- ifelse(two, ev$n1 + ev$n2 - 2, ev$n1 - 1)
    expect_equal(exp(av_logE_z(ev$stat, nEff, ev$g)), ev$ez, tolerance = 1e-10)
    expect_equal(exp(av_logE_t(ev$stat, nEff, nu, ev$g)), ev$et, tolerance = 1e-10)
})

test_that("confidence sequences invert the e-value test", {
    set.seed(5); x <- stats::rnorm(40, 0.3)
    a <- av_analyse(x, NULL, 0, sigma = NA, g = 0.25)
    n <- 30; s <- stats::sd(x[1:n])
    for (lim in c(a$lower[n], a$upper[n]))
        expect_equal(exp(av_logE_t((mean(x[1:n]) - lim) / (s / sqrt(n)), n, n - 1, 0.25)), 20, tolerance = 1e-8)
    z <- av_analyse(x, NULL, 0, sigma = 1, g = 0.25)
    expect_equal(exp(av_logE_z(sqrt(n) * (mean(x[1:n]) - z$upper[n]), n, 0.25)), 20, tolerance = 1e-8)
    expect_true(all(diff(a$p) <= 0))
})

test_that("design: exact n for power agrees with safestats' simulation; single-look n is identical", {
    d <- av_design(TRUE, FALSE, 0.5, 0.05, 0.8)
    expect_lt(abs(d$nPlan - fx$design$z1[["nPlan"]]), 3 * fx$design$z1[["se"]] + 1)
    expect_equal(d$nSingle, fx$design$z1[["batch"]])
    d2 <- av_design(TRUE, TRUE, 0.5, 0.05, 0.8)
    # safestats reports n per group, jmsequential the total
    expect_lt(abs(d2$nPlan / 2 - fx$design$z2[["nPlan"]]), 3 * fx$design$z2[["se"]] + 1)
    expect_equal(d2$nSingle / 2, fx$design$z2[["batch"]])
    dt <- av_design(FALSE, FALSE, 0.5, 0.05, 0.8, nsim = 1000)
    expect_lt(abs(dt$nPlan - fx$design$t1[["nPlan"]]), 3 * fx$design$t1[["se"]] + 3)
    expect_equal(dt$nSingle, fx$design$t1[["batch"]])
    expect_lt(d$alpha, 0.05)
})

test_that("the anytime analysis runs and rejects on the example data", {
    sq <- read.csv(test_path("..", "..", "data", "sprt_sequential.csv"))
    sq$arm <- factor(sq$arm)
    r <- anytime(data = sq, dep = "score", group = "arm")
    m <- r$monitor$asDF
    expect_equal(m$n[m$decision == "Reject H0 (E ≥ 1/α)"], 73)
    r <- anytime(data = data.frame(), sdKnown = TRUE, samples = "one")
    expect_equal(r$summary$asDF$value[4], "54")
})
