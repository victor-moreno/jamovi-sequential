# Smoke test, run by tools/install.sh right after installing the module into
# jamovi desktop and/or the docker container (the runner has put jamovi's base
# library and the installed module first on .libPaths() and attached it).
# Reference values: rpact 4.4.0 getSampleSizeSurvival for the same design.
r <- jmsequential::gsd(data = data.frame(), endpoint = "surv", hr = 0.7, ctrlMedian = 12,
                       accrualTime = 24, followUp = 12, dropoutProb = 0.05, dropoutTime = 12,
                       futility = "spending", futSf = "hsd", futGamma = -2)
ev <- r$summary$getCell(rowKey = "events", "value")$value
tm <- r$timing$asDF$timeH1[1]
stopifnot(identical(ev, "353"), abs(tm - 16.62567) < 1e-4)
cat(sprintf("   smoke test passed: %s events, first look at %.2f\n", ev, tm))

# phase II: Simon optimal design for p0 = 0.2, p1 = 0.4, alpha 0.05, power 0.8 (clinfun::ph2simon)
s <- jmsequential::phase2(data = data.frame())$stages$asDF
stopifnot(identical(as.integer(s$n), c(13L, 43L)), identical(as.integer(s$fut), c(3L, 12L)))
cat("   smoke test passed: Simon optimal 3/13, 12/43\n")

# interim monitoring: stop for efficacy at the second look
dat <- data.frame(n = c(57, 114), est = c(0.5, 0.55), se = 2 / sqrt(c(57, 114)))
m <- jmsequential::gsd(data = dat, lookN = "n", estimate = "est", se = "se")
stopifnot(identical(m$monitor$asDF$decision[2], "Stop: reject H0"), m$final$asDF$p < 0.01)
cat("   smoke test passed: monitoring stops at look 2, p =", signif(m$final$asDF$p, 3), "\n")

# SPRT: exact binomial OC (gsDesign::binomialSPRT, maxn 400) and a t-SPRT on data
p <- jmsequential::sprt(data = data.frame())$oc$asDF$pH1
stopifnot(abs(p[1] - 0.03991942) < 1e-6, abs(p[3] - 0.8238091) < 1e-6)
set.seed(40); invisible(rbinom(80, 1, 0.4)); x <- round(rnorm(80) + 0.6 * (rep(1:2, 40) == 2), 2)
t2 <- jmsequential::sprt(data = data.frame(x = x, g = factor(rep(c("C", "T"), 40))), dep = "x", group = "g",
                         model = "t2", nsim = 200)$test$asDF
stopifnot(t2$n == 41)
cat(sprintf("   smoke test passed: SPRT alpha %.4f, t-SPRT decides at n = %d\n", p[1], t2$n))

# MaxSPRT: Poisson critical value for T = 20 (Sequential::CV.Poisson) and a signal on the example data
sv <- data.frame(cases = c(0, 1, 2, 3, 5, 6, 8, 10, 12, 13, 15, 17, 18, 20), expected = round(0.65 * 1:14, 2))
mx <- jmsequential::maxsprt(data = sv, cases = "cases", expected = "expected")
stopifnot(identical(mx$summary$getCell(rowKey = "cv", "value")$value, "3.628123"),
          which(mx$monitor$asDF$decision == "Signal: reject H0")[1] == 12)
cat("   smoke test passed: MaxSPRT cv 3.628123, signal at look 12\n")

# Bayesian monitoring: sceptical prior stops for efficacy at look 2; exact two-arm type I error
r <- jmsequential::bayesmon(data = data.frame(n = c(57, 114), est = c(0.5, 0.55), se = 2 / sqrt(c(57, 114))),
                            lookN = "n", estimate = "est", se = "se")
b <- jmsequential::bayesmon(data = data.frame(), model = "binary")$oc$asDF$reject[1]
stopifnot(identical(r$monitor$asDF$decision[2], "Stop: efficacy"), abs(b - 0.04835856) < 1e-6)
cat(sprintf("   smoke test passed: Bayesian stop at look 2, two-arm alpha %.4f\n", b))

# anytime-valid: exact n for 80 % power with continuous monitoring (safestats: 55 +- 3 by simulation)
a <- jmsequential::anytime(data = data.frame(), sdKnown = TRUE, samples = "one")$summary$asDF$value[4]
stopifnot(identical(a, "54"))
cat("   smoke test passed: anytime-valid n =", a, "\n")
