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
