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
