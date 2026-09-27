skip_if_not_installed("jmvcore")

run_gsd <- function(...) gsd(data = data.frame(), ...)
cell <- function(tbl, key, col) tbl$getCell(rowKey = key, col)$value

test_that("the default design (means) matches rpact's normal-approximation sample size", {
    r <- run_gsd()
    b <- r$bounds$asDF
    expect_equal(nrow(b), 3)
    expect_equal(b$effZ, c(3.710303, 2.511427, 1.993048), tolerance = 1e-5)
    # rpact getSampleSizeMeans(normalApproximation = TRUE): 170.1114 for delta 0.5, sd 1
    expect_equal(cell(r$summary, "n", "value"), "171")
    expect_equal(nrow(r$oc$asDF), 4)
})

test_that("survival designs fill the timing table and match rpact", {
    r <- run_gsd(endpoint = "surv", hr = 0.7, ctrlMedian = 12, accrualTime = 24, followUp = 12,
                 dropoutProb = 0.05, dropoutTime = 12, futility = "spending", futSf = "hsd", futGamma = -2)
    # rpact getSampleSizeSurvival: 352.3525 events, 550.0681 patients, look 1 at 16.62567
    expect_equal(cell(r$summary, "events", "value"), "353")
    expect_equal(cell(r$summary, "n", "value"), "551")
    tm <- r$timing$asDF
    expect_equal(tm$timeH1[1], 16.62567, tolerance = 1e-5)
    expect_equal(nrow(tm), 3)
})

test_that("every efficacy family and futility type runs, one- and two-sided", {
    for (sided in c("one", "two")) {
        for (ef in c("ldof", "ldpocock", "kd", "hsd", "exponential", "obf", "pocock", "wt", "hp", "pt", "triangular")) {
            r <- run_gsd(sided = sided, alpha = if (sided == "one") 0.025 else 0.05, efficacy = ef)
            expect_true(all(is.finite(r$bounds$asDF$effZ)), label = paste(sided, ef))
        }
        for (fu in c("spending", "fixedz", "cp")) {
            r <- run_gsd(sided = sided, alpha = if (sided == "one") 0.025 else 0.05, futility = fu, binding = TRUE)
            expect_true(all(is.finite(r$bounds$asDF$effZ)), label = paste(sided, fu))
        }
    }
    r <- run_gsd(efficacy = "user", effUser = "0.1, 0.4, 1")
    expect_true(all(is.finite(r$bounds$asDF$effZ)))
})

test_that("unequal timing given as cumulative sizes, and proportions", {
    r <- run_gsd(timing = "100, 250, 400", endpoint = "props", p0 = 0.3, p1 = 0.45)
    expect_equal(r$bounds$asDF$t, c(0.25, 0.625, 1))
})
