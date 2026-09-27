test_that("sample size, events, accrual, follow-up and look times match rpact", {
    cases <- readRDS(test_path("fixtures", "size_oracle.rds"))
    for (cs in cases) {
        d <- gs_design(3, sided = cs$sided, alpha = cs$alpha, eff = list(type = "spending", sf = "LDOF"),
                       fut = list(type = "spending", sf = "HSD", param = -2))
        lab <- paste(cs$source, "sided", cs$sided, "ratio", cs$ratio)
        expect_equal(size_means(d, 5, 10, cs$ratio)$Nmax, cs$means, tolerance = 1e-6, label = lab)
        expect_equal(size_props(d, .45, .3, cs$ratio)$Nmax, cs$props, tolerance = 1e-6, label = lab)
        D <- events_needed(d$eta, .7, cs$ratio)
        expect_equal(D, cs$events, tolerance = 1e-6, label = lab)
        lam0 <- log(2) / 12; gam <- dropout_hazard(.05, 12)
        pl <- surv_plan(D, lam0, .7, gam, cs$ratio, solve = "n", A = 24, F = 12)
        expect_equal(pl$N, cs$nSurv, tolerance = 1e-6, label = lab)
        tm <- surv_look_times(d$t * D, pl$N, 24, 12, lam0, .7, gam, cs$ratio)
        expect_equal(tm, cs$times, tolerance = 1e-6, label = lab)
        expect_equal(surv_plan(D, lam0, .7, gam, cs$ratio, solve = "accrual", rate = 20, F = 12)$A,
                     cs$accrual, tolerance = 1e-6, label = lab)
        expect_equal(surv_plan(D, lam0, .7, gam, cs$ratio, solve = "followup", A = 24, N = 600)$F,
                     cs$followup, tolerance = 1e-6, label = lab)
    }
})
