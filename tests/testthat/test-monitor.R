fx <- readRDS(test_path("fixtures", "monitor_oracle.rds"))

test_that("monitoring matches rpact: repeated CIs, conditional power, stagewise inference", {
    for (cs in fx) {
        d <- do.call(gs_design, cs$design)
        m <- gs_monitor(d, cs$looks, final = cs$final)
        lab <- paste(cs$source, cs$label)
        # repeated CIs use the two-sided test without futility (Jennison & Turnbull 1989);
        # rpact's binding-futility intervals differ by < 1e-4
        tol <- if (cs$binding) 1e-3 else 1e-6
        expect_equal(m$looks$rciLo, cs$rciLo, tolerance = tol, label = lab)
        expect_equal(m$looks$rciHi, cs$rciHi, tolerance = tol, label = lab)
        if (!is.na(cs$cp) && is.na(m$stopAt)) expect_equal(m$looks$cpH1[nrow(cs$looks)], cs$cp, tolerance = 1e-6, label = lab)
        if (!is.null(m$inference)) {
            expect_equal(m$inference$p.value, cs$p, tolerance = 1e-5, label = lab)
            expect_equal(m$inference$mue, cs$mue, tolerance = 1e-5, label = lab)
            expect_equal(c(m$inference$lower, m$inference$upper), c(cs$lower, cs$upper), tolerance = 1e-5, label = lab)
        }
    }
})

test_that("two-sided stagewise inference is symmetric in the direction of the effect", {
    d <- gs_design(3, sided = 2, alpha = .05, beta = .2, eff = list(type = "spending", sf = "LDOF"))
    lk <- data.frame(t = c(1/3, 2/3), z = c(0.3, 3.2), est = c(0.1, 0.8), se = c(0.8, 0.25) / c(1, 1))
    a <- gs_monitor(d, lk)$inference
    lk$z <- -lk$z; lk$est <- -lk$est
    b <- gs_monitor(d, lk)$inference
    expect_equal(a$p.value, b$p.value, tolerance = 1e-10)
    expect_equal(a$mue, -b$mue, tolerance = 1e-8)
    expect_equal(c(a$lower, a$upper), -c(b$upper, b$lower), tolerance = 1e-8)
})

test_that("boundaries at the planned information equal the design boundaries", {
    for (e in list(list(type = "spending", sf = "LDOF"), list(type = "spending", sf = "HSD", param = -2),
                   list(type = "WT", delta = 0))) {
        d <- gs_design(4, alpha = .025, beta = .1, eff = e, fut = list(type = "spending", sf = "LDOF"))
        b <- gs_bounds_at(d, d$t)
        expect_equal(b$up, d$up, tolerance = 1e-7)
        expect_equal(b$lo, d$lo, tolerance = 1e-7)
    }
})
