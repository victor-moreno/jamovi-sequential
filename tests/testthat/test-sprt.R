fx <- readRDS(test_path("fixtures", "sprt_oracle.rds"))

test_that("LLR paths match the SPRT package and sprtt", {
    b <- fx$sprt$binom
    o <- sprt_analyse("binom", b$x, NULL, list(p0 = .2, p1 = .4), .05, .2)
    expect_equal(o$stat[seq_along(b$llr)], b$llr, tolerance = 1e-10)
    expect_equal(if (is.na(b$n)) "continue" else o$n, if (is.na(b$n)) o$decision else b$n)
    n1 <- fx$sprt$norm1
    o <- sprt_analyse("norm1", n1$x, NULL, list(mu0 = 0, mu1 = .5, sigma = 1), .05, .2)
    expect_equal(o$stat[seq_along(n1$llr)], n1$llr, tolerance = 1e-10)
    expect_equal(o$n, n1$n)
    t1 <- fx$sprt$t1
    ours <- sapply(c("one", "two"), function(tl) utils::tail(sprt_llr_path("t1", t1$x, NULL, list(mu0 = 0, d1 = .5), tl), 1))
    expect_equal(unname(ours), unname(t1$llr), tolerance = 1e-8)
    t2 <- fx$sprt$t2
    expect_equal(utils::tail(sprt_llr_path("t2", t2$x, t2$g, list(d1 = .5), "one"), 1), t2$llr, tolerance = 1e-8)
})

test_that("exact binomial OC matches gsDesign::binomialSPRT", {
    par <- list(p0 = .2, p1 = .4)
    tr <- sapply(c(.2, .4), function(p) sprt_binom_exact(p, par, .05, .2, nMax = 60, truncRule = "upper")[1:2])
    expect_equal(c(tr[1, ], tr[2, ]), unname(fx$sprt$binomExact$trunc60), tolerance = 1e-7)
    op <- sapply(c(.2, .4), function(p) sprt_binom_exact(p, par, .05, .2)[1:2])
    expect_equal(c(op[1, ], op[2, ]), unname(fx$sprt$binomExact$open), tolerance = 1e-6)
})

test_that("Wald approximations: alpha and power at the hypotheses, the zero-drift point for the binomial", {
    par <- list(p0 = .2, p1 = .4)
    expect_equal(sprt_wald_oc("binom", .2, par, .05, .2)[["pH1"]], .05, tolerance = 1e-8)
    expect_equal(sprt_wald_oc("binom", .4, par, .05, .2)[["pH1"]], .8, tolerance = 1e-8)
    pn <- list(mu0 = 0, mu1 = .5, sigma = 1)
    expect_equal(sprt_wald_oc("norm1", 0, pn, .05, .2)[["pH1"]], .05, tolerance = 1e-8)
    # the OC is continuous through the zero-drift rate (the slope of the boundary lines)
    s <- -log(.6 / .8) / (log(2) - log(.6 / .8))
    expect_equal(sprt_wald_oc("binom", s, par, .05, .2)[["pH1"]], sprt_wald_oc("binom", s + 1e-6, par, .05, .2)[["pH1"]],
                 tolerance = 1e-4)
})

test_that("exact normal OC (engine) agrees with simulation", {
    par <- list(mu0 = 0, mu1 = .5, sigma = 1)
    for (th in c(0, .5)) {
        e <- sprt_norm_exact(th, "norm1", par, .05, .2)
        s <- sprt_simulate(th, "norm1", par, .05, .2, nsim = 4000, seed = 2)
        expect_lt(abs(e[["pH1"]] - s[["pH1"]]), 0.015)
        expect_lt(abs(e[["asn"]] / s[["asn"]] - 1), 0.04)
    }
})

test_that("the sprt analysis runs in design and data modes", {
    r <- sprt(data = data.frame())
    expect_equal(r$oc$asDF$pH1[1], unname(fx$sprt$binomExact$open[1]), tolerance = 1e-6)
    set.seed(3)
    dat <- data.frame(y = rbinom(60, 1, 0.4), x = rnorm(60, 0.4), g = factor(rep(c("ctl", "trt"), 30)))
    r <- sprt(data = dat, dep = "x", group = "g", model = "t2", d1 = 0.8, nsim = 300)
    expect_equal(nrow(r$test$asDF), 1)
    r <- sprt(data = dat, dep = "y", boundary = "triangular", truncate = TRUE, nMax = 40)
    expect_equal(r$test$asDF$n, 30)
    r <- sprt(data = data.frame(), model = "t1", boundary = "triangular")
    expect_false(r$summary$visible)
})

test_that("MaxSPRT critical values and performance match the Sequential package", {
    for (cs in fx$maxsprt) {
        looks <- if (cs$monitoring == "group") cs$looks else NULL
        size <- if (is.null(cs$size)) max(cs$looks) else cs$size
        d <- maxsprt_design(cs$model, cs$monitoring, 0.05, size, cs$M, looks, cs$z, c(1.5, 2))
        lab <- paste(cs$model, cs$monitoring, size, cs$M, cs$z)
        # Sequential reports the lower end of the interval of equivalent critical values, rounded up
        expect_lt(abs(d$cv - cs$cv), 1e-5, label = lab)
        if (!is.null(cs$alpha)) expect_equal(d$alpha, cs$alpha, tolerance = 1e-6, label = lab)
        expect_equal(d$performance$power, unname(cs$perf[1, ]), tolerance = 1e-6, label = lab)
        expect_equal(d$performance$timeToSignal, unname(cs$perf[2, ]), tolerance = 1e-6, label = lab)
        expect_equal(d$performance$length, unname(cs$perf[3, ]), tolerance = 1e-6, label = lab)
    }
})

test_that("the maxsprt analysis runs and signals", {
    dat <- data.frame(obs = c(1, 3, 5, 9, 12, 16), exp = c(1.5, 3, 4.2, 5.1, 6, 7))
    r <- maxsprt(data = dat, cases = "obs", expected = "exp")
    m <- r$monitor$asDF
    expect_equal(m$decision[6], "Signal: reject H0")
    expect_equal(r$summary$asDF$value[6], "3.628123")
    r <- maxsprt(data = dat, model = "binomial", cases = "obs", expected = "exp")
    expect_false(r$monitor$visible)
})
