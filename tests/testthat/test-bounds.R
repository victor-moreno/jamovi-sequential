test_that("engine reproduces the classic error inflation of repeated tests", {
    t <- (1:5) / 5
    p <- gs_probs(t, up = rep(qnorm(.975), 5), lo = rep(-qnorm(.975), 5))
    expect_equal(sum(p$up + p$lo), 0.1417, tolerance = 1e-4)
})

test_that("classical constants match Jennison & Turnbull (2000) ch. 2", {
    poc <- gs_design(5, sided = 2, alpha = .05, eff = list(type = "WT", delta = .5))
    obf <- gs_design(5, sided = 2, alpha = .05, eff = list(type = "WT", delta = 0))
    expect_equal(poc$C, 2.413, tolerance = 1e-3)
    expect_equal(obf$C, 2.040, tolerance = 1e-3)
    expect_equal(poc$inflation, 1.207, tolerance = 1e-3)
    expect_equal(obf$inflation, 1.026, tolerance = 1e-3)
})

test_that("bounds, futility and inflation match rpact / gsDesign", {
    cases <- readRDS(test_path("fixtures", "bounds_oracle.rds"))
    for (cs in cases) {
        d <- do.call(gs_design, cs$args)
        K <- length(d$up)
        lab <- paste(cs$source, deparse(cs$args, width.cutoff = 500L)[1])
        expect_equal(d$up, cs$up, tolerance = 1e-5, label = lab)
        expect_equal(d$inflation, cs$IF, tolerance = 1e-5, label = lab)
        if (!all(is.na(cs$fut))) {
            f <- if (d$sided == 1) d$lo[-K] else d$wedge[-K]
            f[!is.finite(f) | (d$sided == 2 & f <= 1e-5)] <- NA
            expect_equal(f, cs$fut, tolerance = 1e-5, label = lab)
        }
    }
})

test_that("non-oracle futility designs control alpha and power by construction", {
    d <- gs_design(3, sided = 2, alpha = .05, eff = list(type = "spending", sf = "LDOF"),
                   fut = list(type = "fixedz", z = c(.3, .8)), binding = TRUE)
    expect_equal(d$alphaBinding, .05, tolerance = 1e-7)
    expect_equal(d$power, .9, tolerance = 1e-7)
    d <- gs_design(4, alpha = .025, eff = list(type = "spending", sf = "LDOF"),
                   fut = list(type = "cp", cp = .2, under = "trend"), binding = TRUE)
    expect_equal(d$alphaBinding, .025, tolerance = 1e-7)
    # the binding bound is self-consistent: CP under the trend equals the threshold at it
    t <- d$t[2]; z <- d$lo[2]
    cp <- 1 - pnorm((d$up[4] - z / sqrt(t)) / sqrt(1 - t))
    expect_equal(cp, .2, tolerance = 1e-6)
})

test_that("triangular tests follow Whitehead's formulas and attain alpha and power", {
    d <- tri_design(5, sided = 1, alpha = .025, beta = .1)
    k <- 1 + qnorm(.9) / qnorm(.975)
    expect_equal(unname(d$triangle["a"]), k * log(1 / .05))
    expect_equal(unname(d$triangle["c"]), 1 / (2 * k))
    expect_lt(abs(d$alphaBinding - .025), .003)
    expect_lt(abs(d$power - .9), .01)
    # the corrected boundaries close exactly at the last look, so no look is unreachable
    expect_true(all(d$lo[-5] < d$up[-5]))
    expect_equal(d$lo[5], d$up[5])
    d40 <- tri_design(40, sided = 1, alpha = .025, beta = .1)
    expect_lt(abs(d40$alphaBinding - .025), 5e-4)
    dd <- tri_design(5, sided = 2, alpha = .05, beta = .1)
    expect_lt(abs(dd$alphaBinding - .05), .005)
})

test_that("engine exit probabilities agree with simulation for a two-sided inner-wedge design", {
    d <- gs_design(3, sided = 2, alpha = .05, eff = list(type = "spending", sf = "LDPocock"),
                   fut = list(type = "fixedz", z = c(.4, .9)), binding = TRUE)
    set.seed(11)
    n <- 1e5; t <- d$t; th <- d$eta * .6
    s <- matrix(rnorm(n * 3), n) * rep(sqrt(diff(c(0, t))), each = n) + rep(th * diff(c(0, t)), each = n)
    z <- t(apply(s, 1, cumsum)) / rep(sqrt(t), each = n)
    stopAt <- rep(NA, n); side <- rep(NA, n)
    for (k in 1:3) {
        open <- is.na(stopAt)
        u <- open & z[, k] >= d$up[k]; l <- open & z[, k] <= d$lo[k]; w <- open & abs(z[, k]) < d$wedge[k]
        stopAt[u | l | w] <- k; side[u] <- "up"; side[l] <- "lo"; side[w] <- "wedge"
    }
    p <- gs_probs(t, d$up, d$lo, d$wedge, theta = th)
    sim <- c(mean(side == "up"), mean(side == "wedge"))
    ex <- c(sum(p$up), sum(p$wedge))
    expect_true(all(abs(sim - ex) < 3 * sqrt(ex * (1 - ex) / n) + 1e-4))
})
