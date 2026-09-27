fx <- readRDS(test_path("fixtures", "phase2_oracle.rds"))

test_that("Simon optimal and minimax designs match clinfun::ph2simon", {
    for (cs in fx$simon) {
        s <- ph2_simon(cs$p0, cs$p1, cs$alpha, cs$beta)
        lab <- paste(cs$p0, cs$p1, cs$alpha, cs$beta)
        expect_equal(unlist(s$optimal[c("r1", "n1", "r", "n")], use.names = FALSE), cs$optimal[1:4], label = lab)
        expect_equal(s$optimal$en0, cs$optimal[5], tolerance = 1e-8, label = lab)
        expect_equal(unlist(s$minimax[c("r1", "n1", "r", "n")], use.names = FALSE), cs$minimax[1:4], label = lab)
        # the admissible set contains both, and each minimises q n + (1 - q) EN0 on its q-range
        a <- s$admissible
        expect_true(all(c(s$optimal$n, s$minimax$n) %in% a$n))
    }
})

test_that("two-stage inference: UMVUE and p-value match clinfun; CI reduces to Clopper-Pearson", {
    for (cs in fx$inference) {
        stage <- if (cs$x <= cs$r1) 1 else 2
        r <- ph2_twostage_inference(cs$x, stage, cs$r1, cs$n1, cs$n, cs$p0)
        expect_equal(r$umvue, cs$umvue, tolerance = 1e-10)
        expect_equal(r$p.value, cs$p.value, tolerance = 1e-10)
    }
    for (x in c(0, 3, 7, 20)) {
        r <- ph2_twostage_inference(x, 2, -1, 1, 20, .2)
        cp <- stats::binom.test(x, 20)$conf.int
        expect_equal(c(r$lower, r$upper), as.numeric(cp), tolerance = 1e-6)
    }
})

test_that("Bayesian predictive probability and boundaries match ph2bayes", {
    pp <- sapply(c(5, 10, 14, 17), function(y) pred_prob(y, 23, 40, .6, .4, .6, .9))
    expect_equal(pp, fx$bayes$pp, tolerance = 1e-10)
    b <- ph2_bayes(0.6, 0.8, 0.6, 0.4, 1:40, thetaT = 0.9, thetaL = 0.05, thetaU = 1)
    sb <- fx$bayes$stopbound
    sb <- sb[!is.na(sb$bound) & sb$n < 40, ]
    expect_equal(b$fut[sb$n], sb$bound)
})

test_that("Fleming boundaries reproduce the published examples", {
    f <- ph2_fleming(.2, .05, c(15, 25))
    expect_equal(f$fut, c(3, 8)); expect_equal(f$eff, c(7, 9))
    f <- ph2_fleming(.05, .05, c(25, 50))
    expect_equal(f$fut, c(0, 5)); expect_equal(f$eff, c(5, 6))
})

test_that("multistage operating characteristics match a two-stage closed form", {
    oc <- ms_oc(c(9, 17), c(0, 2), c(10, 3), .05)
    expect_equal(oc$reject, 0.04660496, tolerance = 1e-7)
    expect_equal(oc$en, 11.95800472, tolerance = 1e-7)
})

test_that("the phase2 analysis runs every design, with observed data", {
    skip_if_not_installed("jmvcore")
    r <- phase2(data = data.frame(), candidates = TRUE, analyse = TRUE, obsN = 43, obsX = 14)
    st <- r$stages$asDF
    expect_equal(st$n, c(13, 43))
    expect_equal(st$fut, c(3, 12))
    ob <- r$observed$asDF
    # Koyama & Chen p-value and Jung & Kim UMVUE for 14/43 after the 3/13 first stage
    ref <- ph2_twostage_inference(14, 2, 3, 13, 43, 0.2)
    expect_equal(ob$pvalue, ref$p.value, tolerance = 1e-10)
    expect_equal(ob$umvue, ref$umvue, tolerance = 1e-10)
    for (d in c("single", "fleming", "predictive", "posterior"))
        expect_true(nrow(phase2(data = data.frame(), design = d)$stages$asDF) >= 1, label = d)
})

test_that("invalid rates and observations give notices quickly", {
    skip_if_not_installed("jmvcore")
    t <- system.time(r <- phase2(data = data.frame(), p0 = 0.4, p1 = 0.3))[3]
    expect_lt(t, 5)
    expect_match(r$get("inputProblem")$content, "p1 must be larger")
    r <- phase2(data = data.frame(), analyse = TRUE, obsN = 10, obsX = 12)
    expect_match(r$get("observedProblem")$content, "cannot exceed")
})

test_that("the two-stage UMVUE is exactly unbiased (and the observed rate is not)", {
    r1 <- 3; n1 <- 13; n <- 43; n2 <- n - n1; p <- 0.3
    d1 <- dbinom(0:n1, n1, p); d2 <- dbinom(0:n2, n2, p)
    eu <- sum(d1[1:(r1 + 1)] * (0:r1) / n1); en <- eu
    for (a in (r1 + 1):n1) for (b in 0:n2) {
        w <- d1[a + 1] * d2[b + 1]
        eu <- eu + w * ph2_twostage_inference(a + b, 2, r1, n1, n, 0.2)$umvue
        en <- en + w * (a + b) / n
    }
    expect_equal(eu, p, tolerance = 1e-10)
    expect_lt(en, p - 0.03)
})

test_that("the outcome map describes two-stage designs and is refused for other designs", {
    r <- phase2(data = data.frame(), outcomePlot = TRUE)
    st <- r$outcomePlot$state
    expect_equal(unlist(st[c("n1", "n2", "fut1", "eff2")]), c(n1 = 13, n2 = 30, fut1 = 3, eff2 = 13))
    p <- plot_ph2_outcomes(st, c("blue", "red"), list(ggplot2::theme_minimal()),
                           list(stop0 = "s0", stop1 = "s1", no = "no", yes = "yes", x = "x", y = "y"))
    d <- p$data
    expect_equal(nrow(d), 14 * 31)
    expect_true(all(d$what[d$x1 <= 3] == "s0"))
    expect_true(all((d$what == "yes") == (d$x1 > 3 & d$x1 + d$x2 >= 13)))
    r <- phase2(data = data.frame(), design = "fleming", outcomePlot = TRUE)
    expect_false(r$outcomePlot$visible)
    # the notice sits where the map would be, not at the top
    nm <- r$itemNames
    expect_equal(match("outcomeMap", nm), match("outcomePlot", nm) - 1)
    r <- phase2(data = data.frame(), design = "single", outcomePlot = TRUE)
    expect_gt(match("outcomeMap", r$itemNames), match("ocPlot", r$itemNames))
})

test_that("the stopping-boundary plot shades the decisions per look; trajectories and projection toggle", {
    r <- phase2(data = data.frame(), analyse = TRUE, obsN = 11, obsX = 6)
    st <- r$boundPlot$state
    lab <- list(fut = "no", eff = "yes", cont = "cont", now = "Now: %d of %d (rate %.2f)", reach = "reach",
                proj = "fan %s%% %.2f-%.2f", x = "x", y = "y")
    base <- list(ggplot2::theme_minimal())
    p <- plot_ph2_bounds(st, NULL, base, lab)
    reg <- p$layers[[1]]$data
    # look 1 (13 patients): stop if <= 3 responses, else continue; look 2: promising if >= 13 of 43
    expect_setequal(unique(as.character(reg$what)), c("no", "cont", "yes"))
    expect_true(all(reg$y <= reg$x + 0.5 + 1e-9))
    polys <- function(p) Filter(function(l) inherits(l$geom, "GeomPolygon") && is.null(l$mapping$group), p$layers)
    expect_length(polys(p), 2)
    # the fan spans the slopes of the interval shown in the observed-data table
    fan <- polys(p)[[2]]$data
    ob <- r$observed$asDF
    expect_equal(fan$y[2:3], 6 + c(ob$lower, ob$upper) * (43 - 11), tolerance = 1e-10)
    st$show <- "traj"
    cone <- polys(plot_ph2_bounds(st, NULL, base, lab))
    expect_length(cone, 1)
    expect_equal(max(cone[[1]]$data$y), 6 + 43 - 11)
    st$show <- character(0)
    expect_length(polys(plot_ph2_bounds(st, NULL, base, lab)), 0)
})

test_that("the outcome map marks the current status and fades what can no longer be reached", {
    r <- phase2(data = data.frame(), analyse = TRUE, obsN = 11, obsX = 6, outcomePlot = TRUE)
    lab <- list(stop0 = "s0", stop1 = "s1", no = "no", yes = "yes", now = "now", x = "x", y = "y")
    p <- plot_ph2_outcomes(r$outcomePlot$state, NULL, list(ggplot2::theme_minimal()), lab)
    d <- p$data
    expect_equal(sort(unique(d$x1[d$reach])), 6:8)
    now <- p$layers[[2]]$data
    expect_equal(unlist(now), c(x1 = 6, x2 = 0))
    r <- phase2(data = data.frame(), analyse = TRUE, obsN = 25, obsX = 8, outcomePlot = TRUE)
    p <- plot_ph2_outcomes(r$outcomePlot$state, NULL, list(ggplot2::theme_minimal()), lab)
    now <- p$layers[[2]]$data
    # in stage 2 the split is unknown: every stage-1 count that passed (> 3) and fits 8 of 25
    expect_equal(now$x1, 4:8)
    expect_true(all(now$x1 + now$x2 == 8))
})

test_that("P(ends promising | data) is exact and handles planned looks and the end", {
    st <- list(nk = c(13, 43), fut = c(3, 12), eff = c(14, 13))
    # 6/11: stage 1 ends with 6-8 responses (never <= 3); promising if the total reaches 13
    manual <- function(p) sum(stats::dbinom(0:2, 2, p) * stats::pbinom(13 - (6 + 0:2) - 1, 30, p, lower.tail = FALSE))
    for (p in c(0.2, 0.4, 6 / 11)) expect_equal(ph2_cond_promising(st, 11, 6, p), manual(p), tolerance = 1e-12)
    expect_equal(ph2_cond_promising(st, 13, 3, 0.3), 0)
    expect_equal(ph2_cond_promising(st, 43, 13, 0.3), 1)
    expect_true(is.na(ph2_cond_promising(st, 50, 20, 0.3)))
    r <- phase2(data = data.frame(), analyse = TRUE, obsN = 11, obsX = 6)
    ob <- r$observed$asDF
    expect_equal(c(ob$cond0, ob$cond1), c(manual(0.2), manual(0.4)), tolerance = 1e-12)
})
