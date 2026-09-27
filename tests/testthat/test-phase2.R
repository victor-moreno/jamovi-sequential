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
