fx <- readRDS(test_path("fixtures", "bayes_oracle.rds"))

test_that("the JZS Bayes factor matches BayesFactor and adaptive integration", {
    b <- fx$bf
    two <- !is.na(b$n2)
    nEff <- ifelse(two, b$n1 * b$n2 / (b$n1 + b$n2), b$n1)
    nu <- ifelse(two, b$n1 + b$n2 - 2, b$n1 - 1)
    ours <- vapply(seq_len(nrow(b)), function(i) jzs_logbf(b$t[i], nEff[i], nu[i], b$r[i]), 0)
    # BayesFactor's own approximation is accurate to about 2e-5 on the log scale
    expect_lt(max(abs(ours - b$logbf)), 5e-5)
    ex <- function(t, nEff, nu, r) log(stats::integrate(function(g) exp(-0.5 * log1p(nEff * g) + ((nu + 1) / 2) *
        (log1p(t^2 / nu) - log1p(t^2 / (nu * (1 + nEff * g))))) * r / sqrt(2 * pi) * g^-1.5 * exp(-r^2 / (2 * g)),
        0, Inf, rel.tol = 1e-12, subdivisions = 1000)$value)
    for (i in c(1, 17, 33, 44)) expect_equal(ours[i], ex(b$t[i], nEff[i], nu[i], b$r[i]), tolerance = 1e-9)
})

test_that("two-arm posterior probability equals the exact formula; predictive probability equals enumeration", {
    exact <- function(aC, bC, aT, bT) { i <- 0:(aT - 1)
        1 - sum(exp(lbeta(aC + i, bC + bT) - log(bT + i) - lbeta(1 + i, bT) - lbeta(aC, bC))) }
    for (cs in list(c(3, 19, 12, 10), c(1, 61, 12, 50), c(40, 62, 55, 47), c(2, 2, 2, 2))) {
        # exact() gives P(pC > pT) with the arms swapped, so P(pT - pC > 0) = 1 - P(pC > pT)
        expect_equal(bb2_prob(cs[1], cs[2], cs[3], cs[4], 0), exact(cs[3], cs[4], cs[1], cs[2]), tolerance = 1e-12)
    }
    d <- bb2_design(c(20, 40, 60), 1, 1, 0, 0.975, 0.1)
    s <- d$eff[[3]]
    pc <- bb_pmf(40, 1 + 8, 1 + 12); pt <- bb_pmf(35, 1 + 12, 1 + 13)
    brute <- sum(outer(0:40, 0:35, function(yc, yt) (12 + yt >= s[8 + yc + 1])) * outer(pc, pt))
    expect_equal(bb2_pp(8, 20, 12, 25, 60, 1, 1, s), brute, tolerance = 1e-12)
})

test_that("exact two-arm operating characteristics agree with simulation", {
    d <- bb2_design(c(20, 40, 60), 1, 1, 0, 0.975, 0.1)
    oc <- bb2_oc(d, 0.3, 0.5)
    set.seed(11); nsim <- 20000; rej <- 0
    for (i in seq_len(nsim)) {
        xc <- cumsum(stats::rbinom(3, 20, 0.3)); xt <- cumsum(stats::rbinom(3, 20, 0.5))
        for (k in 1:3) { e <- xt[k] >= d$eff[[k]][xc[k] + 1]; f <- xt[k] <= d$fut[[k]][xc[k] + 1]
            if (e || f) { rej <- rej + e; break } }
    }
    expect_lt(abs(oc$reject - rej / nsim), 0.012)
})

test_that("normal approximation: flat prior reproduces repeated significance testing; thresholds are consistent", {
    d <- bayes_normal_design(c(60, 120, 180), 2, bayes_prior_normal("flat", 0.5), 0, 0.975)
    ref <- sum(gs_probs(c(60, 120, 180) / 4, rep(stats::qnorm(0.975), 3))$up)
    expect_equal(bayes_normal_oc(d, 0)$reject, ref, tolerance = 1e-10)
    pr <- bayes_prior_normal("sceptical", 0.5)
    f <- bayes_fut_est(15, 45, pr, 0.1, 0.95, 0.2)
    expect_equal(bayes_pp_normal(f, 15, 45, pr, 0.1, 0.95), 0.2, tolerance = 1e-12)
    e <- bayes_eff_est(15, pr, 0.1, 0.95)
    post <- bayes_post_normal(e, 15, pr)
    expect_equal(stats::pnorm((post$mean - 0.1) / post$sd), 0.95, tolerance = 1e-12)
})

test_that("the bayesmon analysis runs for the three models", {
    dat <- read.csv(test_path("..", "..", "data", "gsd_interim_means.csv"))
    r <- bayesmon(data = dat, lookN = "n", estimate = "estimate", se = "se")
    m <- r$monitor$asDF
    expect_equal(m$decision, c("Continue", "Stop: efficacy"))
    r <- bayesmon(data = data.frame(), model = "binary")
    expect_equal(r$oc$asDF$reject[1], bb2_oc(bb2_design(c(20, 40, 60), 1, 1, 0, 0.975, 0.1), 0.3, 0.3)$reject, tolerance = 1e-12)
    sq <- read.csv(test_path("..", "..", "data", "sprt_sequential.csv"))
    sq$arm <- factor(sq$arm)
    r <- bayesmon(data = sq, model = "bf", dep = "score", group = "arm", nsim = 200)
    expect_true(any(grepl("Evidence for H1", r$bfTest$asDF$decision)))
})
