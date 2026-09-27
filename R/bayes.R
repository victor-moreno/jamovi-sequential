# Bayesian sequential monitoring.
#
# Normal approximation (estimate and SE per look, normal prior): posterior, posterior
# probability of an effect beyond delta, and predictive probability of success at the final
# analysis (Spiegelhalter, Freedman & Parmar 1994). Every rule is a threshold on the estimate at
# each look, so the frequentist operating characteristics are exact with the gsd engine.
# Two arms, binary outcome: independent Beta posteriors, P(pT - pC > delta) by quadrature,
# predictive probability by beta-binomial sums, exact operating characteristics by dynamic
# programming over the response counts of both arms.
# JZS Bayes factor t-test for raw data (Rouder et al. 2009); Bayes factor design analysis
# (Schönbrodt & Wagenmakers 2018) by seeded simulation.
#
# Effects are oriented so that positive values favour the treatment (dir = -1 flips them).

# ---- normal approximation ------------------------------------------------------------

# Prior (mean, precision) on the oriented scale. Sceptical: centred on 0 with P(theta > theta1)
# = 5 %; enthusiastic: centred on theta1 with P(theta < 0) = 5 % (Spiegelhalter et al. 1994).
bayes_prior_normal <- function(type, theta1, mean = 0, sd = 1) {
    prec1 <- (stats::qnorm(0.95) / theta1)^2
    switch(type,
        flat = c(mean = 0, prec = 0),
        sceptical = c(mean = 0, prec = prec1),
        enthusiastic = c(mean = theta1, prec = prec1),
        custom = c(mean = mean, prec = 1 / sd^2))
}

bayes_post_normal <- function(est, I, prior) {
    prec <- prior[["prec"]] + I
    list(mean = (prior[["prec"]] * prior[["mean"]] + I * est) / prec, sd = 1 / sqrt(prec))
}

# Smallest estimate at information I with P(theta > delta | data) >= gamma.
bayes_eff_est <- function(I, prior, delta, gamma) {
    P <- prior[["prec"]] + I
    (P * (delta + stats::qnorm(gamma) / sqrt(P)) - prior[["prec"]] * prior[["mean"]]) / I
}

# Predictive probability that the final analysis (information IK) meets the efficacy rule:
# the final estimate is normal given the current posterior.
bayes_pp_normal <- function(est, I, IK, prior, delta, gamma) {
    eK <- bayes_eff_est(IK, prior, delta, gamma)
    P <- prior[["prec"]] + I
    R <- IK - I
    m <- (prior[["prec"]] * prior[["mean"]] + I * est) / P
    mu <- (I * est + R * m) / IK
    s <- sqrt(R + R^2 / P) / IK
    ifelse(R <= 1e-12 * IK, as.numeric(est >= eK), stats::pnorm((mu - eK) / pmax(s, 1e-300)))
}

# Largest estimate with predictive probability below gammaF (linear in the estimate).
bayes_fut_est <- function(I, IK, prior, delta, gamma, gammaF) {
    eK <- bayes_eff_est(IK, prior, delta, gamma)
    P <- prior[["prec"]] + I
    R <- IK - I
    s <- sqrt(R + R^2 / P) / IK
    slope <- I * (1 + R / P) / IK
    off <- R * prior[["prec"]] * prior[["mean"]] / (P * IK)
    (eK + stats::qnorm(gammaF) * s - off) / slope
}

# Design: cumulative sizes nk, standard error of one unit sigmaU (SE = sigmaU / sqrt(n)).
# Stop for efficacy if P(theta > delta) >= gammaE (at interim looks if effInterim), for futility
# if the predictive probability of final success is < gammaF (NA: no futility rule).
bayes_normal_design <- function(nk, sigmaU, prior, delta, gammaE, gammaF = NA, effInterim = TRUE) {
    I <- nk / sigmaU^2
    K <- length(I)
    e <- bayes_eff_est(I, prior, delta, gammaE)
    if (!effInterim && K > 1) e[-K] <- Inf
    f <- rep(-Inf, K)
    if (!is.na(gammaF) && K > 1) f[-K] <- bayes_fut_est(I[-K], I[K], prior, delta, gammaE, gammaF)
    f[K] <- e[K]
    f <- pmin(f, e)
    list(nk = nk, I = I, eff = e, fut = f, up = pmin(e * sqrt(I), 40), lo = pmax(f * sqrt(I), -40),
         prior = prior, delta = delta, gammaE = gammaE, gammaF = gammaF, sigmaU = sigmaU)
}

# Frequentist operating characteristics at the (oriented) true effect theta.
bayes_normal_oc <- function(des, theta, futility = TRUE) {
    K <- length(des$I)
    lo <- if (futility) des$lo else c(rep(-40, K - 1), des$lo[K])
    p <- gs_probs(des$I, des$up, lo, theta = theta)
    list(pEff = p$up, pFut = p$lo, reject = sum(p$up), futStop = sum(p$lo[-K]),
         en = sum(des$nk * (p$up + p$lo)))
}

# ---- two arms, binary outcome ---------------------------------------------------------

# Tanh-sinh (double-exponential) nodes and weights on (0, 1): exponential convergence even
# with endpoint singularities (5e-14 against the exact formula for integer parameters).
tanh_sinh <- function(h = 1 / 10, tmax = 3.2) {
    t <- seq(-tmax, tmax, by = h)
    s <- pi / 2 * sinh(t)
    list(x = (1 + tanh(s)) / 2, w = h * pi / 4 * cosh(t) / cosh(s)^2)
}

# P(pT - pC > delta) with pC ~ Beta(aC, bC), pT ~ Beta(aT, bT) (vectorised over aT, bT).
# Integrated over the quantile u of pC, so the integrand is bounded even when the Beta
# density is not (parameters below 1).
bb2_prob <- function(aC, bC, aT, bT, delta, gl = tanh_sinh()) {
    # the integrand is 1 where pC + delta < 0 and 0 where pC + delta > 1: integrate only between
    u0 <- if (delta < 0) stats::pbeta(-delta, aC, bC) else 0
    u1 <- if (delta > 0) stats::pbeta(1 - delta, aC, bC) else 1
    q <- stats::qbeta(u0 + (u1 - u0) * gl$x, aC, bC) + delta
    vapply(seq_along(aT), function(i) u0 + (u1 - u0) * sum(gl$w * stats::pbeta(q, aT[i], bT[i], lower.tail = FALSE)), 0)
}

# Beta-binomial probabilities of 0..m responses among m, Beta(a, b) for the rate.
bb_pmf <- function(m, a, b) {
    y <- 0:m
    exp(lchoose(m, y) + lbeta(y + a, m - y + b) - lbeta(a, b))
}

# For each xC in 0..n, the smallest xT meeting P(pT - pC > delta | data) >= gamma (n + 1: none).
bb2_eff_thresholds <- function(n, a, b, delta, gamma, gl = tanh_sinh()) {
    vapply(0:n, function(xC) {
        f <- function(xT) bb2_prob(a + xC, b + n - xC, a + xT, b + n - xT, delta, gl) >= gamma
        if (!f(n)) return(n + 1)
        lo <- -1; hi <- n                          # f(lo) FALSE by convention, f(hi) TRUE
        while (hi - lo > 1) { mid <- (lo + hi) %/% 2; if (f(mid)) hi <- mid else lo <- mid }
        hi
    }, 0)
}

# Predictive probability of final success given xC of nC and xT of nT; s = final thresholds
# (for 0..N responses in the control arm), N = final size per arm.
bb2_pp <- function(xC, nC, xT, nT, N, a, b, s) {
    mC <- max(0, N - nC); mT <- max(0, N - nT)
    pc <- bb_pmf(mC, a + xC, b + nC - xC)
    pt <- bb_pmf(mT, a + xT, b + nT - xT)
    tail <- c(rev(cumsum(rev(pt))), 0)             # P(yT >= j), j = 0..mT + 1
    need <- s[pmin(xC + 0:mC, N) + 1] - xT         # treatment responses still needed
    sum(pc * tail[pmin(pmax(need, 0), mT + 1) + 1])
}

# For each xC, the largest xT with predictive probability < gammaF (-1: none).
bb2_fut_thresholds <- function(n, N, a, b, s, gammaF) {
    vapply(0:n, function(xC) {
        f <- function(xT) bb2_pp(xC, n, xT, n, N, a, b, s) < gammaF
        if (!f(0)) return(-1)
        lo <- 0; hi <- n + 1                       # f(lo) TRUE, f(hi) FALSE by convention
        while (hi - lo > 1) { mid <- (lo + hi) %/% 2; if (f(mid)) lo <- mid else hi <- mid }
        lo
    }, 0)
}

# Design with nk per arm at each look, common Beta(a, b) prior.
bb2_design <- function(nk, a, b, delta, gammaE, gammaF = NA, effInterim = TRUE) {
    gl <- tanh_sinh()
    K <- length(nk); N <- nk[K]
    s <- bb2_eff_thresholds(N, a, b, delta, gammaE, gl)
    eff <- lapply(seq_len(K), function(k)
        if (k == K) s else if (effInterim) bb2_eff_thresholds(nk[k], a, b, delta, gammaE, gl) else rep(nk[k] + 1, nk[k] + 1))
    fut <- lapply(seq_len(K), function(k)
        if (k == K) s - 1 else if (is.na(gammaF)) rep(-1, nk[k] + 1) else bb2_fut_thresholds(nk[k], N, a, b, s, gammaF))
    list(nk = nk, eff = eff, fut = fut, a = a, b = b, delta = delta, gammaE = gammaE, gammaF = gammaF)
}

# Exact operating characteristics at true rates pC, pT (dynamic programming).
bb2_oc <- function(des, pC, pT, futility = TRUE) {
    K <- length(des$nk)
    dist <- matrix(1, 1, 1)                        # rows: xC = 0.., columns: xT = 0..
    prev <- 0
    pEff <- pFut <- numeric(K)
    trans <- function(n0, n1, p) {
        d <- n1 - n0
        outer(0:n1, 0:n0, function(i, j) stats::dbinom(i - j, d, p))
    }
    for (k in seq_len(K)) {
        n <- des$nk[k]
        dist <- trans(prev, n, pC) %*% dist %*% t(trans(prev, n, pT))
        xT <- matrix(0:n, n + 1, n + 1, byrow = TRUE)
        e <- xT >= des$eff[[k]]                    # the thresholds recycle down the rows (xC)
        f <- if (futility || k == K) xT <= des$fut[[k]] & !e else matrix(FALSE, n + 1, n + 1)
        pEff[k] <- sum(dist[e]); pFut[k] <- sum(dist[f])
        dist[e | f] <- 0
        prev <- n
    }
    list(pEff = pEff, pFut = pFut, reject = sum(pEff), futStop = sum(pFut[-K]),
         en = sum(des$nk * (pEff + pFut)))
}

# ---- JZS Bayes factor for t-tests ----------------------------------------------------------

# log BF10 of the two-sided JZS test (Cauchy prior with scale r on the standardised effect,
# Jeffreys prior on the variance): an integral over g ~ InvGamma(1/2, r^2/2) of the t
# likelihood ratio, on a fixed grid in log g (vectorised over t, nEff, nu).
jzs_logbf <- function(t, nEff, nu, r = sqrt(2) / 2) {
    # trapezoid in log g: the integrand is analytic, step 0.2 is exact to 1e-13 (checked
    # against adaptive integration; BayesFactor::ttest.tstat agrees to 2e-5)
    h <- 0.2
    u <- seq(2 * log(r) - 10, 2 * log(r) + 45, by = h)
    g <- exp(u)
    lp <- log(r / sqrt(2 * pi)) - u / 2 - r^2 / (2 * g)        # log prior density times g
    A <- 1 + outer(nEff, g)                                     # rows: observations
    q <- t^2 / nu
    L <- -0.5 * log(A) + ((nu + 1) / 2) * (log1p(q) - log1p(q / A))
    M <- sweep(L, 2, lp, `+`)
    mx <- apply(M, 1, max)
    mx + log(rowSums(exp(M - mx)) * h)
}

# Running t statistics: one sample against mu0 (g NULL) or two samples (g = 1 / 2, 2 minus 1).
t_path <- function(x, g = NULL, mu0 = 0) {
    if (is.null(g)) {
        n <- seq_along(x)
        m <- cumsum(x) / n
        v <- (cumsum(x^2) - n * m^2) / pmax(n - 1, 1)
        ok <- n >= 2 & v > 0
        return(list(t = ifelse(ok, (m - mu0) / sqrt(v / n), 0), nEff = n, nu = n - 1, ok = ok,
                    n1 = n, n2 = NULL, est = m, d = ifelse(ok, (m - mu0) / sqrt(v), NA)))
    }
    st <- two_sample_stats(x, g)
    nu <- st$n1 + st$n2 - 2
    ok <- st$n1 >= 1 & st$n2 >= 1 & nu >= 1 & !is.na(st$sp) & st$sp > 0
    se <- st$sp * sqrt(1 / st$n1 + 1 / st$n2)
    list(t = ifelse(ok, st$diff / se, 0), nEff = st$n1 * st$n2 / (st$n1 + st$n2), nu = nu, ok = ok,
         n1 = st$n1, n2 = st$n2, est = st$diff, d = ifelse(ok, st$diff / st$sp, NA))
}

# Sequential Bayes factor test on data in arrival order: from nMin observations (per group
# for two samples), stop when BF10 >= k (H1) or <= 1 / k (H0); at nMax (total, NA: none) the
# result is inconclusive.
bf_analyse <- function(x, g = NULL, r = sqrt(2) / 2, k = 10, nMin = 10, nMax = NA, mu0 = 0) {
    p <- t_path(x, g, mu0)
    lbf <- rep(0, length(x))
    if (any(p$ok)) lbf[p$ok] <- jzs_logbf(p$t[p$ok], p$nEff[p$ok], p$nu[p$ok], r)
    elig <- p$ok & (if (is.null(g)) p$n1 >= nMin else p$n1 >= nMin & p$n2 >= nMin)
    N <- length(x)
    hit <- which(elig & abs(lbf) >= log(k))[1]
    if (!is.na(nMax) && !is.na(hit) && hit > nMax) hit <- NA
    end <- if (!is.na(hit)) hit else if (!is.na(nMax) && N >= nMax) nMax else NA
    dec <- if (!is.na(hit)) (if (lbf[hit] > 0) "H1" else "H0") else if (!is.na(end)) "inconclusive" else "continue"
    list(logbf = lbf, eligible = elig, n = if (is.na(end)) N else end, decision = dec, N = N,
         t = p$t, est = p$est, d = p$d, n1 = p$n1, n2 = p$n2)
}

# Bayes factor design analysis: seeded simulation of the sequential Bayes factor test under a
# standardised effect delta (one sample, or two samples adding one observation per group).
bf_simulate <- function(delta, two = FALSE, r = sqrt(2) / 2, k = 10, nMin = 10, nMax = 100,
                        nsim = 1000, seed = 1) {
    set.seed(seed)
    dec <- character(nsim); nAt <- numeric(nsim)
    m <- nMin:nMax                                          # per group (two samples) or total
    for (i in seq_len(nsim)) {
        if (two) {
            xc <- stats::rnorm(nMax); xt <- stats::rnorm(nMax, delta)
            j <- seq_len(nMax)
            m1 <- cumsum(xc) / j; m2 <- cumsum(xt) / j
            ss <- cumsum(xc^2) - j * m1^2 + cumsum(xt^2) - j * m2^2
            sp <- sqrt(pmax(ss, 0) / pmax(2 * j - 2, 1))
            tt <- ((m2 - m1) / (sp * sqrt(2 / j)))[m]; nEff <- m / 2; nu <- 2 * m - 2
        } else {
            x <- stats::rnorm(nMax, delta)
            j <- seq_len(nMax)
            mm <- cumsum(x) / j; v <- (cumsum(x^2) - j * mm^2) / pmax(j - 1, 1)
            tt <- (mm / sqrt(v / j))[m]; nEff <- m; nu <- m - 1
        }
        # evaluate a growing prefix: most paths stop early
        len <- min(length(m), 16)
        repeat {
            lb <- jzs_logbf(tt[seq_len(len)], nEff[seq_len(len)], nu[seq_len(len)], r)
            hit <- which(abs(lb) >= log(k))[1]
            if (!is.na(hit) || len == length(m)) break
            len <- min(length(m), 2 * len)
        }
        dec[i] <- if (is.na(hit)) "inconclusive" else if (lb[hit] > 0) "H1" else "H0"
        nAt[i] <- (if (is.na(hit)) nMax else m[hit]) * (if (two) 2 else 1)
    }
    c(pH1 = mean(dec == "H1"), pH0 = mean(dec == "H0"), pInc = mean(dec == "inconclusive"),
      asn = mean(nAt), median = stats::median(nAt))
}

# ---- orchestration for the jamovi analysis ---------------------------------------------

# Credible interval of pT - pC by inverting P(pT - pC > d) (monotone in d).
bb2_interval <- function(aC, bC, aT, bT, level, gl = tanh_sinh()) {
    f <- function(d, target) bb2_prob(aC, bC, aT, bT, d, gl) - target
    a2 <- (1 - level) / 2
    c(stats::uniroot(f, c(-1, 1), target = 1 - a2, tol = 1e-9)$root,
      stats::uniroot(f, c(-1, 1), target = a2, tol = 1e-9)$root)
}

bayes_compute <- function(o) {
    model <- o$model
    t <- if (model != "bf") gsd_timing(o$timing, o$nLooks) else NULL
    gF <- if (isTRUE(o$futility)) o$gammaF else NA
    if (model == "normal") {
        if (o$theta1 == 0) stop("The effect under H1 must not be 0 (its sign gives the favourable direction).")
        if (o$sigma <= 0) stop("The standard deviation must be positive.")
        dir <- sign(o$theta1)
        th1 <- abs(o$theta1)
        prior <- bayes_prior_normal(o$prior, th1, dir * o$priorMean, o$priorSD)
        sigmaU <- if (o$arms == "two") 2 * o$sigma else o$sigma
        nk <- o$nMax * t
        des <- bayes_normal_design(nk, sigmaU, prior, dir * o$delta, o$gammaE, gF, o$effInterim)
        grid <- c(0, th1 / 2, th1, 1.5 * th1)
        oc <- lapply(grid, function(th) bayes_normal_oc(des, th))
        nb <- bayes_normal_oc(des, 0, futility = FALSE)$reject
        return(list(model = model, dir = dir, theta1 = th1, prior = prior, des = des, t = t,
                    oc = data.frame(theta = dir * grid, reject = vapply(oc, `[[`, 0, "reject"),
                                    futStop = vapply(oc, `[[`, 0, "futStop"), en = vapply(oc, `[[`, 0, "en")),
                    ocLooks = oc, alphaNB = nb))
    }
    if (model == "binary") {
        if (o$pC <= 0 || o$pC >= 1 || o$pT <= 0 || o$pT >= 1 || o$pC == o$pT)
            stop("The response rates must lie strictly between 0 and 1 and differ (their order gives the favourable direction).")
        dir <- sign(o$pT - o$pC)
        nk <- round(o$bNMax * t)
        if (any(nk < 1) || (length(nk) > 1 && any(diff(nk) <= 0)))
            stop("The looks must fall on increasing numbers of patients per arm: use fewer looks or a larger maximum.")
        # lower rates favourable: swap the arms so that the rule reads P(pC - pT > delta)
        des <- bb2_design(nk, o$priorA, o$priorB, o$delta, o$gammaE, gF, o$effInterim)
        rates <- function(pT) if (dir > 0) c(o$pC, pT) else c(pT, o$pC)
        w <- abs(o$pT - o$pC)
        grid <- c(o$pC, o$pC + dir * w / 2, o$pT, min(0.999, max(0.001, o$pT + dir * w / 2)))
        oc <- lapply(grid, function(p) { r <- rates(p); bb2_oc(des, r[1], r[2]) })
        r0 <- rates(o$pC)
        nb <- bb2_oc(des, r0[1], r0[2], futility = FALSE)$reject
        return(list(model = model, dir = dir, des = des, t = t, nk = nk,
                    oc = data.frame(theta = grid, reject = vapply(oc, `[[`, 0, "reject"),
                                    futStop = vapply(oc, `[[`, 0, "futStop"), en = 2 * vapply(oc, `[[`, 0, "en")),
                    ocLooks = oc, alphaNB = nb))
    }
    two <- o$bfSamples == "two"
    if (o$bfNMin < 2 || o$bfNMax <= o$bfNMin) stop("The maximum n must exceed the minimum n (at least 2).")
    sims <- lapply(c(0, o$d1), function(dd) bf_simulate(dd, two, o$rScale, o$bfK, o$bfNMin, o$bfNMax, o$nsim, o$seed))
    nEffMax <- if (two) o$bfNMax / 2 else o$bfNMax
    list(model = model, two = two, sims = sims,
         maxBF01 = exp(-jzs_logbf(0, nEffMax, if (two) 2 * o$bfNMax - 2 else o$bfNMax - 1, o$rScale)))
}

# Monitoring a normal-approximation design: n, est, se (natural scale) per look.
bayes_monitor_normal <- function(res, o, n, est, se, level) {
    ok <- is.finite(n) & is.finite(est) & is.finite(se) & se > 0
    if (!any(ok)) stop("Interim data: no row has the cumulative N, the estimate and its standard error.")
    n <- n[ok]; est <- est[ok]; se <- se[ok]
    if (any(diff(n) <= 0)) stop("Interim data: the cumulative N must increase from row to row.")
    d <- res$des; dir <- res$dir
    e <- dir * est; I <- 1 / se^2
    IK <- I * pmax(o$nMax / n, 1)
    post <- bayes_post_normal(e, I, d$prior)
    z <- stats::qnorm(1 - (1 - level) / 2)
    pDelta <- stats::pnorm((post$mean - d$delta) / post$sd)
    pp <- bayes_pp_normal(e, I, IK, d$prior, d$delta, d$gammaE)
    final <- n >= o$nMax * (1 - 1e-9)
    lims <- cbind(post$mean - z * post$sd, post$mean + z * post$sd) * dir
    data.frame(n = n, est = est, postMean = dir * post$mean, lower = pmin(lims[, 1], lims[, 2]),
               upper = pmax(lims[, 1], lims[, 2]), pBenefit = stats::pnorm(post$mean / post$sd),
               pDelta = pDelta, pp = ifelse(final, NA, pp), final = final,
               decision = bayes_decisions(pDelta >= d$gammaE, if (is.na(d$gammaF)) FALSE else pp < d$gammaF,
                                          final, o$effInterim))
}

# Monitoring a two-arm binary design: cumulative responses and patients per arm.
bayes_monitor_binary <- function(res, o, xC, nC, xT, nT, level) {
    ok <- is.finite(xC) & is.finite(nC) & is.finite(xT) & is.finite(nT)
    if (!any(ok)) stop("Interim data: no row has responses and patients for both arms.")
    xC <- xC[ok]; nC <- nC[ok]; xT <- xT[ok]; nT <- nT[ok]
    if (any(xC < 0 | xT < 0 | xC > nC | xT > nT) || any(abs(c(xC, nC, xT, nT) - round(c(xC, nC, xT, nT))) > 1e-9))
        stop("Interim data: responses and patients must be whole numbers, with responses at most the patients.")
    d <- res$des; a <- d$a; b <- d$b
    N <- utils::tail(d$nk, 1)
    s <- d$eff[[length(d$eff)]]
    # lower rates favourable: the arms enter the rule swapped (see bayes_compute)
    sw <- res$dir < 0
    A <- if (sw) list(x = xT, n = nT) else list(x = xC, n = nC)     # the arm expected lower
    B <- if (sw) list(x = xC, n = nC) else list(x = xT, n = nT)
    gl <- tanh_sinh()
    rows <- lapply(seq_along(xC), function(i) {
        aA <- a + A$x[i]; bA <- b + A$n[i] - A$x[i]; aB <- a + B$x[i]; bB <- b + B$n[i] - B$x[i]
        pDelta <- bb2_prob(aA, bA, aB, bB, d$delta, gl)
        pBen <- bb2_prob(aA, bA, aB, bB, 0, gl)
        final <- min(nC[i], nT[i]) >= N
        pp <- if (final) NA else bb2_pp(A$x[i], A$n[i], B$x[i], B$n[i], N, a, b, s)
        # posterior of pT - pC on the natural scale
        aC <- a + xC[i]; bC <- b + nC[i] - xC[i]; aT <- a + xT[i]; bT <- b + nT[i] - xT[i]
        ci <- bb2_interval(aC, bC, aT, bT, level, gl)
        c(postMean = aT / (aT + bT) - aC / (aC + bC), lower = ci[1], upper = ci[2], pBenefit = pBen,
          pDelta = pDelta, pp = pp, final = final)
    })
    m <- as.data.frame(do.call(rbind, rows))
    m$final <- m$final > 0
    m$decision <- bayes_decisions(m$pDelta >= d$gammaE, if (is.na(d$gammaF)) FALSE else !is.na(m$pp) & m$pp < d$gammaF,
                                  m$final, o$effInterim)
    cbind(data.frame(nC = nC, xC = xC, nT = nT, xT = xT), m)
}

# Decisions look by look: efficacy / futility stops at interim looks, success at the final.
bayes_decisions <- function(eff, fut, final, effInterim) {
    K <- length(eff); fut <- rep_len(fut, K)
    out <- character(K); stopped <- FALSE
    for (k in seq_len(K)) {
        out[k] <- if (stopped) "after" else if (final[k]) (if (eff[k]) "success" else "fail")
                  else if (eff[k] && effInterim) "efficacy" else if (fut[k]) "futility" else "continue"
        if (out[k] %in% c("success", "fail", "efficacy", "futility")) stopped <- TRUE
    }
    out
}
