# Phase II single-arm designs with a binary endpoint: everything exact binomial.
#
# H0: p <= p0 (uninteresting response rate) against H1: p >= p1 (promising).
# A multistage design is described by cumulative sample sizes n_k and, per stage,
#   fut_k: stop and declare "not promising" if cumulative responses <= fut_k (-1 = never)
#   eff_k: stop and declare "promising"     if cumulative responses >= eff_k (n_k + 1 = never)
# with fut_K = eff_K - 1 at the last stage, so every path ends in a decision.

# Exact operating characteristics by recursion over the cumulative responses.
ms_oc <- function(nk, fut, eff, p) {
    K <- length(nk)
    dist <- stats::dbinom(0:nk[1], nk[1], p)
    pEff <- pFut <- numeric(K)
    for (k in seq_len(K)) {
        x <- seq_along(dist) - 1
        pEff[k] <- sum(dist[x >= eff[k]])
        pFut[k] <- sum(dist[x <= fut[k]])
        if (k == K) break
        dist[x >= eff[k] | x <= fut[k]] <- 0
        m <- nk[k + 1] - nk[k]
        dist <- stats::convolve(dist, rev(stats::dbinom(0:m, m, p)), type = "open")
        dist[dist < 0] <- 0
    }
    stop <- pEff + pFut
    list(reject = sum(pEff), pEff = pEff, pFut = pFut,
         pet = sum(stop[-K]), en = sum(nk * stop))
}

# Smallest single-stage exact design: reject H0 if responses > r.
ph2_single <- function(p0, p1, alpha, beta, nmax = 5000) {
    ph2_check_rates(p0, p1)
    for (n in seq_len(nmax)) {
        r <- stats::qbinom(1 - alpha, n, p0)
        if (stats::pbinom(r, n, p0, lower.tail = FALSE) > alpha) r <- r + 1
        if (stats::pbinom(r, n, p1, lower.tail = FALSE) >= 1 - beta)
            return(list(n = n, r = r, alpha = stats::pbinom(r, n, p0, lower.tail = FALSE),
                        power = stats::pbinom(r, n, p1, lower.tail = FALSE)))
    }
    stop("no single-stage design with at most ", nmax, " patients meets these error rates")
}

# Simon (1989) two-stage designs: stop after n1 if responses <= r1; after n, promising if > r.
# For every n, the design with the smallest expected size under p0 among those meeting
# alpha and power; optimal = smallest EN0, minimax = smallest n; admissible designs
# (Jung et al. 2004) minimise q n + (1 - q) EN0 for some weight q in [0, 1].
# Smallest n at which ANY test with at most n patients can reach the power: the
# randomised Neyman-Pearson test on the total responses is the most powerful level-alpha
# test, so no two-stage design can be smaller. (The exact single-stage n is NOT a lower
# bound: p0 = 0.2, p1 = 0.4, 5 %/80 % has a two-stage minimax design with n = 33 below it.)
ph2_check_rates <- function(p0, p1) {
    if (!(p0 > 0 && p1 < 1 && p1 > p0))
        stop("Response rates: the target rate p1 must be larger than the uninteresting rate p0.")
}

ph2_nlow <- function(p0, p1, alpha, beta, nmax = 5000) {
    ph2_check_rates(p0, p1)
    for (n in seq_len(nmax)) {
        r <- stats::qbinom(1 - alpha, n, p0)
        if (stats::pbinom(r, n, p0, lower.tail = FALSE) > alpha) r <- r + 1
        d0 <- stats::dbinom(r, n, p0)
        g <- if (d0 > 0) min(1, max(0, (alpha - stats::pbinom(r, n, p0, lower.tail = FALSE)) / d0)) else 0
        if (stats::pbinom(r, n, p1, lower.tail = FALSE) + g * stats::dbinom(r, n, p1) >= 1 - beta) return(n)
    }
    stop("these error rates need more than ", nmax, " patients.")
}

ph2_simon <- function(p0, p1, alpha, beta, nmax = NULL) {
    nlow <- ph2_nlow(p0, p1, alpha, beta)
    if (nlow > 500) stop("Simon design: at least ", nlow, " patients are needed; the search is limited to 500.")
    if (is.null(nmax) || nmax <= 0) nmax <- max(2 * nlow, nlow + 20)
    rows <- list()
    for (n in max(2, nlow):nmax) {
        best <- NULL
        for (n1 in 1:(n - 1)) {
            n2 <- n - n1
            rej <- function(p) {
                # P(X1 > r1, X1 + X2 > r) for all r1 (rows 0..n1-1) and r (cols 0..n)
                d1 <- stats::dbinom(0:n1, n1, p)
                tail2 <- c(rep(1, n1), stats::pbinom(0:n, n2, p, lower.tail = FALSE))
                idx <- outer(0:n1, 0:n, function(x1, r) r - x1 + n1 + 1)
                M <- d1 * matrix(tail2[idx], nrow = n1 + 1)
                C <- apply(M[(n1 + 1):1, , drop = FALSE], 2, cumsum)[(n1 + 1):1, , drop = FALSE]
                C[2:(n1 + 1), , drop = FALSE]
            }
            R0 <- rej(p0)
            R0[col(R0) - 1 < row(R0) - 1] <- Inf        # require r >= r1
            ok <- R0 <= alpha
            rmin <- apply(ok, 1, function(v) match(TRUE, v))
            has <- !is.na(rmin)
            if (!any(has)) next
            r1s <- which(has) - 1
            rs <- rmin[has] - 1
            R1 <- rej(p1)
            pow <- R1[cbind(r1s + 1, rs + 1)]
            good <- pow >= 1 - beta
            if (!any(good)) next
            r1s <- r1s[good]; rs <- rs[good]
            pet0 <- stats::pbinom(r1s, n1, p0)
            en0 <- n1 + (1 - pet0) * n2
            j <- which.min(en0)
            if (is.null(best) || en0[j] < best$en0 - 1e-12)
                best <- list(n = n, n1 = n1, r1 = r1s[j], r = rs[j], en0 = en0[j], pet0 = pet0[j],
                             alpha = R0[r1s[j] + 1, rs[j] + 1], power = pow[good][j])
        }
        if (!is.null(best)) rows[[length(rows) + 1]] <- as.data.frame(best)
    }
    if (length(rows) == 0) stop("no two-stage design with at most ", nmax, " patients meets these error rates")
    cand <- do.call(rbind, rows)
    cand <- cand[order(cand$n), ]
    # admissible: lower convex hull of (n, EN0), with the range of weights q for which each is best
    q <- seq(0, 1, by = 0.0005)
    loss <- outer(q, cand$n) + outer(1 - q, cand$en0)
    win <- apply(loss, 1, which.min)
    adm <- sort(unique(win))
    cand$qLo <- NA; cand$qHi <- NA
    for (i in adm) { cand$qLo[i] <- min(q[win == i]); cand$qHi[i] <- max(q[win == i]) }
    list(cand = cand, optimal = cand[which.min(cand$en0), ],
         minimax = cand[cand$n == min(cand$n), ][which.min(cand$en0[cand$n == min(cand$n)]), ],
         admissible = cand[!is.na(cand$qLo), ])
}

ph2_simon_pick <- function(s, type, q = 0.5) {
    switch(type, optimal = s$optimal, minimax = s$minimax,
           admissible = { a <- s$admissible; a[which(q >= a$qLo & q <= a$qHi)[1], ] })
}

# Fleming (1982) multistage: cumulative sizes N_g, total N, z = z_{1-alpha}, [x]* nearest integer.
ph2_fleming <- function(p0, alpha, nk) {
    N <- nk[length(nk)]
    z <- stats::qnorm(1 - alpha)
    nint <- function(x) floor(x + 0.5)
    pA <- (sqrt(N * p0) + z * sqrt(1 - p0))^2 / (N + z^2)
    eff <- nint(nk * p0 + z * sqrt(N * p0 * (1 - p0))) + 1
    fut <- nint(nk * pA - z * sqrt(N * pA * (1 - pA)))
    K <- length(nk)
    fut[K] <- eff[K] - 1
    list(nk = nk, fut = pmax(fut, -1), eff = pmin(eff, nk + 1), pA = pA)
}

# Normal-approximation single-stage N, the usual starting point for Fleming's design.
ph2_normal_n <- function(p0, p1, alpha, beta) {
    ceiling(((stats::qnorm(1 - alpha) * sqrt(p0 * (1 - p0)) +
              stats::qnorm(1 - beta) * sqrt(p1 * (1 - p1))) / (p1 - p0))^2)
}

# ---- Bayesian designs (Beta(a, b) prior) -------------------------------------------

post_prob <- function(x, n, a, b, p) stats::pbeta(p, a + x, b + n - x, lower.tail = FALSE)

# Lee & Liu (2008): probability that the trial, continued to nmax, ends with
# P(p > p0 | data) > thetaT, averaging over the beta-binomial predictive distribution.
pred_prob <- function(x, n, nmax, a, b, p0, thetaT) {
    m <- nmax - n
    y <- 0:m
    lw <- lchoose(m, y) + lbeta(a + x + y, b + n - x + m - y) - lbeta(a + x, b + n - x)
    success <- post_prob(x + y, nmax, a, b, p0) > thetaT
    sum(exp(lw)[success])
}

# Boundaries of a Bayesian design at looks nk (the last = nmax).
#   type "predictive": futility if PP < thetaL, efficacy if PP > thetaU (Lee & Liu)
#   type "posterior":  futility if P(p > p1) < thetaL, efficacy if P(p > p0) > thetaU (Thall & Simon style)
#   final: promising if P(p > p0) > thetaT
ph2_bayes <- function(p0, p1, a, b, nk, thetaT, thetaL, thetaU, type = "predictive") {
    K <- length(nk); nmax <- nk[K]
    fut <- rep(-1, K); eff <- nk + 1
    for (k in seq_len(K)) {
        n <- nk[k]; x <- 0:n
        if (k < K) {
            if (type == "predictive") {
                pp <- vapply(x, pred_prob, 0, n = n, nmax = nmax, a = a, b = b, p0 = p0, thetaT = thetaT)
                fcrit <- pp < thetaL; ecrit <- pp > thetaU
            } else {
                fcrit <- post_prob(x, n, a, b, p1) < thetaL
                ecrit <- post_prob(x, n, a, b, p0) > thetaU
            }
            if (any(fcrit)) fut[k] <- max(x[fcrit])
            if (thetaU < 1 && any(ecrit)) eff[k] <- min(x[ecrit])
        } else {
            ok <- post_prob(x, n, a, b, p0) > thetaT
            eff[k] <- if (any(ok)) min(x[ok]) else n + 1
            fut[k] <- eff[k] - 1
        }
    }
    list(nk = nk, fut = fut, eff = pmax(eff, fut + 1))
}

# ---- inference after a two-stage (Simon) design -------------------------------------
# Ordering of Koyama & Chen (2008): stage-1 stops (x1 <= r1) rank below every stage-2
# outcome; within a stage, more responses are more extreme.

ph2_twostage_tail <- function(p, x, stage, r1, n1, n, strict = FALSE) {
    x <- x + strict
    if (stage == 1) return(stats::pbinom(x - 1, n1, p, lower.tail = FALSE))
    x1 <- (r1 + 1):n1
    sum(stats::dbinom(x1, n1, p) * stats::pbinom(x - x1 - 1, n - n1, p, lower.tail = FALSE))
}

ph2_twostage_inference <- function(x, stage, r1, n1, n, p0, level = 0.95) {
    a2 <- (1 - level) / 2
    g <- function(p, strict = FALSE) ph2_twostage_tail(p, x, stage, r1, n1, n, strict)
    lo <- if (x == 0) 0 else stats::uniroot(function(p) g(p) - a2, c(1e-10, 1 - 1e-10), tol = 1e-10)$root
    ntot <- if (stage == 1) n1 else n
    hi <- if (x == ntot) 1 else stats::uniroot(function(p) 1 - g(p, TRUE) - a2, c(1e-10, 1 - 1e-10), tol = 1e-10)$root
    # UMVUE (Jung & Kim 2004): x1 / n1 at stage 1; otherwise a ratio of binomial sums
    umvue <- if (stage == 1) x / n1 else {
        n2 <- n - n1
        x1 <- max(r1 + 1, x - n2):min(n1, x)
        sum(choose(n1 - 1, x1 - 1) * choose(n2, x - x1)) / sum(choose(n1, x1) * choose(n2, x - x1))
    }
    list(p.value = g(p0), lower = lo, upper = hi, umvue = umvue, naive = x / ntot)
}

# ---- from analysis options to results -------------------------------------------------

ph2_prior <- function(o) {
    switch(o$prior,
        meanp0 = c(o$p0, 1 - o$p0), jeffreys = c(0.5, 0.5), uniform = c(1, 1),
        custom = c(o$priorA, o$priorB))
}

ph2_looks <- function(first, every, nmax) {
    if (first >= nmax) return(nmax)
    unique(c(seq(first, nmax - 1, by = max(1, every)), nmax))
}

# Design: stages (nk, fut, eff) plus what the tables report. `simon` is a cached
# ph2_simon() search (it takes seconds), or NULL to run it here.
ph2_compute <- function(o, simon = NULL) {
    ph2_check_rates(o$p0, o$p1)
    beta <- 1 - o$power
    out <- list(design = o$design)
    if (o$design == "simon") {
        s <- if (is.null(simon)) ph2_simon(o$p0, o$p1, o$alpha, beta, o$nMax) else simon
        d <- ph2_simon_pick(s, o$simonType, o$admWeight)
        if (is.null(d) || nrow(d) == 0 || is.na(d$n))
            stop("Simon design: no admissible design for weight q = ", o$admWeight, ".")
        st <- list(nk = c(d$n1, d$n), fut = c(d$r1, d$r), eff = c(d$n1 + 1, d$r + 1))
        out$simon <- s; out$chosen <- d
    } else if (o$design == "single") {
        d <- ph2_single(o$p0, o$p1, o$alpha, beta)
        st <- list(nk = d$n, fut = d$r, eff = d$r + 1)
    } else if (o$design == "fleming") {
        nk <- parse_numbers(o$flemingN, "Cumulative sample sizes")
        if (length(nk) == 0) {
            N <- ph2_normal_n(o$p0, o$p1, o$alpha, beta)
            nk <- unique(round(seq_len(o$flemingStages) / o$flemingStages * N))
        }
        if (any(diff(nk) <= 0) || any(nk <= 0) || any(nk != round(nk)))
            stop("Cumulative sample sizes: whole numbers, positive and increasing (e.g. 15, 25).")
        st <- ph2_fleming(o$p0, o$alpha, nk)
    } else {
        pr <- ph2_prior(o)
        if (any(pr <= 0)) stop("Prior: both Beta parameters must be positive.")
        nk <- ph2_looks(o$bayesFirst, o$bayesEvery, o$bayesNmax)
        st <- ph2_bayes(o$p0, o$p1, pr[1], pr[2], nk, o$thetaT, o$thetaL, o$thetaU, o$design)
        out$prior <- pr
    }
    out$stages <- st
    out$oc0 <- ms_oc(st$nk, st$fut, st$eff, o$p0)
    out$oc1 <- ms_oc(st$nk, st$fut, st$eff, o$p1)
    grid <- c(o$p0, (o$p0 + o$p1) / 2, o$p1, min(0.99, o$p1 + (o$p1 - o$p0) / 2))
    out$oc <- do.call(rbind, lapply(grid, function(p) {
        r <- ms_oc(st$nk, st$fut, st$eff, p)
        data.frame(p = p, reject = r$reject, pet = r$pet, en = r$en)
    }))
    out
}

# Observed data: n patients evaluated, x responses.
# Probability that the trial ends promising given x responses in n patients, if the true rate is
# p: exact over the remaining patients, applying the stopping rules at every later look. A look
# at exactly n has already decided (1 or 0); an unplanned n is assumed to have passed the earlier
# looks. NA beyond the maximum sample size.
ph2_cond_promising <- function(st, n, x, p) {
    K <- length(st$nk)
    k <- match(n, st$nk)
    if (!is.na(k) && (x >= st$eff[k] || x <= st$fut[k] || k == K)) return(as.numeric(x >= st$eff[k]))
    later <- which(st$nk > n)
    if (!length(later)) return(NA_real_)
    dist <- 1; prev <- n; succ <- 0
    for (k in later) {
        d <- st$nk[k] - prev
        # distribution of the extra responses since the current status, still running
        dist <- stats::convolve(dist, rev(stats::dbinom(0:d, d, p)), type = "open")
        dist[dist < 0] <- 0
        y <- x + seq_along(dist) - 1
        e <- y >= st$eff[k]
        succ <- succ + sum(dist[e])
        dist[e | y <= st$fut[k]] <- 0
        prev <- st$nk[k]
    }
    min(1, succ)
}

ph2_observe <- function(o, res, n, x, level = 0.95) {
    st <- res$stages
    if (x > n) stop("Observed data: responses cannot exceed the number of patients.")
    k <- match(n, st$nk)
    decision <- if (is.na(k)) NA_character_ else if (x >= st$eff[k]) "promising" else
        if (x <= st$fut[k]) "not promising" else "continue"
    cp <- if (n > 0) as.numeric(stats::binom.test(x, n, conf.level = level)$conf.int) else c(0, 1)
    out <- list(n = n, x = x, k = k, decision = decision, naive = if (n > 0) x / n else NA,
                lower = cp[1], upper = cp[2],
                p.value = stats::pbinom(x - 1, n, o$p0, lower.tail = FALSE), adjusted = FALSE,
                cond = vapply(c(o$p0, o$p1, if (n > 0) x / n else NA), function(p)
                    if (is.na(p)) NA_real_ else ph2_cond_promising(st, n, x, p), 0))
    if (res$design == "simon" && !is.na(k) && (k == 2 || decision == "not promising")) {
        d <- res$chosen
        inf <- ph2_twostage_inference(x, k, d$r1, d$n1, d$n, o$p0, level)
        out[c("p.value", "lower", "upper")] <- inf[c("p.value", "lower", "upper")]
        out$umvue <- inf$umvue
        out$adjusted <- TRUE
    }
    if (!is.null(res$prior)) {
        a <- res$prior[1] + x; b <- res$prior[2] + n - x
        out$postMean <- a / (a + b)
        out$crLower <- stats::qbeta((1 - level) / 2, a, b)
        out$crUpper <- stats::qbeta(1 - (1 - level) / 2, a, b)
        out$pGt0 <- stats::pbeta(o$p0, a, b, lower.tail = FALSE)
        nmax <- max(st$nk)
        out$pp <- if (n < nmax) pred_prob(x, n, nmax, res$prior[1], res$prior[2], o$p0, o$thetaT) else NA
    }
    out
}
