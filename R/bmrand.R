# Randomized phase II trial stratified by a biomarker (BM+ / BM−), with the decision algorithm of
# Freidlin, McShane, Polley & Korn (2012, J Clin Oncol 30:3304) and its constants editable.
#
# Every effect is on a benefit scale theta (> 0 favours the experimental arm):
#   survival  theta = -log(HR), HR = experimental : control (the gsd convention; HR < 1 is benefit)
#   binary    theta = pE - pC
# In stratum s the statistic is Z_s ~ N(theta_s sqrt(I_s), 1), and the strata are independent.
# Rule at the final analysis:
#   Z+ > z(1 - alpha1) ?  BM− effect from its two-sided `level` CI:
#                           upper < thetaL -> enrichment (BM+ only); lower > thetaU -> all-comers;
#                           otherwise -> biomarker-stratified phase III
#                      :  overall (inverse-variance) Z > z(1 - alpha2) -> all-comers, else stop
# Optional futility look in both strata at information fraction t: Z_s1 <= zf in BM+ stops the
# trial; in BM− it closes the stratum, which then enters the final rule with its interim data.

bm_rec_names <- c("stop", "enrich", "strat", "all")

# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
gauss_legendre <- function(n) {
    i <- seq_len(n - 1)
    J <- matrix(0, n, n)
    J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- i / sqrt(4 * i^2 - 1)
    e <- eigen(J, symmetric = TRUE)
    list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

bvn_gl <- local({
    rules <- NULL
    function() {
        if (is.null(rules)) rules <<- lapply(c(6, 12, 20), function(n) {
            g <- gauss_legendre(n); keep <- g$x > 0
            list(x = g$x[keep], w = g$w[keep])
        })
        rules
    }
})

# P(X > h, Y > k) for a standard bivariate normal with correlation r (scalar), vectorised over
# h and k: Drezner & Wesolowsky (1990) as refined by Genz (2004, Stat Comput 14:251).
bvn_upper <- function(h, k, r) {
    n <- max(length(h), length(k))
    h <- rep_len(h, n); k <- rep_len(k, n)
    out <- numeric(n)
    inf_h <- is.infinite(h); inf_k <- is.infinite(k)
    out[(h == Inf) | (k == Inf)] <- 0
    both <- h == -Inf & k == -Inf; out[both] <- 1
    sel <- h == -Inf & !inf_k; out[sel] <- stats::pnorm(-k[sel])
    sel <- k == -Inf & !inf_h; out[sel] <- stats::pnorm(-h[sel])
    fin <- !inf_h & !inf_k
    if (!any(fin)) return(out)
    h <- h[fin]; k <- k[fin]
    rl <- bvn_gl()
    g <- if (abs(r) < 0.3) rl[[1]] else if (abs(r) < 0.75) rl[[2]] else rl[[3]]
    x <- g$x; w <- g$w
    hk <- h * k
    if (abs(r) < 0.925) {
        hs <- (h^2 + k^2) / 2
        asr <- asin(r) / 2
        bvn <- 0
        for (s in c(-1, 1)) {
            sn <- sin(asr * (1 + s * x))
            bvn <- bvn + colSums(w * exp((outer(sn, hk) - rep(hs, each = length(sn))) / (1 - sn^2)))
        }
        res <- bvn * asr / (2 * pi) + stats::pnorm(-h) * stats::pnorm(-k)
    } else {
        if (r < 0) { k <- -k; hk <- -hk }
        bvn <- numeric(length(h))
        if (abs(r) < 1) {
            as <- (1 - r) * (1 + r); a <- sqrt(as); bs <- (h - k)^2
            cc <- (4 - hk) / 8; d <- (12 - hk) / 16
            asr <- -(bs / as + hk) / 2
            bvn <- ifelse(asr > -100, a * exp(asr) * (1 - cc * (bs - as) * (1 - d * bs / 5) / 3 + cc * d * as^2 / 5), 0)
            b <- sqrt(bs)
            sp <- sqrt(2 * pi) * stats::pnorm(-b / a)
            bvn <- bvn - ifelse(hk > -100, exp(-hk / 2) * sp * b * (1 - cc * bs * (1 - d * bs / 5) / 3), 0)
            a2 <- a / 2
            for (i in seq_along(x)) for (s in c(-1, 1)) {
                xs <- (a2 + a2 * s * x[i])^2
                rs <- sqrt(1 - xs)
                asr <- -(bs / xs + hk) / 2
                sp <- 1 + cc * xs * (1 + d * xs)
                ep <- exp(-hk * xs / (2 * (1 + rs)^2)) / rs
                bvn <- bvn + ifelse(asr > -100, a2 * w[i] * exp(asr) * (ep - sp), 0)
            }
            bvn <- -bvn / (2 * pi)
        }
        res <- if (r > 0) bvn + stats::pnorm(-pmax(h, k))
               else ifelse(h >= k, -bvn,
                           ifelse(h < 0, stats::pnorm(k) - stats::pnorm(h), stats::pnorm(-h) - stats::pnorm(-k)) - bvn)
    }
    out[fin] <- pmin(1, pmax(0, res))
    out
}

# P(X <= h, Y <= k).
bvn_lower <- function(h, k, r) bvn_upper(-h, -k, r)

# Recommendation probabilities for elementwise vectors of true effects and information.
# rule: alpha1, alpha2, level, thetaL, thetaU, interim (logical), frac, zfut.
bm_rand_probs <- function(thP, thN, IP, IN, rule, nq = 96) {
    m <- max(length(thP), length(thN), length(IP), length(IN))
    thP <- rep_len(thP, m); thN <- rep_len(thN, m); IP <- rep_len(IP, m); IN <- rep_len(IN, m)
    c1 <- stats::qnorm(1 - rule$alpha1); c2 <- stats::qnorm(1 - rule$alpha2)
    zg <- stats::qnorm((1 + rule$level) / 2)
    interim <- isTRUE(rule$interim)
    t <- rule$frac; zf <- rule$zfut; rho <- sqrt(t)
    m2P <- thP * sqrt(IP)
    if (interim) {
        m1P <- thP * sqrt(t * IP)
        pFutP <- stats::pnorm(zf - m1P)
        pPass <- bvn_upper(zf - m1P, c1 - m2P, rho)
        dens <- function(z) stats::dnorm(z - m2P) *
            stats::pnorm((m1P + rho * (z - m2P) - zf) / sqrt(1 - rho^2))
    } else {
        pFutP <- rep(0, m)
        pPass <- stats::pnorm(m2P - c1)
        dens <- function(z) stats::dnorm(z - m2P)
    }
    # BM− ends either closed at the interim (information t I−) or at the final analysis
    m2N <- thN * sqrt(IN)
    br <- if (interim) {
        m1N <- thN * sqrt(t * IN)
        pc <- stats::pnorm(zf - m1N)
        list(list(I = t * IN, mass = pc, G = function(x) stats::pnorm(pmin(x, zf) - m1N)),
             list(I = IN, mass = 1 - pc, G = function(x) stats::pnorm(x - m2N) - bvn_lower(zf - m1N, x - m2N, rho)))
    } else list(list(I = IN, mass = rep(1, m), G = function(x) stats::pnorm(x - m2N)))
    # quadrature over Z+ <= c1 for the overall test (the only term coupling the strata)
    gl <- gauss_legendre(nq)
    lo <- pmin(c1, m2P - 9)
    half <- (c1 - lo) / 2
    Z <- outer(half, gl$x) + (c1 + lo) / 2
    W <- outer(half, gl$w) * matrix(dens(Z), m)
    enrich <- strat <- all <- numeric(m)
    for (b in br) {
        a <- rule$thetaL * sqrt(b$I) - zg
        u <- rule$thetaU * sqrt(b$I) + zg
        Ga <- b$G(a); Gu <- b$G(u)
        enrich <- enrich + pPass * Ga
        strat <- strat + pPass * (Gu - Ga)
        all <- all + pPass * (b$mass - Gu)
        hz <- (c2 * sqrt(IP + b$I) - sqrt(IP) * Z) / sqrt(b$I)
        tail <- b$mass - matrix(b$G(as.vector(hz)), m)
        all <- all + rowSums(W * tail)
    }
    res <- cbind(stop = pmax(0, 1 - enrich - strat - all), enrich = enrich, strat = strat, all = all)
    # futility stops: BM+ stops the trial (both strata at the interim); BM− closes its stratum
    pFutN <- if (interim) br[[1]]$mass else rep(0, m)
    attr(res, "fut") <- cbind(pos = pFutP, neg = pFutN)
    res
}

# Information of one stratum. Survival: events d with allocation r (log-rank, Schoenfeld).
# Binary: n patients, rates pC and pC + theta, unpooled variance.
bm_info <- function(endpoint, size, theta = 0, pC = NULL, r = 1) {
    if (endpoint == "surv") size * r / (1 + r)^2
    else {
        pE <- pC + theta
        nE <- size * r / (1 + r); nC <- size / (1 + r)
        1 / (pE * (1 - pE) / nE + pC * (1 - pC) / nC)
    }
}

# Design: sample size per stratum from the BM+ power target, then BM− by prevalence or ratio.
# o: endpoint, alpha1, alpha2, power, level, thetaT (target BM+), thetaM (intermediate), thetaL, thetaU,
#    sizes ("fixed" / "natural"), prev, ratioNeg, alloc (r), interim, frac, zfut;
#    binary: pcPos, pcNeg; survival: medPos, medNeg, accrual, followUp, dropoutProb, dropoutTime.
bm_rand_design <- function(o) {
    r <- o$alloc
    if (!(o$thetaT > 0)) stop("The target effect in BM+ must favour the experimental arm.")
    if (!(o$thetaL < o$thetaU)) stop("The BM− thresholds must be ordered: the 'no benefit' limit below the 'benefit' limit.")
    if (o$sizes == "natural" && !(o$prev > 0 && o$prev < 1)) stop("The prevalence of BM+ must be between 0 and 1.")
    zsum <- stats::qnorm(1 - o$alpha1) + stats::qnorm(o$power)
    d <- list(o = o, r = r)
    if (o$endpoint == "surv") {
        if (!(o$accrual > 0) || o$followUp < 0) stop("Accrual must be positive and follow-up not negative.")
        gam <- dropout_hazard(o$dropoutProb, o$dropoutTime)
        lamP <- log(2) / o$medPos; lamN <- log(2) / o$medNeg
        q <- function(lam, th) surv_events_at(o$accrual + o$followUp, 1, o$accrual, lam, exp(-th), gam, r)
        # event probability per patient at the final analysis, under the design effects
        d$qPos <- q(lamP, o$thetaT); d$qNeg <- q(lamN, o$thetaM)
        dPos <- ceiling(zsum^2 * (1 + r)^2 / (r * o$thetaT^2) - 1e-9)
        nPos <- ceiling(dPos / d$qPos - 1e-9)
    } else {
        pE <- o$pcPos + o$thetaT
        if (!(pE < 1 && o$pcPos > 0)) stop("Response rates must stay between 0 and 1.")
        nPos <- ceiling(zsum^2 * (1 + r) * (pE * (1 - pE) / r + o$pcPos * (1 - o$pcPos)) / o$thetaT^2 - 1e-9)
    }
    if (o$sizes == "natural") {
        N <- ceiling(nPos / o$prev - 1e-9)
        nNeg <- N - nPos
    } else {
        nNeg <- ceiling(o$ratioNeg * nPos - 1e-9)
        N <- nPos + nNeg
    }
    d$N <- N; d$nPos <- nPos; d$nNeg <- nNeg
    if (o$endpoint == "surv") { d$dPos <- dPos; d$dNeg <- nNeg * d$qNeg }
    d
}

# Information of both strata for BM+ sizes nP (vector) out of N (natural) or the fixed sizes.
bm_rand_info <- function(d, thP, thN, nP = d$nPos) {
    o <- d$o
    nN <- if (o$sizes == "natural") d$N - nP else d$nNeg
    if (o$endpoint == "surv") {
        # fixed: BM+ events as designed; natural: events proportional to the realised stratum sizes
        eP <- if (o$sizes == "natural") nP * d$qPos else d$dPos
        eN <- nN * d$qNeg
        list(IP = bm_info("surv", eP, r = d$r), IN = bm_info("surv", eN, r = d$r), sizeP = eP, sizeN = eN)
    } else
        list(IP = bm_info("bin", nP, thP, o$pcPos, d$r), IN = bm_info("bin", nN, thN, o$pcNeg, d$r),
             sizeP = nP, sizeN = nN)
}

# Operating characteristics at true effects (thP[i], thN[i]): recommendation probabilities and the
# expected events (survival) or patients (binary). Natural prevalence: exact mixture over N+.
bm_rand_oc <- function(d, thP, thN, tol = 1e-10) {
    o <- d$o
    rule <- list(alpha1 = o$alpha1, alpha2 = o$alpha2, level = o$level, thetaL = o$thetaL,
                 thetaU = o$thetaU, interim = o$interim, frac = o$frac, zfut = o$zfut)
    one <- function(thp, thn) {
        if (o$sizes == "natural") {
            nP <- 0:d$N
            w <- stats::dbinom(nP, d$N, o$prev)
            keep <- w > tol & nP > 0 & nP < d$N
            nP <- nP[keep]; w <- w[keep] / sum(w[keep])
        } else { nP <- d$nPos; w <- 1 }
        inf <- bm_rand_info(d, thp, thn, nP)
        p <- bm_rand_probs(thp, thn, inf$IP, inf$IN, rule)
        f <- attr(p, "fut")
        t <- if (isTRUE(o$interim)) o$frac else 1
        # BM+ futility stops both strata at the interim; otherwise BM− may close alone
        es <- f[, "pos"] * t * (inf$sizeP + inf$sizeN) +
              (1 - f[, "pos"]) * (inf$sizeP + inf$sizeN * (t * f[, "neg"] + (1 - f[, "neg"])))
        c(colSums(w * p), expected = sum(w * es))
    }
    out <- t(mapply(one, rep_len(thP, max(length(thP), length(thN))), thN))
    as.data.frame(out)
}

# The four scenarios of Freidlin et al. (no effect; BM+ only; equal; BM+ larger), plus a user pair.
bm_rand_scenarios <- function(o, user = NULL) {
    s <- data.frame(key = c("null", "pos", "equal", "larger"),
                    thP = c(0, o$thetaT, o$thetaM, o$thetaT), thN = c(0, 0, o$thetaM, o$thetaM))
    if (!is.null(user)) s <- rbind(s, data.frame(key = "user", thP = user[1], thN = user[2]))
    s
}

# ---- observed data ---------------------------------------------------------------------------

# Log-rank for one stratum: experimental arm observed minus expected events and its variance.
bm_logrank <- function(time, event, exp) {
    tt <- sort(unique(time[event == 1]))
    if (length(tt) == 0) return(c(OmE = 0, V = 0, events = 0))
    ord <- order(time)
    time <- time[ord]; event <- event[ord]; exp <- exp[ord]
    OmE <- 0; V <- 0
    for (u in tt) {
        risk <- time >= u
        n <- sum(risk); n1 <- sum(risk & exp)
        dd <- sum(time == u & event == 1); d1 <- sum(time == u & event == 1 & exp)
        OmE <- OmE + d1 - dd * n1 / n
        if (n > 1) V <- V + dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    }
    c(OmE = OmE, V = V, events = sum(event))
}

# Cox model with one binary covariate (experimental arm), Breslow ties, Newton-Raphson.
bm_cox1 <- function(time, event, exp) {
    x <- as.numeric(exp)
    tt <- sort(unique(time[event == 1]))
    if (length(tt) == 0 || all(x == x[1])) return(c(beta = NA, se = NA))
    dd <- vapply(tt, function(u) sum(time == u & event == 1), 0)
    s1 <- vapply(tt, function(u) sum(x[time == u & event == 1]), 0)
    atrisk <- outer(time, tt, ">=")
    n0 <- colSums(atrisk); n1 <- colSums(atrisk * x)
    b <- 0
    for (it in 1:50) {
        e1 <- n1 * exp(b); e0 <- n0 - n1 + e1
        pr <- e1 / e0
        U <- sum(s1 - dd * pr); H <- sum(dd * pr * (1 - pr))
        if (H <= 0) return(c(beta = NA, se = NA))
        step <- U / H
        b <- b + step
        if (abs(step) < 1e-10) break
        if (abs(b) > 20) return(c(beta = NA, se = NA))
    }
    e1 <- n1 * exp(b); pr <- e1 / (n0 - n1 + e1)
    c(beta = b, se = 1 / sqrt(sum(dd * pr * (1 - pr))))
}

# Per-stratum and overall statistics from patient-level data.
# y: response (0/1) for binary; time + event for survival. exp: logical, experimental arm.
bm_rand_stats <- function(endpoint, pos, exp, y = NULL, time = NULL, event = NULL, level = 0.95) {
    one <- function(sel) {
        e <- exp[sel]
        if (endpoint == "surv") {
            lr <- bm_logrank(time[sel], event[sel], e)
            cx <- bm_cox1(time[sel], event[sel], e)
            list(n = sum(sel), nE = sum(e), nC = sum(!e), eventsE = sum(event[sel][e]),
                 eventsC = sum(event[sel][!e]), OmE = lr[["OmE"]], V = lr[["V"]],
                 z = if (lr[["V"]] > 0) -lr[["OmE"]] / sqrt(lr[["V"]]) else NA,
                 theta = -cx[["beta"]], se = cx[["se"]])
        } else {
            yE <- y[sel][e]; yC <- y[sel][!e]
            pE <- mean(yE); pC <- mean(yC)
            se <- sqrt(pE * (1 - pE) / length(yE) + pC * (1 - pC) / length(yC))
            list(n = sum(sel), nE = length(yE), nC = length(yC), xE = sum(yE), xC = sum(yC),
                 pE = pE, pC = pC, theta = pE - pC, se = se, z = if (se > 0) (pE - pC) / se else NA)
        }
    }
    sp <- one(pos); sn <- one(!pos)
    # overall: stratified log-rank, or the inverse-variance combination of the two differences
    if (endpoint == "surv") {
        V <- sp$V + sn$V
        zAll <- if (V > 0) -(sp$OmE + sn$OmE) / sqrt(V) else NA
        cxAll <- c(theta = NA, se = NA)
        w <- 1 / c(sp$se, sn$se)^2
        thAll <- sum(w * c(sp$theta, sn$theta)) / sum(w); seAll <- 1 / sqrt(sum(w))
    } else {
        w <- 1 / c(sp$se, sn$se)^2
        thAll <- sum(w * c(sp$theta, sn$theta)) / sum(w); seAll <- 1 / sqrt(sum(w))
        zAll <- thAll / seAll
    }
    # interaction on the benefit scale, and Gail & Simon's (1985) test for a qualitative interaction
    dInt <- sp$theta - sn$theta; seInt <- sqrt(sp$se^2 + sn$se^2)
    zz <- c(sp$theta / sp$se, sn$theta / sn$se)
    Q <- min(sum(zz[zz > 0]^2), sum(zz[zz < 0]^2))
    pGS <- if (is.na(Q)) NA else if (Q > 0) 0.5 * stats::pchisq(Q, 1, lower.tail = FALSE) else 1
    list(pos = sp, neg = sn, all = list(theta = thAll, se = seAll, z = zAll),
         inter = list(theta = dInt, se = seInt, z = dInt / seInt, Q = Q, pGS = pGS), level = level)
}

# The decision path for observed statistics (final analysis).
bm_rand_decide <- function(s, o) {
    c1 <- stats::qnorm(1 - o$alpha1); c2 <- stats::qnorm(1 - o$alpha2)
    zg <- stats::qnorm((1 + o$level) / 2)
    step1 <- isTRUE(s$pos$z > c1)
    lo <- s$neg$theta - zg * s$neg$se; hi <- s$neg$theta + zg * s$neg$se
    rec <- if (step1) {
        if (isTRUE(hi < o$thetaL)) "enrich" else if (isTRUE(lo > o$thetaU)) "all" else "strat"
    } else if (isTRUE(s$all$z > c2)) "all" else "stop"
    list(step1 = step1, c1 = c1, c2 = c2, negLo = lo, negHi = hi, overall = isTRUE(s$all$z > c2), rec = rec)
}
