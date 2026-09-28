# Biomarker-Stratified Phase II: options -> design, scenarios and operating characteristics.
# Engines: R/bmrand.R (randomized, Freidlin et al. 2012) and R/bmsingle.R (single arm).

# ---- randomized ----------------------------------------------------------------------------

# Effects on the benefit scale theta (-log HR, or the difference in response rates).
bm_rand_thetas <- function(o) {
    if (o$endpoint == "surv") {
        if (!(o$hrNoBen > o$hrBen))
            stop("The BM− hazard-ratio limits must be ordered: the enrichment limit above the all-comers limit.")
        list(thetaT = -log(o$hrTarget), thetaM = -log(o$hrMid), thetaL = -log(o$hrNoBen), thetaU = -log(o$hrBen))
    } else {
        if (!(o$diffNoBen < o$diffBen))
            stop("The BM− difference limits must be ordered: the enrichment limit below the all-comers limit.")
        list(thetaT = o$diffTarget, thetaM = o$diffMid, thetaL = o$diffNoBen, thetaU = o$diffBen)
    }
}

bm_rand_options <- function(o) {
    nm <- c("endpoint", "alpha1", "alpha2", "power", "level", "sizes", "prev", "ratioNeg", "alloc",
            "interim", "frac", "zfut", "pcPos", "pcNeg", "medPos", "medNeg", "accrual", "followUp",
            "dropoutProb", "dropoutTime")
    c(stats::setNames(lapply(nm, function(n) o[[n]]), nm), bm_rand_thetas(o))
}

bm_rand_compute <- function(o) {
    ro <- bm_rand_options(o)
    d <- bm_rand_design(ro)
    user <- if (isTRUE(o$userScen)) {
        if (o$endpoint == "surv") {
            if (!(o$userPos > 0 && o$userNeg > 0)) stop("Scenario hazard ratios must be positive.")
            c(-log(o$userPos), -log(o$userNeg))
        } else c(o$userPos, o$userNeg)
    } else NULL
    sc <- bm_rand_scenarios(ro, user)
    if (o$endpoint == "binary") {
        pe <- c(o$pcPos + sc$thP, o$pcNeg + sc$thN)
        if (any(pe <= 0 | pe >= 1)) stop("A scenario gives a response rate outside (0, 1) in the experimental arm.")
    }
    list(d = d, sc = sc, oc = bm_rand_oc(d, sc$thP, sc$thN), ro = ro)
}

# Recommendation probabilities on a grid of true effects (for the map), vectorised over the grid
# and, for natural prevalence, over the BM+ sizes.
bm_rand_grid <- function(d, thP, thN, tol = 1e-8) {
    o <- d$o
    rule <- list(alpha1 = o$alpha1, alpha2 = o$alpha2, level = o$level, thetaL = o$thetaL,
                 thetaU = o$thetaU, interim = o$interim, frac = o$frac, zfut = o$zfut)
    g <- expand.grid(thP = thP, thN = thN)
    if (o$sizes == "natural") {
        nP <- 0:d$N; w <- stats::dbinom(nP, d$N, o$prev)
        keep <- w > tol & nP > 0 & nP < d$N; nP <- nP[keep]; w <- w[keep] / sum(w[keep])
    } else { nP <- d$nPos; w <- 1 }
    G <- g[rep(seq_len(nrow(g)), each = length(nP)), ]
    inf <- bm_rand_info(d, G$thP, G$thN, rep(nP, nrow(g)))
    p <- bm_rand_probs(G$thP, G$thN, inf$IP, inf$IN, rule)
    P <- rowsum(p * rep(w, nrow(g)), rep(seq_len(nrow(g)), each = length(nP)))
    cbind(g, P[, bm_rec_names])
}

# ---- single arm ----------------------------------------------------------------------------

# Response rates of the design: binary as given; time to event as landmark rates (exponential).
bm_single_rates <- function(o) {
    r <- if (o$endpoint == "surv")
        c(p0n = ss_landmark(o$m0Neg, o$lmTime), p0p = ss_landmark(o$m0Pos, o$lmTime),
          p1n = ss_landmark(o$m1Neg, o$lmTime), p1p = ss_landmark(o$m1Pos, o$lmTime))
    else c(p0n = o$p0Neg, p0p = o$p0Pos, p1n = o$p1Neg, p1p = o$p1Pos)
    if (!(r[["p1n"]] > r[["p0n"]] && r[["p1p"]] > r[["p0p"]]))
        stop(if (o$endpoint == "surv") "Target medians must be longer than the uninteresting medians in both strata."
             else "Target rates must be above the uninteresting rates in both strata.")
    r
}

bm_user_design <- function(o) {
    d <- list(n1n = o$udN1Neg, n1p = o$udN1Pos, k1n = o$udK1Neg, k1p = o$udK1Pos, nn = o$udNNeg,
              np = o$udNPos, kn = o$udKNeg, kp = o$udKPos, nep = o$udNePos, kep = o$udKePos)
    if (d$nn < d$n1n || d$np < d$n1p || d$nep < d$n1p)
        stop("User-defined design: the totals of each route must be at least the stage-1 sizes.")
    d
}

# Automatic search limits: they contain every optimum of Parashar et al. (2016, Table IV).
bm_search_limits <- function(r, alpha, beta, nMaxNeg = 0, nMaxPos = 0) {
    c(neg = if (nMaxNeg > 0) nMaxNeg else ceiling(2.5 * ph2_nlow(r[["p0n"]], r[["p1n"]], alpha, beta)) + 2,
      pos = if (nMaxPos > 0) nMaxPos else ceiling(1.4 * ph2_nlow(r[["p0p"]], r[["p1p"]], alpha, beta / 4)) + 2)
}

# The design itself (the search is cached by the caller).
bm_single_design <- function(o) {
    r <- bm_single_rates(o)
    alpha <- o$saAlpha; beta <- 1 - o$saPower
    switch(o$saDesign,
        parashar = {
            lim <- bm_search_limits(r, alpha, beta, o$nMaxNeg, o$nMaxPos)
            if (max(lim) > 250) stop("The search is limited to 250 patients per stratum; the rates need more.")
            b <- ss_search_route(r[["p0n"]], r[["p0p"]], r[["p1n"]], r[["p1p"]], alpha, beta,
                                 lim[["neg"]], lim[["pos"]], o$saType)
            list(type = "route", d = b$d, limits = lim)
        },
        jones = list(type = "route", d = bm_user_design(o)),
        isotonic = list(type = "iso", d = bm_user_design(o)),
        parallel = {
            # Bonferroni: each stratum at alpha / 2, so P(any false go) <= alpha
            sn <- ph2_simon_pick(ph2_simon(r[["p0n"]], r[["p1n"]], alpha / 2, beta), o$saType)
            sp <- ph2_simon_pick(ph2_simon(r[["p0p"]], r[["p1p"]], alpha / 2, beta), o$saType)
            list(type = "parallel", dn = ss_from_simon(sn), dp = ss_from_simon(sp))
        })
}

# Operating characteristics at rates (pn, pp) for the design and the stratum-size option.
bm_single_oc1 <- function(des, o, r, pn, pp) {
    natural <- o$sizes == "natural"
    curtail <- o$saDesign == "parashar"
    switch(des$type,
        route = if (natural) ss_route_oc_natural(des$d, pn, pp, o$prev, r[["p0n"]], r[["p0p"]], curtail)
                else ss_route_oc(des$d, pn, pp, curtail),
        iso = ss_iso_oc(des$d, pn, pp),
        parallel = if (natural) ss_parallel_oc_natural(des$dn, des$dp, pn, pp, o$prev, r[["p0n"]], r[["p0p"]])
                   else ss_parallel_oc(des$dn, des$dp, pn, pp))
}

bm_single_scenarios <- function(o, r) {
    s <- data.frame(key = c("null", "pos", "all"), pn = c(r[["p0n"]], r[["p0n"]], r[["p1n"]]),
                    pp = c(r[["p0p"]], r[["p1p"]], r[["p1p"]]))
    if (isTRUE(o$userScen)) {
        u <- if (o$endpoint == "surv") c(ss_landmark(o$userNeg, o$lmTime), ss_landmark(o$userPos, o$lmTime))
             else c(o$userNeg, o$userPos)
        if (any(u <= 0 | u >= 1)) stop("Scenario rates must be between 0 and 1 (medians must be positive).")
        s <- rbind(s, data.frame(key = "user", pn = u[1], pp = u[2]))
    }
    s
}

bm_single_compute <- function(o, des) {
    r <- bm_single_rates(o)
    sc <- bm_single_scenarios(o, r)
    oc <- t(mapply(function(pn, pp) bm_single_oc1(des, o, r, pn, pp)[c("stop", "pos", "neg", "all", "en", "pet")],
                   sc$pn, sc$pp))
    list(des = des, r = r, sc = sc, oc = as.data.frame(oc),
         natural = o$sizes == "natural" && des$type != "iso")
}

# Maximum and stage sizes of a design, per stratum.
bm_single_sizes <- function(des) {
    if (des$type == "parallel")
        return(list(stage1 = c(neg = des$dn$n1, pos = des$dp$n1), max = c(neg = des$dn$n, pos = des$dp$n)))
    d <- des$d
    list(stage1 = c(neg = d$n1n, pos = d$n1p), unsel = c(neg = d$nn, pos = d$np), enrich = c(neg = d$n1n, pos = d$nep),
         max = max(d$nn + d$np, d$n1n + d$nep))
}

# Decision map of a two-stage single-arm design: the route after every stage-1 outcome.
bm_single_routes <- function(des, o) {
    if (des$type == "parallel") {
        g <- expand.grid(xn = 0:des$dn$n1, xp = 0:des$dp$n1)
        cn <- g$xn >= des$dn$k1; cp <- g$xp >= des$dp$k1
        g$route <- ifelse(cn & cp, "both", ifelse(cp, "pos", ifelse(cn, "neg", "stop")))
        return(g)
    }
    d <- des$d
    g <- expand.grid(xn = 0:d$n1n, xp = 0:d$n1p)
    if (des$type == "iso") {
        both <- g$xn >= d$k1n & (g$xn + g$xp) * d$n1n >= d$k1n * (d$n1n + d$n1p)
        enr <- !both & (g$xp >= d$k1p | (g$xn + g$xp) * d$n1p >= d$k1p * (d$n1n + d$n1p))
        g$route <- ifelse(both, "both", ifelse(enr, "pos", "stop"))
    } else {
        cur <- o$saDesign == "parashar"
        g$route <- ifelse(g$xn >= d$k1n, ifelse(cur & g$xn >= d$kn, "goAll", "both"),
                   ifelse(g$xp >= d$k1p, ifelse(cur & g$xp >= d$kep, "goPos", "pos"), "stop"))
    }
    g
}
