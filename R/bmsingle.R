# Single-arm phase II designs stratified by a biomarker (BM− / BM+), binary response, exact binomial.
#
# Route design (Jones & Holmgren 2007; Parashar et al. 2016), cut-offs "X >= k":
#   stage 1: n1n BM− and n1p BM+ patients.
#     X1n >= k1n -> unselected stage 2, up to nn BM− and np BM+ in total:
#                   Xn >= kn -> all-comers; else Xp >= kp -> BM+ only; else stop
#     else X1p >= k1p -> enrichment stage 2, up to nep BM+ in total: Xp >= kep -> BM+ only; else stop
#     else stop
#   curtail: stop at stage 1 when the final cut-off is already reached (X1n >= kn, or X1p >= kep on
#   the enrichment route); it changes the expected size, not the decisions.
# Rejection probabilities (Parashar's notation): R1 all-comers, R2 BM+ only without enrichment,
# R3 BM+ only after enrichment. With X1n, X1p independent:
#   R1 = S, R2 = (C - S) T, R3 = (1 - C) E, where C = P(X1n >= k1n), S = P(X1n >= k1n, Xn >= kn),
#   T = P(Xp >= kp), E = P(X1p >= k1p, X1p + X2 >= kep).
# Parallel designs: an independent two-stage design in each stratum.

# P(X1 >= k1, X1 + X2 >= k), X1 ~ Bin(n1, p), X2 ~ Bin(n - n1, p): vectorised over k.
ss_reject <- function(n1, k1, n, k, p) {
    x1 <- 0:n1
    d1 <- stats::dbinom(x1, n1, p)
    ok <- x1 >= k1
    if (!any(ok)) return(rep(0, length(k)))
    vapply(k, function(kk) sum(d1[ok] * stats::pbinom(kk - x1[ok] - 1, n - n1, p, lower.tail = FALSE)), 0)
}

ss_tail <- function(n, k, p) stats::pbinom(k - 1, n, p, lower.tail = FALSE)

# Exact OC of a route design at true rates (pn, pp). d: n1n, n1p, k1n, k1p, nn, np, kn, kp, nep, kep.
ss_route_oc <- function(d, pn, pp, curtail = TRUE) {
    C <- ss_tail(d$n1n, d$k1n, pn)
    S <- ss_reject(d$n1n, d$k1n, d$nn, d$kn, pn)
    T <- ss_tail(d$np, d$kp, pp)
    E <- ss_reject(d$n1p, d$k1p, d$nep, d$kep, pp)
    R1 <- S; R2 <- (C - S) * T; R3 <- (1 - C) * E
    # stage 2 is enrolled on the unselected route unless curtailed, and on the enrichment route
    goU <- if (curtail) C - ss_tail(d$n1n, max(d$k1n, d$kn), pn) else C
    cP <- ss_tail(d$n1p, d$k1p, pp)
    goE <- (1 - C) * (if (curtail) cP - ss_tail(d$n1p, max(d$k1p, d$kep), pp) else cP)
    en <- d$n1n + d$n1p + (d$nn - d$n1n + d$np - d$n1p) * goU + (d$nep - d$n1p) * goE
    c(stop = 1 - R1 - R2 - R3, pos = R2 + R3, neg = 0, all = R1,
      R1 = R1, R2 = R2, R3 = R3, en = en, pet = 1 - goU - goE)
}

# Cut-off at a realised size n for a cut-off k planned at size n0: the smallest k' whose null tail
# P(Bin(n, p0) >= k') does not exceed the planned one (natural prevalence; our extension).
ss_adapt <- function(k, n0, n, p0) {
    if (k <= 0) return(rep(0L, length(n)))
    target <- ss_tail(n0, k, p0) * (1 + 1e-9)
    vapply(n, function(m) {
        kk <- 0:(m + 1)
        as.integer(kk[which(ss_tail(m, kk, p0) <= target)[1]])
    }, 0L)
}

# Binomial weights of the BM+ count among m all-comers, trimmed below tol.
ss_split <- function(m, prev, tol = 1e-10) {
    x <- 0:m
    w <- stats::dbinom(x, m, prev)
    keep <- w > tol
    list(x = x[keep], w = w[keep] / sum(w[keep]))
}

# Route design with all-comers enrolment: stage 1 of n1n + n1p patients and the unselected stage 2
# of (nn - n1n) + (np - n1p) patients split at random; the enrichment stage 2 enrols BM+ only.
# Cut-offs re-derived at the realised sizes (ss_adapt), OC averaged exactly over the splits.
ss_route_oc_natural <- function(d, pn, pp, prev, p0n, p0p, curtail = TRUE) {
    N1 <- d$n1n + d$n1p; N2 <- d$nn - d$n1n + d$np - d$n1p; e2 <- d$nep - d$n1p
    s1 <- ss_split(N1, prev); s2 <- ss_split(N2, prev)
    acc <- c(stop = 0, pos = 0, neg = 0, all = 0, R1 = 0, R2 = 0, R3 = 0, en = 0, pet = 0)
    for (i in seq_along(s1$x)) {
        m1p <- s1$x[i]; m1n <- N1 - m1p
        k1n <- ss_adapt(d$k1n, d$n1n, m1n, p0n); k1p <- ss_adapt(d$k1p, d$n1p, m1p, p0p)
        kep <- ss_adapt(d$kep, d$nep, m1p + e2, p0p)
        C <- ss_tail(m1n, k1n, pn)
        E <- ss_reject(m1p, k1p, m1p + e2, kep, pp)
        cP <- ss_tail(m1p, k1p, pp)
        goE <- (1 - C) * (if (curtail) cP - ss_tail(m1p, max(k1p, kep), pp) else cP)
        # unselected route: the stage-2 split
        m2p <- s2$x; m2n <- N2 - m2p
        mn <- m1n + m2n; mp <- m1p + m2p
        kn <- ss_adapt(d$kn, d$nn, mn, p0n); kp <- ss_adapt(d$kp, d$np, mp, p0p)
        S <- vapply(seq_along(m2p), function(j) ss_reject(m1n, k1n, mn[j], kn[j], pn), 0)
        T <- ss_tail(mp, kp, pp)
        R1 <- sum(s2$w * S); R2 <- sum(s2$w * (C - S) * T); R3 <- (1 - C) * E
        # curtail only when rejection is certain whatever the split: the cut-off at the largest size
        knc <- ss_adapt(d$kn, d$nn, m1n + N2, p0n)
        goU <- if (curtail) C - ss_tail(m1n, max(k1n, knc), pn) else C
        en <- N1 + N2 * goU + e2 * goE
        acc <- acc + s1$w[i] * c(stop = 1 - R1 - R2 - R3, pos = R2 + R3, neg = 0, all = R1,
                                 R1 = R1, R2 = R2, R3 = R3, en = en, pet = 1 - goU - goE)
    }
    acc
}

# Parallel two-stage designs, one per stratum: stage 1 n1, continue if X1 >= k1, promising if X >= k.
ss_parallel_oc <- function(dn, dp, pn, pp) {
    rn <- ss_reject(dn$n1, dn$k1, dn$n, dn$k, pn); rp <- ss_reject(dp$n1, dp$k1, dp$n, dp$k, pp)
    cn <- ss_tail(dn$n1, dn$k1, pn); cp <- ss_tail(dp$n1, dp$k1, pp)
    c(stop = (1 - rn) * (1 - rp), pos = (1 - rn) * rp, neg = rn * (1 - rp), all = rn * rp,
      en = dn$n1 + (dn$n - dn$n1) * cn + dp$n1 + (dp$n - dp$n1) * cp, pet = (1 - cn) * (1 - cp))
}

# Parallel designs with all-comers enrolment at stage 1 (split at random); stage 2 per stratum as planned.
ss_parallel_oc_natural <- function(dn, dp, pn, pp, prev, p0n, p0p) {
    N1 <- dn$n1 + dp$n1
    s1 <- ss_split(N1, prev)
    acc <- 0
    for (i in seq_along(s1$x)) {
        m1p <- s1$x[i]; m1n <- N1 - m1p
        a <- list(n1 = m1n, k1 = ss_adapt(dn$k1, dn$n1, m1n, p0n), n = m1n + dn$n - dn$n1)
        a$k <- ss_adapt(dn$k, dn$n, a$n, p0n)
        b <- list(n1 = m1p, k1 = ss_adapt(dp$k1, dp$n1, m1p, p0p), n = m1p + dp$n - dp$n1)
        b$k <- ss_adapt(dp$k, dp$n, b$n, p0p)
        acc <- acc + s1$w[i] * ss_parallel_oc(a, b, pn, pp)
    }
    acc
}

# Simon design as a parallel stratum design (continue if X1 >= r1 + 1, promising if X >= r + 1).
ss_from_simon <- function(s) list(n1 = s$n1, k1 = s$r1 + 1, n = s$n, k = s$r + 1)

# Landmark survival rate under exponential survival: S(t) = exp(-log(2) t / median).
ss_landmark <- function(median, time) exp(-log(2) * time / median)

# Landmark status from time and event: 1 alive/free at t, 0 event by t, NA censored before t.
ss_landmark_status <- function(time, event, t) {
    ifelse(time >= t, 1L, ifelse(event == 1, 0L, NA_integer_))
}

# ---- observed data ---------------------------------------------------------------------------

# Clopper-Pearson interval (naive: ignores the two-stage rule).
ss_cp <- function(x, n, level = 0.95) {
    a <- (1 - level) / 2
    if (n == 0) return(c(NA, NA))
    c(if (x == 0) 0 else stats::qbeta(a, x, n - x + 1), if (x == n) 1 else stats::qbeta(1 - a, x + 1, n - x))
}

# Walk patients in arrival order through a route design. pos: logical BM+, y: 0/1 response.
# natural: stage 1 is the first n1n + n1p patients whatever their biomarker, the unselected stage 2
# the next (nn - n1n) + (np - n1p), and cut-offs are re-derived at the realised sizes.
ss_route_observe <- function(d, pos, y, natural = FALSE, p0n = NULL, p0p = NULL, iso = FALSE) {
    # isotonic rule: BM− is judged on min(own, pooled) counts and BM+ on max(own, pooled)
    gn <- function(xn, xp, k, n, m) if (iso) xn >= k && (xn + xp) * n >= k * (n + m) else xn >= k
    gp <- function(xn, xp, k, n, m) if (iso) xp >= k || (xn + xp) * m >= k * (n + m) else xp >= k
    idx <- seq_along(y)
    N1 <- d$n1n + d$n1p
    if (natural) {
        s1 <- idx[idx <= N1]
        s1n <- s1[!pos[s1]]; s1p <- s1[pos[s1]]
        k1n <- ss_adapt(d$k1n, d$n1n, length(s1n), p0n); k1p <- ss_adapt(d$k1p, d$n1p, length(s1p), p0p)
        done1 <- length(y) >= N1
    } else {
        s1n <- utils::head(idx[!pos], d$n1n); s1p <- utils::head(idx[pos], d$n1p)
        k1n <- d$k1n; k1p <- d$k1p
        done1 <- length(s1n) == d$n1n && length(s1p) == d$n1p
    }
    x1n <- sum(y[s1n]); x1p <- sum(y[s1p])
    out <- list(stage1 = c(nn = length(s1n), xn = x1n, np = length(s1p), xp = x1p),
                k1 = c(n = k1n, p = k1p), route = NA, decision = NA, done1 = done1, done2 = FALSE,
                s2 = c(nn = 0, xn = 0, np = 0, xp = 0), k2 = c(n = NA, p = NA))
    if (!done1) { out$decision <- "stage1"; return(out) }
    rest <- idx[idx > max(c(s1n, s1p, 0))]
    if (natural) rest <- idx[idx > N1]
    if (gn(x1n, x1p, k1n, length(s1n), length(s1p))) {
        out$route <- "unselected"
        if (natural) {
            N2 <- d$nn - d$n1n + d$np - d$n1p
            s2 <- utils::head(rest, N2); s2n <- s2[!pos[s2]]; s2p <- s2[pos[s2]]
            done2 <- length(s2) == N2
        } else {
            s2n <- utils::head(setdiff(idx[!pos], s1n), d$nn - d$n1n)
            s2p <- utils::head(setdiff(idx[pos], s1p), d$np - d$n1p)
            done2 <- length(s2n) == d$nn - d$n1n && length(s2p) == d$np - d$n1p
        }
        nN <- length(s1n) + length(s2n); nP <- length(s1p) + length(s2p)
        kn <- if (natural) ss_adapt(d$kn, d$nn, nN, p0n) else d$kn
        kp <- if (natural) ss_adapt(d$kp, d$np, nP, p0p) else d$kp
        xn <- x1n + sum(y[s2n]); xp <- x1p + sum(y[s2p])
        out$s2 <- c(nn = length(s2n), xn = xn - x1n, np = length(s2p), xp = xp - x1p)
        out$k2 <- c(n = kn, p = kp)
        out$done2 <- done2
        out$decision <- if (!done2) "stage2" else if (gn(xn, xp, kn, nN, nP)) "all"
                        else if (gp(xn, xp, kp, nN, nP)) "pos" else "stop"
    } else if (gp(x1n, x1p, k1p, length(s1n), length(s1p))) {
        out$route <- "enrichment"
        e2 <- d$nep - d$n1p
        s2p <- utils::head(setdiff(idx[pos], s1p), e2)
        nP <- length(s1p) + length(s2p)
        kep <- if (natural) ss_adapt(d$kep, d$nep, length(s1p) + e2, p0p) else d$kep
        xp <- x1p + sum(y[s2p])
        out$s2 <- c(nn = 0, xn = 0, np = length(s2p), xp = xp - x1p)
        out$k2 <- c(n = NA, p = kep)
        out$done2 <- length(s2p) == e2
        out$decision <- if (!out$done2) "stage2" else if (xp >= kep) "pos" else "stop"
    } else {
        out$route <- "stop"; out$decision <- "stop"
    }
    out
}

# Parallel designs: each stratum on its own (natural: stage 1 = the first n1 of both strata overall).
ss_parallel_observe <- function(dn, dp, pos, y, natural = FALSE, p0n = NULL, p0p = NULL) {
    idx <- seq_along(y)
    N1 <- dn$n1 + dp$n1
    one <- function(des, sel, s1, p0) {
        m1 <- length(s1)
        k1 <- if (natural) ss_adapt(des$k1, des$n1, m1, p0) else des$k1
        ntot <- m1 + des$n - des$n1
        k <- if (natural) ss_adapt(des$k, des$n, ntot, p0) else des$k
        s2 <- utils::head(setdiff(sel, s1)[setdiff(sel, s1) > max(c(s1, 0))], des$n - des$n1)
        x1 <- sum(y[s1]); x <- x1 + sum(y[s2])
        done1 <- if (natural) length(y) >= N1 else m1 == des$n1
        dec <- if (!done1) "stage1" else if (x1 < k1) "stop"
               else if (length(s2) < des$n - des$n1) "stage2" else if (x >= k) "promising" else "stop"
        list(n1 = m1, x1 = x1, n2 = length(s2), x2 = x - x1, k1 = k1, k = k, decision = dec)
    }
    if (natural) {
        s1 <- idx[idx <= N1]
        a <- one(dn, idx[!pos], s1[!pos[s1]], p0n); b <- one(dp, idx[pos], s1[pos[s1]], p0p)
    } else {
        a <- one(dn, idx[!pos], utils::head(idx[!pos], dn$n1), p0n)
        b <- one(dp, idx[pos], utils::head(idx[pos], dp$n1), p0p)
    }
    final <- all(c(a$decision, b$decision) %in% c("stop", "promising"))
    rec <- if (!final) NA else if (a$decision == "promising" && b$decision == "promising") "all"
           else if (b$decision == "promising") "pos" else if (a$decision == "promising") "neg" else "stop"
    list(neg = a, pos = b, rec = rec)
}

# ---- design search (Parashar et al. 2016 routes) ----------------------------------------------

# P(X1 >= k1, X1 + X2 >= k) for every k1 = 0..n1 (rows) and k = 0..n + 1 (columns).
ss_rejmat <- function(n1, n, p) {
    d1 <- stats::dbinom(0:n1, n1, p)
    tail2 <- stats::pbinom((-n1 - 1):(n + 1) - 1, n - n1, p, lower.tail = FALSE)
    idx <- outer(0:n1, 0:(n + 1), function(x1, k) k - x1 + n1 + 2)
    M <- d1 * matrix(tail2[idx], nrow = n1 + 1)
    apply(M[(n1 + 1):1, , drop = FALSE], 2, cumsum)[(n1 + 1):1, , drop = FALSE]
}

# Cut-offs searched for n patients: a band around the expected counts under p0 and p1
# (as in ph2hetero, widened by `band`).
ss_band <- function(n, p0, p1, band, lo = 0, hi = n + 1) {
    a <- max(lo, round(n * p0) - band); b <- min(hi, round(n * p1) + band)
    if (a > b) integer(0) else a:b
}

# Pareto filter: keep rows not dominated when every column is to be minimised.
ss_pareto <- function(M) {
    if (nrow(M) <= 1) return(rep(TRUE, nrow(M)))
    o <- do.call(order, as.data.frame(M))
    keep <- logical(nrow(M))
    kept <- M[0, , drop = FALSE]
    for (i in o) {
        x <- M[i, ]
        dom <- nrow(kept) > 0 && any(colSums(t(kept) <= x) == ncol(M))
        if (!dom) { keep[i] <- TRUE; kept <- rbind(kept, x) }
    }
    keep
}

# Optimal (smallest EN under H0) or minimax (smallest maximum N, then EN) route design subject to
# R123(p0-, p0+) <= alpha, R1(p1-) >= 1 - beta and R23(p0-, p1+) >= 1 - beta (Parashar et al. 2016),
# with early stops (curtailment). EN = N1 + N2 P(k1- <= X1- < k-) + N2E P(X1- < k1-) P(k1+ <= X1+ < kE+).
# The constraints decompose as R1 = S, R23 = (C - S) T + (1 - C) E, so BM− designs (S, C), BM+
# enrichment designs (E) and BM+ unselected tests (T) are enumerated separately, Pareto-filtered and
# combined by branch and bound on EN.
ss_search_route <- function(p0n, p0p, p1n, p1p, alpha, beta, nmaxN, nmaxP, type = "optimal", band = 2) {
    # BM− part: stage 1 (N1n, k1n) and the unselected total (Nn, kn)
    rows <- list()
    for (n1 in 1:(nmaxN - 1)) for (n in (n1 + 1):nmaxN) {
        k1s <- ss_band(n1, p0n, p1n, band, lo = 1, hi = n1); ks <- ss_band(n, p0n, p1n, band, lo = 1)
        if (!length(k1s) || !length(ks)) next
        R1 <- ss_rejmat(n1, n, p1n)[k1s + 1, ks + 1, drop = FALSE]
        ok <- R1 >= 1 - beta
        if (!any(ok)) next
        R0 <- ss_rejmat(n1, n, p0n)[k1s + 1, ks + 1, drop = FALSE]
        ok <- ok & R0 <= alpha & outer(k1s, ks, "<=")
        if (!any(ok)) next
        w <- which(ok, arr.ind = TRUE)
        rows[[length(rows) + 1]] <- cbind(n1 = n1, k1 = k1s[w[, 1]], n = n, k = ks[w[, 2]], S0 = R0[w], S1 = R1[w])
    }
    if (!length(rows)) stop("No BM− design with at most ", nmaxN, " patients reaches the power; increase the maximum.")
    N <- as.data.frame(do.call(rbind, rows))
    N$C0 <- ss_tail(N$n1, N$k1, p0n)
    N$u0 <- N$C0 - ss_tail(N$n1, pmax(N$k1, N$k), p0n)
    N$cost <- N$n1 + (N$n - N$n1) * N$u0
    grp <- split(seq_len(nrow(N)), paste(N$n1, N$k1))
    keep <- unlist(lapply(grp, function(i) i[ss_pareto(cbind(N$S0[i], N$u0[i], N$n[i]))]))
    N <- N[keep, ]
    N <- if (type == "optimal") N[order(N$cost), ] else N[order(N$n, N$cost), ]
    Nn1 <- N$n1; Nk1 <- N$k1; Nn <- N$n; Nk <- N$k; NS0 <- N$S0; NC0 <- N$C0; Nu0 <- N$u0; Ncost <- N$cost
    best <- list(en = Inf, maxn = Inf)
    better <- function(en, mx) if (type == "minimax") mx < best$maxn || (mx == best$maxn && en < best$en - 1e-9) else en < best$en - 1e-9
    for (n1p in 1:(nmaxP - 1)) {
        if (type == "optimal" && min(Ncost) + n1p >= best$en) break
        # BM+ enrichment designs (k1p, Nep, kep) and unselected BM+ tests (Np, kp)
        er <- list()
        for (ne in (n1p + 1):nmaxP) {
            k1s <- ss_band(n1p, p0p, p1p, band, lo = 1, hi = n1p); ks <- ss_band(ne, p0p, p1p, band, lo = 1)
            if (!length(k1s) || !length(ks)) next
            E0 <- ss_rejmat(n1p, ne, p0p)[k1s + 1, ks + 1, drop = FALSE]
            E1 <- ss_rejmat(n1p, ne, p1p)[k1s + 1, ks + 1, drop = FALSE]
            ok <- outer(k1s, ks, "<=") & E0 <= alpha
            w <- which(ok, arr.ind = TRUE)
            if (nrow(w)) er[[length(er) + 1]] <- cbind(ne = ne, k1 = k1s[w[, 1]], k = ks[w[, 2]], E0 = E0[w], E1 = E1[w])
        }
        if (!length(er)) next
        E <- as.data.frame(do.call(rbind, er))
        E$v0 <- ss_tail(n1p, E$k1, p0p) - ss_tail(n1p, pmax(E$k1, E$k), p0p)
        E$cost <- (E$ne - n1p) * E$v0
        E <- E[ss_pareto(cbind(E$E0, -E$E1, E$cost, E$ne)), ]
        E0 <- E$E0; E1 <- E$E1; Ecost <- E$cost; Ene <- E$ne; Ek1 <- E$k1; Ek <- E$k
        tr <- lapply(n1p:nmaxP, function(np) {
            ks <- c(ss_band(np, p0p, p1p, band, lo = 1), np + 1)
            cbind(np = np, k = ks, T0 = ss_tail(np, ks, p0p), T1 = ss_tail(np, ks, p1p))
        })
        T <- as.data.frame(do.call(rbind, tr))
        T <- T[ss_pareto(cbind(T$T0, -T$T1, T$np)), ]
        T <- T[order(T$np), ]
        Tnp <- T$np; T0 <- T$T0; T1 <- T$T1; Tk <- T$k
        cand <- if (type == "optimal") which(Ncost + n1p < best$en) else which(pmax(Nn, Nn1 + n1p + 1) <= best$maxn)
        for (i in cand) {
            B0 <- 1 - NC0[i]; A0 <- NC0[i] - NS0[i]; u0 <- Nu0[i]
            base <- Ncost[i] + n1p
            # E rows that can still meet both constraints (T0 >= 0, T1 <= 1) and improve on the best
            ok <- B0 * E0 <= alpha - NS0[i] & B0 * E1 + A0 >= 1 - beta
            if (type == "optimal") ok <- ok & base + B0 * Ecost < best$en
            else ok <- ok & Nn1[i] + Ene <= best$maxn
            if (!any(ok)) next
            e <- which(ok)
            tsel <- if (type == "optimal") which(base + u0 * (Tnp - n1p) + B0 * min(Ecost[e]) < best$en)
                    else which(Nn[i] + Tnp <= best$maxn)
            if (!length(tsel)) next
            tau0 <- (alpha - NS0[i] - B0 * E0[e]) / A0
            tau1 <- (1 - beta - B0 * E1[e]) / A0
            F <- outer(tau0, T0[tsel], ">=") & outer(tau1, T1[tsel], "<=")
            if (!any(F)) next
            en <- base + outer(B0 * Ecost[e], u0 * (Tnp[tsel] - n1p), "+")
            en[!F] <- Inf
            if (type == "minimax") {
                mx <- outer(Nn1[i] + Ene[e], Nn[i] + Tnp[tsel], pmax)
                mx[!F] <- Inf
                m <- min(mx); en[mx > m] <- Inf
            }
            j <- which.min(en)
            ee <- e[(j - 1) %% length(e) + 1]; tt <- tsel[(j - 1) %/% length(e) + 1]
            mxj <- max(Nn1[i] + Ene[ee], Nn[i] + Tnp[tt])
            if (better(en[j], mxj))
                best <- list(en = en[j], maxn = mxj, d = list(n1n = Nn1[i], n1p = n1p, k1n = Nk1[i], k1p = Ek1[ee],
                             nn = Nn[i], np = Tnp[tt], kn = Nk[i], kp = Tk[tt], nep = Ene[ee], kep = Ek[ee]))
        }
    }
    if (is.null(best$d)) stop("No design found within the maximum sizes; increase them.")
    best
}

# Isotonic stratified rule (Li & Ivanova 2024, rules as stated in the paper): the strata are compared
# with pooled counts when the observed rates violate p- <= p+ (PAVA), X~- = N- min(p^-, pbar) and
# X~+ = N+ max(p^+, pbar). Stage 1: X~1- >= k1- -> both continue; else X~1+ >= k1+ -> enrich; else stop.
# Final: X~- >= k- -> all-comers; else X~+ >= k+ -> BM+ only. Enriched: X+ >= kE+. No early stops.
# Exact: stage 2 by two-dimensional convolution (Toeplitz matrix products).
ss_iso_oc <- function(d, pn, pp) {
    geq_n <- function(xn, xp, k, n, m) xn >= k & (xn + xp) * n >= k * (n + m)
    geq_p <- function(xn, xp, k, n, m) xp >= k | (xn + xp) * m >= k * (n + m)
    d1n <- stats::dbinom(0:d$n1n, d$n1n, pn); d1p <- stats::dbinom(0:d$n1p, d$n1p, pp)
    P1 <- outer(d1n, d1p)
    X1n <- row(P1) - 1; X1p <- col(P1) - 1
    both <- geq_n(X1n, X1p, d$k1n, d$n1n, d$n1p)
    enr <- !both & geq_p(X1n, X1p, d$k1p, d$n1n, d$n1p)
    toep <- function(n1, n, p) {
        d2 <- stats::dbinom(0:(n - n1), n - n1, p)
        A <- matrix(0, n + 1, n1 + 1)
        for (j in 0:n1) A[j + 1:length(d2), j + 1] <- d2
        A
    }
    F <- toep(d$n1n, d$nn, pn) %*% (P1 * both) %*% t(toep(d$n1p, d$np, pp))
    Xn <- row(F) - 1; Xp <- col(F) - 1
    r1 <- geq_n(Xn, Xp, d$kn, d$nn, d$np)
    R1 <- sum(F[r1]); R2 <- sum(F[!r1 & geq_p(Xn, Xp, d$kp, d$nn, d$np)])
    e1 <- colSums(P1 * enr)
    R3 <- sum(e1 * ss_tail(d$nep - d$n1p, d$kep - (0:d$n1p), pp))
    pB <- sum(P1[both]); pE <- sum(e1)
    c(stop = 1 - R1 - R2 - R3, pos = R2 + R3, neg = 0, all = R1, R1 = R1, R2 = R2, R3 = R3,
      en = d$n1n + d$n1p + (d$nn - d$n1n + d$np - d$n1p) * pB + (d$nep - d$n1p) * pE, pet = 1 - pB - pE)
}
