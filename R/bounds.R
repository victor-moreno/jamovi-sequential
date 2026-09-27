# Group sequential boundaries on the Z scale.
#
# Designs are solved on the standardised scale: information fractions t
# (t_K = 1 at the design stage) and drift eta = theta1 * sqrt(I_max), so that
# Z_k ~ N(theta/theta1 * eta * sqrt(t_k), 1). Any endpoint then gets
# I_max = (eta / theta1)^2 (endpoints.R); the inflation factor is
# (eta / eta_fixed)^2.
#
# spec (see gs_design):
#   sided   1 or 2 (two-sided = symmetric; alpha is the total two-sided level)
#   eff     list(type = "spending", sf, param)   error spending (spending.R)
#           list(type = "WT", delta)             Wang & Tsiatis (1987); 0.5 Pocock, 0 O'Brien-Fleming
#           list(type = "HP", z)                 Haybittle (1971) - Peto et al. (1976)
#           list(type = "PT", delta, delta0)     Pampallona & Tsiatis (1994), with its own futility
#   fut     list(type = "none")
#           list(type = "spending", sf, param)   beta spending (Pampallona, Tsiatis & Kim 2001)
#           list(type = "fixedz", z)             user z-bounds at the interim looks (|Z| for two-sided)
#           list(type = "cp", cp, under)         conditional-power threshold, under "trend" or "H1"
#   binding TRUE: efficacy bounds account for futility stopping; FALSE: they ignore it

# Efficacy bound at one look spending `target` under H0 (both tails if two-sided).
gs_solve_eff <- function(state, I, target, sided) {
    if (!(target > 1e-15)) return(Inf)
    f <- if (sided == 1) function(b) gs_upper(state, I, b, 0) - target
         else function(b) gs_upper(state, I, b, 0) + gs_lower(state, I, -b, 0) - target
    lo <- if (sided == 1) -10 else 0
    if (f(lo) <= 0) return(lo)
    stats::uniroot(f, c(lo, 40), tol = 1e-12)$root
}

# Futility bound spending `target` under the alternative (drift eta), capped at `cap`.
gs_solve_fut <- function(state, I, target, cap, sided, eta) {
    if (sided == 1) {
        if (!(target > 1e-15)) return(-Inf)
        f <- function(a) gs_lower(state, I, a, eta) - target
        hi <- min(cap, 40)
        if (f(hi) <= 0) return(hi)
        return(stats::uniroot(f, c(-40, hi), tol = 1e-12)$root)
    }
    if (!(target > 1e-15)) return(0)
    f <- function(a) gs_lower(state, I, a, eta) - gs_lower(state, I, -a, eta) - target
    hi <- min(cap, 40)
    if (f(hi) <= 0) return(hi)
    stats::uniroot(f, c(0, hi), tol = 1e-12)$root
}

# Futility z-bound from a conditional-power threshold, assuming only the final
# analysis remains (B-value algebra, Lan & Wittes 1988). under = "trend" uses the
# current estimate as the drift; "H1" uses the design drift eta.
gs_cp_bound <- function(t, cp, bK, eta, under) {
    q <- stats::qnorm(1 - cp)
    if (under == "trend") sqrt(t) * (bK - sqrt(1 - t) * q)
    else (bK - eta * (1 - t) - sqrt(1 - t) * q) / sqrt(t)
}

# One pass over the looks for a given drift eta (and constant C for WT/PT).
gs_seq <- function(t, eta, spec, C = NA, bK = NA) {
    K <- length(t); sided <- spec$sided; eff <- spec$eff; fut <- spec$fut
    a1 <- spec$alpha / sided
    cumA <- switch(eff$type,
        spending = c(sided * sf_value(eff$sf, a1, t[-K], eff$param), spec$alpha),
        HP = c(rep(NA, K - 1), spec$alpha),
        rep(NA, K))
    bfix <- switch(eff$type,
        spending = rep(NA, K),
        HP = c(rep(eff$z, K - 1), NA),
        C * t^(eff$delta - 0.5))
    cumB <- if (fut$type == "spending") {
        if (!is.null(fut$cumB)) fut$cumB
        else c(sf_value(fut$sf, spec$beta, t[-K], fut$param), spec$beta)
    } else NULL

    up <- rep(Inf, K); lo <- rep(-Inf, K); wedge <- rep(0, K)
    s0 <- gs_start(); s1 <- gs_start()
    spent0 <- 0; spentB <- 0; aK <- NA
    for (k in seq_len(K)) {
        b <- bfix[k]
        if (is.na(b)) b <- gs_solve_eff(s0, t[k], cumA[k] - spent0, sided)
        spent0 <- spent0 + gs_upper(s0, t[k], b, 0) +
            (if (sided == 2) gs_lower(s0, t[k], -b, 0) else 0)
        up[k] <- b
        if (sided == 2) lo[k] <- -b

        a <- NA
        if (fut$type == "spending") {
            # the final look is not capped, so the outer search can match a_K to b_K
            a <- if (k < K && isTRUE(fut$skip[k])) 0
                 else gs_solve_fut(s1, t[k], cumB[k] - spentB, if (k < K) b else Inf, sided, eta)
            if (k == K) aK <- a
        } else if (k < K) {
            a <- switch(fut$type,
                none = NA,
                fixedz = fut$z[k],
                cp = gs_cp_bound(t[k], fut$cp, bK, eta, fut$under),
                pt = eta * sqrt(t[k]) - (eta - C) * t[k]^(eff$delta0 - 0.5))
            if (!is.na(a)) a <- min(a, b)
        }
        if (k < K && !is.na(a)) {
            if (sided == 1) lo[k] <- a else wedge[k] <- max(0, a)
        }
        if (fut$type == "spending")
            spentB <- spentB + (if (sided == 1) gs_lower(s1, t[k], min(a, b), eta)
                                else gs_lower(s1, t[k], min(a, b), eta) - gs_lower(s1, t[k], -min(a, b), eta))
        if (k < K) {
            r0 <- if (isTRUE(spec$binding)) gs_region(lo[k], up[k], wedge[k])
                  else gs_region(if (sided == 2) -b else -Inf, b, 0)
            s0 <- gs_advance(s0, t[k], r0, 0)
            s1 <- gs_advance(s1, t[k], gs_region(lo[k], up[k], wedge[k]), eta)
        }
    }
    # final look: every remaining path stops; decide by b_K
    if (sided == 1) lo[K] <- up[K] else wedge[K] <- up[K]
    list(up = up, lo = lo, wedge = wedge, alphaH0 = spent0, aK = aK, bK = up[K])
}

# Two-sided beta spending follows rpact (>= 3.3.0): the one-sided design at
# alpha/2 is solved first, and looks where its futility bound would be negative
# get no inner wedge. Beta spending is then rescaled to start at the first look
# that has one (rpact's default betaAdjustment = TRUE).
gs_twosided_beta <- function(t, alpha, beta, eff, fut, binding) {
    K <- length(t)
    aux <- gs_design(timing = t, sided = 1, alpha = alpha / 2, beta = beta,
                     eff = eff, fut = fut, binding = binding)
    f1 <- aux$lo[-K]
    cumB <- c(sf_value(fut$sf, beta, t[-K], fut$param), beta)
    kMin <- which(f1 > 0)[1]
    if (!is.na(kMin) && kMin > 1) {
        b0 <- cumB[kMin - 1]
        adj <- (cumB - b0) / (cumB[K] - b0) * cumB[K]
        cumB <- c(rep(0, kMin - 1), adj[kMin:K])
    }
    fut$skip <- f1 < 0
    fut$cumB <- cumB
    fut
}

gs_power <- function(d, eta) sum(gs_probs(d$t, d$up, d$lo, d$wedge, theta = eta)$up)

# Solve a design: bounds and the drift eta that gives power 1 - beta.
gs_design <- function(K = 3, timing = NULL, sided = 1, alpha = 0.025, beta = 0.1,
                      eff = list(type = "spending", sf = "LDOF"),
                      fut = list(type = "none"), binding = FALSE) {
    t <- if (is.null(timing)) seq_len(K) / K else timing / timing[length(timing)]
    K <- length(t)
    if (K > 1 && any(diff(t) <= 0)) stop("information fractions must increase")
    # keep the futility rule as given: monitoring re-derives the timing-dependent parts
    fut0 <- fut
    if (eff$type == "PT") fut <- list(type = "pt")
    if (sided == 2 && fut$type == "spending") fut <- gs_twosided_beta(t, alpha, beta, eff, fut, binding)
    spec <- list(sided = sided, alpha = alpha, beta = beta, eff = eff, fut = fut, fut0 = fut0,
                 binding = binding && fut$type != "none")
    etaFix <- stats::qnorm(1 - alpha / sided) + stats::qnorm(1 - beta)

    # constant C of the classical WT / PT families, solved for the alpha level
    needC <- eff$type %in% c("WT", "PT")
    solveC <- function(eta, bK) {
        if (!needC) return(NA)
        stats::uniroot(function(C) gs_seq(t, eta, spec, C, bK)$alphaH0 - alpha,
                       c(0.1, 20), tol = 1e-12)$root
    }
    # C depends on eta only when futility is binding
    C0 <- if (needC && !spec$binding) solveC(etaFix, NA) else NA

    run <- function(eta) {
        C <- if (needC && spec$binding) NA else C0
        bK <- NA
        # CP futility needs b_K, which (when binding) depends on the futility bounds
        for (it in 1:30) {
            if (needC && spec$binding) C <- solveC(eta, bK)
            s <- gs_seq(t, eta, spec, C, bK)
            if (fut$type != "cp" || (!is.na(bK) && abs(s$bK - bK) < 1e-9)) break
            bK <- s$bK
        }
        s$C <- C; s$t <- t
        s
    }

    if (fut$type == "spending") {
        eta <- stats::uniroot(function(e) { s <- run(e); s$aK - s$bK },
                              c(0.8, 1.6) * etaFix, extendInt = "upX", tol = 1e-10)$root
    } else {
        eta <- stats::uniroot(function(e) gs_power(run(e), e) - (1 - beta),
                              c(0.8, 1.6) * etaFix, extendInt = "upX", tol = 1e-10)$root
    }
    d <- run(eta)
    gs_finish(d, eta, etaFix, spec)
}

# Attach operating characteristics to a solved design.
gs_finish <- function(d, eta, etaFix, spec) {
    t <- d$t; K <- length(t)
    p0 <- gs_probs(t, d$up, d$lo, d$wedge, theta = 0)
    p1 <- gs_probs(t, d$up, d$lo, d$wedge, theta = eta)
    # H0 efficacy exits ignoring futility (the alpha a non-binding design controls)
    effLo <- if (spec$sided == 2) -d$up else c(rep(-Inf, K - 1), d$up[K])
    p0nb <- gs_probs(t, d$up, effLo, rep(0, K), theta = 0)
    cumAlpha <- cumsum(p0nb$up + (if (spec$sided == 2) p0nb$lo else 0))
    futH1 <- if (spec$sided == 1) p1$lo else p1$wedge
    futH1[K] <- 0
    list(t = t, up = d$up, lo = d$lo, wedge = d$wedge, sided = spec$sided,
         eta = eta, etaFix = etaFix, inflation = (eta / etaFix)^2, C = d$C,
         H0 = p0, H1 = p1, cumAlpha = cumAlpha, cumBetaEarly = cumsum(futH1),
         power = sum(p1$up), alpha = sum(p0nb$up + (if (spec$sided == 2) p0nb$lo else 0)),
         alphaBinding = sum(p0$up + (if (spec$sided == 2) p0$lo else 0)),
         EtH0 = sum(t * rowSums(p0)), EtH1 = sum(t * rowSums(p1)), spec = spec)
}

# Whitehead's triangular test (one-sided) and double triangular test (two-sided).
# On the (V, Z_score) plane, with theta1 = 1 so that V is in units of 1/theta1^2:
#   k = 1 + z_{1-beta} / z_{1-alpha1},  a = k ln(1/(2 alpha1)),  c = 1/(2k),  V_max = a/c
#   reject H0: Z_score >= a + cV;  accept H0: Z_score <= -a + 3cV
# (alpha1 = alpha one-sided, alpha/2 per tail two-sided). Discrete looks use the
# "Christmas-tree" correction a - 0.583 sqrt(V_k - V_{k-1}) (Whitehead 1997), and
# V_max is chosen so that the corrected lines meet at the final look, where the
# test rejects when Z_score >= 2cV. Attained alpha and power are computed
# exactly by gs_finish, not taken from the approximations.
tri_design <- function(K = 5, timing = NULL, sided = 1, alpha = 0.025, beta = 0.1,
                       correction = TRUE) {
    t <- if (is.null(timing)) seq_len(K) / K else timing / timing[length(timing)]
    K <- length(t)
    a1 <- alpha / sided
    k <- 1 + stats::qnorm(1 - beta) / stats::qnorm(1 - a1)
    a <- k * log(1 / (2 * a1)); cc <- 1 / (2 * k)
    # the corrected boundaries must meet at the last look (as in SAS PROC SEQDESIGN),
    # so every planned look is reachable: solve c V_K = a - 0.583 sqrt(V_K - V_{K-1})
    dK <- t[K] - (if (K > 1) t[K - 1] else 0)
    Vmax <- if (correction)
        stats::uniroot(function(v) cc * v - a + 0.583 * sqrt(v * dK), c(0, a / cc), tol = 1e-12)$root
    else a / cc
    spec <- list(sided = sided, alpha = alpha, beta = beta, eff = list(type = "triangular"),
                 fut = list(type = "triangular"), binding = TRUE)
    d <- c(list(t = t, C = NA), tri_bounds(t, Vmax, a, cc, sided, correction))
    etaFix <- stats::qnorm(1 - a1) + stats::qnorm(1 - beta)
    out <- gs_finish(d, sqrt(Vmax), etaFix, spec)
    out$triangle <- c(a = a, c = cc, Vmax = Vmax, Vapex = a / cc, correction = correction)
    out
}

# Triangular boundaries on the Z scale at looks t (fractions of Vmax; the last is final).
tri_bounds <- function(t, Vmax, a, cc, sided, correction = TRUE) {
    K <- length(t)
    V <- t * Vmax
    ak <- if (correction) a - 0.583 * sqrt(diff(c(0, V))) else rep(a, K)
    up <- (ak + cc * V) / sqrt(V)
    inner <- (-ak + 3 * cc * V) / sqrt(V)
    up[K] <- 2 * cc * sqrt(V[K])
    inner <- pmin(inner, up)
    lo <- rep(-Inf, K); wedge <- rep(0, K)
    if (sided == 1) {
        lo <- inner; lo[K] <- up[K]
    } else {
        lo <- -up; wedge <- pmax(0, inner); wedge[K] <- up[K]
    }
    list(up = up, lo = lo, wedge = wedge)
}
