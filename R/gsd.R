# Group sequential design: from the analysis options (a plain list, as in
# jamovi's self$options) to every number the results show. No jmvcore here, so
# it is testable without jamovi and reused by the tutorial.

# jamovi shows NA in a number column as NaN; an empty string leaves the cell blank.
blank_na <- function(values) lapply(values, function(x) if (is.numeric(x) && length(x) == 1 && is.na(x)) "" else x)

# "0.3, 0.6, 1" or cumulative sizes "100, 200, 300" -> information fractions.
parse_numbers <- function(x, what = "") {
    x <- trimws(unlist(strsplit(as.character(x), "[,; ]+")))
    x <- x[nzchar(x)]
    v <- suppressWarnings(as.numeric(x))
    if (any(is.na(v)))
        stop(what, ": could not read '", paste(x, collapse = ", "), "'; use numbers separated by commas.")
    v
}

gsd_timing <- function(timing, K) {
    v <- parse_numbers(timing, "Information fractions")
    if (length(v) == 0) return(seq_len(K) / K)
    if (any(v <= 0) || (length(v) > 1 && any(diff(v) <= 0)))
        stop("Information fractions: the values must be positive and increasing ",
             "(e.g. 0.3, 0.6, 1, or cumulative sizes such as 100, 200, 300).")
    v / v[length(v)]
}

# Efficacy and futility specifications for gs_design().
gsd_eff <- function(o, K) {
    switch(o$efficacy,
        ldof = list(type = "spending", sf = "LDOF"),
        ldpocock = list(type = "spending", sf = "LDPocock"),
        kd = list(type = "spending", sf = "KD", param = o$effRho),
        hsd = list(type = "spending", sf = "HSD", param = o$effGamma),
        exponential = list(type = "spending", sf = "exponential", param = o$effNu),
        user = {
            p <- parse_numbers(o$effUser, "Cumulative proportions of α")
            if (length(p) == 0)
                stop("User-defined spending: enter the cumulative proportions of α, one per analysis and ",
                     "ending in 1, in 'Cumulative proportions of α' (for ", K, " analyses, e.g. ",
                     paste(format(round(seq_len(K)^2 / K^2, 2)), collapse = ", "), ").")
            if (length(p) != K)
                stop("User-defined spending: ", length(p), " proportions were given for ", K,
                     " analyses; give one per analysis.")
            if (any(diff(c(0, p)) < 0) || abs(p[K] - 1) > 1e-8)
                stop("User-defined spending: the cumulative proportions of α must increase and end at 1.")
            list(type = "spending", sf = "user", param = list(t = NULL, p = p))
        },
        obf = list(type = "WT", delta = 0),
        pocock = list(type = "WT", delta = 0.5),
        wt = list(type = "WT", delta = o$wtDelta),
        hp = list(type = "HP", z = o$hpZ),
        pt = list(type = "PT", delta = o$ptDelta1, delta0 = o$ptDelta0))
}

gsd_fut <- function(o, K) {
    switch(o$futility,
        none = list(type = "none"),
        spending = switch(o$futSf,
            ldof = list(type = "spending", sf = "LDOF"),
            ldpocock = list(type = "spending", sf = "LDPocock"),
            kd = list(type = "spending", sf = "KD", param = o$futRho),
            hsd = list(type = "spending", sf = "HSD", param = o$futGamma),
            exponential = list(type = "spending", sf = "exponential", param = o$futNu)),
        fixedz = {
            z <- parse_numbers(o$futZ, "Futility z-values")
            if (length(z) == 0)
                stop("Fixed z-values: enter the futility z-values, one for all interim analyses or one for each.")
            list(type = "fixedz", z = rep_len(z, max(K - 1, 1)))
        },
        cp = list(type = "cp", cp = o$futCp, under = o$futCpUnder))
}

gsd_compute <- function(o) {
    tt <- gsd_timing(o$timing, o$nLooks)
    K <- length(tt)
    sided <- if (identical(o$sided, "two")) 2 else 1
    beta <- 1 - o$power
    if (o$efficacy == "triangular") {
        d <- tri_design(timing = tt, sided = sided, alpha = o$alpha, beta = beta)
    } else {
        eff <- gsd_eff(o, K)
        if (identical(eff$sf, "user")) eff$param$t <- tt
        d <- tryCatch(gs_design(timing = tt, sided = sided, alpha = o$alpha, beta = beta,
                                eff = eff, fut = gsd_fut(o, K), binding = isTRUE(o$binding)),
                      error = function(e) stop("The boundaries cannot be computed for these settings (",
                                               conditionMessage(e), "). Check α, power and the boundary parameters."))
    }
    d$K <- K
    r <- o$ratio

    # sample size and information per look
    sz <- switch(o$endpoint,
        means = size_means(d, o$meanDiff, o$sd, r),
        props = size_props(d, o$p1, o$p0, r),
        generic = size_generic(d, o$genEffect, o$genSd),
        surv = {
            D <- events_needed(d$eta, o$hr, r, o$eventsMethod)
            Dfix <- events_needed(d$etaFix, o$hr, r, o$eventsMethod)
            lam0 <- surv_hazard(o$ctrlSpec, median = o$ctrlMedian, hazard = o$ctrlHazard,
                                time = o$survTime, surv = o$survProb)
            gam <- dropout_hazard(o$dropoutProb, o$dropoutTime)
            pl <- surv_plan(D, lam0, o$hr, gam, r, solve = o$solveFor,
                            A = o$accrualTime, F = o$followUp, rate = o$accrualRate, N = o$nTotal)
            list(Nmax = D, Nfix = Dfix, Imax = D * r / (1 + r)^2, plan = pl, lam0 = lam0, gam = gam)
        })
    theta1 <- switch(o$endpoint,
        means = o$meanDiff, props = o$p1 - o$p0, generic = o$genEffect, surv = log(o$hr))

    Ik <- tt * sz$Imax
    toEffect <- function(z) {
        e <- sign(theta1) * z / sqrt(Ik)
        if (o$endpoint == "surv") exp(e) else e
    }
    fut <- if (sided == 1) d$lo else d$wedge
    fut[K] <- NA
    hasFut <- any(is.finite(fut) & (sided == 1 | fut > 1e-5), na.rm = TRUE)
    fut[!is.finite(fut) | (sided == 2 & fut <= 1e-5)] <- NA

    effStop0 <- d$H0$up + (if (sided == 2) d$H0$lo else 0)
    futStop0 <- if (sided == 1) d$H0$lo else d$H0$wedge
    effStop1 <- d$H1$up + (if (sided == 2) d$H1$lo else 0)
    futStop1 <- if (sided == 1) d$H1$lo else d$H1$wedge
    looks <- data.frame(
        look = seq_len(K), t = tt, n = tt * sz$Nmax,
        effZ = d$up,
        effP = (if (sided == 2) 2 else 1) * stats::pnorm(d$up, lower.tail = FALSE),
        effEst = toEffect(d$up),
        futZ = fut, futEst = toEffect(fut),
        cumAlpha = d$cumAlpha,
        stopH0 = effStop0 + futStop0, stopH1 = effStop1 + futStop1)

    # operating characteristics over a grid of true effects
    grid <- c(0, 0.5, 1, 1.5)
    oc <- do.call(rbind, lapply(grid, function(g) {
        p <- gs_probs(tt, d$up, d$lo, d$wedge, theta = g * d$eta)
        rej <- p$up + (if (sided == 2) p$lo else 0)
        fstop <- if (sided == 1) p$lo else p$wedge
        data.frame(g = g,
                   effect = if (o$endpoint == "surv") exp(g * theta1) else g * theta1,
                   reject = sum(rej),
                   earlyEff = sum(rej[-K]),
                   earlyFut = sum(fstop[-K]),
                   en = sum(tt * rowSums(p)) * sz$Nmax)
    }))

    out <- list(d = d, size = sz, looks = looks, oc = oc, hasFut = hasFut, sided = sided,
                theta1 = theta1, endpoint = o$endpoint)

    if (o$endpoint == "surv") {
        pl <- sz$plan
        tH1 <- surv_look_times(tt * sz$Nmax, pl$N, pl$A, pl$F, sz$lam0, o$hr, sz$gam, r)
        tH0 <- surv_look_times(tt * sz$Nmax, pl$N, pl$A, pl$F, sz$lam0, 1, sz$gam, r)
        out$timing <- data.frame(look = seq_len(K), events = tt * sz$Nmax, timeH1 = tH1, timeH0 = tH0,
                                 enrolled = pl$N * pmin(tH1, pl$A) / pl$A)
        stop1 <- rowSums(d$H1); stop0 <- rowSums(d$H0)
        out$duration <- c(max = pl$A + pl$F, H1 = sum(tH1 * stop1), H0 = sum(tH0 * stop0, na.rm = TRUE))
    }
    out
}
