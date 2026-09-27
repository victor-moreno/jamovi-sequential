# Recursive numerical integration for group sequential tests.
#
# Implements the grid and Simpson-rule recursion of Jennison & Turnbull (2000,
# ch. 19), after Armitage, McPherson & Rowe (1969). Z_k ~ N(theta * sqrt(I_k), 1)
# with independent score increments; the only special functions needed are
# pnorm/dnorm, which keeps the module pure R (no mvtnorm, no compiled code).
#
# A look is described by three stopping rules on the Z scale:
#   up    stop above:  Z >= up          (efficacy, upper)
#   lo    stop below:  Z <= lo          (futility one-sided, or efficacy lower two-sided)
#   wedge stop inside: |Z| < wedge      (two-sided futility "inner wedge"; 0 = none)
# The continuation region is (lo, up) minus [-wedge, wedge], i.e. a union of at
# most two intervals, so one code path serves one- and two-sided designs.

# Standard grid points for a N(mean, sd^2) density (J&T 19.2, r = 18 as gsDesign).
gs_grid_points <- function(mean, sd, r = 18) {
    i <- seq_len(6 * r - 1)
    x <- ifelse(i <= r - 1, -3 - 4 * log(r / i),
         ifelse(i <= 5 * r, -3 + 3 * (i - r) / (2 * r),
                3 + 4 * log(r / (6 * r - i))))
    mean + sd * x
}

# Grid restricted to [a, b] plus Simpson weights (odd points = grid, even = midpoints).
# hmax: largest gap allowed between grid points. The J&T grid is built for a few looks; with
# many looks (an SPRT looks after every observation) the transition density between looks is
# narrower than the grid, above all in the tails, and the recursion blows up. Gaps wider than
# hmax (a third of the transition sd, see gs_advance) are filled with equally spaced points.
gs_simpson <- function(x, a, b, hmax = Inf) {
    lo <- max(a, min(x)); hi <- min(b, max(x))
    if (!(hi > lo)) return(NULL)
    p <- sort(unique(c(lo, x[x > lo & x < hi], hi)))
    d <- diff(p)
    if (any(d > hmax)) {
        extra <- unlist(lapply(which(d > hmax), function(i) {
            k <- ceiling(d[i] / hmax)
            p[i] + d[i] * seq_len(k - 1) / k
        }))
        p <- sort(c(p, extra))
    }
    m <- length(p)
    d <- diff(p)
    z <- numeric(2 * m - 1)
    w <- numeric(2 * m - 1)
    z[seq(1, 2 * m - 1, 2)] <- p
    z[seq(2, 2 * m - 2, 2)] <- (p[-1] + p[-m]) / 2
    w[seq(2, 2 * m - 2, 2)] <- 4 * d / 6
    w[seq(1, 2 * m - 1, 2)] <- c(d, 0) / 6 + c(0, d) / 6
    list(z = z, w = w)
}

# Continuation intervals of one look as a 2-column matrix (possibly 0 rows).
gs_region <- function(lo, up, wedge = 0) {
    if (is.na(wedge) || wedge <= 0) {
        iv <- rbind(c(lo, up))
    } else {
        iv <- rbind(c(lo, min(up, -wedge)), c(max(lo, wedge), up))
    }
    iv[iv[, 2] > iv[, 1], , drop = FALSE]
}

# State of the recursion after a look: sub-density on a grid of the Z scale,
# already multiplied by the Simpson weights. The start state is a point mass
# at S_0 = 0 with I_0 = 0, or at a conditioning value (z0, I0) for
# conditional power, so the first look needs no special case.
gs_start <- function(z0 = 0, I0 = 0) {
    list(z = z0, wh = 1, I = I0, z0 = z0, I0 = I0)
}

# Mean and sd of the increment-driven transition to information I.
gs_trans <- function(state, I, theta) {
    dI <- I - state$I
    if (dI <= 0) stop("information must increase between looks")
    list(mu = state$z * sqrt(state$I) + theta * dI, sd = sqrt(dI))
}

# P(Z at information I >= c, and continuation so far) under theta.
gs_upper <- function(state, I, c, theta) {
    if (c == Inf) return(0)
    tr <- gs_trans(state, I, theta)
    sum(state$wh * stats::pnorm((c * sqrt(I) - tr$mu) / tr$sd, lower.tail = FALSE))
}

# P(Z at information I <= c, and continuation so far) under theta.
gs_lower <- function(state, I, c, theta) {
    if (c == -Inf) return(0)
    tr <- gs_trans(state, I, theta)
    sum(state$wh * stats::pnorm((c * sqrt(I) - tr$mu) / tr$sd))
}

# Advance the recursion to information I over the continuation region.
gs_advance <- function(state, I, region, theta, r = 18) {
    tr <- gs_trans(state, I, theta)
    # the grid follows the unconditional distribution of Z at I given the start
    m <- (state$z0 * sqrt(state$I0) + theta * (I - state$I0)) / sqrt(I)
    s <- sqrt((I - state$I0) / I)
    x <- gs_grid_points(m, s, r)
    # sd of the transition on the Z scale: when it is narrow (many looks, e.g. an SPRT), the grid
    # must resolve it; designs with a few looks keep the validated J&T grid unchanged
    kappa <- tr$sd / sqrt(I)
    hmax <- if (kappa < 0.25) kappa / 3 else Inf
    zz <- numeric(0); ww <- numeric(0)
    for (j in seq_len(nrow(region))) {
        g <- gs_simpson(x, region[j, 1], region[j, 2], hmax)
        if (!is.null(g)) { zz <- c(zz, g$z); ww <- c(ww, g$w) }
    }
    if (length(zz) == 0)
        return(list(z = 0, wh = 0, I = I, z0 = state$z0, I0 = state$I0))
    # density of Z at I: sqrt(I)/sd * phi((z sqrt(I) - mu) / sd), summed over the old grid
    dens <- stats::dnorm(outer(zz * sqrt(I), tr$mu, "-") / tr$sd) %*% state$wh
    list(z = zz, wh = ww * as.vector(dens) * sqrt(I) / tr$sd,
         I = I, z0 = state$z0, I0 = state$I0)
}

# Exit probabilities for a whole design under theta.
#   info   information at each look (any scale; only ratios and theta*sqrt(info) matter)
#   up, lo, wedge  per-look stopping rules (see header)
#   start  optional list(z, I) to condition on an observed statistic
# Returns a data frame with P(stop above), P(stop below), P(stop in wedge) per look.
gs_probs <- function(info, up, lo = rep(-Inf, length(info)), wedge = rep(0, length(info)),
                     theta = 0, start = NULL, r = 18) {
    K <- length(info)
    st <- if (is.null(start)) gs_start() else gs_start(start$z, start$I)
    pu <- pl <- pw <- numeric(K)
    for (k in seq_len(K)) {
        u <- up[k]; l <- min(lo[k], up[k]); w <- wedge[k]
        pu[k] <- gs_upper(st, info[k], u, theta)
        pl[k] <- gs_lower(st, info[k], l, theta)
        if (!is.na(w) && w > 0) {
            wl <- max(l, -w); wu <- min(u, w)
            if (wu > wl)
                pw[k] <- gs_lower(st, info[k], wu, theta) - gs_lower(st, info[k], wl, theta)
        }
        if (k < K) st <- gs_advance(st, info[k], gs_region(l, u, w), theta, r)
    }
    data.frame(up = pu, lo = pl, wedge = pw)
}
