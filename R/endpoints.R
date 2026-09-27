# Sample size of a solved design (bounds.R) for each endpoint.
#
# A design gives the drift eta (= theta1 * sqrt(I_max)) and the inflation factor
# R = (eta / etaFix)^2. The endpoint turns information into observations.
# Conventions follow rpact, which is the test oracle:
#   means        normal approximation: N = (eta / delta)^2 sigma^2 (1 + r)^2 / r
#   proportions  fixed-sample N from the pooled-H0 / unpooled-H1 formula, times R
#   generic      effect theta1 and the SD of one observation's contribution
# r is the allocation ratio treatment : control; N is the total over both groups.

size_means <- function(d, delta, sd, r = 1) {
    k <- sd^2 * (1 + r)^2 / r / delta^2
    list(Nmax = d$eta^2 * k, Nfix = d$etaFix^2 * k, Imax = (d$eta / delta)^2)
}

size_props <- function(d, p1, p0, r = 1) {
    za <- stats::qnorm(1 - d$spec$alpha / d$spec$sided)
    zb <- stats::qnorm(1 - d$spec$beta)
    pbar <- (r * p1 + p0) / (1 + r)
    Nfix <- (za * sqrt(pbar * (1 - pbar) * (1 + r)^2 / r) +
             zb * sqrt((1 + r) * (p1 * (1 - p1) / r + p0 * (1 - p0))))^2 / (p1 - p0)^2
    Nmax <- Nfix * d$inflation
    # information at N under the alternative (unpooled), for effect-scale boundaries
    list(Nmax = Nmax, Nfix = Nfix,
         Imax = Nmax / ((1 + r) * (p1 * (1 - p1) / r + p0 * (1 - p0))))
}

size_generic <- function(d, theta, sd) {
    list(Nmax = (d$eta * sd / theta)^2, Nfix = (d$etaFix * sd / theta)^2, Imax = (d$eta / theta)^2)
}

# Events for a log-rank / Cox comparison with hazard ratio hr (treatment : control).
events_needed <- function(eta, hr, r = 1, method = "schoenfeld") {
    if (method == "freedman") eta^2 * (1 + r * hr)^2 / (r * (1 - hr)^2)
    else eta^2 * (1 + r)^2 / (r * log(hr)^2)
}

# Boundaries on the effect scale: theta at which Z_k equals z at information I_k.
effect_at <- function(z, Ik) z / sqrt(Ik)
