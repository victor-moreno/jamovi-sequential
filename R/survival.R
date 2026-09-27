# Survival designs: patients, accrual, follow-up and the calendar time of looks.
#
# Exponential survival (control hazard lam0, treatment lam0 * hr), exponential
# dropout (hazard gam, both arms), uniform accrual over [0, A] of N patients,
# minimum follow-up F, so the final analysis is at A + F (Lachin & Foulkes 1986).
# Expected events by calendar time tau, for an arm with N_j patients:
#   E_j(tau) = N_j / A * lam / mu * [ m - (exp(-mu (tau - m)) - exp(-mu tau)) / mu ],
#   mu = lam + gam, m = min(tau, A).
# Look k is at the tau_k where total expected events reach t_k * D (Kim & Tsiatis 1990).

surv_hazard <- function(spec, median = NULL, hazard = NULL, time = NULL, surv = NULL) {
    switch(spec,
        median = log(2) / median,
        hazard = hazard,
        survival = -log(surv) / time)
}

# Dropout hazard from the probability of dropping out by `time`.
dropout_hazard <- function(prob, time) if (prob <= 0) 0 else -log(1 - prob) / time

surv_events_at <- function(tau, N, A, lam0, hr, gam, r = 1) {
    arm <- function(Nj, lam) {
        mu <- lam + gam
        m <- pmin(tau, A)
        Nj / A * lam / mu * (m - (exp(-mu * (tau - m)) - exp(-mu * tau)) / mu)
    }
    arm(N * r / (1 + r), lam0 * hr) + arm(N / (1 + r), lam0)
}

# Solve the survival design for D events.
#   solve = "n":        A and F given, find N
#   solve = "accrual":  accrual rate and F given, find A (N = rate * A)
#   solve = "followup": N and A given, find F
surv_plan <- function(D, lam0, hr, gam = 0, r = 1, solve = "n",
                      A = NULL, F = NULL, rate = NULL, N = NULL) {
    ev <- function(N, A, F) surv_events_at(A + F, N, A, lam0, hr, gam, r)
    if (solve == "n") {
        N <- D / (ev(1, A, F))
    } else if (solve == "accrual") {
        f <- function(a) ev(rate * a, a, F) - D
        hi <- 1
        while (f(hi) < 0) { hi <- hi * 2; if (hi > 1e6) stop("accrual rate too low for the required events") }
        A <- stats::uniroot(f, c(1e-9, hi), tol = 1e-10)$root
        N <- rate * A
    } else {
        f <- function(fu) ev(N, A, fu) - D
        if (f(1e6) < 0) stop("sample size too small: the required events can never be observed")
        hi <- 1
        while (f(hi) < 0) hi <- hi * 2
        F <- if (f(0) >= 0) 0 else stats::uniroot(f, c(0, hi), tol = 1e-10)$root
    }
    list(N = N, A = A, F = F, D = D)
}

# Calendar time at which expected events reach each target, under hazard ratio hr.
surv_look_times <- function(targets, N, A, F, lam0, hr, gam = 0, r = 1) {
    vapply(targets, function(dk) {
        # the target may be reached after A + F (e.g. under H0 when hr > 1): search on
        f <- function(tau) surv_events_at(tau, N, A, lam0, hr, gam, r) - dk
        hi <- A + F
        while (f(hi) < 0) { hi <- hi * 2; if (hi > 1e6 * (A + F)) return(NA_real_) }
        stats::uniroot(f, c(1e-9, hi), tol = 1e-10)$root
    }, numeric(1))
}
