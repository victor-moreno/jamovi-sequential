# Error-spending functions: cumulative error spent by information fraction t.
#
# `a` is the level spent in ONE tail: the one-sided alpha, or alpha/2 per tail
# of a symmetric two-sided design (the gsDesign convention, which reproduces the
# classic Lan-DeMets tables). Beta spending uses the same functions with a = beta.
#
#   LDOF        Lan & DeMets (1983), O'Brien-Fleming type: 2 - 2 Phi(z_{1-a/2} / sqrt(t))
#   LDPocock    Lan & DeMets (1983), Pocock type:          a log(1 + (e - 1) t)
#   KD          Kim & DeMets (1987) power family:           a t^rho           (param = rho)
#   HSD         Hwang, Shih & DeCani (1990):                a (1 - e^(-g t)) / (1 - e^(-g))  (param = gamma)
#   exponential Anderson & Clark (2010):                    a^(t^(-nu))       (param = nu)
#   user        cumulative proportions at the planned looks, linearly interpolated
#               (param = list(t = planned fractions, p = cumulative proportions))

sf_names <- c("LDOF", "LDPocock", "KD", "HSD", "exponential", "user")

sf_value <- function(sf, a, t, param = NULL) {
    t <- pmin(pmax(t, 0), 1)
    v <- switch(sf,
        LDOF = ifelse(t > 0, 2 - 2 * stats::pnorm(stats::qnorm(1 - a / 2) / sqrt(pmax(t, 1e-300))), 0),
        LDPocock = a * log(1 + (exp(1) - 1) * t),
        KD = a * t^param,
        HSD = if (abs(param) < 1e-10) a * t else a * (1 - exp(-param * t)) / (1 - exp(-param)),
        exponential = ifelse(t > 0, a^(pmax(t, 1e-300)^(-param)), 0),
        user = a * stats::approx(c(0, param$t), c(0, param$p), xout = t, rule = 2)$y,
        stop("unknown spending function: ", sf))
    pmin(pmax(v, 0), a)
}
