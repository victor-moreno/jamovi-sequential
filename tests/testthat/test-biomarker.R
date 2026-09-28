# Biomarker-Stratified Phase II. Reference values: Freidlin et al. (2012, J Clin Oncol) Table 2 and
# the NCI program RP2BM; Parashar et al. (2016, Pharm Stat) Tables III and IV; Li & Ivanova (2024)
# designs under the rules stated in their paper; closed forms and seeded simulation.

bm <- function(...) do.call(biomarker, utils::modifyList(list(bmLevel = NULL, armLevel = NULL, respLevel = NULL,
                                                            eventLevel = NULL), list(...), keep.null = TRUE))
# absolute agreement, for references printed to a few decimals or simulated
expect_near <- function(x, y, tol) expect_lte(max(abs(unname(x) - unname(y))), tol)
route <- function(k1n, k1p, n1n, n1p, kep, nep, kn, kp, nn, np)
    list(n1n = n1n, n1p = n1p, k1n = k1n, k1p = k1p, nn = nn, np = np, kn = kn, kp = kp, nep = nep, kep = kep)

test_that("the bivariate normal CDF matches the closed form at the origin and its margins", {
    for (r in c(-0.97, -0.6, -0.2, 0, 0.3, 0.7071, 0.93, 0.99))
        expect_equal(bvn_lower(0, 0, r), 0.25 + asin(r) / (2 * pi), tolerance = 1e-14)
    expect_equal(bvn_lower(c(-1, 0.5, 2), 30, 0.6), pnorm(c(-1, 0.5, 2)), tolerance = 1e-12)
    expect_equal(bvn_upper(0.3, -0.4, 0.8), bvn_lower(-0.3, 0.4, 0.8), tolerance = 1e-15)
    skip_if_not_installed("mvtnorm")
    h <- c(-1.3, 0.2, 1.1); k <- c(0.4, -0.8, 1.9)
    for (r in c(-0.95, 0.5, 0.95)) {
        ref <- mapply(function(a, b) mvtnorm::pmvnorm(upper = c(a, b), corr = matrix(c(1, r, r, 1), 2))[1], h, k)
        expect_equal(bvn_lower(h, k, r), ref, tolerance = 1e-10)
    }
})

test_that("the randomized design reproduces Freidlin et al. (2012) Table 2 (prevalence 33%)", {
    # BM+ sized for 90% power at HR 2 (control / experimental), stratum sizes at the natural prevalence
    o <- list(endpoint = "surv", alpha1 = 0.10, alpha2 = 0.05, power = 0.90, level = 0.80,
              thetaT = log(2), thetaM = log(1.5), thetaL = log(1.3), thetaU = log(1.5), sizes = "natural",
              prev = 1 / 3, ratioNeg = 2, alloc = 1, interim = FALSE, frac = 0.5, zfut = 0, medPos = 6, medNeg = 6,
              accrual = 24, followUp = 12, dropoutProb = 0, dropoutTime = 12)
    oc <- bm_rand_oc(bm_rand_design(o), c(0, log(2)), c(0, 0))
    # columns stop / enrich / stratified / all-comers; the paper simulated 50,000 trials
    expect_near(unlist(oc[1, bm_rec_names]), c(0.87, 0.06, 0.05, 0.03), 0.011)
    expect_near(unlist(oc[2, bm_rec_names]), c(0.09, 0.50, 0.40, 0.01), 0.011)
})

test_that("exact recommendation probabilities agree with simulation, with an interim futility look", {
    rule <- list(alpha1 = 0.10, alpha2 = 0.05, level = 0.80, thetaL = log(1.3), thetaU = log(1.5),
                 interim = TRUE, frac = 0.5, zfut = 0)
    IP <- 14; IN <- 28; thP <- log(1.5); thN <- log(1.5)
    ex <- bm_rand_probs(thP, thN, IP, IN, rule)
    set.seed(7); n <- 2e5; t <- 0.5
    z1P <- rnorm(n, thP * sqrt(t * IP)); z2P <- sqrt(t) * z1P + sqrt(1 - t) * rnorm(n, thP * sqrt((1 - t) * IP))
    z1N <- rnorm(n, thN * sqrt(t * IN)); z2N <- sqrt(t) * z1N + sqrt(1 - t) * rnorm(n, thN * sqrt((1 - t) * IN))
    cN <- z1N <= 0
    ZN <- ifelse(cN, z1N, z2N); I2 <- ifelse(cN, t * IN, IN)
    est <- ZN / sqrt(I2); hw <- qnorm(0.9) / sqrt(I2)
    Zall <- (sqrt(IP) * z2P + sqrt(I2) * ZN) / sqrt(IP + I2)
    rec <- ifelse(z2P > qnorm(0.9), ifelse(est + hw < log(1.3), "enrich", ifelse(est - hw > log(1.5), "all", "strat")),
                  ifelse(Zall > qnorm(0.95), "all", "stop"))
    rec[z1P <= 0] <- "stop"
    mc <- as.vector(table(factor(rec, bm_rec_names)) / n)
    expect_near(unname(ex[1, ]), mc, 0.005)
})

test_that("the randomized sample size follows Schoenfeld and the prevalence", {
    o <- list(endpoint = "surv", alpha1 = 0.10, alpha2 = 0.05, power = 0.90, level = 0.80,
              thetaT = log(2), thetaM = log(1.5), thetaL = log(1.3), thetaU = log(1.5), sizes = "natural",
              prev = 1 / 3, ratioNeg = 2, alloc = 1, interim = FALSE, frac = 0.5, zfut = 0, medPos = 6, medNeg = 6,
              accrual = 24, followUp = 12, dropoutProb = 0, dropoutTime = 12)
    d <- bm_rand_design(o)
    expect_equal(d$dPos, ceiling((qnorm(0.9) + qnorm(0.9))^2 * 4 / log(2)^2))
    expect_equal(d$N, ceiling(d$nPos * 3))
    oc <- bm_rand_oc(d, c(0, log(2)), c(0, 0))
    expect_equal(rowSums(oc[, bm_rec_names]), c(1, 1), tolerance = 1e-12)
    expect_gt(oc$enrich[2], 0.45)
})

test_that("observed log-rank, Cox and stratified statistics match the survival package", {
    skip_if_not_installed("survival")
    set.seed(11); n <- 240
    pos <- runif(n) < 0.4; ex <- runif(n) < 0.5
    lam <- ifelse(pos, 0.1, 0.07) * ifelse(ex, ifelse(pos, 0.5, 0.9), 1)
    tt <- rexp(n, lam); cens <- runif(n, 5, 30)
    time <- round(pmin(tt, cens), 1); event <- as.integer(tt <= cens)
    s <- bm_rand_stats("surv", pos, ex, time = time, event = event)
    for (g in c(TRUE, FALSE)) {
        st <- if (g) s$pos else s$neg
        sd <- survival::survdiff(survival::Surv(time, event) ~ ex, subset = pos == g)
        cx <- survival::coxph(survival::Surv(time, event) ~ ex, subset = pos == g, ties = "breslow")
        expect_equal(st$z^2, sd$chisq, tolerance = 1e-8)
        expect_equal(-st$theta, unname(coef(cx)), tolerance = 1e-7)
        expect_equal(st$se, sqrt(vcov(cx)[1]), tolerance = 1e-7)
    }
    sa <- survival::survdiff(survival::Surv(time, event) ~ ex + survival::strata(pos))
    expect_equal(s$all$z^2, sa$chisq, tolerance = 1e-8)
})

test_that("route designs reproduce Parashar et al. (2016) Table IV exactly", {
    tab <- list(
        list(.10, .10, c(.048, .800, .800, .623, 110.2), route(3, 2, 44, 34, 7, 104, 9, 4, 135, 53)),
        list(.10, .15, c(.049, .801, .801, .653, 77.9), route(2, 2, 32, 21, 6, 67, 7, 3, 106, 29)),
        list(.10, .25, c(.050, .800, .800, .571, 60.0), route(2, 1, 34, 8, 4, 29, 6, 2, 87, 9)),
        list(.15, .15, c(.050, .802, .801, .611, 46.9), route(2, 1, 20, 12, 4, 43, 6, 2, 66, 21)),
        list(.15, .25, c(.046, .803, .802, .561, 32.5), route(1, 1, 12, 7, 4, 28, 4, 2, 43, 11)),
        list(.15, .35, c(.045, .801, .800, .615, 27.8), route(1, 1, 11, 5, 3, 15, 4, 2, 47, 7)),
        list(.25, .25, c(.045, .802, .801, .695, 18.5), route(1, 1, 6, 6, 3, 24, 3, 2, 23, 13)),
        list(.25, .40, c(.038, .802, .801, .742, 13.5), route(1, 1, 6, 4, 2, 9, 3, 2, 23, 5)))
    for (r in tab) {
        a <- ss_route_oc(r[[4]], .03, .03); b <- ss_route_oc(r[[4]], r[[1]], r[[1]]); c <- ss_route_oc(r[[4]], .03, r[[2]])
        expect_near(c(1 - a[["stop"]], b[["R1"]], c[["R2"]] + c[["R3"]], a[["pet"]]), r[[3]][1:4], 0.0006)
        expect_near(a[["en"]], r[[3]][5], 0.051)
    }
})

test_that("the search finds Parashar's optimal designs", {
    for (r in list(list(.25, .40, 30, 30, route(1, 1, 6, 4, 2, 9, 3, 2, 23, 5)),
                   list(.25, .25, 40, 40, route(1, 1, 6, 6, 3, 24, 3, 2, 23, 13)),
                   list(.15, .35, 60, 30, route(1, 1, 11, 5, 3, 15, 4, 2, 47, 7)))) {
        b <- ss_search_route(.03, .03, r[[1]], r[[2]], 0.05, 0.2, r[[3]], r[[4]])
        expect_equal(b$d[names(r[[5]])], r[[5]], ignore_attr = TRUE)
    }
})

test_that("the Jones & Holmgren design reproduces Parashar et al. (2016) Table III", {
    d <- route(2, 1, 34, 14, 5, 50, 4, 4, 53, 27)
    expect_near(1 - ss_route_oc(d, .03, .03, curtail = FALSE)[["stop"]], 0.079, 0.0006)
    for (s in list(c(.03, .03, .067, .012, 65.79), c(.03, .15, .067, .720, 80.21), c(.10, .15, .755, NA, 80.03),
                   c(.15, .30, .952, NA, 80.10))) {
        x <- ss_route_oc(d, s[1], s[2], curtail = FALSE)
        expect_near(x[["R1"]], s[3], 0.0006)
        expect_near(x[["en"]], s[5], 0.02)
        if (!is.na(s[4])) expect_near(x[["R2"]] + x[["R3"]], s[4], 0.0006)
    }
    # Parashar's design without early stops (Li & Ivanova's comparator), p0 0.20, p1 0.40
    e <- route(6, 4, 19, 14, 15, 47, 19, 9, 60, 27)
    e00 <- ss_route_oc(e, .2, .2, FALSE); e01 <- ss_route_oc(e, .2, .4, FALSE); e11 <- ss_route_oc(e, .4, .4, FALSE)
    expect_near(c(1 - e00[["stop"]], e00[["pet"]], e00[["en"]], 0.5 * e00[["en"]] + 0.25 * e01[["en"]] + 0.25 * e11[["en"]]),
                c(.0497, .584, 50.14, 62.30), 0.006)
})

test_that("the isotonic rules give the type I errors of Li & Ivanova's designs under their stated rules", {
    tab <- list(list(.15, .0543, route(2, 6, 14, 28, 16, 70, 12, 14, 50, 60)),
                list(.20, .0592, route(2, 6, 10, 20, 17, 55, 11, 12, 34, 39)),
                list(.25, .0599, route(2, 7, 10, 20, 22, 62, 16, 15, 42, 41)),
                list(.30, .0541, route(5, 8, 15, 20, 26, 63, 18, 18, 42, 44)))
    for (r in tab) expect_near(1 - ss_iso_oc(r[[3]], r[[1]], r[[1]])[["stop"]], r[[2]], 6e-5)
    d <- tab[[2]][[3]]
    set.seed(4); n <- 1e5; pn <- .25; pp <- .4
    X1n <- rbinom(n, d$n1n, pn); X1p <- rbinom(n, d$n1p, pp)
    gn <- function(a, b, k, u, v) a >= k & (a + b) * u >= k * (u + v)
    gp <- function(a, b, k, u, v) b >= k | (a + b) * v >= k * (u + v)
    both <- gn(X1n, X1p, d$k1n, d$n1n, d$n1p); enr <- !both & gp(X1n, X1p, d$k1p, d$n1n, d$n1p)
    Xn <- X1n + rbinom(n, d$nn - d$n1n, pn); Xp <- X1p + rbinom(n, d$np - d$n1p, pp); Xe <- X1p + rbinom(n, d$nep - d$n1p, pp)
    rec <- ifelse(both, ifelse(gn(Xn, Xp, d$kn, d$nn, d$np), "all", ifelse(gp(Xn, Xp, d$kp, d$nn, d$np), "pos", "stop")),
                  ifelse(enr & Xe >= d$kep, "pos", "stop"))
    ex <- ss_iso_oc(d, pn, pp)
    expect_near(unname(ex[c("stop", "pos", "all")]), as.vector(table(factor(rec, c("stop", "pos", "all"))) / n), 0.005)
})

test_that("natural prevalence: exact averaging over the stratum sizes agrees with simulation", {
    d <- route(2, 2, 32, 21, 6, 67, 7, 3, 106, 29)
    ex <- ss_route_oc_natural(d, .10, .10, prev = 21 / 53, p0n = .03, p0p = .03)
    set.seed(2); n <- 2e4
    N1 <- 53; N2 <- 82; e2 <- 46
    out <- character(n)
    for (i in seq_len(n)) {
        m1p <- rbinom(1, N1, 21 / 53); m1n <- N1 - m1p
        x1n <- rbinom(1, m1n, .10); x1p <- rbinom(1, m1p, .10)
        out[i] <- if (x1n >= ss_adapt(2, 32, m1n, .03)) {
            m2p <- rbinom(1, N2, 21 / 53); m2n <- N2 - m2p
            xn <- x1n + rbinom(1, m2n, .10); xp <- x1p + rbinom(1, m2p, .10)
            if (xn >= ss_adapt(7, 106, m1n + m2n, .03)) "all" else if (xp >= ss_adapt(3, 29, m1p + m2p, .03)) "pos" else "stop"
        } else if (x1p >= ss_adapt(2, 21, m1p, .03)) {
            if (x1p + rbinom(1, e2, .10) >= ss_adapt(6, 67, m1p + e2, .03)) "pos" else "stop"
        } else "stop"
    }
    mc <- as.vector(table(factor(out, c("stop", "pos", "all"))) / n)
    expect_near(unname(ex[c("stop", "pos", "all")]), mc, 0.012)
    expect_true(all(ss_adapt(7, 106, 60:140, .03) >= 0) && !is.unsorted(ss_adapt(7, 106, 60:140, .03)))
})

test_that("landmark rates and statuses", {
    expect_equal(ss_landmark(3, 6), 0.25)
    expect_equal(ss_landmark_status(c(7, 2, 3), c(0, 1, 0), 6), c(1L, 0L, NA))
})

test_that("the analysis runs in every mode and analyses the example datasets", {
    r <- bm(data = data.frame())
    oc <- r$oc$asDF
    expect_equal(oc$stop + oc$enrich + oc$strat + oc$all, rep(1, 4), tolerance = 1e-9)
    expect_equal(r$sizes$asDF$events[1], 55)
    rn <- bm(data = data.frame(), sizes = "natural", endpoint = "binary", interim = TRUE, userScen = TRUE,
             userPos = 0.25, userNeg = 0.05)
    expect_equal(nrow(rn$oc$asDF), 5)
    rd <- read.csv(test_path("..", "..", "data", "biomarker_randomized.csv"), stringsAsFactors = TRUE)
    a <- bm(data = rd, biomarker = "biomarker", bmLevel = "Positive", arm = "arm", armLevel = "Experimental",
            time = "pfs", event = "progression", eventLevel = "Yes")
    ob <- a$observed$asDF
    expect_equal(ob$est[1], 0.4428487, tolerance = 1e-6)
    expect_equal(ob$z[1], 2.8636932, tolerance = 1e-6)
    expect_equal(a$decision$asDF$result[4], "Biomarker-stratified phase III")
    sa <- read.csv(test_path("..", "..", "data", "biomarker_singlearm.csv"), stringsAsFactors = TRUE)
    s <- bm(data = sa, trial = "single", endpoint = "binary", saDesign = "jones", biomarker = "biomarker",
            bmLevel = "Positive", response = "response", respLevel = "Yes")
    expect_equal(s$decision$asDF$result[3], "Promising in BM+ only: phase III in BM+ (enrichment)")
    expect_equal(s$observed$asDF$x2[2], 17)
    p <- bm(data = data.frame(), trial = "single", endpoint = "binary", saDesign = "parallel")
    expect_equal(nrow(p$oc$asDF), 3)
    iso <- bm(data = data.frame(), trial = "single", endpoint = "binary", saDesign = "isotonic", sizes = "natural")
    expect_equal(nrow(iso$oc$asDF), 3)
})
