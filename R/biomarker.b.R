biomarkerClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "biomarkerClass",
    inherit = biomarkerBase,
    private = list(

        .init = function() {
            o <- self$options
            r <- self$results
            rand <- o$trial == "randomized"
            for (key in private$.summaryKeys())
                r$summary$addRow(rowKey = key, values = list(name = private$.summaryLabel(key)))
            sz <- r$sizes
            strata <- if (rand) c(pos = .("BM+"), neg = .("BM−"), total = .("Total"))
                      else c(neg = .("BM−"), pos = .("BM+"), total = .("Total"))
            for (key in names(strata)) sz$addRow(rowKey = key, values = list(stratum = strata[[key]]))
            par <- !rand && o$saDesign == "parallel"
            for (cn in c("patients", "exp", "ctrl")) sz$getColumn(cn)$setVisible(rand)
            sz$getColumn("events")$setVisible(rand && o$endpoint == "surv")
            for (cn in c("stage1", "total")) sz$getColumn(cn)$setVisible(!rand)
            for (cn in c("unsel", "enrich")) sz$getColumn(cn)$setVisible(!rand && !par)
            for (k in seq_len(private$.nRules())) r$rules$addRow(rowKey = k, values = list())
            for (key in private$.scenarioKeys())
                r$oc$addRow(rowKey = key, values = list(scenario = private$.scenarioLabel(key)))
            oc <- r$oc
            oc$getColumn("strat")$setVisible(rand)
            oc$getColumn("neg")$setVisible(par)
            oc$getColumn("pet")$setVisible(!rand)
            oc$getColumn("expected")$setTitle(if (rand && o$endpoint == "surv") .("Expected events") else .("Expected N"))
            if (!is.null(o$biomarker)) {
                ob <- r$observed
                keys <- if (rand) c(pos = .("BM+"), neg = .("BM−"), all = .("Overall (stratified)"), inter = .("Interaction (BM+ vs BM−)"))
                        else c(neg = .("BM−"), pos = .("BM+"))
                for (key in names(keys)) ob$addRow(rowKey = key, values = list(stratum = keys[[key]]))
                for (cn in c("nE", "nC", "xE", "xC", "z", "p")) ob$getColumn(cn)$setVisible(rand)
                for (cn in c("n1", "x1", "n2", "x2")) ob$getColumn(cn)$setVisible(!rand)
                ob$getColumn("est")$setTitle(if (rand) (if (o$endpoint == "surv") .("HR") else .("Difference")) else .("Rate"))
                ob$getColumn("xE")$setSuperTitle(if (o$endpoint == "surv") .("Events") else .("Responses"))
                ob$getColumn("xC")$setSuperTitle(if (o$endpoint == "surv") .("Events") else .("Responses"))
                for (key in if (rand) c("s1", "s2", "s3", "rec") else c("s1", "s2", "rec"))
                    r$decision$addRow(rowKey = key, values = list())
            }
        },

        .run = function() {
            o <- self$options
            if (o$trial == "randomized") private$.runRandomized() else private$.runSingle()
        },

        # ---- randomized ----------------------------------------------------------------------

        .runRandomized = function() {
            o <- self$options
            r <- self$results
            res <- tryCatch(bm_rand_compute(o), error = function(e) e)
            if (inherits(res, "error")) return(private$.fail(conditionMessage(res)))
            d <- res$d; ro <- res$ro
            surv <- o$endpoint == "surv"
            private$.fillSummary(private$.randSummary(res))
            nat <- o$sizes == "natural"
            r1 <- d$r / (1 + d$r)
            sizes <- list(pos = if (nat) d$N * o$prev else d$nPos, neg = if (nat) d$N * (1 - o$prev) else d$nNeg)
            sizes$total <- d$N
            ev <- list(pos = if (surv) (if (nat) d$N * o$prev * d$qPos else d$dPos) else NA,
                       neg = if (surv) sizes$neg * d$qNeg else NA)
            ev$total <- if (surv) ev$pos + ev$neg else NA
            for (key in c("pos", "neg", "total"))
                r$sizes$setRow(rowKey = key, values = blank_na(list(patients = sizes[[key]], exp = sizes[[key]] * r1,
                                                                    ctrl = sizes[[key]] * (1 - r1), events = ev[[key]])))
            note <- if (nat) sprintf(.("Natural prevalence: %d patients enrolled whatever their biomarker; the stratum sizes shown are expected values (BM+ sized for the power target at the expected size)."), d$N)
                    else sprintf(.("Fixed stratum sizes: about %d patients are screened to fill both strata (prevalence %s)."),
                                 ceiling(max(d$nPos / o$prev, d$nNeg / (1 - o$prev))), o$prev)
            if (surv) note <- paste(note, sprintf(.("Events expected at the final analysis (accrual %s + follow-up %s) under the target HR in BM+ and the intermediate HR in BM−; the operating characteristics treat these numbers of events as fixed."),
                                                  o$accrual, o$followUp))
            r$sizes$setNote("n", note, init = FALSE)
            private$.fillRules(private$.randRules(res))
            if (o$ocTable) {
                for (i in seq_len(nrow(res$sc))) {
                    key <- res$sc$key[i]
                    r$oc$setRow(rowKey = key, values = list(effPos = private$.effText(res$sc$thP[i]),
                        effNeg = private$.effText(res$sc$thN[i]), stop = res$oc$stop[i], enrich = res$oc$enrich[i],
                        strat = res$oc$strat[i], all = res$oc$all[i], neg = "", pet = "", expected = res$oc$expected[i]))
                }
                note <- .("Exact probabilities (normal approximation of the test statistics; independent strata).")
                if (isTRUE(o$interim)) note <- paste(note, .("The interim futility look is our extension: Freidlin et al. recommend it but did not specify or simulate it."))
                r$oc$setNote("m", note, init = FALSE)
            }
            obs <- private$.observeRandomized(res)
            img <- r$regionPlot
            if (o$regionPlot && is.null(img$state)) img$setState(private$.randRegionState(res, obs))
            img <- r$ocMap
            if (o$ocMap && is.null(img$state)) {
                th <- if (surv) -log(exp(seq(log(0.3), log(1.2), length.out = 25))) else seq(-0.1, min(0.5, 0.95 - max(o$pcPos, o$pcNeg)), length.out = 25)
                g <- bm_rand_grid(d, th, th)
                img$setState(list(kind = "rand", endpoint = o$endpoint, g = g,
                                  sc = res$sc[, c("key", "thP", "thN")]))
            }
        },

        .observeRandomized = function(res) {
            o <- self$options
            r <- self$results
            if (is.null(o$biomarker)) return(NULL)
            surv <- o$endpoint == "surv"
            need <- c(o$biomarker, o$arm, if (surv) c(o$time, o$event) else o$response)
            lev <- c(o$bmLevel, o$armLevel, if (surv) o$eventLevel else o$respLevel)
            if (is.null(o$arm) || (surv && (is.null(o$time) || is.null(o$event))) || (!surv && is.null(o$response)) ||
                length(lev) < 3 || any(!nzchar(lev))) {
                private$.notice(if (surv) .("Assign the biomarker, the arm, the time and the event, and choose their levels (BM+, experimental arm, event).")
                                else .("Assign the biomarker, the arm and the response, and choose their levels (BM+, experimental arm, response)."),
                                "dataNeeded", jmvcore::NoticeType$INFO, at = match("observed", self$results$itemNames))
                r$observed$setVisible(FALSE); r$decision$setVisible(FALSE)
                return(NULL)
            }
            dat <- self$data
            pos <- as.character(dat[[o$biomarker]]) == o$bmLevel
            ex <- as.character(dat[[o$arm]]) == o$armLevel
            ok <- !is.na(pos) & !is.na(ex)
            if (surv) {
                tt <- jmvcore::toNumeric(dat[[o$time]]); ev <- as.character(dat[[o$event]]) == o$eventLevel
                ok <- ok & !is.na(tt) & !is.na(ev)
            } else {
                y <- as.character(dat[[o$response]]) == o$respLevel
                ok <- ok & !is.na(y)
            }
            s <- tryCatch({
                if (sum(ok & pos) < 2 || sum(ok & !pos) < 2) stop(.("Each stratum needs patients in both arms."))
                if (any(tapply(ex[ok], pos[ok], function(v) length(unique(v))) < 2)) stop(.("Each stratum needs patients in both arms."))
                if (surv) bm_rand_stats("surv", pos[ok], ex[ok], time = tt[ok], event = as.integer(ev[ok]), level = o$ciLevel)
                else bm_rand_stats("binary", pos[ok], ex[ok], y = as.integer(y[ok]), level = o$ciLevel)
            }, error = function(e) e)
            if (inherits(s, "error")) {
                private$.notice(conditionMessage(s), "dataProblem", jmvcore::NoticeType$WARNING,
                                at = match("observed", self$results$itemNames))
                r$observed$setVisible(FALSE); r$decision$setVisible(FALSE)
                return(NULL)
            }
            zc <- stats::qnorm((1 + o$ciLevel) / 2)
            row <- function(st, events = TRUE) {
                eff <- private$.effPoint(st$theta)
                lo <- private$.effPoint(st$theta - zc * st$se); hi <- private$.effPoint(st$theta + zc * st$se)
                if (surv) { t <- lo; lo <- hi; hi <- t }
                list(est = eff, lower = lo, upper = hi, z = st$z,
                     p = if (is.na(st$z)) NA else stats::pnorm(st$z, lower.tail = FALSE))
            }
            for (key in c("pos", "neg")) {
                st <- s[[key]]
                v <- c(list(n = st$n, nE = st$nE, nC = st$nC,
                            xE = if (surv) st$eventsE else st$xE, xC = if (surv) st$eventsC else st$xC), row(st))
                r$observed$setRow(rowKey = key, values = blank_na(c(v, list(n1 = NA, x1 = NA, n2 = NA, x2 = NA))))
            }
            a <- s$all
            r$observed$setRow(rowKey = "all", values = blank_na(c(list(n = s$pos$n + s$neg$n, nE = s$pos$nE + s$neg$nE,
                nC = s$pos$nC + s$neg$nC, xE = if (surv) s$pos$eventsE + s$neg$eventsE else s$pos$xE + s$neg$xE,
                xC = if (surv) s$pos$eventsC + s$neg$eventsC else s$pos$xC + s$neg$xC, n1 = NA, x1 = NA, n2 = NA, x2 = NA),
                row(a))))
            it <- s$inter
            r$observed$setRow(rowKey = "inter", values = blank_na(c(list(n = NA, nE = NA, nC = NA, xE = NA, xC = NA,
                n1 = NA, x1 = NA, n2 = NA, x2 = NA), row(it))))
            note <- if (surv) .("HR experimental / control: Cox model per stratum (Breslow ties); Z and p from the log-rank test (overall: stratified log-rank). Interaction: ratio HR(BM+) / HR(BM−).")
                    else .("Difference experimental − control with unpooled standard errors; overall: inverse-variance combination of the strata. Interaction: difference BM+ − BM−.")
            note <- paste(note, sprintf(.("Gail & Simon test for a qualitative interaction (effects in opposite directions): Q = %s, p = %s."),
                                        format(round(it$Q, 3)), format.pval(it$pGS, digits = 3)))
            r$observed$setNote("m", note, init = FALSE)
            dd <- bm_rand_decide(s, res$ro)
            private$.fillRandDecision(s, dd, res)
            list(s = s, dd = dd)
        },

        .fillRandDecision = function(s, dd, res) {
            o <- self$options
            r <- self$results
            f2 <- function(x) formatC(x, format = "f", digits = 2)
            lv <- paste0(format(round(100 * o$level)), "%")
            ci <- private$.effRange(dd$negLo, dd$negHi)
            r$decision$setRow(rowKey = "s1", values = list(step = .("1. BM+ test"),
                result = sprintf(.("Z = %s %s %s: %s"), f2(s$pos$z), if (dd$step1) ">" else "≤", f2(dd$c1),
                                 if (dd$step1) .("positive") else .("not positive"))))
            if (dd$step1) {
                r$decision$setRow(rowKey = "s2", values = list(step = .("2. BM− effect"),
                    result = sprintf(.("%s CI %s"), lv, ci)))
            } else {
                r$decision$setRow(rowKey = "s2", values = list(step = .("2. Overall test"),
                    result = sprintf(.("Z = %s %s %s: %s"), f2(s$all$z), if (dd$overall) ">" else "≤", f2(dd$c2),
                                     if (dd$overall) .("positive") else .("not positive"))))
            }
            d <- res$d
            info <- if (o$endpoint == "surv")
                sprintf(.("BM+ %d of %s planned events, BM− %d of %s"), s$pos$eventsE + s$pos$eventsC, format(round(d$dPos)),
                        s$neg$eventsE + s$neg$eventsC, format(round(d$dNeg)))
            else sprintf(.("BM+ %d of %s planned patients, BM− %d of %s"), s$pos$n, format(round(d$nPos)), s$neg$n, format(round(d$nNeg)))
            if (isTRUE(o$interim)) {
                cl <- function(z) if (isTRUE(z <= o$zfut)) .("close") else .("continue")
                info <- paste0(info, "; ", sprintf(.("interim futility (Z ≤ %s): BM+ %s, BM− %s"), o$zfut, cl(s$pos$z), cl(s$neg$z)))
            }
            r$decision$setRow(rowKey = "s3", values = list(step = .("Information"), result = info))
            r$decision$setRow(rowKey = "rec", values = list(step = .("Recommendation"), result = private$.recText(dd$rec)))
            r$decision$setNote("m", .("The final rule applied to the data as they are; with less information than planned the error rates differ from the design."), init = FALSE)
        },

        .randRegionState = function(res, obs) {
            o <- self$options
            d <- res$d; ro <- res$ro
            nP <- if (o$sizes == "natural") d$N * o$prev else d$nPos
            inf <- bm_rand_info(d, 0, 0, nP)
            ob <- if (is.null(obs)) NULL else c(thP = obs$s$pos$theta, seP = obs$s$pos$se, thN = obs$s$neg$theta, seN = obs$s$neg$se)
            if (!is.null(ob) && any(!is.finite(ob))) ob <- NULL
            list(kind = "rand", endpoint = o$endpoint, IP = inf$IP, IN = inf$IN, c1 = stats::qnorm(1 - ro$alpha1),
                 c2 = stats::qnorm(1 - ro$alpha2), zg = stats::qnorm((1 + ro$level) / 2), thetaL = ro$thetaL,
                 thetaU = ro$thetaU, obs = ob, ciz = stats::qnorm((1 + o$ciLevel) / 2))
        },

        .randSummary = function(res) {
            o <- self$options
            d <- res$d
            surv <- o$endpoint == "surv"
            lv <- paste0(format(round(100 * o$level)), "%")
            list(design = .("Randomized, stratified by biomarker (Freidlin et al. 2012)"),
                 endpoint = if (surv) .("Time to event: hazard ratio experimental / control (log-rank)")
                            else .("Response: difference in rates, experimental − control"),
                 power = sprintf(.("%s at %s in BM+ (one-sided α = %s)"), o$power, private$.effText(res$ro$thetaT), o$alpha1),
                 sizes = if (o$sizes == "natural") sprintf(.("Natural prevalence of BM+: %s"), o$prev)
                         else sprintf(.("Fixed: %s BM− per BM+ patient"), o$ratioNeg),
                 interim = if (isTRUE(o$interim)) sprintf(.("Futility at %s of the information in each stratum: Z ≤ %s closes the stratum; BM+ closing stops the trial"),
                                                           o$frac, o$zfut) else .("None"),
                 alloc = sprintf(.("%s experimental per control"), o$alloc))
        },

        .randRules = function(res) {
            o <- self$options
            surv <- o$endpoint == "surv"
            lv <- paste0(format(round(100 * o$level)), "%")
            what <- if (surv) .("HR") else .("difference")
            rows <- list(
                list(step = "1", condition = sprintf(.("BM+: one-sided test at α = %s is positive"), o$alpha1), action = .("go to 3")),
                list(step = "2", condition = sprintf(.("Otherwise: overall test (stratified) at α = %s is positive"), o$alpha2),
                     action = paste0(.("All-comers phase III"), "; ", .("otherwise stop"))),
                list(step = "3a", condition = sprintf(.("%s CI of the BM− %s entirely %s %s"), lv, what, if (surv) .("above") else .("below"),
                                                      if (surv) o$hrNoBen else o$diffNoBen), action = .("BM+ only (enrichment) phase III")),
                list(step = "3b", condition = sprintf(.("%s CI of the BM− %s entirely %s %s"), lv, what, if (surv) .("below") else .("above"),
                                                      if (surv) o$hrBen else o$diffBen), action = .("All-comers phase III")),
                list(step = "3c", condition = .("Otherwise"), action = .("Biomarker-stratified phase III")))
            if (isTRUE(o$interim))
                rows <- c(list(list(step = "0", condition = sprintf(.("Interim, at %s of the information: Z ≤ %s in a stratum"), o$frac, o$zfut),
                                    action = .("BM+: stop the trial; BM−: close it (it enters step 2 or 3 with its interim data)"))), rows)
            rows
        },

        # ---- single arm ----------------------------------------------------------------------

        .runSingle = function() {
            o <- self$options
            r <- self$results
            des <- r$searchCache$state
            if (is.null(des) || o$saDesign %in% c("jones", "isotonic")) {
                des <- tryCatch(bm_single_design(o), error = function(e) e)
                if (inherits(des, "error")) return(private$.fail(conditionMessage(des)))
                if (!o$saDesign %in% c("jones", "isotonic")) r$searchCache$setState(des)
            }
            res <- tryCatch(bm_single_compute(o, des), error = function(e) e)
            if (inherits(res, "error")) return(private$.fail(conditionMessage(res)))
            if (o$sizes == "natural" && des$type == "iso")
                private$.notice(.("Natural prevalence is not available for the isotonic rules: the results assume fixed stratum sizes."),
                                "isoNatural", jmvcore::NoticeType$INFO)
            private$.fillSummary(private$.singleSummary(res))
            sz <- bm_single_sizes(des)
            if (des$type == "parallel") {
                for (key in c("neg", "pos")) r$sizes$setRow(rowKey = key, values = list(stage1 = sz$stage1[[key]], total = sz$max[[key]], unsel = "", enrich = ""))
                r$sizes$setRow(rowKey = "total", values = list(stage1 = sum(sz$stage1), total = sum(sz$max), unsel = "", enrich = ""))
            } else {
                for (key in c("neg", "pos")) r$sizes$setRow(rowKey = key, values = list(stage1 = sz$stage1[[key]],
                    unsel = sz$unsel[[key]], enrich = if (key == "pos") sz$enrich[[key]] else sz$stage1[[key]], total = ""))
                r$sizes$setRow(rowKey = "total", values = list(stage1 = sum(sz$stage1), unsel = sum(sz$unsel),
                                                               enrich = sum(sz$enrich), total = sz$max))
            }
            note <- if (res$natural) sprintf(.("Natural prevalence (BM+ %s): stage 1 and the stage 2 of both strata enrol all-comers, so the stratum sizes vary; cut-offs are re-derived at the realised sizes, keeping their planned null tail probability (our extension). The sizes shown are the planned ones."), o$prev)
                    else sprintf(.("Fixed stratum sizes: about %d patients are screened to fill stage 1 (prevalence %s)."),
                                 ceiling(max(sz$stage1[["pos"]] / o$prev, sz$stage1[["neg"]] / (1 - o$prev))), o$prev)
            r$sizes$setNote("n", note, init = FALSE)
            private$.fillRules(private$.singleRules(des))
            if (o$ocTable) {
                for (i in seq_len(nrow(res$sc))) {
                    key <- res$sc$key[i]
                    r$oc$setRow(rowKey = key, values = list(effPos = private$.rateText(res$sc$pp[i]), effNeg = private$.rateText(res$sc$pn[i]),
                        stop = res$oc$stop[i], enrich = res$oc$pos[i], strat = "", neg = res$oc$neg[i], all = res$oc$all[i],
                        expected = res$oc$en[i], pet = res$oc$pet[i]))
                }
                note <- .("Exact binomial probabilities.")
                if (o$endpoint == "surv") note <- paste(note, sprintf(.("Response = alive and event-free at the landmark (%s); rates from the medians assuming exponential survival."), o$lmTime))
                if (des$type == "parallel") note <- paste(note, .("Parallel designs: each stratum at α/2, so the probability of any false go is at most α."))
                r$oc$setNote("m", note, init = FALSE)
            }
            obs <- private$.observeSingle(res)
            img <- r$regionPlot
            if (o$regionPlot && is.null(img$state))
                img$setState(list(kind = "single", g = bm_single_routes(des, o), type = des$type,
                                  n1n = if (des$type == "parallel") des$dn$n1 else des$d$n1n,
                                  n1p = if (des$type == "parallel") des$dp$n1 else des$d$n1p,
                                  obs = if (is.null(obs)) NULL else obs$stage1))
            img <- r$ocMap
            if (o$ocMap && is.null(img$state)) {
                k <- if (res$natural) 11 else 25
                pg <- seq(0.02, 0.98, length.out = k)
                pr <- res$r
                lo <- max(0.01, min(pr) - 0.15); hi <- min(0.99, max(pr) + 0.2)
                pg <- seq(lo, hi, length.out = k)
                g <- expand.grid(pn = pg, pp = pg)
                P <- t(mapply(function(a, b) bm_single_oc1(des, o, pr, a, b)[c("stop", "pos", "neg", "all")], g$pn, g$pp))
                img$setState(list(kind = "single", g = cbind(g, P), sc = res$sc))
            }
        },

        .observeSingle = function(res) {
            o <- self$options
            r <- self$results
            if (is.null(o$biomarker)) return(NULL)
            surv <- o$endpoint == "surv"
            if ((surv && (is.null(o$time) || is.null(o$event) || !nzchar(o$eventLevel %||% ""))) ||
                (!surv && (is.null(o$response) || !nzchar(o$respLevel %||% ""))) || !nzchar(o$bmLevel %||% "")) {
                private$.notice(if (surv) .("Assign the biomarker, the time and the event, and choose their levels (BM+, event). Patients must be in arrival order.")
                                else .("Assign the biomarker and the response, and choose their levels (BM+, response). Patients must be in arrival order."),
                                "dataNeeded", jmvcore::NoticeType$INFO, at = match("observed", self$results$itemNames))
                r$observed$setVisible(FALSE); r$decision$setVisible(FALSE)
                return(NULL)
            }
            dat <- self$data
            pos <- as.character(dat[[o$biomarker]]) == o$bmLevel
            if (surv) {
                tt <- jmvcore::toNumeric(dat[[o$time]]); ev <- as.integer(as.character(dat[[o$event]]) == o$eventLevel)
                y <- ss_landmark_status(tt, ev, o$lmTime)
                nCens <- sum(!is.na(tt) & !is.na(ev) & is.na(y))
            } else {
                y <- as.integer(as.character(dat[[o$response]]) == o$respLevel)
                nCens <- 0
            }
            ok <- !is.na(pos) & !is.na(y)
            pos <- pos[ok]; y <- y[ok]
            des <- res$des
            nat <- res$natural
            if (des$type == "parallel") {
                a <- ss_parallel_observe(des$dn, des$dp, pos, y, nat, res$r[["p0n"]], res$r[["p0p"]])
                st <- list(neg = a$neg, pos = a$pos)
                for (key in c("neg", "pos")) {
                    x <- st[[key]]
                    ci <- ss_cp(x$x1 + x$x2, x$n1 + x$n2, o$ciLevel)
                    r$observed$setRow(rowKey = key, values = blank_na(list(n1 = x$n1, x1 = x$x1, n2 = x$n2, x2 = x$x2,
                        n = x$n1 + x$n2, est = if (x$n1 + x$n2 > 0) (x$x1 + x$x2) / (x$n1 + x$n2) else NA, lower = ci[1], upper = ci[2],
                        nE = NA, nC = NA, xE = NA, xC = NA, z = NA, p = NA)))
                }
                ds <- function(x, lab) sprintf(.("%s: %d/%d responses at stage 1 (continue if ≥ %d); %d/%d in total (promising if ≥ %d): %s"),
                                               lab, x$x1, x$n1, x$k1, x$x1 + x$x2, x$n1 + x$n2, x$k, private$.stageText(x$decision))
                r$decision$setRow(rowKey = "s1", values = list(step = .("BM−"), result = ds(a$neg, .("BM−"))))
                r$decision$setRow(rowKey = "s2", values = list(step = .("BM+"), result = ds(a$pos, .("BM+"))))
                r$decision$setRow(rowKey = "rec", values = list(step = .("Recommendation"),
                    result = if (is.na(a$rec)) .("Not reached yet") else private$.recText(a$rec)))
                stage1 <- c(xn = a$neg$x1, xp = a$pos$x1)
            } else {
                a <- ss_route_observe(des$d, pos, y, nat, res$r[["p0n"]], res$r[["p0p"]], iso = des$type == "iso")
                s1 <- a$stage1; s2 <- a$s2
                for (key in c("neg", "pos")) {
                    n1 <- s1[[paste0("n", substr(key, 1, 1))]]; x1 <- s1[[paste0("x", substr(key, 1, 1))]]
                    n2 <- s2[[paste0("n", substr(key, 1, 1))]]; x2 <- s2[[paste0("x", substr(key, 1, 1))]]
                    ci <- ss_cp(x1 + x2, n1 + n2, o$ciLevel)
                    r$observed$setRow(rowKey = key, values = blank_na(list(n1 = n1, x1 = x1, n2 = n2, x2 = x2, n = n1 + n2,
                        est = if (n1 + n2 > 0) (x1 + x2) / (n1 + n2) else NA, lower = ci[1], upper = ci[2],
                        nE = NA, nC = NA, xE = NA, xC = NA, z = NA, p = NA)))
                }
                route <- switch(if (is.na(a$route)) "none" else a$route,
                    unselected = .("stage 2 in both strata"), enrichment = .("stage 2 in BM+ only (enrichment)"),
                    stop = .("stop"), none = .("stage 1 not complete"))
                r$decision$setRow(rowKey = "s1", values = list(step = .("Stage 1"),
                    result = sprintf(.("BM− %d/%d (cut-off %d), BM+ %d/%d (cut-off %d): %s"), s1[["xn"]], s1[["nn"]], a$k1[["n"]],
                                     s1[["xp"]], s1[["np"]], a$k1[["p"]], route)))
                s2txt <- if (!isTRUE(a$done1) || identical(a$route, "stop")) ""
                    else if (a$route == "unselected")
                        sprintf(.("BM− %d/%d (all-comers if ≥ %d), BM+ %d/%d (BM+ only if ≥ %d)"), s1[["xn"]] + s2[["xn"]],
                                s1[["nn"]] + s2[["nn"]], a$k2[["n"]], s1[["xp"]] + s2[["xp"]], s1[["np"]] + s2[["np"]], a$k2[["p"]])
                    else sprintf(.("BM+ %d/%d (BM+ only if ≥ %d)"), s1[["xp"]] + s2[["xp"]], s1[["np"]] + s2[["np"]], a$k2[["p"]])
                r$decision$setRow(rowKey = "s2", values = list(step = .("Stage 2"), result = s2txt))
                r$decision$setRow(rowKey = "rec", values = list(step = .("Recommendation"), result = switch(a$decision,
                    stage1 = .("Stage 1 in progress"), stage2 = .("Stage 2 in progress"), private$.recText(a$decision))))
                stage1 <- c(xn = s1[["xn"]], xp = s1[["xp"]])
            }
            note <- .("Patients in arrival order; stage 1 is the first patients of each stratum (natural prevalence: the first patients overall). Rates with exact (Clopper-Pearson) intervals that ignore the two-stage rule.")
            if (des$type == "iso") note <- paste(note, .("Isotonic rules: BM− is judged on the smaller and BM+ on the larger of its own and the pooled response rate."))
            if (surv) note <- paste(note, sprintf(.("Response = event-free at the landmark (%s); %d patients censored before it are not evaluable."), o$lmTime, nCens))
            r$observed$setNote("m", note, init = FALSE)
            list(stage1 = stage1)
        },

        .singleSummary = function(res) {
            o <- self$options
            des <- res$des; r <- res$r
            f3 <- function(x) formatC(x, format = "f", digits = 3)
            oc <- function(pn, pp) bm_single_oc1(des, o, r, pn, pp)
            n0 <- oc(r[["p0n"]], r[["p0p"]]); nAll <- oc(r[["p1n"]], r[["p1p"]]); nPos <- oc(r[["p0n"]], r[["p1p"]])
            name <- switch(o$saDesign,
                parashar = paste0(.("Optimal stratified two-stage (Parashar et al. 2016)"), ", ",
                                  if (o$saType == "minimax") .("minimax") else .("optimal")),
                jones = .("User-defined, Jones & Holmgren (2007) rules"),
                isotonic = .("User-defined, isotonic rules (Li & Ivanova 2024)"),
                parallel = paste0(.("Parallel Simon designs"), ", ", if (o$saType == "minimax") .("minimax") else .("optimal")))
            v <- list(design = name,
                      endpoint = if (o$endpoint == "surv") sprintf(.("Event-free at %s (landmark, exponential survival)"), o$lmTime)
                                 else .("Response (binary)"),
                      hyp = sprintf(.("BM−: %s vs %s; BM+: %s vs %s"), f3(r[["p0n"]]), f3(r[["p1n"]]), f3(r[["p0p"]]), f3(r[["p1p"]])),
                      alpha = sprintf(.("%s (target %s)"), f3(1 - n0[["stop"]]), o$saAlpha),
                      powAll = sprintf(.("%s (target %s)"), f3(nAll[["all"]]), o$saPower),
                      powPos = sprintf(.("%s (target %s)"), f3(nPos[["pos"]]), o$saPower),
                      en0 = formatC(n0[["en"]], format = "f", digits = 2),
                      search = if (!is.null(des$limits)) sprintf(.("Up to %d BM− and %d BM+ patients; cut-offs near the expected counts"),
                                                                 des$limits[["neg"]], des$limits[["pos"]]) else "")
            v[private$.summaryKeys()]
        },

        .singleRules = function(des) {
            o <- self$options
            if (des$type == "parallel") {
                a <- des$dn; b <- des$dp
                return(list(
                    list(step = .("BM−"), condition = sprintf(.("Stage 1: %d patients, continue if ≥ %d responses; promising if ≥ %d of %d"), a$n1, a$k1, a$k, a$n), action = ""),
                    list(step = .("BM+"), condition = sprintf(.("Stage 1: %d patients, continue if ≥ %d responses; promising if ≥ %d of %d"), b$n1, b$k1, b$k, b$n), action = ""),
                    list(step = .("Result"), condition = .("Both promising / BM+ only / BM− only / neither"),
                         action = .("All-comers / BM+ only / BM− only / stop"))))
            }
            d <- des$d
            iso <- des$type == "iso"
            pooled <- function(k, n, m, who) if (iso) sprintf(.(" (and the pooled rate × %d ≥ %d)"), n, k) else ""
            list(
                list(step = "1", condition = sprintf(.("Stage 1: %d BM− and %d BM+ patients"), d$n1n, d$n1p), action = ""),
                list(step = "1a", condition = sprintf(.("BM− responses ≥ %d%s"), d$k1n, if (iso) .(" and the pooled count test") else ""),
                     action = sprintf(.("Stage 2 in both strata: up to %d BM− and %d BM+"), d$nn, d$np)),
                list(step = "1b", condition = sprintf(.("Otherwise, BM+ responses ≥ %d%s"), d$k1p, if (iso) .(" or the pooled count test") else ""),
                     action = sprintf(.("Stage 2 in BM+ only: up to %d BM+"), d$nep)),
                list(step = "1c", condition = .("Otherwise"), action = .("Stop")),
                list(step = "2a", condition = sprintf(.("Both strata: BM− responses ≥ %d of %d"), d$kn, d$nn),
                     action = .("Promising in all-comers")),
                list(step = "2b", condition = sprintf(.("Otherwise, BM+ responses ≥ %d of %d"), d$kp, d$np),
                     action = .("Promising in BM+ only")),
                list(step = "2c", condition = sprintf(.("BM+ only: responses ≥ %d of %d"), d$kep, d$nep),
                     action = .("Promising in BM+ only")))
        },

        # ---- shared --------------------------------------------------------------------------

        .fillSummary = function(vals) {
            for (key in names(vals)) if (key %in% private$.summaryKeys())
                self$results$summary$setRow(rowKey = key, values = list(value = vals[[key]]))
        },

        .fillRules = function(rows) {
            tb <- self$results$rules
            for (k in seq_along(rows)) {
                if (k > tb$rowCount) tb$addRow(rowKey = k, values = list())
                tb$setRow(rowKey = k, values = rows[[k]])
            }
        },

        .nRules = function() {
            o <- self$options
            if (o$trial == "randomized") 5 + isTRUE(o$interim) else if (o$saDesign == "parallel") 3 else 7
        },

        .summaryKeys = function() {
            o <- self$options
            if (o$trial == "randomized") c("design", "endpoint", "power", "sizes", "alloc", "interim")
            else c("design", "endpoint", "hyp", "alpha", "powAll", "powPos", "en0", if (o$saDesign == "parashar") "search")
        },

        .summaryLabel = function(key) {
            switch(key, design = .("Design"), endpoint = .("Endpoint"), power = .("Power"), sizes = .("Stratum sizes"),
                   alloc = .("Allocation"), interim = .("Interim analysis"), hyp = .("Uninteresting vs target rates"),
                   alpha = .("Type I error: any go under H0 (attained)"), powAll = .("Power: all-comers when both strata respond"),
                   powPos = .("Power: BM+ only when only BM+ responds"), en0 = .("Expected N under H0"),
                   search = .("Search"))
        },

        .scenarioKeys = function() {
            o <- self$options
            k <- if (o$trial == "randomized") c("null", "pos", "equal", "larger") else c("null", "pos", "all")
            if (isTRUE(o$userScen)) c(k, "user") else k
        },

        .scenarioLabel = function(key) {
            switch(key, null = .("No effect"), pos = .("BM+ only"), equal = .("Equal (intermediate)"),
                   larger = .("BM+ larger"), all = .("Both strata"), user = .("User"))
        },

        .recText = function(rec) {
            switch(rec, stop = .("Stop: no further testing"), enrich = .("Phase III in BM+ only (enrichment)"),
                   pos = .("Promising in BM+ only: phase III in BM+ (enrichment)"),
                   strat = .("Biomarker-stratified phase III"), all = .("Phase III in all-comers"),
                   neg = .("Promising in BM− only (unexpected)"), rec)
        },

        .stageText = function(dec) {
            switch(dec, stage1 = .("stage 1 in progress"), stage2 = .("stage 2 in progress"),
                   stop = .("not promising"), promising = .("promising"), dec)
        },

        # effect display: HR = exp(-theta) for survival, the difference otherwise
        .effPoint = function(theta) if (self$options$endpoint == "surv") exp(-theta) else theta,

        .effText = function(theta) {
            if (self$options$endpoint == "surv") sprintf("HR %s", formatC(exp(-theta), format = "f", digits = 2))
            else sprintf("%+.2f", theta)
        },

        .effRange = function(lo, hi) {
            if (self$options$endpoint == "surv")
                sprintf("%s – %s", formatC(exp(-hi), format = "f", digits = 2), formatC(exp(-lo), format = "f", digits = 2))
            else sprintf("%s – %s", formatC(lo, format = "f", digits = 3), formatC(hi, format = "f", digits = 3))
        },

        .rateText = function(p) formatC(p, format = "f", digits = 3),

        .fail = function(text) {
            private$.notice(text, "inputProblem", jmvcore::NoticeType$ERROR)
            for (nm in c("summary", "sizes", "rules", "oc", "observed", "decision", "regionPlot", "ocMap"))
                self$results[[nm]]$setVisible(FALSE)
            invisible(NULL)
        },

        # Notice, not `type: Notice` in the yaml, which does not compile (jamovi skill)
        .notice = function(text, name, type, at = 1) {
            n <- jmvcore::Notice$new(options = self$options, name = name, type = type)
            n$setContent(text)
            self$results$insert(at, n)
        },

        .plotRegions = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- private$.recLabels()
            if (st$kind == "rand") {
                surv <- st$endpoint == "surv"
                labels$x <- if (surv) .("Estimated HR in BM+") else .("Estimated difference in BM+")
                labels$y <- if (surv) .("Estimated HR in BM−") else .("Estimated difference in BM−")
                labels$obs <- .("Observed (with interval)")
                print(plot_bm_rand_regions(st, list(theme_only(ggtheme)), labels))
            } else {
                labels$x <- sprintf(.("BM− responses at stage 1 (of %d)"), st$n1n)
                labels$y <- sprintf(.("BM+ responses at stage 1 (of %d)"), st$n1p)
                labels$both <- .("Stage 2, both strata"); labels$posRoute <- .("Stage 2, BM+ only")
                labels$negRoute <- .("Stage 2, BM− only"); labels$stopRoute <- .("Stop")
                labels$goAll <- .("All-comers already reached"); labels$goPos <- .("BM+ already reached")
                labels$now <- .("Observed at stage 1")
                print(plot_bm_routes(st, list(theme_only(ggtheme)), labels))
            }
            TRUE
        },

        .plotOCMap = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- private$.recLabels()
            if (st$kind == "rand") {
                surv <- st$endpoint == "surv"
                labels$x <- if (surv) .("True HR in BM+") else .("True difference in BM+")
                labels$y <- if (surv) .("True HR in BM−") else .("True difference in BM−")
            } else {
                labels$x <- .("True rate in BM−"); labels$y <- .("True rate in BM+")
            }
            labels$prob <- .("Probability")
            labels$sc <- c(null = .("No effect"), pos = .("BM+ only"), equal = .("Equal"), larger = .("BM+ larger"),
                           all = .("Both"), user = .("User"))
            print(plot_bm_ocmap(st, list(theme_only(ggtheme)), labels))
            TRUE
        },

        .recLabels = function() list(stop = .("Stop"), enrich = .("BM+ only"), strat = .("Stratified"),
                                      all = .("All-comers"), neg = .("BM− only")))
)
