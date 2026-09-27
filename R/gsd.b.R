
gsdClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "gsdClass",
    inherit = gsdBase,
    private = list(

        # Rows exist by the end of .init (jamovi restores results into them);
        # .run only writes values. Everything is recomputed per run: the design
        # takes milliseconds, and the plots keep their state until cleared.
        .init = function() {
            o <- self$options
            r <- self$results
            K <- tryCatch(length(gsd_timing(o$timing, o$nLooks)), error = function(e) o$nLooks)
            surv <- o$endpoint == "surv"

            for (key in private$.summaryKeys())
                r$summary$addRow(rowKey = key, values = list(name = private$.summaryLabel(key)))

            b <- r$bounds
            for (k in seq_len(K)) b$addRow(rowKey = k, values = list(look = k))
            b$getColumn("n")$setTitle(if (surv) .("Events") else .("N"))
            effTitle <- switch(o$endpoint, surv = .("HR"), generic = .("Effect"), .("Difference"))
            b$getColumn("effEst")$setTitle(effTitle)
            b$getColumn("futEst")$setTitle(effTitle)
            fut <- o$efficacy %in% c("pt", "triangular") || o$futility != "none"
            b$getColumn("futZ")$setVisible(fut)
            b$getColumn("futEst")$setVisible(fut)

            oc <- r$oc
            for (g in c("0", "0.5", "1", "1.5"))
                oc$addRow(rowKey = g, values = list(hyp = private$.ocLabel(g)))
            oc$getColumn("effect")$setTitle(effTitle)
            oc$getColumn("en")$setTitle(if (surv) .("Expected events") else .("Expected N"))
            oc$getColumn("earlyFut")$setVisible(fut)

            if (surv)
                for (k in seq_len(K)) r$timing$addRow(rowKey = k, values = list(look = k))

            # interim data: one row per analysis (.init sees only the header, so read the rows)
            if (!is.null(o$lookN)) {
                dat <- if (nrow(self$data) > 0) self$data else self$readDataset()
                nr <- sum(!is.na(jmvcore::toNumeric(dat[[o$lookN]])))
                m <- r$monitor
                for (k in seq_len(nr)) m$addRow(rowKey = k, values = list(look = k))
                m$getColumn("n")$setTitle(if (surv) .("Events") else .("N"))
                m$getColumn("est")$setTitle(effTitle)
                m$getColumn("fut")$setVisible(fut)
                r$final$addRow(rowKey = 1, values = list())
            }
        },

        .run = function() {
            o <- self$options
            res <- tryCatch(gsd_compute(private$.optionList()), error = function(e) e)
            # an incomplete or inconsistent input is a message for the user, not an
            # analysis error with a stack trace: say what to fix and show nothing else
            if (inherits(res, "error")) {
                private$.notice(conditionMessage(res))
                for (nm in c("summary", "bounds", "oc", "timing", "boundPlot", "spendPlot", "powerPlot"))
                    self$results[[nm]]$setVisible(FALSE)
                return()
            }
            d <- res$d
            K <- d$K
            surv <- o$endpoint == "surv"

            # summary
            vals <- private$.summaryValues(res)
            for (key in names(vals))
                self$results$summary$setRow(rowKey = key, values = list(value = vals[[key]]))

            # boundaries
            b <- self$results$bounds
            lk <- res$looks
            for (k in seq_len(K)) {
                row <- as.list(lk[k, ])
                row$n <- ceiling(row$n - 1e-9)
                b$setRow(rowKey = k, values = blank_na(row))
            }
            private$.boundNotes(res)

            # operating characteristics
            if (o$ocTable) {
                for (i in seq_len(nrow(res$oc))) {
                    x <- res$oc[i, ]
                    self$results$oc$setRow(rowKey = as.character(x$g),
                        values = blank_na(list(effect = x$effect, reject = x$reject, earlyEff = x$earlyEff,
                                               earlyFut = x$earlyFut, en = x$en)))
                }
            }

            # survival timing
            if (surv) {
                for (k in seq_len(K)) {
                    x <- res$timing[k, ]
                    self$results$timing$setRow(rowKey = k, values = blank_na(list(
                        events = ceiling(x$events - 1e-9), timeH1 = x$timeH1, timeH0 = x$timeH0,
                        enrolled = ceiling(x$enrolled - 1e-9))))
                }
            }

            # plot data (state is set only when absent, so unrelated clicks do not redraw)
            img <- self$results$boundPlot
            if (is.null(img$state))
                img$setState(list(looks = lk, sided = res$sided, endpoint = o$endpoint))
            img <- self$results$spendPlot
            if (is.null(img$state)) {
                cumBeta <- if (res$hasFut) c(d$cumBetaEarly[-K], d$spec$beta) else NULL
                img$setState(list(t = d$t, cumAlpha = d$cumAlpha, alpha = d$spec$alpha,
                                  cumBeta = cumBeta, beta = d$spec$beta))
            }
            if (!is.null(o$lookN)) {
                mon <- tryCatch(private$.monitor(res), error = function(e) e)
                if (inherits(mon, "error")) {
                    private$.notice(conditionMessage(mon), "monitorProblem", jmvcore::NoticeType$WARNING)
                    for (nm in c("monitor", "final", "monitorPlot")) self$results[[nm]]$setVisible(FALSE)
                }
            }

            img <- self$results$powerPlot
            if (is.null(img$state) && o$powerPlot) {
                g <- seq(if (res$sided == 2) -1.5 else -0.5, 1.5, by = 0.05)
                pr <- lapply(g, function(x) gs_probs(d$t, d$up, d$lo, d$wedge, theta = x * d$eta))
                img$setState(list(
                    g = g,
                    reject = vapply(pr, function(p) sum(p$up) + (if (res$sided == 2) sum(p$lo) else 0), 0),
                    en = vapply(pr, function(p) sum(d$t * rowSums(p)) * res$size$Nmax, 0)))
            }
        },

        # Interim data: boundaries at the observed information, decisions, repeated CIs,
        # conditional and predictive power, and stagewise inference once the trial stops.
        .monitor = function(res) {
            o <- self$options
            dat <- self$data
            num <- function(v) if (is.null(v)) NULL else jmvcore::toNumeric(dat[[v]])
            dir <- if (res$theta1 < 0) -1 else 1
            looks <- monitor_looks(n = num(o$lookN), est = num(o$estimate), se = num(o$se), z = num(o$zStat),
                                   frac = num(o$infoFrac), Nmax = res$size$Nmax, Imax = res$size$Imax, dir = dir)
            m <- gs_monitor(res$d, looks, final = o$finalLook, dir = dir)
            nat <- if (o$endpoint == "surv") exp else identity
            tb <- self$results$monitor
            lk <- m$looks
            for (k in seq_len(nrow(lk))) {
                if (k > tb$rowCount) tb$addRow(rowKey = k, values = list(look = k))
                x <- lk[k, ]
                tb$setRow(rowKey = k, values = blank_na(list(
                    n = looks$n[k], t = x$t, est = nat(looks$est[k]), z = x$z, bval = x$b, up = x$up, fut = x$fut,
                    decision = private$.decisionLabel(x$decision, k == nrow(lk) && m$final),
                    rciLo = nat(x$rciLo), rciHi = nat(x$rciHi), cpH1 = x$cpH1, cpTrend = x$cpTrend, pp = x$pp)))
            }
            tb$setNote("scale", sprintf(.("Boundaries recomputed at the observed information. Z is oriented so that positive values favour H1. Repeated CIs at %s%%; conditional power under the design effect (H1) and the current trend; predictive power with a flat prior."),
                                        format(round(100 * m$level, 1))), init = FALSE)
            ft <- self$results$final
            if (is.null(m$inference)) {
                ft$setVisible(FALSE)
            } else {
                fi <- m$inference
                lo <- nat(fi$lower); hi <- nat(fi$upper)
                ft$setRow(rowKey = 1, values = blank_na(list(
                    look = fi$stage, naive = nat(fi$naive), naiveLo = nat(fi$naiveLower), naiveHi = nat(fi$naiveUpper),
                    mue = nat(fi$mue), lower = min(lo, hi), upper = max(lo, hi), p = fi$p.value)))
                ft$setNote("ord", sprintf(.("Stagewise ordering (Tsiatis, Rosner & Mehta 1984): %s p-value, median-unbiased estimate and %s%% interval that account for the stopping rule."),
                                          if (res$sided == 2) .("two-sided") else .("one-sided"),
                                          format(round(100 * m$level, 1))), init = FALSE)
            }
            img <- self$results$monitorPlot
            if (is.null(img$state)) {
                b <- m$bounds
                futb <- if (res$sided == 2) b$wedge else b$lo
                # the final look closes the continuation region; it is not a futility boundary
                futb[length(futb)] <- NA
                img$setState(list(
                    looks = data.frame(n = m$tt * res$size$Nmax, t = m$tt, effZ = b$up, effP = NA, effEst = NA,
                                       futZ = ifelse(is.finite(futb) & (res$sided == 1 | futb > 1e-5), futb, NA), futEst = NA),
                    obs = data.frame(n = looks$n, t = looks$t, z = looks$z),
                    sided = res$sided, endpoint = o$endpoint))
            }
            m
        },

        .decisionLabel = function(dec, last) {
            switch(dec,
                efficacy = if (last) .("Final: reject H0") else .("Stop: reject H0"),
                accept = .("Final: H0 not rejected"),
                futility = .("Stop: futility"),
                futilityNB = .("Futility boundary crossed (non-binding)"),
                continue = .("Continue"),
                after = .("(after stopping)"))
        },

        # Notice, not `type: Notice` in the yaml, which does not compile (jamovi skill)
        .notice = function(text, name = "inputProblem", type = jmvcore::NoticeType$ERROR) {
            n <- jmvcore::Notice$new(options = self$options, name = name, type = type)
            n$setContent(text)
            self$results$insert(1, n)
        },

        .optionList = function() {
            nm <- c("endpoint", "sided", "alpha", "power", "nLooks", "timing", "ratio", "meanDiff", "sd",
                    "p0", "p1", "genEffect", "genSd", "hr", "ctrlSpec", "ctrlMedian", "ctrlHazard",
                    "survTime", "survProb", "eventsMethod", "solveFor", "accrualTime", "accrualRate",
                    "followUp", "nTotal", "dropoutProb", "dropoutTime", "efficacy", "effRho", "effGamma",
                    "effNu", "effUser", "wtDelta", "hpZ", "ptDelta1", "ptDelta0", "futility", "futSf",
                    "futRho", "futGamma", "futNu", "futZ", "futCp", "futCpUnder", "binding")
            stats::setNames(lapply(nm, function(n) self$options[[n]]), nm)
        },

        .summaryKeys = function() {
            o <- self$options
            keys <- c("test", "looks", "efficacy", "futility", "alpha", "power", "inflation")
            if (o$futility != "none" && !o$binding && !o$efficacy %in% c("pt", "triangular"))
                keys <- c(keys, "alphaObeyed")
            if (o$endpoint == "surv")
                c(keys, "events", "eventsFix", "eventsH0", "eventsH1", "n", "accrual", "followup",
                  "duration", "durationH1", "durationH0")
            else c(keys, "n", "nGroups", "nFix", "nH0", "nH1")
        },

        .summaryLabel = function(key) {
            switch(key,
                test = .("Test"), looks = .("Analyses"), efficacy = .("Efficacy boundary"),
                futility = .("Futility boundary"), alpha = .("Type I error"),
                alphaObeyed = .("Type I error if futility stops are obeyed"),
                power = .("Power"), inflation = .("Inflation factor"),
                events = .("Maximum events"), eventsFix = .("Fixed-design events"),
                eventsH0 = .("Expected events under H0"), eventsH1 = .("Expected events under H1"),
                n = .("Maximum sample size"), accrual = .("Accrual duration"),
                followup = .("Minimum follow-up"), duration = .("Maximum study duration"),
                durationH1 = .("Expected study duration under H1"),
                durationH0 = .("Expected study duration under H0"),
                nGroups = .("Maximum per group (treatment / control)"),
                nFix = .("Fixed-design sample size"), nH0 = .("Expected sample size under H0"),
                nH1 = .("Expected sample size under H1"))
        },

        .summaryValues = function(res) {
            o <- self$options
            d <- res$d
            f3 <- function(x) formatC(x, format = "f", digits = 3)
            f1 <- function(x) formatC(x, format = "f", digits = 1)
            up <- function(x) format(ceiling(x - 1e-9))
            effLab <- private$.effLabel()
            futLab <- if (o$efficacy == "pt") .("Pampallona-Tsiatis") else if (o$efficacy == "triangular")
                .("Triangular (lower boundary)") else private$.futLabel()
            binding <- o$efficacy == "triangular" || ((o$futility != "none" || o$efficacy == "pt") && o$binding)
            if (res$hasFut) futLab <- paste0(futLab, ", ", if (binding) .("binding") else .("non-binding"))
            v <- list(
                test = paste0(if (res$sided == 2) .("Two-sided") else .("One-sided"), ", α = ", format(o$alpha)),
                looks = paste0(d$K, " (t = ", paste(formatC(d$t, format = "f", digits = 3), collapse = ", "), ")"),
                efficacy = effLab, futility = futLab,
                alpha = f3(if (binding) d$alphaBinding else d$alpha),
                alphaObeyed = f3(d$alphaBinding),
                power = f3(d$power), inflation = f3(d$inflation))
            s <- res$size
            if (o$endpoint == "surv") {
                pl <- s$plan
                v <- c(v, list(events = up(s$Nmax), eventsFix = up(s$Nfix),
                               eventsH0 = f1(d$EtH0 * s$Nmax), eventsH1 = f1(d$EtH1 * s$Nmax),
                               n = up(pl$N), accrual = f3(pl$A), followup = f3(pl$F),
                               duration = f3(res$duration[["max"]]), durationH1 = f3(res$duration[["H1"]]),
                               durationH0 = f3(res$duration[["H0"]])))
            } else {
                r <- o$ratio
                v <- c(v, list(n = up(s$Nmax),
                               nGroups = paste(up(s$Nmax * r / (1 + r)), "/", up(s$Nmax / (1 + r))),
                               nFix = up(s$Nfix), nH0 = f1(d$EtH0 * s$Nmax), nH1 = f1(d$EtH1 * s$Nmax)))
            }
            v[private$.summaryKeys()]
        },

        .effLabel = function() {
            o <- self$options
            switch(o$efficacy,
                ldof = .("Lan-DeMets O'Brien-Fleming spending"),
                ldpocock = .("Lan-DeMets Pocock spending"),
                kd = paste0(.("Kim-DeMets spending"), ", ρ = ", o$effRho),
                hsd = paste0(.("Hwang-Shih-DeCani spending"), ", γ = ", o$effGamma),
                exponential = paste0(.("Exponential spending"), ", ν = ", o$effNu),
                user = paste0(.("User-defined spending"), " (", o$effUser, ")"),
                obf = .("O'Brien-Fleming"), pocock = .("Pocock"),
                wt = paste0(.("Wang-Tsiatis"), ", Δ = ", o$wtDelta),
                hp = paste0(.("Haybittle-Peto"), ", z = ", o$hpZ),
                pt = paste0(.("Pampallona-Tsiatis"), ", Δ = ", o$ptDelta1, " / ", o$ptDelta0),
                triangular = .("Triangular test (Whitehead)"))
        },

        .futLabel = function() {
            o <- self$options
            switch(o$futility,
                none = .("None"),
                spending = paste0(.("Beta spending"), ": ", switch(o$futSf,
                    ldof = .("Lan-DeMets O'Brien-Fleming"), ldpocock = .("Lan-DeMets Pocock"),
                    kd = paste0(.("Kim-DeMets"), ", ρ = ", o$futRho),
                    hsd = paste0(.("Hwang-Shih-DeCani"), ", γ = ", o$futGamma),
                    exponential = paste0(.("Exponential"), ", ν = ", o$futNu))),
                fixedz = paste0(.("Fixed z-values"), " (", o$futZ, ")"),
                cp = paste0(.("Conditional power"), " < ", o$futCp, " (",
                            if (o$futCpUnder == "trend") .("current trend") else .("design effect"), ")"))
        },

        .ocLabel = function(g) {
            switch(g, "0" = .("No effect (H0)"), "0.5" = .("Half the design effect"),
                   "1" = .("Design effect (H1)"), "1.5" = .("1.5 × design effect"))
        },

        .boundNotes = function(res) {
            o <- self$options
            b <- self$results$bounds
            note <- .("Information: fraction of the maximum. Effect: the estimate at which Z reaches the boundary.")
            if (res$sided == 2)
                note <- paste(note, .("Two-sided: efficacy if |Z| reaches the boundary; futility if |Z| is below the futility value."))
            b$setNote("scale", note, init = FALSE)
            if (o$efficacy == "triangular")
                b$setNote("tri", .("Whitehead's triangular test with the Christmas-tree correction; the attained α and power are exact."), init = FALSE)
        },

        .plotBounds = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            surv <- st$endpoint == "surv"
            scale <- self$options$boundScale
            xvar <- self$options$boundX
            labels <- list(eff = .("Efficacy"), fut = .("Futility"),
                           x = if (xvar == "t") .("Information fraction") else if (surv) .("Events") else .("Sample size"),
                           y = switch(scale, z = .("Z boundary"), b = .("B-value (Z·√t)"), p = .("Nominal p-value"),
                                      effect = if (surv) .("Hazard ratio") else .("Effect")))
            print(plot_bounds(st, scale, theme$color, list(theme_only(ggtheme)), labels, xvar))
            TRUE
        },

        .plotMonitor = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            surv <- st$endpoint == "surv"
            xvar <- self$options$boundX
            labels <- list(eff = .("Efficacy"), fut = .("Futility"), obs = .("Observed"),
                           x = if (xvar == "t") .("Information fraction") else if (surv) .("Events") else .("Sample size"),
                           y = .("Z"))
            print(plot_monitor(st, theme$color, list(theme_only(ggtheme)), labels, xvar))
            TRUE
        },

        .plotSpend = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(alpha = .("Type I error (α)"), beta = .("Type II error (β)"),
                           x = .("Information fraction"), y = .("Proportion of error spent"))
            print(plot_spend(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        },

        .plotPower = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            surv <- self$options$endpoint == "surv"
            labels <- list(power = .("Probability of rejecting H0"),
                           en = if (surv) .("Expected events") else .("Expected sample size"),
                           x = .("True effect / design effect (log hazard ratio scale for survival)"))
            if (!surv) labels$x <- .("True effect / design effect")
            print(plot_power(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        })
)
