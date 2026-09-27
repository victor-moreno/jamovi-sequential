
bayesmonClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "bayesmonClass",
    inherit = bayesmonBase,
    private = list(

        .init = function() {
            o <- self$options
            r <- self$results
            for (key in private$.summaryKeys())
                r$summary$addRow(rowKey = key, values = list(name = private$.summaryLabel(key)))
            if (o$model != "bf") {
                K <- tryCatch(length(gsd_timing(o$timing, o$nLooks)), error = function(e) o$nLooks)
                for (k in seq_len(K)) r$looks$addRow(rowKey = k, values = list(look = k))
                labs <- if (o$model == "normal") c(.("No effect"), .("Half the effect"), .("H1"), .("1.5 × H1"))
                        else c(.("No difference"), .("Midway"), .("H1"), .("Beyond H1"))
                for (i in 1:4) r$oc$addRow(rowKey = i, values = list(label = labs[i]))
                bin <- o$model == "binary"
                for (cn in c("eff", "fut")) r$looks$getColumn(cn)$setVisible(!bin)
                r$looks$getColumn("n")$setTitle(if (bin) .("N per arm") else .("N"))
                r$oc$getColumn("theta")$setTitle(if (bin) .("Treatment rate") else .("True effect"))
            } else {
                r$bfDesign$addRow(rowKey = 1L, values = list(label = .("No effect")))
                r$bfDesign$addRow(rowKey = 2L, values = list(label = .("Design effect")))
            }
            hasN <- o$model == "normal" && !is.null(o$lookN) && !is.null(o$estimate) && !is.null(o$se)
            hasB <- o$model == "binary" && !is.null(o$xC) && !is.null(o$nC) && !is.null(o$xT) && !is.null(o$nT)
            hasF <- o$model == "bf" && !is.null(o$dep)
            r$monitor$setVisible(hasN); r$monitorBin$setVisible(hasB); r$bfTest$setVisible(hasF)
            if (!(hasN || hasB)) r$probPlot$setVisible(FALSE)
            if (!hasF) r$bfPlot$setVisible(FALSE)
            if (hasN || hasB) {
                dat <- if (nrow(self$data) > 0) self$data else self$readDataset()
                vars <- if (hasN) c(o$lookN, o$estimate, o$se) else c(o$xC, o$nC, o$xT, o$nT)
                n <- sum(stats::complete.cases(dat[vars]))
                tb <- if (hasN) r$monitor else r$monitorBin
                for (i in seq_len(n)) tb$addRow(rowKey = i, values = list(look = i))
            }
        },

        .run = function() {
            o <- self$options
            r <- self$results
            # the design (and the simulation for the Bayes factor) is cached on the summary
            res <- r$summary$state
            if (is.null(res)) {
                res <- tryCatch(bayes_compute(private$.optionList()), error = function(e) e)
                if (inherits(res, "error")) return(private$.fail(conditionMessage(res)))
                r$summary$setState(res)
            }
            vals <- private$.summaryValues(res)
            for (key in names(vals)) r$summary$setRow(rowKey = key, values = list(value = vals[[key]]))

            if (o$model == "bf") return(private$.runBF(res))

            f3 <- function(x) formatC(x, format = "f", digits = 3)
            dir <- res$dir
            K <- length(res$t)
            for (k in seq_len(K)) {
                v <- list(look = k, n = if (o$model == "normal") res$des$nk[k] else res$nk[k],
                          stop0 = res$ocLooks[[1]]$pEff[k] + res$ocLooks[[1]]$pFut[k],
                          stop1 = res$ocLooks[[3]]$pEff[k] + res$ocLooks[[3]]$pFut[k])
                if (o$model == "normal") {
                    e <- res$des$eff[k]; f <- res$des$fut[k]
                    v$eff <- if (!is.finite(e)) "—" else paste(if (dir > 0) "≥" else "≤", f3(dir * e))
                    v$fut <- if (k == K || !is.finite(f)) "—" else paste(if (dir > 0) "≤" else "≥", f3(dir * f))
                }
                r$looks$setRow(rowKey = k, values = v)
            }
            if (o$model == "binary")
                r$looks$setNote("bin", .("The thresholds depend on the responses in both arms; the probabilities are exact."), init = FALSE)
            if (o$ocTable) {
                for (i in seq_len(nrow(res$oc)))
                    r$oc$setRow(rowKey = i, values = as.list(res$oc[i, c("theta", "reject", "futStop", "en")]))
                r$oc$setNote("exact", if (o$model == "normal")
                    .("Frequentist probabilities, exact (the rules are thresholds on the estimate at each look).")
                    else .("Frequentist probabilities, exact (dynamic programming over the responses in both arms); expected N counts both arms."),
                    init = FALSE)
            }

            mon <- NULL
            if (r$monitor$visible || r$monitorBin$visible) {
                mon <- tryCatch(private$.monitor(res), error = function(e) e)
                if (inherits(mon, "error")) {
                    private$.notice(conditionMessage(mon), "monitorProblem", jmvcore::NoticeType$WARNING)
                    for (nm in c("monitor", "monitorBin", "probPlot")) r[[nm]]$setVisible(FALSE)
                    mon <- NULL
                }
            }

            img <- r$boundPlot
            if (o$model == "normal" && o$boundPlot && is.null(img$state)) {
                d <- res$des
                df <- data.frame(x = d$nk, y = dir * d$eff, series = .("Success"), type = "line")
                if (!is.na(d$gammaF) && K > 1)
                    df <- rbind(df, data.frame(x = d$nk[-K], y = dir * d$fut[-K], series = .("Futility"), type = "line"))
                df <- df[is.finite(df$y), ]
                if (!is.null(mon)) df <- rbind(df, data.frame(x = mon$n, y = mon$est, series = .("Observed"), type = "point"))
                img$setState(list(df = df, points = TRUE, h = data.frame(y = 0), inkSeries = .("Observed")))
            }
            img <- r$probPlot
            if (!is.null(mon) && o$probPlot && is.null(img$state)) {
                x <- if (o$model == "normal") mon$n else (mon$nC + mon$nT) / 2
                df <- data.frame(x = x, y = mon$pDelta, series = .("P(beyond δ)"), type = "line")
                if (any(!is.na(mon$pp))) df <- rbind(df, data.frame(x = x[!is.na(mon$pp)], y = mon$pp[!is.na(mon$pp)],
                                                                  series = .("Predictive probability"), type = "line"))
                h <- data.frame(y = c(o$gammaE, if (isTRUE(o$futility)) o$gammaF))
                img$setState(list(df = df, points = TRUE, h = h, ylim = c(0, 1)))
            }
            img <- r$densityPlot
            if (o$densityPlot && is.null(img$state)) img$setState(private$.densityState(res, mon))
            img <- r$ocPlot
            if (o$ocPlot && is.null(img$state)) img$setState(private$.ocCurve(res))
        },

        .monitor = function(res) {
            o <- self$options
            dat <- self$data
            num <- function(v) jmvcore::toNumeric(dat[[v]])
            if (o$model == "normal") {
                m <- bayes_monitor_normal(res, o, num(o$lookN), num(o$estimate), num(o$se), o$level)
                tb <- self$results$monitor
                for (i in seq_len(nrow(m))) {
                    if (i > tb$rowCount) tb$addRow(rowKey = i, values = list(look = i))
                    tb$setRow(rowKey = i, values = blank_na(list(n = m$n[i], est = m$est[i], postMean = m$postMean[i],
                        lower = m$lower[i], upper = m$upper[i], pBenefit = m$pBenefit[i], pDelta = m$pDelta[i],
                        pp = m$pp[i], decision = private$.decisionLabel(m$decision[i]))))
                }
            } else {
                m <- bayes_monitor_binary(res, o, num(o$xC), num(o$nC), num(o$xT), num(o$nT), o$level)
                tb <- self$results$monitorBin
                for (i in seq_len(nrow(m))) {
                    if (i > tb$rowCount) tb$addRow(rowKey = i, values = list(look = i))
                    tb$setRow(rowKey = i, values = blank_na(list(xC = m$xC[i], nC = m$nC[i], xT = m$xT[i], nT = m$nT[i],
                        postMean = m$postMean[i], lower = m$lower[i], upper = m$upper[i], pBenefit = m$pBenefit[i],
                        pDelta = m$pDelta[i], pp = m$pp[i], decision = private$.decisionLabel(m$decision[i]))))
                }
            }
            tb$setNote("pp", sprintf(.("Posterior intervals at level %s. The predictive probability is that of success at the planned maximum."),
                                     o$level), init = FALSE)
            m
        },

        .runBF = function(res) {
            o <- self$options
            r <- self$results
            if (o$ocTable) {
                for (i in 1:2) {
                    s <- res$sims[[i]]
                    r$bfDesign$setRow(rowKey = i, values = list(d = c(0, o$d1)[i], pH1 = s[["pH1"]], pH0 = s[["pH0"]],
                        pInc = s[["pInc"]], asn = s[["asn"]], median = s[["median"]]))
                }
                note <- sprintf(.("Simulated with %d runs per row (seed %d). N counts all observations."), o$nsim, o$seed)
                if (res$maxBF01 < o$bfK)
                    note <- paste(note, sprintf(.("With at most %d observations (per group) BF₀₁ cannot exceed %s, so evidence for H0 at k = %s is impossible."),
                                                o$bfNMax, formatC(res$maxBF01, format = "f", digits = 2), o$bfK))
                r$bfDesign$setNote("sim", note, init = FALSE)
            }
            if (!r$bfTest$visible) return(invisible(NULL))
            a <- tryCatch(private$.analyseBF(), error = function(e) e)
            if (inherits(a, "error")) {
                private$.notice(conditionMessage(a), "dataProblem", jmvcore::NoticeType$WARNING)
                for (nm in c("bfTest", "bfPlot")) r[[nm]]$setVisible(FALSE)
                return(invisible(NULL))
            }
            img <- r$bfPlot
            if (o$probPlot && is.null(img$state)) {
                k <- seq_len(a$n)
                keep <- a$eligible[k]
                df <- data.frame(x = k[keep], y = exp(a$logbf[k][keep]), series = .("BF₁₀"), type = "line")
                st <- list(df = df, points = length(k) <= 200, logy = TRUE, h = data.frame(y = c(o$bfK, 1 / o$bfK, 1)))
                if (a$decision %in% c("H1", "H0")) st$mark <- c(x = a$n, y = exp(a$logbf[a$n]))
                img$setState(st)
            }
        },

        .analyseBF = function() {
            o <- self$options
            dat <- self$data
            x <- jmvcore::toNumeric(dat[[o$dep]])
            g <- NULL
            if (o$bfSamples == "two") {
                if (is.null(o$group)) stop(.("A two-sample test needs a grouping variable with two levels."))
                gf <- droplevels(as.factor(dat[[o$group]]))
                if (nlevels(gf) != 2) stop(.("The grouping variable must have exactly two levels."))
                g <- as.integer(gf)
                keep <- !is.na(x) & !is.na(g)
                g <- g[keep]
                levs <- levels(gf)
            } else keep <- !is.na(x)
            x <- x[keep]
            nMax <- o$bfNMax * (if (is.null(g)) 1 else 2)
            a <- bf_analyse(x, g, o$rScale, o$bfK, o$bfNMin, nMax, o$mu0)
            tb <- self$results$bfTest
            shown <- which(a$eligible[seq_len(a$n)])
            if (length(shown) > 20) shown <- unique(c(shown[round(seq(1, length(shown), length.out = 20))], a$n))
            for (i in shown) {
                dec <- if (i == a$n) switch(a$decision, H1 = .("Evidence for H1: stop"), H0 = .("Evidence for H0: stop"),
                                            inconclusive = .("Maximum reached: inconclusive"), .("Continue")) else ""
                tb$addRow(rowKey = i, values = blank_na(list(n = i, bf = exp(a$logbf[i]), d = a$d[i], decision = dec)))
            }
            note <- sprintf(.("JZS Bayes factor, Cauchy prior scale r = %s, two-sided; tests from %d observations%s."),
                            o$rScale, o$bfNMin, if (is.null(g)) "" else .(" per group"))
            if (!is.null(g)) note <- paste(note, sprintf(.("Effect d = (%s − %s) / pooled SD."), levs[2], levs[1]))
            tb$setNote("bf", note, init = FALSE)
            a
        },

        .densityState = function(res, mon) {
            o <- self$options
            if (o$model == "normal") {
                pr <- res$prior; dir <- res$dir
                last <- if (!is.null(mon)) nrow(mon) else 0
                sdL <- if (last) (mon$upper[last] - mon$lower[last]) / (2 * stats::qnorm(1 - (1 - o$level) / 2)) else NA
                # each curve over its own mean ± 4 SD, with 0 and the H1 effect in view
                lo <- c(0, dir * res$theta1); hi <- lo
                if (pr[["prec"]] > 0) { lo <- c(lo, dir * pr[["mean"]] - 4 / sqrt(pr[["prec"]])); hi <- c(hi, dir * pr[["mean"]] + 4 / sqrt(pr[["prec"]])) }
                if (last) { lo <- c(lo, mon$postMean[last] - 4 * sdL); hi <- c(hi, mon$postMean[last] + 4 * sdL) }
                x <- seq(min(lo), max(hi), length.out = 301)
                df <- NULL
                if (pr[["prec"]] > 0)
                    df <- data.frame(x = x, y = stats::dnorm(x, dir * pr[["mean"]], 1 / sqrt(pr[["prec"]])), series = .("Prior"), type = "line")
                if (last) df <- rbind(df, data.frame(x = x, y = stats::dnorm(x, mon$postMean[last], sdL),
                                                     series = sprintf(.("Posterior at look %d"), last), type = "line"))
                if (is.null(df)) return(NULL)
                return(list(df = df, points = FALSE, h = NULL, xlab = .("Effect")))
            }
            a <- o$priorA; b <- o$priorB
            last <- if (!is.null(mon)) nrow(mon) else 0
            x <- seq(0.001, 0.999, length.out = 400)
            if (last) {
                df <- rbind(data.frame(x = x, y = stats::dbeta(x, a + mon$xC[last], b + mon$nC[last] - mon$xC[last]),
                                       series = .("Control"), type = "line"),
                            data.frame(x = x, y = stats::dbeta(x, a + mon$xT[last], b + mon$nT[last] - mon$xT[last]),
                                       series = .("Treatment"), type = "line"))
            } else df <- data.frame(x = x, y = stats::dbeta(x, a, b), series = .("Prior (each arm)"), type = "line")
            list(df = df, points = FALSE, h = NULL, xlab = .("Response rate"))
        },

        .ocCurve = function(res) {
            o <- self$options
            if (o$model == "normal") {
                g <- seq(0, 1.5 * res$theta1, length.out = 31)
                oc <- lapply(g, function(th) bayes_normal_oc(res$des, th))
                return(list(p = res$dir * g, reject = vapply(oc, `[[`, 0, "reject"), en = vapply(oc, `[[`, 0, "en"),
                            p0 = 0, p1 = res$dir * res$theta1))
            }
            dir <- res$dir; w <- abs(o$pT - o$pC)
            g <- seq(o$pC, min(0.99, max(0.01, o$pT + dir * w / 2)), length.out = 21)
            oc <- lapply(g, function(p) if (dir > 0) bb2_oc(res$des, o$pC, p) else bb2_oc(res$des, p, o$pC))
            list(p = g, reject = vapply(oc, `[[`, 0, "reject"), en = 2 * vapply(oc, `[[`, 0, "en"), p0 = o$pC, p1 = o$pT)
        },

        .decisionLabel = function(d) {
            switch(d, efficacy = .("Stop: efficacy"), futility = .("Stop: futility"), success = .("Final: success"),
                   fail = .("Final: not successful"), after = .("After the stop"), .("Continue"))
        },

        .fail = function(text) {
            private$.notice(text, "inputProblem", jmvcore::NoticeType$ERROR)
            for (nm in c("summary", "looks", "oc", "bfDesign", "monitor", "monitorBin", "bfTest", "boundPlot",
                         "probPlot", "bfPlot", "densityPlot", "ocPlot"))
                self$results[[nm]]$setVisible(FALSE)
            invisible(NULL)
        },

        # Notice, not `type: Notice` in the yaml, which does not compile (jamovi skill)
        .notice = function(text, name, type) {
            n <- jmvcore::Notice$new(options = self$options, name = name, type = type)
            n$setContent(text)
            self$results$insert(1, n)
        },

        .optionList = function() {
            nm <- c("model", "nMax", "nLooks", "timing", "arms", "sigma", "theta1", "prior", "priorMean", "priorSD",
                    "pC", "pT", "bNMax", "priorA", "priorB", "delta", "gammaE", "effInterim", "futility", "gammaF",
                    "bfSamples", "rScale", "bfK", "bfNMin", "bfNMax", "d1", "nsim", "seed")
            stats::setNames(lapply(nm, function(n) self$options[[n]]), nm)
        },

        .summaryKeys = function() {
            if (self$options$model == "bf") c("model", "prior", "rule", "alpha", "power", "en0", "en1")
            else c("model", "prior", "rule", "fut", "looks", "alpha", "alphaNB", "power", "en0", "en1")
        },

        .summaryLabel = function(key) {
            bf <- self$options$model == "bf"
            switch(key, model = .("Model"), prior = .("Prior"), rule = .("Success rule"), fut = .("Futility rule"),
                   looks = .("Analyses at N"),
                   alpha = if (bf) .("P(evidence for H1) with no effect") else .("P(success) with no effect"),
                   alphaNB = .("P(success) with no effect, futility ignored"),
                   power = if (bf) .("P(evidence for H1) at the design effect") else .("P(success) under H1"),
                   en0 = .("Expected N with no effect"), en1 = if (bf) .("Expected N at the design effect") else .("Expected N under H1"))
        },

        .summaryValues = function(res) {
            o <- self$options
            f3 <- function(x) formatC(x, format = "f", digits = 3)
            f1 <- function(x) formatC(x, format = "f", digits = 1)
            if (o$model == "bf") {
                s0 <- res$sims[[1]]; s1 <- res$sims[[2]]
                return(list(model = if (res$two) .("Two-sample JZS Bayes factor t-test") else .("One-sample JZS Bayes factor t-test"),
                            prior = sprintf(.("Cauchy(0, %s) on the standardised effect"), o$rScale),
                            rule = sprintf(.("Stop when BF₁₀ ≥ %s (H1) or ≤ 1/%s (H0), from %d to %d observations%s"),
                                           o$bfK, o$bfK, o$bfNMin, o$bfNMax, if (res$two) .(" per group") else ""),
                            alpha = f3(s0[["pH1"]]), power = f3(s1[["pH1"]]), en0 = f1(s0[["asn"]]), en1 = f1(s1[["asn"]])))
            }
            dir <- res$dir
            if (o$model == "normal") {
                pr <- res$prior
                prior <- switch(o$prior, flat = .("Flat"),
                    sprintf(.("Normal, mean %s, SD %s (worth %s patients)"), format(signif(dir * pr[["mean"]], 3)),
                            format(signif(1 / sqrt(pr[["prec"]]), 3)), f1(pr[["prec"]] * res$des$sigmaU^2)))
                rule <- sprintf(.("P(effect %s %s | data) ≥ %s"), if (dir > 0) ">" else "<", o$delta, o$gammaE)
                looks <- paste(format(signif(res$des$nk, 4)), collapse = ", ")
                model <- if (o$arms == "two") .("Normal approximation, two arms (N total)") else .("Normal approximation, one arm")
            } else {
                prior <- sprintf(.("Beta(%s, %s) for each arm"), o$priorA, o$priorB)
                rule <- sprintf(if (dir > 0) .("P(pT − pC > %s | data) ≥ %s") else .("P(pC − pT > %s | data) ≥ %s"), o$delta, o$gammaE)
                looks <- paste0(paste(res$nk, collapse = ", "), " ", .("per arm"))
                model <- .("Two arms, binary outcome (beta-binomial)")
            }
            if (!o$effInterim) rule <- paste(rule, .("(final analysis only)"))
            fut <- if (isTRUE(o$futility)) sprintf(.("Stop if the predictive probability of success < %s"), o$gammaF) else .("None")
            oc <- res$oc
            list(model = model, prior = prior, rule = rule, fut = fut, looks = looks,
                 alpha = f3(oc$reject[1]), alphaNB = f3(res$alphaNB), power = f3(oc$reject[3]),
                 en0 = f1(oc$en[1]), en1 = f1(oc$en[3]))
        },

        .plotBounds = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(plot_series(st, theme$color, list(theme_only(ggtheme)), list(x = .("N"), y = .("Estimate"))))
            TRUE
        },

        .plotProb = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            x <- if (self$options$model == "binary") .("N per arm") else .("N")
            print(plot_series(st, theme$color, list(theme_only(ggtheme)), list(x = x, y = .("Probability"))))
            TRUE
        },

        .plotBF = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(plot_series(st, theme$color, list(theme_only(ggtheme)), list(x = .("Observations n"), y = .("BF₁₀ (log scale)"))))
            TRUE
        },

        .plotDensity = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(plot_series(st, theme$color, list(theme_only(ggtheme)), list(x = st$xlab, y = .("Density"))))
            TRUE
        },

        .plotOC = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(reject = .("P(success)"), en = .("Expected N"),
                           x = if (self$options$model == "binary") .("Treatment rate") else .("True effect"))
            print(plot_ph2_oc(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        })
)
