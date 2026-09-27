
maxsprtClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "maxsprtClass",
    inherit = maxsprtBase,
    private = list(

        .init = function() {
            o <- self$options
            r <- self$results
            for (key in private$.summaryKeys())
                r$summary$addRow(rowKey = key, values = list(name = private$.summaryLabel(key)))
            pois <- o$model == "poisson"
            r$thresholds$getColumn("look")$setTitle(if (pois && o$monitoring == "continuous") .("Expected cases (up to)")
                                                    else if (pois) .("Expected cases") else .("All cases"))
            r$monitor$getColumn("exp")$setTitle(if (pois) .("Expected") else .("All cases"))
            if (!is.null(o$cases) && !is.null(o$expected)) {
                dat <- if (nrow(self$data) > 0) self$data else self$readDataset()
                n <- sum(!is.na(dat[[o$cases]]) & !is.na(dat[[o$expected]]))
                for (i in seq_len(n)) r$monitor$addRow(rowKey = i, values = list(look = i))
            }
        },

        .run = function() {
            o <- self$options
            r <- self$results
            # the exact critical value takes seconds for long surveillance: cache it
            d <- r$summary$state
            if (is.null(d)) {
                d <- tryCatch(maxsprt_compute(private$.optionList()), error = function(e) e)
                if (inherits(d, "error")) return(private$.fail(conditionMessage(d)))
                r$summary$setState(d)
            }
            vals <- private$.summaryValues(d)
            for (key in names(vals)) r$summary$setRow(rowKey = key, values = list(value = vals[[key]]))
            if (d$alpha <= 0)
                private$.notice(.("No count can signal at this α with this surveillance length and minimum number of cases: lengthen the surveillance or lower the minimum."),
                                "noSignal", jmvcore::NoticeType$WARNING)

            pf <- d$performance
            for (i in seq_len(nrow(pf)))
                r$performance$addRow(rowKey = i, values = blank_na(as.list(pf[i, ])))
            unit <- if (d$model == "poisson") .("expected cases under H0") else .("cases")
            r$performance$setNote("unit", sprintf(.("Times in %s. The time to signal is the mean among the surveillances that signal; the length of surveillance counts those that do not as %s."),
                                                  unit, format(signif(d$size, 4))), init = FALSE)

            if (o$thresholdTable) {
                th <- d$thresholds
                tb <- r$thresholds
                if (d$model == "poisson" && d$monitoring == "continuous") {
                    for (i in seq_len(nrow(th)))
                        tb$addRow(rowKey = i, values = list(look = th$mu[i], signal = th$events[i]))
                    tb$setNote("rule", .("Signal as soon as the cases reach the threshold while the expected count is at most the value shown."), init = FALSE)
                } else {
                    for (i in seq_len(nrow(th)))
                        tb$addRow(rowKey = i, values = blank_na(list(look = th$look[i],
                            signal = if (is.finite(th$threshold[i])) th$threshold[i] else NA)))
                    if (any(!is.finite(th$threshold)))
                        tb$setNote("blank", .("Blank: no count can signal at that look."), init = FALSE)
                }
            }

            m <- NULL
            if (!is.null(o$cases) && !is.null(o$expected)) {
                m <- tryCatch(private$.monitor(d), error = function(e) e)
                if (inherits(m, "error")) {
                    private$.notice(conditionMessage(m), "dataProblem", jmvcore::NoticeType$WARNING)
                    r$monitor$setVisible(FALSE)
                    m <- NULL
                }
            }

            img <- r$boundaryPlot
            if (is.null(img$state) && o$boundaryPlot) {
                st <- list(boundary = maxsprt_boundary(d), continuous = d$monitoring == "continuous",
                           slope = if (d$model == "poisson") 1 else d$p0, poisson = d$model == "poisson")
                if (!is.null(m)) {
                    st$path <- data.frame(x = c(0, m$exp), y = c(0, m$obs))
                    k <- which(m$signal)[1]
                    if (!is.na(k)) st$signal <- c(x = m$exp[k], y = m$obs[k])
                }
                img$setState(st)
            }
            img <- r$powerPlot
            if (is.null(img$state) && o$powerPlot) {
                g <- seq(1, max(2, 1.25 * max(d$performance$rr)), length.out = 21)
                pf <- maxsprt_design_perf(d, g)
                img$setState(list(p = g, reject = pf$power, en = pf$timeToSignal, p0 = 1, p1 = 1))
            }
        },

        .monitor = function(d) {
            o <- self$options
            dat <- self$data
            m <- maxsprt_monitor(d, jmvcore::toNumeric(dat[[o$cases]]), jmvcore::toNumeric(dat[[o$expected]]))
            tb <- self$results$monitor
            for (i in seq_len(nrow(m))) {
                dec <- if (m$after[i]) .("After the signal") else if (m$signal[i]) .("Signal: reject H0")
                       else if (m$obs[i] < d$M) .("No signal (too few cases)")
                       else if (m$beyond[i]) .("No signal (beyond the planned length)")
                       else if (m$exp[i] >= d$size * (1 - 1e-9)) .("No signal: end of surveillance")
                       else .("No signal")
                if (i > tb$rowCount) tb$addRow(rowKey = i, values = list(look = i))
                tb$setRow(rowKey = i, values = blank_na(list(obs = m$obs[i], exp = m$exp[i],
                    rr = if (is.finite(m$rr[i])) m$rr[i] else NA, llr = m$llr[i], cv = d$cv, decision = dec)))
            }
            note <- if (d$model == "poisson") .("Relative risk = observed / expected cases.")
                    else sprintf(.("Relative risk = z·c / (n − c) with c risk-window cases of n and matching ratio z = %s."), o$z)
            if (d$monitoring == "group")
                note <- paste(note, .("The critical value controls α for the planned looks; testing at other times changes the attained α."))
            tb$setNote("rr", note, init = FALSE)
            m
        },

        .fail = function(text) {
            private$.notice(text, "inputProblem", jmvcore::NoticeType$ERROR)
            for (nm in c("summary", "performance", "thresholds", "monitor", "boundaryPlot", "powerPlot"))
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
            nm <- c("model", "alpha", "survLength", "N", "z", "M", "monitoring", "nLooks", "looks", "rr")
            stats::setNames(lapply(nm, function(n) self$options[[n]]), nm)
        },

        .summaryKeys = function() c("model", "hyp", "length", "M", "monitoring", "cv", "alpha"),

        .summaryLabel = function(key) {
            switch(key, model = .("Model"), hyp = .("Hypotheses"), length = .("Surveillance length"),
                   M = .("Minimum cases for a signal"), monitoring = .("Monitoring"),
                   cv = .("Critical value (LLR)"), alpha = .("Type I error (attained)"))
        },

        .summaryValues = function(d) {
            o <- self$options
            pois <- d$model == "poisson"
            mon <- if (d$monitoring == "continuous") .("Continuous") else
                sprintf(.("Group, looks at %s"), paste(format(signif(d$looks, 4)), collapse = ", "))
            list(model = if (pois) .("Poisson") else sprintf(.("Binomial, matching ratio z = %s"), o$z),
                 hyp = .("H0: RR = 1  vs  H1: RR > 1"),
                 length = if (pois) sprintf(.("%s expected cases"), format(signif(d$size, 4)))
                          else sprintf(.("%d cases"), as.integer(d$size)),
                 M = format(d$M), monitoring = mon,
                 cv = formatC(d$cv, format = "f", digits = 6),
                 alpha = sprintf("%s (%s %s)", formatC(d$alpha, format = "f", digits = 4), .("target"), o$alpha))
        },

        .plotBoundary = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(signal = .("Cases needed for a signal"), signalled = .("Signal"),
                           x = if (st$poisson) .("Expected cases under H0") else .("All cases"),
                           y = if (st$poisson) .("Observed cases") else .("Cases in the risk window"))
            print(plot_maxsprt(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        },

        .plotPower = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(reject = .("Power"), en = .("Expected time to signal"), x = .("Relative risk"))
            print(plot_ph2_oc(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        })
)
