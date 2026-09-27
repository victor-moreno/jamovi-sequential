
anytimeClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "anytimeClass",
    inherit = anytimeBase,
    private = list(

        .init = function() {
            o <- self$options
            r <- self$results
            for (key in private$.summaryKeys())
                r$summary$addRow(rowKey = key, values = list(name = private$.summaryLabel(key)))
            r$monitor$getColumn("stat")$setTitle(if (o$sdKnown) "z" else "t")
            if (is.null(o$dep)) for (nm in c("ePlot", "csPlot")) r[[nm]]$setVisible(FALSE)
        },

        .run = function() {
            o <- self$options
            r <- self$results
            d <- r$summary$state
            if (is.null(d)) {
                d <- tryCatch(av_compute(private$.optionList()), error = function(e) e)
                if (inherits(d, "error")) return(private$.fail(conditionMessage(d)))
                r$summary$setState(d)
            }
            vals <- private$.summaryValues(d)
            for (key in names(vals)) r$summary$setRow(rowKey = key, values = list(value = vals[[key]]))
            r$summary$setNote("method", if (d$method == "exact")
                .("Exact probabilities (numerical integration over a test after every observation, or every pair for two samples).")
                else sprintf(.("Simulated with %d runs (seed %d), testing after every observation (every pair for two samples)."),
                             o$nsim, o$seed), init = FALSE)

            img <- r$curvePlot
            if (o$curvePlot && is.null(img$state)) {
                cv <- d$curve
                df <- data.frame(x = cv$n, y = cv$reject, series = .("Anytime-valid, every observation"), type = "line")
                img$setState(list(df = df, points = FALSE, h = data.frame(y = o$power), ylim = c(0, 1),
                                  mark = if (!is.na(d$nPlan)) c(x = d$nPlan, y = o$power) else NULL))
            }

            if (is.null(o$dep)) return(invisible(NULL))
            a <- tryCatch(private$.analyse(d), error = function(e) e)
            if (inherits(a, "error")) {
                private$.notice(conditionMessage(a), "dataProblem", jmvcore::NoticeType$WARNING)
                for (nm in c("monitor", "ePlot", "csPlot")) r[[nm]]$setVisible(FALSE)
                return(invisible(NULL))
            }
            ok <- which(a$ok)
            img <- r$ePlot
            if (o$ePlot && is.null(img$state)) {
                df <- data.frame(x = a$n[ok], y = a$e[ok], series = .("E-value"), type = "line")
                img$setState(list(df = df, points = length(ok) <= 200, logy = TRUE,
                                  h = data.frame(y = c(1 / o$alpha, 1)),
                                  mark = if (a$reject) c(x = a$at, y = a$e[a$at]) else NULL))
            }
            img <- r$csPlot
            if (o$csPlot && is.null(img$state)) {
                fin <- ok[is.finite(a$lower[ok]) & is.finite(a$upper[ok])]
                if (length(fin)) {
                    df <- data.frame(n = a$n[fin], est = a$est[fin], lower = a$lower[fin], upper = a$upper[fin])
                    # clip the first, very wide intervals so that the later ones stay readable
                    span <- stats::quantile(c(df$lower, df$upper), c(0.1, 0.9))
                    w <- diff(span)
                    img$setState(list(df = df, null = if (o$samples == "one") o$mu0 else 0,
                                      ylim = c(span[[1]] - w, span[[2]] + w)))
                }
            }
        },

        .analyse = function(d) {
            o <- self$options
            dat <- self$data
            x <- jmvcore::toNumeric(dat[[o$dep]])
            g <- NULL
            if (o$samples == "two") {
                if (is.null(o$group)) stop(.("A two-sample test needs a grouping variable with two levels."))
                gf <- droplevels(as.factor(dat[[o$group]]))
                if (nlevels(gf) != 2) stop(.("The grouping variable must have exactly two levels."))
                g <- as.integer(gf)
                keep <- !is.na(x) & !is.na(g)
                g <- g[keep]
                levs <- levels(gf)
            } else keep <- !is.na(x)
            x <- x[keep]
            if (length(x) < 2) stop(.("The outcome needs at least two complete observations."))
            a <- av_analyse(x, g, o$mu0, if (o$sdKnown) o$sigma else NA, d$g, o$alpha)
            tb <- self$results$monitor
            for (i in av_rows(length(x), a$at, o$every)) {
                if (!a$ok[i]) next
                dec <- if (a$reject && i == a$at) .("Reject H0 (E ≥ 1/α)") else if (a$reject && i > a$at) .("After rejection")
                       else .("Continue")
                tb$addRow(rowKey = i, values = blank_na(list(n = i, est = a$est[i], stat = a$stat[i], e = a$e[i], p = a$p[i],
                    lower = if (is.finite(a$lower[i])) a$lower[i] else NA,
                    upper = if (is.finite(a$upper[i])) a$upper[i] else NA, decision = dec)))
            }
            note <- sprintf(.("Two-sided test and %s%% confidence sequence, valid at every n; the always-valid p is the smallest 1/E so far."),
                            format(100 * (1 - o$alpha)))
            if (!is.null(g)) note <- paste(note, sprintf(.("Estimate = %s − %s."), levs[2], levs[1]))
            tb$setNote("av", note, init = FALSE)
            a
        },

        .fail = function(text) {
            private$.notice(text, "inputProblem", jmvcore::NoticeType$ERROR)
            for (nm in c("summary", "monitor", "curvePlot", "ePlot", "csPlot")) self$results[[nm]]$setVisible(FALSE)
            invisible(NULL)
        },

        # Notice, not `type: Notice` in the yaml, which does not compile (jamovi skill)
        .notice = function(text, name, type) {
            n <- jmvcore::Notice$new(options = self$options, name = name, type = type)
            n$setContent(text)
            self$results$insert(1, n)
        },

        .optionList = function() {
            nm <- c("samples", "sdKnown", "sigma", "alpha", "delta", "power", "mixing", "g", "nMax", "nsim", "seed")
            stats::setNames(lapply(nm, function(n) self$options[[n]]), nm)
        },

        .summaryKeys = function() c("test", "evalue", "rule", "nPlan", "asn", "median", "power", "alpha", "nSingle", "nFix"),

        .summaryLabel = function(key) {
            p <- self$options$power
            switch(key, test = .("Test"), evalue = .("E-value"), rule = .("Rule"),
                   nPlan = sprintf(.("N for power %s, monitoring after every observation"), p),
                   asn = .("Expected N at the design effect (to the horizon)"),
                   median = .("Median and 90th percentile of the stopping N"),
                   power = .("Power by the horizon"), alpha = .("Type I error by the horizon"),
                   nSingle = sprintf(.("N for power %s, one analysis with the e-value"), p),
                   nFix = sprintf(.("N for power %s, fixed-sample test"), p))
        },

        .summaryValues = function(d) {
            o <- self$options
            two <- o$samples == "two"
            f3 <- function(x) formatC(x, format = "f", digits = 3)
            tot <- if (two) paste0(" ", .("(total)")) else ""
            nn <- function(x) if (is.na(x)) sprintf(.("more than %d"), d$nMax) else paste0(format(x), tot)
            test <- paste0(if (two) .("Difference in means") else .("Mean"), ", ",
                           if (o$sdKnown) sprintf(.("known SD %s (mixture SPRT)"), o$sigma) else .("unknown SD (safe t-test)"))
            list(test = test,
                 evalue = sprintf(.("Gaussian mixture, g = %s (standardised effect scale)"), format(signif(d$g, 4))),
                 rule = sprintf(.("Reject H0: effect = 0 the first time E ≥ 1/α = %s"), format(signif(1 / o$alpha, 4))),
                 nPlan = nn(d$nPlan), asn = paste0(formatC(d$asn, format = "f", digits = 1), tot),
                 median = paste0(if (is.na(d$median)) "—" else d$median, " / ", if (is.na(d$q90)) "—" else d$q90, tot),
                 power = if (two) sprintf(.("%s (horizon %d in total)"), f3(d$power), d$nMax)
                         else sprintf(.("%s (horizon %d)"), f3(d$power), d$nMax),
                 alpha = f3(d$alpha), nSingle = nn(d$nSingle), nFix = nn(d$nFix))
        },

        .plotCurve = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(plot_series(st, theme$color, list(theme_only(ggtheme)),
                              list(x = .("N"), y = .("P(reject H0 by N) at the design effect"))))
            TRUE
        },

        .plotE = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(plot_series(st, theme$color, list(theme_only(ggtheme)), list(x = .("Observations n"), y = .("E-value (log scale)"))))
            TRUE
        },

        .plotCS = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(plot_cs(st, theme$color, list(theme_only(ggtheme)), list(x = .("Observations n"), y = .("Effect"))))
            TRUE
        })
)
