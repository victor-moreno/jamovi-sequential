
sprtClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "sprtClass",
    inherit = sprtBase,
    private = list(

        .init = function() {
            o <- self$options
            r <- self$results
            for (key in private$.summaryKeys())
                r$summary$addRow(rowKey = key, values = list(name = private$.summaryLabel(key)))
            labs <- c(.("H0"), .("Midway"), .("H1"), .("Beyond H1"))
            for (i in 1:4) r$oc$addRow(rowKey = i, values = list(label = labs[i]))
            sim <- o$model %in% c("t1", "t2")
            for (cn in c("pH1", "asn")) r$oc$getColumn(cn)$setSuperTitle(if (sim) .("Simulated") else .("Exact"))
            wald <- o$boundary == "wald" && !sim
            for (cn in c("waldPH1", "waldASN")) r$oc$getColumn(cn)$setVisible(wald)
            r$oc$getColumn("theta")$setTitle(private$.thetaTitle())
            statTitle <- if (o$boundary == "triangular") .("Z") else .("LLR")
            r$test$getColumn("stat")$setTitle(statTitle)
            r$path$getColumn("stat")$setTitle(statTitle)
            r$test$getColumn("estimate")$setTitle(switch(o$model, binom = .("Proportion"),
                norm1 = , t1 = .("Mean"), .("Difference")))
            if (!is.null(o$dep)) r$test$addRow(rowKey = 1, values = list())
            # a visibility expression cannot combine a variable with &&
            else r$path$setVisible(FALSE)
        },

        .thetaTitle = function() {
            switch(self$options$model, binom = .("True proportion"), norm1 = .("True mean"),
                   norm2 = .("True difference"), .("True effect (d)"))
        },

        .run = function() {
            o <- self$options
            r <- self$results
            # the t-tests simulate: cache the design on the summary (its clearWith lists the design)
            d <- r$summary$state
            if (is.null(d)) {
                d <- tryCatch(sprt_compute(private$.optionList(), o$nsim, o$seed), error = function(e) e)
                if (inherits(d, "error")) return(private$.fail(conditionMessage(d)))
                r$summary$setState(d)
            }
            vals <- private$.summaryValues(d)
            for (key in names(vals)) r$summary$setRow(rowKey = key, values = list(value = vals[[key]]))

            if (o$ocTable) {
                for (i in seq_len(nrow(d$oc))) {
                    x <- d$oc[i, ]
                    r$oc$setRow(rowKey = i, values = blank_na(list(theta = x$theta, pH1 = x$pH1, asn = x$asn,
                        waldPH1 = x$waldPH1, waldASN = x$waldASN)))
                }
                notes <- character(0)
                if (d$method == "simulation")
                    notes <- sprintf(.("Simulated with %d runs per row (seed %d); an untruncated test is stopped at 1000 observations."),
                                     o$nsim, o$seed)
                else if (d$boundary == "wald")
                    notes <- .("Exact values by dynamic programming (proportions) or numerical integration (means). Wald's approximations ignore the overshoot of the thresholds and any truncation.")
                else
                    notes <- .("Whitehead's boundaries are approximations for data observed one at a time; the exact error rates are shown.")
                r$oc$setNote("method", notes, init = FALSE)
            }

            a <- NULL
            if (!is.null(o$dep)) {
                a <- tryCatch(private$.analyse(d), error = function(e) e)
                if (inherits(a, "error")) {
                    private$.notice(conditionMessage(a), "dataProblem", jmvcore::NoticeType$WARNING)
                    for (nm in c("test", "path")) r[[nm]]$setVisible(FALSE)
                    a <- NULL
                }
            }

            img <- r$pathPlot
            if (is.null(img$state) && o$pathPlot) img$setState(sprt_plot_state(d, a))
            img <- r$ocPlot
            if (is.null(img$state) && o$ocPlot) img$setState(sprt_oc_curve(d))
        },

        .analyse = function(d) {
            o <- self$options
            dat <- self$data
            x <- jmvcore::toNumeric(dat[[o$dep]])
            g <- NULL
            two <- o$model %in% c("norm2", "t2")
            if (two) {
                if (is.null(o$group)) stop(.("A two-sample test needs a grouping variable with two levels."))
                gf <- droplevels(as.factor(dat[[o$group]]))
                if (nlevels(gf) != 2) stop(.("The grouping variable must have exactly two levels."))
                g <- as.integer(gf)
                keep <- !is.na(x) & !is.na(g)
                g <- g[keep]
                levs <- levels(gf)
            } else keep <- !is.na(x)
            x <- x[keep]
            if (length(x) < 1) stop(.("The outcome has no complete observations."))
            a <- sprt_data(d, x, g)
            k <- a$n
            dec <- switch(a$decision, H1 = .("Accept H1 (reject H0)"), H0 = .("Accept H0"), .("Continue sampling"))
            if (a$truncated) dec <- paste(dec, .("at the maximum sample size"))
            self$results$test$setRow(rowKey = 1, values = blank_na(list(
                N = a$N, n = if (a$decision == "continue") NA else k, decision = dec,
                stat = a$stat[k], lower = a$lower[k], upper = a$upper[k], estimate = a$estimate)))
            note <- .("The estimate at the decision is biased: tests that stop early for H1 tend to overestimate the effect.")
            if (two) note <- paste(note, sprintf(.("Difference = %s − %s."), levs[2], levs[1]))
            self$results$test$setNote("bias", note, init = FALSE)
            if (o$pathTable) {
                tb <- self$results$path
                for (i in seq_len(if (a$decision == "continue") a$N else k))
                    tb$addRow(rowKey = i, values = blank_na(list(n = i, V = if (is.null(a$V)) NA else a$V[i],
                        stat = a$stat[i], lower = a$lower[i], upper = a$upper[i])))
            }
            a
        },

        .fail = function(text) {
            private$.notice(text, "inputProblem", jmvcore::NoticeType$ERROR)
            for (nm in c("summary", "oc", "test", "path", "pathPlot", "ocPlot"))
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
            nm <- c("model", "p0", "p1", "mu0", "mu1", "delta1", "sigma", "d1", "tails", "alpha", "power",
                    "boundary", "truncate", "nMax", "truncRule")
            stats::setNames(lapply(nm, function(n) self$options[[n]]), nm)
        },

        .summaryKeys = function() c("test", "hyp", "bounds", "trunc", "alpha", "power", "asn0", "asn1", "nmax", "fixed"),

        .summaryLabel = function(key) {
            switch(key, test = .("Test"), hyp = .("Hypotheses"), bounds = .("Boundaries"),
                   trunc = .("Truncation"), alpha = .("Type I error"), power = .("Power"),
                   asn0 = .("Expected n under H0"), asn1 = .("Expected n under H1"),
                   nmax = .("Maximum n"), fixed = .("Fixed-sample n (normal approximation)"))
        },

        .summaryValues = function(d) {
            o <- self$options
            f3 <- function(x) formatC(x, format = "f", digits = 3)
            f2 <- function(x) formatC(x, format = "f", digits = 2)
            model <- switch(o$model, binom = .("Proportion, one sample"), norm1 = .("Mean, known SD, one sample"),
                            norm2 = .("Difference in means, known SD"), t1 = .("One-sample t-test"),
                            t2 = .("Two-sample t-test"))
            bnd <- if (d$boundary == "wald") .("Wald SPRT") else .("Triangular test (Whitehead)")
            par <- d$par
            hyp <- switch(o$model,
                binom = sprintf("H0: p = %s  vs  H1: p = %s", par$p0, par$p1),
                norm1 = sprintf("H0: μ = %s  vs  H1: μ = %s (σ = %s)", par$mu0, par$mu1, par$sigma),
                norm2 = sprintf("H0: δ = 0  vs  H1: δ = %s (σ = %s)", par$delta1, par$sigma),
                if (d$tails == "two") sprintf("H0: d = 0  vs  H1: |d| = %s", par$d1)
                else sprintf("H0: d = 0  vs  H1: d = %s", par$d1))
            if (o$model == "t1" && par$mu0 != 0) hyp <- paste0(hyp, sprintf(.(", d = (μ − %s) / σ"), par$mu0))
            bounds <- if (d$boundary == "wald") {
                th <- d$thresholds
                sprintf(.("Accept H1 if LLR ≥ %s; accept H0 if LLR ≤ %s"), f3(th[["upper"]]), f3(th[["lower"]]))
            } else {
                tp <- d$triangle
                sprintf(.("Accept H1 if Z ≥ a + cV; accept H0 if Z ≤ −a + 3cV (a = %s, c = %s; closes by V = %s)"),
                        f3(tp$a), f3(tp$c), f2(tp$a / tp$c))
            }
            trunc <- if (is.na(d$nMax)) .("None") else if (d$boundary == "triangular")
                sprintf(.("At n = %d: accept H1 if Z ≥ 2cV"), d$nMax)
                else sprintf(if (d$truncRule == "sign") .("At n = %d: accept H1 if the LLR is positive")
                             else .("At n = %d: accept H1 only at the upper threshold"), d$nMax)
            oc <- d$oc
            nmax <- if (!is.na(d$nMax)) format(d$nMax) else if (d$boundary == "triangular") format(max(oc$nEnd))
                    else .("None (open test)")
            two <- o$model %in% c("norm2", "t2")
            list(test = paste0(model, ", ", bnd), hyp = hyp, bounds = bounds, trunc = trunc,
                 alpha = sprintf("%s (%s %s)", f3(oc$pH1[1]), .("target"), o$alpha),
                 power = sprintf("%s (%s %s)", f3(oc$pH1[3]), .("target"), o$power),
                 asn0 = f2(oc$asn[1]), asn1 = f2(oc$asn[3]), nmax = nmax,
                 fixed = if (two) sprintf(.("%d (total, both groups)"), d$fixedN) else format(d$fixedN))
        },

        .plotPath = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(upper = .("Accept H1"), lower = .("Accept H0"),
                           x = if (st$triangular) .("Information V") else .("Observations n"),
                           y = if (st$triangular) .("Efficient score Z") else .("Log-likelihood ratio"))
            print(plot_sprt_path(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        },

        .plotOC = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(reject = .("P(accept H1)"), en = .("Expected sample size"), x = private$.thetaTitle())
            print(plot_ph2_oc(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        })
)
