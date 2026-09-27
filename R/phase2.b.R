
phase2Class <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "phase2Class",
    inherit = phase2Base,
    private = list(

        .init = function() {
            o <- self$options
            r <- self$results
            for (key in private$.summaryKeys())
                r$summary$addRow(rowKey = key, values = list(name = private$.summaryLabel(key)))
            K <- tryCatch(private$.nLooks(), error = function(e) 1)
            for (k in seq_len(K)) r$stages$addRow(rowKey = k, values = list(stage = k))
            labs <- c(.("p0"), .("Midway"), .("p1"), .("Beyond p1"))
            for (i in 1:4) r$oc$addRow(rowKey = i, values = list(label = labs[i]))
            if (o$analyse) r$observed$addRow(rowKey = 1, values = list())
            bayes <- o$design %in% c("predictive", "posterior")
            ob <- r$observed
            ob$getColumn("umvue")$setVisible(o$design == "simon")
            for (cn in c("postMean", "crLower", "crUpper", "pGt0", "pp")) ob$getColumn(cn)$setVisible(bayes)
        },

        .nLooks = function() {
            o <- self$options
            switch(o$design, simon = 2, single = 1,
                fleming = { v <- parse_numbers(o$flemingN, "Cumulative sample sizes")
                            if (length(v)) length(v) else o$flemingStages },
                length(ph2_looks(o$bayesFirst, o$bayesEvery, o$bayesNmax)))
        },

        .run = function() {
            o <- self$options
            opts <- private$.optionList()
            simon <- NULL
            if (o$design == "simon") {
                # the search takes seconds: cache it on the element whose clearWith
                # lists exactly the options it depends on
                simon <- self$results$candidates$state
                if (is.null(simon)) {
                    simon <- tryCatch(ph2_simon(o$p0, o$p1, o$alpha, 1 - o$power, o$nMax), error = function(e) e)
                    if (!inherits(simon, "error")) self$results$candidates$setState(simon)
                }
                if (inherits(simon, "error")) return(private$.fail(conditionMessage(simon)))
            }
            res <- tryCatch(ph2_compute(opts, simon), error = function(e) e)
            if (inherits(res, "error")) return(private$.fail(conditionMessage(res)))

            st <- res$stages
            K <- length(st$nk)
            vals <- private$.summaryValues(res)
            for (key in names(vals))
                self$results$summary$setRow(rowKey = key, values = list(value = vals[[key]]))

            tb <- self$results$stages
            for (k in seq_len(K)) {
                if (k > tb$rowCount) tb$addRow(rowKey = k, values = list(stage = k))
                tb$setRow(rowKey = k, values = blank_na(list(
                    n = st$nk[k],
                    fut = if (st$fut[k] >= 0) st$fut[k] else NA,
                    eff = if (st$eff[k] <= st$nk[k]) st$eff[k] else NA,
                    stop0 = res$oc0$pEff[k] + res$oc0$pFut[k],
                    stop1 = res$oc1$pEff[k] + res$oc1$pFut[k])))
            }
            if (o$design %in% c("predictive", "posterior"))
                tb$setNote("rule", private$.bayesRule(), init = FALSE)

            if (o$ocTable)
                for (i in seq_len(nrow(res$oc)))
                    self$results$oc$setRow(rowKey = i, values = as.list(res$oc[i, ]))

            if (o$design == "simon" && o$candidates) {
                ct <- self$results$candidates
                cand <- res$simon$cand
                for (i in seq_len(nrow(cand))) {
                    x <- cand[i, ]
                    ct$addRow(rowKey = i, values = blank_na(list(r1 = x$r1, n1 = x$n1, r = x$r, n = x$n,
                        en0 = x$en0, pet0 = x$pet0, alpha = x$alpha, power = x$power, qLo = x$qLo, qHi = x$qHi)))
                }
                ct$setNote("adm", .("Best design for each maximum N. Admissible designs minimise q·N + (1 − q)·EN(p0) for q in the range shown (Jung et al. 2004)."), init = FALSE)
            }

            obs <- NULL
            if (o$analyse) {
                obs <- tryCatch(ph2_observe(opts, res, o$obsN, o$obsX, o$ciLevel), error = function(e) e)
                if (inherits(obs, "error")) {
                    private$.notice(conditionMessage(obs), "observedProblem", jmvcore::NoticeType$WARNING)
                    self$results$observed$setVisible(FALSE)
                    obs <- NULL
                } else private$.fillObserved(obs, res)
            }

            img <- self$results$ocPlot
            if (is.null(img$state) && o$ocPlot) {
                pg <- seq(0.005, 0.995, by = 0.005)
                oc <- lapply(pg, function(p) ms_oc(st$nk, st$fut, st$eff, p))
                img$setState(list(p = pg, reject = vapply(oc, `[[`, 0, "reject"),
                                  en = vapply(oc, `[[`, 0, "en"), p0 = o$p0, p1 = o$p1))
            }
            img <- self$results$boundPlot
            if (is.null(img$state)) {
                # the interval of the observed rate is drawn as slopes from the current status
                show <- c("traj", "proj", "hyp")[c(o$boundTraj, o$boundProj, o$boundHyp)]
                img$setState(list(stages = st, obs = if (is.null(obs)) NULL else c(n = obs$n, x = obs$x),
                                  rate = if (is.null(obs)) NULL else obs$naive,
                                  ci = if (is.null(obs)) NULL else c(obs$lower, obs$upper),
                                  level = o$ciLevel, p0 = o$p0, p1 = o$p1, show = show))
            }
            # outcome map: every (stage-1, stage-2) outcome of a two-stage design and its decision
            img <- self$results$outcomePlot
            if (o$outcomePlot) {
                if (K != 2) {
                    # shown where the map would be, not at the top of the results
                    private$.notice(.("The outcome map needs a two-stage design (Simon, or two stages or looks)."),
                                    "outcomeMap", jmvcore::NoticeType$INFO,
                                    at = match("outcomePlot", self$results$itemNames))
                    img$setVisible(FALSE)
                } else if (is.null(img$state))
                    img$setState(list(n1 = st$nk[1], n2 = st$nk[2] - st$nk[1], fut1 = st$fut[1],
                                      eff1 = st$eff[1], eff2 = st$eff[2],
                                      obs = if (is.null(obs)) NULL else c(n = obs$n, x = obs$x)))
            }
        },

        .fail = function(text) {
            private$.notice(text, "inputProblem", jmvcore::NoticeType$ERROR)
            for (nm in c("summary", "stages", "candidates", "oc", "observed", "ocPlot", "boundPlot", "outcomePlot"))
                self$results[[nm]]$setVisible(FALSE)
            invisible(NULL)
        },

        # Notice, not `type: Notice` in the yaml, which does not compile (jamovi skill)
        .notice = function(text, name, type, at = 1) {
            n <- jmvcore::Notice$new(options = self$options, name = name, type = type)
            n$setContent(text)
            self$results$insert(at, n)
        },

        .fillObserved = function(obs, res) {
            o <- self$options
            dec <- switch(if (is.na(obs$decision)) "none" else obs$decision,
                promising = .("Promising: reject H0"), `not promising` = .("Not promising: stop"),
                continue = .("Continue"), none = .("Not a planned analysis"))
            self$results$observed$setRow(rowKey = 1, values = blank_na(list(
                n = obs$n, x = obs$x, decision = dec, naive = obs$naive,
                umvue = if (is.null(obs$umvue)) NA else obs$umvue,
                lower = obs$lower, upper = obs$upper, pvalue = obs$p.value,
                postMean = if (is.null(obs$postMean)) NA else obs$postMean,
                crLower = if (is.null(obs$crLower)) NA else obs$crLower,
                crUpper = if (is.null(obs$crUpper)) NA else obs$crUpper,
                pGt0 = if (is.null(obs$pGt0)) NA else obs$pGt0,
                pp = if (is.null(obs$pp)) NA else obs$pp)))
            note <- if (obs$adjusted)
                .("The p-value and the exact interval account for the two-stage stopping rule (Koyama & Chen 2008); UMVUE of Jung & Kim (2004).")
            else .("Exact binomial p-value and Clopper-Pearson interval; they do not account for interim looks.")
            if (!is.null(res$prior))
                note <- paste(note, sprintf(.("Posterior from the Beta(%s, %s) prior; predictive probability of a positive final analysis."),
                                            format(signif(res$prior[1], 3)), format(signif(res$prior[2], 3))))
            self$results$observed$setNote("inf", note, init = FALSE)
        },

        .bayesRule = function() {
            o <- self$options
            fut <- if (o$design == "predictive") sprintf(.("predictive probability of success < %s"), o$thetaL)
                   else sprintf(.("P(p > p1 | data) < %s"), o$thetaL)
            eff <- if (o$thetaU >= 1) .("no early efficacy stop") else if (o$design == "predictive")
                sprintf(.("efficacy if predictive probability > %s"), o$thetaU)
                else sprintf(.("efficacy if P(p > p0 | data) > %s"), o$thetaU)
            paste0(.("Futility if "), fut, "; ", eff, "; ",
                   sprintf(.("final analysis promising if P(p > p0 | data) > %s."), o$thetaT))
        },

        .optionList = function() {
            nm <- c("design", "p0", "p1", "alpha", "power", "simonType", "admWeight", "nMax", "flemingStages",
                    "flemingN", "prior", "priorA", "priorB", "bayesNmax", "bayesFirst", "bayesEvery",
                    "thetaT", "thetaL", "thetaU")
            stats::setNames(lapply(nm, function(n) self$options[[n]]), nm)
        },

        .summaryKeys = function() {
            o <- self$options
            keys <- c("design", "hyp", "rule", "nmax", "alpha", "power", "en0", "en1", "pet0", "pet1")
            if (o$design %in% c("predictive", "posterior")) keys <- c(keys, "prior")
            keys
        },

        .summaryLabel = function(key) {
            switch(key, design = .("Design"), hyp = .("Hypotheses"), rule = .("Rule"),
                   nmax = .("Maximum sample size"), alpha = .("Type I error (attained)"),
                   power = .("Power (attained)"), en0 = .("Expected N under p0"),
                   en1 = .("Expected N under p1"), pet0 = .("P(early stop) under p0"),
                   pet1 = .("P(early stop) under p1"), prior = .("Prior"))
        },

        .summaryValues = function(res) {
            o <- self$options
            st <- res$stages
            f3 <- function(x) formatC(x, format = "f", digits = 3)
            f2 <- function(x) formatC(x, format = "f", digits = 2)
            des <- switch(o$design,
                simon = paste0(.("Simon two-stage"), ", ", switch(o$simonType, optimal = .("optimal"),
                    minimax = .("minimax"), admissible = paste0(.("admissible"), ", q = ", o$admWeight))),
                single = .("Single-stage exact binomial"),
                fleming = .("Fleming multistage"),
                predictive = .("Bayesian predictive probability (Lee & Liu)"),
                posterior = .("Bayesian posterior probability"))
            rule <- if (o$design == "simon") {
                d <- res$chosen
                sprintf(.("Stop after %d patients if ≤ %d responses; promising if > %d responses in %d"),
                        d$n1, d$r1, d$r, d$n)
            } else if (o$design == "single") {
                sprintf(.("Promising if ≥ %d responses in %d"), st$eff[1], st$nk[1])
            } else sprintf(.("Looks after %s patients (see Stopping rules)"), paste(st$nk, collapse = ", "))
            v <- list(design = des,
                      hyp = sprintf("H0: p ≤ %s  vs  H1: p ≥ %s", o$p0, o$p1),
                      rule = rule, nmax = format(max(st$nk)),
                      alpha = sprintf("%s (%s %s)", f3(res$oc0$reject), .("target"), o$alpha),
                      power = sprintf("%s (%s %s)", f3(res$oc1$reject), .("target"), o$power),
                      en0 = f2(res$oc0$en), en1 = f2(res$oc1$en),
                      pet0 = f3(res$oc0$pet), pet1 = f3(res$oc1$pet))
            if (!is.null(res$prior))
                v$prior <- sprintf("Beta(%s, %s)", format(signif(res$prior[1], 3)), format(signif(res$prior[2], 3)))
            v[private$.summaryKeys()]
        },

        .plotOC = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(reject = .("P(declared promising)"), en = .("Expected sample size"),
                           x = .("True response rate"))
            print(plot_ph2_oc(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        },

        .plotOutcomes = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(stop0 = .("Stopped after stage 1: not promising"), stop1 = .("Stopped after stage 1: promising"),
                           no = .("Not promising"), yes = .("Promising"), now = .("Current status"),
                           x = sprintf(.("Responses in stage 1 (of %d)"), st$n1),
                           y = sprintf(.("Responses in stage 2 (of %d)"), st$n2))
            print(plot_ph2_outcomes(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        },

        .plotBounds = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            labels <- list(fut = .("Not promising"), eff = .("Promising"), cont = .("Continue"),
                           now = .("Now: %d of %d (rate %.2f)"), reach = .("Possible trajectories"),
                           proj = .("Red fan: projection with the %s%% interval %.2f–%.2f"),
                           x = .("Patients"), y = .("Responses"))
            print(plot_ph2_bounds(st, theme$color, list(theme_only(ggtheme)), labels))
            TRUE
        })
)
