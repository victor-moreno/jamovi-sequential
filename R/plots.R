# Plot builders: plain data in, ggplot out. Colours and the base theme come from
# the caller (jamovi's palette and ggtheme), so these stay usable outside jamovi.

# jamovi's ggtheme is a list (theme + palette scales); keep only the theme so the
# plots' own scales are not replaced (see the jamovi skill, lifecycle-and-results).
theme_only <- function(ggtheme) {
    if (is.null(ggtheme)) return(NULL)
    if (inherits(ggtheme, "theme")) return(ggtheme)
    parts <- Filter(function(x) inherits(x, "theme"), ggtheme)
    if (length(parts) == 0L) NULL else parts[[1L]]
}

# Size (mm) of text drawn inside a plot, as a fraction of the theme's text size: fixed sizes
# that suit a 12 pt theme came out at 60 % of the axis text in jamovi, whose theme uses 16 pt.
ann_size <- function(base, rel = 0.8) {
    sz <- tryCatch(ggplot2::calc_element("text", base[[1]])$size, error = function(e) NULL)
    if (!is.numeric(sz)) sz <- 12
    rel * sz / ggplot2::.pt
}

# scale: z, b (B-value Z*sqrt(t): Whitehead's score up to a constant, so triangular
# boundaries are straight lines), p or effect; xvar: "n" (sample size / events) or "t".
plot_bounds <- function(st, scale, cols, base, labels, xvar = "n") {
    K <- nrow(st$looks)
    lk <- st$looks
    lk$x <- if (xvar == "t") lk$t else lk$n
    y_eff <- switch(scale, z = lk$effZ, b = lk$effZ * sqrt(lk$t), p = lk$effP, effect = lk$effEst)
    y_fut <- switch(scale,
        z = lk$futZ, b = lk$futZ * sqrt(lk$t),
        p = (if (st$sided == 2) 2 * stats::pnorm(lk$futZ, lower.tail = FALSE) else stats::pnorm(lk$futZ, lower.tail = FALSE)),
        effect = lk$futEst)
    mirror <- st$sided == 2 && scale %in% c("z", "b")
    df <- data.frame(x = lk$x, y = y_eff, what = labels$eff, grp = "e1")
    if (mirror)
        df <- rbind(df, data.frame(x = lk$x, y = -y_eff, what = labels$eff, grp = "e2"))
    ok <- !is.na(y_fut)
    if (any(ok)) {
        # the futility line ends where it meets the efficacy boundary at the last look
        fx <- c(lk$x[ok], lk$x[K]); fy <- c(y_fut[ok], y_eff[K])
        df <- rbind(df, data.frame(x = fx, y = fy, what = labels$fut, grp = "f1"))
        if (mirror)
            df <- rbind(df, data.frame(x = fx, y = -fy, what = labels$fut, grp = "f2"))
    }
    df$what <- factor(df$what, c(labels$eff, labels$fut))
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = what, linetype = what, group = grp)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 2.8) +
        ggplot2::scale_colour_manual(values = stats::setNames(cols[1:2], c(labels$eff, labels$fut)), drop = TRUE) +
        ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22"), c(labels$eff, labels$fut)), drop = TRUE) +
        ggplot2::labs(x = labels$x, y = labels$y, colour = NULL, linetype = NULL)
    if (scale %in% c("z", "b")) p <- p + ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60")
    if (xvar == "t") p <- p + ggplot2::scale_x_continuous(limits = c(0, 1))
    if (scale == "p" || (scale == "effect" && st$endpoint == "surv")) p <- p + ggplot2::scale_y_log10()
    p + base + ggplot2::theme(legend.position = "top")
}

plot_spend <- function(st, cols, base, labels) {
    df <- data.frame(t = c(0, st$t), y = c(0, st$cumAlpha / st$alpha), what = labels$alpha)
    if (!is.null(st$cumBeta))
        df <- rbind(df, data.frame(t = c(0, st$t), y = c(0, st$cumBeta / st$beta), what = labels$beta))
    df$what <- factor(df$what, unique(df$what))
    ggplot2::ggplot(df, ggplot2::aes(t, y, colour = what, linetype = what)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey60", linetype = "dotted") +
        ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2.6) +
        ggplot2::scale_colour_manual(values = cols[seq_along(levels(df$what))]) +
        ggplot2::scale_linetype_manual(values = c("solid", "22")[seq_along(levels(df$what))]) +
        # clip visually, not by scale limits: a final value of 1 + 1e-7 would be dropped
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = labels$x, y = labels$y, colour = NULL, linetype = NULL) +
        base + ggplot2::theme(legend.position = "top")
}

# Power and expected sample size have different scales: two panels, never two axes.
plot_power <- function(st, cols, base, labels) {
    df <- rbind(data.frame(g = st$g, y = st$reject, panel = labels$power),
                data.frame(g = st$g, y = st$en, panel = labels$en))
    df$panel <- factor(df$panel, c(labels$power, labels$en))
    ggplot2::ggplot(df, ggplot2::aes(g, y)) +
        ggplot2::geom_vline(xintercept = c(0, 1), linewidth = 0.3, colour = "grey60", linetype = "dotted") +
        ggplot2::geom_line(linewidth = 0.9, colour = cols[1]) +
        ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
        ggplot2::labs(x = labels$x, y = NULL) + base
}

# Phase II: probability of a "promising" conclusion and expected size against the true
# response rate (two panels, never two axes), with p0 and p1 marked.
plot_ph2_oc <- function(st, cols, base, labels) {
    df <- rbind(data.frame(p = st$p, y = st$reject, panel = labels$reject),
                data.frame(p = st$p, y = st$en, panel = labels$en))
    df$panel <- factor(df$panel, c(labels$reject, labels$en))
    ggplot2::ggplot(df, ggplot2::aes(p, y)) +
        ggplot2::geom_vline(xintercept = c(st$p0, st$p1), linewidth = 0.3, colour = "grey60", linetype = "dotted") +
        ggplot2::geom_line(linewidth = 0.9, colour = cols[1]) +
        ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
        ggplot2::labs(x = labels$x, y = NULL) + base
}

# Phase II plots share fixed colours with a meaning (the theme palette starts with black):
# promising green, not promising blue, continue grey, the current status red.
ph2_palette <- function() {
    c(yes = "#8FD18E", yesEarly = "#CFEBCD", no = "#7FA3E0", stop0 = "#D5E1F6", cont = "#BDBDBD",
      yesInk = "#2E8B3A", noInk = "#3B6FC4", now = "#D62728")
}

# Outcome map for a two-stage design: every (stage-1, stage-2) response count and the decision it
# leads to; columns stopped after stage 1 are drawn whole. With observed data (st$obs: n, x) the
# cells matching the current status are outlined in red and the final outcomes that can no longer
# be reached are faded.
plot_ph2_outcomes <- function(st, cols, base, labels) {
    pal <- ph2_palette()
    g <- expand.grid(x1 = 0:st$n1, x2 = 0:st$n2)
    what <- ifelse(g$x1 <= st$fut1, "stop0", ifelse(g$x1 >= st$eff1, "stop1",
                   ifelse(g$x1 + g$x2 >= st$eff2, "yes", "no")))
    lev <- c(stop0 = labels$stop0, no = labels$no, yes = labels$yes, stop1 = labels$stop1)
    keep <- names(lev)[names(lev) %in% what]
    g$what <- factor(lev[what], lev[keep])
    fills <- c(stop0 = pal[["stop0"]], no = pal[["no"]], yes = pal[["yes"]], stop1 = pal[["yesEarly"]])
    g$reach <- TRUE
    now <- NULL
    if (!is.null(st$obs)) {
        n <- st$obs[["n"]]; x <- st$obs[["x"]]; N <- st$n1 + st$n2
        if (n <= st$n1) {
            # still in stage 1: x responses so far, n1 - n patients to go in stage 1
            now <- data.frame(x1 = x, x2 = 0)
            g$reach <- g$x1 >= x & g$x1 <= x + st$n1 - n
        } else if (n <= N) {
            # in stage 2: the split of the x responses between the stages is not known
            x1 <- max(st$fut1 + 1, x - (n - st$n1)):min(st$n1, x)
            now <- data.frame(x1 = x1, x2 = x - x1)
            g$reach <- g$x1 %in% x1 & g$x2 >= x - g$x1 & g$x2 <= x - g$x1 + N - n
        }
    }
    step <- function(n) max(1, ceiling(n / 10))
    p <- ggplot2::ggplot(g, ggplot2::aes(x1, x2)) +
        ggplot2::geom_tile(ggplot2::aes(fill = what, alpha = reach), width = 0.9, height = 0.9) +
        ggplot2::scale_fill_manual(values = stats::setNames(fills[keep], lev[keep])) +
        ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none")
    if (!is.null(now))
        p <- p + ggplot2::geom_tile(data = now, ggplot2::aes(colour = labels$now), fill = NA,
                                    width = 0.9, height = 0.9, linewidth = 1.2) +
            ggplot2::scale_colour_manual(values = stats::setNames(pal[["now"]], labels$now))
    p + ggplot2::scale_x_continuous(breaks = seq(0, st$n1, by = step(st$n1)), expand = c(0, 0.3)) +
        ggplot2::scale_y_continuous(breaks = seq(0, st$n2, by = step(st$n2)), expand = c(0, 0.3)) +
        ggplot2::labs(x = labels$x, y = labels$y, fill = NULL, colour = NULL) + base +
        # the current-status key on a line of its own: beside the decisions it overflowed jamovi's plot
        ggplot2::theme(legend.position = "top", legend.box = "vertical", panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 1),
                        colour = ggplot2::guide_legend(order = 2))
}

# Phase II stopping boundaries: responses against patients. Shaded regions give the decision that
# the next look would take (for patients in (n[k-1], n[k]], the rule at look k), clipped so that
# responses never exceed patients; region edges lie half-way between response counts. Layers
# (st$show): "hyp" expected responses under p0 and p1 (slopes from the origin); with observed data
# (st$obs: n, x), "traj" the triangle of trajectories still possible and their range at each later
# look, "proj" the projection at the observed rate with its interval (st$rate, st$ci) as slopes.
plot_ph2_bounds <- function(st, cols, base, labels) {
    pal <- ph2_palette()
    ink <- tryCatch(ggplot2::calc_element("text", base[[1]])$colour, error = function(e) NULL)
    if (is.null(ink)) ink <- "grey30"
    show <- if (is.null(st$show)) c("traj", "proj", "hyp") else st$show
    s <- st$stages
    K <- length(s$nk); N <- s$nk[K]
    clip <- function(px, py, f) {                  # keep f(x, y) <= 0 (one half-plane)
        ox <- oy <- numeric(0)
        for (i in seq_along(px)) {
            j <- if (i == length(px)) 1 else i + 1
            a <- f(px[i], py[i]); b <- f(px[j], py[j])
            if (a <= 0) { ox <- c(ox, px[i]); oy <- c(oy, py[i]) }
            if ((a <= 0) != (b <= 0)) {
                t <- a / (a - b)
                ox <- c(ox, px[i] + t * (px[j] - px[i])); oy <- c(oy, py[i] + t * (py[j] - py[i]))
            }
        }
        list(x = ox, y = oy)
    }
    lev <- c(fut = labels$fut, cont = labels$cont, eff = labels$eff)
    reg <- list(); id <- 0
    add <- function(x0, x1, y0, y1, what) {
        if (y1 <= y0) return(invisible())
        p <- clip(c(x0, x1, x1, x0), c(y0, y0, y1, y1), function(x, y) y - x - 0.5)
        if (length(p$x) < 3) return(invisible())
        id <<- id + 1
        reg[[id]] <<- data.frame(x = p$x, y = p$y, id = id, what = lev[[what]])
    }
    for (k in seq_len(K)) {
        x0 <- if (k == 1) 0 else s$nk[k - 1]; x1 <- s$nk[k]
        f <- s$fut[k]; e <- s$eff[k]
        if (f >= 0) add(x0, x1, -0.5, f + 0.5, "fut")
        if (e <= x1) add(x0, x1, e - 0.5, x1 + 0.5, "eff")
        add(x0, x1, max(f, -1) + 0.5, min(e, x1 + 1) - 0.5, "cont")
    }
    reg <- do.call(rbind, reg)
    present <- lev[lev %in% unique(reg$what)]
    reg$what <- factor(reg$what, present)
    fills <- stats::setNames(c(pal[["no"]], pal[["cont"]], pal[["yes"]]), lev)[present]
    ok_f <- s$fut >= 0
    ok_e <- s$eff <= s$nk
    pts <- rbind(data.frame(n = s$nk[ok_f], y = s$fut[ok_f], what = "fut"),
                 data.frame(n = s$nk[ok_e], y = s$eff[ok_e], what = "eff"))
    p <- ggplot2::ggplot() +
        ggplot2::geom_polygon(data = reg, ggplot2::aes(x, y, group = id, fill = what), alpha = 0.45, colour = NA) +
        ggplot2::scale_fill_manual(values = fills, name = NULL) +
        ggplot2::geom_vline(xintercept = s$nk, colour = "grey55", linewidth = 0.3, linetype = "dotted") +
        ggplot2::geom_point(data = pts, ggplot2::aes(n, y), size = 2.6,
                            colour = ifelse(pts$what == "fut", pal[["noInk"]], pal[["yesInk"]]))
    lab_at_end <- function(y, text, colour, vjust = -0.5) ggplot2::annotate("text", x = N, y = y, label = text, hjust = 1.03,
                                                               vjust = vjust, size = ann_size(base, 0.75), colour = colour)
    if ("hyp" %in% show && !is.null(st$p0)) {
        p <- p + ggplot2::annotate("segment", x = 0, xend = N, y = 0, yend = st$p0 * N, colour = pal[["noInk"]],
                                   linetype = "42", linewidth = 0.5) +
            ggplot2::annotate("segment", x = 0, xend = N, y = 0, yend = st$p1 * N, colour = pal[["yesInk"]],
                              linetype = "42", linewidth = 0.5) +
            lab_at_end(st$p0 * N, sprintf("p0 = %s", st$p0), pal[["noInk"]]) +
            lab_at_end(st$p1 * N, sprintf("p1 = %s", st$p1), pal[["yesInk"]])
    }
    if (!is.null(st$obs)) {
        n0 <- st$obs[["n"]]; x0 <- st$obs[["x"]]
        if (n0 < N && "traj" %in% show) {
            later <- s$nk[s$nk > n0]
            p <- p + ggplot2::annotate("polygon", x = c(n0, N, N), y = c(x0, x0, x0 + N - n0),
                                       fill = ink, alpha = 0.05, colour = ink, linetype = "22", linewidth = 0.4) +
                ggplot2::annotate("segment", x = later, xend = later, y = x0, yend = x0 + later - n0,
                                  colour = ink, linewidth = 1) +
                # below the apex: above it, the label ran into the status text in wider fonts
                lab_at_end(x0 + N - n0, labels$reach, ink, vjust = 1.6)
        }
        # the rate is the slope: project the current status with the observed rate and the limits
        # of its interval, a fan inside the triangle of possible trajectories
        status <- sprintf(labels$now, x0, n0, st$rate)
        if (n0 < N && "proj" %in% show && !is.null(st$ci)) {
            lo <- st$ci[1]; hi <- st$ci[2]
            p <- p + ggplot2::annotate("polygon", x = c(n0, N, N), y = c(x0, x0 + lo * (N - n0), x0 + hi * (N - n0)),
                                       fill = pal[["now"]], alpha = 0.10, colour = NA) +
                ggplot2::annotate("segment", x = n0, xend = N, y = x0, yend = x0 + c(lo, hi) * (N - n0),
                                  colour = pal[["now"]], linewidth = 0.45) +
                ggplot2::annotate("segment", x = n0, xend = N, y = x0, yend = x0 + st$rate * (N - n0),
                                  colour = pal[["now"]], linetype = "42", linewidth = 0.6)
            status <- paste0(status, "\n", sprintf(labels$proj, format(100 * st$level), lo, hi))
        }
        # responses cannot exceed patients, so the top-left corner is always empty: the status goes there
        p <- p + ggplot2::annotate("point", x = n0, y = x0, shape = 4, size = 4, stroke = 1.5, colour = pal[["now"]]) +
            ggplot2::annotate("text", x = 0.02 * N, y = 0.98 * N, label = status, hjust = 0, vjust = 1, size = ann_size(base),
                              colour = pal[["now"]], lineheight = 1.1)
    }
    p + ggplot2::coord_cartesian(xlim = c(0, N + max(0.8, N / 60)), ylim = c(-0.8, N + 0.8), expand = FALSE) +
        ggplot2::labs(x = labels$x, y = labels$y) + base + ggplot2::theme(legend.position = "top")
}

# Interim monitoring: boundaries at the observed (and remaining planned) looks, with the
# observed Z path on top.
plot_monitor <- function(st, cols, base, labels, xvar = "n") {
    ob <- st$obs
    ob$x <- if (xvar == "t") ob$t else ob$n
    obsCol <- if (length(cols) >= 3) cols[3] else "grey30"
    # geom_text ignores the theme: take the theme's own text colour (dark mode)
    ink <- tryCatch(ggplot2::calc_element("text", base[[1]])$colour, error = function(e) NULL)
    if (is.null(ink)) ink <- "grey30"
    last <- ob[nrow(ob), ]
    plot_bounds(st, "z", cols, base, labels, xvar) +
        ggplot2::geom_line(data = ob, ggplot2::aes(x, z), inherit.aes = FALSE, colour = obsCol, linewidth = 0.8) +
        ggplot2::geom_point(data = ob, ggplot2::aes(x, z), inherit.aes = FALSE, colour = obsCol, size = 3.2, shape = 17) +
        ggplot2::annotate("text", x = last$x, y = last$z, label = labels$obs, vjust = -1.3, size = ann_size(base), colour = ink)
}

# SPRT: boundaries (LLR against n, or Z against V for the triangular test) with the observed
# path. st: bx, bup, blo (boundary lines), xT / yT (truncation segment), path (x, y), stop.
plot_sprt_path <- function(st, cols, base, labels) {
    ink <- tryCatch(ggplot2::calc_element("text", base[[1]])$colour, error = function(e) NULL)
    if (is.null(ink)) ink <- "grey30"
    b <- rbind(data.frame(x = st$bx, y = st$bup, what = labels$upper),
               data.frame(x = st$bx, y = st$blo, what = labels$lower))
    b$what <- factor(b$what, c(labels$upper, labels$lower))
    p <- ggplot2::ggplot() +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60", linetype = "dotted") +
        ggplot2::geom_line(data = b, ggplot2::aes(x, y, colour = what), linewidth = 0.9) +
        ggplot2::scale_colour_manual(values = stats::setNames(cols[1:2], c(labels$upper, labels$lower)))
    if (!is.null(st$xT) && !is.na(st$xT))
        p <- p + ggplot2::annotate("segment", x = st$xT, xend = st$xT, y = st$yT[1], yend = st$yT[2],
                                   colour = ink, linetype = "22", linewidth = 0.6)
    if (!is.null(st$path)) {
        pc <- if (length(cols) >= 3) cols[3] else ink
        p <- p + ggplot2::geom_path(data = st$path, ggplot2::aes(x, y), colour = pc, linewidth = 0.6)
        if (nrow(st$path) <= 200)
            p <- p + ggplot2::geom_point(data = st$path, ggplot2::aes(x, y), colour = pc, size = 1.1)
        if (!is.null(st$stop))
            p <- p + ggplot2::annotate("point", x = st$stop[["x"]], y = st$stop[["y"]], shape = 21, size = 3.6,
                                       stroke = 1.2, colour = pc, fill = NA)
    }
    p + ggplot2::expand_limits(x = 0) +
        ggplot2::labs(x = labels$x, y = labels$y, colour = NULL) + base + ggplot2::theme(legend.position = "top")
}

# MaxSPRT: cases needed for a signal against the expected count (Poisson) or all cases
# (binomial), the null expectation, and the observed cumulative counts.
plot_maxsprt <- function(st, cols, base, labels) {
    ink <- tryCatch(ggplot2::calc_element("text", base[[1]])$colour, error = function(e) NULL)
    if (is.null(ink)) ink <- "grey30"
    b <- st$boundary
    xmax <- max(b$x, if (!is.null(st$path)) st$path$x)
    p <- ggplot2::ggplot() +
        ggplot2::annotate("segment", x = 0, xend = xmax, y = 0, yend = st$slope * xmax,
                          colour = "grey60", linetype = "dotted", linewidth = 0.5) +
        # group monitoring tests only at the looks: points, not a step line
        (if (st$continuous) ggplot2::geom_step(data = b, ggplot2::aes(x, y, colour = labels$signal), linewidth = 0.9, direction = "hv")
         else ggplot2::geom_point(data = b, ggplot2::aes(x, y, colour = labels$signal), size = 2.8, shape = 15))
    p <- p + ggplot2::scale_colour_manual(values = stats::setNames(cols[1], labels$signal))
    if (!is.null(st$path)) {
        pc <- if (length(cols) >= 3) cols[3] else ink
        p <- p + ggplot2::geom_path(data = st$path, ggplot2::aes(x, y), colour = pc, linewidth = 0.6) +
            ggplot2::geom_point(data = st$path, ggplot2::aes(x, y), colour = pc, size = 1.6)
        if (!is.null(st$signal))
            p <- p + ggplot2::annotate("point", x = st$signal[["x"]], y = st$signal[["y"]], shape = 21, size = 3.8,
                                       stroke = 1.2, colour = pc, fill = NA) +
                ggplot2::annotate("text", x = st$signal[["x"]], y = st$signal[["y"]], label = labels$signalled,
                                  vjust = -1.3, size = ann_size(base), colour = ink)
    }
    p + ggplot2::expand_limits(x = 0, y = 0) +
        ggplot2::labs(x = labels$x, y = labels$y, colour = NULL) + base + ggplot2::theme(legend.position = "top")
}

# Generic series plot. st$df: x, y, series, type ("line": line and points if st$points,
# "point": points only); st$h: reference lines (y, label); st$logy: log10 y axis; st$mark:
# a highlighted point (x, y).
plot_series <- function(st, cols, base, labels) {
    ink <- tryCatch(ggplot2::calc_element("text", base[[1]])$colour, error = function(e) NULL)
    if (is.null(ink)) ink <- "grey30"
    df <- st$df
    lev <- unique(df$series)
    df$series <- factor(df$series, lev)
    pal <- stats::setNames(rep_len(if (length(cols)) cols else "grey30", length(lev)), lev)
    if (!is.null(st$inkSeries)) pal[st$inkSeries] <- ink
    ln <- df[df$type != "point", ]
    pts <- if (isTRUE(st$points)) df else df[df$type == "point", ]
    p <- ggplot2::ggplot()
    if (!is.null(st$h) && nrow(st$h))
        p <- p + ggplot2::geom_hline(data = st$h, ggplot2::aes(yintercept = y), linetype = "22",
                                     colour = "grey55", linewidth = 0.45)
    if (nrow(ln)) p <- p + ggplot2::geom_line(data = ln, ggplot2::aes(x, y, colour = series), linewidth = 0.9)
    if (nrow(pts)) p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x, y, colour = series), size = 2.2)
    if (!is.null(st$mark))
        p <- p + ggplot2::annotate("point", x = st$mark[["x"]], y = st$mark[["y"]], shape = 21, size = 3.8,
                                   stroke = 1.2, colour = ink, fill = NA)
    if (isTRUE(st$logy)) p <- p + ggplot2::scale_y_log10()
    if (!is.null(st$ylim)) p <- p + ggplot2::coord_cartesian(ylim = st$ylim)
    p + ggplot2::scale_colour_manual(values = pal) +
        ggplot2::labs(x = labels$x, y = labels$y, colour = NULL) + base + ggplot2::theme(legend.position = "top")
}

# Confidence sequence: band of the anytime-valid interval against n, with the estimate and the
# null value.
plot_cs <- function(st, cols, base, labels) {
    ink <- tryCatch(ggplot2::calc_element("text", base[[1]])$colour, error = function(e) NULL)
    if (is.null(ink)) ink <- "grey30"
    df <- st$df
    ggplot2::ggplot(df, ggplot2::aes(n)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = cols[1], alpha = 0.25) +
        ggplot2::geom_hline(yintercept = st$null, linetype = "22", colour = "grey55", linewidth = 0.45) +
        ggplot2::geom_line(ggplot2::aes(y = est), colour = ink, linewidth = 0.7) +
        ggplot2::coord_cartesian(ylim = st$ylim) +
        ggplot2::labs(x = labels$x, y = labels$y) + base
}
