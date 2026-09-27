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

# Phase II stopping boundaries: responses against patients at each look.
plot_ph2_bounds <- function(st, cols, base, labels) {
    s <- st$stages
    ok_f <- s$fut >= 0
    ok_e <- s$eff <= s$nk
    df <- rbind(data.frame(n = s$nk[ok_f], y = s$fut[ok_f], what = labels$fut),
                data.frame(n = s$nk[ok_e], y = s$eff[ok_e], what = labels$eff))
    df$what <- factor(df$what, c(labels$eff, labels$fut))
    p <- ggplot2::ggplot(df, ggplot2::aes(n, y, colour = what, linetype = what)) +
        ggplot2::geom_step(linewidth = 0.8, direction = "hv") +
        ggplot2::geom_point(size = 2.8) +
        ggplot2::scale_colour_manual(values = stats::setNames(cols[1:2], c(labels$eff, labels$fut))) +
        ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22"), c(labels$eff, labels$fut))) +
        ggplot2::expand_limits(x = 0, y = 0) +
        ggplot2::labs(x = labels$x, y = labels$y, colour = NULL, linetype = NULL)
    if (!is.null(st$obs)) {
        # geom_text ignores the theme: take the theme's own text colour (dark mode)
        ink <- tryCatch(ggplot2::calc_element("text", base[[1]])$colour, error = function(e) NULL)
        if (is.null(ink)) ink <- "grey30"
        p <- p + ggplot2::annotate("point", x = st$obs[["n"]], y = st$obs[["x"]], shape = 4, size = 4, stroke = 1.3,
                                   colour = if (length(cols) >= 3) cols[3] else ink) +
            ggplot2::annotate("text", x = st$obs[["n"]], y = st$obs[["x"]], label = labels$obs, vjust = -1.2,
                              size = 3.5, colour = ink)
    }
    p + base + ggplot2::theme(legend.position = "top")
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
        ggplot2::annotate("text", x = last$x, y = last$z, label = labels$obs, vjust = -1.3, size = 3.5, colour = ink)
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
                                  vjust = -1.3, size = 3.5, colour = ink)
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

# Phase II outcome map for a two-stage design: every (stage-1, stage-2) response count and the
# decision it leads to. Columns stopped after stage 1 are drawn whole (stage 2 never happens).
plot_ph2_outcomes <- function(st, cols, base, labels) {
    g <- expand.grid(x1 = 0:st$n1, x2 = 0:st$n2)
    what <- ifelse(g$x1 <= st$fut1, "stop0", ifelse(g$x1 >= st$eff1, "stop1",
                   ifelse(g$x1 + g$x2 >= st$eff2, "yes", "no")))
    lev <- c(stop0 = labels$stop0, no = labels$no, yes = labels$yes, stop1 = labels$stop1)
    keep <- names(lev)[names(lev) %in% what]
    g$what <- factor(lev[what], lev[keep])
    fills <- c(stop0 = grDevices::adjustcolor(cols[2], alpha.f = 0.3), no = cols[2], yes = cols[1],
               stop1 = grDevices::adjustcolor(cols[1], alpha.f = 0.45))
    step <- function(n) max(1, ceiling(n / 10))
    ggplot2::ggplot(g, ggplot2::aes(x1, x2, fill = what)) +
        ggplot2::geom_tile(width = 0.9, height = 0.9) +
        ggplot2::scale_fill_manual(values = stats::setNames(fills[keep], lev[keep])) +
        ggplot2::scale_x_continuous(breaks = seq(0, st$n1, by = step(st$n1)), expand = c(0, 0.3)) +
        ggplot2::scale_y_continuous(breaks = seq(0, st$n2, by = step(st$n2)), expand = c(0, 0.3)) +
        ggplot2::labs(x = labels$x, y = labels$y, fill = NULL) + base +
        ggplot2::theme(legend.position = "top", panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}
