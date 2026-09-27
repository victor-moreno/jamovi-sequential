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
