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

plot_bounds <- function(st, scale, cols, base, labels) {
    K <- nrow(st$looks)
    lk <- st$looks
    y_eff <- switch(scale, z = lk$effZ, p = lk$effP, effect = lk$effEst)
    y_fut <- switch(scale,
        z = lk$futZ,
        p = (if (st$sided == 2) 2 * stats::pnorm(lk$futZ, lower.tail = FALSE) else stats::pnorm(lk$futZ, lower.tail = FALSE)),
        effect = lk$futEst)
    df <- data.frame(x = lk$n, y = y_eff, what = labels$eff, grp = "e1")
    if (st$sided == 2 && scale == "z")
        df <- rbind(df, data.frame(x = lk$n, y = -lk$effZ, what = labels$eff, grp = "e2"))
    ok <- !is.na(y_fut)
    if (any(ok)) {
        # the futility line ends where it meets the efficacy boundary at the last look
        fx <- c(lk$n[ok], lk$n[K]); fy <- c(y_fut[ok], y_eff[K])
        df <- rbind(df, data.frame(x = fx, y = fy, what = labels$fut, grp = "f1"))
        if (st$sided == 2 && scale == "z")
            df <- rbind(df, data.frame(x = fx, y = -fy, what = labels$fut, grp = "f2"))
    }
    df$what <- factor(df$what, c(labels$eff, labels$fut))
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = what, linetype = what, group = grp)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 2.8) +
        ggplot2::scale_colour_manual(values = stats::setNames(cols[1:2], c(labels$eff, labels$fut)), drop = TRUE) +
        ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22"), c(labels$eff, labels$fut)), drop = TRUE) +
        ggplot2::labs(x = labels$x, y = labels$y, colour = NULL, linetype = NULL)
    if (scale == "z") p <- p + ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60")
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
