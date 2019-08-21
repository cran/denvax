suppressPackageStartupMessages({
  require(data.table)
  require(jsonlite)
  require(ggplot2)
  require(cowplot)
})

#' @example
#' .args <- c("results/fig-scales.rda", "results/scan-rois.rds", "results/scan-params.rds", "results/fig-scan.png")
.args <- commandArgs(trailingOnly = T)

# load pre-computed shared figure elements; e.g., scales, helper functions, labels
load(.args[1])

refpars <- readRDS(.args[3])[, seropos := 1 - seroneg9 ]

reducedrois <- readRDS(.args[2])[mechanism == "binary"][
  refpars, on = .(foi, disparity)
][, roi,
  by = .(foi, disparity, seropos, p_H, log10OR, nu, tau, A, L)
]

rename <- function(res, ref) { names(res) <- ref; return(res) }

foiref <- refpars[, seropos[1], by = .(ref = foi)][, rename(V1, ref)]
pHref <- refpars[, p_H[1], by = .(ref = disparity)][, rename(V1, ref)]
ORref <- refpars[, log10OR[1], by = .(ref = disparity)][, rename(V1, ref)]

disref <- list(
  md = quote(Disparity ~ (p[H] * ", " * log[10] * OR)),
  ll = quote(phantom("Disparity" ~ (p[H] * ", " * log[10] * OR))),
  ul = quote(phantom("Disparity" ~ (p[H] * ", " * log[10] * OR)))
)

plot.dt <- reducedrois[
  disparity == "ll" & foi == "ul" &
  (nu - 0.5) ^ 2 < 0.0001 &
  log10(tau) %in% seq(-2, -1, by = .5)
]

scale_start_age <- scale_x_continuous(
  "Initial Age for Routine Testing",
  breaks = 5:10,
  expand = expand_scale(add = 0.5)
)

scale_num_tests   <- scale_y_continuous(
  "Maximum # of Tests",
  breaks = function(lims) seq(round(lims[1]), round(lims[2]), by = 1),
  expand = expand_scale(add = 0.5)
)

p <- ggplot(plot.dt) + aes(
  A, L,
  fill = roi, z = roi, color = c("neg", "pos")[(roi > 0) + 1]
) + facet_grid(
  disparity + nu + foi ~ tau,
  labeller = label_bquote(
    rows = atop(
      nu * " = " * .(sprintf("%0.1f", nu)) * ", " * .(
        c(ll = "Low", md = "Mid", ul = "High")[disparity]
      ) * " Disparity",
      .(
        c(ll = 50, md = 70, ul = 90)[foi]
      ) * "% " * S ^ "+" * " in 9-year-olds"
    ),
    cols = tau * " = " * .(sprintf("%0.2f", tau))
  )
) +
  geom_heat +
  coord_cartesian(xlim = c(5, 10), ylim = c(1, 10)) +
  scale_num_tests + scale_start_age + scale_roi +
  scale_contour + thm + theme(
    panel.spacing.x = unit(12, "pt")
  )

save_plot(
  tail(.args, 1), p,
  ncol = plot.dt[, length(unique(tau))], nrow = 1,
  base_height = 3.5, base_width = 4
)
