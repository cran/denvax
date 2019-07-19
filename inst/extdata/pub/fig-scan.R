suppressPackageStartupMessages({
  require(data.table)
  require(jsonlite)
  require(ggplot2)
  require(cowplot)
})

#' @example
#' .args <- c("results/fig-scales.rda", "results/scan-rois.rds", "results/scan-params.rds", "results/fig-scan.png")
.args <- commandArgs(trailingOnly = T)

load(.args[1])

refpars <- readRDS(.args[3])[, seropos := 1 - seroneg9 ]

reducedrois <- readRDS(.args[2])[mechanism == "binary"][
  refpars, on = .(foi, disparity)
][, .(roi = max(roi)),
  by = .(foi, disparity, seropos, p_H, log10OR, nu, tau)
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

p <- ggplot(reducedrois) + aes(
  log10(tau), nu,
  fill = roi, z = roi, color = c("neg", "pos")[(roi > 0) + 1]
) +
  facet_grid(
    disparity ~ foi,
    labeller = label_bquote(
      cols = .(
        sprintf("%s\n%d%%",
          ifelse(foi == "md", "% Seropositive 9-year-olds", ""),
          foiref[foi] * 100)
      ),
      rows = atop(
        .(
          disref[[disparity]]),
          (.(pHref[disparity]) * ", " * .(ORref[disparity]))
        )
    ), scales = "free"
  ) +
  geom_heat +
  thm + scale_contour + scale_roi + scale_nu + scale_tau

save_plot(tail(.args, 1), p, ncol = 3, nrow = 3, base_height = 4 * 2 / 3)
