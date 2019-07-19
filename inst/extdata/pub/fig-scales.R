suppressPackageStartupMessages({
  require(ggplot2)
})

.args <- commandArgs(trailingOnly = T)

scale_start_age <- scale_x_continuous(
  "Initial Age for Routine Testing",
  expand = expand_scale(add = 0.1)
)

scale_end_age   <- scale_y_continuous(
  "Maximum Age for Testing",
  expand = expand_scale(add = 0.1)
)

scale_num_tests   <- scale_y_continuous(
  "Maximum # of Tests",
  breaks = function(lims) seq(round(lims[1]), round(lims[2]), by = 1),
  expand = expand_scale(add = 0.1)
)

scale_tau <- scale_x_continuous(
  expression(tau * ", Test Cost Fraction (log scale)"),
  labels = function(b) sprintf("%0.2f", 10 ^ b)
)

scale_nu <- scale_y_continuous(
  expression(nu * ", Vaccine Cost Fraction"),
  labels = function(b) sprintf("%0.1f", b)
)

scale_roi       <- scale_fill_gradient2(
  "ROI",
  limits = c(-1, 2), breaks = seq(-1, 2, by = 0.5),
  guide = guide_colorbar(order = 1)
)

pm <- function(f) sprintf("\U00B1%0.2f", f)

scale_contour <- list(
  scale_linetype_manual(
    breaks = c("zero", "c25", "c50"),
    labels = c(zero = "0.0", c25 = pm(0.25), c50 = pm(0.5)),
    values = c(zero = "solid", c25 = "dashed", c50 = "dotted"),
    guide = guide_legend(
      title = NULL,
      override.aes = list(color = rep("black", 3)),
      order = 5
    )
  ),
  scale_color_manual(
    values = c(zero = "black", pos = "dodgerblue", neg = "firebrick"),
    guide = "none"
  ),
  scale_alpha_manual(
    values = c(zero = 0.5, c25 = 0.5, c50 = 1),
    guide = "none"
  )
)

thm <- theme_minimal()

geom_heat <- list(
  geom_tile(aes(color = NULL)),
  geom_contour(
    aes(linetype = "zero", color = "zero", alpha = "zero"),
    breaks = 0
  ),
  geom_contour(
    aes(linetype = "c25", alpha = "c25"),
    breaks = c(-0.25, 0.25)
  ),
  geom_contour(
    aes(linetype = "c50", alpha = "c50"),
    breaks = c(-0.5, 0.5)
  )
)

save(list = ls(), file = tail(.args, 1))
