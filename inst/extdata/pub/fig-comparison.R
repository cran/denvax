suppressPackageStartupMessages({
  require(data.table)
  require(jsonlite)
  require(ggplot2)
})

#' @example
#' .args <- c("results/fig-scales.rda", "results/example-rois.rds", "pub/malaysia_v_peru.json", "results/fig-compare.png")
.args <- commandArgs(trailingOnly = T)

load(.args[1])
rois <- readRDS(.args[2])
refnutau <- jsonlite::read_json(.args[3])

pars <- lapply(
  refnutau,        # change the input parameters...
  function(pars) list( # into their dimensionless versions
    nu  =  pars$V / pars$S,
    tau = pars$T / pars$S
  )
)

# X=\U1D708
# Y=\U1D70F

country_labels <- sprintf("%s: \U1D708 = %0.2f, Y = %0.2f",
  names(pars), c(pars[[1]]$nu, pars[[2]]$nu), c(pars[[1]]$tau, pars[[2]]$tau)
)
names(country_labels) <- names(pars)

rois[, lbl := sprintf(
  "'%s: '*nu*' = %0.2f, '*tau*' = %0.2f'", Country, nu, tau
)]

p <- ggplot(rois) +
  aes(
    x = A, y = A + L - 1,
    fill = roi, z = roi, color = c("neg", "pos")[(roi > 0) + 1]
  ) +
  facet_grid(mechanism ~ lbl, labeller = labeller(
    mechanism = c(ordinal = "Ordinal Test", binary = "Binary Test"),
    lbl   = label_parsed
  )) +
  geom_heat +
  thm + scale_contour + scale_roi + scale_start_age + scale_end_age

ggsave(tail(.args, 1), p, width = 6, height = 5, dpi = 600)
