# start to finish analysis

require(data.table)
require(denvax)
require(ggplot2)
require(cowplot)

.args <- c("inputs/serodata-morrison2010.csv", "results/ROI-morrison.png")
# your source serodata, in this example the Morrison 2010 data set
# and the plot target
.args <- commandArgs(trailingOnly = T) # or use from the command line

# read in the data
serodata <- fread(.args[1])

# fit the data
fit <- with(serodata, serofit(sero = Seropositive, N = Number, age.min = Age))

# simulate life histories
lh <- synthetic.pop(fit)

# estimate population probabilites for various life trajectories
probs <- nPxA(lh)

# compute the ROI equation coefficients
rc <- ROIcoeffs(probs, Ls = 10)

# use those coefficients to determine ROI surfaces
res <- ROI(rc,
  nus = seq(0.3, 0.9, by = 0.3),
  taus = 10 ^ seq(-1.5, -.5, by = .5)
)

# visualize results
p <- ggplot(res[L %in% c(1, 5, 10)]) +
  aes(
    x = A, y = roi,
    linetype = mechanism, alpha = factor(L),
    group = interaction(mechanism, L)
  ) +
  facet_grid(nu ~ tau, labeller = labeller(
    nu = function(nu) {
      paste0(c("\n", "\u03BD = \n", "\n"), sprintf("%.2g", as.numeric(nu)))
    },
    tau = function(tau) {
      paste0(c("\n", "\u03C4 = \n", "\n"), sprintf("%.2g", as.numeric(tau)))
    }
  )) +
  geom_line() + theme_minimal() +
  scale_y_continuous("ROI") +
  scale_x_continuous("Initial Testing Age") +
  scale_alpha_discrete("Max. # of Tests", breaks = c(1, 5, 10)) +
  scale_linetype_discrete("Vaccine Mechanism") +
  coord_cartesian(ylim = c(-1, 2))

save_plot(tail(.args, 1), p, ncol = 3, nrow = 3)
