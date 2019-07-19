suppressPackageStartupMessages({
  require(denvax); require(data.table)
})

.args <- c("results/scan-params.rds", "results/scan-rois.rds")
.args <- commandArgs(trailingOnly = T)

parcombns  <- readRDS(.args[1])
outputfile   <- tail(.args, 1)

props <- parcombns[,
  nPxA(synthetic.pop(.BY, rngseed = 5)),
  by = .(p_H, f_L, f_H, foi, disparity)
]

rcs <- props[,
  denvax::ROIcoeffs(.SD),
  by = .(foi, disparity)
]

# using those coefficients and the proposed cost fractions
# compute the ROI surfaces
rois <- rcs[,
  denvax::ROI(.SD), by = .(foi, disparity)
]

saveRDS(rois, file = tail(.args, 1))
