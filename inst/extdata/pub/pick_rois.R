suppressPackageStartupMessages({
  require(denvax); require(data.table); require(jsonlite)
})

#' expected to be called with three arguments
#' first: reference rds with digested life history proportions
#' second: json file with
#'   (1) keys for each of the countries to compare
#'   (2) that point to objects with costs for
#'   secondary infection (S), testing (T), a vaccination (V)
#' last: target rds file name to store resulting ROIs
#' command line invocation from the published usage gives:
#' @example
#' .args <- c(
#'   "results/proportions.rds",
#'   "pub/malaysia_v_peru.json",
#'   "results/example-rois.rds"
#' )
.args <- commandArgs(trailingOnly = T)

proportions  <- readRDS(.args[1])
tarCountries <- jsonlite::read_json(.args[2])
outputfile   <- tail(.args, 1)

pars <- lapply(
  tarCountries,        # change the input parameters...
  function(pars) list( # into their dimensionless versions
    nu =  pars$V / pars$S,
    tau = pars$T / pars$S
  )
)

rcs <- proportions[
  Country %in% names(pars), # for only the countries of interest...
  denvax::ROIcoeffs(.SD),   # compute the ROI surface coefficients
  by = Country                # n.b., this uses the default age, test limits
]

# using those coefficients and the proposed cost fractions
# compute the ROI surfaces
rois <- rcs[,
  denvax::ROI(.SD,
    nus  = pars[[Country[1]]]$nu,
    taus = pars[[Country[1]]]$tau
  ), by = Country
]

saveRDS(rois, file = tail(.args, 1))
