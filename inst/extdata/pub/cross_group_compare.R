suppressPackageStartupMessages({
  require(denvax); require(data.table); require(jsonlite)
})

#' expected to be called with at least two arguments
#' first: json file(s) with model parameters
#' last: target rds file name to store synthesized pop digests
#' command line invocation from the published usage gives:
#' @example
#' .args <- c(
#'   "results/fit-morrison2010.json",
#'   "results/fit-lazou2016.json",
#'   "results/proportions.rds"
#' )
.args <- commandArgs(trailingOnly = T)

jsonfiles <- head(.args, -1); outputfile <- tail(.args, 1)

allpars <- data.table::rbindlist( # bind together all the results from ...
  lapply(
    jsonfiles,           # for all the json files
    jsonlite::read_json, # read them in
    simplifyVector = T   # using the read_json simplification scheme
  )
)

# then for all the parameter combinations
# generate the same synthetic population (rngseed = 1)
# and determine the life history proportions
results <- allpars[,
  denvax::nPxA(denvax::synthetic.pop(.SD, rngseed = 1)),
  by = Country
]
# this will take a bit to run...

saveRDS(results, file = outputfile)
