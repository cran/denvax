suppressPackageStartupMessages({
  require(denvax); require(data.table); require(jsonlite)
})

#' expected to be called with at least two arguments
#' first: json file(s) with model parameters
#' last: target rds file name to store synthesized pop digests
#' command line invocation from the published usage gives:
#' .args <- c(
#'   "results/fit-morrison2010.json",
#'   "results/fit-lazou2016.json",
#'   "results/ranges.rds"
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

survp <- function(f, A) (1 - f) ^ A
failp <- function(f, A) 1 - survp(f, A)

allpars[, seropos9 := p_H * failp(f_H, 9) + (1 - p_H) * failp(f_L, 9) ]
allpars[, seroneg9 := 1 - seropos9 ]
allpars[, OR9 := (
  failp(f_H, 9) * survp(f_L, 9)
) / (
  failp(f_L, 9) * survp(f_H, 9)
) ]

scan_labels <- c("ll", "md", "ul")

serorange <- allpars[, {
  ll <- min(floor(seropos9 * 10)) / 10
  ul <- max(ceiling(seropos9 * 10)) / 10
  .(seroneg9 = 1 - seq(ll, ul, length.out = 3), foi = scan_labels)
} ]

disprange <- allpars[, {
  llfrac <- min(floor(p_H * 10) / 10)
  ulfrac <- max(ceiling(p_H * 10) / 10)
  llOR <- min(floor(log10(OR9)))
  ulOR <- max(ceiling(log10(OR9)))
  .(
    p_H = seq(llfrac, ulfrac, length.out = 3),
    log10OR = seq(ulOR, llOR, length.out = 3),
    # note the opposite direction here - corresponds to observed pattern in data
    disparity = scan_labels
  )
}]

result <- serorange[, as.data.table(disprange), by = .(seroneg9, foi)]
foicalc <- function(seroneg9, log10OR9, p_H) {
  OR9 <- 10 ^ log10OR9
  A <- -(1 - OR9) * (1 - p_H)
  B <- (1 - OR9) * (seroneg9 - p_H) - OR9
  C <- OR9 * seroneg9
  slA <- -(B + sqrt(B ^ 2 - 4 * A * C)) / (2 * A)
  f_L <- 1 - slA ^ (1 / 9)
  shA <- pmax( (seroneg9 - (1 - p_H) * slA) / p_H, 0)
  f_H <- 1 - shA ^ (1 / 9)
  list(f_L = f_L, f_H = f_H)
}

result[, c("f_L", "f_H") := foicalc(seroneg9, log10OR, p_H)]

saveRDS(result, file = outputfile)
