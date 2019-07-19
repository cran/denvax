require(denvax)
require(data.table)

# expected to be called with two arguments
# first: an rds containing a matrix. each row should correspond to a life history
#  columns correspond to a year in life.
#
# second: target rds file to store synthesized populations
.args <- commandArgs(trailingOnly = T)

# load the lifehistory
lifehistories <- readRDS(.args[1])

# compute the proportion coefficients
probs <- denvax::nPxA(lifehistories)

saveRDS(probs, file = tail(.args, 1))
