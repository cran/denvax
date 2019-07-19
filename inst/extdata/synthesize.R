require(denvax)
require(jsonlite)

# expected to be called with two arguments
# first: json with model parameters
# second: target rds file to store synthesized populations
.args <- commandArgs(trailingOnly = T)

# load the model parameters; considers first entry in json by default
if (length(.args) > 2) {
  tar <- as.integer(.args[2])
} else {
  tar <- 1L
}
pars <- as.list(
  jsonlite::read_json(.args[1], simplifyVector = T)[tar, c("f_H", "f_L", "p_H")]
)

# create a population with the sample defaults
warning("creating synthetic population; this will take some time to complete.")
results <- denvax::synthetic.pop(pars)

saveRDS(results, file = tail(.args, 1))
