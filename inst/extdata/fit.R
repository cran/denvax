require(denvax)
require(jsonlite)
require(data.table)

# expected to be called with two arguments
# first: csv with sero data
# second: target json file to store fit data
.args <- commandArgs(trailingOnly = T)

# this example assumes csv has headers Seropositve, Number, and Age (among any others)
# and that Age is either a year (like the Morrison 2010 data), or a range (like the L'Azou 2016 data)
serodata <- data.table::fread(.args[1])
if (class(serodata$Age) == "integer") {
  serodata[, lower := Age ]
  serodata[, upper := Age ]
} else {
  serodata[, lower := as.integer(gsub("(\\d+)-.+", "\\1", Age)) ]
  serodata[, upper := as.integer(gsub(".+-(\\d+)", "\\1", Age)) ]
}

# this  assumes that all of the data is grouped by Country (like Morrison 2010 and L'Azou 2016),
# rather than some other grouping
# if there are different headers, or no headers, substitute the correct columns into
# denvax::serofit invocation, by= group specification
result <- serodata[,
  denvax::serofit(Seropositive, Number, lower, upper),
  by = Country
]

cat(
  jsonlite::prettify(jsonlite::toJSON(result, auto_unbox = T)),
  file = tail(.args, 1)
)
