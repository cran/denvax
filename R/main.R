#' denvax: Simple Dengue Test and Vaccinate Cost Thresholds
#'
#' Provides the mathematical model described by "Serostatus Testing & Dengue Vaccine Cost-Benefit Thresholds"
#' in \href{https://doi.org/10.1098/rsif.2019.0234}{<doi:10.1098/rsif.2019.0234>}.  Using the functions in the package,
#' that analysis can be repeated using sample life histories, either synthesized from local seroprevalence data
#' using other functions in this package (as in the manuscript) or from some other source.
#' The package provides a vignette which walks through the analysis in the publication, as well as a function
#' to generate a project skeleton for such an analysis.
#'
#' @section \code{denvax} functions:
#' \describe{
#'   \item{\code{\link{serofit}}}{fits serosurvey data against two-risk model}
#'   \item{\code{\link{synthetic.pop}}}{using parameter fit data, generate a sample population}
#'   \item{\code{\link{nPxA}}}{estimate dengue infection outcome probabilities based on a synthetic population}
#'   \item{\code{\link{ROIcoeffs}}}{using population outcome probabilities, compute ROI equation coefficients}
#'   \item{\code{\link{ROI}}}{compute ROIs from setting coefficients and cost scenarios}
#'   \item{\code{\link{build.project}}}{create a template project for estimating ROIs}
#' }
#'
#' @docType package
#' @name denvax
NULL

# NOT EXPORTED: constructor, reductor, use.default, inv.logit, prob, llp

####################################################################
# functions that swap behaviors based on what packages are available

constructor <- function(...) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::data.table(...)
  } else data.frame(...)
}

reductor <- function(...) { if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::rbindlist(...)
} else Reduce(rbind, ...) }

# TODO: parallel apply for prune.exposure? Rcpp?

####################################################################
# utility functions

# use as with default function
use.default <- function(val, arg, note=warning) {
  note(sprintf("using default for %s...\n", arg))
  return(val)
}

####################################################################
# functions perform maximum likelihood related computations

# inverse logit
inv.logit <- function(l) 1 / (1 + exp(-l))

# probability of the two-risk model with logit scale arguments
prob <- function(A, lgtf_H, lgtRR, lgtp_H) {
  f_H <- inv.logit(lgtf_H)
  f_L <- f_H * inv.logit(lgtRR)
  p_H <- inv.logit(lgtp_H)
  return(p_H * (1 - (1 - f_H) ^ A) + (1 - p_H) * (1 - (1 - f_L) ^ A))
}

# full likelihood computation for all the measurements
llp <- function(pars, Amins, Amaxs, pos, N, pfunc) {
  avep <- mapply(
    function(mn, mx, ...) mean(pfunc(mn:mx, ...)),
    Amins, Amaxs, MoreArgs = as.list(pars)
  )
  if (any(avep <= 0) | any(avep >= 1)) {
    return(1000)
  } else return(-sum(stats::dbinom(pos, N, avep, log = T)))
}

#' Model fitting for serological data
#'
#' @param seropositive, the number of seropositive samples for each age group; length(seropositive) must be at least 3
#' @param N, the total number of samples for each age group; length(N) must equal length(seropositive)
#' @param age.min, the low age in age groups; defaults to `1:length(seropositive)`, i.e. assumes
#' the seropositive data corresponds to yearly cohorts starting at age 1.
#' @param age.max, the upper age in age groups; defaults to `age.min`, i.e. assumes each category corresponds
#' to a single year
#'
#' @details Fits a constant force of infection, two-risk category model using
#'   seroprevalence survey data. i.e.:
#'
#'   $$
#'   P_+(A) = p_H * (1-(1-f_H)^A) + (1-p_H) * (1-(1-f_L)^A)
#'   $$
#'
#'   This probability is fit to the seroprevalence by age category data,
#'   using maximum likelihood and \code{\link[stats]{optim}}.
#'
#' @return a list of best-fit parameters, all numeric values:
#' \describe{
#'   \item{f_H}{force of infection, for the high risk group}
#'   \item{f_L}{force of infection, for the low risk group}
#'   \item{p_H}{the proportion of the population at high risk}
#' }
#'
#' @examples
#' require(denvax);
#' data(morrison2010) # has counts by age
#' fit <- with(morrison2010, serofit(sero=Seropositive, N=Number, age.min=Age))
#' if (requireNamespace("data.table", quietly = TRUE)) {
#' data(lazou2016) # has counts by age range, instead of counts for every year
#' # this example uses `data.table`` functions to simplify processing
#' # several groups at once
#'   lazou2016[,{
#'     agerange <- data.table::tstrsplit(Age, "-")
#'     serofit(
#'       sero     = Seropositive,
#'       N        = Number,
#'       age.min  = as.integer(agerange[[1]]),
#'       age.max  = as.integer(agerange[[2]])
#'     )
#'   }, by = Country]
#' }
#' @export
serofit <- function(
  seropositive, N,
  age.min = use.default(1L:length(seropositive), "age.min"),
  age.max = use.default(age.min, "age.max")
) {
stopifnot(length(seropositive) >= 3, length(seropositive) == length(N))
# TODO adjust to take seropositive percents

  # fit model parameters in logit space
  pars <- c(lgtf_H = 0, lgtRR = 0, lgtp_H = 0)
  ret <- stats::optim(pars, llp,
    Amins = age.min, Amaxs = age.max,
    pos = seropositive, N = N,
    pfunc = prob,
    method = "L-BFGS-B" # selected because it does not require limits
  )
  newpar <- inv.logit(ret$par)
  newpar["lgtRR"] <- newpar["lgtRR"] * newpar["lgtf_H"]
  names(newpar) <- c("f_H", "f_L", "p_H")
  return(as.list(newpar))
}

# sample H vs L fois
rfoi <- function(n, f_H, f_L, p_H) sample(
  c(f_H, f_L), n, replace = T, prob = c(p_H, 1 - p_H)
)
# equivalent, but this version is slower once n large:
# rfoi.alt <- function(n, f_H, f_L, p_H) ifelse(runif(n) < p_H, f_H, f_L)

seros <- function(nSeries, toAge, popsize, ntypes=4L) matrix(
  sample.int(ntypes, size = toAge * nSeries, replace = T),
  nrow = popsize * nSeries,
  ncol = toAge,
  byrow = T
)

expfoi <- function(foi, nSeries, toAge) matrix(
  rep(foi, each = nSeries * toAge),
  ncol = toAge, byrow = T
)

#' @importFrom utils tail head
NULL
# TODO need to check handling of NAs
prune.exposures <- function(outcomes) t(apply(outcomes, 1, function(thisrow) {
  if (any(thisrow > 0, na.rm = T)) {
    tar <- which.max(thisrow > 0)
    maxage <- length(thisrow) # TODO accommodate NAs
    while (tar & (tar != maxage)) {
      thisrow[tar + 1] <- 0 # prevent any infection in next year
      tsero <- thisrow[tar] # get this infection type
      # remove all future infections of this type
      rmvs <- tail(which(thisrow == tsero), -1)
      thisrow[rmvs] <- 0
      # know 1 is excluded, since it's zero'd;
      # also tar != maxage, so no empty sets
      ntar <- which.max(thisrow[-(1:tar)] > 0)
      tar <- ifelse(ntar == 1, 0, tar + ntar) # update target
    }
  }
  thisrow
}))

#   Future: in addition to 0-4, `NA_integer` after end of individual series lifetime. Not
#   currently supported when generating series, but can be accommodated by the digesting
#   function in \link{nPxA}.

#' Compute synthetic population trajectories from model parameters
#'
#' @param pars, a list with elements `f_H`, `f_L`, and `p_H`; see \link{serofit} return value
#' @param runs, the number of different serotype timelines to simulate
#' @param popsize, the number of individuals to sample for each run
#' @param maxAge, the length of each lifetime
#' @param rngseed, an optional seed for the random number generator
#'
#' @details Using fitted parameters for a two-risk, constant force of infection model,
#'   simulate a dengue annual exposures model for the requested number of serotype series (`runs`)
#'   and individuals (`popsize`).  The resulting matrix is a collection of integers, 0-4.
#'   0 indicates no infection, 1-4 infection by the corresponding serotype.
#'
#' @return a matrix of integers 0-4, rows `runs*popsize` x columns `maxAge`
#'
#' @examples
#' require(denvax);
#' data(morrison2010) # has counts by age
#' fit <- with(morrison2010, serofit(sero=Seropositive, N=Number, age.min=Age))
#' m2010pop <- synthetic.pop(fit, runs = 10, popsize = 10) # small sample size for example run time
#' head(m2010pop)
#'
#' @export
synthetic.pop <- function(
  pars, runs=100, popsize=1000, maxAge=70, rngseed=NULL
) with(pars, {
  if (!is.null(rngseed)) set.seed(rngseed)

  # each individuals foi
  foi <- rfoi(popsize, f_H, f_L, p_H)
  # make all the seroseries;
  # pattern is circulations 1:n, 1:n, etc for each individual
  outcomes <- seros(runs, maxAge, popsize)

  # create a lifetime by replicating each individual foi to maxage,
  # then comparing random series
  exposure <- runif(popsize * maxAge * runs) < expfoi(foi, runs, maxAge)

  # everywhere exposures are NOT successful, set to exposing pathogen to 0
  outcomes[!exposure] <- 0
  # apply infection rules:
  #  - no re-infection with same sero, and
  #  - no infections immediately following an infection
  res <- prune.exposures(outcomes)
  dimnames(res) <- list(sample = NULL, age = NULL)
  return(res)
})

# helper function for nPxA
crunch <- function(
  qualifying_sumlh, popAtAge
) colSums(qualifying_sumlh, na.rm = T) / popAtAge

#' Compute the nPx(A), C(A) proportions from a population of life histories
#'
#' @param lifehistory, a matrix with rows (sample individuals) and columns (outcome in year of life); see \link{synthetic.pop} return value
#'
#' @details computes the relevant nPx(A) and C(A): the probabilities of the various life trajectories, by age.
#'   See \href{https://doi.org/10.1098/rsif.2019.0234}{<doi:10.1098/rsif.2019.0234>}, SI section II.A (Cost Benefit Equations: Definitions)
#'
#' @return a \code{\link[base]{data.frame}} (\code{\link[data.table]{data.table}}, if available) with columns
#' \describe{
#'   \item{A}{integer; the reference year of life, from 1 to \code{dim(lifehistory)[2]}}
#'   \item{p_0}{numeric; probability of 0 lifetime infections}
#'   \item{p_1p}{numeric; probability of 1 or more lifetime infections}
#'   \item{p_2p}{numeric; probability of 2 or more lifetime infections}
#'   \item{p0_1}{numeric; probability of 0 infections at age A, and 1 lifetime infection}
#'   \item{p0_1p}{numeric; probability of 0 infections at age A, and 1 or more lifetime infections}
#'   \item{p0_2p}{numeric; probability of 0 infections at age A, and 2 or more lifetime infections}
#'   \item{p1_A}{numeric; probability of 1 infection at age A, and 1 or more lifetime infections}
#'   \item{p1_2p}{numeric; probability of 1 infection at age A, and 2 or more lifetime infections}
#'   \item{p1p_A}{numeric; probability of 1 or more infections at age A, and 1 or more lifetime infections}
#'   \item{p2p_A}{numeric; probability of 2 or more infections at age A, and 2 or more lifetime infections}
#'   \item{CA}{numeric; probability of converting from seronegative to seropositive between age A and A+1}
#' }
#'
#' @examples
#' require(denvax);
#' data(morrison2010) # has counts by age
#' fit <- with(morrison2010, serofit(sero=Seropositive, N=Number, age.min=Age))
#' m2010pop <- synthetic.pop(fit, runs = 10, popsize = 10) # small sample size for example run time
#' m2010lh <- nPxA(m2010pop)
#' m2010lh
#' with(m2010lh,
#'   plot(A, p0_2p*100, type="l",
#'     xlab="Age", ylab="%", ylim = c(0, 100),
#'     main="Individuals w/ No Infections,\nbut that will have 2"
#'   )
#' )
#'
#' @export
nPxA <- function(lifehistory) {

  # total population - num of circulation series x num people exposed to them
  tpop <- colSums(!is.na(lifehistory))
  maxAge <- dim(lifehistory)[2]

  # make the by-age infection counts
  suminfs <- t(apply(lifehistory > 0, 1, cumsum))
  lifetimeinfs <- apply(suminfs, 1, max, na.rm = T)

  twoplus <- lifetimeinfs > 1
  oneplus <- lifetimeinfs > 0
  justone <- xor(oneplus, twoplus)
  zerothensome  <- crunch(suminfs[oneplus, , drop = F] == 0, tpop)
  zerothenone   <- crunch(suminfs[justone, , drop = F] == 0, tpop)
  zerothen2plus <- crunch(suminfs[twoplus, , drop = F] == 0, tpop)
  onethenmore   <- crunch(suminfs[twoplus, , drop = F] == 1, tpop)
  onenow        <- crunch(suminfs == 1, tpop)
  somenow       <- crunch(suminfs > 0,  tpop)
  twoormorenow  <- crunch(suminfs > 1,  tpop)
  l_none        <- crunch(!is.na(suminfs[!oneplus, , drop = F]), tpop)
  l_onep        <- crunch(!is.na(suminfs[oneplus, , drop = F]), tpop)
  l_twop        <- crunch(!is.na(suminfs[twoplus, , drop = F]), tpop)

  # TODO move to unit testing
  # average # of tests for starting at age 5, testing 11 times (i.e., to age 15)
  # mtests <- mean(apply(suminfs[, 5:15] != 0, 1, function(rw) ifelse(any(rw),which.max(rw),length(rw))))
  # should be identical to ave.tests calc

  # TODO move to unit testing
  # dagvfrac <- mean(apply(suminfs[, 5:15] != 0, 1, any))
  # privfrac <- mean(apply(suminfs[, 5:15] != 0, 1, any)) - mean(suminfs[,5]>1)
  # should be identical vac.frac calc
  # dag vfrac = 1 - unlist(probabilities[sl, "p0_1p"]+probabilities[sl, "p_0"])
  # pri vfrac = 1 - unlist(probabilities[sl, "p0_1p"]+probabilities[sl, "p_0"]) - unlist(probabilities[A, "p2p_A"])
  # vs old calc
  # coeffs <- unlist(probabilities[ssl, "p0_1p"])*unlist(probabilities[ssl, "CA"])
  # refv <- cumsum(coeffs)
  # baseres$Vmul <- c(0, refv)


  seroconvs <- t(apply(cbind(0, suminfs), 1, diff))[lifetimeinfs != 0,, drop = F]
  seroconvbyage <- rle(sort(apply(seroconvs, 1, which.max)) - 1)
  ref <- rep(0, maxAge)
  ref[seroconvbyage$values + 1] <- seroconvbyage$lengths
  totalseroconverted <- cumsum(ref)
  tot <- tail(totalseroconverted, 1)
  remainingtoconvert <- c(tot, head(tot - totalseroconverted, -1))
  seroconvA <- ref / remainingtoconvert
  seroconvA[is.na(seroconvA)] <- 0

  return(constructor(
    A     = 1L:maxAge,
    p_0   = l_none,
    p_1p  = l_onep,
    p_2p  = l_twop,
    p0_1  = zerothenone,
    p0_1p = zerothensome,
    p0_2p = zerothen2plus,
    p1_A  = onenow,
    p1_2p = onethenmore,
    p1p_A = somenow,
    p2p_A = twoormorenow,
    CA    = c(tail(seroconvA, -1), NaN)
  ))


}


#' Compute the Return on Investment (ROI) surface coefficients from population probabilities
#'
#' @param probabilities, a \code{\link[base]{data.frame}} (or \code{\link[data.table]{data.table}}) with the probabilities resulting from \link{nPxA}. Rows must correspond to ages, starting with age 1
#' @param As, the starting age(s) to consider
#' @param Ls, the maximum number of tests for each age; should either be an integer per age or a single integer for all ages.
#'   The default behavior computes the number of tests (for each age) that makes the maximum of `As` the maximum testing age
#'   Note: results will also be provided for shorter testing intervals, as the intermediate coefficients are calculated as part
#'   of computing the value at the maximum \code{L}
#'
#' @details computes the coefficients for the economic calculations
#' @return a \code{\link[base]{data.frame}} (\code{\link[data.table]{data.table}}, if available) with columns:
#' \describe{
#'   \item{A}{integer; the age when routine test-then-vaccinate strategy starts (from \code{As})}
#'   \item{L}{integer; the maximum number of tests for routine test-then-vaccinate strategy (from \code{Ls})}
#'   \item{vacfrac}{numeric; the fraction of individuals participating in this strategy that get vaccinated}
#'   \item{pri.offset}{numeric; the (additive) reduction in \code{vacfrac} if using the ordinal test}
#'   \item{Sfrac}{numeric; the proportion experiencing second infection costs}
#'   \item{Fresp}{numeric; the F/S cost fraction term, when comparing vaccination with and without testing}
#'   \item{Sgain}{numeric; the S term, when comparing vaccination with and without testing}
#' }
#'
#' @examples
#' require(denvax);
#' data(morrison2010) # has counts by age
#' fit <- with(morrison2010, serofit(sero=Seropositive, N=Number, age.min=Age))
#' m2010pop <- synthetic.pop(fit, runs = 10, popsize = 10) # small sample size for example run time
#' m2010lh <- nPxA(m2010pop)
#' rc <- ROIcoeffs(m2010lh, As=5:10, Ls=5)
#' pp <- par()
#' par(mfrow=c(1, 2))
#' rcs <- subset(rc, A==10 & L < 11)
#' with(rcs, plot(
#'   L, aveTests, type="l",
#'   xlab="Max # of Tests Allowed",
#'   ylab="Ave # of Tests Administered",
#'   main="Starting @ Age 10",
#'   ylim=c(1, 3)
#' ))
#' rcs <- subset(rc, A==5 & L < 11)
#' with(rcs, plot(
#'   L, aveTests, type="l",
#'   xlab="Max # of Tests Allowed",
#'   ylab="",
#'   main="Starting @ Age 5",
#'   ylim=c(1, 3)
#' ))
#' par(pp)
#'
#' @export
ROIcoeffs <- function(
  probabilities,
  As = 5:20,
  Ls = (diff(range(As))+1):1
) reductor(mapply(function(A,L) {
  Lref <- A+L-1 # the last age of testing
  sl <- A:Lref # the range of ages for the intervention
  baseres <- constructor(
    A=A, L=1:L,
    vacfrac = 1 - with(probabilities[sl,], p0_1p+p_0),
    pri.offset = probabilities[A,]$p2p_A,
    Sfrac = with(probabilities[A,], p_2p-p2p_A) - probabilities[sl,]$p0_2p,
    Fresp = probabilities[A,]$p0_1p,
    Sgain = probabilities[A,]$p0_1p - probabilities[sl,]$p0_2p,
    aveTests = 1
  )

  if (L!=1) { # update the average number of tests
    ssl <- A:(Lref-1)
    coeffs <- with(probabilities[ssl,], p0_1p*CA)
    refc <- cumsum(coeffs*(1:(L-1)))
    baseres$aveTests <- 1 + with(probabilities[Lref,], p_0+p0_1p)*((1:L)-1) + c(0, refc)
  }
  return(baseres)
}, A=As, L=Ls, SIMPLIFY = F))

# addresses check notes for CRAN; these variables are actually available from within(rcoeffs, ...)
aveTests <- vacfrac <- pri.offset <- Sfrac <- Fresp <- Sgain <- NULL

# TODO fix documentation
#' Compute the ROI surfaces given test and vaccine cost fractions.
#'
#' @param rcoeffs, a data.frame with the ROI surface coefficients from \link{ROIcoeffs}
#' @param nus, the series of normalized vaccine costs to use for ROI calcs
#' @param taus, the series of normalized test costs to use for ROI calcs
#'
#' @details tabulates ROI
#' @return a `data.frame` (`data.table`, if available) with columns:
#' \describe{
#'   \item{nu}{numeric, the normalized vaccine cost used}
#'   \item{tau}{numeric, the normalized test cost used}
#'   \item{mechanism}{character, either "ordinal" or "binary" corresponding to the type of test}
#'   \item{A}{integer; the age when routine test-then-vaccinate strategy starts (from \code{As})}
#'   \item{L}{integer; the maximum number of tests for routine test-then-vaccinate strategy (from \code{Ls})}
#'   \item{cost}{numeric; the intervention cost (as a fraction of second infection cost)}
#'   \item{benefit}{
#'     numeric; the difference in health outcome cost (as a fraction of second infection cost) minus `cost`;
#'     positive values indicate positive net benefit
#'   }
#'   \item{roi}{numeric; return on investment: `benefit` over `cost`}
#' }
#'
#' @examples
#' require(denvax);
#' data(morrison2010) # has counts by age
#' fit <- with(morrison2010, serofit(sero=Seropositive, N=Number, age.min=Age))
#' m2010pop <- synthetic.pop(fit, runs = 10, popsize = 10) # small sample size for example run time
#' m2010lh <- nPxA(m2010pop)
#' L <- 5
#' rc <- ROIcoeffs(m2010lh, As=5:10, Ls=L)
#' rois <- ROI(rc, nus = 0.5, taus = 0.01)
#' srois <- subset(rois, mechanism == "binary")
#' mrois <- matrix(srois$roi, nrow = L)
#' contour(x=unique(srois$L), y=unique(srois$A), z=mrois,
#'   xlab = "Max # of Tests", ylab = "Initial Age", main="ROI Contour"
#' )
#'
#' @export
ROI <- function(
  rcoeffs,
  nus=seq(0.3, 0.9, by = 0.05),
  taus = 10 ^ seq(-2, -0.5, by = 0.05)
) {
  grd <- expand.grid(tau = taus, nu = nus)
  res <- mapply(function(nu, tau) with(rcoeffs, {
    cost <- c(aveTests * tau + (vacfrac - pri.offset) * nu, aveTests * tau + vacfrac * nu)
    benefit <- Sfrac - cost
    roi <- benefit / cost
    mechanism <- c(rep("ordinal", length(aveTests)), rep("binary", length(aveTests)))
#    fd.net <- interceptf + nu # this associated with the with vs without testing gain
    constructor(
      nu = nu, tau = tau, mechanism = mechanism,
      A = rep(A, times=2), L = rep(L, times = 2),
      cost = cost, benefit = benefit, roi = roi
    )
  }), nu = grd$nu, tau = grd$tau, SIMPLIFY = F)
  ret <- reductor(res)
  return(ret[, c("nu", "tau", "mechanism", "A", "L", "cost", "benefit", "roi")])
}



#' Creates an ROI estimation project.
#'
#' @param targetdir, path to enclosing directory. If this directory does not exist, will attempt to create it (recursively)
#' @param overwrite, overwrite existing files corresponding to the project skeleton elements?
#' @param copy_pub, copy the `/pub` folder (which contains the analyses from the publication)
#'
#' @details This function sets up the skeleton of an analysis to go from seroprevalence data to the ROI estimation surface.
#'   That skeleton uses a series of separate scripts for each analytical step (fitting, simulation, analysis, and application),
#'   connected via the command line build tool \code{make}.  This approach allows clean substitution for various stages (e.g.,
#'   using a different model to generate life histories). The following files are created:
#' \describe{
#'   \item{Makefile}{the dependencies for various analysis stages}
#'   \item{README.md}{brief notes about project parts}
#'   \item{fit.R}{script for fitting seroprevalence data}
#'   \item{synthesize.R}{script for generating synthetic populations}
#'   \item{digest.R}{script for converting life histories into probability coefficients for ROI calculation}
#'   \item{simple.R}{a quick example start-to-finish analysis}
#' }
#'
#' @return logical, indicating error free creation of the project skeleton; there may still be other warnings
#'
#' @examples
#' require(denvax)
#' tardir <- tempdir() # replace with desired target
#' build.project(tardir)
#' list.files(tardir, recursive = TRUE)
#'
#' @export
build.project <- function(
  targetdir,
  overwrite = FALSE,
  copy_pub  = TRUE
) {
  if (missing(targetdir)) {
    warning("No target directory supplied.")
    return(FALSE)
  }
  if (!length(targetdir)) {
    warning("Did not understand targetdir argument; expecting a character vector (of length 1).")
    return(FALSE)
  }
  if (length(targetdir) != 1) warning("length(targetdir) > 1; only using first element.")
  if (!dir.exists(targetdir[1]) && !dir.create(targetdir[1], recursive = T)) {
    warning(sprintf("'%s' does not exist and could not be created.", targetdir[1]))
    return(FALSE)
  }
  srcdir   <- system.file("extdata", package = "denvax")
  srcfiles <- list.files(srcdir, pattern = "*.R", full.names = T, recursive = F, include.dirs = F)

  res <- file.copy(from      = srcfiles,
                   to        = targetdir[1],
                   overwrite = overwrite,
                   recursive = TRUE)

  if (copy_pub){
    srcdir <- system.file("extdata/pub", package = "denvax")
    srcfiles <-  list.files(srcdir, pattern = "*.R", full.names = T, recursive = F, include.dirs = F)
    targetdir <- paste(targetdir[1], "pub", sep = .Platform$file.sep)
    dir.create(targetdir)
    res <- c(res, file.copy(from      = srcfiles,
                     to        = targetdir,
                     overwrite = overwrite,
                     recursive = TRUE))
  }


  if (!all(res)) warning("Something may have gone wrong with copying.")
  return(any(res))
}
