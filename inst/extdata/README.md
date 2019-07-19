# Project Template

The `denvax::build.project(...)` function creates a skeleton project your analysis.

The `Makefile` provides the relationships between the various single step scripts, and a way to run them as a workflow. Alternatively, they can be run manually from the command line or in an application like Rstudio.

The `simple.R` script demonstrates a simple start-to-finish analysis (using the Morrison 2010 data included in the package).

In either approach, the steps are:

1. develop fit(s) for local serological data (e.g., use the `fit.R` script in Rstudio, write your own script with the `denvax::serofit` function, or use make with the target `%-fit.json` in Makefile).
2. use a fit to synthesize a population (e.g., `synthesize.R` script, the `denvax::synthetic.pop` function, or the target `%-lh.rds` in Makefile)
3. use the synthetic population to estimate the lifetime outcome probabilities (e.g., `digest.R` script, the `denvax::npxa` function, or the target `%-npxa.rds` in Makefile)
4. combine the probabilities with the desired boundaries for interventions (initial testing age and maximum number of tests) and to compute the ROI equation coefficients (`denvax::ROIcoeffs` function)
5. finally, combine the coefficients with desired cost fraction ranges (i.e., the cost of vaccination and testing as fractions of the estimated burden of secondary infection, or $\nu$ and $\tau$ respectively), to determine the threshold region for positive ROI (`denvax::ROI` function).

The `pub/` folder includes all the scripts necessary to generate the figures associated with the original publication of the analysis.
