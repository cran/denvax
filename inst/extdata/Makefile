#########################################################
# MAKEFILE ORGANIZATION
# 1. DEFINITIONS
# 2. GENERIC PROJECT TARGETS
#########################################################

R = Rscript $^ $@

# default locations for inputs and outputs
OUTPUTDIR ?= results
INPUTDIR ?= inputs

dirs: $(OUTPUTDIR) $(INPUTDIR)

$(OUTPUTDIR) $(INPUTDIR):
	mkdir -p $@

.PHONY: cleanfits cleansims

cleanfits: | $(OUTPUTDIR)
	rm -f $|/fit-*.json

cleanlh: | $(OUTPUTDIR)
	rm -f $|/lh-*.rds

cleannpxa: | $(OUTPUTDIR)
	rm -f $|/npxa-*.rds

cleanreports: | $(OUTPUTDIR)
	rm -f $|/*.pdf

.PRECIOUS: $(OUTPUTDIR)/fit-%.json $(OUTPUTDIR)/lh-%.rds $(OUTPUTDIR)/npxa-%.rds

#########################################################



############## GENERIC PROJECT TARGETS ##################

# if using age-seroprevalence data
$(OUTPUTDIR)/fit-%.json: fit.R $(INPUTDIR)/serodata-%.csv | $(OUTPUTDIR)
	$(R)

# create life histories with two-risk constant foi model;
# this step could be replaced by an alternative model for creating life histories
$(OUTPUTDIR)/lh-%.rds: synthesize.R $(OUTPUTDIR)/fit-%.json | $(OUTPUTDIR)
	$(R)

# create probability distribution coefficients
$(OUTPUTDIR)/npxa-%.rds: digest.R $(OUTPUTDIR)/lh-%.rds | $(OUTPUTDIR)
	$(R)

