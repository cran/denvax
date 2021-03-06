include ../Makefile

INDIR = ../$(INPUTDIR)
OUTDIR = ../$(OUTPUTDIR)

################## PUBLICATION CODE #####################

# example targets for conducting the published analysis steps

# make reference csv (serodata-morrison2010.csv and serodata-lazou2016.csv)
$(INDIR)/serodata-%.csv: | $(INDIR)
	R -e "data.table::fwrite(denvax::$*, '$@')"

.PRECIOUS: $(INDIR)/serodata-morrison2010.csv $(INDIR)/serodata-lazou2016.csv

# compute nPxA for all groups
$(OUTDIR)/proportions.rds: cross_group_compare.R $(OUTDIR)/fit-morrison2010.json $(OUTDIR)/fit-lazou2016.json | $(OUTDIR)
	$(R)

# compute ROIs for Malaysia & Peru
$(OUTDIR)/example-rois.rds: pick_rois.R $(OUTDIR)/proportions.rds malaysia_v_peru.json | $(OUTDIR)
	$(R)

# compute bounding parameters based on serodata fits
$(OUTDIR)/scan-params.rds: cross_group_range.R $(OUTDIR)/fit-morrison2010.json $(OUTDIR)/fit-lazou2016.json | $(OUTDIR)
	$(R)

# compute ROIs for bounding parameters
$(OUTDIR)/scan-rois.rds: scan_rois.R $(OUTDIR)/scan-params.rds | $(OUTDIR)
	$(R)

$(addprefix $(OUTDIR)/fig-compare.,png pdf tiff jpg): fig-comparison.R $(addprefix $(OUTDIR)/,fig-scales.rda example-rois.rds) malaysia_v_peru.json | $(OUTDIR)
	$(R)

$(addprefix $(OUTDIR)/fig-scan.,png pdf tiff jpg): fig-scan.R $(addprefix $(OUTDIR)/,fig-scales.rda scan-rois.rds scan-params.rds) | $(OUTDIR)
	$(R)

$(addprefix $(OUTDIR)/fig-ntests.,png pdf tiff jpg): fig-ntests.R $(addprefix $(OUTDIR)/,fig-scales.rda scan-rois.rds scan-params.rds) | $(OUTDIR)
	$(R)

$(OUTDIR)/fig-scales.rda: fig-scales.R | $(OUTDIR)
	$(R)

.PHONY: manuscript

manuscript: $(addprefix $(OUTDIR)/,fig-compare.pdf fig-scan.pdf fig-ntests.pdf)

cleanmanuscript:
	cd $(OUTDIR); rm -f *.rds *.json *.rda
	cd $(INDIR); rm -f serodata-lazou2016.csv serodata-morrison2010.csv
