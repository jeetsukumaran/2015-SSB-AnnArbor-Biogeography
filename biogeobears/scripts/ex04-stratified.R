#!/usr/bin/env Rscript

## Load BioGeoBEARS and apply patches
source("scripts/biogeobears-utilities.R")
load.biogeobears()

###################################################################
## Set up the phylogeny and geography
phylogeny.path = "data/Psychotria_5.2.newick"
psychotria.tree = read.tree(phylogeny.path)
geography.path = "data/Psychotria_geog.data"
disp.strat.run = define_BioGeoBEARS_run()
disp.strat.run$trfn = phylogeny.path
disp.strat.run$geogfn = geography.path


###################################################################
## Set up the stratification
disp.strat.run$timesfn = "data/timeperiods.txt"
disp.strat.run$dispersal_multipliers_fn = "data/Psychotria-stratified-dispersal.txt"
# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
disp.strat.run = readfiles_BioGeoBEARS_run(disp.strat.run)
# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
disp.strat.run = section_the_tree(
        inputs=disp.strat.run,
        make_master_table=TRUE,
        plot_pieces=FALSE)
# The stratified tree is described in this table:
disp.strat.run$master_table

# Configure the run
disp.strat.run = configure.standard.biogeobears.run(disp.strat.run)

## Check
check_BioGeoBEARS_run(disp.strat.run)

## Execute the analysis
disp.strat.results = bears_optim_run(disp.strat.run)
save(disp.strat.results, file="results/Psychotria.DEC.disp.stratified.Rdata")
