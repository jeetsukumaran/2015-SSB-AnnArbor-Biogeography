#!/usr/bin/env Rscript

## Load BioGeoBEARS and apply patches
source("scripts/biogeobears-utilities.R")
load.biogeobears()

###################################################################
## Set up the phylogeny and geography
phylogeny.path = "data/Psychotria_5.2.newick"
psychotria.tree = read.tree(phylogeny.path)
geography.path = "data/Psychotria_geog.data"
dec.aa.run = define_BioGeoBEARS_run()
dec.aa.run$trfn = phylogeny.path
dec.aa.run$geogfn = geography.path


###################################################################
## Set up the stratification
dec.aa.run$timesfn = "data/timeperiods.txt"
dec.aa.run$areas_allowed_fn = "data/Psychotria-areas-allowed.txt"
dec.aa.run = readfiles_BioGeoBEARS_run(dec.aa.run)
# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
dec.aa.run = section_the_tree(
        inputs=dec.aa.run,
        make_master_table=TRUE,
        plot_pieces=FALSE)
# The stratified tree is described in this table:
dec.aa.run$master_table

# Configure the run
dec.aa.run = configure.standard.biogeobears.run(dec.aa.run)

## Check
check_BioGeoBEARS_run(dec.aa.run)

## Execute the analysis
dec.aa.results = bears_optim_run(dec.aa.run)
save(dec.aa.results, file="results/Psychotria.DEC.dec.aa.Rdata")
