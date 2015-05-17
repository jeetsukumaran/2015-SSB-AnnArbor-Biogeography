#!/usr/bin/env Rscript

## Load BioGeoBEARS and apply patches
source("scripts/load_biogeobears.R")

## Create the run object for the DEC analysis
dec.run = define_BioGeoBEARS_run()

## Inspect the phylogeny
phylogeny.path = "data/Psychotria_5.2.newick"
moref(phylogeny.path)
psychotria.tree = read.tree(phylogeny.path)
plot(psychotria.tree)

## Set the phylogeny
dec.run$trfn = phylogeny.path
dec.run$trfn

## Inspect the geography
geography.path = "data/Psychotria_geog.data"
moref(geography.path)
getranges_from_LagrangePHYLIP(lgdata_fn=geography.path)

## Set the geography
dec.run$geogfn = geography.path
dec.run$geogfn

## Set a standard run configuration
## Note, due to the way R passes objects, you need to
## rebind the name "dec.run" to the returned object
dec.run = configure.standard.biogeobears.run(dec.run)