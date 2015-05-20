#!/usr/bin/env Rscript

## Load BioGeoBEARS and apply patches
source("scripts/biogeogbears-utilities.R")
load.biogeobears()

## Create the run object for the DEC analysis
dec.run = define_BioGeoBEARS_run()

## Inspect the phylogeny
phylogeny.path = "data/Psychotria_5.2.newick"
moref(phylogeny.path)
psychotria.tree = read.tree(phylogeny.path)
plot(psychotria.tree)
prt(psychotria.tree, get_tipnames=TRUE)

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

## Check
check_BioGeoBEARS_run(dec.run)

## Execute the analysis
dec.results = bears_optim_run(dec.run)
save(dec.results, file="results/Psychotria.DEC.Rdata")

## View the results by ranges
plot.biogeobears.results.ranges(
        results.object=dec.results,
        plot.type="text",
        analysis.title="DEC"
        )

## View the results by areas
plot.biogeobears.results.areas(
        results.object=dec.results,
        analysis.title="DEC"
        )

## Get the results as a table
dec.range.table = get.biogeobears.results.by.range.table(dec.results)
write.table(dec.range.table, "results/Psychotria.DEC.ranges.tsv", sep="\t", row.names=F)

## Get the area results as a table
dec.area.table = get.biogeobears.results.by.area.table(dec.results)
write.table(dec.area.table, "results/Psychotria.DEC.areas.tsv", sep="\t", row.names=F)

## Look up a node
ndi = getMRCA(psychotria.tree, c("P_hawaiiensis_Makaopuhi","P_wawraeDL7428"))
dec.range.table[ndi,]
