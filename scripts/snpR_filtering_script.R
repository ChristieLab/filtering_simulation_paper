# run this also for the perfect data for the proposal
library(snpR); library(ggplot2); source("/home/hemstrow/filtering_simulation_paper/scripts/filt_func.R")

args <- commandArgs(TRUE)
i <- as.numeric(args[1])
tdir <- as.character(args[2])
parmfile <- as.character(args[3])
outfile <- as.character(args[4])
ne_estimator_path <- as.character(args[5])

d <- readRDS("../data/monarch_nomaf.RDS")
d <- d[pop = c("NAM", "HAW", "GUA", "ROT", "SAM", "NOR")]

# read in the parameter file:
## 1: filter facet
## 2: analysis facet
## 3: chr
## 4: maf
## 5: mac
## 6: hwe
## 7: hwe_facet
## 8: min_ind
## 9: min_loci
## 10: seed
## 11: par
## 12: step
## 13: sigma
parmfile <- read.table(parmfile, header = FALSE)[i,]

filter_facet <- as.character(parmfile[1])
analysis_facet <-  as.character(parmfile[2])
chr <-  as.character(parmfile[3])
subset_seed <- as.numeric(parmfile[10])

par <- parmfile[11]

results <- flt_func(d = d,
                     maf = as.numeric(parmfile[4]),
                     mac = as.numeric(parmfile[5]),
                     analysis_facet = analysis_facet,
                     filter_facet = filter_facet,
                     hwe = as.numeric(parmfile[6]),
                     hwe_facet = as.character(parmfile[7]),
                     min_ind = as.numeric(parmfile[8]),
                     min_loci = as.numeric(parmfile[9]),
                     chr = chr,
                     subset_seed = subset_seed,
                     step = as.numeric(parmfile[12]),
                     sigma = as.numeric(parmfile[13]),
                     par = as.numeric(par),
                     tdir = tdir,
                     ne_estimator_path = ne_estimator_path)

saveRDS(list(res = results, parms = parmfile), outfile)
