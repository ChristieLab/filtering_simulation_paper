## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------
library(snpR); library(data.table); library(coala);
# remotes::install_github("hemstrow/snpR")
args <- commandArgs(TRUE)
parmfile <- data.table::fread(as.character(args[1]), header = F)
iter <- as.numeric(args[2])
msms_path <- as.character(args[3])
java_path <- as.character(args[4])

parmfile <- unlist(parmfile[iter,])
sf <- as.numeric(parmfile[1])
ts <- as.numeric(parmfile[2])
strn <- as.numeric(parmfile[3])
r <- as.numeric(parmfile[4])
out <- as.character(parmfile[5])

## -----------------------------------------------------------------------------------------------------------------------------
N0 <- 10000 # starting population size (at present)
mu <- 1e-8 # per-base mutation rate
len <- 10000000 # chromosome length
r <- r # average number of recombination events per chr
g <- 1 # years per gen
t1 <- ts # time of ingroup split
t2 <- 1000 # time of outgroup split
samples <- 50 # sample size
time_selection <- ts # time where selection starts
starting_selection_freq <- sf # starting frequency of allele under selection; hard 1/(2N), soft .05
selection_strength <- 2*N0*(strn) # need to adjust after author feedback

output_file_name <- out # no file extension


## -----------------------------------------------------------------------------------------------------------------------------
# for the one chr under selection
model_basic_sel <- coal_model(sample_size = c(samples, samples, samples), loci_number = 1, ploidy = 2, loci_length = len) +
  feat_pop_merge(par_expr(t1/(4*N0*g)), 3, 2) + #pop. 3 merges into pop. 2, 0.5 coal. units in the past
  feat_pop_merge(par_expr(t2/(4*N0*g)), 2, 1) + #pops. 2+3 merge into pop. 1, 0.3 coal. units before 2+3 divergence
  feat_mutation(rate = par_expr(4*N0*mu*len)) + # 4*N*mu, where mu is mut. rate per locus
  feat_selection(strength_A = selection_strength, 
                 time = par_expr(time_selection/(4*N0*g)), 
                 population = 3, 
                 start = TRUE,
                 start_frequency = starting_selection_freq,
                 Ne = N0,
                 position = 0.25,
                 force_keep = TRUE,
                 locus_group = "all") +
  feat_recombination(rate = par_expr(4*N0*r)) +
  sumstat_seg_sites() + # generate segregating snps
  par_named("N0") +
  par_named("mu") +
  par_named("r") +
  par_named("g") +
  par_named("t1") +
  par_named("t2") +
  par_named("time_selection")


activate_msms(msms_path, java_path)
sim <- simulate(model_basic_sel, seed = 123 + iter,  pars = c(N0 = N0, len = len,
                                                              r = r, g = g, t1 = t1, t2 = t2,  
                                                              mu = mu,
                                                              time_selection = time_selection))


## -----------------------------------------------------------------------------------------------------------------------------
saveRDS(sim, paste0(output_file_name, ".RDS"))
