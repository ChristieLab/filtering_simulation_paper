library(coala); library(foreach); library(doParallel); library(parallel)
# comp <- data.table::fread("../aus_monarchs/data/genotypes.geno")

lens <- read.table("/scratch/bell/whemstro/filtering_simulation_paper/data/mon_chr_lengths.txt")


coala::activate_scrm()

mu <- 2.9e-9
N0 <- 1e4
m <- 10/N0
bottle <- 0.01

ntasks <- 30

cl <- parallel::makePSOCKcluster(ntasks)
doParallel::registerDoParallel(cl)

#loop through each set
out <- foreach::foreach(q = 1:ntasks, .packages = "coala") %dopar% {
  tlen <- lens$V2[q]
  d <- coal_model(c(100, 100), loci_number = 1, loci_length = tlen, ploidy = 2) +
    sumstat_seg_sites() +
    feat_mutation(4*N0*mu*tlen) + 
    feat_migration(4*N0*m, symmetric = TRUE) +
    feat_pop_merge(.8, pop_source = 2, 1) +
    feat_size_change(1, population = 2, time = .02) +
    feat_size_change(bottle, population = 2, time = .01) +
    feat_size_change(1, population = 2, time = 0) +
    feat_recombination(4*N0, locus_group = 1)
  
  res <- simulate(d)
  res$seg_sites
}

saveRDS(out, "/scratch/bell/whemstro/filtering_simulation_paper/data/raw_mig_sim_dat.RDS")
