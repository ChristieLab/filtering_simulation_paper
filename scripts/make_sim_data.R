library(coala); library(foreach); library(doParallel); library(parallel)
# comp <- data.table::fread("../aus_monarchs/data/genotypes.geno")

lens <- read.table("data/mon_chr_lengths.txt")


coala::activate_scrm()

mu <- 2.9e-9
N0 <- 1e4
m <- .05

ntasks <- 30

cl <- parallel::makePSOCKcluster(ntasks)
doParallel::registerDoParallel(cl)

#loop through each set
out <- foreach::foreach(q = 1:ntasks) %dopar% {
  tlen <- lens$V2[q]
  d <- coal_model(c(100, 100), loci_number = 1, loci_length = tlen, ploidy = 2) +
    sumstat_seg_sites() +
    feat_mutation(4*N0*mu*tlen) + 
    feat_migration(4*N0*m, symmetric = TRUE) + 
    feat_size_change(.1, population = 2, time = .3) +
    feat_size_change(1, population = 2, time = .5) +
    feat_recombination(4*N0, locus_group = 1)
  
  res <- simulate(d)
  res$seg_sites
}

saveRDS(out, "mig_sim_dat.RDS")
