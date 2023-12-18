library(snpR)
#=============parse command args================
args <- commandArgs(TRUE)
iter <- as.numeric(args[1])
infile <- as.character(args[2])
parmfile <- as.character(args[3])
outfile_prefix <- as.character(args[4])
return_all_fst <- as.logical(args[5])
return_sfs <- as.logical(args[6])
ne_estimator_path <- as.character(args[7])
par <- as.numeric(args[8])
temp_path <- as.character(args[9])
if(par == 0){par <- FALSE}

#============parse parmfile args===============

parmfile <- data.table::fread(parmfile, header = TRUE)
parmfile <- as.data.frame(parmfile)

i <- iter
maf <- as.numeric(parmfile[i,1])
mgc <- as.numeric(parmfile[i,2])
hwe <- as.numeric(parmfile[i,3])
min_ind <- as.numeric(parmfile[i,4])
min_loci <- as.numeric(parmfile[i,5])
LD <- as.logical(parmfile[i,6])
LD_sigma <- as.numeric(parmfile[i,7])
LD_cut <- as.numeric(parmfile[i,8])

#================filter============
cat("Filtering.\n")

# read in and change to temp dir
f <- readRDS(infile)
cdir <- getwd()
tempdir <- tempfile("Rtemp", temp_path)
dir.create(tempdir)
setwd(tempdir)

# filter
if(maf > 0 | mgc > 0 | hwe > 0 | min_ind > 0 | min_loci > 0 | LD){
  f <- filter_snps(f, maf = maf, mgc = mgc, maf_facets = "pop", 
                   hwe = hwe, hwe_facets = "pop", min_ind = min_ind,
                   min_loci = min_loci, verbose = TRUE, LD_prune_sigma = ifelse(LD, LD_sigma, FALSE),
                   LD_prune_r = LD_cut, LD_prune_facet = "pop.chr", LD_prune_par = par)
}

#===============statistics=========
cat("Calculating statistics.\n")
f <- calc_he(f, "pop")
f <- calc_ho(f, "pop")
f <- calc_pi(f, "pop")
if(length(summarize_facets(f, "pop")$pop) > 1){
  f <- calc_pairwise_fst(f, "pop")
}
f <- calc_fis(f, "pop")
f <- calc_tajimas_d(f, "pop", global = TRUE)
f <- calc_seg_sites(f, "pop")
f <- calc_private(f, "pop")

if(nrow(f) >= 30000){
  f <- calc_ne(f, "pop", chr = "chr", pcrit = 0, nsnps = 30000, NeEstimator_path = ne_estimator_path)
}
if(nrow(f) < 30000){
  f <- calc_ne(f, "pop", chr = "chr", pcrit = 0, NeEstimator_path = ne_estimator_path)
}

#===============graphics============
cat("Generating pca/sfs.\n")
pca <- plot_clusters(f, "pop", simplify_output = TRUE)


cat("SFS.\n")
upops <- summarize_facets(f, "pop")$pop
upops <- names(upops)
sfs <- data.table(pop = character(), num_min_alleles = numeric(), count = numeric(), proj_prop = numeric(), proj = numeric())


if(return_sfs){
 # make an sfs for each pop, picking the best projection from 50% to 90% of the
  # max possible number of gene copies. Checking 50, 70, and 90%
  for(j in 1:length(upops)){
    proj <- summarize_facets(f, "pop")[[1]]
    proj <- floor(proj[match(upops[j], names(proj))]*c(.5, .7, .9))*2
    
    for(k in 1:length(proj)){
      snpR:::.make_it_quiet(tsfs <- try(calc_sfs(f, "pop", upops[j], proj[k], fold = TRUE), silent = TRUE))
      if(!is(tsfs, "try-error")){
        sfs <- rbind(sfs, data.table(pop = upops[j], num_min_alleles = 1:length(tsfs), count = tsfs, proj_prop = c(.5, .7, .9)[k],
                     proj = proj[k]))
      }
    }
  }
}

#==============save and return============
setwd(cdir)

cat("Done. Returning.\n")
# set the outfile prefix
outfile_prefix <- paste0(outfile_prefix, "_r", iter, "_")
cat("Outfile prefix: ", outfile_prefix, "\n")

# get the filt summary tab
filt_sum <- data.table(maf = maf, mgc = mgc, hwe = hwe, min_ind = min_ind, min_loci = min_loci, LD = LD)

# basic stats
r1 <- get.snpR.stats(f, "pop", c("pi", "ho", "he", "seg_sites", "fis", "private", "ne", "tajimas_d"))
ne <- r1$pop
r1 <- r1$weighted.means

r1 <- merge(r1, ne[c("pop", "LDNe_0")], by.x = "subfacet", by.y = "pop")
r1$facet <- r1$snp.facet <- r1$snp.subfacet <- r1$global_num_seg <- NULL
r1 <- cbind(r1, filt_sum)

data.table::fwrite(r1, paste0(outfile_prefix, "stats.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, scipen = 999)

# fst
if(length(summarize_facets(f, "pop")$pop) > 1){
  fst <- get.snpR.stats(f, "pop", "fst")
  if(return_all_fst){
    fst_all <- cbind(fst$pairwise, filt_sum)
    fst_all$facet <- NULL
    data.table::fwrite(fst_all, paste0(outfile_prefix, "fst_all.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, scipen = 999)
  }
  fst <- fst$weighted.means
  fst$snp.facet <- fst$facet <- fst$snp.subfacet <- NULL
  fst <- cbind(fst, filt_sum)
  data.table::fwrite(fst, paste0(outfile_prefix, "fst.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, scipen = 999)
}

# sfs
if(return_sfs){
  sfs <- cbind(sfs, filt_sum)
  data.table::fwrite(sfs, paste0(outfile_prefix, "sfs.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, scipen = 999)
}

# pca
pca <- cbind(pca$pca$data[,c("pop", "PC1", "PC2")], filt_sum)
data.table::fwrite(pca, paste0(outfile_prefix, "pca.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, scipen = 999)

