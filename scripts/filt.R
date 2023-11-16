library(snpR)
#=============parse command args================
args <- commandArgs(TRUE)
i <- as.numeric(args[1])
parmfile <- as.character(args[2])
outfile_prefix <- as.character(args[3])
ne_estimator_path <- as.character(args[4])
par <- as.numeric(args[5])
temp_path <- as.character(args[6])
if(par == 0){par <- FALSE}

#============parse parmfile args===============

parmfile <- data.table::fread(parmfile)

maf <- parmfile[i,1]
hwe <- parmfile[i,2]
min_ind <- parmfile[i,3]
min_loci <- parmfile[i,4]
LD <- parmfile[i,5]
LD_sigma <- parmfile[i,6]
LD_cut <- parmfile[i,7]
do_sel <- parmfile[i,8]
sel_win_size <- parmfile[i,9]
infile <- parmfile[i,10]

#================filter============
cat("Filtering.\n")

# read in and change to temp dir
f <- readRDS(infile)
cdir <- getwd()
tempdir <- tempfile("Rtemp", temp_path)
dir.create(tempdir)
setwd(tempdir)

# filter
f <- filter_snps(f, maf = maf, maf_facets = "pop", 
                 hwe = hwe, hwe_facets = "pop", min_ind = min_ind,
                 min_loci = min_loci, verbose = TRUE, LD_prune_sigma = ifelse(LD, LD_sigma, FALSE),
                 LD_prune_r = LD_cut, LD_prune_facet = "pop.chr", LD_prune_par = par)

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

if(do_sel){
  f <- calc_smoothed_averages(f, "pop.chr", sigma = sel_win_size, triple_sigma = FALSE, gaussian = FALSE, par = par)
}

#===============graphics============
cat("Generating pca/sfs.\n")
pca <- plot_clusters(f, "pop")


cat("SFS.\n")
upops <- summarize_facets(f, "pop")$pop
upops <- names(upops)
sfs <- data.table(pop = numeric(), num_min_alleles = numeric(), count = numeric(), proj_prop = numeric(), proj = numeric())

# make an sfs for each pop, picking the best projection from 50% to 90% of the
# max possible number of gene copies. Checking 50, 70, and 90%
for(i in 1:length(upops)){
  proj <- summarize_facets(f, "pop")[[1]]
  proj <- floor(proj[match(upops[i], names(proj))]*c(.5, .7, .9))*2
  
  for(j in 1:length(proj)){
    snpR:::.make_it_quiet(tsfs <- calc_sfs(f, "pop", upops[i], proj[j], fold = TRUE))
    sfs <- rbind(sfs, data.table(pop = upops[i], num_min_alleles = 1:length(tsfs), count = tsfs, proj_prop = c(.5, .7, .9)[j],
                 proj = proj[j]))
  }
}

#==============save and return
setwd(cdir)

cat("Done. Returning.\n")
# set the outfile prefix
outfile_prefix <- paste0(outfile_prefix, "_r", i, "_")

# get the filt summary tab
filt_sum <- data.table(maf = maf, mac = mac, hwe = hwe, min_ind = min_ind, min_loci = min_loci, LD = LD)

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
  fst <- get.snpR.stats(f, "pop", "fst")$weighted.means
  fst$snp.facet <- fst$facet <- fst$snp.subfacet <- NULL
  fst <- cbind(fst, filt_sum)
  data.table::fwrite(fst, paste0(outfile_prefix, "fst.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, scipen = 999)
}

# sfs
sfs <- cbind(sfs, filt_sum)
data.table::fwrite(sfs, paste0(outfile_prefix, "sfs.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, scipen = 999)

# pca
pca <- cbind(pca$pca$data[,c("pop", "PC1", "PC2")], filt_sum)
data.table::fwrite(pca, paste0(outfile_prefix, "pca.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, scipen = 999)
