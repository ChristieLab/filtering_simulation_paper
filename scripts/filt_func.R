flt_func <- function(d, maf = FALSE, mac = 0, analysis_facet, filter_facet = NULL, 
                     hwe, hwe_facet, min_ind, min_loci,
                     chr, subset_seed, step, sigma, par = FALSE, tdir,
                     ne_estimator_path){
  
  # check 0s/FALSEs
  
  if(is.numeric(maf) & maf == FALSE){
    maf <- FALSE
  }
  
  if(is.numeric(par) & par == FALSE){
    par <- FALSE
  }
  
  cat("Filtering.\n")
  f <- filter_snps(d, maf = maf, mac = mac, maf_facets = filter_facet, 
                   hwe = hwe, hwe_facets = hwe_facet, min_ind = min_ind,
                   min_loci = min_loci, verbose = TRUE)
  
  cat("Basic stats.\n")
  f <- calc_pi(f, analysis_facet)
  # f <- calc_he(f, analysis_facet)
  f <- calc_ho(f, analysis_facet)
  f <- calc_fis(f, analysis_facet)
  # f <- calc_prop_poly(f, analysis_facet)
  
  cat("Fst and/or Tajima's D.\n")
  if(analysis_facet != ".base"){
    f <- calc_pairwise_fst(f, analysis_facet)
    f <- calc_tajimas_d(f, paste0(analysis_facet, ".", chr), sigma, step, triple_sigma = FALSE, par = par, verbose = TRUE)
    set.seed(subset_seed)
    if(nrow(f) > 1000){
      ld <- calc_pairwise_ld(f, analysis_facet, ss = 1000, par = par, verbose = TRUE)
    }
    else{
      ld <- calc_pairwise_ld(f, analysis_facet, par = par, verbose = TRUE)
    }
    
  }
  else{
    f <- calc_tajimas_d(f, chr, sigma, step, triple_sigma = FALSE, par = par, verbose = TRUE)
    set.seed(subset_seed)
    if(nrow(f) > 1000){
      ld <- calc_pairwise_ld(f, ss = 1000, par = par, verbose = TRUE)
    }
    else{
      ld <- calc_pairwise_ld(f, par = par, verbose = TRUE)
    }
  }
  
  if(nrow(f) > 30000){
    set.seed(subset_seed)
    sub.f <- f[sample(nrow(f), 30000, replace = FALSE),]
  }
  else{
    sub.f <- f
  }
  
  owd <- getwd()
  setwd(tdir)
  
  cat("Ne.\n")
  sub.f <- calc_ne(sub.f, analysis_facet, chr = chr, pcrit = 0, NeEstimator_path = ne_estimator_path)
  
  setwd(owd)  
  
  cat("PCA.\n")
  if(analysis_facet != ".base"){
    pca <- plot_clusters(f, analysis_facet)
  }
  else{
    pca <- list(plots = list(pca = NULL))
  }
  
  cat("SFS.\n")
  if(analysis_facet != ".base"){
    options <- summarize_facets(f, analysis_facet)
    options <- unlist(options)
    sfs <- vector("list", length(options))
    names(sfs) <- options
    for(i in 1:length(options)){
      keep.samps <- which(sample.meta(f)[[analysis_facet]] == options[i])
      tf <- f[,keep.samps]
      tf <- filter_snps(f, non_poly = TRUE)
      if(nrow(tf) > 5 & ncol(tf) > 100){
        sfs[[i]] <- calc_sfs(tf, projection = nsamps(tf), fold = TRUE)
        sfs[[i]] <- data.frame(mac = 1:length(sfs[[i]]), count = sfs[[i]], maf = maf, mac = mac)
      }
      else{
        sfs[i] <- list(NULL)
      }
    }
    
    sfs <- data.table::rbindlist(sfs, idcol = "pop")
  }
  else{
    sfs <- calc_sfs(f, projection = nsamps(f), fold = TRUE)
    sfs <- data.frame(mac = 1:length(sfs), count = sfs, maf = maf, mac = mac)
  }
  
  
  cat("Done. Returning.\n")
  
  
  
  if(analysis_facet != ".base"){
    r1 <- get.snpR.stats(f, analysis_facet, c("pi", "ho", "fis", "fst"))$weighted.means
    r2 <- get.snpR.stats(f, paste0(analysis_facet, ".", chr), "tajimas_d")$weighted.means[
      which(get.snpR.stats(f, paste0(analysis_facet, ".", chr), "tajimas_d")$weighted.means$snp.subfacet == ".OVERALL_MEAN"),
    ]
    ld <- get.snpR.stats(ld, analysis_facet, "ld")$LD$prox
    ld <- ld[,c("proximity", "sample.subfacet", "sample.facet", "CLD")]
    ld$maf <- maf
    ld$mac <- mac
  }
  else{
    r1 <- get.snpR.stats(f, analysis_facet, c("pi", "ho", "fis"))$weighted.means
    r2 <- get.snpR.stats(f, chr, "tajimas_d")$weighted.means[
      which(get.snpR.stats(f, chr, "tajimas_d")$weighted.means$snp.subfacet == ".OVERALL_MEAN"),
    ]
    ld <- get.snpR.stats(ld, stats = "ld")$LD$prox
  }
  
  
  
  r <- merge(r1, r2, by = c("facet", "subfacet"), all = TRUE)
  r <- r[,-grep("\\.x$", colnames(r))]
  r <- r[,-grep("facet\\.y$", colnames(r))]
  colnames(r) <- gsub("\\.y$", "", colnames(r))
  
  r3 <- get.snpR.stats(sub.f, analysis_facet, "ne")$pop
  colnames(r3)[which(colnames(r3) == "pop")] <- "subfacet"
  r3 <- r3[,-grep("CI", colnames(r3))]
  r <- merge(r, r3, by = c("facet", "subfacet"), all = TRUE)
  r$snp.facet <- NULL
  r$snp.subfacet <- NULL
  r$nsnps <- nrow(f)
  r <- reshape2::melt(r, id.vars = c("facet", "subfacet"))
  
  r <- na.omit(r)
  r$maf <- maf
  r$mac <- mac
  
  return(list(stats = r,
              pca = pca$plots$pca,
              sfs = sfs,
              ld = ld))
}
