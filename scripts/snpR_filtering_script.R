# run this also for the perfect data for the proposal

library(snpR); library(ggplot2);

d <- readRDS("data/simulated_coalescent.RDS")

filter_facet <- ".base"
analysis_facet <- ".base"
chr <- "group"

mafs <- seq(0, .1, by = .01)
macs <- c(1, 3, 5)
hwe <- 1e-6
hwe_facet <- ".base"
step <- 50
sigma <- 100
subset_seed <- 1234

results <- vector("list", length = length(c(mafs, macs)))

prog <- 1

ufs <- unique(c(filter_facet, analysis_facet, hwe_facet))
ufs <- ufs[-grep("\\.base", ufs)]

# pre-compute maf for efficiency at filter facets (doesn't need to re-run each time the filtering happens)
if(length(ufs) != 0){
  d <- calc_maf(d, ufs)
}

flt_func <- function(d, maf = FALSE, mac = 0, analysis_facet, filter_facet = NULL, hwe, hwe_facet, chr, subset_seed, step, sigma, par = FALSE){
  
  cat("Filtering.\n")
  f <- filter_snps(d, maf = maf, mac = mac, maf_facets = filter_facet, hwe = hwe, hwe_facets = hwe_facet, verbose = FALSE)
  
  cat("Basic stats.\n")
  f <- calc_pi(f, analysis_facet)
  f <- calc_he(f, analysis_facet)
  f <- calc_ho(f, analysis_facet)
  f <- calc_fis(f, analysis_facet)
  f <- calc_prop_poly(f, analysis_facet)
  
  cat("Fst and/or Tajima's D.\n")
  if(analysis_facet != ".base"){
    f <- calc_pairwise_fst(f, analysis_facet)
    f <- calc_tajimas_d(f, paste0(analysis_facet, ".", chr), sigma, step, triple_sigma = FALSE)
  }
  else{
    f <- calc_tajimas_d(f, chr, sigma, step, triple_sigma = FALSE, par = par)
  }
  
  if(nrow(f) > 3000){
    set.seed(subset_seed)
    sub.f <- f[sample(nrow(f), 3000, replace = FALSE),]
  }
  else{
    sub.f <- f
  }
  
  cat("Ne.\n")
  sub.f <- calc_ne(sub.f, analysis_facet, chr = chr, pcrit = 0)
  
  cat("PCA.\n")
  if(analysis_facet != ".base"){
    pca <- plot_clusters(f, analysis_facet)
  }
  else{
    pca <- list(plots = list(pca = NULL))
  }
  
  cat("SFS.\n")
  sfs <- calc_sfs(f, projection = nsamps(f), fold = TRUE)
  sfs <- data.frame(mac = 1:length(sfs), count = sfs, maf = maf, mac = mac)
  
  cat("Done. Returning.\n")
  
  
 
  if(analysis_facet != ".base"){
    r1 <- get.snpR.stats(f, analysis_facet, c("pi", "he", "ho", "fis", "fst", "prop_poly"))$weighted.means
    r2 <- get.snpR.stats(f, paste0(analysis_facet, ".", chr), "tajimas_d")$weighted.means[
      which(get.snpR.stats(f, paste0(analysis_facet, ".", chr), "tajimas_d")$weighted.means$snp.subfacet == ".OVERALL_MEAN"),
    ]
  }
  else{
    r1 <- get.snpR.stats(f, analysis_facet, c("pi", "he", "ho", "fis", "prop_poly"))$weighted.means
    r2 <- get.snpR.stats(f, chr, "tajimas_d")$weighted.means[
      which(get.snpR.stats(f, chr, "tajimas_d")$weighted.means$snp.subfacet == ".OVERALL_MEAN"),
    ]
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
              sfs = sfs))
}

for(i in 1:length(mafs)){
  print(i)
  results[[prog]] <- flt_func(d, mafs[i], mac = 0, analysis_facet, filter_facet, hwe, hwe_facet, chr, subset_seed, step, sigma)
  prog <- prog + 1
}

for(i in 1:length(macs)){
  results[[prog]] <- flt_func(d, maf = FALSE, mac = macs[i], analysis_facet, filter_facet, hwe, hwe_facet, chr, subset_seed, step, sigma)
  prog <- prog + 1
}



maf_res <- results[1:length(mafs)]
mac_res <- results[(length(mafs) + 1):length(results)]
rm(results)

bs <- purrr::map(maf_res, "stats")
bs <- dplyr::bind_rows(bs)
ggplot(bs, aes(x = maf, y = value)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y")


theta_scew <- bs[grep("theta", bs$variable),]
theta_scew <- as.data.table(theta_scew)
theta_scew <- data.table::dcast(theta_scew, facet + subfacet + maf + mac ~ variable, value.var = "value")
theta_scew[,scew := (weighted_mean_ts.theta - weighted_mean_ws.theta)/((weighted_mean_ts.theta + weighted_mean_ws.theta)/2)]
theta_scew <- merge(theta_scew, bs[bs$variable == "weighted_mean_D",], by = c("facet", "subfacet", "maf", "mac"))
theta_scew <- theta_scew[,-8]
colnames(theta_scew)[8] <- "Tajima's D"

theta_scew_labs <- theta_scew
theta_scew_labs$weighted_mean_ts.theta <- theta_scew_labs$weighted_mean_ts.theta - 5
theta_scew_labs$weighted_mean_ws.theta <- theta_scew_labs$weighted_mean_ws.theta - 5

pad.x <- 2
pad.y <- -5
titles <- 16
text <- 14
ggplot(theta_scew, aes(x = weighted_mean_ts.theta, y = weighted_mean_ws.theta, color = maf)) +
  geom_point(size = 10) + geom_abline(slope = 1, intercept = 0) +
  geom_segment(aes(x = weighted_mean_ts.theta,
                   y = weighted_mean_ws.theta,
                   xend = weighted_mean_ts.theta + pad.x, 
                   yend = weighted_mean_ws.theta + pad.y)) +
  geom_label(aes(x = weighted_mean_ts.theta + pad.x, 
                 y = weighted_mean_ws.theta + pad.y, 
                 label = round(`Tajima's D`, 2)), size = 5.5) +
  scale_color_viridis_c() +
  theme_bw() +
  xlab(expression("Tajima\'s" ~ theta)) +
  ylab(expression("Watterson\'s" ~ theta)) +
  scale_y_continuous(limits = c(15, 90)) +
  scale_x_continuous(limits = c(50, 90)) +
  theme(axis.text = element_text(size = text),
        axis.title = element_text(size = titles),
        legend.text = element_text(size = text),
        legend.title = element_text(size = titles))


tsm <- data.table::melt(theta_scew, id.vars = c("facet", "subfacet", "maf", "mac", "scew", "Tajima's D"))
ggplot(tsm, aes(x = maf, y = value, shape = variable, color = `Tajima's D`)) +
  geom_line(size = 2) +
  geom_point(size = 8) +
  scale_color_viridis_c() +
  theme_bw() +
  ylab(expression(theta)) +
  theme(axis.text = element_text(size = text),
        axis.title = element_text(size = titles),
        legend.text = element_text(size = text),
        legend.title = element_text(size = titles)) +
  scale_shape_manual(labels = list(bquote(`Tajima's`~theta), bquote(`Watterson's`~theta)),
                     values = c(16, 17)) +
  xlab("MAF filter")


pc <- purrr::map(maf_res, "sfs")
pc <- dplyr::bind_rows(pc)
ggplot(pc, aes(x = mac, y = log10(count), color = maf, group = maf)) +
  geom_line() + theme_bw() +
  ylim(c(2, 5))


# bs <- purrr::map(maf_res, "stats")
# bs <- purrr::map(bs, "single")
# names(bs) <- mafs
# bs <- dplyr::bind_rows(bs, .id = "maf")
# bs$maf <- as.numeric(bs$maf)
# bs <- reshape2::melt(bs[,-c(2, 4, 5, 6)], id = c("maf", "subfacet"))


# bs <- purrr::map(maf_res, "stats")
# bs <- purrr::map(bs, "weighted.means")
# names(bs) <- mafs
# bs <- dplyr::bind_rows(bs, .id = "maf")
# bs$maf <- as.numeric(bs$maf)
# bs <- reshape2::melt(bs[,-c(2, 4, 5)], id = c("maf", "subfacet"))
# bs <- na.omit(bs)
# bs2 <- purrr::map(maf_res, "tsd")
# bs2 <- purrr::map(bs2, "weighted.means")
# names(bs2) <- mafs
# bs2 <- dplyr::bind_rows(bs2, .id = "maf")
# bs2 <- bs2[which(bs2$snp.subfacet == ".OVERALL_MEAN"),]
# bs2$maf <- as.numeric(bs2$maf)
# bs2 <- reshape2::melt(bs2[,-c(2, 4, 5)], id.vars = c("maf", "subfacet"))
# bs <- rbind(bs, bs2)
# bs$variable <- gsub("weighted_mean_", "", bs$variable)
# bs <- bs[-which(bs$variable == "pi"),]
# bs3 <- purrr::map(maf_res, "ne")
# bs3 <- purrr::map(bs3, "pop")
# names(bs3) <- mafs
# bs3 <- as.data.frame(dplyr::bind_rows(bs3, .id = "maf"))
# bs3$maf <- as.numeric(bs3$maf)
# bs3 <- bs3[,c(1,3:4)]
# bs3$variable <- "Ne"
# bs3 <- bs3[,c(1,2,4,3)]
# colnames(bs3) <- c("maf", "subfacet", "variable", "value")
# bs <- rbind(bs, bs3)
# bs <- bs[-which(bs$variable == "num_seg"),]
# bs$variable <- factor(bs$variable, levels = c("ho", "he", "fis", "prop_poly", "Ne", "fst", "D", "ts.theta", "ws.theta"))
# 
# library(ggplot2)
# ggplot(bs, aes(x = as.factor(maf), y = value, color = as.factor(maf), shape = as.factor(maf))) +
#   geom_point(size = 4) + theme_bw() +
#   facet_wrap(~variable, scales = "free_y") +
#   scale_color_viridis_d()
