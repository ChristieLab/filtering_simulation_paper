library(data.table); library(ggplot2); library(mgcv); library(tidygam); library(gridExtra)

fl <- list.files("results/", "\\.RDS", full.names = TRUE)
r <- vector("list", length(fl))

bind_parms <- function(y, parms, i){
  y$maf <- parms[1,4]
  y$mac <- ifelse(parms[1,5] == 0, FALSE, parms[1,5])
  y$hwe <- parms[1,6]
  y$miss <- parms[1,8]
  y$i <- i
  return(y)
}

for(i in 1:length(fl)){
  print(i)
  
  y <- readRDS(fl[i])
  
  uf <- unique(y$res$ld$sample.subfacet)
  nd <- vector("list", length(uf))
  for(j in 1:length(uf)){
    if(nrow(y$res$ld[y$res$ld$sample.subfacet == uf[j],]) > 2){
      try(ldm <- gam(formula = CLD ~ s(proximity, bs = "cs"), data = y$res$ld[y$res$ld$sample.subfacet == uf[j],]), silent = TRUE)
      try(nd[[j]] <- tidygam::predict_gam(ldm, 100), silent = TRUE)
      nd[[j]]$sample.subfacet <- uf[j]
    }
  }
  nd <- rbindlist(nd)
  
  r[[i]]$ld <- bind_parms(nd, y$parms, i)
  r[[i]]$stats <- bind_parms(y$res$stats, y$parms, i)
  y$res$sfs$mac.1 <- NULL
  colnames(y$res$sfs)[which(colnames(y$res$sfs) == "mac")] <- "min_allele_count"
  r[[i]]$sfs <- bind_parms(y$res$sfs, y$parms, i)
  if(!is.null(y$res$pca)){
    r[[i]]$pca <- bind_parms(as.data.frame(y$res$pca$data[,1:5]), y$parms, i)
  }
  
  names(r)[i] <- fl[i]
  gc();gc()
}

saveRDS(r, "data/processed/filt_results.RDS")

stats <- purrr::map(r, "stats")
pca <- purrr::map(r, "pca")
sfs <- purrr::map(r, "sfs")
ld <- purrr::map(r, "ld")

stats <- data.table::rbindlist(stats)
pca <- data.table::rbindlist(pca, fill = TRUE)
sfs <- data.table::rbindlist(sfs, fill = TRUE)
ld <- data.table::rbindlist(ld)

# stats$subfacet[stats$subfacet == ".base"] <- "NUB"
# sfs$pop[is.na(sfs$pop)] <- "NUB"

stats$variable <- gsub("weighted_mean_", "", stats$variable)

multigrepl <- function(pattern, x){
  return(Reduce(f = "|", lapply(c(pattern), function(y) grepl(y, x))))
}

stats$dataset <- ifelse(multigrepl(c("IMC", "GBT", "NUB"), stats$subfacet), "Yellow Perch", 
                        ifelse(multigrepl(c("WAL", "QUI", "ROC"), stats$subfacet),
                               "Cutthroat",
                               "Monarchs"))

num_seg <- stats[stats$variable == "prop_poly",]
num_seg <- merge(num_seg, unique(stats[stats$variable == "nsnps",c("value", "i")]),
                 by = "i")
num_seg[, value := value.x * value.y]
num_seg[,variable := "num_seg"]
num_seg <- num_seg[,c("facet", "subfacet", "variable", "value", "maf", "mac", "hwe", "miss", "i", "dataset")]

stats <- stats[-which(stats$variable == "num_seg"),]
stats <- rbind(stats, num_seg)
stats <- stats[-which(stats$variable == "prop_poly"),]
# stats <- stats[-which(stats$dataset == "Monarchs" &
#                         !stats$subfacet %in% c("GUA", "NAM", "HAW",
#                                                "GUA~NAM", "GUA~HAW", "HAW~NAM")),]

stats$mnseg <- 0
stats[,mnseg := max(as.numeric(value), na.rm = TRUE), by = .(dataset, subfacet, variable)]
stats$mnseg[stats$variable != "num_seg"] <- NA
stats$value[stats$variable == "num_seg"] <- stats$value[stats$variable == "num_seg"]/stats$mnseg[stats$variable == "num_seg"]
stats$mnseg <- NULL
var_tab <- data.frame(ref = unique(stats$variable),
                      var_p = c("pi", "H[O]", "F[IS]", "F[ST]", "Tajima*'\\''*'s'*' D'", "Tajima*'\\''*'s'*' 'theta",
                                "Watterson*'\\''*'s'*' 'theta", "N[e]", "N[snp]", "P[seg]"))

stats$variable <- var_tab$var_p[match(stats$variable, var_tab$ref)]


plot_stats <- c("Tajima*'\\''*'s'*' D'", "H[O]", "P[seg]", "F[ST]")


# maf
pa <- ggplot(stats[stats$mac == FALSE  & stats$variable %in% plot_stats &
                     stats$hwe == 1e-6 & stats$miss == .5,], 
             aes(x = maf, y = value, color = dataset, group = subfacet)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~variable, scales = "free_y", strip.position = "left", labeller = label_parsed, ncol = 1) +
  theme_bw() +
  khroma::scale_color_highcontrast() +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)) +
  xlab("MAF filter")

# miss
pb <- ggplot(stats[stats$mac == 1 & stats$variable %in% plot_stats &
               stats$hwe == 1e-6,], 
       aes(x = 1-miss, y = value, color = dataset, group = subfacet)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~variable, scales = "free_y", strip.position = "left", labeller = label_parsed, ncol = 1) +
  scale_x_reverse() +
  theme_bw() +
  khroma::scale_color_highcontrast() +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18)) +
  xlab("% Missing Data")





pl1 <- ggpubr::get_legend(pa)
# pl2 <- ggpubr::get_legend(pa2)


# grid.arrange(pa + guides(color = "none"), 
#              pb + guides(color = "none"), 
#              pa2 + guides(color = "none"), 
#              pb2 + guides(color = "none"),
#              pl1,
#              layout_matrix = matrix(c(1,2,3,4,5,5), nrow = 2),
#              widths = c(1,  .33333, .2))



# load in the structure data

str <- readRDS("data/processed/monarch_structure/plot.RDS")

# plot everything
g1 <- arrangeGrob(pa + guides(color = "none") + ggtitle("A") + 
                    theme(plot.title = element_text(size = 18, hjust = -.16)), 
                  pb + guides(color = "none") + ggtitle("B")+ 
                    theme(plot.title = element_text(size = 18, hjust = -.16)),
                  pl1, 
                  layout_matrix = matrix(c(1,2,3), nrow = 1), 
                  widths = c(1, 1, .3))

gb <- arrangeGrob(g1, str$p + guides(color = "none", fill = "none") +
                    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                          strip.text = element_text(angle = 0),
                          plot.title = element_text(size = 18, hjust = -.02)) +
                    ylab("Cluster")  +
                    ggtitle("C"),
                  heights = c(1, .6))

ggsave("data/processed/box2.pdf", gb, "pdf", height = 11, width = 15)

# ggplot(ld[ld$mac == FALSE &
#                   ld$hwe == 1e-6 & ld$miss == .5,], 
#              aes(x = proximity, y = CLD, color = maf, ymin = lower_ci, ymax = upper_ci, group = maf)) +
#   geom_line(size = 1) +
#   facet_wrap(~sample.subfacet, strip.position = "left", scales = "free") +
#   theme_bw() +
#   khroma::scale_color_batlow(name = "MAF filter") +
#   theme(axis.title.y = element_blank(),
#         strip.background = element_blank(), 
#         strip.placement = "outside",
#         strip.text = element_text(size = 18),
#         axis.title.x = element_text(size = 18),
#         axis.text = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 18)) +
#   xlab("Proximity (bp)") +
#   ylab("LD (r2)")
# 
# 
# ggplot(sfs[sfs$mac == FALSE & sfs$hwe == 1e-6 & sfs$miss == .5,], 
#        aes(x = min_allele_count, y = log10(count), color = maf, group = maf)) +
#   geom_line(size = 1) +
#   theme_bw() +
#   facet_wrap(~pop, scale="free") +
#   khroma::scale_color_batlow(name = "MAF filter") +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14)) +
#   xlab("Minor Allele Count") +
#   ylab("log10(Number of SNPs)")


# # hwe
# ggplot(stats[!grepl("~", stats$subfacet) & stats$variable %in% plot_stats & 
#                stats$mac == 1 & stats$maf == FALSE & stats$miss == .5,],
#        aes(x = log10(hwe), y = value, color = subfacet)) +
#   geom_line(size = 1) +
#   facet_wrap(~variable, scales = "free_y", strip.position = "left") +
#   theme_bw() +
#   scale_color_brewer(palette = "Set3", name = "Population") +
#   theme(axis.title.y = element_blank(),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         strip.text = element_text(size = 18),
#         axis.title.x = element_text(size = 18),
#         axis.text = element_text(size = 14)) +
#   xlab("HWP p-value alpha")

# # fst maf
# pa2 <- ggplot(stats[stats$mac == FALSE & grepl("~", stats$subfacet) & stats$variable == "F[ST]" &
#                stats$hwe == 1e-6 & stats$miss == .5,], 
#        aes(x = maf, y = value, color = dataset, group = subfacet)) +
#   geom_line(size = 1) +
#   geom_point(size = 3) +
#   facet_wrap(~variable, scales = "free_y", strip.position = "left", labeller = label_parsed) +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         strip.background = element_blank(), 
#         strip.placement = "outside",
#         strip.text = element_text(size = 18),
#         axis.title.x = element_text(size = 18),
#         axis.text = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 18)) +
#   khroma::scale_color_highcontrast() +
#   xlab("MAF filter") +
#   ylab("Fst")
# 
# 
# # fst miss
# pb2 <- ggplot(stats[stats$mac == 1 & grepl("~", stats$subfacet) & stats$variable == "F[ST]" &
#                stats$hwe == 1e-6,], 
#        aes(x = 1-miss, y = value, color = dataset, group = subfacet)) +
#   geom_line(size = 1) +
#   geom_point(size = 3) +
#   facet_wrap(~variable, scales = "free_y", strip.position = "left", labeller = label_parsed) +
#   scale_x_reverse() +
#   theme_bw() +
#   khroma::scale_color_highcontrast() +
#   theme(axis.title.y = element_blank(),
#         strip.background = element_blank(), 
#         strip.placement = "outside",
#         strip.text = element_text(size = 18),
#         axis.title.x = element_text(size = 18),
#         axis.text = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 18)) +
#   xlab("Maximum proportion missing data")

# # fst hwe
# ggplot(stats[grepl("~", stats$subfacet) &
#                stats$mac == 1 & stats$maf == FALSE & stats$miss == .5,],
#        aes(x = log10(hwe), y = value, color = subfacet)) +
#   geom_line(size = 1) +
#   facet_wrap(~variable, scales = "free_y", strip.position = "left") +
#   theme_bw() +
#   scale_color_brewer(palette = "Set3", name = "Population") +
#   theme(axis.title.y = element_blank(),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         strip.text = element_text(size = 18),
#         axis.title.x = element_text(size = 18),
#         axis.text = element_text(size = 14)) +
#   xlab("HWP p-value alpha")
# 
# 
# # pca maf
# ggplot(pca[pca$mac == FALSE & 
#              pca$hwe == 1e-6 & pca$miss == .5,], 
#        aes(x = PC1, y = PC2, color = pop)) +
#   geom_point() +
#   facet_wrap(~maf) +
#   theme_bw() +
#   scale_color_brewer(palette = "Set3", name = "Population") +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 18),
#         axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14))
# 
# # pca hwe
# ggplot(pca[pca$mac == 1 & pca$miss == .5,], 
#        aes(x = PC1, y = PC2, color = pop)) +
#   geom_point() +
#   facet_wrap(~hwe) +
#   theme_bw() +
#   scale_color_brewer(palette = "Set3", name = "Population") +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 18),
#         axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14))
# 
# 
# # pca hwe
# ggplot(pca[pca$mac == 1 & pca$hwe == 1e-6,], 
#        aes(x = PC1, y = PC2, color = pop)) +
#   geom_point() +
#   facet_wrap(~miss) +
#   theme_bw() +
#   scale_color_brewer(palette = "Set3", name = "Population") +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 18),
#         axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14))
# 
# 
# # sfs maf
# ggplot(sfs[sfs$mac == FALSE & sfs$hwe == 1e-6 & sfs$miss == .5,], 
#        aes(x = min_allele_count, y = log10(count), color = maf, group = maf)) +
#   geom_line(size = 1) +
#   theme_bw() +
#   facet_wrap(~pop, scale="free") +
#   khroma::scale_color_batlow(name = "MAF filter") +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14)) +
#   xlab("Minor Allele Count") +
#   ylab("log10(Number of SNPs)")
# 
# 
# # sfs hwe
# ggplot(sfs[sfs$mac == 1 & sfs$miss == .5 & sfs$maf == FALSE,], 
#        aes(x = min_allele_count, y = log10(count), color = log10(hwe), group = hwe)) +
#   geom_line(size = 1) +
#   theme_bw() +
#   facet_wrap(~pop, scale="free") +
#   khroma::scale_color_batlow(name = "MAF filter") +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14)) +
#   xlab("Minor Allele Count") +
#   ylab("log10(Number of SNPs)")
# 
# 
# # sfs miss
# ggplot(sfs[sfs$mac == 1 & sfs$hwe == 1e-6 & sfs$maf == FALSE,], 
#        aes(x = min_allele_count, y = log10(count), color = miss, group = miss)) +
#   geom_line(size = 1) +
#   theme_bw() +
#   facet_wrap(~pop, scale="free") +
#   khroma::scale_color_batlow(name = "MAF filter") +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14)) +
#   xlab("Minor Allele Count") +
#   ylab("log10(Number of SNPs)")
