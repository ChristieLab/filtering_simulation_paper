library(data.table); library(ggplot2)

fl <- list.files("results/", "maf_res", full.names = TRUE)
r <- vector("list", length(fl))

for(i in 1:length(fl)){
  print(i)
  y <- readRDS(fl[i])
  y$pca <- y$pca$data[,c("samp", "pop", "PC1", "PC2")]
  y$pca$mac <- y$stats$mac[1]
  y$pca$maf <- y$stats$maf[1]
  colnames(y$sfs) <- c("n_min_al", "count", "maf", "mac")
  r[[i]] <- y
}


stats <- purrr::map(r, "stats")
pca <- purrr::map(r, "pca")
sfs <- purrr::map(r, "sfs")

stats <- data.table::rbindlist(stats)
pca <- data.table::rbindlist(pca)
sfs <- data.table::rbindlist(sfs)

stats$variable <- gsub("weighted_mean_", "", stats$variable)

ggplot(stats[stats$mac == 0 & stats$variable != "fst" & !grepl("~", stats$subfacet) & stats$variable != "nsnps",], 
       aes(x = maf, y = value, color = subfacet)) +
  geom_line(size = 1) +
  facet_wrap(~variable, scales = "free_y", strip.position = "left") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2", name = "Population") +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  xlab("MAF filter")

ggplot(stats[stats$mac == 0 & stats$variable == "fst",], 
       aes(x = maf, y = value, color = subfacet)) +
  geom_line(size = 1) +
  theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2", name = "Population") +
  xlab("MAF filter") +
  ylab("Fst")
  
ggplot(pca[pca$mac == 0,], 
       aes(x = PC1, y = PC2, color = pop)) +
  geom_point() +
  facet_wrap(~maf) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2", name = "Population") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

ggplot(sfs[sfs$mac == 0,], 
       aes(x = n_min_al, y = log10(count), color = maf, group = maf)) +
  geom_line(size = 1) +
  theme_bw() +
  khroma::scale_color_batlow(name = "MAF filter") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  xlab("Minor Allele Count") +
  ylab("log10(Number of SNPs)")
