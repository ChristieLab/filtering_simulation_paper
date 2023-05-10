library(data.table); library(ggplot2)

fl <- list.files("results/", "\\.RDS", full.names = TRUE)
r <- vector("list", length(fl))

for(i in 1:length(fl)){
  print(i)
  
  y <- readRDS(fl[i])
  y$ld <- NULL # skip for now, too big
  if(!is.null(y$PCA)){
    y$pca <- y$pca$data[,c("samp", "pop", "PC1", "PC2")]
    y$pca$mac <- y$stats$mac[1]
    y$pca$maf <- y$stats$maf[1]
  }
  if(ncol(y$sfs) == 5){
    colnames(y$sfs) <- c("pop", "n_min_al", "count", "maf", "mac")
  }
  else{
    colnames(y$sfs) <- c("n_min_al", "count", "maf", "mac")
    
  }
  r[[i]] <- y
  names(r)[i] <- fl[i]
}

lseq <- function(from, to, length.out ) exp(seq(log(from), log(to), length.out = length.out))
mafs <- seq(0, .1, by = .01)
macs <- c(1, 3, 5)
hwes <- lseq(1e-6, .05, 10)

stats <- purrr::map(r, "stats")
pca <- purrr::map(r, "pca")
sfs <- purrr::map(r, "sfs")

iter <- gsub(".+res_r", "", fl)
iter <- gsub("\\.RDS", "", iter)
iter <- as.numeric(iter)

names(stats) <- iter
names(pca) <- iter
names(sfs) <- iter
stats <- data.table::rbindlist(stats, idcol = "iter")
pca <- data.table::rbindlist(pca, idcol = "iter")
sfs <- data.table::rbindlist(sfs, idcol = "iter", fill = TRUE)

stats$subfacet[stats$subfacet == ".base"] <- "NUB"
sfs$pop[is.na(sfs$pop)] <- "NUB"

hwe_tab <- data.frame(iter = 1:sum(length(mafs), length(macs), length(hwes)))
hwe_tab$hwe[hwe_tab$iter <= length(macs) + length(mafs)] <- 1e-6
hwe_tab$hwe[(length(macs) + length(mafs) + 1):sum(length(mafs), length(macs), length(hwes))] <- hwes

stats$iter <- as.numeric(stats$iter)
stats$hwe <- hwe_tab$hwe[match(stats$iter, hwe_tab$iter)]
sfs$hwe <- as.numeric(sfs$iter)
sfs$hwe <- hwe_tab$hwe[match(sfs$iter, hwe_tab$iter)]


stats$variable <- gsub("weighted_mean_", "", stats$variable)
plot_stats <- c("D", "ho", "num_seg")

ggplot(stats[stats$mac == 0 & !grepl("~", stats$subfacet) & stats$variable %in% plot_stats &
               stats$hwe == 1e-6,], 
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

ggplot(stats[!grepl("~", stats$subfacet) & stats$variable %in% plot_stats &
               stats$iter >= 15,], 
       aes(x = log10(hwe), y = value, color = subfacet)) +
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
  xlab("HWP p-value alpha")

ggplot(stats[stats$mac == 0 & stats$variable == "fst"&
               stats$hwe == 1e-6,], 
       aes(x = maf, y = value, color = subfacet)) +
  geom_line(size = 1) +
  theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2", name = "Population") +
  xlab("MAF filter") +
  ylab("Fst")

ggplot(stats[stats$variable == "fst"&
               stats$iter >= 15,], 
       aes(x = log10(hwe), y = value, color = subfacet)) +
  geom_line(size = 1) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2", name = "Population") +
  xlab("HWP p-value alpha") +
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

ggplot(sfs[sfs$mac == 0 & sfs$hwe == 1e-6,], 
       aes(x = n_min_al, y = log10(count), color = maf, group = maf)) +
  geom_line(size = 1) +
  theme_bw() +
  facet_wrap(~pop, scale="free") +
  khroma::scale_color_batlow(name = "MAF filter") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  xlab("Minor Allele Count") +
  ylab("log10(Number of SNPs)")


ggplot(sfs[sfs$iter >= 15,], 
       aes(x = n_min_al, y = log10(count), color = log10(hwe), group = iter)) +
  geom_line(size = 1) +
  theme_bw() +
  facet_wrap(~pop, scale="free") +
  khroma::scale_color_batlow(name = "HWE filter") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  xlab("Minor Allele Count") +
  ylab("log10(Number of SNPs)")

