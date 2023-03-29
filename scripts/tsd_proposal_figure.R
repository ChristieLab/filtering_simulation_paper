library(snpR); library(ggplot2)
d <- readRDS("data/monarch_nomaf.RDS")
d <- d[pop = "NAM"]

d1 <- filter_snps(d, hwe = 1e-6)
d1 <- calc_tajimas_d(d1, facets = "group", sigma = 50, triple_sigma = FALSE, par = 4)
d1 <- get.snpR.stats(d1, "group", "tajimas_d")

d2 <- filter_snps(d, hwe = 1e-6, maf = 0.02)
d2 <- calc_tajimas_d(d2, facets = "group", sigma = 50, triple_sigma = FALSE, par = 4)
d2 <- get.snpR.stats(d2, "group", "tajimas_d")

d3 <- filter_snps(d, hwe = 1e-6, maf = 0.05)
d3 <- calc_tajimas_d(d3, facets = "group", sigma = 50, triple_sigma = FALSE, par = 4)
d3 <- get.snpR.stats(d3, "group", "tajimas_d")

com <- dplyr::bind_rows(list("0" = d1$single.window,
                             "0.02" = d2$single.window,
                             "0.05" = d3$single.window), .id = "maf")


com$maf <- as.numeric(com$maf)
p <- plot_manhattan(com, plot_var = "D", chr = "snp.subfacet")$plot

dat <- p$data
dat <- as.data.table(dat)
dat[,D_mean := mean(pvar), by = as.factor(maf)]
dat[,chr_maf := paste0(maf, "_", chr)]

ggplot(dat, aes(x = cum.bp, y = pvar, color = as.factor(maf))) +
  geom_line(aes(group = chr_maf), alpha = 0.4) + 
  # geom_point(alpha = 0.3) +
  geom_hline(aes(yintercept = D_mean, color = as.factor(maf)), size = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 26),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 26)) +
  scico::scale_color_scico_d(palette = "batlow", end = .7) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 5e7, 1e8, 1.5e8),
                     labels = c("0", "50", "100", "150")) +
  guides(color=guide_legend(title = "MAF filter")) +
  ylab("Tajima's D") +
  xlab("Cumulative Genomic Position (mb)")
