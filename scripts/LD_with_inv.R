library(snpR); library(ggplot2)
stick <- data.table::fread("data/stickleback.txt", colClasses = "character")
stick$position <- as.numeric(stick$position)
stick <- stick[-grep("scaffold", stick$group),]
stick <- import.snpR.data(stick[,-c(1:3)], snp.meta = stick[,1:3], 
                          sample.meta = data.frame(pop = substr(colnames(stick)[-c(1:3)], 1, 3)),
                          mDat = "0000")

format_snps(stick, "vcf", outfile = "data/all_STB_chr.vcf", chr = "group")
format_snps(stick[group = "groupIX"], "vcf", outfile = "data/inv_STB_chr.vcf", chr = "group")
format_snps(stick[group = "groupXIX"], "vcf", outfile = "data/sex_STB_chr.vcf", chr = "group")


stick <- stick[group = "groupIX"]
stick <- filter_snps(stick, mgc = 1, hwe = 1e-6, hwe_facets = "pop", min_loci = .5, min_ind = .5)


stick <- calc_pairwise_ld(stick, "pop.group", par = 4, verbose = TRUE)

prox <- get.snpR.stats(stick, "pop.group", "LD")$prox

ggplot(prox[prox$sample.subfacet == "ASP",], aes(x = proximity, y = CLD)) + 
  geom_point(alpha = .5) +
  theme_bw()

