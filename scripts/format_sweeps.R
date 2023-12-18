library(snpR)

args <- commandArgs(TRUE)
file <- as.character(args[1])
sim <- readRDS(file)

len <- 10000000

genos <- sim$seg_sites

i <- 1
meta <- data.table(chr = i, position = floor(as.numeric(genos[[i]]$position)*len))
genos[[i]] <- genos[[i]]$snps[seq(1, nrow(genos[[i]]$snps), 2),] +
  genos[[i]]$snps[seq(2, nrow(genos[[i]]$snps), 2),]
genos[[i]] <- t(genos[[i]])
genos<- as.data.table(genos[[i]])

keeps <- c(sample(51:100, 30, FALSE), sample(101:150, 30, FALSE))
keeps <- sort(keeps)
genos <- as.data.frame(genos)

srd <- import.snpR.data(genos[,keeps], meta, data.frame(pop = rep(c("B", "C"), each = 30)))

srd <- srd[pop = c("B", "C")]
srd <- calc_maf(srd, "pop")

ofile <- tools::file_path_sans_ext(file)
ofile <- file.path("snpR_RDS", paste0(ofile, ".RDS"))

saveRDS(srd, ofile)
