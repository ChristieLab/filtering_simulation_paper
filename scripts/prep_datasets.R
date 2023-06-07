library(snpR); library(data.table)

#===========populus dataset============
populus <- readLines("data/populus/populus.txt")
populus <- strsplit(populus, "\t")
populus <- lapply(populus, function(x) gsub("\\|.+", "", x))
populus <- as.data.frame(populus)
populus <- as.data.table(t(populus))
colnames(populus) <- c("snpID", unlist(populus[1,-1]))
populus <- populus[-1,]
populus[populus == ""] <- "NN"
snp_meta <- populus[,1]
snp_meta[, c("p1", "p2", "position") := tstrsplit(snpID, "_", fixed = TRUE)]
snp_meta[,scaffold := paste0(p1, "_", p2)]

sample_meta <- openxlsx::read.xlsx("data/populus/evo12497-sup-0002-tables1.xlsx")
colnames(sample_meta) <- sample_meta[1,]
sample_meta <- sample_meta[-1,]
sample_meta <- sample_meta[match(colnames(populus)[-1], sample_meta$Accession),]
sample_meta <- sample_meta[,c(1, 3)]
colnames(sample_meta) <- c("sampleID", "pop")
sample_meta$pop <- gsub(" ", "_", sample_meta$pop)

pg <- import.snpR.data(populus[,-1], snp.meta = snp_meta[,c("scaffold", "position")], sample.meta = sample_meta)
saveRDS(pg, "data/populus/populus.RDS") # note, use Skeena, Vancouver_Island, and Quesnel drainages probably (large sample sizes, should have population structure)
# p <- plot_clusters(pg, "pop")
# p2 <- plot_clusters(pg[pop = c("Vancouver_Island", "Skeena", "Quesnel")], "pop")
# plotly::ggplotly(p$plot$pca)
# plotly::ggplotly(p2$plot$pca)


#===========mallards=======
mall <- read_structure("data/mallards/Supplementary Structure input file.str", header_cols = 2, 
                       rows_per_individual = 1, marker_names = TRUE, mDat = "0")
sample_meta <- sample.meta(mall)
colnames(sample_meta)[1:2] <- c("sampleID", "pop")
sample_meta$pop <- as.character(sample_meta$pop)
sample.meta(mall) <- sample_meta
# p2 <- plot_clusters(mall[pop = c("32", "12", "13")], "pop")$plot
# p <- plot_clusters(mall, "pop")$plot
saveRDS(mall, "data/mallards/mallards.RDS") #note, pops to use: 32, 12, 13 (three clusters, large sizes)


#============mosquitoes=======
mosquitos <- fread("data/mosquitos/UnfilteredSNPs_FL_SCA.ped")
map <- fread("data/mosquitos/UnfilteredSNPs_FL_SCA.map")
sample_meta <- mosquitos[,1:2]
colnames(sample_meta) <- c("pop", "sample_number")
mosquitos <- mosquitos[,-c(1:6)]
mosquitos <- as.matrix(mosquitos)
mosquitos[mosquitos == 0] <- "N"
dm <- dim(mosquitos)
pm <- paste0(mosquitos[,seq(1,ncol(mosquitos), by = 2)],
                    mosquitos[,seq(2,ncol(mosquitos), by = 2)])
pm <- matrix(pm, nrow = dm[1], ncol = dm[2]/2)
pm <- t(pm)
rm(mosquitos)
map <- map[,c(1,4)]
colnames(map) <- c("scaffold", "position")

pm <- import.snpR.data(pm, snp.meta = map, sample.meta = sample_meta)
p <- plot_clusters(pm, "pop")
plotly::ggplotly(p$plots$pca)
p2 <- plot_clusters(pm[pop = c("FM", "GG", "MV")], "pop")
plotly::ggplotly(p2$plots$pca)
saveRDS(pm, "data/mosquitos/UnfilteredSNPs_FL_SCA.RDS") # note: use MV, GG, and FM


#================jewelweed==============
jewelweed <- read_vcf("data/jewelweed/Dryad-populations.snps-stacks2.1.m3max3M4n4-filtered.vcf")
key <- fread("data/jewelweed/sample_name_key.csv")
sample.meta(jewelweed) <- as.data.frame(key[match(sample.meta(jewelweed)$sampID, key$Sample_name),])

p <- plot_clusters(jewelweed, "Population")
plotly::ggplotly(p$plots$pca)
p2 <- plot_clusters(jewelweed[Population = c("4033", "4079", "4057")], "Population")
plotly::ggplotly(p2$plots$pca)

saveRDS(jewelweed, "data/jewelweed/jewelweed.RDS")


