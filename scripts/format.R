library(snpR)

# read in
args <- commandArgs(TRUE)
file <- as.character(args[1])
sample_meta <- as.character(args[2])
outfile <- as.character(args[3])


d <- read_vcf(file)
sample_meta <- read.table(sample_meta, header = F, sep = "\t")

# prepare a new snp metadata with "chr" instead of "CHROM" and no "."
snp_meta <- snp.meta(d)
snp_meta <- snp_meta[,c("CHROM", "position")]
snp_meta$chr <- gsub("\\.", "_", snp_meta$CHROM)
snp_meta <- snp_meta[,c("chr", "position")]

# read in sample metadata and set populations
sample.meta(d) <- data.frame(sampID = colnames(d),
                             pop = sample_meta[match(colnames(d), unlist(sample_meta[,2])),1])

# update with the new snp meta, done here for speed
snp.meta(d) <- snp_meta

# pre-tabulate pop
d <- calc_maf(d, "pop")

# save
saveRDS(d, outfile)
