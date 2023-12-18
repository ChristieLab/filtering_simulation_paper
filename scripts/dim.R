library(snpR)
#=============parse command args================
args <- commandArgs(TRUE)
infile <- as.character(args[1])

#================filter============
cat("Filtering.\n")

# read in and change to temp dir
f <- readRDS(infile)
f <- filter_snps(f, non_poly = TRUE)

cat("file: ", infile, "\n")
cat("SNPs: ", nrow(f), "\nSamples: ", ncol(f), "\n")
