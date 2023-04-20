library(Biostrings); library(data.table)
offset_overlapping_positions <- function(meta){
  meta[,2] <- round(meta[,2])
  meta$ord <- 1:nrow(meta)
  
  # function that checks for duplicates, returns dups, sorted meta, and pos_id
  dup_check <- function(meta){
    pos_id <- paste0(meta[,1], "_", meta[,2])
    meta <- meta[order(pos_id),]
    pos_id <- sort(pos_id)
    dups <- duplicated(pos_id) | rev(duplicated(rev(pos_id)))
    return(list(dups = dups, meta = meta, pos_id = pos_id))
  }
  
  # function that takes a dup check result and adjusts any duplicate sites
  adjust_pos <- function(dc){
    # unpack dc
    meta <- dc$meta
    dups <- dc$dups
    pos_id <- dc$pos_id
    
    # get a dup count table
    tab <- table(pos_id[dups])
    
    # make a new vector of positions for each duplicate site and overwrite
    ## prepare vectors
    cent <- ceiling(tab/2)
    current_pos <- meta[match(names(tab), pos_id), 2]
    lower <- current_pos - (cent - 1)
    upper <- ifelse(tab %% 2 == 0, current_pos + cent, current_pos + cent - 1)
    vec <- apply(matrix(c(lower, upper), ncol = 2), 1, function(x){x[1]:x[2]})
    
    ## overwrite and return
    meta[dups,2] <- unlist(vec)
    
    return(meta)
  }
  
  # run adjustment and re-dup check until there are no duplicated sites
  dc <- dup_check(meta)
  while(any(dc$dups)){
    meta <- adjust_pos(dc) # adjust
    dc <- dup_check(meta) # dup check
  }
  
  meta <- meta[order(meta$ord),]
  meta$ord <- NULL
  
  return(meta)
}

# read in simulated data
d <- readRDS("data/raw_mig_sim_dat.RDS")
d <- unlist(d, recursive = FALSE)

# read in and sort monarch genome
mon_ref <- readDNAStringSet("data/ncbi-genomes-2023-03-29/GCF_009731565.1_Dplex_v4_genomic.fna.gz")
mon_ref <- mon_ref[grep("NC_", names(mon_ref))]

ord <- gsub("^.+chromosome ", "", names(mon_ref))
ord <- gsub(",.+", "", ord)
ord <- gsub("Z", "1", ord)
ord <- as.numeric(ord)
ord <- order(ord)

mon_ref <- mon_ref[ord]

# get positions of mutations from sim data
conv <- purrr::map(d, "position")
for(i in 1:length(conv)){
  conv[[i]] <- data.frame(chr = names(mon_ref[i]), position = unlist(conv[[i]]))
  conv[[i]]$position <- floor(conv[[i]]$position*length(mon_ref[[i]]))
}
conv <- dplyr::bind_rows(conv)
conv <- offset_overlapping_positions(conv)

# grab the seg sites
d <- purrr::map(d, "snps")
names(d) <- names(mon_ref)
d <- lapply(d, function(x) as.data.frame(t(x)))
d <- data.table::rbindlist(d, idcol = "chr")

# make an alternative reference with all of the non-ref alleles
# alt_ref <- mon_ref
# for(i in 1:length(alt_ref)){
#   ref_seq <- strsplit(as.character(mon_ref[i][[1]][conv[conv$chr == names(mon_ref)[i],2]]), "")[[1]]
#   if(any(ref_seq == "N")){
#     mon_ref[i][[1]][conv[conv$chr == names(mon_ref)[i],2]][which(ref_seq == "N")] <-
#       Biostrings::DNAString(paste0(sample(c("A", "T", "C", "G"), sum(ref_seq == "N"), replace = TRUE), collapse = ""))
#   }
#   tref <- alt_ref[i][[1]][conv[conv$chr == names(mon_ref)[i],2]]
#   alt <- sample(c("A", "T", "C", "G"), length(tref), replace = TRUE)
#   matches <- which(unlist(strsplit(as.character(tref), "")) == alt)
#   while(length(matches) > 0){
#     alt[matches] <- sample(c("A", "T", "C", "G"), length(matches), replace = TRUE)
#     matches <- which(unlist(strsplit(as.character(tref), "")) == alt)
#   }
#   alt_ref[i][[1]][conv[conv$chr == names(mon_ref)[i],2]] <- Biostrings::DNAString(paste0(alt, collapse = ""))
# }

# figure out ref and alt alleles for each seg site
ref_alt <- data.table(chr = character(),
                      position = numeric(),
                      ref = character(),
                      alt = character())

for(i in 1:length(mon_ref)){
  ref_seq <- unlist(strsplit(as.character(mon_ref[i][[1]]), ""))[conv$position[conv$chr == names(mon_ref)[i]]]
  if(any(ref_seq == "N")){
    fill.pos <- which(ref_seq == "N")
    fill <- sample(c("A", "T", "C", "G"), sum(ref_seq == "N"), replace = TRUE)
    ref_seq[fill.pos] <-
      fill
    mon_ref[i][[1]][conv$position[conv$chr == names(mon_ref)[i]]][fill.pos] <- Biostrings::DNAString(paste0(fill, collapse = ""))
  }
  
  alt_seq <- sample(c("A", "T", "C", "G"), length(ref_seq), replace = TRUE)
  matches <- which(ref_seq == alt_seq)
  while(length(matches) > 0){
    alt_seq[matches] <- sample(c("A", "T", "C", "G"), length(matches), replace = TRUE)
    matches <- which(ref_seq == alt_seq)
  }
  ref_alt <- rbind(ref_alt, 
                   data.table(chr = names(mon_ref)[i],
                              position = conv$position[conv$chr == names(mon_ref)[i]],
                              ref = ref_seq,
                              alt = alt_seq))
}

# make a VCF with the ref and alt
p1 <- seq(1, 400, by = 2)
p2 <- seq(2, 400, by = 2)
vcf <- paste0(unlist(d[,-1][,..p1]), "|", unlist(d[,-1][,..p2]))
vcf <- matrix(vcf, nrow(d), 200)
vcf <- as.data.table(vcf)
colnames(vcf) <- paste0(rep(c("pop_A_", "pop_B_"), each = 100), rep(1:100, length.out = 200))
write("##fileformat=VCFv4.0\n##fileDate=20230330\n##source=scrm\n##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele Count\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "data/monarch_sim_gene_copies/variants.vcf")
ac <- rowSums(d[,-1])
chr.sim.names <- gsub("^.+chromosome ", "", ref_alt$chr)
chr.sim.names <- gsub(",.+", "", chr.sim.names)
chr.sim.names <- gsub("Z", "1", chr.sim.names)
chr.sim.names <- paste0("chr_", chr.sim.names)
vcf <- cbind(data.table(`#CHROM` = chr.sim.names,
                        POS = ref_alt$position,
                        ID = paste0("SNP_", 1:nrow(d)),
                        REF = ref_alt$ref,
                        ALT = ref_alt$alt,
                        QUAL = ".",
                        FILTER = "PASS",
                        INFO = paste0("NS=200;AC=",ac),
                        FORMAT = "GT"),
             vcf)

fwrite(vcf, "data/monarch_sim_gene_copies/variants.vcf", 
            quote = F,
            sep = "\t", 
            append = TRUE)
Biostrings::writeXStringSet(x = mon_ref,
                            filepath = "data/monarch_sim_gene_copies/N_filled_ref.fasta")




# make fastas for each gene copy
# for(i in 1:400){
#   cat(i, "\n")
#   this_gc <- mon_ref
#   tc <- i + 1
#   for(j in 1:length(mon_ref)){
#     flips <- conv$position[which(d$chr == names(mon_ref[j]) & d[[tc]] == 1)]
#     this_gc[[j]][flips] <- alt_ref[[j]][flips]
#   }
#   
#   of <- paste0("data/monarch_sim_gene_copies/gc_", i, ".fasta")
#   Biostrings::writeXStringSet(x = this_gc,
#                               filepath = of)
# }
