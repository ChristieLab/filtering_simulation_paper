library(snpR); library(foreach); library(ggplot2); library(gridExtra)

#==========================with the satrapa data===========

# d <- read_structure("data/satrapa_structure/str_in/satrapa_unlinked.str", header_cols = 1, rows_per_individual = 2, mDat = "-9", marker_names = FALSE)
# sample.meta(d)$pop <- gsub("_", "", stringi::stri_extract(sample.meta(d)$V1, regex = "_[A-Z]{2}_"))
# # d <- d[pop = c("Chuckwall", "Dean", "Fraser", "Burke", "Lillooet")]
# d <- filter_snps(d, non_poly = TRUE)
# 
# setwd("data/processed/sat_structure")
# 
# macs <- 0:11
# res <- vector("list", 12)
# pal <- khroma::color("batlow")
# 
# for(i in macs){
#   dir.create(paste0("structure_mac_", i))
#   setwd(paste0("structure_mac_", i))
#   print(i)
#   if(i != 0){
#     td <- filter_snps(d, mac = i)
#   }
#   else{
#     td <- d
#   }
#   
#   
#   # if(nrow(td) > 5000){
#   #   td <- td[sort(sample(nrow(td), size = 5000, replace = FALSE)),]
#   # }
#   
#   res[[i+1]] <- plot_structure(td, facet = "pop", reps = 1, k = 3,
#                  burnin = 10000,
#                  iterations = 100000,
#                  structure_path = "C://usr/bin/structure.exe",
#                  method = "structure",
#                  alt.palette = pal(3),
#                  cleanup = FALSE)
#   setwd("..")
#   
#   # res[[i]] <- plot_clusters(td, "pop")
# }
# 
# for(i in macs){
#   setwd(paste0("structure_mac_", i))
#   res[[i + 1]] <- plot_structure("structure_outfile", k = 3,
#                                   alt.palette = pal(3))
#     
#     #              facet = sample.meta(d)$pop,
#     # #              facet.order = c("AZ", "BC", "CA", "CO", "ID", "WA", "WV", "NB", "MX"))$plot +
#     # theme(axis.title = element_blank(),
#     #       axis.text = element_blank(),
#     #       axis.ticks = element_blank(),
#     #       strip.text = element_blank(),
#     #       axis.text.x = element_blank()) +
#     # guides(color = "none", fill = "none")
#   setwd("..")
# }
# 
# 
# fix_clust <- function(x){
#   
#   if(length(x) == 1){
#     return(x)
#   }
#   
#   #loop through each q object
#   for (i in 2:length(x)){
#     #see which columns in the previous run are the most similar to each column
#     
#     #initialize mapping df
#     mdf <- data.frame(tcol = 1:ncol(x[[i]]), pcol = numeric(ncol(x[[i]])),
#                       ed = numeric(ncol(x[[i]])))
#     
#     #loop through each column and find where to map it.
#     for (j in 1:ncol(x[[i]])){
#       
#       #intialize euc distance vector
#       elist <- numeric(ncol(x[[i - 1]]))
#       
#       #compare to each other col.
#       for(tk in 1:ncol(x[[i-1]])){
#         #save euclidean dist
#         elist[tk] <- sum((x[[i]][,j] - x[[i-1]][,tk])^2, na.rm = T)
#       }
#       
#       #save results
#       mdf[j,2] <- which.min(elist)
#       mdf[j,3] <- min(elist)
#     }
#     
#     #reassign clusters in this qdf
#     ##which is the new cluster? Probably that with the most distance to any original clusters.
#     dups <- duplicated(mdf[,2]) | duplicated(mdf[,2], fromLast = T)
#     nc <- which.max(mdf[dups,3])
#     mdf[dups,2][nc] <- nrow(mdf)
#     mdf <- mdf[order(mdf[,2]),]
#     
#     ##reasign clusters
#     tdf <- x[[i]]
#     tdf <- tdf[,mdf[,1]]
#     
#     ##replace object in x with the re-arranged qfile.
#     colnames(tdf) <- colnames(x[[i]])
#     x[[i]] <- tdf
#   }
#   
#   return(x)
# }
# 
# qls <- purrr::map(res, "data")
# qls <- purrr::map(qls, "K_3")
# qls <- purrr::map(qls, "r_1")
# 
# fqls <- fix_clust(qls)
# for(i in macs){
#   fqls[[i + 1]]$mac <- i
#   fqls[[i + 1]]$pop <- sample.meta(d)$pop
#   fqls[[i + 1]]$ID <- 1:nrow(fqls[[i + 1]])
# }
# 
# fqls <- rbindlist(fqls)
# fqls$pop <- factor(fqls$pop, c("AZ", "BC", "CA", "CO", "ID", "WA", "WV", "NB", "MX"))
# fqls$mac <- factor(fqls$mac, 11:0)
# mfqls <- melt(fqls, id.vars = c("mac", "pop", "ID"))
# 
# 
# p <- ggplot2::ggplot(mfqls, ggplot2::aes(ID, value, color = variable, fill = variable)) +
#   ggplot2::facet_wrap(~mac, ncol = 1, strip.position = "right") +
#   ggplot2::theme_bw() +
#   ggplot2::geom_bar(stat = "identity", width = .9) +
#   ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5,0.75)) +
#   ggplot2::ylab("Cluster Membership Proportion") +
#   ggplot2::theme(axis.text.x = element_blank(),
#                  strip.text = ggplot2::element_text(size = 14),
#                  axis.title.x = element_blank(),
#                  axis.title.y = element_text(size = 18),
#                  strip.background = ggplot2::element_blank(),
#                  panel.grid = ggplot2::element_blank(),
#                  panel.spacing = ggplot2::unit(0.1, "lines"),
#                  panel.background = element_rect(fill = "black")) +
#   scale_x_discrete(expand = c(0,0))  +
#   khroma::scale_color_batlow(discrete = TRUE) +
#   khroma::scale_fill_batlow(discrete = TRUE) +
#   guides(fill = "none", color = "none")
# p
# 
# saveRDS(list(d = qls, p = p), "plot.RDS")


#==========with the monarch data=================
setwd("../monarch_structure")

d <- readRDS("../../monarch_nomaf.RDS") # subset this down to the same size as the satrapa data for consistancy
d <- d[pop = c("GUA", "ROT", "HAW", "NAM")]
d <- filter_snps(d, non_poly = TRUE)
d <- d[sample(nrow(d), 3819, FALSE)]

for(i in 2:11){
  dir.create(paste0("structure_mac_", i))
  setwd(paste0("structure_mac_", i))
  print(i)
  if(i != 0){
    td <- filter_snps(d, mac = i)
  }
  else{
    td <- d
  }
  
  res[[i+1]] <- plot_structure(td, facet = "pop", reps = 1, k = 3,
                               burnin = 10000,
                               iterations = 100000,
                               structure_path = "C://usr/bin/structure.exe",
                               method = "structure",
                               alt.palette = pal(3),
                               cleanup = FALSE)
  setwd("..")
  
  # res[[i]] <- plot_clusters(td, "pop")
}

for(i in macs){
  setwd(paste0("structure_mac_", i))
  res[[i + 1]] <- plot_structure("structure_outfile", k = 3,
                                 alt.palette = pal(3),
                                 facet = sample.meta(d)$pop,
                                 facet.order = c("NAM", "HAW", "GUA", "ROT"))
  
  #              facet = sample.meta(d)$pop,
  # #              facet.order = c("AZ", "BC", "CA", "CO", "ID", "WA", "WV", "NB", "MX"))$plot +
  # theme(axis.title = element_blank(),
  #       axis.text = element_blank(),
  #       axis.ticks = element_blank(),
  #       strip.text = element_blank(),
  #       axis.text.x = element_blank()) +
  # guides(color = "none", fill = "none")
  setwd("..")
}



qls <- purrr::map(res, "data")
qls <- purrr::map(qls, "K_3")
qls <- purrr::map(qls, "r_1")

fqls <- fix_clust(qls)
for(i in macs){
  fqls[[i + 1]]$mac <- i
  fqls[[i + 1]]$pop <- sample.meta(d)$pop
  fqls[[i + 1]]$ID <- 1:nrow(fqls[[i + 1]])
}

fqls <- rbindlist(fqls)
fqls$pop <- factor(fqls$pop, c("NAM", "HAW", "GUA", "ROT"))
fqls$mac <- factor(fqls$mac, 11:0)
mfqls <- melt(fqls, id.vars = c("mac", "pop", "ID"))
mfqls$variable <- gsub("Cluster", "", mfqls$variable)

breaks <- table(mfqls$pop)/(3*12)
breaks <- cumsum(breaks)
breaks <- c(0, breaks)
breaks <- rowMeans(cbind(breaks[-length(breaks)], breaks[-1]))

p <- ggplot2::ggplot(mfqls, ggplot2::aes(ID, value, color = variable, fill = variable)) +
  ggplot2::facet_wrap(~mac, ncol = 1, strip.position = "right") +
  ggplot2::theme_bw() +
  ggplot2::geom_bar(stat = "identity", width = .75) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5,0.75)) +
  ggplot2::ylab("Cluster Membership Proportion") +
  ggplot2::theme(strip.text.y = ggplot2::element_text(size = 10, angle = 0),
                 axis.title = element_text(size = 18),
                 axis.text.y = element_text(size = 14),
                 axis.text.x = element_text(size = 14, angle = 90),
                 strip.background = ggplot2::element_blank(),
                 axis.ticks.x = element_blank(),
                 panel.grid = ggplot2::element_blank(),
                 panel.spacing = ggplot2::unit(0.1, "lines"),
                 panel.background = element_rect(fill = "black")) +
  scale_x_continuous(expand = c(0,0), labels = c("NAM", "HAW", "GUA", "ROT"),
                   breaks = breaks - 1) +
  khroma::scale_color_vibrant() +
  khroma::scale_fill_vibrant() +
  guides(color = "none", fill = "none") +
  xlab("Population")
  
saveRDS(list(d = mfqls, p = p, br = breaks), "plot.RDS")


