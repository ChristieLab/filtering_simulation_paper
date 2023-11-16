library(data.table)

sf <- c(.05, 1/(2*10000))
ts <- c(50, 100, 200)
r <- c(seq(.1, 1, by = .1), seq(2, 10, by = 1))
strn <- .2

d <- expand.grid(sf, ts, strn, r)
d <- unique(d)
options(scipen = 999)
d$out <- paste0(d[,1], "_", d[,2], "_", d[,3], "_", d[,4])


data.table::fwrite(d, "scripts/sweep_sim_parms.txt", quote = F, row.names = F, col.names = F, sep = "\t")
