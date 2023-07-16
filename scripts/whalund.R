library(ggplot2); library(snpR)
alpha <- .35
mon <- readRDS("data/monarch_nomaf.RDS")
mon <- mon[pop = c("GUA", "ROT")]
mon <- filter_snps(mon, non_poly = TRUE)

# mon <- calc_ho(mon)
# mon <- calc_ho(mon, "pop")
# mon <- calc_he(mon)
# mon <- calc_he(mon, "pop")
mon <- calc_fis(mon)
mon <- calc_fis(mon, "pop")
mon <- calc_pairwise_fst(mon, "pop")
mon <- calc_hwe(mon)
mon <- calc_hwe(mon, "pop")

rp <- get.snpR.stats(mon, "pop", c("fst"))
rs <- get.snpR.stats(mon, stats = c("ho", "he", "fis", "hwe"))

pd <- merge(rp$pairwise, rs$single, by = c("group", "position"))
pd$fst[pd$fst < 0] <- 0
# ggplot(pd, aes(x = fst, y = fis)) + geom_point(alpha = 0.5) +
#   geom_smooth() +
#   theme_bw() +
#   ylab(bquote(F[IS])) +
#   xlab(bquote(F[ST])) +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 18))
# 
# ggplot(pd, aes(x = fst, y = -log10(pHWE))) + geom_point(alpha = alpha) +
#   geom_smooth() +
#   geom_hline(yintercept = -log10(1e-6), mapping = aes(xmin = 0, xmax = 1), color = "red") +
#   theme_bw() +
#   ylab(bquote(-log[10](p[HWE]))) +
#   xlab(bquote(F[ST])) +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 18))

rsp <- get.snpR.stats(mon, "pop", c("fis", "hwe"))$single
rsp <- merge(rsp[which(rsp$subfacet == "GUA"),], rsp[which(rsp$subfacet == "ROT"),], by = c("group", "position"), 
             suffixes = c("_GUA", "_ROT"))
rsp <- as.data.table(rsp)
rsp[,mfis := rowMeans(cbind(fis_GUA, fis_ROT), na.rm = TRUE)]
rsp[,mHWE := rowMeans(cbind(pHWE_ROT, pHWE_GUA), na.rm = TRUE)]

pd2 <- merge(pd, rsp, by = c("group", "position"))
pd2 <- as.data.table(pd2)
pd2[,bad_reject := ((pHWE_ROT >= 1e-6 & pHWE_ROT >= 1e-6) & pHWE <= 1e-6)]

# ggplot(pd2, aes(x = mfis, y = fis, color = fst)) +
#   geom_point(alpha = alpha) +
#   geom_abline(slope = 1, intercept = 0, color = "red") +
#   khroma::scale_color_batlow() +
#   theme_bw()  +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 18)) +
#   xlab(bquote("Mean" ~ F[IS])) +
#   ylab(bquote("Combined" ~ F[IS]))
# 
# ggplot(pd2, aes(x = fst, y = fis - mfis)) +
#   geom_point(alpha = alpha) +
#   khroma::scale_color_batlow() +
#   theme_bw() +
#   geom_smooth() +
#   xlab(bquote(F[ST])) +
#   ylab(bquote("Combined" ~ - "Mean" ~ F[IS])) +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 18))
# 
# ggplot(pd2, aes(x = fst, y = pHWE - mHWE)) +
#   geom_point(alpha = alpha, aes(color = bad_reject)) +
#   scale_color_manual(values = c("black", "red")) +
#   theme_bw() +
#   geom_smooth() +
#   xlab(bquote(F[ST])) +
#   ylab(bquote("Combined" ~ - "Mean" ~ p[HWE])) +
#   guides(color = guide_legend(title = "Falsely Rejected")) +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 18))


pd3_1 <- melt(pd2[,c("group", "position", "fis", "fis_ROT", "fis_GUA")], id.vars = c("group", "position"))
pd3_1$level <- ifelse(grepl("ROT", pd3_1$variable), "ROT", ifelse(grepl("GUA", pd3_1$variable), "GUA", "Combined"))
colnames(pd3_1) <- c("group", "position", "variable", "fis", "level")
pd3_1$variable <- NULL

pd3_2 <- melt(pd2[,c("group", "position", "pHWE", "pHWE_ROT", "pHWE_GUA")], id.vars = c("group", "position"))
pd3_2$level <- ifelse(grepl("ROT", pd3_2$variable), "ROT", ifelse(grepl("GUA", pd3_2$variable), "GUA", "Combined"))
colnames(pd3_2) <- c("group", "position", "variable", "pHWE", "level")
pd3_2$variable <- NULL

pd3_3 <- pd2[,c("group", "position", "fst")]

pd3 <- merge(pd3_1, pd3_2, by = c("group", "position", "level"))
pd3 <- merge(pd3, pd3_3, by = c("group", "position"))

pd3[, fst_rank := rank(fst, ties.method = "last"), by = level]
pd3[, fis_rank := rank(fis, ties.method = "last"), by = level]

pd3$level <- as.factor(pd3$level)

ldm <- mgcv::gam(formula = fis ~ s(fst, bs = "cs") + 
                   # s(level, bs = "re", k = 3) +
                   s(level, fst, bs = "re"),  data = pd3, silent = TRUE)
nd <- tidygam::predict_gam(ldm, 1000)

pal <- khroma::color("batlow")

ggplot(nd, aes(x = fst, y = fis, color = level)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = level), alpha = .5, color = NA) +
  scale_color_manual(values = pal(4)[1:3]) +
  scale_fill_manual(values = pal(4)[1:3]) +
  theme_bw() +
  ylab(bquote(F[IS])) +
  xlab(bquote("Pairwise " ~ F[ST])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)) +
  guides(color = "none", fill = guide_legend(title = "Population"))


p <- ggplot(pd3[pd3$pHWE > 0.05,], aes(x = fst, y = fis)) +
  geom_point(alpha = alpha*2, color = "lightgray") +
  geom_point(data = pd3[pd3$pHWE <= 0.05,], aes(color = -log10(pHWE)), alpha = alpha) +
  geom_smooth(data = pd3, color = "red") +
  khroma::scale_color_batlow(range = c(.2, 1)) +
  theme_bw() +
  facet_wrap(~level) +
  ylab(bquote(F[IS])) +
  xlab(bquote(F[ST])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

ggsave("data/processed/whalund_plot.pdf", p, "pdf", width = 15, height = 11)



ggplot(mapping = aes(x = fst, y = fis)) +
  geom_point(data = pd3[which(pd3$pHWE > 0.05),], color = "gray", alpha = alpha/3) +
  geom_point(data = pd3[which(pd3$pHWE <= 0.05),], aes(color = -log10(pHWE)), alpha = alpha) +
  geom_smooth(data = pd3, color = "red", method = "lm") +
  khroma::scale_color_batlow(range = c(.2, 1)) +
  theme_bw() +
  facet_wrap(~level) +
  ylab(bquote(F[IS])) +
  xlab(bquote(F[ST])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

