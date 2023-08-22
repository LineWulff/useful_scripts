

setwd("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R5/clustering/MNP/")

Idents(MNP_LP) <- MNP_LP@meta.data$res0.5_ident
levels(MNP_LP)
levels(MNP_LP) <- c("mono","mono int","mac","act_1","act_2","cDC2B","cDC2","cDC1","prol DC","mig cDC2_1","mono int_1","pDC","CD16 mono","debris","mig cDC2_2")

Idents(MNP_LP) <- MNP_LP@meta.data$integrated_snn_res.0.5
levels(MNP_LP)
levels(MNP_LP) <- c("2","0","3","5","7","1","4","6","8","9","10","11","12","13","14")

Idents(MNP_LP) <- MNP_LP@meta.data$integrated_snn_res.0.5
markers_res0.5 <- FindAllMarkers(MNP_LP,
                                 only.pos = T, 
                                 test.use = "wilcox",
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)
col_HM <- c("#C99800","#F8766D","#A3A500","#00BA38","#00C0AF","#E58700","#6BB100","#00BF7D","#00BCD8","#00B0F6","#619CFF","#B983FF","#E76BF3","#FD61D1","#FF67A4")

round="MNP"
project="MNP_LP"
top10 <- as.data.frame(markers_res0.5 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC))
png(paste(dato,project,round,"R2_res0.5_top10DEG_ordered.png",sep="_"),height = 2000, width = 2000, res = 150)
DoHeatmap(MNP_LP, features = top10$gene, group.colors = col_HM) + NoLegend()
dev.off()


