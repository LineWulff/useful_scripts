##### Example of plotting with
##### DoMultipleBarHeatmap 

visu_obj <- readRDS(file="/Volumes/Mucosal-immunology/WA group/Line/10x_other/Kristina/analysis/multiple_merge/MON_WT-DEL-79/R1/Rdata/R2_MON_traj_v1000_t100_v2.rds")

Used_markers <- c("Cx3cr1","Ly6c1", "Ly6c2","Ly6a", "Cxcl1", "Cxcr3","Mafb")

#run####
# Setting groups of metdatato plot from and their colours (won't show up in Legend, so important to know which are which)
Idents(visu_obj) <- visu_obj@meta.data$traj_res.0.03
cols.use <- list(traj_res.0.03=hue_pal()(length(levels(visu_obj@meta.data$traj_res.0.03))), model=hue_pal()(6))
names(cols.use[['traj_res.0.03']]) <- level(visu_obj@meta.data$traj_res.0.03)
names(cols.use[["model"]]) <- c( "MON_WT_NI", "MON_WT_INF", "MON_79_NI" , "MON_79_INF", "MON_DEL_NI", "MON_DEL_INF")

# You can add additional groups but groups will be ordered by the "group.by" first and foremost
# You can change order of top level groups and genes by script "ordering_HM_&genes.R"
DoMultiBarHeatmap(subset(visu_obj, downsample=500), draw.lines= T , lines.width = 10,
                  features = Used_markers, slot="data", assay = "integrated",
                  disp.min = 0, disp.max = 4,
                  group.by = "traj_res.0.03",
                  additional.group.by = c("model"),
                  cols.use = cols.use)+ #NoLegend()+
  scale_fill_gradientn(colors = c("blue", "yellow", "red"), na.value = "white")
ggsave("20_11_04_SILI_HM_UsedMarkers.pdf", height = 15, width = 15)
