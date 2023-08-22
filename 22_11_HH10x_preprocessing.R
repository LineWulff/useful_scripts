#################################################################################
###### ----------------- pat3 ILF analysis - initial ----------------- ##########
#################################################################################

###### Libraries ######
library(gplots)
library(dplyr)
library(stringr)
library(scales)
library(Rmagic)
library(ggplot2)
library(viridis)
library(clustree)
library(ccRemover)
library(rgl)
#library(Seurat, lib.loc = "/Library/Frameworks/R.framework/Versions/3.5/Resources/other_libs")
library(Seurat)

###### variables used through script ######
rm(list=ls())
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
# project
project <- "HH26-MLN-CD45neg"

# HH23-MLN-CD45neg
# HH23-PP-INF-CD45neg
# HH23-SI-INF-MILF-CD45neg
# HH23-SILP-HEA-CD45neg
# HH23-SILP-INF-CD45neg

#"HH26-SILP-HEA-CD45neg" #"HH26-SILP-INF-CD45neg" #"HH23-SI-INF-MILF-CD45pos"#"HH23-SILP-HEA-HLADR"

setwd(paste("/Volumes/Mucosal-immunology/WA group/Helmsley_Data/22_07_01/",project,sep="/"))

#### read in data ####
pat_data <- Read10X(data.dir = paste(getwd(),"filtered_feature_bc_matrix/", sep="/"))

###### As Seurat objects ######
pat_data <- CreateSeuratObject(counts = pat_data, project = project, min.cells = 3, min.features = 100)

###### ----------------- Quality Control part 1 ----------------- ######
pat_data[["percent.mt"]] <- PercentageFeatureSet(pat_data, pattern = "^MT-")

#### cell cycle 
data("cc.genes") ## Tirosh et al. 2015 HUMAN GENES
sgenes <- cc.genes$s.genes
g2mgenes <- cc.genes$g2m.genes

pat_data <- CellCycleScoring(pat_data, s.features = sgenes, g2m.features = g2mgenes, set.ident = FALSE)


###### ----------------- QC plots ----------------- ######
setwd(paste(getwd(),"QC/",sep="/"))
# nFeature = number of genes per barcode
# nCount = number of reads per barcode = count depth

### sample
png(paste(dato,project,"QC_measures1.png",sep="_"),width = 1200, height = 1000, res = 150)
VlnPlot(pat_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
dev.off()

pat_data_QC <- as.data.frame(pat_data@meta.data)
pat_data_QC <- pat_data_QC[order(pat_data_QC$nCount_RNA),]
pat_data_QC <- cbind(pat_data_QC, rank = c(1:length(pat_data_QC$nCount_RNA)))

## Count depth
#red line = zoom for next plot
png(paste(dato,project,"QC_CountDepthHisto_1.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=nCount_RNA))+geom_histogram(fill="grey",bins = 100)+geom_vline(xintercept = 6000, colour = "red")+
  xlab("Count Depth")+ggtitle(paste(project,"_ILF",sep=""))
dev.off()

# change count depth cut off suggestion
png(paste(dato,project,"QC_CountDepthHisto_2.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC[pat_data_QC$nCount_RNA<6000,], aes(x=nCount_RNA))+geom_histogram(fill="grey", bins = 30)+
  geom_vline(xintercept = 2500, colour = "red")+xlab("Count Depth")+ggtitle(paste(project,"_ILF",sep=""))
dev.off()

## Number of genes
png(paste(dato,project,"QC_NumberGeneshHisto_1.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=nFeature_RNA))+geom_histogram(fill="grey",bins = 100)+
  xlab("Number of genes")+ggtitle(paste(project,"_ILF",sep=""))+geom_vline(xintercept = c(800,4500), colour = "red")
dev.off()

png(paste(dato,project,"QC_NumberGeneshHisto_2.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC[pat_data_QC$nFeature_RNA<2000,], aes(x=nFeature_RNA))+geom_histogram(fill="grey",bins = 100)+
  xlab("Number of genes")+ggtitle(paste(project,"_ILF",sep=""))+geom_vline(xintercept = 800, colour = "red")
dev.off()

png(paste(dato,project,"QC_CountDepthRanked.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=-rank, y=nCount_RNA))+geom_point()+
  xlab("Barcode rank")+ylab("Count Depth")+ggtitle(paste(project,sep=""))+
  geom_hline(yintercept = 2500, colour = "red")
dev.off()

png(paste(dato,project,"QC_CountDepth_GeneCount_1.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=nCount_RNA, y=nFeature_RNA, colour=percent.mt))+geom_point()+
  scale_colour_viridis()+xlab("Count Depth")+ylab("Number of genes")+ggtitle(paste(project,sep=""))+
  geom_vline(xintercept = c(2500,25000), colour = "red")+geom_hline(yintercept = c(1000,4500),colour="red")
dev.off()

png(paste(dato,project,"QC_CountDepth_GeneCount_2.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC[pat_data_QC$nCount_RNA<6000,], aes(x=nCount_RNA, y=nFeature_RNA, colour=percent.mt))+geom_point(size=0.2)+
  scale_colour_viridis()+xlab("Count Depth")+ylab("Number of genes")+ggtitle(paste(project,sep=""))+
  geom_hline(yintercept = 600, colour = "red")
dev.off()

#Set cut off slightly lower than the lines (can be adjusted later)
### Set thresholds and subset data
pat_data <- subset(pat_data, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500 & percent.mt < 10 )#& nCount_RNA < 25000)

### Normalize data
pat_data <- NormalizeData(pat_data)

### save objects separately
pat_data <- FindVariableFeatures(pat_data, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(pat_data), 10)

pat_data <- RenameCells(pat_data, add.cell.id = project)

setwd(paste("/Volumes/Mucosal-immunology/WA group/Helmsley_Data/22_07_01/",project,sep="/"))
saveRDS(object=pat_data, file=paste(project,".rds",sep=""))

plot1 <- VariableFeaturePlot(pat_data)
#png(paste(dato,project,"TopVariableeGenes.png",sep = "_"),height = 1000,width = 1000, res =150)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
#dev.off()

## Add prefix to separate samples from each otehr in integrated analysis

dims = 60

## Human cell cycle genes from ccremover
cell_cycle_indices <- gene_indexer(rownames(pat_data@assays$RNA), species = "human", 
                                   name_type = "symbols" )
human_cell_cycle <- rownames(pat_data@assays$RNA)[cell_cycle_indices]

pat_data[["cc.exp"]] <- colSums(pat_data@assays$RNA@counts[human_cell_cycle, , drop = FALSE])/pat_data[[paste0("nCount_","RNA")]]*100

pat_data <- ScaleData(pat_data, verbose = TRUE, vars.to.regress = c("nCount_RNA", "percent.mt","cc.exp","G2M.Score","S.Score"), do.scale = TRUE)
pat_data <- RunPCA(pat_data, npcs = 60, verbose = FALSE)

print(pat_data[["pca"]], dims = 1:5, nfeatures = 10)

### Deciding PCs to use 
# VizDimLoadings(pat_data, dims = 1:2, reduction = "pca")                               
# ggsave(paste(dato,"VizPCA.png",sep="_"), width = 20, height = 10, units = "cm", dpi = 600, scale=1.5, limitsize = FALSE)
# 
# DimPlot(pat_data, reduction = "pca", group.by = "orig.ident")                              # PCA-plot grouped by sample
# ggsave(paste(dato,"PCA_plot_60PCs.png",sep=""), width = 10, height = 10, units = "cm", dpi = 600, scale=1.7, limitsize = FALSE)
# 
# png(paste(dato,"PCA_heatmaps_60PCs.png",sep = "_"), width = 20, height = 110, units = "cm", res = 600)
# par(mfrow = c(1, 1))
# DimHeatmap(pat_data, dims = 1:dims, cells = 500, balanced = TRUE)                        # PC heatmaps (1-50 PCs)
# dev.off()

# JackStraw analysis
# pat_data <- JackStraw(pat_data, num.replicate = 100, dims=dims)
# pat_data <- ScoreJackStraw(pat_data, dims = 1:dims)
# JackStrawPlot(pat_data, dims = 1:dims)#+theme(legend.position = "none")
# ggsave(paste(dato,"JackStraw_60PCs.png"), width = 30, height = 15, units = "cm", dpi = 600, scale=1.5, limitsize = FALSE)

# Elbow-plot - this can evaluate up to 60PCs in seurat 3
ElbowPlot(pat_data, ndims = 60)
ggsave(paste(dato,"ElbowPlot_60PCs.png",sep = "_"), width = 20, height = 15, units = "cm", dpi = 600, scale=1.5, limitsize = FALSE)

### Select PC for integrated analysis
dims = 20

###### ----------------- Clustering and UMAP dim. reduction ----------------- ######
pat_data <- RunUMAP(pat_data, dims = 1:dims, reduction = "pca")

#png(paste(dato,round,project,"UMAP_nCountRNA.png",sep = "_"),height = 800, width = 800, res =150)
FeaturePlot(pat_data, features = c("nCount_RNA","nFeature_RNA"), reduction = "umap", cols = c("blue","grey","red"))
#dev.off()

## looking at cc influence on clustering
# png(paste(dato,round,project,"UMAP_CCmeasures.png",sep = "_"),height = 2*800, width = 2*800, res =150)
# FeaturePlot(pat_data, features = c("S.Score","G2M.Score","cc.exp"), reduction = "umap", cols = c("grey","red"))
# dev.off()
# png(paste(dato,round,project,"UMAP_CCphase.png",sep = "_"),height = 800, width = 800, res =150)
# DimPlot(pat_data, group.by = "Phase", reduction = "umap")
# dev.off()

#### Identity markers ####
#setwd("/Volumes/sund$/WA group/Tom and Line data/patient_3/08/analysis/ID_markers/")
### MAC markers
MAC_genes <- c("MERTK","CD163","CD14")
#png(paste(dato,round,project,"MACgenes.png",sep = "_"),height = 1000, width = 1000, res =150)
FeaturePlot(pat_data, features = MAC_genes, reduction = "umap", cols = c("grey","blue"),ncol = 2)
#dev.off()

### MO markers (CD14+)
MO_genes <- c("S100A9","FCN1","CD14","VCAN")
#png(paste(dato,round,project,"MOgenes.png",sep = "_"),height = 800, width = 800, res =150)
FeaturePlot(pat_data, features = MO_genes, reduction = "umap", cols = c("grey","blue"),ncol = 2)
#dev.off()

### DC2 genes
DC2_genes <- c("CD1C","CD207","FLT3","IRF4")
#png(paste(dato,round,project,"DC2genes.png",sep = "_"),height = 1000, width = 1000, res =150)
FeaturePlot(pat_data, features = DC2_genes, reduction = "umap", cols = c("grey","blue"),ncol = 2)
#dev.off()

### DC1 genes
DC1_genes <- c("XCR1","CADM1","CLEC9A")
#png(paste(dato,round,project,"DC1genes.png",sep = "_"),height = 1000, width = 1000, res =150)
FeaturePlot(pat_data, features = DC1_genes, reduction = "umap", cols = c("grey","blue"),ncol = 2)
#dev.off()

#### contaminants ####
### B cells 
Bcell_genes <- c("CD79A","CD79B","IGHA1","IGHM")
#png(paste(dato,round,project,"Bcellgenes.png",sep = "_"),height = 1000, width = 1000, res =150)
FeaturePlot(pat_data, features = Bcell_genes, reduction = "umap", cols = c("grey","blue"),ncol = 2)
#dev.off()

### ILC
#png(paste(dato,round,project,"ILCgenes.png",sep = "_"),height = 500, width = 500, res =150)
FeaturePlot(pat_data, features = "KIT", reduction = "umap", cols = c("grey","blue"))
#dev.off()

### Endothelial
End_genes <- c("VWF","SPARC")
#png(paste(dato,round,project,"Endogenes.png",sep = "_"),height = 500, width = 500, res =150)
FeaturePlot(pat_data, features = End_genes, reduction = "umap", cols = c("grey","blue"),ncol = 2)
#dev.off()

### T cells
T_genes <- c("CD3E","CD3D")
#png(paste(dato,round,project,"Endogenes.png",sep = "_"),height = 500, width = 500, res =150)
FeaturePlot(pat_data, features = T_genes, reduction = "umap", cols = c("grey","blue"),ncol = 2)



