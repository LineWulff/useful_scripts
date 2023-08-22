
library(ggplot2); theme_set(theme_classic())

mlCols <- readRDS(url("http://web.stanford.edu/group/butcherlab/mlColors.rds"))
mlColsFunc <- function(x) mlCols[1:x]

scalar1 <- function(x) {x / sqrt(sum(x^2))}

#select 3d projection columns and copy with column names from JMP
# select the three umap columns from JMP and copy them
original_data <- clipr::read_clip_tbl(header=F) 
#select label column in JMP and copy (without column names)
#labs <- clipr::read_clip() 
#cellID <- clipr::read_clip() 

#change to fit script from JMP 
#red arrow copy script take out rotation part
rot <- c(  -146.056926551631, -62.9614994993169, 23.9266478358373 )

rot <- rev(rot) * (-pi/180)
x_axis <- c(1,0,0)
y_axis <- c(0,1,0)
ref.plane <- c(0,0,1)
rot.mat <- matrix(c(cos(rot[1])*cos(rot[2]),
                    sin(rot[1])*cos(rot[2]),
                    -sin(rot[2]),
                    cos(rot[1])*sin(rot[2])*sin(rot[3])-sin(rot[1])*cos(rot[3]),
                    sin(rot[1])*sin(rot[2])*sin(rot[3])+cos(rot[1])*cos(rot[3]),
                    cos(rot[2])*sin(rot[3]),
                    cos(rot[1])*sin(rot[2])*cos(rot[3])+sin(rot[1])*sin(rot[3]),
                    sin(rot[1])*sin(rot[2])*cos(rot[3])-cos(rot[1])*sin(rot[3]),
                    cos(rot[2])*cos(rot[3])), nrow = 3)

doRotation <- function(x) c(rot.mat %*% x)
x_axis.rot <- scalar1(doRotation(x_axis))
y_axis.rot <- scalar1(doRotation(y_axis))
ref.plane.rot <- scalar1(doRotation(ref.plane))

projected_data <- t(apply(original_data, 1, function(x) {
  c(x_axis.rot %*% (x - ref.plane.rot), y_axis.rot %*% (x - ref.plane.rot))
}))



dat.toplot <- as.data.frame(projected_data)
#add to embeddings in seurat object

# dat.toplot$labs <- labs
# rownames(dat.toplot) <- cellID
# ggplot2::ggplot(dat.toplot, aes(V1, V2, color = labs)) + 
#   ggplot2::geom_point(size = 0.5)
# 
# write.csv(dat.toplot, file="DCtraj_Flattended3Dcoord.csv")
# 
# visu_obj@reductions$flat_3D <- visu_obj@reductions$umap
# visu_obj@reductions$flat_3D@cell.embeddings <- as.matrix(dat.toplot[colnames(visu_obj),1:2])
# colnames(visu_obj@reductions$flat_3D@cell.embeddings) <- c("UMAP_1","UMAP_2")
# DimPlot(visu_obj, group.by = "integrated_snn_res.0.7", reduction = "flat_3D", pt.size = 1)+
#   theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
# 
# save
# 


