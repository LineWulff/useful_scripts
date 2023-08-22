
## the function takes a prefix of a gene name or a whole gene name and looks for matches in all teh seurat object assays
## Always read function into environment first, by running from first line of function (here, line 7)
## Example: 
## My seurat object is called MNP and I want to look for gene names starting with CD10
## genenamesinlist("CD10",MNP)
genenameinlist <- function(gene,seu_obj){
  assays <- Assays(seu_obj)
  gene_list <- c(rep(NA, length(assays)))
  names(gene_list) <- assays
  for (ass in assays){
    DefaultAssay(seu_obj) <- ass
    genes <- rownames(seu_obj)[startsWith(rownames(seu_obj),gene)]
    if (length(genes)>0){
      gene_list[ass] <- list(genes) }}
  return(gene_list)
}


