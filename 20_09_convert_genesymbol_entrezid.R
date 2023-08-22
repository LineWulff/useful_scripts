########################################################################################################
######################### Converting gene symbols to entrez IDs ########################################
########################################################################################################
#load your DEG list of interest
DEG <- read.csv("/FileWGeneList.csv")
#library for musmusculus
library(org.Mm.eg.db)

#Picking out the genelist from the inputted DEG list
my_fave_genes <- DEG$X
#how many genes were in your DEG list?
length(my_fave_genes)

#conversion function
my_fave_genes_conv <- select(org.Mm.eg.db, my_fave_genes, columns = "ENTREZID", "SYMBOL")
#which genes can be converted - and how many is that?
entrezids <- my_fave_genes_conv[!is.na(my_fave_genes_conv$ENTREZID),]
entrezids$SYMBOL
#how many can be converted?
dim(entrezids)[1]

#which genes cannot be converted?
no_entrezid <- my_fave_genes_conv[is.na(my_fave_genes_conv$ENTREZID),]$SYMBOL
no_entrezid
#how many can not be converted?
length(no_entrezid)

#save those which can be converted with their associated logfoldchange in a csv file
DEGwEntrezids <- DEG[which(DEG$X %in% entrezids$SYMBOL),]
DEGwEntrezids <- cbind(DEGwEntrezids, entrezids)
#DEGwEntrezids <- DEGwEntrezids[,c("ENTREZID","avg_logFC", "X")] #remove hashtag from this line if you want to see associated gene name and not just entrezID
DEGwEntrezids <- DEGwEntrezids[,c("ENTREZID","avg_logFC")] #if you remove hashtag fron above put one in the start of this line

#save the file
write.csv(DEGwEntrezids, file="/FileWGeneList_EntrezIDconverted.csv")





