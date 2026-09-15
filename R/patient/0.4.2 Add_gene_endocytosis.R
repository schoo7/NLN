rm(list=ls())

library(Seurat)
library(dplyr)

# seurat.primary <- readRDS("/work/SCCC/s228508/choushi/choushi_sc.log_primary.rds")
# all.genes = rownames(seurat.primary)
# write.table(all.genes, "/work/SCCC/s228508/4.1 Add_gene_signature/all_genes_26460.csv",sep = ",")

# get gene set
geneset <- read.table("/work/SCCC/s228508/choushi/Endocytosis.csv", sep = ",")
colname<- geneset[1,]
colnames(geneset) <- colname
geneset <- geneset[-1,]

# function
#-----------------
remove_empty <- function(x) x[x != ""]
remove_space <- function(x) gsub(" ", "", x, fixed=TRUE)
rename_gene <- function(x, name_correction){
  for(d in 1 : length(x)){
    corrected_name = name_correction[[x[d]]]
    if(! is.null(corrected_name)){
      x[d] = corrected_name
    }else if(toupper(x[d]) %in% all.genes){
      x[d] = toupper(x[d])
    }
  }
  return(x)
}
toupper_gene <- function(geneset){
  for (i in colnames(geneset)) {
    geneset$i <- toupper(geneset$i)
  }
  return(geneset)
}

#-----------------
geneset <- lapply(geneset, remove_empty)
geneset <- lapply(geneset, remove_space)
geneset <- lapply(geneset, toupper_gene)

geneset

saveRDS(geneset, "/work/SCCC/s228508/choushi/Endocytosis.rds")










