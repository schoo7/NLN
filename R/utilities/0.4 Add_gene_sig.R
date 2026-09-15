rm(list=ls())

library(Seurat)
library(dplyr)

seurat.primary <- readRDS("/work/SCCC/s228508/choushi/choushi_sc.log_primary.rds")
all.genes = rownames(seurat.primary)
# write.table(all.genes, "/work/SCCC/s228508/4.1 Add_gene_signature/all_genes_26460.csv",sep = ",")

gene_list <- read.table("/work/SCCC/s228508/choushi/NLN KOvsWT 5000.csv", sep = ",")
colnames(gene_list) <- "genes"
gene_list$genes <- toupper(gene_list$genes)

# top_2000 <- head(gene_list$genes , 2000)
# 
# elements_not_in_allgenes <- setdiff(top_2000, all.genes)
# elements_not_in_allgenes
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

gene_list <- lapply(gene_list, remove_empty)
gene_list <- lapply(gene_list, remove_space)
#---------

# # unmatch items
# #----------
# names <- names(T_markers)
# unmatch_items <- c()
# 
# for (i in names) {
#   a <- setdiff(T_markers[[i]], all.genes)
#   unmatch_items <- append(unmatch_items,a)
# }
# 
# unmatch_items <- unique(unmatch_items)
# #----------------

# find the genes not in all_genes
#---------
elements_not_in_allgenes <- list()

for (i in names(gene_list)) {
  a <- gene_list[[i]]
  elements_not_in_allgenes[[i]] <- setdiff(a, all.genes)
}
elements_not_in_allgenes
#---------

# correct list1
# ------------------------
name_correction = list(
  "C1orf228"="ARMH1",
  "C6orf48"="SNHG32",  
  "AES"="TLE5",    
  "SEPT6"="SEPTIN6", 
  "FAM60A"="SINHCAF",  
  "C17orf51"="LINC02693", 
  "FAM46C"="TENT5C",
  "FAM129A"="NIBAN1",
  "TMEM2"="CEMIP2",
  "CCL3L3"="", 
  "RARRES3"="PLAAT4",
  "ATP5E"="ATP5F1E",
  "TMEM2"="CEMIP2",
  "KIAA1551"="RESF1",
  "ATP5L"="ATP5MG",
  "PLA2G16"="PLAAT3",
  "SEPT1"="SEPTIN1",
  "SLC1A7"="",
  "CXCR1"="",     
  "KIR3DX1"="",
  "C10orf128"="TMEM273", 
  "FAM46A"="TENT5A",
  "FAM19A1"="TAFA1",
  "TMEM2"="CEMIP2",
  "SEPT7"="SEPTIN7",
  "FAM173A"="ANTKMT",
  "FAM26F"="CALHM6",
  "PLA2G16"="PLAAT3",
  "SEPT11"="SEPTIN11",
  "ATP5E"="ATP5F1E",   
  "ATP5G3"="ATP5MC3",  
  "USMG5"="",
  "ATP5C1"="ATP5F1C",  
  "ATP5B"="ATP5F1B",   
  "ATP5J"="ATP5PF",   
  "ATP5J2"="ATP5MF",  
  "MINOS1"="MICOS10",
  "ATP5F1"="ATP5PB", 
  "SEPT11"="SEPTIN11",
  "ATPIF1"="ATP5IF1",
  "UNQ6494"="",
  "C16orf45"="BMERB1",
  "LHFP"="LHFPL6",
  "ATP5G3"="ATP5MC3",
  "SEPT2"="SEPTIN2",
  "ASNA1"="GET3",
  "ATP5G1"="ATP5MC1",
  "ATP5A1"="ATP5F1A",
  "ATP5G2"="ATP5MC2",   
  "ATP5D"="ATP5F1D",
  "GGT1"="GGT1.1",    
  "RTFDC1"="RTF2",
  "C19orf24"="FAM174C",
  "FAM96A"="CIAO2A",
  "C17orf62"="CYBC1"
)
gene_list <- lapply(gene_list, rename_gene, name_correction=name_correction)
gene_list <- lapply(gene_list, remove_empty)
gene_list <- lapply(gene_list, remove_space)
#--------------


# # correct list5
# #--------
# name_correction = list(
#   "FAM96B"="",
#   "C20orf24"="",  
# )
# 
# gene_list <- lapply(gene_list, rename_gene, name_correction=name_correction)
# gene_list <- lapply(gene_list, remove_empty)
# gene_list <- lapply(gene_list, remove_space)
# #--------


# remove all unmatched genes directly
#----------
for (j in names(gene_list)) {
  gene_list[[j]] <- gene_list[[j]][!gene_list[[j]] %in% elements_not_in_allgenes[[j]]]
}

# check again
#----------
elements_not_in_allgenes <- list()

for (i in names(gene_list)) {
  a <- gene_list[[i]]
  elements_not_in_allgenes[[i]] <- setdiff(a, all.genes)
}

duplicates <- names(gene_list)[duplicated(names(gene_list))]
duplicates



#--------------------
top_2000 <- head(gene_list$genes , 2000)

top_2000 <- data.frame(top_2000)
colnames(top_2000) <- "genes"
# top_2000 <- lapply(top_2000, remove_empty)
# top_2000 <- lapply(top_2000, remove_space)

write.table(top_2000, "/work/SCCC/s228508/choushi/top2000.csv")










