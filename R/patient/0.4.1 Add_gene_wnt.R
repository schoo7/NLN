rm(list=ls())

library(Seurat)
library(dplyr)
# BiocManager::install("GSEABase")
library(GSEABase)

gset1 <- getGmt("/work/SCCC/s228508/choushi/HALLMARK_WNT_BETA_CATENIN_SIGNALING.v2023.2.Hs.gmt")
gset2 <- getGmt("/work/SCCC/s228508/choushi/GOBP_EXOCYTOSIS.v2023.2.Hs.gmt")
gset3 <- getGmt("/work/SCCC/s228508/choushi/GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS.v2023.2.Hs.gmt")
gset4 <- getGmt("/work/SCCC/s228508/choushi/HALLMARK_MITOTIC_SPINDLE.v2023.2.Hs.gmt")

gset1 = gset1[["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]]@geneIds
gset2 = gset2[["GOBP_EXOCYTOSIS"]]@geneIds
gset3 = gset3[["GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS"]]@geneIds
gset4 = gset4[["HALLMARK_MITOTIC_SPINDLE"]]@geneIds

gene_list <- list()
gene_list[["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]] <- gset1
gene_list[["GOBP_EXOCYTOSIS"]] <- gset2
gene_list[["GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS"]] <- gset3
gene_list[["HALLMARK_MITOTIC_SPINDLE"]] <- gset4

#-------------------------------
seurat.primary <- readRDS("/work/SCCC/s228508/choushi/choushi_sc.log_primary.rds")
all.genes = rownames(seurat.primary)
# write.table(all.genes, "/work/SCCC/s228508/4.1 Add_gene_signature/all_genes_26460.csv",sep = ",")

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


saveRDS(gene_list, "/work/SCCC/s228508/choushi/gene_list_wnt.rds")



a <- read.table("/work/SCCC/s228508/choushi/Dotplot_Grey.gmt.txt", sep = '')



