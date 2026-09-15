library(Seurat)
endocytosis <- readRDS("/work/SCCC/s228508/choushi/Endocytosis.rds")

selected <- c("GOCC_CLATHRIN_SCULPTED_MONOAMINE_TRANSPORT_VESICLE",
              "GOBP_VESICLE_TRANSPORT_ALONG_MICROTUBULE",
              "GOBP_EXOCYTIC_PROCESS",
              "GOBP_POSITIVE_REGULATION_OF_EXOSOMAL_SECRETION",
              "GOBP_CLATHRIN_COATED_VESICLE_CARGO_LOADING",
              "GOCC_COPI_VESICLE_COAT",
              "REACTOME_COPII_MEDIATED_VESICLE_TRANSPORT")

endocytosis <- endocytosis[selected]

write.gmt <- function(gs, file){
  sink(file)
  lapply(names(gs), function(i){
    cat(paste(c(i,"tmp",gs[[i]]), collapse = "\t"))
    cat("\n")
  })
  sink()
}

write.gmt(geneSets, "/gpfs/gibbs/project/mu/qx75/Grey/5.2.1 CD4/HALLMARK.gmt")

a = read.gmt("/gpfs/gibbs/project/mu/qx75/Grey/5.2.1 CD4/HALLMARK.gmt")


