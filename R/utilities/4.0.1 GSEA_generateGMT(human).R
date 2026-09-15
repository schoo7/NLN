HALLMARK <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/7.0 Myeloid/HALLMARK(240816).rds")
for (i in names(HALLMARK)) {
  HALLMARK[[i]] <- unique(HALLMARK[[i]])
}

kegg <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/7.0 Myeloid/KEGG(240929).rds")
for (i in names(kegg)) {
  kegg[[i]] <- unique(kegg[[i]])
}

reactomes <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/7.0 Myeloid/reactomes(240914).rds")
for (i in names(reactomes)) {
  reactomes[[i]] <- unique(reactomes[[i]])
}

GOALL <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/7.0 Myeloid/GOALL_genesets(240814).rds")
for (i in names(GOALL)) {
  GOALL[[i]] <- unique(GOALL[[i]])
}

combined_list <- c(HALLMARK, kegg, reactomes, GOALL)  

write.gmt <- function(gs, file){
  sink(file)
  lapply(names(gs), function(i){
    cat(paste(c(i,"tmp",gs[[i]]), collapse = "\t"))
    cat("\n")
  })
  sink()
}

write.gmt(combined_list, "/gpfs/gibbs/project/mu/qx75/Grey/choushi/pathways.gmt")

a = read.gmt("/gpfs/gibbs/project/mu/qx75/Grey/choushi/pathways.gmt")
