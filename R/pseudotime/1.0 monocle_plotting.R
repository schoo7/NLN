rm(list=ls())

library('monocle3')
library('ggplot2')
library(devtools)
library('dplyr')
library('ggplot2')
library('circlize')
library('RColorBrewer')
library('ComplexHeatmap')

clip = function(array_x, pmin, pmax) {
  pmin = quantile(array_x, pmin)
  pmax = quantile(array_x, pmax)
  array_x[array_x < pmin] = pmin
  array_x[array_x > pmax] = pmax
  return (array_x)
}

setwd('/gpfs/gibbs/project/mu/qx75/Grey/10.0 pseudotime_trajectory/')
# expr = read.table('expr.txt',sep='\t', row.names = 1, header=T)
meta = read.table('cell_meta.txt', sep = '\t', row.names = 1, header=T)
cds = readRDS('monocle.Rdata')
# expr = as.matrix(exprs(cds))
colData(cds)['phase'] = meta$phase
colData(cds)['Mutation_counts'] = clip(meta$Mutation_counts,0.05,0.95)
colData(cds)['ITH'] = clip(meta$ITH,0.05,0.95)

pdf(file='pseudotime.pdf')
plot_cells(
  cds, color_cells_by = "pseudotime", show_trajectory_graph = F,
  label_cell_groups=FALSE,
  label_leaves=FALSE,cell_size = 1)
dev.off()

pdf(file='Phase.pdf')
plot_cells(
  cds, color_cells_by = "phase", show_trajectory_graph = F,
  label_cell_groups=FALSE,
  label_leaves=FALSE,cell_size = 1)
dev.off()

pdf(file='Mutation_counts.pdf')
plot_cells(
  cds, color_cells_by = "Mutation_counts", show_trajectory_graph = F,
  label_cell_groups=FALSE,
  label_leaves=FALSE,cell_size = 1.5, alpha=0.75)
dev.off()

pdf(file='ITH.pdf')
plot_cells(
  cds, color_cells_by = "ITH", show_trajectory_graph = F,
  label_cell_groups=FALSE,
  label_leaves=FALSE,cell_size = 1, alpha=0.75)
dev.off()

pdf(file='Monocle_cluster.pdf')
plot_cells(
  cds, color_cells_by = "leiden", show_trajectory_graph = F,
  label_cell_groups=FALSE,
  label_leaves=FALSE,cell_size = 1, alpha=0.75)
dev.off()

########## custom plots ##########

expr = as.matrix(exprs(cds))
for (i in 1:2) {
  if (i == 1) {
    gl_table = read.table('PCA scores alternative.csv',sep=',',header = T)
    gl_table = gl_table[,1:6]
  } else {
    gl_table = read.table('downstream_genes.csv',sep=',', header=T)
    gl_table = gl_table[,-6]
  }
  cluster_colors = setNames(
    brewer.pal(name="Set1", n=8), sort(unique(meta$leiden))
  )
  pseudotime_colors = colorRamp2(
    seq(0, 15, length = 9), rev(brewer.pal(name="YlGnBu", n=9)))
  topanno = HeatmapAnnotation(
    Pseudotime = pseudotime(cds)[order(pseudotime(cds))],
    Cluster = meta$leiden[order(pseudotime(cds))],
    col = list(Pseudotime=pseudotime_colors,Cluster=cluster_colors)
  )
  for (gl in colnames(gl_table)) {
    genes = unique(gl_table[gl][,1])
    genes = genes[genes %in% row.names(expr)]
    if (length(genes) == 0) {next}
    plot_data = expr[genes, order(pseudotime(cds))]
    plot_data = t(apply(plot_data,1,function(x){smooth.spline(x,df=3)$y}))
    plot_data = t(apply(plot_data,1,function(x){(x-mean(x))/sd(x)}))
    psut = pseudotime(cds)/max(pseudotime(cds))
    pdf(paste('pseudotime_', gl, '_top.pdf',sep=''), height=length(genes)/4)
    g = Heatmap(
      plot_data,
      col = colorRamp2(seq(from=-2,to=2,length=9),rev(brewer.pal(9, "RdBu"))),
      show_row_names = TRUE,
      show_column_names  = FALSE,
      row_title_rot = 0,
      top_annotation  = topanno,
      cluster_row_slices = FALSE,
      cluster_columns = FALSE,
      cluster_row = T)
    print(g)
    dev.off()
  }
}


#--------
all_genes = row.names(cds)
invalid_genes = grepl('^MT\\.', all_genes) | grepl('^RPS', all_genes) | grepl('^RPL', all_genes)
invalid_genes = all_genes[invalid_genes]
valid_genes = all_genes[!all_genes %in% invalid_genes]

corr_with_pt = cor(pseudo_time, as.matrix(t(expr[valid_genes,])))
top50 = order(corr_with_pt, decreasing = T)[1:50]
top50 = valid_genes[top50]
bot50 = order(corr_with_pt)[1:50]
bot50 = valid_genes[bot50]

genes = c(bot50,top50)
plot_data = expr[genes, order(pseudotime(cds))]
plot_data = t(apply(plot_data,1,function(x){smooth.spline(x,df=3)$y}))
plot_data = t(apply(plot_data,1,function(x){(x-mean(x))/sd(x)}))
psut = pseudotime(cds)/max(pseudotime(cds))
topanno = HeatmapAnnotation(Pseudotime=pseudotime(cds)[order(pseudotime(cds))], Cluster=meta$leiden, Sample=meta$Sample)
pdf('pseudotime_top.pdf',height = 15)
Heatmap(
  plot_data,
  col= colorRamp2(seq(from=-2,to=2,length=9),rev(brewer.pal(9, "YlGnBu"))),
  show_row_names  = TRUE,
  show_column_names = FALSE,
  row_title_rot = 0,
  top_annotation  = topanno,
  cluster_row_slices  = FALSE,
  cluster_columns = FALSE,)
dev.off()
