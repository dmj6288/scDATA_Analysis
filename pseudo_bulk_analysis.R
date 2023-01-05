library(tidyverse)
library(edgeR)
library(DESeq2)
library(Seurat)
library(EnhancedVolcano)
library(ggplot2)
library(readr)

GSE127774_ACC_seurat <- readRDS("D:/Yi/WINTER 2022/Pseudobulk Analyses/Krameeva Human ACC region pseudobulk analysis scripts/data/GSE127774_ACC_seurat.rds")

GSE127774_ACC_seurat_human <- subset(x = GSE127774_ACC_seurat, subset = orig.ident == "H")

GSE127774_ACC_seurat_human <- FindNeighbors(GSE127774_ACC_seurat_human, dims = 1:10)
GSE127774_ACC_seurat_human <- FindClusters(GSE127774_ACC_seurat_human, resolution = 0.5)

pbmc.markers <- FindAllMarkers(GSE127774_ACC_seurat_human, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.15)
new.cluster.ids <- c("Excitatory Neurons", "Astrocyte", "Inhibitory Neurons", "Inhibitory Neurons", "Inhibitory Neurons", "Microglia", "Oligodendrocytes", "OPC", "Excitatory Neurons", "Inhibitory Neurons", "Inhibitory Neurons", "Astrocyte")
names(new.cluster.ids) <- levels(GSE127774_ACC_seurat_human)
GSE127774_ACC_seurat_human <- RenameIdents(GSE127774_ACC_seurat_human, new.cluster.ids)
DimPlot(GSE127774_ACC_seurat_human, reduction = 'tsne', label = TRUE, pt.size = 1, repel = TRUE, label.size = 6)


location_vector = c()
index = 0
for (tube_number in 1:5){
  for (barcode in GSE127774_ACC_seurat_human@assays$RNA@counts@Dimnames[[2]]){
    if (grepl(as.character(tube_number), barcode, fixed = TRUE) == TRUE){
      index = index + 1
    }
  }
  location_vector <- c(location_vector, index)
}

require(gtools)
permutations(n = 9, r = 3, v = 1:9)

replicate1 = GSE127774_ACC_seurat_human@assays$RNA@data[, 1:location_vector[1]]
replicate2 = GSE127774_ACC_seurat_human@assays$RNA@data[, (location_vector[1]+1):location_vector[3]]
replicate3 = GSE127774_ACC_seurat_human@assays$RNA@data[, (location_vector[3]+1):location_vector[5]]

gene_list <- read_csv("myFile.csv", col_names = FALSE)

rownames(replicate1) <- gene_list[[2]]
rownames(replicate2) <- gene_list[[2]]
rownames(replicate3) <- gene_list[[2]]

cell_identity_replicate1 = GSE127774_ACC_seurat_human@active.ident[1:location_vector[1]]
cell_identity_replicate2 = GSE127774_ACC_seurat_human@active.ident[(location_vector[1]+1):location_vector[3]]
cell_identity_replicate3 = GSE127774_ACC_seurat_human@active.ident[(location_vector[3]+1):location_vector[5]]

cell_names <- levels(cell_identity_replicate1)
comparison_vector <- t(permutations(n = length(cell_names), r = 2, v = cell_names))
#comparison_vector <- combn(cell_names, 2)

levels(cell_identity_replicate1) <- paste(levels(GSE127774_ACC_seurat_human@active.ident[1:location_vector[1]]), c('-1', '-1', '-1', '-1', '-1', '-1'))
levels(cell_identity_replicate2) <- paste(levels(GSE127774_ACC_seurat_human@active.ident[(location_vector[1]+1):location_vector[3]]), c('-2', '-2', '-2', '-2', '-2', '-2'))
levels(cell_identity_replicate3) <- paste(levels(GSE127774_ACC_seurat_human@active.ident[(location_vector[3]+1):location_vector[5]]), c('-3', '-3', '-3', '-3', '-3', '-3'))

colnames(replicate1) <- cell_identity_replicate1
colnames(replicate2) <- cell_identity_replicate2
colnames(replicate3) <- cell_identity_replicate3

agg1 <- t(rowsum(t(replicate1), group = colnames(replicate1), na.rm = T))
agg2 <- t(rowsum(t(replicate2), group = colnames(replicate2), na.rm = T))
agg3 <- t(rowsum(t(replicate3), group = colnames(replicate3), na.rm = T))

for(i in 1:2) {       # for-loop over columns
  column_names <- c(paste(comparison_vector[, i][1], "-1"), paste(comparison_vector[, i][1], "-2"), paste(comparison_vector[, i][1], "-3"), paste(comparison_vector[, i][2], "-1"), paste(comparison_vector[, i][2], "-2"), paste(comparison_vector[, i][2], "-3"))
  count_matrix <- cbind(agg1[, column_names[1]], agg2[, column_names[2]], agg3[, column_names[3]], agg1[, column_names[4]], agg2[, column_names[5]], agg3[, column_names[6]])
  colnames(count_matrix) <- column_names
  coldata <- cbind(column_names, c(comparison_vector[,i][1], comparison_vector[,i][1], comparison_vector[,i][1], comparison_vector[,i][2], comparison_vector[,i][2], comparison_vector[,i][2]))
  colnames(coldata) <- c('replicates', 'cell_type')
  
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix), colData = coldata, design = ~ cell_type)
  dds <- DESeq(dds)
  res <- na.omit(results(dds))
  
  print(res)
  
  p <- EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                    title = paste(comparison_vector[, i][2], " vs ",comparison_vector[, i][1]),
                  pCutoff = 10e-4,
                  FCcutoff = 5,
                  pointSize = 1.5,
                  labSize = 6.0,
                  shape = 8,
                  colAlpha = 1)
  
  p + xlim(-8, 8) + ylim(0, 45)
  
  ggsave(filename = paste(comparison_vector[, i][2], " vs ",comparison_vector[, i][1], ".png"),
         width = 10, height = 10, dpi = 300, units = "in", device='png')
  
}

