library(Seurat)
library(SeuratObject)
library(dplyr)
library(biomaRt)  
library(Matrix)   
library(Matrix.utils) 

rna_assay <- homo[["RNA"]]

layers <- Layers(rna_assay)
print(layers)

counts_layers <- grep("^counts", layers, value = TRUE)
print(counts_layers)

gene_ensembl <- rownames(rna_assay)
mart_mouse<- useDataset("mmusculus_gene_ensembl", mart) 
mart_human<- useDataset("hsapiens_gene_ensembl", mart)

human_genes_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), ## mgi_symbol
                       filters = "ensembl_gene_id",
                       values = gene_ensembl,
                       mart = mart_human)

mouse_genes_mapping <- getBM(attributes = c("mgi_symbol", "external_gene_name"),
                             filters = "external_gene_name",
                             values = external_gene_name,
                             mart = mart_mouse)

gene_map <- setNames(genes_annot$mgi_symbol, genes_annot$external_gene_name)
gene_map<-setNames(mouse_genes_mapping$mgi_symbol,mouse_genes_mapping$external_gene_name)
gene_map<-gene_id_map
counts_list <- list()

for (layer_name in counts_layers) {
  cat("processing layer：", layer_name, "\n")
  
  counts_layer <- GetAssayData(object = rna_assay, layer = layer_name)
  
  new_gene_names <- gene_map[rownames(counts_layer)]
  
  new_gene_names[is.na(new_gene_names) | new_gene_names == ""] <- rownames(counts_layer)[is.na(new_gene_names) | new_gene_names == ""]
  
  rownames(counts_layer) <- new_gene_names
  
  counts_layer_aggregated <- aggregate.Matrix(counts_layer, groupings = rownames(counts_layer), fun = "sum")
  
  colnames(counts_layer_aggregated) <- paste(layer_name, colnames(counts_layer_aggregated), sep = "_")
  
  counts_list[[layer_name]] <- counts_layer_aggregated
}

counts_combined <- Reduce(Matrix::cbind2, counts_list)

new_seurat <- CreateSeuratObject(counts = counts_combined, meta.data = filtered_seurat@meta.data)

#standardization
new_seurat <- NormalizeData(new_seurat)
new_seurat <- FindVariableFeatures(new_seurat, selection.method = "vst", nfeatures = 3000)
new_seurat <- ScaleData(new_seurat)
# PCA, Top 20 principal components
new_seurat <- RunPCA(new_seurat)
head(rownames(new_seurat))
saveRDS(new_seurat, file = "datasets_seurat.rds")

que<-new_seurat
ref<-readRDS("Oligodendrocyte.rds")
DefaultAssay(ref) <- "RNA"
anchors <- FindTransferAnchors(reference = ref, query = que,
                               dims = 1:30, reference.reduction = "pca")

# get predicted ID
predictions <- TransferData(anchorset = anchors, refdata = ref$new_celltype,
                            dims = 1:15)
# add metadata
que<- AddMetaData(que,metadata = predictions)

# UMAP projection
ref <- RunUMAP(ref, dims = 1:15, reduction = "pca", return.model = TRUE)

small <- MapQuery(anchorset = anchors, reference = ref, query = que,
                  #refdata = list(seurat_clusters = "seurat_clusters"),  #可以映射多个标签
                  refdata = list(seurat_clusters = "seurat_clusters", celltype="new_celltype"), 
                  reference.reduction = "pca", 
                  reduction.model = "umap")
# MapQuery() is a wrapper around three functions: TransferData(), IntegrateEmbeddings(), and ProjectUMAP().

p1 <- DimPlot(ref, reduction = "umap", group.by = "new_celltype", label = FALSE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference")+xlim(-7,12)+ylim(-7,15)+theme(
                plot.title = element_text(face = "bold", size = 14), 
                axis.title = element_text(face = "bold", size = 12), 
                axis.text = element_text(face = "bold", size = 10)
              )
p2 <- DimPlot(small, reduction = "ref.umap", group.by = "predicted.id", raster=FALSE,label =FALSE,label.size = 3, repel = TRUE) +  
  ggtitle("A53T_mapping") +
  xlim(-7, 12) +
  ylim(-7, 15) +
  theme(
    plot.title = element_text(face = "bold", size = 14),  
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 10) 
  )
p1 + p2
