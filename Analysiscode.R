#download (install.packages('')) and load (library ('') following packages
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(sctransform)
library(limma)
library(VennDiagram)


#load matrix files created through aggr in cell ranger without depth normalization with file output called filtered_feature_bc_matrix
matrix_dir = "/Users/rajabn/Desktop/SSAnalysis/Merged_Aggr_Nonorm_Filteredfeature/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2 #ensure feature names V refers to the right V column so that gene IDs are picked up instead of ensembl IDs (important for removing mitochondrial genes later on)

#create a seurat object. .
seurat_object <- CreateSeuratObject(counts = mat, project = "singlecell", min.cells = 3, min.features = 200)
seurat_object #to see what this consists on
dense.size <- object.size(as.matrix(mat))
dense.size
sparse.size <- object.size(mat)
sparse.size
dense.size/sparse.size
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-") #lower case mt if for mouse mitochondrial genes. capital for human.
head(seurat_object@meta.data, 20) #shows QC metrics for the first 20 cells
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#remove unwanted cells. below are default settings but you can modify these
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#now you have removed unwanted cells, it is time to normalize the data. By default, Seurat employs a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
#you can alternatively input seurat_object <- NormalizeData(seurat_object) instead of above.
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#we now identify highly variable features in order to determine a subset of features that exhibit high cell-to-cell variation in the dataset.
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0, max.overlaps=20)
plot2
#now we apply a linear transformation (scaling) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
#we can visualise both cells and features that define the PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
print(seurat_object[["pca"]], dims = 1:8, nfeatures = 2)
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
DimPlot(seurat_object, reduction = "pca")
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE)
#in order to determine the dimensionality of the dataset, we can first visualise the distribution of pvalues for each PC with a uniform distribution. Significant PCs will show a strong enrichment of features with low pvalues.
#this can be visualised using an Elbow Plot
ElbowPlot(seurat_object)
#to cluster cells, I have used the number of sig. PCs that I observed in the above plots. The findneighbors function is for constuction of a KNN graph based on euclidean distance in PCA space and refines edge weights between any two cells based on the shared overlap in their local neighborhoods (jaccard similarity). It uses the input of previously defined dimensionality of the dataset.
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
#now to actually cluster the cells, we apply modularity optimisation techniques (default is Louvain algorithm). The findclusters function contains a resolution parameter which sets the granularity of downstream clustering. Settings are recommended between 0.4-1.2, but may need to be increased for larger datasets.
seurat_object <- FindClusters(seurat_object, resolution =0.90)
head(Idents(seurat_object), 5) #to have a look at cluster IDs of first 5 cells.

#run non-linear dimensional reduction (UMAP/tSNE)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "umap")
cell.num2 <- table(seurat_object@active.ident)
#To colour with different colours:
p<-DimPlot(seurat_object, reduction = "umap")
p<-p + scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))
p 

#to seperate out barcode.names file to identify samples
library(tidyr)
seurat_object@meta.data$names <- row.names(seurat_object@meta.data) #makes another column that has the row names of the data frame
groups <- separate(seurat_object@meta.data, col = names, into = c("a", "b"), sep = "-")
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
seurat_object@meta.data$group = groups$b #added group column to seurat_object meta.data
DimPlot(seurat_object, reduction = "umap", group.by = "group")


#colour UMAP by library size or my gene expession
FeaturePlot(seurat_object, features = "nCount_RNA")
FeaturePlot(seurat_object, features = "MKI67")


# Rename groups
sample_names = c("PSCM control", "PSCM acute", "PSCM 18hr", "PSCM re-stimulation", "MDM control", "MDM acute", "MDM 18hr", "MDM re-stimulation")
sample_labels = as.factor(seurat_object@meta.data$group)
levels(sample_labels) = sample_names
seurat_object@meta.data$sample_labels = sample_labels
cell.num <- table(seurat_object@meta.data$group)
cell.num2 <- table(seurat_object)

#FindMarkers

#find all markers distinguishing cluster 8 from cluster 11
seurat_object.PSCMcontrolsubsets <- FindMarkers(seurat_object, ident.1= 8, ident.2 = 11, min.pct=0.25)
head(seurat_object.PSCMcontrolsubsets, n=5)


#find all markers distinguishing cluster 0 from 10
seurat_object.PSCMcontrolsubsets <- FindMarkers(seurat_object, ident.1= 0, ident.2 = 10, min.pct=0.25)
head(seurat_object.PSCMcontrolsubsets, n=5)
write.table(seurat_object.PSCMcontrolsubsets, file = "C:/Users/rajabn/Desktop/SSAnalysis/PSCM_acutesubsets.tsv", quote = FALSE, sep = "\t" )

# find markers for every cluster compared to all remaining cells, reporting only the positive ones (number of genes shown can be changed)
seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.30, logfc.threshold = 0.25)
seurat_object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)


# Population and gene expression visualisation
no <- DotPlot(object = seurat_object,
              features = c("TNF", "IL1B", "CCL4L2", "CCL3"),
              group.by = 'group', scale = FALSE)


ggplot(no$data, aes(fill = avg.exp.scaled, y = pct.exp, x = id)) +
  geom_bar(stat = "identity") + facet_wrap(~features.plot) +
  theme_bw() + labs(y =  "Percentage of cells (%)", x = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradient(low = "lightgrey", high = "blue") + scale_x_discrete(labels = sample_names)


# Reactome IFN signalling pathway
library(msigdbr)
library(magrittr)

msigdb <- msigdbr(species = "Homo sapiens")

msigdb.split <- msigdb %>%
  split(f = paste(.$gs_cat, .$gs_subcat, sep = "_")) %>%
  set_names(sub("_$", "", names(.)))

# The following splits each of the MSigDB gene set collections into list by their Entrez gene ID
# If your gene annotation is Ensembl, you could change the x argument in split to:
# x = collection$ensembl_gene
# Or, if your annotation is by gene symbols:
# x = collection$gene_symbol

msigdb.lists <- lapply(msigdb.split, function(collection) split(x = collection$gene_symbol, f = collection$gs_name))

isg.I <- unique(msigdb.lists$H$HALLMARK_INTERFERON_ALPHA_RESPONSE)
isg.II <- unique(msigdb.lists$H$HALLMARK_INTERFERON_GAMMA_RESPONSE)

isg.reactome <- unique(msigdb.lists$`C2_CP:REACTOME`$REACTOME_INTERFERON_SIGNALING)

seurat_reactome <- AddModuleScore(seurat_object, features = list(isg.reactome), name = "Scores_Reactome_")

boxplot(seurat_reactome@meta.data$Scores_Reactome_1 ~ seurat_object@meta.data$group,
        main = "Reactome IFN signalling", ylab = "Module scores", xlab = "", names = levels(sample_labels),
        las = 2)


#heatmaps
group_cluster_table <- table(sample_labels, seurat_object$seurat_clusters)
pheatmap::pheatmap(group_cluster_table, display_numbers = group_cluster_table,
                   annotation_col = data.frame(Cluster = colnames(group_cluster_table), row.names = colnames(group_cluster_table)),
                   annotation_legend = FALSE,
                   annotation_colors = list(Cluster = setNames(c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                                                                 '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
                                                                 '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'),
                                                               colnames(group_cluster_table))),
                                          cluster_rows = FALSE,
                            color = c("#FFFFFF", colorRampPalette(RColorBrewer::brewer.pal(5, "Reds")[-5])(99)))

seurat_object$clustered_groups <- factor(seurat_object$seurat_clusters, levels = levels(seurat_object$seurat_clusters)[new.order])

heat <- DoHeatmap(
  seurat_object,
  features = top2$gene,
  cells = sample(rownames(seurat_object@meta.data), size = 10000, replace = FALSE),
  group.bar = TRUE,
  group.by = "clustered_groups",
  group.colors =  c('#E6194B', '#3CB44B', '#FFE119', '#4363D8', '#F58231', '#911EB4', '#46F0F0', '#F032E6', '#BCF60C', '#FABEBE', '#008080', '#E6BEFF', '#9A6324', '#FFFAC8', '#800000', '#AAFFC3', '#808000', '#FFD8B1', '#000075', '#808080', '#FFFFFF', '#000000')[new.order],
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)+ scale_fill_gradientn(colors = c("blue", "white", "red"))

 #Average Expression of components of TLR pathway
gene_list_TLR4_Dep_Indp <- c("TRAM1", "TICAM1","TRAF3" etc.)
gene_avg_expression <- AverageExpression(seurat_object, slot= "data", features = gene_list_TLR4_Dep_Indp, group.by = "seurat_clusters")


#aggregation based on cluster name
df <- data.frame(row.names = 0:33537)
for (i_cluster in unique(seurat_object@meta.data$sample_labels)) 
{
  message(sprintf("Current cluster:%s\n", i_cluster))
  seurat_object@meta.data[seurat_object@meta.data[,"sample_labels"] == i_cluster,]
  temp <- mat[,row.names(seurat_object@meta.data[seurat_object@meta.data[,"sample_labels"] == i_cluster,])]
  print(dim(temp))
  colSums(temp!=0)
  colnames(temp) <- sample.int(floor(ncol(temp)/100), ncol(temp), replace = TRUE)
  print(colnames(temp))
  for (i_subset in unique(colnames(temp))) 
  {temp_subcluster <- temp[,colnames(temp)==i_subset]
  i_name<-paste(i_cluster, i_subset, sep= "_")
  print(i_name)
  df[i_name] <- rowSums(temp_subcluster)
  }  
}
rownames(df)<-feature.names$V1
colSums(df)
#export file now for projection on atlas
write.table(df, file = "C:/Users/rajabn/Desktop/SSAnalysis/merged_samplenames_100_new.tsv", quote = FALSE, sep = "\t" )
#make sample description file
df_samples <- data.frame(row.names = colnames(df))
df_samples$names <- row.names(df_samples)
df_samples <- separate(df_samples, col = names, into = c("cluster", "sub_cluster"), sep = "_")


#setting up a loop for aggregation for projection
#we want to take for each cluster within the curly braces we want to select out the count matrix that corresponds to that count cluster (everything you want to do is in the braces)
df <- data.frame(row.names = 0:33537)
for (i_cluster in unique(seurat_object@meta.data$seurat_clusters)) 
{
  message(sprintf("Current cluster:%s\n", i_cluster))
  seurat_object@meta.data[seurat_object@meta.data[,"seurat_clusters"] == i_cluster,]
  temp <- mat[,row.names(seurat_object@meta.data[seurat_object@meta.data[,"seurat_clusters"] == i_cluster,])]
  print(dim(temp))
  colSums(temp!=0)
  colnames(temp) <- sample.int(floor(ncol(temp)/150), ncol(temp), replace = TRUE)
  print(colnames(temp))
  for (i_subset in unique(colnames(temp))) 
  {temp_subcluster <- temp[,colnames(temp)==i_subset]
  i_name<-paste(i_cluster, i_subset, sep= "_")
  print(i_name)
  df[i_name] <- rowSums(temp_subcluster)
  }  
}
rownames(df)<-feature.names$V1
colSums(df)
#export file now for projection on atlas
write.table(df, file = "C:/Users/nraja/Desktop/Single Cell processing via Seurat/df.tsv", quote = FALSE, sep = "\t" )
#make sample description file
df_samples <- data.frame(row.names = colnames(df))
df_samples$names <- row.names(df_samples)
df_samples <- separate(df_samples, col = names, into = c("cluster", "sub_cluster"), sep = "_")


#aggregation based on cluster name
df <- data.frame(row.names = 0:33537)
for (i_cluster in unique(seurat_object@meta.data$group)) 
{
  message(sprintf("Current cluster:%s\n", i_cluster))
  seurat_object@meta.data[seurat_object@meta.data[,"group"] == i_cluster,]
  temp <- mat[,row.names(seurat_object@meta.data[seurat_object@meta.data[,"group"] == i_cluster,])]
  print(dim(temp))
  colSums(temp!=0)
  colnames(temp) <- sample.int(floor(ncol(temp)/100), ncol(temp), replace = TRUE)
  print(colnames(temp))
  for (i_subset in unique(colnames(temp))) 
  {temp_subcluster <- temp[,colnames(temp)==i_subset]
  i_name<-paste(i_cluster, i_subset, sep= "_")
  print(i_name)
  df[i_name] <- rowSums(temp_subcluster)
  }  
}
rownames(df)<-feature.names$V1
colSums(df)

#make sample description file
df_samples <- data.frame(row.names = colnames(df))
df_samples$names <- row.names(df_samples)
df_samples <- separate(df_samples, col = names, into = c("cluster", "sub_cluster"), sep = "_")

