
library(tidyverse)

library(Seurat)

# load metadata and count matrix
hum_meta <- read.csv("AIBS_human_meta_mini.csv", row.names = 1)
hum_counts <- read.csv("AIBS_human_counts_mini.csv", row.names = 1)
row.names(hum_meta) <- hum_meta$sample_name #necessary for matching between counts and metadata in Seurat object

# create Seurat object
Seu_hum <- CreateSeuratObject(counts = t(hum_counts), meta.data = hum_meta) # t(hum_counts) transposes so that we have rows as genes and columns as samples as they should be for Seurat

rm(hum_counts, hum_meta) # removing large matrices now that we have Seurat obj

### How to accses seurat's identity class
Idents(Seu_hum) #interacts with active.ident part of seurat obj, by default without providing orig.ident, pulls string before _
Seu_hum@active.ident #same as above

### Processing Seurat object

#### Normalization
Seu_hum <- NormalizeData(Seu_hum, normalization.method = "LogNormalize", scale.factor = 1000000) #changing scale.factor to mil so we get CPM

#look at normalized data
Seu_hum@assays$RNA@data 

#### Find variable features
Seu_hum <- FindVariableFeatures(Seu_hum, selection.method = "vst", nfeatures = 2000) #should see effect of changing nfeatures

#look at most variable features
Seu_hum@assays$RNA@var.features 

#### Scale data
Seu_hum <- ScaleData(Seu_hum, verbose = FALSE)

#### Run Principal Component Analysis (PCA)
Seu_hum <- RunPCA(Seu_hum, npcs = 50, verbose = FALSE) #50 is default, we could choose something smaller based on ElbowPlot below

ElbowPlot(Seu_hum, ndims=50) #see SD of each PC, shows how much explained, use to see how many PC needed to best explain data
#cut at the elbow (can argue where cutoff is, might choose 7 or 20)

#### Find neighbors
Seu_hum <- FindNeighbors(Seu_hum, reduction = "pca", dims = 1:20)

#### Find clusters
Seu_hum <- FindClusters(Seu_hum, resolution = 0.5) #nm.method and annoy.metric have drastic effects on cluster creation

table(Seu_hum$seurat_clusters) #tells you number of cells in each cluster
table(Seu_hum$seurat_clusters, Seu_hum$class_label) #number of cells per class per cluster
table(Seu_hum$seurat_clusters, Seu_hum$subclass_label) #number of cells per subclass per cluster

#### Run UMAP
Seu_hum <- RunUMAP(Seu_hum, reduction = "pca", dims = 1:20)

# visualizing clusters
DimPlot(Seu_hum, reduction = "umap", group.by = "subclass_label", label=TRUE)
DimPlot(Seu_hum, reduction = "umap", group.by = "seurat_clusters", label=TRUE, repel=TRUE)

### Differential expression analysis
#### Finding cluster marker genes
Idents(Seu_hum) <- "seurat_clusters" #set active identity to our newly defined clusters, these are what we want to compare in the following steps

# get genes to distinguish cluster 8 from 3
cluster_8_v_3 <- FindMarkers(Seu_hum, ident.1 = 8, ident.2 = 3, logfc.threshold = log(2), min.pct = 0.50)

# genes to distinguish cluster 8, 2, or 13 from all other clusters
cluster_8_v_all <- FindMarkers(Seu_hum, ident.1 = 8, logfc.threshold = log(2), min.pct = 0.50)
cluster_2_v_all <- FindMarkers(Seu_hum, ident.1 = 2, logfc.threshold = log(2), min.pct = 0.50)
cluster_13_v_all <- FindMarkers(Seu_hum, ident.1 = 2, logfc.threshold = log(2), min.pct = 0.50)

# finding marker genes for all clusters, alternative to running one-by-one as done for 8, 2, 13 (takes a while to run)
all_clusters <- FindAllMarkers(Seu_hum, logfc.threshold = log(2), min.pct = 0.50)

# taking a look at most significantly DE marker genes
cluster_8_v_3  %>% 
  arrange(p_val_adj)
cluster_8_v_all %>% 
  arrange(p_val_adj)
cluster_2_v_all %>% 
  arrange(p_val_adj)

# two different ways of looking for a specific gene in the output dataframes, where genes are rownames
cluster_2_v_all %>% 
  filter(row.names(.)=="GRIN1")

cluster_2_v_all %>% 
  rownames_to_column("gene_name") %>% 
  filter(gene_name=="GRIN1")

#### Using metadata for other comparisons
# making new metadata column to see if any sex differences for the same cell type
Seu_hum@meta.data <- Seu_hum@meta.data %>% 
  mutate(sex_subclass = paste(donor_sex_label, subclass_label, sep="_"))

table(Seu_hum$sex_subclass) #see how many cells we have for each sex+subclass combo

Idents(Seu_hum) <- "sex_subclass" #setting this new column as our active identity

unique(Idents(Seu_hum)) #seeing what our options are for making comparisons

#finding genes that are DE in female-derived microglia vs male-derived
F_microglia_vs_M_microglia <- FindMarkers(Seu_hum, ident.1 = "F_Microglia", ident.2 = "M_Microglia", 
                                          logfc.threshold = log(2), min.pct = 0.25)

F_microglia_vs_M_microglia %>% 
  arrange(p_val_adj)

Seu_hum@assays$RNA@counts["GFAP",] #checking if a gene is present in our counts matrix

### Visualizations
#getting top 6 marker genes for distinguishing cluster 8 cells, saving to plot below
features <- cluster_8_v_all %>% 
  arrange(p_val_adj) %>% 
  head(n=6) %>% 
  row.names()

Idents(Seu_hum) <- "seurat_clusters" #setting our active identity back to our clusters

#### Violin Plot
VlnPlot(Seu_hum, features = features)

#### Feature Plot
FeaturePlot(Seu_hum, features = features)

#### Dot Plot
DotPlot(Seu_hum, features = features) + RotatedAxis()

#### Heat Map
DoHeatmap(subset(Seu_hum, downsample = 100), features = features, size = 3, slot="data") +
  scale_fill_viridis_c()