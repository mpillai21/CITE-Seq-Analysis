#load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)
library(dsb)
library(magrittr)
library(ComplexHeatmap)

#Set working directory
setwd("/Users/Meenakshi Pillai/Desktop/CITE_Seq/Results")

# read raw data using the Seurat function "Read10X" 
raw1 <- Seurat::Read10X("/Users/Meenakshi Pillai/Desktop/Macaque Reference/GSM6674201_BM22-LN1_raw")
cells1 <- Seurat::Read10X("/Users/Meenakshi Pillai/Desktop/Macaque Reference/GSM6674201_BM22-LN1_filtered")

# define cell-containing barcodes and separate cells and empty drops
stained_cells1 = colnames(cells1$`Gene Expression`)
background1 = setdiff(colnames(raw1$`Gene Expression`), stained_cells1)

# split the data into separate matrices for RNA and ADT
prot1 = raw1$`Antibody Capture`
rna1 = raw1$`Gene Expression`

# create metadata of droplet QC stats used in standard scRNAseq processing
mtgene1 = grep(pattern = "^MT-", rownames(rna1), value = TRUE) # used below

md1 = data.frame(
  rna.size = log10(Matrix::colSums(rna1)), 
  prot.size = log10(Matrix::colSums(prot1)), 
  n.gene = Matrix::colSums(rna1 > 0), 
  mt.prop = Matrix::colSums(rna1[mtgene1, ]) / Matrix::colSums(rna1)
)
# add indicator for barcodes Cell Ranger called as cells
md1$drop.class = ifelse(rownames(md1) %in% stained_cells1, 'cell', 'background')

# remove barcodes with no evidence of capture in the experiment
md1 = md1[md1$rna.size > 0 & md1$prot.size > 0, ]

background_drops1 = rownames(
  md1[ md1$prot.size > 1.5 & 
        md1$prot.size < 3 & 
        md1$rna.size < 2.5, ]
) 
background.adt.mtx1 = as.matrix(prot1[ , background_drops1])

# calculate statistical thresholds for droplet filtering.
cellmd1 = md1[md1$drop.class == 'cell', ]

# filter drops with + / - 3 median absolute deviations from the median library size
rna.mult1 = (3*mad(cellmd1$rna.size))
prot.mult1 = (3*mad(cellmd1$prot.size))
rna.lower1 = median(cellmd1$rna.size) - rna.mult1
rna.upper1 = median(cellmd1$rna.size) + rna.mult1
prot.lower1 = median(cellmd1$prot.size) - prot.mult1
prot.upper1 = median(cellmd1$prot.size) + prot.mult1

# filter rows based on droplet quality control metrics
qc_cells1 = rownames(
  cellmd1[cellmd1$prot.size > prot.lower1 & 
           cellmd1$prot.size < prot.upper1 & 
           cellmd1$rna.size > rna.lower1 & 
           cellmd1$rna.size < rna.upper1 & 
           cellmd1$mt.prop < 0.14, ]
)
cell.adt.raw1 = as.matrix(prot1[ , qc_cells1])
cell.rna.raw1 = rna1[ ,qc_cells1]
cellmd1 = cellmd1[qc_cells1, ]

# normalize and denoise with dsb with 
cells.dsb.norm1 = DSBNormalizeProtein(
  cell_protein_matrix = cell.adt.raw1, 
  empty_drop_matrix = background.adt.mtx1, 
  denoise.counts = TRUE, 
  use.isotype.control = FALSE
)

# integrating with Seurat
stopifnot(isTRUE(all.equal(rownames(cellmd1), colnames(cell.adt.raw1))))
stopifnot(isTRUE(all.equal(rownames(cellmd1), colnames(cell.rna.raw1))))

# create Seurat object note: min.cells is a gene filter, not a cell filter
s = Seurat::CreateSeuratObject(counts = cell.rna.raw1, 
                               meta.data = cellmd1,
                               assay = "RNA", 
                               min.cells = 3)

# add dsb normalized matrix "cell.adt.dsb" to the "CITE" data (not counts!) slot
s[["CITE"]] = Seurat::CreateAssayObject(data = cells.dsb.norm1)


# define proteins to use in clustering (non-isptype controls)
prots1 = rownames(s@assays$CITE@data)[1:100]

# cluster and run umap 
s = Seurat::FindNeighbors(object = s, dims = NULL,assay = 'CITE', 
                          features = prots1, k.param = 30, 
                          verbose = FALSE)

# direct graph clustering 
s = Seurat::FindClusters(object = s, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8), 
                         algorithm = 3, 
                         graph.name = 'CITE_snn', 
                         verbose = FALSE)
# umap (optional)
s = RunUMAP(object = s, assay = "CITE", features = prots1,
                    seed.use = 1990, min.dist = 0.6, n.neighbors = 30,
                    verbose = FALSE)
Idents(s) = "CITE_snn_res.1.8"

DefaultAssay(s) <- "CITE"
DimPlot(s, reduction = "umap")

VariableFeatures(s) = prots1
#PCA technique #1
s = s %>% 
  ScaleData() %>% 
  RunPCA(reduction.name = 'apca', verbose = FALSE)
DimPlot(s, reduction = "apca")
ElbowPlot(s, ndims = 20, reduction = "apca")

#PCA #2
s <- RunPCA(s, features = rownames(s), reduction.name = "pca_adt", reduction.key = "pca_adt_", verbose = FALSE)
DimPlot(s, reduction = "pca_adt")
ElbowPlot(s, ndims = 20, reduction = "pca_adt")

# make results dataframe 
d = cbind(s@meta.data, 
          as.data.frame(t(s@assays$CITE@data),
          s@reductions$umap@cell.embeddings)
)


# use pearson residuals as normalized values for pca 
DefaultAssay(s) = "RNA"
s = NormalizeData(s, verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = 'vst', verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
ElbowPlot(s, ndims = 30, reduction = "pca")

#s = Seurat::RunUMAP(object = s, assay = "RNA", features = prots1,
                    #seed.use = 1990, min.dist = 0.2, n.neighbors = 30,
                    #verbose = FALSE)



# run WNN 
s = FindMultiModalNeighbors(
  s, reduction.list = list("pca", "apca"), 
  dims.list = list(1:10, 1:10), 
  modality.weight.name = "RNA.weight", 
  verbose = FALSE
)

# cluster 
s <- FindClusters(s, graph.name = "wsnn", 
                  algorithm = 3, 
                  resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8), 
                  verbose = FALSE, 
                  random.seed = 1990)

# create multimodal heatmap 
vf = VariableFeatures(s,assay = "RNA")

# find marker genes for the joint clusters 
Idents(s) = "wsnn_res.1.8"
DefaultAssay(s)  = "RNA"
rnade = FindAllMarkers(s, features = vf, only.pos = TRUE, verbose = FALSE)
gene_plot = rnade %>% 
  dplyr::filter(avg_log2FC > 1 ) %>%  
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(5) %$% gene %>% unique 
markers_rna <- rnade %>%
  group_by(cluster) %>%
  slice_max(n = 6, order_by = avg_log2FC)

cite_data = GetAssayData(s,slot = 'data',assay = 'CITE') %>% t()
rna_subset = GetAssayData(s,assay = 'RNA',slot = 'data')[gene_plot, ] %>%
  as.data.frame() %>% 
  t() %>% 
  as.matrix()

# combine into dataframe 
d = cbind(s@meta.data, cite_data, rna_subset) 

# calculate the median protein expression per cluster
dat_plot = d %>% 
  dplyr::group_by(wsnn_res.0.2) %>% 
  dplyr::summarize_at(.vars = c(prots1, gene_plot), .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("wsnn_res.0.2") 

prot_col = circlize::colorRamp2(breaks = seq(-1,25, by = 1), 
                                colors = viridis::viridis(n = 27, option = "B"))
p1 = Heatmap(t(dat_plot)[prots1, ], 
             name = "protein", 
             col = prot_col, 
             use_raster = T,
             row_names_gp = gpar(color = "black", fontsize = 5)
)


# mRNA heatmap 
mrna = t(dat_plot)[gene_plot, ]
rna_col = circlize::colorRamp2(breaks = c(-2,-1,0,1,2), 
                               colors = colorspace::diverge_hsv(n = 5))
p2 = Heatmap(t(scale(t(mrna))), 
             name = "mRNA", 
             col = rna_col,
             use_raster = T, 
             clustering_method_columns = 'average',
             column_names_gp = gpar(color = "black", fontsize = 7), 
             row_names_gp = gpar(color = "black", fontsize = 5))


# combine heatmaps 
d=dist(p1@matrix)
d[is.na(d)]=10^50 #this has to be larger than max(dist)
p1.1 <- pheatmap(p1@matrix[1:20,], clustering_distance_rows = d, clustering_method = "single", name = "protein")
d2 <- dist(p2@matrix)
d2[is.na(d2)]=10^50
p2.1 <- pheatmap(p2@matrix, clustering_distance_rows = d2, clustering_method = "single", name= "mRNA")
ht_list = p1 %v% p2
draw(ht_list)

DefaultAssay(s) <- "CITE"
#Visualize the proteins on the clusters
FeaturePlot(s, reduction = "umap", features = rownames(p1@matrix)[13:24], min.cutoff = "q05", max.cutoff = "q95", ncol = 4)


CITE_Markers <- FindAllMarkers(s, assay = "CITE",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CITE_Markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
RNA_Markers <- FindAllMarkers(s, assay = "RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
RNA_Markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)

#Annotation
new.cluster.ids <- c("B_Naive", "B_gc_light_zone", "T_CD8", "Unknown", "Unknown", "T_CD4_Tfh", "T_CD4_Tfh", "pDC", "B_gc_light_zone", "NK", "T_CD4_Naive", "T_CD8", "T_CD8", "T_CD4_Tfh", "Unknown", "B_marginal_zone", "B_Naive", "B_gc_light_zone", "T_CD4_Naive", "pDC", "Macrophage", "B_Memory", "pDC", "Macrophage", "pDC"  )
names(new.cluster.ids) <- levels(s)
s <- RenameIdents(s, new.cluster.ids)
s <- AddMetaData(s, metadata = as.data.frame(Idents(s)))


#Save the final object
saveRDS(s, "Final_CITE_Seq.rds")

DimPlot(s, reduction = "umap")
