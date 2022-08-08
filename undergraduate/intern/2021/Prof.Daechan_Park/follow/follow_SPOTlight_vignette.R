install.packages("devtools") #Install it if you need
devtools::install_github("https://github.com/MarcElosua/SPOTlight/tree/spotlight-0.1.7") #Download the sample_data
                                                                                         #sample_data will saved in "D:\R\R-4.2.1\library\SPOTlight" by default setting 
#Bring the packages
library(Matrix)
library(data.table)
library(Seurat) #Single cell data package
library(SeuratData)
library(dplyr) #data structure package
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(NMF)
library(magrittr)

path_to_data <- system.file(package = "SPOTlight") #Set the path to "SPOTlight"
cortex_sc <- readRDS(glue::glue("{path_to_data}/allen_cortex_dwn.rds")) #Load the single cell sequencing data

# Download the data from SeuratData
if (!("stxBrain.SeuratData" %in% rownames(InstalledData()))) {
  InstallData("stxBrain")
}

# Load SeuratData
anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")

#Create seurat object 
set.seed(123)
cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
  Seurat::RunPCA(., verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

#cell type plot
Seurat::DimPlot(cortex_sc,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()

#cell type dataset
cortex_sc@meta.data %>%
  dplyr::count(subclass) %>%
  gt::gt(.[-1, ]) %>%
  gt::tab_header(
    title = "Cell types present in the reference dataset",
  ) %>%
  gt::cols_label(
    subclass = gt::html("Cell Type")
  )

#Compute marker genes
Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$subclass
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc,  
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)

saveRDS(object = cluster_markers_all,
        file = here::here("inst/markers_sc.RDS"))
#set directory_220712
cluster_markers_all<-readRDS('markers_sc.RDS')

#SPOTlight Decomposition
set.seed(123)

spotlight_ls <- spotlight_deconvolution(
  se_sc = cortex_sc,
  counts_spatial = anterior@assays$Spatial@counts,
  clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 3000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

saveRDS(object = spotlight_ls, file = here::here("inst/spotlight_ls.rds"))

#set directory_220712
spotlight_ls <- readRDS(file = "spotlight_ls.rds")

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

#Assess deconvolution
h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod[[2]])

#topic profile cell type plot
topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90), 
  axis.text = ggplot2::element_text(size = 12))

#topic proportion within cell type
topic_profile_plts[[1]] + theme(axis.text.x = element_text(angle = 90), 
                                axis.text = element_text(size = 12))

#take a look at which genes are the most important for each topic
basis_spotlight <- data.frame(NMF::basis(nmf_mod[[1]]))

colnames(basis_spotlight) <- unique(stringr::str_wrap(nmf_mod[[2]], width = 30))

basis_spotlight %>%
  dplyr::arrange(desc(Astro)) %>%
  round(., 5) %>% 
  DT::datatable(., filter = "top")

#Visualization
# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(anterior)

decon_df <- decon_mtrx %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

anterior@meta.data <- anterior@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

#Specific cell-types
Seurat::SpatialFeaturePlot(
  object = anterior,
  features = c("L2.3.IT", "L6b", "Meis2", "Oligo"),
  alpha = c(0.1, 1))

#Spatial scatterpies
cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
#integration of cell type image
SPOTlight::spatial_scatterpie(se_obj = anterior,
                              cell_types_all = cell_types_all,
                              img_path = "tissue_lowres_image.png", #this is real file name, / means that the path of file lol
                              pie_scale = 0.4)

#Plot spot composition of spots containing cell-types of interest
SPOTlight::spatial_scatterpie(se_obj = anterior,
                              cell_types_all = cell_types_all,
                              img_path = "tissue_lowres_image.png",
                              cell_types_interest = "L6b",
                              pie_scale = 0.8)

#composition of spots containing another cell-type
SPOTlight::spatial_scatterpie(se_obj = anterior,
                              cell_types_all = cell_types_all,
                              img_path = "tissue_lowres_image.png",
                              cell_types_interest = "VLMC",
                              pie_scale = 0.8)
#Spatial interaction graph
#bring the function for interaction graph
graph_ntw <- SPOTlight::get_spatial_interaction_graph(decon_mtrx = decon_mtrx[, cell_types_all])

#graph drawing protocol
deg <- degree(graph_ntw, mode="all")
# Get color palette for difusion
edge_importance <- E(graph_ntw)$importance
# Select a continuous palette
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
# Create a color palette
getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
# Get how many values we need
grad_edge <- seq(0, max(edge_importance), 0.1)
# Generate extended gradient palette dataframe
graph_col_df <- data.frame(value = as.character(grad_edge),
                           color = getPalette(length(grad_edge)),
                           stringsAsFactors = FALSE)
# Assign color to each edge
color_edge <- data.frame(value = as.character(round(edge_importance, 1)), stringsAsFactors = FALSE) %>%
  dplyr::left_join(graph_col_df, by = "value") %>%
  dplyr::pull(color)
# Open a pdf file
plot(graph_ntw,
     # Size of the edge
     edge.width = edge_importance,
     edge.color = color_edge,
     # Size of the buble
     vertex.size = deg,
     vertex.color = "#cde394",
     vertex.frame.color = "white",
     vertex.label.color = "black",
     vertex.label.family = "Ubuntu", # Font family of the label (e.g.??Times??, ??Helvetica??)
     layout = layout.circle)

#compute the correlation
#Remove cell types not predicted to be on the tissue
decon_mtrx_sub <- decon_mtrx[, cell_types_all]
decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]
# Compute correlation
decon_cor <- cor(decon_mtrx_sub)
# Compute correlation P-value
p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)
# Visualize
#install.packages('ggcorrplot')
ggcorrplot::ggcorrplot(
  corr = decon_cor,
  p.mat = p.mat[[1]],
  hc.order = TRUE,
  type = "full",
  insig = "blank",
  lab = TRUE,
  outline.col = "lightgrey",
  method = "square",
  # colors = c("#4477AA", "white", "#BB4444"))
  colors = c("#6D9EC1", "white", "#E46726"),
  title = "Predicted cell-cell proportion correlation",
  legend.title = "Correlation\n(Pearson)") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = 90),
    axis.text = ggplot2::element_text(size = 18, vjust = 0.5))

#downsample data
#Downsample scRNAseq to select gene set and number of cells to train the model
se_sc_down <- downsample_se_obj(se_obj = cortex_sc,
                                clust_vr = "subclass",
                                cluster_markers = cluster_markers_all,
                                cl_n = 100,
                                hvg = 3000)
#Train NMF model
start_time <- Sys.time()
nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers_all, 
                        se_sc = se_sc_down, 
                        mtrx_spatial = anterior@assays$Spatial@counts,
                        clust_vr = "subclass",
                        ntop = NULL,
                        hvg = 3000,
                        transf = "uv",
                        method = "nsNMF")
nmf_mod <- nmf_mod_ls[[1]]

#get basis matrix W
w <- basis(nmf_mod)
#w
dim(w)

# get coefficient matrix H
h <- coef(nmf_mod)
#h
dim(h)

#cell-type specific topic profile
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod_ls[[2]]
)

topic_profile_plts[[2]] + theme(axis.text.x = element_text(angle = 90))

#Spot Deconvolution
#Extract count matrix
spot_counts <- anterior@assays$Spatial@counts

# Subset to genes used to train the model
spot_counts <- spot_counts[rownames(spot_counts) %in% rownames(w), ]

ct_topic_profiles <- topic_profile_per_cluster_nmf(h = h,
                                                   train_cell_clust = nmf_mod_ls[[2]])

decon_mtrx <- mixture_deconvolution_nmf(nmf_mod = nmf_mod,
                                        mixture_transcriptome = spot_counts,
                                        transf = "uv",
                                        reference_profiles = ct_topic_profiles, 
                                        min_cont = 0.09)
