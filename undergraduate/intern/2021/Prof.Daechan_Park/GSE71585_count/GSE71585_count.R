##github. (2020). GitHub. Retrieved from https://github.com/MarcElosua/SPOTlight/tree/spotlight-0.1.7

##Sometimes in GSE, the single cell data saved as CSV format
#Download the csv file from GSE
#Bring the packages
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(devtools)
library(Matrix)
library(data.table)
library(Seurat) #sc data package
library(SeuratData)
library(dplyr) #data structure
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(NMF)
library(magrittr)

#Generate single-cell RNA seq data
count <- read.csv("D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/GSE71585_RefSeq_counts.csv", header=T,row.names = 1) #Gene name is in rows, cell number is in column, counts in grid
temp1 <- as.matrix(count) #Convert it as matrix
#Save sparse matrix(sparse matrix is the most of the value is 0)
sparse.gbm <- Matrix(temp1 , sparse = T )
#Create matrix file
writeMM(obj = sparse.gbm, file="matrix_n.mtx")
#Create features and barcode file
write(x = rownames(sparse.gbm), file = "features_n.tsv")
write(x = colnames(sparse.gbm), file = "barcodes_n.tsv")
#Create Seuratboject
getwd()
data_dir <- 'D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv'
test1<-  Read10X(data.dir = data_dir,gene.column = 1, unique.features = TRUE)
test_seurat <- CreateSeuratObject(counts = test1, project = "seurat", min.cells = 3, min.features = 200)

#Customaize the format
test_seurat <- FindVariableFeatures(test_seurat, selection.method = "vst", nfeatures = 2000)

#I have to run RunPCA, RunUMAP to create object for spotlight
set.seed(0812)
test_seurat <- Seurat::SCTransform(test_seurat, verbose = FALSE) %>%
  Seurat::RunPCA(., verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

#I have to run FindNeighbors, FindClusters to create object for spotlight
test_seurat <- FindNeighbors(test_seurat, dims=1:10)
test_seurat <- FindClusters(test_seurat,resolution=0.5)

#Use levels to find out how many clusters in Seurat object
d<-levels(test_seurat)
d
#Use for loop to name each cluster easily
dl<-length(d)-1 #Cluseter starts at 0
nx<-c()
for (i in 0:dl){
  cc<-paste0('Ty',i)
  nx<-c(nx,cc)
} 

#In this case, cell type will be just number, but if you can give cell name

#Name the levels with cell type
new.cluster.ids <- nx
names(new.cluster.ids) <- levels(test_seurat) #Match the levels with cell type
test_seurat <- RenameIdents(test_seurat, new.cluster.ids) #Change the level of  each ident to cell type

#Add the subclass chapter in meta.data
test_seurat[["subclass"]] <- as.character(Idents(object = test_seurat))

saveRDS(object = test_seurat,
        file = here::here("D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/test_seurat_n.rds"))

#Find markers
cluster_markers_all <- Seurat::FindAllMarkers(object = test_seurat, 
                                                assay = "SCT",
                                                slot = "data",
                                                verbose = TRUE, 
                                                only.pos = TRUE)

saveRDS(object = cluster_markers_all,
        file = here::here("D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/markers_sc_n.rds"))

#====================================================

#Now apply the SeuratObject in the SPotlight vignette

#install.packages("devtools") #Install it if you need
#devtools::install_github("https://github.com/MarcElosua/SPOTlight/tree/spotlight-0.1.7") #Download the sample_data
#sample_data will saved in "D:\R\R-4.2.1\library\SPOTlight" by default setting 

#Load data
path_to_data <- system.file(package = "SPOTlight")
test_seurat <- readRDS(glue::glue("D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/test_seurat_n.rds"))

if (! "stxBrain" %in% SeuratData::AvailableData()[, "Dataset"]) {
  # If dataset not downloaded proceed to download it
  SeuratData::InstallData("stxBrain")
}

# Load data
anterior_a <- SeuratData::LoadData("stxBrain", type = "anterior1")

#check the cluster
Seurat::DimPlot(test_seurat,
                group.by = "subclass",
                label = TRUE) + Seurat::NoLegend()

#check the cell type
test_seurat@meta.data %>%
  dplyr::count(subclass) %>%
  gt::gt(.[-1, ]) %>%
  gt::tab_header(
    title = "Cell types present in the reference dataset",
  ) %>%
  gt::cols_label(
    subclass = gt::html("Cell Type")
  )

#I already found the markers at line 77
#Load markers
cluster_markers_all<- readRDS(file = "D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/markers_sc_n.rds")

#SPOTlight Decomposition
set.seed(0812)
spotlight_ls <- spotlight_deconvolution(
  se_sc = test_seurat,
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

saveRDS(object = spotlight_ls, file = here::here("D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/spotlight_ls_n.rds"))

#Read RDS object
spotlight_ls <- readRDS(file = here::here("D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/spotlight_ls_n.rds"))

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

#Assess deconvolution
h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts_a <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod[[2]])

topic_profile_plts_a[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90), 
  axis.text = ggplot2::element_text(size = 12))

#Look at the individual topic profiles
topic_profile_plts_a[[1]] + theme(axis.text.x = element_text(angle = 90), 
                                  axis.text = element_text(size = 12))


#Look at which genes are the most important
basis_spotlight_a <- data.frame(NMF::basis(nmf_mod[[1]]))

colnames(basis_spotlight_a) <- unique(stringr::str_wrap(nmf_mod[[2]], width = 30))

basis_spotlight_a %>%
  dplyr::arrange(desc(Ty0)) %>%
  round(., 5) %>% 
  DT::datatable(., filter = "top")

#Visualization
# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"] #Exclude res_ss column
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0 #Convert the value that smaller that 0.08 to o
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"]) #Add res_ss column again
rownames(decon_mtrx) <- colnames(anterior_a) #Give the spot information

#I guess, with single cell data, we can classify the Ty9, but there is no such cell type in tissue which is used in reference
#I used different single cell data in this project, so probably it can happen

decon_df_a <- decon_mtrx %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

anterior_a@meta.data <- anterior_a@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df_a, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

#See the Specific cell-types plot
#View the whole cell type
nx
#Select cell types
nxs<-nx[1:4] #In this case, Ty1 to Ty4

#Draw the plot
Seurat::SpatialFeaturePlot(
  object = anterior_a,
  features = c(nxs),
  ncol = 2, #Set the display
  alpha = c(0.1, 1)) #Set the opacity

#Spatial_scatterpie plot will show spots containing at least one of cell types we select
#In this case, Spatial scatterpies for all cell types
cell_types_all_a <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
SPOTlight::spatial_scatterpie(se_obj = anterior_a,
                              cell_types_all = cell_types_all_a, 
                              img_path = "D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/tissue_lowres_image.png",
                              pie_scale = 0.4) 

#Select cell types
nxs1<-nx[1:2] #In this case, Ty1 to Ty4

#Spatial scatterpies for specific cell type
SPOTlight::spatial_scatterpie(se_obj = anterior_a,
                              cell_types_all = cell_types_all_a,
                              img_path = "D:/pankyung/intern_back_up/sub/GSE71585_RefSeq_counts (1).csv/tissue_lowres_image.png",
                              cell_types_interest = nxs1, #Show spots which contain at least Ty0 or Ty1
                              pie_scale = 0.8)