
######################################################
######################U-Map Analysis##################
######################################################
library(ghibli)
ghibli_palettes$PonyoLight

library(dplyr)
library(patchwork)
install.packages("Seurat")
library(Seurat)
library(sceptre)
rm(list = ls())
#Set working directory
setwd("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/Data_Outputs/Seurat_Outputs/SeuratObjects/UMAP")

#load your own data, it is already SCT transformed, and NormalizeData(), ScaleData(), and FindVariableFeatures() has been run!
CD_CROPseq <- readRDS("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/Data_Outputs/Seurat_Outputs/SeuratObjects/A02_gex_sct_integrated.rds")

# # Check the unique batches, ideally only analyze one batch at a time for UMAPs, but only if not previously integrated
# selected_batch <- unique(CD_CROPseq@meta.data$batch_num)[2]  # Change this to the desired batch
# 
# # Subset the Seurat object based on the selected batch
# CD_CROPseq_subset <- subset(CD_CROPseq, subset = batch_num == selected_batch)

# Run UMAP on the integrated data
CD_CROPseq_integrated <- RunUMAP(CD_CROPseq, dims = 1:20)
DimPlot(CD_CROPseq_integrated, reduction = "umap")

# # Run UMAP on the subset
# CD_CROPseq_subset <- RunUMAP(CD_CROPseq_subset, dims = 1:20)
# DimPlot(CD_CROPseq_subset, reduction = "umap")

#Cluster based on similarities between individual cell expression for integrated Seurat object
CD_CROPseq_integrated_cluster <- FindNeighbors(CD_CROPseq_integrated, dims = 1:10)
CD_CROPseq_integrated_cluster <- FindClusters(CD_CROPseq_integrated_cluster, resolution = 0.5)

#Cluster based on similarities between individual cell expression for subset (1 batch) Seurat object
# CD_CROPseq_subset_cluster <- FindNeighbors(CD_CROPseq_subset, dims = 1:10)
# CD_CROPseq_subset_cluster <- FindClusters(CD_CROPseq_subset_cluster, resolution = 0.5)

#Visualize the clusters for the integrated Seurat object
CD_CROPseq_integrated_cluster <- RunUMAP(CD_CROPseq_integrated_cluster, dims = 1:10)
DimPlot(CD_CROPseq_integrated_cluster, reduction = "umap")

#Visualize the clusters (this is the U-MAP) for separate batches
# CD_CROPseq_subset_cluster <- RunUMAP(CD_CROPseq_subset_cluster, dims = 1:10)
# DimPlot(CD_CROPseq_subset_cluster, reduction = "umap")

#Save clustered data for integrated Seurat Object
saveRDS(CD_CROPseq_integrated_cluster, file = "CD_CROPseq_cluster.rds")

#Save clustered data for subset
#saveRDS(CD_CROPseq_subset_cluster, file = paste0("CD_CROPseq_cluster_", selected_batch, ".rds"))

#########################Load from here#########################################
rm(list = ls())

#Load the Cluster file
# CD_CROPseq_integrated <- readRDS("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/Data_Outputs/Seurat_Outputs/SeuratObjects/UMAP/CD_CROPseq_cluster.rds")
 CD_CROPseq_batch1 <- readRDS("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/Data_Outputs/Seurat_Outputs/SeuratObjects/UMAP/CD_CROPseq_cluster_batch1.rds")
# CD_CROPseq_batch2 <- readRDS("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/Data_Outputs/Seurat_Outputs/SeuratObjects/UMAP/CD_CROPseq_cluster_batch2.rds")

#Now let's aim for the markers we care about (Serotonergic & Dopaminergic)
# Find markers for every cluster (calculates differential marker gene expression per cluster)
# Find markers for every cluster compared to all remaining cells, report only the positive ones
CD_CROPseq_cluster <- CD_CROPseq_batch1 #so we don't have to reload the clusters everytime it fails

#Check dopamine marker and opioid marker
interest_markers <- c("FOS", "HTR2A")
FeaturePlot(CD_CROPseq_cluster, 
            features = interest_markers, 
            cols = c("lightgray", "blue"), alpha = 0.3)

#Check STAG1, SETD1A, GRIN2A (MPRA-CROP-seq markers)
MPRACROP_markers <- c("STAG1", "SETD1A", "GRIN2A", "PRPH") 
FeaturePlot(CD_CROPseq_cluster, 
            features = MPRACROP_markers, 
            cols = c("gray", "blue"), alpha = 0.3)

#Check proliferative markers for NPCs and iPSCs
iPSC_markers <- c("POU5F1", "SOX2", "NANOG") #OCT4 -> POU5F1
#OCT4 (Octamer-binding transcription factor 4) is a key transcription factor involved in maintaining pluripotency.
#SOX2 works in conjunction with OCT4 to regulate the expression of genes involved in maintaining stem cell identity.
#NANOG is often used as a marker for undifferentiated pluripotent stem cells.
FeaturePlot(CD_CROPseq_cluster, 
            features = iPSC_markers, 
            cols = c("gray", "red"), alpha = 0.3)

NPC_markers <- c("SOX1", "PAX6", "NES")
#SOX1 is often used as a marker for early neural commitment.
#PAX6 is expressed in neural progenitor cells and is important for the regulation of neural cell fate.
#NES (Nestin) is an intermediate filament protein expressed in neural stem cells and neural progenitor cells.
FeaturePlot(CD_CROPseq_cluster, 
            features = NPC_markers, 
            cols = c("gray", "lightblue"), alpha = 0.3)

EarlyNeuron_markers <- c("DCX", "TUBB3", "NEUROD1") #TUJ1 -> TUBB3
#DCX is often considered a marker for early stages of neuronal development, especially in migrating neuroblasts.
#TUJ1 is a marker for immature neurons and is commonly used to identify early neuronal differentiation.
#NEUROD1 is a transcription factor expressed in neuronal progenitor cells and is involved in the determination and differentiation of neurons. 
FeaturePlot(CD_CROPseq_cluster, 
            features = EarlyNeuron_markers, 
            cols = c("gray", "orange"), alpha = 0.3)

LateNeuron_markers <- c("MAP2", "RBFOX3", "SYP")
#MAP2 is a cytoskeletal protein expressed in the dendrites of mature neurons.
#NEUN, also known as RBFOX3, is a nuclear protein expressed in postmitotic neurons. 
#SYP is a synaptic vesicle protein that plays a role in regulating neurotransmitter release, identifies functional synapses.
FeaturePlot(CD_CROPseq_cluster, 
            features = LateNeuron_markers, 
            cols = c("gray", "lightgreen"), alpha = 0.3)

GlutNeuron_markers <- c("GRIN2A", "GRIN2B", "SLC17A6", "GRIK1")
#GRIN2A and GRIN2B are both glutamate receptors, with GRIN2B (in humans) typically being silenced after birth.
#SLC17A6 is otherwise known as VGLUT2, involved in the transport of glutamate.
#GRIK1 encodes a subunit of kainate receptors, another type of glutamate receptor.
FeaturePlot(CD_CROPseq_cluster, 
            features = GlutNeuron_markers, 
            cols = c("gray", "blue"), alpha = 0.3)

DimPlot(CD_CROPseq_cluster, reduction = "umap", group.by = "batch_num", alpha = 0.3)

# # Check markers for dopaminergic and serotonergic marker receptors
# dopa_receptors <- c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5")
# sero_receptors <- c("HTR1A", "HTR1B", "HTR2A", "HTR2B", "HTR2C")
# marker_genes <- c(dopa_receptors, sero_receptors)
# #marker_genes <- c("DRD1","DRD2","DRD3","DRD4","DRD5","HTR2A","HTR2C","HTR3A")
# 
# # Check markers for dopaminergic and serotonergic marker receptors
# dopa_markers <- c("TH", "DDC", "SLC6A3")
# sero_markers <- c("TPH2", "SLC6A4")
# glut_markers <- c("SLC17A7", "SLC17A6", "CAMK2A")
# gaba_markers <- c("GAD1", "GAD2", "SLC32A1")
# marker_genes <- c(dopa_markers, sero_markers, glut_markers, gaba_markers)
# 
# ampa_receptors <- c("GRIA1", "GRIA2", "GRIA3", "GRIA4")
# kainate_receptors <- c("GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5")
# nmda_receptors <- c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D")
# glut_marker_receptors <- c(ampa_receptors, kainate_receptors, nmda_receptors)
# 
# # Visualize the UMAP plot colored by marker groups
# pdf("UMAP_result.pdf", width = 15, height = 8)
# FeaturePlot(CD_CROPseq_cluster, features = marker_genes)
# dev.off()
# 
# #Plot UMAP for batches in one plot, orig.ident is your batch
# DimPlot(CD_CROPseq_cluster, reduction = "umap", group.by = c("orig.ident"), alpha = 0.3)
# 
# #Plot UMAP for batches separately
# DimPlot(gex_sct, reduction = "umap.dr", split.by = c("batch_num"), alpha = 1)

# pdf("UMAP_glut_result.pdf", width = 20, height = 15)
# FeaturePlot(CD_CROPseq_cluster, 
#             features = ampa_receptors, 
#             cols = c("gray", "red"))
# FeaturePlot(CD_CROPseq_cluster, 
#             features = kainate_receptors, 
#             cols = c("gray", "blue"))
# FeaturePlot(CD_CROPseq_cluster, 
#             features = nmda_receptors, 
#             cols = c("gray", "green"))
# FeaturePlot(CD_CROPseq_cluster, 
#             features = glut_markers, 
#             cols = c("gray", "pink"))
# dev.off()
# # 