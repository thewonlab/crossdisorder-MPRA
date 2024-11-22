################################################
#############SEURAT ANALYSIS####################
################################################
#New version (v5) of Seurat
#install.packages('Seurat')

#Supplementary packages to help Seurat run faster, part of v5
#setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
#install.packages(c("BPCells", "presto", "glmGamPoi"))

options(stringsAsFactors = F)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(BPCells)
#library(presto)
library(glmGamPoi)
library(progress)

rm(list = ls())


############################# Input directory paths and sample info before running the whole script

# Directory to save data
savedir <- "/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Seurat_Outputs"

# Directory of the cellranger outs folder
# Use filtered_feature_bc_matrix folder. Here, barcodes from empty GEM (background noise) are excluded.
# For more details about files in the filtered_feature_bc_matrix folder, see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
#For UMAP
cellrangeroutsdir <-  paste0("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/CellRanger_Outputs/filtered_feature_bc_matrix_In", 1:5)


############################# Setup the seurat object for gene x cell matrix and guide x cell matrix

# Set directory to save data
setwd(savedir)

# Load data from cellranger outs folder
# If you ran cellranger with a feature reference csv file which indicates the data has multiple data types (e.g., Gene Expression + CRISPR Guide Capture), a list containing a sparse matrix of the data from each type will be returned. Otherwise a sparse matrix containing the gene expression data will be returned.
# Matrix means UMI count matrix (row=features, cols=cells)
dat <- Read10X(data.dir = cellrangeroutsdir)

# You can check the data types included (e.g., Gene Expression + CRISPR Guide Capture) and their matrix dimensions using this command
str(dat)

# To access each data type matrix, use $
dim(dat$`Gene Expression`)
dim(dat$`CRISPR Guide Capture`)

# Note that the column dimensions are same for Gene Expression and CRISPR Guide Capture matrices
# Double check that the column names (cells) are identical between Gene Expression and CRISPR Guide Capture matrices
identical(colnames(dat$`CRISPR Guide Capture`), colnames(dat$`Gene Expression`))

# Create a seurat objects that contains both gene matrix & guide matrix (like a multimodal object)
# https://satijalab.org/seurat/articles/multimodal_vignette.html
# Default 'assay' argument for gene matrix = "RNA"
gex <- CreateSeuratObject(counts = dat$`Gene Expression`, project = "cropseq")
gex # Check dimension

# Add guide matrix to assay named 'crispr'
gex[["crispr"]] <- CreateAssayObject(counts = dat$`CRISPR Guide Capture`)
# Check that the "crispr" assay is added
gex

# Use str() command to check the structure of the seurat object and identifiers to access each component
str(gex)
# For example, to access the gene count matrix from the seurat object
gex@assays$RNA$counts
# To access gene names
gex@assays$RNA$counts@Dimnames[[1]]
# To access cell IDs
gex@assays$RNA$counts@Dimnames[[2]]
# To access guide count matrix
gex@assays$crispr$counts
# To access guide names
gex@assays$crispr$counts@Dimnames[[1]]


# Add mitochondrial gene percentage for each cell. Note: The '[[' operator can add columns to object metadata
# NOTE: for mouse genome, use "^mt-" instead of "^MT"
gex[["percent.mt"]] <- PercentageFeatureSet(gex, pattern = "^MT-")
# Check whether mt percentage values look normal (i.e., not all 0)
gex[["percent.mt"]]

# Extract batch numbers from rownames
batch_numbers <- as.integer(sub("^([0-9]+)_.*", "\\1", rownames(gex@meta.data)))

# Add batch_num column to metadata
gex@meta.data$batch_num <- batch_numbers

# Prune cells based on QC metrics
# Refer to selection criteria in other papers, but mainly decide cutoffs values based on how your actual data distribution looks like
# Iterate analysis with different parameters (start from less stringent cutoffs)
# Example: gex_qc <- subset(gex, subset = percent.mt < 10 & nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 40000)
gex_qc <- subset(gex, subset = percent.mt < 15) # Note that guide matrix will be simultaneously pruned
gex_qc # Check dimension
gex_qc@assays$crispr # Check guide matrix dimension

# Filter out cells with less than 1250 unique genes
gex_filtered <- subset(gex_qc, subset = nFeature_RNA >= 1250)

# Check the dimension of the filtered object
gex_filtered

#Save Seurat object with all metadata
saveRDS(gex_filtered, file = "/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Seurat_Outputs/gex_meta_mitofiltered_genefiltered.rds")
#######################################################################
################### Quality control in gene matrix#####################

# View QC metrics
pdf("01_before_QC.pdf")
VlnPlot(gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(gex, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(gex, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# # Prune cells based on QC metrics
# # Refer to selection criteria in other papers, but mainly decide cutoffs values based on how your actual data distribution looks like
# # Iterate analysis with different parameters (start from less stringent cutoffs)
# # Example: gex_qc <- subset(gex, subset = percent.mt < 10 & nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 40000)
# gex_qc <- subset(gex, subset = percent.mt < 15) # Note that guide matrix will be simultaneously pruned
# gex_qc # Check dimension
# gex_qc@assays$crispr # Check guide matrix dimension

# View plots after QC
pdf("02_after_QC.pdf")
VlnPlot(gex_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(gex_qc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(gex_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


############################# Quality control in guide matrix

# Take the guide count matrix
guidemat <- data.frame(gex_qc@assays$crispr@counts)

# Total count of all guides per cell
table1 <- data.frame(colSums(guidemat))
colnames(table1) <- "Total_count_of_all_guides_per_cell_before_removing_0_guide_cells"

# Number of different guides per cell
table1$Number_of_different_guides_per_cell_before_removing_0_guide_cells <- colSums(guidemat != 0)

# Check the number and percentage of cells with 0 guide expression
a <- nrow(table1[table1$Total_count_of_all_guides_per_cell_before_removing_0_guide_cells!=0,])
b <- nrow(table1[table1$Total_count_of_all_guides_per_cell_before_removing_0_guide_cells==0,])
percent_cells_zero_guide <- b/(a+b)*100

# Remove cells with 0 guide expression
guidemat2 <- guidemat[, colSums(guidemat)!=0]

# Check the number and percentage of cells expressing each guide
table2 <- data.frame(rowSums(guidemat2 != 0))
colnames(table2) <- "Number_of_cells_per_guide_after_removing_0_guide_cells"
table2$Total_cell_number_after_removing_0_guide_cells <- ncol(guidemat2)
table2$Percent_of_cells_per_guide_after_removing_0_guide_cells <- table2$Number_of_cells_per_guide_after_removing_0_guide_cells/ncol(guidemat2)*100

# Total count of each guide across cells
table2$Total_count_of_each_guide_after_removing_0_guide_cells <- rowSums(guidemat2)

# Plot distributions
pdf("03_guide_distributions.pdf")

ggplot(table1, aes(x=Total_count_of_all_guides_per_cell_before_removing_0_guide_cells)) + geom_histogram(binwidth=1, color="black", fill='grey')

ggplot(table1, aes(x=Number_of_different_guides_per_cell_before_removing_0_guide_cells)) + geom_histogram(binwidth=1, color="black", fill='grey')

ggplot(table2, aes(x=Number_of_cells_per_guide_after_removing_0_guide_cells)) + geom_histogram(binwidth=50, color="black", fill='grey')

ggplot(table2, aes(x=Percent_of_cells_per_guide_after_removing_0_guide_cells)) + geom_histogram(binwidth=1, color="black", fill='grey')

dev.off()

# Save tables
write.csv(table1, file = "03_guide_table1.csv")
write.csv(table2, file = "03_guide_table2.csv")

# Prune the gex_qc object
# Note that in guidemat2, dash is converted to dot in cell IDs. Convert back to dash so it matches with the cell IDs in gex_qc object.
colnames(guidemat2) <- gsub('[.]', '-', colnames(guidemat2)) # You must bracket dot because dot means 'any character' for gsub
colnames(guidemat2) <- sub("^X", "", colnames(guidemat2)) #There is also an X in front that needs to be removed
gex_qc2 <- gex_qc[, colnames(guidemat2)] # Note that guide matrix will be simultaneously pruned
gex_qc2 # Check dimension

# View QC values again
pdf("04_after_QC2.pdf")
VlnPlot(gex_qc2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(gex_qc2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(gex_qc2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


##################  Normalization & save

# SCT normalization
gex_qc2_sct <- SCTransform(gex_qc2, vars.to.regress = "percent.mt")

# Save pruned and normalized object for next step
#saveRDS(gex_qc2_sct, file = "04_gex_qc2_sct.rds")
