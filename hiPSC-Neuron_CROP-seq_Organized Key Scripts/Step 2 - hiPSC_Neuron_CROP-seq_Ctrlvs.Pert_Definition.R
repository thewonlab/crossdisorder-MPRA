options(stringsAsFactors = F)
library(dplyr)
library(MASS)
library(reshape2)
library(Seurat)

rm(list = ls())

# Directory to save results
#savedir <- "/proj/hyejunglab/cropseq/Jiseok/analysis/"
#setwd(savedir)

# Load in data
obj <- readRDS("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Seurat_Outputs/gex_meta_mitofiltered_genefiltered.rds")
obj
head(obj@meta.data)

genemat <- obj@assays$RNA$counts # Take raw counts
dim(genemat)
head(rownames(genemat))
head(colnames(genemat))
"DCC" %in% rownames(genemat)
"RERE" %in% rownames(genemat)

grnamat <- obj@assays$crispr$counts # Take raw counts
dim(grnamat)
rownames(grnamat)

meta <- obj@meta.data

guides <- data.frame(rownames(grnamat))
colnames(guides) <- c("id")
guides$target_gene <- NA
guides$target_var <- NA
for(i in 1:nrow(guides)){
  
  id <- guides[i, "id"]
    
  if(grepl("DCC", id)){
    guides[i, "target_gene"] <- "DCC"
  } else if(grepl("RERE", id)){
    guides[i, "target_gene"] <- "RERE"
  } else if (grepl("non-targeting", id)){
    guides[i, "target_gene"] <- "NT"
  }
  
  if(grepl("rs", id)){
    rsnum <- unlist(strsplit(id, "rs"))[2]
    guides[i, "target_var"] <- paste0("rs", rsnum)
  } else{
    guides[i, "target_var"] <- "none"
  }
  
}

# Specify target gene, target variant (variant-of-interest, VOI), and perturbation threshold (this can be changed)
targetgene <- "RERE"
targetvar <- "rs301804"
perturb_threshold <- 3 # Minimum counts required to be defined 'perturbed'

######### Control cell definition

# All gRNAs targeting the target gene = 0
grna_exclude <- guides[guides$target_gene==targetgene, "id"]
grna_exclude
grnamat2 <- grnamat[rownames(grnamat) %in% grna_exclude,]
dim(grnamat2)
grnamat3 <- grnamat2[, colSums(grnamat2)==0]
dim(grnamat3)
grnamat3
ctrlcells <- colnames(grnamat3)


######### Perturb cell definition

grna_include <- guides[guides$target_var==targetvar, "id"]
grna_include
grnamat4 <- grnamat[rownames(grnamat) %in% grna_include,]
dim(grnamat4)

grnamat5 <- grnamat4[, apply(grnamat4 >= perturb_threshold, 2, any)] # At least one out of all considered gRNAs should have minimum 3 counts
dim(grnamat5)
grnamat5
pertcells <- colnames(grnamat5)

# The big input table including all genes surrounding DCC and RERE
inputdf <- meta
genes2test <- c("DCC", "MBD2", "POLI", "C18orf54", "RERE", "ENO1", "H6PD", "SPSB1", "SLC45A1", "ERRFI1", "PARK7")
genemat_genes2test <- genemat[rownames(genemat) %in% genes2test,]
dim(genemat_genes2test)
# inputdf <- cbind(inputdf, t(genemat_genes2test))
# 
# var2test <- c("rs4513167", "rs6508210", "rs4614799", "rs8089270", "rs301804")
# for(i in var2test){
#   inputdf[[i]] <- NA
# }

inputdf[rownames(inputdf) %in% ctrlcells, targetvar] <- "ctrl"
inputdf[rownames(inputdf) %in% pertcells, targetvar] <- "pert"

saveRDS(inputdf, file = "/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Glm+Permutation_Analysis_Results/Publication_FinalPermutations/inputdf.rds")

inputdf <- readRDS("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Glm+Permutation_Analysis_Results/Publication_FinalPermutations/inputdf.rds")
#####NT guide portion#####
# Specify target variant (variant-of-interest, VOI), target gene, and perturbation threshold
# NOTE: If targetvar & targetgene is given from command line input, DO NOT specify here
targetvar <- "none"
targetgene <- "RERE"
perturb_threshold <- 3 # Minimum counts required to be defined 'perturbed'

######### Control cell definition

# All gRNAs targeting the target gene = 0
grna_exclude <- guides[guides$target_gene==targetgene, "id"]
grna_exclude
grnamat2 <- grnamat[rownames(grnamat) %in% grna_exclude,]
dim(grnamat2)
grnamat3 <- grnamat2[, colSums(grnamat2)==0]
dim(grnamat3)
grnamat3
ctrlcellsNT <- colnames(grnamat3)


#Filter further for no NT guide-cells
grnamat3.1 <- grnamat[,ctrlcellsNT]
dim(grnamat3.1)
NTgrnas <- guides[guides$target_var==targetvar,] %>% pull(id)
NTgrnas
grnamat3.2 <- grnamat3.1[NTgrnas,]
dim(grnamat3.2)
grnamat3.3 <- grnamat3.2[,colSums(grnamat3.2)==0]
dim(grnamat3.3)
ctrlcellsNT <- colnames(grnamat3.3)
length(ctrlcellsNT)

######### Perturb cell definition

grna_include <- guides[guides$target_var==targetvar, "id"]
grna_include
grnamat4 <- grnamat[rownames(grnamat) %in% grna_include,]
dim(grnamat4)
grnamat5 <- grnamat4[, apply(grnamat4 >= perturb_threshold, 2, any)] # At least one out of all gRNAs in grna_include should have minimum 3 counts
dim(grnamat5)
grnamat5
pertcellsNT <- colnames(grnamat5)
length(pertcellsNT)

#Filter further for other gRNA targets to 'targetgene'
grnamat5.1 <- grnamat[,pertcellsNT]
dim(grnamat5.1)
targetgrnas <- guides[guides$target_gene==targetgene,] %>% pull(id)
targetgrnas
grnamat5.2 <- grnamat5.1[targetgrnas,]
dim(grnamat5.2)
rownames(grnamat5.2)
grnamat5.3 <- grnamat5.2[,colSums(grnamat5.2)==0]
dim(grnamat5.3)
pertcellsNT <- colnames(grnamat5.3)
length(pertcellsNT)

# NT2test <- c("DCC_NT", "RERE_NT")
# for(i in NT2test){
#   inputdf[[i]] <- NA
# }

inputdf[rownames(inputdf) %in% ctrlcellsNT, "DCC_NT"] <- "ctrl"
inputdf[rownames(inputdf) %in% pertcellsNT, "DCC_NT"] <- "pert"

saveRDS(inputdf, file = "/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Glm+Permutation_Analysis_Results/Publication_FinalPermutations/inputdf_final.rds")
        