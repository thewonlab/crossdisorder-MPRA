####Lollipop Plot for rsid's targeted with guides:
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(ghibli)
ghibli_palettes$PonyoLight

rm(list = ls())

setwd("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Cell_Rebuttal+Publication_Analyses/Publication_FinalPermutations/Plots")

#Load in your data 
cropseqdat <- read_csv("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Cell_Rebuttal+Publication_Analyses/Publication_FinalPermutations/Data/result_with_perm.csv") #VOIs have to be >3, but other VOIs will be included in pertcells

#Remove the first column
cropseqdat <- cropseqdat[, -1]

#Add MPRA type information to result_dataframe
# Define the active variants
active_variants <- c('rs301804', 'rs1934138', 'rs4513167', 'rs6681362', 'rs4614799')

# Add the new column 'MPRA_Status' to 'result_dataframe'
cropseqdat$MPRA_Status <- ifelse(cropseqdat$Variant %in% active_variants, 'Active', 'Allelic')

# Renaming the columns for plot
names(cropseqdat)[names(cropseqdat) == "Variant"] <- "rsid" 
names(cropseqdat)[names(cropseqdat) == "coef"] <- "glm_coef" 
names(cropseqdat)[names(cropseqdat) == "p-value"] <- "glm_p" 
names(cropseqdat)[names(cropseqdat) == "Gene"] <- "gene" 

# Ensure columns of data values are numeric
cropseqdat$glm_p <- as.numeric(cropseqdat$glm_p)
cropseqdat$glm_coef <- as.numeric(cropseqdat$glm_coef)

# Given this is CRISPRi, do a one-sided analysis
# Make a copy of the dataframe for manipulation
cropseqdat_adj <- cropseqdat

# Adjust p-values for one-sided test for negative coefficients
cropseqdat_adj$glm_p_left_sided <- ifelse(cropseqdat_adj$glm_coef < 0, cropseqdat_adj$glm_p / 2, 1-cropseqdat_adj$glm_p / 2)

# Create the signed glm p-value column
cropseqdat_adj$signed_glm_p <- -1*log10(cropseqdat_adj$glm_p_left_sided) * sign(cropseqdat_adj$glm_coef)

# Apply FDR correction using the Benjamini-Hochberg method
cropseqdat_adj$glm_FDR <- p.adjust(cropseqdat_adj$glm_p_left_sided, method = "BH")

# Create a new column 'Significance'
cropseqdat_adj$Significance <- ifelse(cropseqdat_adj$glm_FDR < 0.05, TRUE, FALSE) #Use if using permuted dataframe glm p

#Subset rsids to plot and save as csv
#cropseqdat_subset <- subset(cropseqdat_adj, rsid %in% c('rs4513167', 'rs4614799', 'rs6508210', 'rs8089270', 'rs301804'))

# Load gene annotations and information file and modify to filter for protein-coding and make TSS dependent on strand (based on hg19)
load("/proj/hyejunglab/chr/geneAnno_allgenes.rda")
geneAnno = geneAnno1[geneAnno1$gene_biotype=="protein_coding",]
geneAnno$start = ifelse(geneAnno$strand==1, geneAnno$start_position, geneAnno$end_position)

#Include the 'start' values from geneAnno onto cropseqdat based on gene name/symbol
cropseqdat_genepos <- left_join(cropseqdat_adj, geneAnno %>% 
                                  dplyr::select(hgnc_symbol, ensembl_gene_id, start), 
                                by = c("gene" = "hgnc_symbol"))

################################
##Gene-Centered Plot####
# Filter your gene of interest
gene_of_interest <- "RERE" #Change as needed for your gene of interest one at a time (RERE/DCC)
dcc_VOIs <- c("rs4513167", "rs4614799", "rs6508210", "rs8089270") #Change as needed for your variants of interest

#rere_VOIs <- c("rs301804", "rs6681362", "rs1809332", "rs1463052") #Change as needed for your variants of interest
#rere_VOIs <- c("rs301804", "rs1463052") #Change as needed for your variants of interest
rere_VOIs <- c("rs301804") #Change as needed for your variants of interest

# Select the appropriate vector of variants based on gene_of_interest
if (gene_of_interest == "RERE") {
  variants_of_interest <- rere_VOIs
} else if (gene_of_interest == "DCC") {
  variants_of_interest <- dcc_VOIs
} else {
  stop("Invalid gene_of_interest value.")
}


#Now filter for only the gene of interest
genecentered_data <- cropseqdat_adj[cropseqdat_adj$gene == gene_of_interest & cropseqdat_adj$rsid %in% variants_of_interest, ]

# Color mapping for MPRA_status
color_mapping <- c("Allelic" = "gray", "Active" = "#F4ADB3FF")
shape_mapping <- c("TRUE" = 8, "FALSE" = 16)

# Create lollipop plot with color based on MPRA_status and annotations for significance
rsid_plot <- ggplot(genecentered_data, aes(x = rsid, y = signed_glm_p, color = MPRA_Status)) +
  geom_segment(aes(xend = rsid, yend = 0), color = "grey") +  # Add vertical lines
  geom_point(aes(color = MPRA_Status, shape = Significance), size = 3) +  # Include color and shape mappings here
  scale_color_manual(values = color_mapping, name = "MPRA Status", labels = c("Active", "Allelic")) + # Specify colors and add legend
  #scale_shape_manual(values = shape_mapping, name = "Permuted P<0.05") +
  scale_shape_manual(values = shape_mapping, name = "FDR<0.05") +
  theme_minimal() +
  labs(title = paste("GLM Coefficients for variants in", gene_of_interest, "locus"),
       x = "Variant Target",
       y = "Signed GLM P-Value") +  # Update y-axis label
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print the plot
print(rsid_plot)

################################
##Variant-Centered Plot####


###Lollipop Plot for downregulation of nearby genes for the top-performing variant rsid########
# Replace 'your_top_rsid' with the actual top-performing rsid
top_variant_rsid <- c("rs301804") #you can get this from the lollipop plots above for RERE (rs301796,rs301804) or DCC (rs4372757,rs4513167,rs8089270)
top_variant_target <- "RERE"
gene_data_subset <- cropseqdat_genepos[cropseqdat_genepos$rsid %in% top_variant_rsid, ]

#Save as PDF in your working directory
#pdf(paste0("Cropseq_genesnearby_LollipopPlot_", top_variant_rsid, "(", top_variant_target, ")", ".pdf"))

# Omit NAs in the gene data
# gene_data_subset <- na.omit(gene_data_subset)

# Order the data by genomic coordinates
gene_data_subset$start <- as.factor(gene_data_subset$start)
gene_data_subset <- gene_data_subset[order(gene_data_subset$start), ]
gene_data_subset$gene <- factor(
  gene_data_subset$gene, 
  levels = unique(gene_data_subset$gene[order(gene_data_subset$start)])
)

# Create a lollipop plot for downregulation of nearby genes with dashed threshold line
gene_plot <- ggplot(gene_data_subset, aes(x = gene, y = signed_glm_p)) +
  geom_segment(aes(x = gene, xend = gene, y = 0, yend = signed_glm_p), color = "grey") +  # Set lines to gray
  geom_point(aes(color = MPRA_Status, shape = Significance), size = 3) +  # Include color and shape mappings here
  scale_color_manual(values = color_mapping, name = "MPRA Status", labels = c("Active", "Allelic")) + 
  scale_shape_manual(values = shape_mapping, name = "FDR<0.05") +
  theme_minimal() +
  labs(title = paste("Gene Up/Downregulation near", top_variant_rsid[1], "(", top_variant_target, ")"),
       x = "Gene ordered by TSS location (Decreasing)",
       y = "Signed -log10P") +
  coord_flip()

# Show the plot
print(gene_plot)

