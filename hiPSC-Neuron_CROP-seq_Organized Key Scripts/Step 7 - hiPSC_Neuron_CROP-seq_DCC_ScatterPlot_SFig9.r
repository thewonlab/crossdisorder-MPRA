library(readr)
library(dplyr)
library(tibble)
library(ghibli)
ghibli_palettes$PonyoLight
#install.packages("readxl")
library(readxl)
library(ggplot2)


rm(list = ls())

#MPRA allelic color: #F4E3D3FF
#emVar color: #F4ADB3FF

setwd("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Cell_Rebuttal+Publication_Analyses/ScatterPlot/Plots")

#load MPRA element data
file_path <- "/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Cell_Rebuttal+Publication_Analyses/ScatterPlot/Data/STable1.xlsx"
mpra_elem <- read_excel(file_path, sheet = 2)

cropseqdat <- read_csv("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Cell_Rebuttal+Publication_Analyses/Publication_FinalPermutations/Data/result_with_perm.csv") #VOIs have to be >3, but other VOIs will be included in pertcells

#Remove the first column
cropseqdat <- cropseqdat[, -1]

#Add MPRA type information to result_dataframe
# Define the active variants
active_variants <- c('rs1934138', 'rs4513167', 'rs6681362', 'rs4614799')

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

# Merge the dataframes on 'rsid'
merged_df <- merge(mpra_elem, cropseqdat, by.x = "RSID", by.y = "rsid") #use this one if wanting to plot glm coefficient
#merged_df <- merge(mpra_elem, cropseqdat_adj, by.x = "RSID", by.y = "rsid") #use this one if wanting to plot p values

# Filter the merged dataframe to include only 'DCC' in the 'gene' column
merged_df <- merged_df %>% filter(gene %in% c("DCC"))

# Define the rsid values to be labeled and colored
label_rsid <- c("rs4513167", "rs6508210", "rs4614799", "rs8089270") #rs301804
emVar <- "#F4ADB3FF"
allelic <- "#F4E3D3FF"

# Create a new column for colors based on RSID
merged_df$color <- ifelse(merged_df$RSID %in% c("rs4513167", "rs4614799"), emVar, 
                          ifelse(merged_df$RSID %in% c("rs6508210", "rs8089270"), allelic, "black"))

# Create the scatter plot
p <- ggplot(merged_df, aes(x = MPRA_logFC, y = glm_coef)) +
  geom_point(aes(color = color)) +  # Use the color column for coloring points
  geom_text(data = subset(merged_df, RSID %in% label_rsid), 
            aes(label = RSID), 
            vjust = -1, 
            color = "black") +  # Label only the specified RSID values
  scale_color_identity() +  # Use the colors as they are
  theme_minimal() +
  labs(title = "Scatterplot of MPRA_logFC vs glm_coef",
       x = "MPRA_logFC",
       y = "glm_coef")

# Print the plot
print(p)
