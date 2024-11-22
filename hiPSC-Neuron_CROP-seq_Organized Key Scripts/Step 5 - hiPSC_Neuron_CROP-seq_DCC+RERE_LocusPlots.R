#Customized by Alejandro Gomez from Jiseok Lee & Jessica McAfee's LocusPlot codes 6/25/24
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)

options(stringsAsFactors=F)
library(rtracklayer)
library(data.table)
library(readr)
library(dplyr)
library(GenomicRanges)
library(plotgardener)
library(tibble)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ghibli)
ghibli_palettes$PonyoLight

rm(list = ls())

# mpra negative color: #ECD89DFF
# mpra positive color: #ADB7C0FF
# mpra active color: #F4ADB3FF

# Set a directory to save file
setwd("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Cell_Rebuttal+Publication_Analyses/Locus_Plots/Plots")

# Load gene annotations (based on hg19)
load("/proj/hyejunglab/chr/geneAnno_allgenes.rda")
load("/proj/hyejunglab/MPRA/RNAseq/PSYCH/PJD/Jessica/mpra_abc.Rda")

mpra_abc_hgnc = mpra_abc %>%
  left_join(geneAnno1[,1:2], by=c("ensg"="ensembl_gene_id"))

########ATTENTION#######: First define locus to plot & axis scales
locusname <- "DCC" # or RERE
locus <- mpra_abc_hgnc %>%
  filter(hgnc_symbol==locusname)
locus_chrom <- geneAnno1[geneAnno1$hgnc_symbol==locusname, "chromosome_name"]
locus_chrom <- paste0("chr", locus_chrom)
locus_pos <- geneAnno1[geneAnno1$hgnc_symbol==locusname, "start_position"]
locus_range <- 10^6/4
chrstart <- locus_pos - locus_range
chrend <- geneAnno1[geneAnno1$hgnc_symbol==locusname, "end_position"] + locus_range #locus_pos + locus_range

######ATTENTION########: Second define variants to highlight depending on locusname
#crophits_act = c("rs301804") #RERE MPRA-active variant
#crophits_pos = c() #RERE MPRA-allelic variants

crophits_act = c("rs4513167", "rs4614799") #DCC MPRA-active variants
crophits_pos = c("rs6508210","rs8089270") #DCC MPRA-allelic variants




#Load and process all files needed for plotting
#GWAS data
load("/proj/hyejunglab/disorder/crossdisorder_philee/cdg3/finaldat/crossdis_GRanges_colocinput.rda")
gwasdat <- as_tibble(gwasranges)
gwasdat <- gwasdat[, c("seqnames", "start", "p.gwas", "rsid")]
colnames(gwasdat) <- c("chrom", "pos", "p", "snp")
gwasdat$chrom <- paste0("chr", gwasdat$chrom)

#MPRA data: allelic
combined_df <- read.table("/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/mpra_result_updated.tsv", header = TRUE, sep = "\t")
mpradat <- as_tibble(combined_df)
mpradat <- mpradat[, c("chr", "bp", "MPRA_P", "rsid")] # Use un-adjusted MPRA p-value
colnames(mpradat) <- c("chrom", "pos", "p", "snp")
mpradat$group = "MPRA_negative"
mpradat[mpradat$snp %in% crophits_act, "group"] <- "MPRA_active"
mpradat[mpradat$snp %in% crophits_pos, "group"] <- "MPRA_positive"
mpradat$group <- factor(mpradat$group, levels = c("MPRA_active", "MPRA_positive", "MPRA_negative"))
mpradat$chrom <- paste0("chr", mpradat$chrom)

#MPRA data: activity (FDR)
# load("/proj/hyejunglab/MPRA/RNAseq/PSYCH/PJD/Hyejung/MPRAactive_allelic_combination.rda")
# mpract = data.frame(mprange)[,c("seqnames","start","active_fdr","rsid")]
# colnames(mpract) <- c("chrom", "pos", "p", "snp")
# mpract$chrom = as.character(mpract$chrom)
# mpract$group = "MPRA_negative"
# mpract[mpract$snp %in% crophits_act, "group"] <- "MPRA_active"
# mpract[mpract$snp %in% crophits_pos, "group"] <- "MPRA_positive"
# mpract$group <- factor(mpract$group, levels = c("MPRA_active", "MPRA_positive", "MPRA_negative"))

#MPRA data: activity (logFC)
load("/proj/hyejunglab/MPRA/RNAseq/PSYCH/PJD/Hyejung/MPRAactive_allelic_combination.rda")
mpract = data.frame(mprange)[,c("seqnames","start","MPRA_logFC","rsid")]
colnames(mpract) <- c("chrom", "pos", "p", "snp")
mpract$chrom = as.character(mpract$chrom)
mpract$group = "MPRA_negative"
mpract[mpract$snp %in% crophits_act, "group"] <- "MPRA_active"
mpract[mpract$snp %in% crophits_pos, "group"] <- "MPRA_positive"
mpract$group <- factor(mpract$group, levels = c("MPRA_active", "MPRA_positive", "MPRA_negative"))

# Convert MPRA_logFC column to absolute values
mpract$p <- abs(mpract$p)

#Set the range values
mpra_p_max = as.integer(max(abs(-log10(locus$MPRA_P)))) + 4
#mpra_act_max = as.integer(max(abs(-log10(mpract[mpract$group=="MPRA_active","p"])))) + 3
mpra_act_max = as.integer(max(abs(locus$MPRA_logFC))) + 1
#mpra_act_max = as.integer(max(abs(-log10(mpract$p)))) + 2
gwas_max = as.integer(max(abs(-log10(locus$GWAS_P)))) + 2

#Merge gwasdat and mpradat to get MPRA 'group' to include active variants in the gwas plot
gwasdat$group <- "MPRA_negative"
gwasdat[gwasdat$p<5e-8,"group"] <- "GWS"
gwasdat[gwasdat$snp %in% crophits_act, "group"] <- "MPRA_active"
gwasdat[gwasdat$snp %in% crophits_pos, "group"] <- "MPRA_positive"
gwasdat$group <- factor(gwasdat$group, levels = c("MPRA_active", "MPRA_positive", "MPRA_negative", "GWS"))

# Plotgardener
pdf(paste0("alleliclogFC_cropseq_analysis_locusplot_7-1-24_", locusname,".pdf"), width=8.5, height=15)
pageCreate(width = 8.5, height = 15, default.units = "inches", showGuides = FALSE)
panelwidth <- 7
panelheight <- 0.7
dotsize <- 0.3

########## 2) GWAS Manhattan plot
p1 <- plotManhattan(
  data = gwasdat, 
  chrom = locus_chrom,
  chromstart = chrstart,
  chromend = chrend,
  assembly = "hg19",
  fill = colorby("group",
                 palette = colorRampPalette(c("#F4ADB3FF", "#F4E3D3FF", "#ADB7C0FF", "#AADAC5"))),
  cex = dotsize,
  sigLine = FALSE,
  range = c(0 ,gwas_max),
  x = 1, y = 1, width = panelwidth, height = panelheight, just = c("left", "top"), default.units = "inches"
)

annoYaxis(
  plot = p1, 
  at = c(0, gwas_max/2, gwas_max),
  axisLine = TRUE, 
  fontsize = 8,
)

plotText(
  label = "GWAS\n-log10(p)", 
  x = 0.5, y = 1 + panelheight/2, rot = 90, fontsize = 8, just = "center", default.units = "inches"
)

########## 2) MPRA Manhattan plot

# Plot all groups
p2 <- plotManhattan(
  data = mpradat, 
  chrom = locus_chrom,
  chromstart = chrstart,
  chromend = chrend,
  assembly = "hg19",
  fill = colorby("group",
                 palette = colorRampPalette(c("#F4ADB3FF", "#F4E3D3FF", "#ADB7C0FF"))),
  cex = dotsize,
  sigLine = FALSE,
  #sigCol = "#F4ADB3FF",
  range = c(0 , mpra_p_max),
  x = 1, y = 2, width = panelwidth, height = panelheight, just = c("left", "top"), default.units = "inches"
)

annoYaxis(
  plot = p2, 
  at = c(0, mpra_p_max/2, mpra_p_max),
  axisLine = TRUE, 
  fontsize = 8,
)

plotText(
  label = "MPRA allelic\n-log10(p)", 
  x = 0.5, y = 2 + panelheight/2, rot = 90, fontsize = 8, just = "center", default.units = "inches"
)

########## 3) MPRA logFC Manhattan plot
# Plot all groups
p3 <- plotManhattan(
  data = mpract, 
  chrom = locus_chrom,
  chromstart = chrstart,
  chromend = chrend,
  assembly = "hg19",
  fill = colorby("group",
                 palette = colorRampPalette(c("#F4ADB3FF", "#F4E3D3FF", "#ADB7C0FF"))),
  cex = dotsize,
  sigLine = FALSE,
  range = c(0 , mpra_act_max),
  x = 1, y = 3, width = panelwidth, height = panelheight, just = c("left", "top"), default.units = "inches",
  trans=""
)

annoYaxis(
  plot = p3, 
  at = c(0, mpra_act_max/2, mpra_act_max),
  axisLine = TRUE, 
  fontsize = 8,
)

plotText(
  label = "|MPRA allelic\nlogFC|", 
  x = 0.5, y = 3 + panelheight/2, rot = 90, fontsize = 8, just = "center", default.units = "inches"
)

########## 4) Gene track
highlightGene <- data.frame(matrix(nrow = 0, ncol = 2))
highlightGene <- rbind(highlightGene, c(locusname, "#F4ADB3FF"))
colnames(highlightGene) <- c("gene", "color") # To highlight the gene of interest

p4 <- plotGenes(
  chrom = locus_chrom,
  chromstart = chrstart,
  chromend = chrend,
  assembly = "hg19",
  fontsize = 6,
  geneHighlights = highlightGene,
  x = 1, y = 4, width = panelwidth, height = panelheight, just = c("left", "top"), default.units = "inches"
)

plotGenomeLabel(
  chrom = locus_chrom, chromstart = chrstart, chromend = chrend,
  assembly = "hg19",
  x = 1, y = 5, length = panelwidth, height = panelheight, scale = "Mb",
  just = c("left", "top"), default.units = "inches"
)

dev.off()