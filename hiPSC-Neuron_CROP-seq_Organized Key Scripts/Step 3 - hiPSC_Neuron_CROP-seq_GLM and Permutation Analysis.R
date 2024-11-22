options(stringsAsFactors = F)
library(dplyr)
library(ggplot2)
library(MASS)
library(reshape2)
# library(Seurat)
library(parallel)
library(Matrix)

working_dir = "/proj/hyejunglab/cropseq/Alejandro/CD_CROP-Seq_AG/CD_CROP-Seq_REAL_4-4-24/Glm+Permutation_Analysis_Results/Publication_FinalPermutations"
setwd(working_dir)
n_perm = 2000
data = readRDS("Data/inputdf_final.rds")


Y_target1 = c("DCC","MBD2","POLI","C18orf54")
Y_target2 = c("PARK7","ERRFI1","SLC45A1","RERE","ENO1","H6PD","SPSB1")
Y_target3 = c(Y_target1,Y_target2)

X_target1 = c("rs4513167","rs6508210","rs4614799","rs8089270")
X_target2 = c( "rs301804")
X_target3 = c("DCC_NT","RERE_NT")

xy_comb <- expand.grid(Y_target1, X_target1)
xy_comb = rbind(xy_comb,expand.grid(Y_target2, X_target2))
xy_comb = rbind(xy_comb,expand.grid(Y_target3, X_target3))
colnames(xy_comb) <- c("Y_target", "X_target")
xy_comb[,1] = as.vector(xy_comb[,1])
xy_comb[,2] = as.vector(xy_comb[,2])
write.csv(xy_comb,"Data/perm_table/xy_comb.csv")

fit_glm_nb = function(data,x,y){
    data[[x]] = sample(data[[x]])
    formula = as.formula(paste0(y,"~",x,"+nCount_RNA+nCount_crispr+percent.mt++factor(batch_num)"))
    glm_result = glm.nb(formula=formula,data = data)
    glm_coef <- coef(summary(glm_result))[2,1]
    glm_p <- coef(summary(glm_result))[2,4]
    glm_converged = glm_result$converged
    return(c(glm_coef,glm_p,glm_converged))
}

for (i in 1:dim(xy_comb)[1]){
    result <- data.frame(Gene = character(), Variant = character(), coef = numeric(), p_value = numeric(), perm_round = integer(), stringsAsFactors = FALSE)
    x = xy_comb[i,2]
    y = xy_comb[i,1]

    lines = mclapply(1:n_perm,function(k) fit_glm_nb(data,x,y),mc.cores=128)
    for (k in seq_along(lines)) {
      line <- c(y, x, lines[[k]], k)
      result <- rbind(result, data.frame(Gene = line[1], Variant = line[2], coef = as.numeric(line[3]), p_value = as.numeric(line[4]), converged = line[5], perm_round = as.integer(line[6]), stringsAsFactors = FALSE))
    }
    write.csv(result,paste0("Data/perm_table/glm_nb_results",i,".csv"),row.names = FALSE)
}
