#!/usr/bin/Rscript

### 2020.2.19

combined_matrix <- readRDS("input_files/tf_allgenes_normalized_matrix_2020-01-10.rds")
# spars_gene <- readRDS("input_files/PGP1_QCed_raw_counts_sparse.rds")
# matrix_gene <- readRDS("input_files/PGP1_QCed_raw_counts_df.rds")
# cor_matrix <- readRDS("input_files/pearson_cor_pgp1_2019-12-26.rds")



#index <- c(10660, 22977, 10745, 29985, 26162, 1938)
#index <- 41:33734
index <- readRDS("SCT_genes_index.rds")


for (i in index)  #41:33734
{
  mat_nor <- t(combined_matrix[c(1:40, i),])  
  col_genename <- colnames(mat_nor)[41]
  colnames(mat_nor)[41] <- "gene1" # RP11-34P13.3 
  mat_nor <- as.data.frame(mat_nor)
  
  
  fit4 <- lm(gene1 ~ (tf_BHLHE41 +tf_CEBPA+tf_BATF+tf_TAL1+tf_BATF3+
                        tf_FOS+tf_MAFF+tf_MNDA+tf_SPI1+tf_CEBPB+
                        tf_MAFB+tf_HCLS1+tf_IRF8+tf_EGR2+tf_IRF1+
                        tf_CIITA+tf_IFI16+tf_MEF2C+tf_JUN+tf_IRF5+
                        tf_ZFP36+tf_JUNB)^4 , data = mat_nor)
  
  
  gene1_lm_res_fit4 <- data.frame(coefficients=summary(fit4)$coefficients[,1], pvalue= summary(fit4)$coefficients[,4])
  gene1_lm_res_sel_fit4 <- subset(gene1_lm_res_fit4, pvalue < 0.01)
  gene1_lm_res_single_fit4 <- gene1_lm_res_fit4[1:23,]
  write.table(gene1_lm_res_sel_fit4, file = paste("output_files_3000/", col_genename, "_lm_res_fit4.txt", sep = ""), quote = F, sep = "\t", row.names = T)
  write.table(gene1_lm_res_single_fit4, file = paste("output_files_3000/", col_genename, "_lm_res_singletf_fit4.txt", sep = ""), quote = F, sep = "\t", row.names = T)
  
  
  
  fit3 <- lm(gene1 ~ (tf_BHLHE41 +tf_CEBPA+tf_BATF+tf_TAL1+tf_BATF3+
                        tf_FOS+tf_MAFF+tf_MNDA+tf_SPI1+tf_CEBPB+
                        tf_MAFB+tf_HCLS1+tf_IRF8+tf_EGR2+tf_IRF1+
                        tf_CIITA+tf_IFI16+tf_MEF2C+tf_JUN+tf_IRF5+
                        tf_ZFP36+tf_JUNB)^3 , data = mat_nor)
  
  
  gene1_lm_res_fit3 <- data.frame(coefficients=summary(fit3)$coefficients[,1], pvalue= summary(fit3)$coefficients[,4])
  gene1_lm_res_sel_fit3 <- subset(gene1_lm_res_fit3, pvalue < 0.01)
  gene1_lm_res_single_fit3 <- gene1_lm_res_fit3[1:23,]
  write.table(gene1_lm_res_sel_fit3, file = paste("output_files_3000/", col_genename, "_lm_res_fit3.txt", sep = ""), quote = F, sep = "\t", row.names = T)
  write.table(gene1_lm_res_single_fit3, file = paste("output_files_3000/", col_genename, "_lm_res_singletf_fit3.txt", sep = ""), quote = F, sep = "\t", row.names = T)
  
  
  fit2 <- lm(gene1 ~ (tf_BHLHE41 +tf_CEBPA+tf_BATF+tf_TAL1+tf_BATF3+
                        tf_FOS+tf_MAFF+tf_MNDA+tf_SPI1+tf_CEBPB+
                        tf_MAFB+tf_HCLS1+tf_IRF8+tf_EGR2+tf_IRF1+
                        tf_CIITA+tf_IFI16+tf_MEF2C+tf_JUN+tf_IRF5+
                        tf_ZFP36+tf_JUNB)^2 , data = mat_nor)
  
  
  gene1_lm_res_fit2 <- data.frame(coefficients=summary(fit2)$coefficients[,1], pvalue= summary(fit2)$coefficients[,4])
  gene1_lm_res_sel_fit2 <- subset(gene1_lm_res_fit2, pvalue < 0.01)
  gene1_lm_res_single_fit2 <- gene1_lm_res_fit2[1:23,]
  write.table(gene1_lm_res_sel_fit2, file = paste("output_files_3000/", col_genename, "_lm_res_fit2.txt", sep = ""), quote = F, sep = "\t", row.names = T)
  write.table(gene1_lm_res_single_fit2, file = paste("output_files_3000/", col_genename, "_lm_res_singletf_fit2.txt", sep = ""), quote = F, sep = "\t", row.names = T)
  
  
  
  fit1 <- lm(gene1 ~ tf_BHLHE41 +tf_CEBPA+tf_BATF+tf_TAL1+tf_BATF3+
               tf_FOS+tf_MAFF+tf_MNDA+tf_SPI1+tf_CEBPB+
               tf_MAFB+tf_HCLS1+tf_IRF8+tf_EGR2+tf_IRF1+
               tf_CIITA+tf_IFI16+tf_MEF2C+tf_JUN+tf_IRF5+
               tf_ZFP36+tf_JUNB , data = mat_nor)
  
  
  gene1_lm_res_fit1 <- data.frame(coefficients=summary(fit1)$coefficients[,1], pvalue= summary(fit1)$coefficients[,4])
  gene1_lm_res_sel_fit1 <- subset(gene1_lm_res_fit1, pvalue < 0.01)
  gene1_lm_res_single_fit1 <- gene1_lm_res_fit1[1:23,]
  write.table(gene1_lm_res_sel_fit1, file = paste("output_files_3000/", col_genename, "_lm_res_fit1.txt", sep = ""), quote = F, sep = "\t", row.names = T)
  write.table(gene1_lm_res_single_fit1, file = paste("output_files_3000/", col_genename, "_lm_res_singletf_fit1.txt", sep = ""), quote = F, sep = "\t", row.names = T)
  
  
}




