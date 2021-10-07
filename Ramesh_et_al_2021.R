# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------                                    This is the R code used in Ramesh et al. 2021                                           ---------------
# ---------------                       Included here are both proteomic analysis, as well as small RNA analysis                             ---------------
# ---------------                 In all datasets, D2 equals 3 g/L of sugar, D4 equals 30 g/L and D6 equals 300 g/L                          ---------------
# ---------------               If you have any questions or concern, you can ge in touch with me at signe.skog@liu.se                       ---------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------






# ---------------                                            1.         PROTEOMICS                                                           ---------------


#load libraries
library(ROTS)
library(limma)
library(readxl)
library(xlsx)


#Design
diet<-data.frame(row.names = c("D2", "D2_2","D2_3", "D4", "D4_2", "D4_3", "D6", "D6_2", "D6_3"),  Diet= rep(c("D2", "D4", "D6"), each=3))
design<-model.matrix(~ diet$Diet)
# relevel "Diet" to use D4 (30g/L) as intercept for model
diet$Diet<-relevel(diet$Diet, "D4")


#This is the desired Design
# # (Intercept) diet$DietD2 diet$DietD6
# 1           1           1           0
# 2           1           1           0
# 3           1           1           0
# 4           1           0           0
# 5           1           0           0
# 6           1           0           0
# 7           1           0           1
# 8           1           0           1
# 9           1           0           1
# attr(,"assign")
# [1] 0 1 1
# attr(,"contrasts")
# attr(,"contrasts")$`diet$Diet`
# [1] "contr.treatment"


#Clean up data
df<-read_excel("~/Rashmi correction factor_total counts.xlsx", sheet = 1)
df<-df[,2:10]

#Run model - we used Limma with Bonferroni corrected adjusted p values
fit<-lmFit(df, design)
fit<-eBayes(fit)
top<-topTable(fit, number=542)
topd2<-topTable(fit, coef="diet$DietD2", number=542)
topd6<-topTable(fit, coef="diet$DietD6", number=542)

topd2<-topd2[match(rownames(topd6), rownames(topd2)),]
top<-top[match(rownames(topd6), rownames(top)),]

#Save results
write.xlsx(x=top, file="C:/Users/sigsk47/Desktop/Total_top.xlsx")
write.xlsx(x=topd2, file="C:/Users/sigsk47/Desktop/D2_top.xlsx")
write.xlsx(x=topd6, file="C:/Users/sigsk47/Desktop/D6_top.xlsx")



# ---------------                                             2.      small RNA analysis                                                     ---------------







#  2.1 Lane merging, Adaptor trimming, Quality control                  ----------------        
#  2.2 Read Counting and Phenotypic information                         ----------------      
#  2.3 Annotation                                                       ----------------     
#  2.4 Data Filtering and Normalization                                 ----------------      

#HUMAN
pac_h<-PAC_filter(pac_h, size=c(18,50), threshold=10, coverage = 60, anno_target = list("genome", "mis0"))
pac_h<-PAC_norm(pac_h)
pac_h<-PAC_filter(pac_h, norm="cpm", threshold=10, coverage=25)

#FLY

pac_d<-PAC_filter(pac_d, size=c(18,50), threshold=10, coverage = 60, anno_target = list("genome", "mis0"))
pac_d<-PAC_norm(pac_d)
pac_d<-PAC_filter(pac_d, norm="cpm", threshold=10, coverage=25)

pac_d<-PAC_filter(pac_d, anno_target = list("Biotypes_mis0", c("lncRNA","miRNA","Mt_tRNA","no_anno","other","piRNA","protein_coding","snoRNA","snRNA","tRNA")))

pac_d<-PAC_summary(pac_d, norm="cpm", type="means", pheno_target = list("DI"))
pac_d<-PAC_summary(pac_d, norm="cpm", type="log2FC", pheno_target = list("DI"))
pac_d<-PAC_summary(pac_d, norm="cpm", type="log2FC", pheno_target = list("DI"), rev=TRUE)


#  2.5 miRNA analysis                                                   ----------------   


pac_miRNA_plot <- PAC_trna(pac_miRNA, norm="cpm", filter = NULL,
                        join = TRUE, top = 10, log2fc = TRUE,
                        pheno_target = list("DI", c("D4_C", "D6_C")),
                        anno_target_1 = list("mir"),
                        anno_target_2 = list("mir"))

pac_miRNA_plot_all <- PAC_trna(pac_miRNA, norm="cpm", filter = NULL,
                        join = TRUE, top = 41, log2fc = TRUE,
                        pheno_target = list("DI", c("D4_C", "D6_C")),
                        anno_target_1 = list("mir"),
                        anno_target_2 = list("mir"))

#  2.5.1 Negative Binomial distribution of miRNA

library(MASS)

df<-data.frame(pac_mir$norm$cpm, tr=pac_mir$Anno$mir)
df_l<-tidyr::gather(df, key="sample", value="cpm", -tr)

df_l$Diet<-rep(c("D4", "D6"), each=1170)
df_l$Int<-rep(c("C", "NAC", "C", "NAC"), each=585)
df_l$DI <- paste0(df_l$Diet, sep="-", df_l$Int)

#nb_func creates a data frame with negative binomial statistics for miRNA comparing D4 C vs D6 C

nb_func <- function(z){
  library(lme4)
  res_list <- list(NULL)
  res_short_list <- list(NULL)
  df <- data.frame(matrix(NA, ncol=7, nrow=length(unique(z$tr))))
  colnames(df) <- c("hypothesis","variable", "estimate", "se", "z", "p", "mt") 
  for(i in 1:length(unique(z$tr))){
    nam <- as.character(unique(z$tr)[i])
    data <- z[z$tr == as.character(nam),] 
    y <- data$cpm
    res <- summary(nb_res <-glm.nb(as.integer(y) ~ DI, data=data))
    res_list[[i]] <- res
    res_short_list[[i]] <- res$coefficient
    df$hypothesis[i] <- rownames(res$coefficient)[3]
    df$variable[i] <- nam
    df$estimate[i] <- res$coefficient[3,1]
    df$se[i] <- res$coefficient[3,2]
    df$z[i] <- res$coefficient[3,3]
    df$p[i] <- res$coefficient[3,4]
    names(res_list)[i] <- paste(nam, "nbmix", sep="_")
    names(res_short_list)[i] <- paste(nam, "nbmix", sep="_")
  }
  fin_list <- list(summary=df, short_result=res_short_list, long_result=res_list)
  return(fin_list)
}



#nb_func_NAC creates a data frame with negative binomial statistics for miRNA comparing D4 C vs D4 NAC

nb_func_NAC <- function(z){
  library(lme4)
  res_list <- list(NULL)
  res_short_list <- list(NULL)
  df <- data.frame(matrix(NA, ncol=7, nrow=length(unique(z$tr))))
  colnames(df) <- c("hypothesis","variable", "estimate", "se", "z", "p", "mt") 
  for(i in 1:length(unique(z$tr))){
    nam <- as.character(unique(z$tr)[i])
    data <- z[z$tr == as.character(nam),] 
    y <- data$cpm
    res <- summary(nb_res <-glm.nb(as.integer(y) ~ DI, data=data))
    res_list[[i]] <- res
    res_short_list[[i]] <- res$coefficient
    df$hypothesis[i] <- rownames(res$coefficient)[3]
    df$variable[i] <- nam
    df$estimate[i] <- res$coefficient[2,1]
    df$se[i] <- res$coefficient[2,2]
    df$z[i] <- res$coefficient[2,3]
    df$p[i] <- res$coefficient[2,4]
    names(res_list)[i] <- paste(nam, "nbmix", sep="_")
    names(res_short_list)[i] <- paste(nam, "nbmix", sep="_")
  }
  fin_list <- list(summary=df, short_result=res_short_list, long_result=res_list)
  return(fin_list)
}


nb_res_list <- nb_func(df_l)
nb_res_list_NAC <- nb_func_NAC(df_l)



#  2.6 tsRNA analysis                                                   ----------------   


#  2.6.2 Negative Binomial Model using functions from section 2.5.2

df<-data.frame(pac_mir$norm$cpm, tr=pac_mir$Anno$mir)
df_l<-tidyr::gather(df, key="sample", value="cpm", -tr)

df_l$Diet<-rep(c("D4", "D6"), each=1170)
df_l$Int<-rep(c("C", "NAC", "C", "NAC"), each=585)
df_l$DI <- paste0(df_l$Diet, sep="-", df_l$Int)

nb_res_list <- nb_func(df_l)
nb_res_list_NAC <- nb_func_NAC(df_l)

#  2.7 Data Visulaization Figure 3                                      ----------------   

Library(patchwork)
#color code-------------------------
col45<-c("#d56d64",
         "#e19355",
         "#f4dd91",
         "#4c9051",
         "#b1cecf",
         "#a790cb",
         "#6282be")
scales::show_col(col45)

#  Fig 3.B-pies

#  Fig 3.C-sizedist

#  Fig 3.D - heidi

#  Fig 3.E
pac_miRNA_plot$plots$Expression_Anno_1$Grand_means + pac_miRNA_plot$plots$Log2FC_Anno_1

#  Fig 3.F - miRNA

#  Fig 3.G - covplot mir-10

#  2.8 Data Visualization Figure 4                                      ---------------- 


#  Fig 4 B - piecharts mito vs nuc

#  Fig 4 C - heatmap

#  Fig 4 D - scatterplots

#  2.9 Data Visualization Supplementary Figures                         ---------------- 

#  Sup. Fig. S3 A - jitter diet 

#  Sup. Fig. S3 B - jitter diet vs nac 

#  Sup. Fig. S3 C - all miRNA 
pac_miRNA_plot_all$plots$Expression_Anno_1$Grand_means + pac_miRNA_plot_all$plots$Log2FC_Anno_1

#  Sup. Fig. S4 A -

#  Sup. Fig. S4 B -

#  Sup. Fig. S4 C -

#  2.9.
