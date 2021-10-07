# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------                                    This is the R code used in Ramesh et al. 2021                                           ---------------
# ---------------                       Included here are both proteomic analysis, as well as small RNA analysis                             ---------------
# ---------------                 In all datasets, D2 equals 3 g/L of sugar, D4 equals 30 g/L and D6 equals 300 g/L                          ---------------
# ---------------               If you have any questions or concern, you can ge in touch with me at signe.skog@liu.se                       ---------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------






# ---------------                                            1.         PROTEOMICS                                                           ---------------







#----------------------------------------------------------------
#final notes
#For Rashmi data, we used limma with the model
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
#on her data which she corrected (working name was Rashmi corrected values)
#from this, we picked out the coef we wanted, D2 and D6, and used BF as adj p val


library(ROTS)
library(limma)
library(readxl)
library(xlsx)

diet<-data.frame(row.names = c("D2", "D2_2","D2_3", "D4", "D4_2", "D4_3", "D6", "D6_2", "D6_3"),
                 Diet= rep(c("D2", "D4", "D6"), each=3),
                 Batch=c(1,2,2,3,4,4,5,6,6))
design<-model.matrix(~ diet$Diet)
#we use d4 as intercept
diet$Diet<-relevel(diet$Diet, "D4")



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

#color code-------------------------
col45<-c("#d56d64",
         "#e19355",
         "#f4dd91",
         "#4c9051",
         "#b1cecf",
         "#a790cb",
         "#6282be")
scales::show_col(col45)
