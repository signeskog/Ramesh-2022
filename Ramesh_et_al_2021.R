# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------                                    This is the R code used in Ramesh et al. 2021                                           ---------------
# ---------------                       Included here are both proteomic analysis, as well as small RNA analysis                             ---------------
# ---------------                 In all datasets, D2 equals 3 g/L of sugar, D4 equals 30 g/L and D6 equals 300 g/L                          ---------------
# ---------------       If you have any questions or concern, you can ge in touch with me at signe.skog@liu.se or here at GitHub             ---------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------



# ---------------                                            1.         PROTEOMICS                                                           ---------------


#load libraries
library(ROTS)
library(limma)
library(readxl)
library(xlsx)


#Design
diet<-data.frame(row.names = c("D2", "D2_2","D2_3", "D4", "D4_2", "D4_3", "D6", "D6_2", "D6_3"),  Diet= rep(c("D2", "D4", "D6"), each=3))
# relevel "Diet" to use D4 (30g/L) as intercept for model
diet$Diet<-relevel(diet$Diet, "D4")
design<-model.matrix(~ diet$Diet)


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
df<-read_excel("/home/sigsk47/Documents/Data/proteomics/Rashmi correction factor_total counts.xlsx", sheet = 1)
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



# ---------------                                            2.      small RNA analysis                                                     ---------------

library(seqpac)
library(patchwork)
library(ggplot2)

#  2.1 Adaptor trimming, Quality control and Read Counting              ---------------


#  2.1.1 Drosphila PAC file                                             ---------------


counts_d<-make_counts(input = "/run/media/sigsk47/Elements/fasta/NAC",
                      trimming="seqpac",
                      threads=10,
                      plot=TRUE,
                      parse="default_neb", 
                      save_temp=FALSE)



pheno_d<-data.frame(row.names = colnames(counts_d$counts), names=colnames(counts_d$counts))
pac_d<-make_PAC(pheno = pheno_d, counts = counts_d, anno = NULL)


#  2.1.2 Human PAC file                                                   ------------


counts_h<-make_counts(input = "/run/media/sigsk47/Elements/hs2/Merged_fastq",
                      trimming="seqpac",
                      threads=1,
                      plot=TRUE,
                      parse="default_neb", 
                      save_temp=TRUE)
pheno_h<-data.frame(row.names = colnames(counts_h$counts), names=colnames(counts_h$counts))
pac_h<-make_PAC(pheno = pheno_h, counts = counts_h, anno = NULL)


#  2.3 Annotation                                                       ---------------     


#  2.3.1 Genomic annotation fly                                           ---------------


map_reanno(PAC=pac_d, ref_paths=list(genome="/run/media/sigsk47/Elements1/PhD/R/Genomes/Drosophila_melanogaster_ref/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome.fa"),
                                       #genome="/home/sigsk47/Documents/Genomes/fly/genome/EnsGenome/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"), 
           output_path="home/sigsk47/Documents/Genomes/out", 
           type="internal",
           mismatches=0, 
           import="genome", 
           threads=1, 
           keep_temp=TRUE)
reanno_genome <- make_reanno(reanno_path = "home/sigsk47/Documents/Genomes/out", 
                             PAC= pac_d, 
                             mis_fasta_check = TRUE)
pac_d <- add_reanno(reanno_genome, 
                    type="genome", 
                    genome_max="all", 
                    mismatches=0, 
                    merge_pac = pac_d)
pac_d <- PAC_filter(pac_d, 
                    anno_target = list("genome", c("mis0")))


#  2.3.2 Biotype annotation fly                                         ---------------

ref_b=list(Ensenmbl="/run/media/sigsk47/Elements1/PhD/R/Genomes/Drosophila_melanogaster/Drosophila_melanogaster/Ensembl/Drosophila_melanogaster.BDGP6.ncrna.fa",
           tRNA="/run/media/sigsk47/Elements1/Documents/PhD/R/Genomes/Drosophila_melanogaster/Drosophila_melanogaster/tRNA_reanno/tRNA_mature.fa",
           rRNA="/run/media/sigsk47/Elements1/PhD/R/Genomes/Drosophila_melanogaster/Drosophila_melanogaster/rRNA_reanno/drosophila_rRNA_all.fa",
           piRNA="/run/media/sigsk47/Elements1/PhD/R/Genomes/Drosophila_melanogaster/Drosophila_melanogaster/piRNA_piRBase/piR_dme.fa",
           miRNA="/run/media/sigsk47/Elements1/PhD/R/Genomes/Drosophila_melanogaster/Drosophila_melanogaster/miRNA/miRBase_21-dme.fa",
           protein_coding="/home/sigsk47/Documents/Genomes/fly/PC/DM_ens_BDGP6.32.cdna.fa")
map_reanno(PAC=pac_d, 
           ref_paths=ref_b, 
           output_path="/home/sigsk47/Documents/Genomes/out", 
           type="internal", 
           mismatches=0, 
           import="biotype", 
           threads=1, 
           keep_temp=F)
reanno_bio<- make_reanno(reanno_path="/home/sigsk47/Documents/Genomes/out", 
                         PAC=pac_d, 
                         mis_fasta_check = TRUE)
bio_search<- list(Ensenmbl=c("lncRNA", "lincRNA", "miRNA", "pre_miRNA", "rRNA", "snoRNA", 
                             "snRNA", "tRNA", "lncRNA","misc","ribozyme","scaRNA","scRNA","snoRNA","snRNA" 
                             ,"sRNA", "Mt_tRNA", "lncRNA","Mt_rRNA", "rRNA",
                             "5SrRNA", "28SrRNA", "18SrRNA", "CR41606", "CR40582", "5.8SrRNA",
                             "CR40611", "CR40679", "pre-rRNA", "2SrRNA", "CR40594", "CR40621"),
                  rRNA=c("18S_rRNA", "16S_rRNA", "12S_rRNA", "28S_rRNA", "Other_rRNA", "5S_rRNA",
                         "pre_S45_rRNA", "5.8S_rRNA"),
                  miRNA="dme-",
                  tRNA =c("tRNA", "mt-tRNA"),
                  piRNA="piR",
                  protein_coding=c("FB"))
pac_d <- add_reanno(reanno=reanno_bio, 
                   bio_search=bio_search, 
                   type="biotype", 
                   bio_perfect=TRUE, 
                   mismatches = 0, 
                   merge_pac=pac_d)
hierarchy <- list(rRNA="18S_rRNA|16S_rRNA|12S_rRNA|28S_rRNA|Other_rRNA|5S_rRNA|pre_S45_rRNA|5.8S_rRNA|Ensenmbl_rRNA",
                  Mt_tRNA="tRNA_mt-tRNA|Ensenmbl_Mt_tRNA",
                  tRNA="Ensenmbl_tRNA|tRNA_tRNA",
                  miRNA ="^miRNA|Ensenmbl_miRNA|Ensenmbl_pre_miRNA",
                  lncRNA="Ensenmbl_lincRNA|Ensenmbl_lncRNA",
                  piRNA="piRNA",
                  protein_coding="protein_coding")
pac_d <- simplify_reanno(input=pac_d, 
                         hierarchy=hierarchy, 
                         mismatches=0,
                         bio_name="Biotypes_mis0",
                         merge_pac = TRUE)


#  2.3.3 Genomic annotation Human                                         ---------------


map_reanno(PAC=pac_h, 
           ref_paths=list(genome="/home/sigsk47/Documents/Genomes/human/genome/hsref_GRCh38.fa"), 
           output_path="home/sigsk47/Documents/Genomes/out", 
           type="internal",
           mismatches=0, 
           import="genome", 
           threads=1, 
           keep_temp=TRUE)
reanno_genome<- make_reanno(reanno_path = "home/sigsk47/Documents/Genomes/out",  
                            PAC=pac_h, 
                            mis_fasta_check = TRUE,
                            output = "list")
pac_h <- add_reanno(reanno_genome, 
                    type="genome", 
                    genome_max="all", 
                    mismatches=0, 
                    merge_pac = pac_h)
pac_h<-PAC_filter(pac_h, anno_target = list("genome", c("mis0")))


#  2.3.4 Biotype annotation human                                        ---------------


ref_b=list(Ensenmbl="/home/sigsk47/Documents/Genomes/human/ncRNA/HS.GRCh38.ncRNA.fa",
           tRNA="/home/sigsk47/Documents/Genomes/human/tRNA_r/gtrnadb_hs_19.fa",
           rRNA="/home/sigsk47/Documents/Genomes/human/rRNA/rnacentral_all_rrna.fa",
           piRNA="/home/sigsk47/Documents/Genomes/human/piRNA/hsa_piRbase_2.0.fa",
           miRNA="/home/sigsk47/Documents/Genomes/human/miRNA/miRbase_hs.fa",
           protein_coding="/home/sigsk47/Documents/Genomes/human/pc/GRCh38_latest_rna.fa")
map_reanno(PAC=pac_h, 
           ref_paths=ref_b, 
           output_path="/home/sigsk47/Documents/Genomes/out", 
           type="internal", 
           mismatches=0, 
           import="biotype", 
           threads=1, 
           keep_temp=F)
reanno_bio<- make_reanno(reanno_path="/home/sigsk47/Documents/Genomes/out", 
                         PAC=pac_h, 
                         mis_fasta_check = TRUE,
                         output="list")
bio_search<- list(Ensenmbl=c("lncRNA", "lincRNA", "miRNA", "pre_miRNA", "rRNA", "snoRNA", 
                             "snRNA", "tRNA", "lncRNA","misc","ribozyme","scaRNA","scRNA","snoRNA","snRNA" 
                             ,"sRNA", "Mt_tRNA", "lncRNA","Mt_rRNA", "rRNA",
                             "5SrRNA", "28SrRNA", "18SrRNA", "CR41606", "CR40582", "5.8SrRNA",
                             "CR40611", "CR40679", "pre-rRNA", "2SrRNA", "CR40594", "CR40621"),
                  rRNA=c("rRNA", "ribosomal_RNA", "5S", "human"),
                  miRNA="hsa-",
                  tRNA =c("tRNA", "mt-tRNA"),
                  piRNA="piR",
                  protein_coding=c("NM", "NR", "XM", "XR"))
anno <- add_reanno(reanno_bio, 
                   bio_search=bio_search, 
                   type="biotype", 
                   bio_perfect=T, 
                   mismatches = 0)
hierarchy <- list(rRNA="Ensenmbl_rRNA|Ensenmbl_Mt_rRNA|^rRNA",
                  Mt_tRNA="tRNA_mt-tRNA|Ensenmbl_Mt_tRNA",
                  tRNA="Ensenmbl_tRNA|tRNA_tRNA",
                  miRNA ="^miRNA|Ensenmbl_miRNA|Ensenmbl_pre_miRNA",
                  lncRNA="Ensenmbl_lincRNA|Ensenmbl_lncRNA",
                  piRNA="piRNA",
                  protein_coding="protein_coding")
pac_h <- simplify_reanno(input=anno, 
                         hierarchy=hierarchy, 
                         mismatches=0, 
                         bio_name="Biotypes_mis0", 
                         merge_pac = pac_h)

save(pac_h, file="~/Human_may21.Rdata")


#  2.4 Data Filtering and Normalization                                 ----------------      


#  2.4.1 Drosophila data filtering and normalization


pac_d<-PAC_filter(pac_d, size=c(18,50), threshold=10, coverage = 60, anno_target = list("genome", "mis0"))
pac_d<-PAC_norm(pac_d)
pac_d<-PAC_filter(pac_d, norm="cpm", threshold=10, coverage=25)

pac_d<-PAC_filter(pac_d, anno_target = list("Biotypes_mis0", c("lncRNA","miRNA","Mt_tRNA","no_anno","other","piRNA","protein_coding","tRNA")))

pac_d<-PAC_summary(pac_d, norm="cpm", type="means", pheno_target = list("DI"))
pac_d<-PAC_summary(pac_d, norm="cpm", type="log2FC", pheno_target = list("DI"))
pac_d<-PAC_summary(pac_d, norm="cpm", type="log2FC", pheno_target = list("DI"), rev=TRUE)


#  2.4.1 Human data filtering                                           ----------------


pac_h<-PAC_filter(pac_h, size=c(18,50), threshold=10, coverage = 60, anno_target = list("genome", "mis0"))
pac_h<-PAC_norm(pac_h)
pac_h<-PAC_filter(pac_h, norm="cpm", threshold=10, coverage=25)


#  2.5 miRNA analysis                                                   ----------------   


#  2.5.1 miRNA mapping and data preparation                             ----------------


pac_miRNA<-PAC_filter(pac_d, 
                      anno_target=list("Biotypes_mis0",c("miRNA")))
map_object<-PAC_mapper(pac_miRNA, 
                       ref="/home/sigsk47/Documents/Genomes/fly/miRNA/miRBase_21-dme.fa", 
                       mismatches=0,
                       threads=1,
                       report_string = T)
map_object <-  map_object[!unlist(lapply(map_object, function(x){x[[2]][1,1] == "no_hits"}))]


miRNA_class<-function (PAC, map, terminal = 5) 
{
  logi_no_hits <- unlist(lapply(map, function(x) {
    x[[2]][1, 1] == "no_hits"
  }))
  if (any(logi_no_hits) == TRUE) {
    cat("The map object contains references without hits.", 
        "\nThese will be removed from output.")
    map <- map[!unlist(lapply(map, function(x) {
      x[[2]][1, 1] == "no_hits"
    }))]
    stopifnot(any(unlist(lapply(map, function(x) {
      x[[2]][1, 1] == "no_hits"
    }))) == FALSE)
  }
  type_vector <- lapply(map, function(x) {
    align <- x$Alignments
    reference_length <- x$Ref_seq@ranges@width
    ref_name <- x$Ref_seq@ranges@NAMES
    # terminal_type <- ifelse(align$Align_start <= terminal, 
    #                         "5'", ifelse(reference_length - align$Align_end <= 
    #                                        terminal, "3'", "i'"))
    # half_type <- ifelse(align$type_start_loop2 == TRUE | 
    #                       align$type_end_loop2 == TRUE, "half", "tRF")
    return(data.frame(miRNA_ref = ref_name, seq = rownames(align))) 
           #class = paste(terminal_type, half_type, sep = "-"), 
           #decoder = align$decoder, acceptor = align$acceptor))
  })
    type_vector <- do.call("rbind", type_vector)
    rownames(type_vector) <- NULL
    type_vector <- split(type_vector, type_vector$seq)
    finished <- lapply(type_vector, function(x) {
      df <- data.frame(seq = paste0(sort(unique(x$seq)), collapse = ";"), 
                       miRNA_ref = paste0(sort(unique(x$miRNA_ref)), collapse = ";"))
      return(df)
    })
    finished <- do.call("rbind", finished)
    PAC@Anno$seq <- rownames(PAC@Anno)
    pac_trna <- PAC_filter(PAC, anno_target = list("seq", finished$seq), 
                           subset_only = TRUE)
    stopifnot(identical(rownames(pac_trna@Anno), as.character(finished$seq)))
    pac_trna@Anno <- cbind(pac_trna@Anno[, 1, drop = FALSE], 
                           finished[, -1])
    return(pac_trna)
}
pac_mirna<-miRNA_class(pac_miRNA, map=map_object)

df<-pac_mirna@Anno

df<-df[match(rownames(pac_miRNA@Anno), rownames(df)),]
identical(rownames(df), rownames(pac_miRNA@Anno))

pac_miRNA@Anno$mir<-df$`finished[, -1]`
pac_miRNA@Anno$mir<-sub(" .*", "\\1", pac_miRNA@Anno$mir)

# pac_miRNA<-list(pheno=pac_d@Pheno, anno=as.data.frame(Book2$`Refined annotation`), norm=as.data.frame(Book2[,6:25]), summary=as.data.frame(Book2[,26:29]))
# 
# pac_miRNA<-asS4(pac_miRNA)
# pac_miRNA$Pheno<-pac_miRNA$pheno
# pac_miRNA$Anno<-pac_miRNA$anno
# pac_miRNA$norm$cpm<-pac_miRNA$norm

pac_miRNA_plot <- PAC_trna(pac_miRNA, norm="cpm", filter = NULL,
                           join = FALSE, top = 10, log2fc = TRUE,
                           pheno_target = list("DI", c("D4_C", "D6_C", "D4_NAC", "D6_NAC")),
                           anno_target_1 = list("mir"),
                           anno_target_2 = list("mir"))

pac_miRNA_plot_all <- PAC_trna(pac_miRNA, norm="cpm", filter = NULL,
                               join = FALSE, top = 41, log2fc = TRUE,
                               pheno_target = list("DI", c("D4_C", "D6_C", "D4_NAC", "D6_NAC")),
                               anno_target_1 = list("mir"),
                               anno_target_2 = list("mir"))


df_t<-data.frame(pac_miRNA@norm$cpm, pac_miRNA@Anno$mir)
df_t2<-aggregate(df_t, by=list(df_t$pac_miRNA.Anno.mir), FUN=mean)
df_t3<-as.data.frame(df_t2[,1:6])
df_t3$d4n<-rowMeans(df_t2[,7:11])
df_t3$d6<-rowMeans(df_t2[,12:16])
df_fc_d6<-

df_d4nac<-aggregate(pac_miRNA_plot_all$data$anno_target_1$D4_NAC$value, by=list(pac_miRNA_plot_all$data$anno_target_1$D4_NAC$Group.1), FUN=mean)


#  2.5.2 Negative Binomial distribution of miRNA


library(MASS)

df<-data.frame(pac_miRNA@norm$cpm, tr=pac_miRNA@Anno$mir)
df_l<-tidyr::gather(df, key="sample", value="cpm", -tr)

df_l$Diet<-sub("_.*", "\\1", df_l$sample)
df_l$DI <- sub('^([^_]+_[^_]+).*', '\\1', df_l$sample)


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
    df$hypothesis[i] <- rownames(res$coefficient)[2]
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


#  2.6.1 tsRNA mapping and data preparation                             ----------------

pac_tr<-PAC_filter(pac_d, 
                   anno_target=list("Biotypes_mis0", c("tRNA", "Mt_tRNA")))
map_object_tRNA <- PAC_mapper(pac_tr, 
                              ref="/home/sigsk47/Documents/Genomes/fly/tRNA/tRNA.fa", 
                              mismatches=0, 
                              threads=1, 
                              report_string = T)
map_object<-map_rangetype(map_object_tRNA, type="percent")

map_object_ss <- map_rangetype(map_object_tRNA, min_loop_width = 4,
                               type="ss", 
                               ss="/home/sigsk47/Documents/Genomes/fly/tRNA/Finished_tRNA_Drosophila_loop_anno.ss")       
map_object_ss <-  map_object_ss[!unlist(lapply(map_object_ss, function(x){x[[2]][1,1] == "no_hits"}))]

pac_tRNA <- tRNA_class(pac_tr, 
                       map=map_object_ss, 
                       terminal=5)

map_object <-  map_object[!unlist(lapply(map_object, function(x){x[[2]][1,1] == "no_hits"}))]
pac_tRNA <- tRNA_class(pac_tr, 
                       map=map_object, 
                       terminal=5)

pac_tRNA@Anno$type <- paste0(pac_tRNA@Anno$decoder, sep="-", pac_tRNA@Anno$acceptor)
pac_tRNA@Anno$type2 <- paste0(pac_tRNA@Anno$type, sep="-", pac_tRNA@Anno$class)


#  2.6.2 Mitochondrial vs nuclear tsRNA analysis                        ----------------


pac_tRNA@Anno$genome<-c("nuc")
pac_comb<-pac_tr@Anno[match(rownames(pac_tRNA@Anno), rownames(pac_tr@Anno)),]
pac_tRNA@Anno[grepl("chrM", pac_comb$mis0_genome),]$genome<-"mito"

df<-data.frame(pac_tRNA@summary$cpmMeans_DI, 
               type=paste0(pac_tRNA@Anno$type, sep="-", pac_tRNA@Anno$genome), 
               gen=pac_tRNA@Anno$genome, 
               type=pac_tRNA@Anno$class,
               genome=paste0(pac_tRNA@Anno$type, sep="-", pac_tRNA@Anno$class, sep="-", pac_tRNA@Anno$genome))
df_l<-tidyr::gather(df, key="sample", value="cpm", -c(type, gen, type.1, genome))
df_s<-aggregate(df_l$cpm, by=list(df_l$sample), FUN=mean)
#for heat map
df_hm<-df[,c(1:5)]
df_hm<-aggregate(df_hm[,1:4], by=list(df_hm$type), FUN=mean)
rownames(df_hm)<-df_hm$Group.1


out<-pheatmap::pheatmap(df_hm[,c(2,4,3,5)],
                        scale="row", 
                        cluster_col=FALSE, 
                        show_rownames = TRUE, 
                        cutree_rows = 4,
                        color = colorRampPalette(c("#0417D6", "#F9F4D1", "#D04116"))(50))

clstrs<-sort(cutree(out$tree_row, k=4))
clstrs<-as.data.frame(clstrs)


#  2.6.3 Negative Binomial Model using functions from section 2.5.2    ----------------


pac_tRNA@Anno$type3<-paste0(pac_tRNA@Anno$type, sep="-", pac_tRNA@Anno$genome)
df<-data.frame(pac_tRNA@norm$cpm, tr=pac_tRNA@Anno$type3)
df_l<-tidyr::gather(df, key="sample", value="cpm", -tr)

df_l$Diet<-sub("_.*", "\\1", df_l$sample)
df_l$DI <- sub('^([^_]+_[^_]+).*', '\\1', df_l$sample)

nb_res_list <- nb_func(df_l)
nb_res_list_NAC <- nb_func_NAC(df_l)

#  2.7 Data Visulaization Figure 3                                      ----------------   

Library(patchwork)


#  2.8 Color code                                                       ---------------- 


col45<-c("#d56d64",
         "#e19355",
         "#f4dd91",
         "#4c9051",
         "#b1cecf",
         "#a790cb",
         "#6282be")
scales::show_col(col45)

col_pie<-c("#a790cb", #lncrna
           "#f4dd91", #mirna
           "#b1cecf", #mt trna
           "#e19355", #pc
           "#d56d64", #piRNA
           "#4c9051", #other
           "#6282be") #tRNA)
scales::show_col(col_pie)


c45_10<-c("#d56d64",
          "#e19355",
          "#f4dd91",
          "#4c9051",
          "#b1cecf",
          "#a790cb",
          "#6282be",
          "#b77356",
          "#9c8c7c",
          "#8eaaa2")
scales::show_col(c45_10)


#  Fig 3.B


pie_h<-PAC_pie(pac_h, 
               anno_target = list("Biotypes_mis0",  c("protein_coding" , "piRNA","Mt_tRNA", "miRNA", "lncRNA", "tRNA")),
               summary="all", colors = rev(col_pie), angle=5)
pie_h

pie_d<-PAC_pie(pac_d, 
               anno_target = list("Biotypes_mis0", c("piRNA","protein_coding" , "Mt_tRNA", "miRNA", "lncRNA", "other", "tRNA")),
               summary="all", colors=rev(col_pie), angle=-100)

pie_d

pie_h + pie_d


#  Fig 3.C-sizedist


psd<-PAC_sizedist(pac_d, 
                  norm="cpm",
                  anno_target = list("Biotypes_mis0", c("lncRNA",  "miRNA", "Mt_tRNA", "other", "protein_coding", "piRNA", "tRNA")), 
                  summary_target = list("cpmMeans_DI"), 
                  colors = c45_10)

(psd$Histograms$D4_C+psd$Histograms$D4_NAC+psd$Histograms$D6_C+psd$Histograms$D6_NAC)+plot_layout(guides = 'collect')


#  Fig 3.D - heidi


df_plot<-data.frame(d4vsd6=pac_d@summary$Log2FC_DI$D4_C_vs_D6_C, 
                    d4vsnac=pac_d@summary$Log2FC_DI_rev$D4_NAC_vs_D4_C,
                    bio=pac_d@Anno$Biotypes_mis0)

#Organize biotypes
df_plot$bio<-factor(df_plot$bio, levels=c("miRNA","Mt_tRNA","tRNA", "other"))
df_plot$bio[is.na(df_plot$bio)]<-"other"
df_plot<-df_plot[order(df_plot$bio, decreasing=T),]

#plot graph
logp<-ggplot(df_plot, aes(x=d4vsnac, y=d4vsd6, color=bio))+
  geom_jitter(size=7, alpha=0.8)+
  theme_classic()+
  geom_hline(yintercept = 0, alpha=0.4)+
  geom_vline(xintercept= 0, alpha=0.4)+
  scale_color_manual(values = c( "#e19355",
                                 "#f4dd91",
                                 "#6282be", "#BFBFBF"))+
  ylab("logFC 30 g/L / 300 g/L ")+
  xlab("logFC C / NAC")+
  theme(legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size=18),
        axis.title = element_text(size=20))+
  annotate(geom="text", x= -1, y= 1.7, label="high 30g/L, low NAC", size=5,color="black")+
  annotate(geom="text", x= 2, y= -1.5, label="low 30g/L, high NAC", size=5,color="black")

logp


#  Fig 3.E


pac_miRNA_plot$plots$Expression_Anno_1$Grand_means + pac_miRNA_plot$plots$Log2FC_Anno_1


#  Fig 3.F - miRNA


seqs<-rownames(pac_miRNA@Anno[pac_miRNA@Anno$mir %in% "dme-mir-10",])
df_m10<-data.frame(pac_miRNA@norm$cpm[rownames(pac_miRNA@norm$cpm) %in% seqs,])
df_m10<-as.data.frame(colMeans(df_m10))
colnames(df_m10)<-"cpm"
df_m10$sample<-rownames(df_m10)
#df_m10_l<-gather(df_m10, key="sample", value="cpm")
df_m10$DI<-sub('^([^_]+_[^_]+).*', '\\1', df_m10$sample)
df_m10$DI<-factor(df_m10$DI, levels=c("D4_c", "D6_c", "D4_NAc", "D6_NAc"))

mir10plot<-ggplot(df_m10, aes(x=DI, y=cpm, fill=DI))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = c("#e67823","#ca9e7d","#d98f59","#cab09e"))

mir10plot


#  Fig 3.G - covplot mir-10


cp_mir10<-PAC_covplot(pac_miRNA, map=map_object, summary_target = list("cpmMeans_DI"),
                      map_target="dme-mir-10") 
cp_mir10$`dme-mir-10`
  
  
#  2.8 Data Visualization Figure 4                                      ---------------- 


#  Fig 4 B - piecharts mito vs nuc


df_nuc<-df_l[df_l$gen %in% "nuc",]
df2<-aggregate(df_nuc$cpm, by=list(df_nuc$type.1), FUN=mean)
df2$perc<-df2$x/sum(df2$x)*100

pie_nuc<-ggplot(df2, aes(x="", y=perc, fill=Group.1))+
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("#84b67d", "#c2d7a3", "#5167a4", "#ffd96b"))

df_mit<-df_l[df_l$gen %in% "mito",]
df22<-aggregate(df_mit$cpm, by=list(df_mit$type.1), FUN=mean)
df22$perc<-df22$x/sum(df22$x)*100

pie_mit<-ggplot(df22, aes(x="", y=perc, fill=Group.1))+
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("#84b67d", "#c2d7a3", "#5167a4", "#ffd96b", "#f3ebca"))

pie_nuc / pie_mit


#  Fig 4 C - heatmap


df_trna_hm<-data.frame(pac_trna$summary$cpmMeans_DI, type=pac_trna$Anno$type2)
df_trna_hm<-aggregate(.~type, data=df_trna_hm, mean)
rownames(df_trna_hm)<-df_trna_hm$type
df_trna_hm<-t(df_trna_hm[,2:5])
pheatmap(df_trna_hm,
         cluster_rows = F, cluster_cols = T,
         legend = T,
         show_rownames = T,
         gaps_row = c(1,2,3),
         scale="column",
         angle_col=c("45"),
         main="individual transcripts (all tRNA, n=229)",
         fontsize=13,
         color = colorRampPalette(c("#497279", "white", "#A44A3F"))(50))


#  Fig 4 D - coveragegraphs


pac_tRNA@summary$cpmMeans_DI<-pac_tRNA@summary$cpmMeans_DI[,c(3,4,2,1)]
cp_mit_tRNA<-PAC_covplot(pac_tRNA, map=map_object_tRNA, summary_target = list("cpmMeans_DI"),
                         style = "solid",
                         colors = c(  "#828c9e", "#aeb8c7", "#b8b8b8", "#f3dc8f"))
cp_mit_tRNA$`mt-tRNA-Gly-TCC` / cp_mit_tRNA$`mt-tRNA-Ile-GAT` / cp_mit_tRNA$`mt-tRNA-Met-CAT`


#  2.9 Data Visualization Supplementary Figures                         ---------------- 


#  Sup. Fig. S3 A - jitter diet 


pac_d@Anno$gen<-c("nuc")
pac_d@Anno[grepl("chrM", pac_d@Anno$mis0_genome),]$gen<-"mito"
pac_d@Anno$genbio<-paste0(pac_d@Anno$Biotypes_mis0, sep=" ", pac_d@Anno$gen)

jp_S3A<-PAC_jitter(pac_d, summary_target = list("Log2FC_DI"), anno_target = list("genbio", 
  c("lncRNA nuc", "miRNA nuc", "tRNA nuc", "Mt_tRNA mito", "protein_coding nuc", 
    "protein_coding mito", "piRNA nuc", "piRNA mito", "other nuc")),
     colors = c45_10)

jp_S3A$D4_C_vs_D6_C


#  Sup. Fig. S3 B - jitter diet vs nac 


jp_S3A$D4_C_vs_D4_NAC


#  Sup. Fig. S3 C - all miRNA 


pac_miRNA_plot_all$plots$Expression_Anno_1$Grand_means + pac_miRNA_plot_all$plots$Log2FC_Anno_1


#  Sup. Fig. S4 A -


df_fin<-aggregate(df_l$cpm, by=list(df_l$sample, df_l$type.1, df_l$gen), FUN=mean)
df_fin_all<-aggregate(df_l$cpm, by=list(df_l$sample, df_l$type.1), FUN=mean)
df_fin_all$Group.1<-factor(df_fin_all$Group.1, levels = c("D4_C", "D6_C", "D4_NAC", "D6_NAC"))

all_cpm<-ggplot(df_fin_all, aes(x=Group.1, y=x, fill=Group.2))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=c("#84b67d", "#c2d7a3", "#5167a4", "#ffd96b", "#f3ebca"))


#  Sup. Fig. S4 B -


clstrs$type<-rownames(clstrs)
df_tot<-merge(df_l, clstrs, by.x= "type", by.y = "type", all.x=TRUE)
df_fin<-aggregate(df_tot$cpm, by=list(df_tot$sample, df_tot$type.1, df_tot$clstrs), FUN=mean)
df_fin$Group.1<-factor(df_fin$Group.1, levels = c("D4_C", "D6_C", "D4_NAC", "D6_NAC"))
df_fin$Group.3<-factor(df_fin$Group.3, levels = c(3,4,2,1))

ggplot(df_fin, aes(x=Group.1, y=x, fill=Group.2))+
  geom_bar(position="stack", stat="identity")+
  facet_grid(.~Group.3)+
  scale_fill_manual(values=c("#84b67d", "#c2d7a3", "#5167a4", "#ffd96b", "#f3ebca"))


#  Sup. Fig. S4 C -

ggplot(df_fin, aes(x=Group.1, y=x, fill=Group.2))+
  geom_bar(position="fill", stat="identity")+
  facet_grid(.~Group.3)+
  scale_fill_manual(values=c("#84b67d", "#c2d7a3", "#5167a4", "#ffd96b", "#f3ebca"))

