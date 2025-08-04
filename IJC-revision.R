# Revision 

#CD274 samples: 36/501

rm(list=ls())
myinf1 = "/FILEPATH/PubDat/Dataset/Firehose/CNV/processed/LUSC_CNV_Symbol.rda"

fra.thr = 0.1	## in more than 10% samples
amp.thr = log2(2.8/2)
del.thr = log(1.2/2)




load(myinf1)
data = mydata
dim(data)
xx = as.numeric(substr(colnames(data), 14, 15))
se = which(xx==1)
data = data[,se]
dim(data)
colnames(data) = substr(colnames(data), 1, 12)


xx = apply(data<del.thr, 1, sum)
se = which(xx>ncol(data)*fra.thr)
length(se)
dat2 = data[se,]
dat2=data

dat2 = t(dat2)
for(k in 1:ncol(dat2))
{
  dat2[,k] = ifelse(dat2[,k]<del.thr, 1, 0)
}
tmp = apply(dat2, 2, sum)
dat2 = dat2[, order(tmp, decreasing=T)]

apply(dat2, 2, sum)

# CDKN2A  FOXP1   PCM1 
# 140     55     55 

xx= cor(dat2)
for(k in 1:ncol(xx))
{
  xx[k,k] = 0
}
# apply(xx, 1, max)
# CDKN2A       PCM1        WRN 
# 0.04718381 0.89279783 0.89279783 




cd274 = data['CD274',]
cdkn2a =data['CDKN2A',]
plot(cdkn2a, cd274, main = "CD274 vs CDKN2A CNV",
     xlab = "CDKN2A CNV", ylab = "CD274 CNV")
cor.test(cdkn2a, cd274)  
#Pearson's product-moment correlation
# data:  cdkn2a and cd274
# t = 6.3714, df = 499, p-value = 4.263e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1912761 0.3533947
# sample estimates:
#   cor 
# 0.2742831 


cdkn2a_del = ifelse(data["CDKN2A", ] < del.thr, 1, 0)
cd274_del  = ifelse(data["CD274", ] < del.thr, 1, 0)

common = intersect(names(cdkn2a_del), names(cd274_del))
cdkn2a_del = cdkn2a_del[common]
cd274_del = cd274_del[common]
combo_table = table(CDKN2A = cdkn2a_del, CD274 = cd274_del)
print(combo_table)

total_samples = length(common)
both_del = combo_table["1", "1"]
cdkn2a_only = combo_table["1", "0"]
cd274_only = combo_table["0", "1"]

#Number of mutual coverage
cat("\nOut of", total_samples, "samples:\n")
cat("- CDKN2A deletion only:", cdkn2a_only, "(", round(100 * cdkn2a_only / total_samples, 2), "% )\n")
cat("- CD274 deletion only:", cd274_only, "(", round(100 * cd274_only / total_samples, 2), "% )\n")
cat("- Both deletions:", both_del, "(", round(100 * both_del / total_samples, 2), "% )\n")


####################
rm(list=ls())
myinf1 = "/FILEPATH/WorSpa/m1_cancer/Lung/IntGen/result/TCGA-LUSC_GenomicEvent_iRAS_SQC.txt"
myinf2 = "/FILEPATH/WorSpa/m1_cancer/Lung/IntGen/data/TCGA_LUSC_Freq_GenomicEvents_Gene.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("_up\\.ES", "", tmp)
colnames(data) = tmp
se = grep("uni.adj", colnames(data))
data = data[,se]
colnames(data) = gsub(".up.ES", "", colnames(data))
colnames(data) = gsub("uni.adj__", "", colnames(data))
cdkn2a_del_sig = data$CDKN2A__DEL #deletion signature
names(cdkn2a_del_sig) = rownames(data)
cd274 #expression


x <- cdkn2a_del_sig
y <- cd274

plot(x, y,
     xlab = "CDKN2A deletion signature",
     ylab = "CD274 expression",
     main = "CD274 vs CDKN2A-Del signature",
     pch = 16, col = "darkblue")

abline(lm(y ~ x), col = "red", lwd = 2)
cor.test(x, y)
#Pearson's product-moment correlation

# data:  x and y
# t = -0.31454, df = 499, p-value = 0.7532
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.10155705  0.07361439
# sample estimates:
#   cor 
# cdkn2a 

cd274_del  = ifelse(data["CD274", ] < del.thr, 1, 0)


group = rep(NA, length(cdkn2a_del_sig))
names(group) = names(cdkn2a_del_sig)

for (i in 1:length(cdkn2a_del_sig)) {
  if (cdkn2a_del[i] == 1 & cd274_del[i] == 0) {
    group[i] = "CDKN2A-del / CD274-wt"
  } else if (cdkn2a_del[i] == 0 & cd274_del[i] == 1) {
    group[i] = "CDKN2A-wt / CD274-del"
  } else if (cdkn2a_del[i] == 1 & cd274_del[i] == 1) {
    group[i] = "CDKN2A-del / CD274-del"
  } else {
    group[i] = "CDKN2A-wt / CD274-wt"
  }
}

dat = data.frame(
  group = factor(group),
  score = as.numeric(cdkn2a_del_sig)
)

dat$group = factor(dat$group, levels = c(
  "CDKN2A-wt / CD274-wt",
  "CDKN2A-del / CD274-wt",
  "CDKN2A-wt / CD274-del",
  "CDKN2A-del / CD274-del"
))


pal <- c("#4DBBD5FF", "#E64B35FF", "#00A087FF", "#3C5488FF")
names(pal) <- levels(dat$group)

my_comparisons <- list(
  c("CDKN2A-wt / CD274-wt", "CDKN2A-del / CD274-wt"),
  c("CDKN2A-wt / CD274-wt", "CDKN2A-wt / CD274-del"),
  c("CDKN2A-wt / CD274-wt", "CDKN2A-del / CD274-del"),
  c("CDKN2A-del / CD274-wt", "CDKN2A-del / CD274-del")
)

g <- ggplot(dat, aes(x = group, y = score, fill = group)) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.4) +
  scale_fill_manual(values = pal) + ylab("CDKN2A-DEL  Score") +
  labs(title = "CDKN2A and CD274 CNV Status") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(20, 22, 24,26)) +
  theme(axis.text.x=element_text(angle = 30,size=16),
            axis.text.y=element_text(size=16),
            axis.title = element_text(size=16),
            plot.title = element_text(size=16, face="bold"),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16), 
            panel.background = element_rect(fill='transparent'), 
            plot.background = element_rect(fill='transparent', color=NA), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.background = element_rect(fill='transparent'), 
            legend.box.background = element_rect(fill='transparent'),
            panel.border = element_rect(color = "black",
                                        fill = NA,
                                        size = 1),
           aspect.ratio = 1,legend.position = "none") 

print(g)





#### Whether cdkn2aDEL primary in immune cold (?) 
## we can use the immune marker genes to separate hot vs. cold samples
# Can compare signature score variaiton between the two groups which dataset to use
# rm(list=ls())
library(reshape)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(stringr)
library(ComplexHeatmap) 
library(ggsci)
library(circlize)
myinf1 = "/FILEPATH/PubDat/Dataset/Firehose/RNAseqv2_Hiseq/processed/LUSC_RNAseqv2_Tumor_Symbol.rda" #tcga-RNAseq
load(file= myinf1)
data = mydata	
xx = apply(data>0, 1, sum)
se = which(xx>=10)## >0 in at least 10 samples
data = data[se,]
dim(data)
data = log10(data+1)
data= as.data.frame((data))
gene_set<-read.table("/mount/ictr1/chenglab/xwang/RCC_data/Figure/immune_genesets.txt",sep="\t",header=T,row.names=1)
set_list<-c()

for (i in row.names(gene_set)){
  if (i %in% c("Translation","Proliferation") ){
    next
  }
  set_list<- c(set_list,str_trim(unlist(strsplit(gene_set[i,],","))))
} 

library(ComplexHeatmap) 
genes<-intersect(set_list,row.names(data))
immune_expression<-data[genes,]
mat<-as.matrix(immune_expression)
scaled_mat = t(scale(t(mat)))
scaled_mat=mat

gene_ss<-data.frame(marker=rep(0,length(genes)))
rownames(gene_ss)<-genes
for ( i in genes){
  for (j in row.names(gene_set)){
    gene_ss[i,"marker"]=ifelse(i %in% str_trim(unlist(strsplit(gene_set[j,],","))),j,0)
    if (gene_ss[i,"marker"] !=0){
      break
    }
  } 
}


set.seed(123) 
gene_ss$marker<-as.factor(gene_ss$marker)
col_fun = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
lgd=Legend(col_fun = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")))
pal<-pal_frontiers("default")(8)
names(pal)<-levels(gene_ss$marker)
row_ha = rowAnnotation( df=gene_ss,col=list(marker=pal),show_annotation_name=FALSE)

ht<-Heatmap(scaled_mat, right_annotation=row_ha,
            show_column_names = FALSE, cluster_rows = FALSE,name="expression",
            show_row_dend = FALSE,column_km=2, column_gap = unit(5,"mm"),
            show_row_names=FALSE,column_title = c("immune hot", "immune cold"),
            column_title_gp = gpar(fill = c("#E64B35", "#3C5488"),col="white"))
draw(ht)

a<-column_dend(ht)
b<-column_order(ht)
clust <- lapply(names(b), function(i){
  out <- data.frame(sampleID = colnames(scaled_mat[,b[[i]]]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)

table(clust)

#---
# Compare Immune cold whether CDKN2A-Del enriched
clust$cdkn2a_del <- cdkn2a_del[clust$sampleID]

table_cluster = table(clust$Cluster, clust$cdkn2a_del)
print(table_cluster)

clust$in_cluster2 <- clust$Cluster == "cluster2" #IS this immune cold??
table_2x2 <- table(clust$in_cluster2, clust$cdkn2a_del)
print(table_2x2)

fisher.test(table_2x2)
library(ggplot2)
ggplot(na.omit(clust), aes(x = Cluster, fill = as.factor(cdkn2a_del))) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("0" = "#4DBBD5FF", "1" = "#E64B35FF"),
                    name = "CDKN2A Deletion") +
  ylab("Proportion of Samples") +
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title = element_text(size=16),
        plot.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16), 
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill='transparent'), 
        legend.box.background = element_rect(fill='transparent'),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        aspect.ratio = 1) +
  ggtitle("CDKN2A Deletion Frequency Across Clusters")




#######
rm(list=ls())
myinf1 = "/FILEPATH/WorSpa/Lab/T35_MPM_ITHsig/result/TCGA_MESO_ORACLE_PRACME_scores.txt"
myinf2 = "/FILEPATH/PubDat/Dataset/Firehose/done/Thorsson_2018_TCGA_immunelandscape.csv"

score = as.data.frame(cdkn2a_del_sig)

imm= read.table(myinf2, sep=",", header=T, row.names=1, stringsAsFactors=F)
se = which(imm$TCGA.Study=="LUSC")
data = imm[se,]

comxx = intersect(row.names(data), row.names(score))
data = data[comxx,]
score = score[comxx,]
dim(data)		

data = data[, c(4:31, 36:54)]
cor.df=(cor(data, score, method="s", use="pair"))
cor.df=as.data.frame(cor.df)
write.table(cor.df,'/mount/amos1/home/ylin/lung-prelesion/output/cor_df.txt' , sep="\t", quote=F)


#JAK2 co-loss
rm(list=ls())
myinf1 = "/FILEPATH/PubDat/Dataset/Firehose/CNV/processed/LUSC_CNV_Symbol.rda"

fra.thr = 0.1	## in more than 10% samples
amp.thr = log2(2.8/2)
del.thr = log(1.2/2)


load(myinf1)
data = mydata
dim(data)
xx = as.numeric(substr(colnames(data), 14, 15))
se = which(xx==1)
data = data[,se]
dim(data)
colnames(data) = substr(colnames(data), 1, 12)


xx = apply(data<del.thr, 1, sum)
se = which(xx>ncol(data)*fra.thr)
length(se)
dat2 = data[se,]
dat2=data

dat2 = t(dat2)
for(k in 1:ncol(dat2))
{
  dat2[,k] = ifelse(dat2[,k]<del.thr, 1, 0)
}
tmp = apply(dat2, 2, sum)
dat2 = dat2[, order(tmp, decreasing=T)]

apply(dat2, 2, sum)

# CDKN2A  FOXP1   PCM1 
# 140     55     55 

# xx= cor(dat2)
# for(k in 1:ncol(xx))
# {
#   xx[k,k] = 0
# }
# apply(xx, 1, max)
# CDKN2A       PCM1        WRN 
# 0.04718381 0.89279783 0.89279783 


jak2 = data['JAK2',]
cdkn2a =data['CDKN2A',]
jak2_del  = ifelse(data["JAK2", ] < del.thr, 1, 0)
cdkn2a_del  = ifelse(data["CDKN2A", ] < del.thr, 1, 0)
plot(cdkn2a, jak2, main = "JAK2 vs CDKN2A CNV",
     xlab = "CDKN2A CNV", ylab = "JAK2 CNV")
cor.test(cdkn2a, jak2)  
#Pearson's product-moment correlation
# data:  cdkn2a and jak2
# t = 6.8368, df = 499, p-value = 2.373e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.2104489 0.3707542
# sample estimates:
#   cor 
# 0.2926565
myinf1 = "/FILEPATH/WorSpa/m1_cancer/Lung/IntGen/result/TCGA-LUSC_GenomicEvent_iRAS_SQC.txt"
myinf2 = "/FILEPATH/WorSpa/m1_cancer/Lung/IntGen/data/TCGA_LUSC_Freq_GenomicEvents_Gene.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("_up\\.ES", "", tmp)
colnames(data) = tmp
se = grep("uni.adj", colnames(data))
data = data[,se]
colnames(data) = gsub(".up.ES", "", colnames(data))
colnames(data) = gsub("uni.adj__", "", colnames(data))
cdkn2a_del_sig = data$CDKN2A__DEL #deletion signature
names(cdkn2a_del_sig) = rownames(data)

x <- cdkn2a_del_sig
y <- jak2

plot(x, y,
     pch = 16, col = "darkblue")

abline(lm(y ~ x), col = "red", lwd = 2)
cor.test(x, y) # not significant

group = rep(NA, length(cdkn2a_del_sig))
names(group) = names(cdkn2a_del_sig)

for (i in 1:length(cdkn2a_del_sig)) {
  if (cdkn2a_del[i] == 1 & jak2_del[i] == 0) {
    group[i] = "CDKN2A-del / jak2-wt"
  } else if (cdkn2a_del[i] == 0 & jak2_del[i] == 1) {
    group[i] = "CDKN2A-wt / jak2-del"
  } else if (cdkn2a_del[i] == 1 & jak2_del[i] == 1) {
    group[i] = "CDKN2A-del / jak2-del"
  } else {
    group[i] = "CDKN2A-wt / jak2-wt"
  }
}

dat = data.frame(
  group = factor(group),
  score = as.numeric(cdkn2a_del_sig)
)

dat$group = factor(dat$group, levels = c(
  "CDKN2A-wt / jak2-wt",
  "CDKN2A-del / jak2-wt",
  "CDKN2A-wt / jak2-del",
  "CDKN2A-del / jak2-del"
))


pal <- c("#4DBBD5FF", "#E64B35FF", "#00A087FF", "#3C5488FF")
names(pal) <- levels(dat$group)

my_comparisons <- list(
  c("CDKN2A-wt / jak2-wt", "CDKN2A-del / jak2-wt"),
  c("CDKN2A-wt / jak2-wt", "CDKN2A-wt / jak2-del"),
  c("CDKN2A-wt / jak2-wt", "CDKN2A-del / jak2-del"),
  c("CDKN2A-del / jak2-wt", "CDKN2A-del / jak2-del")
)

g <- ggplot(dat, aes(x = group, y = score, fill = group)) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.4) +
  scale_fill_manual(values = pal) + ylab("CDKN2A-DEL  Score") +
  labs(title = "CDKN2A and jak2 CNV Status") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(20, 22, 24,26)) +
  theme(axis.text.x=element_text(angle = 30,size=16),
        axis.text.y=element_text(size=16),
        axis.title = element_text(size=16),
        plot.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16), 
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill='transparent'), 
        legend.box.background = element_rect(fill='transparent'),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        aspect.ratio = 1,legend.position = "none") 

print(g)





#### Whether co-loss (cdkn2a and jak2) primary in immune cold
## we can use the immune marker genes to separate hot vs. cold samples
# Can compare signature score variaiton between the two groups which dataset to use
# rm(list=ls())
library(reshape)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(stringr)
library(ComplexHeatmap) 
library(ggsci)
library(circlize)
myinf1 = "/FILEPATH/PubDat/Dataset/Firehose/RNAseqv2_Hiseq/processed/LUSC_RNAseqv2_Tumor_Symbol.rda" #tcga-RNAseq
load(file= myinf1)
data = mydata	
xx = apply(data>0, 1, sum)
se = which(xx>=10)## >0 in at least 10 samples
data = data[se,]
dim(data)
data = log10(data+1)
data= as.data.frame((data))
gene_set<-read.table("/mount/ictr1/chenglab/xwang/RCC_data/Figure/immune_genesets.txt",sep="\t",header=T,row.names=1)
set_list<-c()

for (i in row.names(gene_set)){
  if (i %in% c("Translation","Proliferation") ){
    next
  }
  set_list<- c(set_list,str_trim(unlist(strsplit(gene_set[i,],","))))
} 

library(ComplexHeatmap) 
genes<-intersect(set_list,row.names(data))
immune_expression<-data[genes,]
mat<-as.matrix(immune_expression)
scaled_mat = t(scale(t(mat)))
scaled_mat=mat

gene_ss<-data.frame(marker=rep(0,length(genes)))
rownames(gene_ss)<-genes
for ( i in genes){
  for (j in row.names(gene_set)){
    gene_ss[i,"marker"]=ifelse(i %in% str_trim(unlist(strsplit(gene_set[j,],","))),j,0)
    if (gene_ss[i,"marker"] !=0){
      break
    }
  } 
}


set.seed(123) 
gene_ss$marker<-as.factor(gene_ss$marker)
col_fun = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
lgd=Legend(col_fun = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")))
pal<-pal_frontiers("default")(8)
names(pal)<-levels(gene_ss$marker)
row_ha = rowAnnotation( df=gene_ss,col=list(marker=pal),show_annotation_name=FALSE)

ht<-Heatmap(scaled_mat, right_annotation=row_ha,
            show_column_names = FALSE, cluster_rows = FALSE,name="expression",
            show_row_dend = FALSE,column_km=2, column_gap = unit(5,"mm"),
            show_row_names=FALSE,column_title = c("immune cold", "immune hot"),
            column_title_gp = gpar(fill = c("#3C5488","#E64B35"),col="white"))
draw(ht)

a<-column_dend(ht)
b<-column_order(ht)
clust <- lapply(names(b), function(i){
  out <- data.frame(sampleID = colnames(scaled_mat[,b[[i]]]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)

table(clust)

#---
# Compare Immune cold whether CDKN2A-Del enriched
clust$cdkn2a_del <- cdkn2a_del[clust$sampleID]

table_cluster = table(clust$Cluster, clust$cdkn2a_del)
print(table_cluster)

clust$in_cluster2 <- clust$Cluster == "cluster1" #IS this immune cold??
table_2x2 <- table(clust$in_cluster2, clust$cdkn2a_del)
print(table_2x2)

fisher.test(table_2x2)

del_status_df <- data.frame(
  sampleID = names(group),
  CNV_group = group,
  stringsAsFactors = FALSE
)

df_merge <- merge(clust, del_status_df, by = "sampleID")
table(df_merge$Cluster, df_merge$CNV_group)
prop.table(table(df_merge$Cluster, df_merge$CNV_group), margin = 1) * 100  # Row-wise percentage
library(ggplot2)
se= which(df_merge$CNV_group =="CDKN2A-wt / jak2-wt" ) #358
ggplot(df_merge[-se,], aes(x = Cluster, fill = CNV_group)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Percentage", x = "Immune Cluster", fill = "CNV Group") +
  scale_fill_manual(values = c(
    "CDKN2A-wt / jak2-wt" = "#999999",
    "CDKN2A-del / jak2-wt" = "#E64B35",
    "CDKN2A-wt / jak2-del" = "#00A087",
    "CDKN2A-del / jak2-del" = "#3C5488"
  )) +
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title = element_text(size=16),
        plot.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16), 
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill='transparent'), 
        legend.box.background = element_rect(fill='transparent'),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        aspect.ratio = 1) +
  ggtitle("CNV Group Distribution")
df_merge$coloss <- ifelse(df_merge$CNV_group == "CDKN2A-del / jak2-del", 1, 0)
table_test <- table(df_merge$Cluster, df_merge$coloss)
fisher.test(table_test)


#################
# [Ref3.2]
# we can use the immune marker genes to separate hot vs. cold samples
# Can compare signature score variaiton between the two groups which dataset to use
rm(list=ls())
library(reshape)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(stringr)
library(ComplexHeatmap) 
library(ggsci)
library(circlize)
myinf1= "/FILEPATH/PubDat/Cancer/Lung/PreCancer/Mascaux_GSE33479/Mascaux_GSE33479_Symbol.txt"
data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)

gene_set<-read.table("/mount/ictr1/chenglab/xwang/RCC_data/Figure/immune_genesets.txt",sep="\t",header=T,row.names=1)
set_list<-c()

for (i in row.names(gene_set)){
  if (i %in% c("Translation","Proliferation") ){
    next
  }
  set_list<- c(set_list,str_trim(unlist(strsplit(gene_set[i,],","))))
} 

library(ComplexHeatmap) 
genes<-intersect(set_list,row.names(data))
immune_expression<-data[genes,]
mat<-as.matrix(immune_expression)
scaled_mat = t(scale(t(mat)))
scaled_mat=mat

gene_ss<-data.frame(marker=rep(0,length(genes)))
rownames(gene_ss)<-genes
for ( i in genes){
  for (j in row.names(gene_set)){
    gene_ss[i,"marker"]=ifelse(i %in% str_trim(unlist(strsplit(gene_set[j,],","))),j,0)
    if (gene_ss[i,"marker"] !=0){
      break
    }
  } 
}


set.seed(123) 
gene_ss$marker<-as.factor(gene_ss$marker)
col_fun = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
lgd=Legend(col_fun = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")))
pal<-pal_frontiers("default")(8)
names(pal)<-levels(gene_ss$marker)
row_ha = rowAnnotation( df=gene_ss,col=list(marker=pal),show_annotation_name=FALSE)

ht<-Heatmap(scaled_mat, right_annotation=row_ha,
            show_column_names = FALSE, cluster_rows = FALSE,name="expression",
            show_row_dend = FALSE,column_km=2, column_gap = unit(5,"mm"),
            show_row_names=FALSE,column_title = c("immune hot", "immune cold"),
            column_title_gp = gpar(fill = c("#E64B35", "#3C5488"),col="white"))
draw(ht)
# pdf("/mount/ictr1/chenglab/xwang/RCC_data/Figure/Fig3d.pdf",width=8.8,height=2.16)

a<-column_dend(ht)
b<-column_order(ht)
clust <- lapply(names(b), function(i){
  out <- data.frame(sampleID = colnames(scaled_mat[,b[[i]]]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)

# write.table(clust, file= "/mount/ictr1/chenglab/xwang/RCC_data/immune_cluster.txt", sep="\t", quote=F, row.names=FALSE)

library(reshape)
myinf1 = "/FILEPATH/WorSpa/m1_cancer/Lung/IntGen/result/Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf3 = "/FILEPATH/PubDat/Cancer/Lung/PreCancer/Mascaux_GSE33479/Clinical_info.txt"
info = read.table(myinf3, sep="\t", header=T, quote="", row.names=1)

data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("\\.up\\.ES", "", colnames(data))
se = grep("uni.adj__", colnames(data))
data = data[,se]
colnames(data) = gsub("uni.adj__", "", colnames(data))

#-----
str_sub(colnames(data),start=-5, end=-4) <- "-" 
df= read.csv("/mount/amos1/home/ylin/lung-prelesion/output/amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")
se1=which(colnames(data)%in%df$rep)
se2=which(!colnames(data)%in%df$entrezgene_accession)
data = data[,c(se1,se2)]
Sig = data

se= which(clust$Cluster=="cluster1")
c1 = clust$sampleID[which(clust$Cluster=="cluster1")]#88hot vs 34cold
se= which(rownames(Sig)%in%c1) #cold

res=NULL
for (i in 1:ncol(Sig)){
  xx=t.test(Sig[se,i],Sig[-se,i])
  tstat=mean(Sig[se,i])-mean(Sig[-se,i])
  xx=xx$p.value
  tmp = cbind(tstat,xx)
  res = rbind(res,tmp)
}
rownames(res) = colnames(Sig)
colnames(res) = c("Difference","pval")
res=as.data.frame(res)
length(which(res$pval <0.05))# two-sided: 37 ;greater: 32

#------------------------------------------
res=NULL

for (i in 1:ncol(Sig)){
  xx=t.test(Sig[se,i],Sig[-se,i])
  xx=xx$p.value
  res = rbind(res,xx)
}
rownames(res) = colnames(Sig)
length(which(res<0.05))
# [1] 37
res
# str_sub(colnames(Sig),start=-5, end=-4) <- "-" 
gene = "CDKN2A-DEL"
i=which(colnames(Sig)==gene)
# i=which(colnames(data)=="NF1-MUT")
score = Sig[,i]
xx = (score- mean(score)) / sd(score)
xx=as.data.frame(xx)
rownames(xx)=rownames(data)
se= which(rownames(Sig)%in%c1)
xx$group = "-"
xx$group[-se] ="Immune Hot"
xx$group[se]="Immune Cold"
dat=xx
mean(xx[se,1])
mean(xx[-se,1])
colnames(dat) = c("score","group")
ttest= t.test(xx$xx[se],xx$xx[-se])
ttest$p.value
library(ggpubr)
library(rstatix)
my_comparisons <- list( c("Immune Hot" , "Immune Cold"))
stat.test <- xx %>% 
  wilcox_test(xx ~ group, comparisons = my_comparisons) %>%
  add_xy_position(x = "group") %>%
  mutate(myformatted.p = paste0("p = ", ifelse(p < 1e-5, format(p, scientific = T), signif(p, digits = 1))))
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
pal=c( "#3C5488","#E64B35")
# levels(dat$group)= c(paste0(mut," (n=",table(dat$group)[1], ")"), paste0("WO (n=",table(dat$group)[2],")"))
levels(dat$group)= c(paste0("Immune \n Cold"), paste0("Immune \n Hot"))

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA,width=0.5/length(unique(dat$group)))
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Signature score")
g=g+labs(title = gene)
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 4/3)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g+stat_compare_means(comparisons = my_comparisons,label.y = c(3)) # Add pairwise comparisons p-value
print(g)

dev.off()


#################################################################################
#CD274
rm(list=ls())
#---------
myinf1 = "/FILEPATH/WorSpa/m1_cancer/Lung/IntGen/result/Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf2 = "/FILEPATH/PubDat/Cancer/Lung/PreCancer/Mascaux_GSE33479/Mascaux_GSE33479_Symbol.txt"
myinf3 = "/FILEPATH/PubDat/Cancer/Lung/PreCancer/Mascaux_GSE33479/Clinical_info.txt"
info = read.table(myinf3, sep="\t", header=T, quote="", row.names=1)
data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("\\.up\\.ES", "", colnames(data))
se = grep("uni.adj__", colnames(data))
data = data[,se]
colnames(data) = gsub("uni.adj__", "", colnames(data))
Sig = data

data = read.table(myinf2, sep="\t", header=T, quote="", row.names=1)
Imm = t(data)

exp = Imm[,which(colnames(Imm)=="CD274")]
#---------------
mycat = table(info$histology)
mycat = mycat[c(6, 2, 3, 4, 5, 7, 1, 8)]

res = matrix(0, ncol(Sig), length(mycat))
row.names(res) = colnames(Sig)
colnames(res) = names(mycat)
res = as.data.frame(res)

res.p = matrix(0, ncol(Sig), length(mycat))
row.names(res.p) = colnames(Sig)
colnames(res.p) = names(mycat)
res.p = as.data.frame(res.p)
calculate_correlation_p_value <- function(col, vector) {
  cor_test_result <- cor.test(col, vector, method = "s")
  return(cor_test_result$p.value)
}

exp= t(exp)
for(k in 1:ncol(res))
{
  se = which(info$histology==names(mycat)[k]) #stage
  scc = cor(Sig[se,], exp[se], method="s")
  res[,k] = scc
  pval<- sapply(Sig[se,], calculate_correlation_p_value, vector = exp[se])
  res.p[,k]=pval
}


res = round(res, 3)

xx = t(res)
se = which(colnames(xx)%in%c("CDKN2A__DEL"))
xx = xx[,se]
cor= as.data.frame(xx)
cor$stage = row.names(cor)
melted_cor <- melt(cor, id.vars = "stage")
names(melted_cor) =c("stage","sig","cor")

xx = t(round(res.p,3))
se = which(colnames(xx)%in%c("CDKN2A__DEL" ))
xx = xx[,se]
p=as.data.frame(xx)

p$stage = row.names(p)
melted_p <- melt(p, id.vars = "stage")
names(melted_p) =c("stage","sig","pval")


df <- merge(melted_cor, melted_p,  by = c("stage","sig"), all = FALSE)
custom_order <- c("normal","hyperplasia","metaplasia", "mild dysplasia","moderate dysplasia", "severe dysplasia", "carcinoma in situ","squamous cell carcinoma")

#Remove certain stages
se=which(df$stage %in% c("normal","moderate dysplasia","severe dysplasia"))
df= df[-se,]
# se = which(custom_order %in%c("normal","hyperplasia","moderate dysplasia","severe dysplasia") )
# custom_order = custom_order[-se]

# Order the "Category" column with the custom order
df$stage <- factor(df$stage, levels = custom_order)
df$Signifcant = "-"

df$Signifcant[which(df$cor<0)] = "Negative"
df$Signifcant[which(df$cor>0)] = "Positive"
df$Signifcant[which(df$pval>0.05)] = "Not significant"
df$Signifcant[which(df$cor>0.3)] = "Positive"
df$Signifcant[which(df$cor< (-0.3))] = "Negative"
table(df$Signifcant)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

ggplot(data = df, aes(x=sig, y = stage)) +
  geom_point(aes( color = Signifcant,size = abs(cor)),shape=16) +
  scale_size_continuous(range = c(0.1, 13), limits=c(0.1,1), name = 'Correlation')+
  scale_color_manual(values=alpha(c('darkblue',"white",'darkred')))+
  theme(axis.text.x=element_text(angle = 90,size=16),
        axis.text.y=element_text(size=16),
        axis.title = element_text(size=16),
        plot.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16), panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) #transparent legend pan)


dev.off()



cdkn2a_sig = Sig[, "CDKN2A__DEL"]
median_val = median(cdkn2a_sig, na.rm = TRUE)

high_samples = rownames(Sig)[which(cdkn2a_sig > median_val)]
low_samples  = rownames(Sig)[which(cdkn2a_sig <= median_val)]
run_correlation_by_group <- function(sample_ids, Sig, exp, info, mycat) {
  res = matrix(0, ncol(Sig), length(mycat))
  row.names(res) = colnames(Sig)
  colnames(res) = names(mycat)
  res = as.data.frame(res)
  
  res.p = matrix(0, ncol(Sig), length(mycat))
  row.names(res.p) = colnames(Sig)
  colnames(res.p) = names(mycat)
  res.p = as.data.frame(res.p)
  
  exp_group = exp[sample_ids]
  Sig_group = Sig[sample_ids, ]
  info_group = info[sample_ids, ]
  
  for(k in 1:ncol(res)) {
    se = which(info_group$histology == names(mycat)[k])
    if (length(se) > 3) {
      scc = cor(Sig_group[se,], exp_group[se], method="s")
      res[,k] = scc
      pval <- sapply(Sig_group[se,], calculate_correlation_p_value, vector = exp_group[se])
      res.p[,k] = pval
    } else {
      res[,k] = NA
      res.p[,k] = NA
    }
  }
  
  list(cor = round(t(res), 3), p = round(t(res.p), 3))
}

res_high = run_correlation_by_group(high_samples, Sig, exp, info, mycat)
res_low  = run_correlation_by_group(low_samples, Sig, exp, info, mycat)

get_long_df <- function(res, sig_name) {
  cor = as.data.frame(res$cor)
  cor$stage = rownames(cor)
  melted_cor = melt(cor, id.vars = "stage")
  colnames(melted_cor) = c("stage", "sig", "cor")
  
  p = as.data.frame(res$p)
  p$stage = rownames(p)
  melted_p = melt(p, id.vars = "stage")
  colnames(melted_p) = c("stage", "sig", "pval")
  
  df = merge(melted_cor, melted_p, by = c("stage", "sig"))
  df = df[df$sig == sig_name, ]
  df$group = ifelse(sig_name == "CDKN2A__DEL", "CDKN2A__DEL", sig_name)
  df
}

df_high = get_long_df(res_high, "CDKN2A__DEL")
df_low  = get_long_df(res_low, "CDKN2A__DEL")

df_high$SigGroup = "High CDKN2A-del score"
df_low$SigGroup  = "Low CDKN2A-del score"

df_all = rbind(df_high, df_low)

df_all$Signifcant = "Not significant"
# df_all$Signifcant[df_all$pval < 0.05 & df_all$cor > 0.3] = "Positive"
# df_all$Signifcant[df_all$pval < 0.05 & df_all$cor < -0.3] = "Negative"

df_all$Signifcant[df_all$pval < 0.05 ] = "Positive"
df_all$Signifcant[df_all$pval < 0.05 ] = "Negative"


ggplot(df_all, aes(x = sig, y = stage)) +
  geom_point(aes(color = Signifcant, size = abs(cor)), shape = 16) +
  facet_wrap(~SigGroup) +
  scale_size_continuous(range = c(0.1, 13), limits=c(0.1,1), name = 'Correlation') +
  scale_color_manual(values = alpha(c('darkblue', "white", 'darkred'))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "Correlation with CD274 Expression", y = "Stage", x = "")

#------cd8t
library(reshape)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

myinf1 = "/FILEPATH/WorSpa/m1_cancer/Lung/IntGen/result/Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf2 = "/FILEPATH/WorSpa/m1_cancer/Lung/IntGen/result/Mascaux_GSE33479_ImmGen_representative_TILs.txt"
myinf3 = "/FILEPATH/PubDat/Cancer/Lung/PreCancer/Mascaux_GSE33479/Clinical_info.txt"

info = read.table(myinf3, sep = "\t", header = TRUE, quote = "", row.names = 1)

data = read.table(myinf1, sep = "\t", header = TRUE, quote = "", row.names = 1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("\\.up\\.ES", "", colnames(data))
se = grep("uni.adj__", colnames(data))
data = data[, se]
colnames(data) = gsub("uni.adj__", "", colnames(data))
Sig = data

data = read.table(myinf2, sep = "\t", header = TRUE, quote = "", row.names = 1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("_up\\.ES", "", colnames(data))
Imm = data

common_samples = Reduce(intersect, list(rownames(Sig), rownames(Imm), rownames(info)))
Sig = Sig[common_samples, ]
Imm = Imm[common_samples, ]
info = info[common_samples, ]


mycat = table(info$histology)
mycat = mycat[c(6, 2, 3, 4, 5, 7, 1, 8)]

cdkn2a_scores = Sig[, "CDKN2A__DEL"]
threshold = median(cdkn2a_scores, na.rm = TRUE)
high_samples = rownames(Sig)[cdkn2a_scores > threshold]
low_samples  = rownames(Sig)[cdkn2a_scores <= threshold]

exp = as.data.frame(Imm[,which(colnames(Imm)=="Monocyte")])
rownames(exp)=rownames(Imm)
run_correlation_by_group <- function(sample_ids, Sig, exp, info, mycat) {
  res = matrix(0, ncol(Sig), length(mycat))
  row.names(res) = colnames(Sig)
  colnames(res) = names(mycat)
  res = as.data.frame(res)
  
  res.p = matrix(0, ncol(Sig), length(mycat))
  row.names(res.p) = colnames(Sig)
  colnames(res.p) = names(mycat)
  res.p = as.data.frame(res.p)
  
  exp_group = exp[sample_ids,1]
  Sig_group = Sig[sample_ids, ]
  info_group = info[sample_ids, ]
  
  for(k in 1:ncol(res)) {
    se = which(info_group$histology == names(mycat)[k])
    if (length(se) > 3) {
      scc = cor(Sig_group[se,], exp_group[se], method="s")
      res[,k] = scc
      pval <- sapply(Sig_group[se,], calculate_correlation_p_value, vector = exp_group[se])
      res.p[,k] = pval
    } else {
      res[,k] = NA
      res.p[,k] = NA
    }
  }
  
  list(cor = round(t(res), 3), p = round(t(res.p), 3))
}

res_high = run_correlation_by_group(high_samples, Sig, exp, info, mycat)
res_low  = run_correlation_by_group(low_samples, Sig, exp, info, mycat)

get_long_df <- function(res, sig_name) {
  cor = as.data.frame(res$cor)
  cor$stage = rownames(cor)
  melted_cor = melt(cor, id.vars = "stage")
  colnames(melted_cor) = c("stage", "sig", "cor")
  
  p = as.data.frame(res$p)
  p$stage = rownames(p)
  melted_p = melt(p, id.vars = "stage")
  colnames(melted_p) = c("stage", "sig", "pval")
  
  df = merge(melted_cor, melted_p, by = c("stage", "sig"))
  df = df[df$sig == sig_name, ]
  df$group = ifelse(sig_name == "CDKN2A__DEL", "CDKN2A__DEL", sig_name)
  df
}

df_high = get_long_df(res_high, "CDKN2A__DEL")
df_low  = get_long_df(res_low, "CDKN2A__DEL")

df_high$SigGroup = "High CDKN2A-del score"
df_low$SigGroup  = "Low CDKN2A-del score"

df_all = rbind(df_high, df_low)

df_all$Signifcant = "Not significant"
# df_all$Signifcant[df_all$pval < 0.05 & df_all$cor > 0.3] = "Positive"
# df_all$Signifcant[df_all$pval < 0.05 & df_all$cor < -0.3] = "Negative"

df_all$Signifcant[df_all$pval < 0.05 ] = "Positive"
df_all$Signifcant[df_all$pval < 0.05 ] = "Negative"



ggplot(df_all, aes(x = sig, y = stage)) +
  geom_point(aes(color = Signifcant, size = abs(cor)), shape = 16) +
  facet_wrap(~SigGroup) +
  scale_size_continuous(range = c(0.1, 13), limits = c(0.1, 1), name = 'Correlation') +
  scale_color_manual(values = c("Negative" = "darkblue", "Positive" = "darkred", "Not significant" = "grey")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "CDKN2A-del Signature Correlation with CD8T by Signature Group", y = "Stage", x = "")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.1] define weight profiles of chrosome bands
rm(list=ls())
myinf1 = "/FILEPATHPubDat/Dataset/Firehose/RNAseqv2_Hiseq/processed/LUSC_RNAseqv2_Tumor_Symbol.rda"
myinf2 = "/FILEPATHPubDat/Dataset/Firehose/Clinical/processed2/LUSC_Clincial_info.txt"
myinf3A = "/FILEPATHPubDat/Dataset/Firehose/done/TCGA_Firehose_CNV_Chrband_All.txt"
myinf3B = "/FILEPATHWorSpa/m1_cancer/Lung/IntGen/data/TCGA_LUSC_Freq_GeneCNV.txt"
myoutf1 = "/FILEPATH/lung-prelesion/output/TCGA-LUNG-GenomicEvent-tempProfile_4base.txt"

#----------------------------
data = read.table(myinf3A, sep="\t", header=T, row.names=1, quote="")
se = grep("LUSC__", row.names(data))
data = data[se,]
dim(data)
row.names(data) = gsub("LUSC__", "", row.names(data))
xx = as.numeric(substr(row.names(data), 14, 15))
se = which(xx==1)
data = data[se,]
dim(data)
row.names(data) = substr(row.names(data), 1, 12)
for(k in 1:ncol(data))
{
  data[,k] = ifelse(data[,k]==-1, 1, 0)
}
som1 = data


soma = read.table(myinf3B, sep="\t", header=T, row.names=1, quote="")
for(k in 1:ncol(soma))
{
  soma[,k] = ifelse(soma[,k]>0, 1, 0)
}
som2 = soma

##----------------
comxx = intersect(row.names(som1), row.names(som2))
som1 = som1[comxx,]
som2 = som2[comxx,]

tmp = data.frame(CDKN2A__DEL=som2[, "CDKN2A__DEL"], chr9_p21 = som1[, "chr9_p21"])
row.names(tmp) = row.names(som1)
soma = tmp
for(k in 1:ncol(soma))
{
  soma[,k] = ifelse(soma[,k]>0, 1, 0)
}
apply(soma, 2, sum)
sum(soma[,1]*soma[,2])		## 36

#---------------------------
load(file= myinf1)
data = mydata	## 262 samples
xx = apply(data>0, 1, sum)
se = which(xx>=10)## >0 in at least 10 samples
data = data[se,]
dim(data)
data = log10(data+1)

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
se = c("age_at_initial_pathologic_diagnosis", "gender", "stage_event.pathologic_stage", "number_pack_years_smoked")
info = info[,se]
colnames(info) = c("age", "gender", "stage", "smoker")
info$gender = as.character(info$gender)
info$age = as.numeric(info$age)
info$smoker = as.numeric(as.character(info$smoker))
xx = as.character(info$stage)
xx = gsub("stage ", "", xx)
se = grep("^iv", xx)
xx[se] = "IV"
se = grep("^iii", xx)
xx[se] = "III"
se = grep("^ii", xx)
xx[se] = "II"
se = grep("^i", xx)
xx[se] = "I"
info$stage = xx

comSam = intersect(colnames(data), row.names(info))
data = data[, comSam]
info = info[comSam,]
raw.info = info
raw.data = data
#~~~~~~~~~~~~~~~~~~~~
tmp = matrix(0, nrow(data), ncol(soma))
row.names(tmp) = row.names(data)
colnames(tmp) = colnames(soma)
p1 = p2 = p3 = p4  =  tmp
## uni.noj --> p1:   univariate using each genomic feature without adjust by clinical variables
## uni.adj --> p2:   univariate using each genomic feature adjusted by clinical variables
## mul.noj --> p3:   mulvariate using all genomic features without adjust by clinical variables
## mul.adj --> p4:   mulvariate using all genomic features adjusted by clinical variables


for(s in 1:ncol(soma))
{
  mut = row.names(soma)[soma[,s]>0]
  wt =  row.names(soma)[soma[,s]==0]	
  mut = mut[mut%in%colnames(data)]
  wt = wt[wt%in%colnames(data)]
  xx = c(mut, wt)
  info = raw.info[xx, ]
  data = raw.data[, xx]
  mut = ifelse(row.names(info)%in%wt, 0, 1)
  info = cbind(mut, info)
  
  pval1 = beta1 = rep(0, nrow(data))
  pval2 = beta2 = rep(0, nrow(data))
  
  for(k in 1:nrow(data))
  {
    cat("\r\r\r", s, "-->", k)
    mytf = as.numeric(data[k,])
    xx = cbind(mytf, info)
    mylog <- lm(mytf ~ mut, data = xx)
    mylog = summary(mylog)
    pval1[k] =  mylog[["coefficients"]][2, 4]
    beta1[k] =  mylog[["coefficients"]][2, 1]
    
    mylog <- lm(mytf ~ mut + age + gender + stage + smoker, data = xx, family = "binomial")
    mylog = summary(mylog)
    pval2[k] =  mylog[["coefficients"]][2, 4]
    beta2[k] =  mylog[["coefficients"]][2, 1]
  }
  xx = -log10(pval1)
  xx = ifelse(beta1>0, xx, -xx)
  p1[,s] = xx
  xx = -log10(pval2)
  xx = ifelse(beta2>0, xx, -xx)
  p2[,s] = xx		
}


#--------------------------
res = p1
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp1 = ifelse(tmp>0, tmp, 0)
  tmp1[tmp1>10] = 10
  res1[,k] = tmp1
  tmp2 = ifelse(tmp<0, -tmp, 0)
  tmp2[tmp2>10] = 10
  res2[,k] = tmp2
}
colnames(res1)= paste(colnames(res1), ".up", sep="")
colnames(res2)= paste(colnames(res2), ".dn", sep="")
myres = cbind(res1, res2)
row.names(myres) = row.names(res)
colnames(myres) = paste("uni.noj", colnames(myres), sep="__")
minv = min(myres)
maxv = max(myres)
myres = (myres-minv)/(maxv-minv)
prof1 = myres


res = p2
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp1 = ifelse(tmp>0, tmp, 0)
  tmp1[tmp1>10] = 10
  res1[,k] = tmp1
  tmp2 = ifelse(tmp<0, -tmp, 0)
  tmp2[tmp2>10] = 10
  res2[,k] = tmp2
}
colnames(res1)= paste(colnames(res1), ".up", sep="")
colnames(res2)= paste(colnames(res2), ".dn", sep="")
myres = cbind(res1, res2)
row.names(myres) = row.names(res)
colnames(myres) = paste("uni.adj", colnames(myres), sep="__")
minv = min(myres)
maxv = max(myres)
myres = (myres-minv)/(maxv-minv)
prof2 = myres


#------------------------------------------------------------
profile = cbind(prof1, prof2)
se1 = grep("\\.up", colnames(profile))
se2 = grep("\\.dn", colnames(profile))
profile = profile[, c(se1, se2)]
dim(profile)	## 20500

write.table(profile, myoutf1, sep="\t", quote=F)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1.2] Similarity between CDNA2A-DEL and chr9_p21__DEL
rm(list=ls())
myinf1 = "/FILEPATH/lung-prelesion/output/TCGA-LUNG-GenomicEvent-tempProfile_4base.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
nn = ncol(data)/2
tmp = data[, (1:nn)] - data[, (nn+1):(nn*2)]
se = grep("uni.adj__", colnames(tmp))
tmp = tmp[,se]

cor(tmp, method="s")

uni.adj__CDKN2A__DEL.up uni.adj__chr9_p21.up
uni.adj__CDKN2A__DEL.up               1.0000000            0.6334436
uni.adj__chr9_p21.up                  0.6334436            1.0000000

cor(tmp, method="p")




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [2] calculate score
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1]Teixeira_GSE108124
rm(list=ls())
myinf1 = "/FILEPATHPubDat/Cancer/Lung/PreCancer/Teixeira_GSE108124/Teixeira_GSE94611_Symbol_expr.txt"
myinf2 = "/FILEPATH/lung-prelesion/output/TCGA-LUNG-GenomicEvent-tempProfile_4base.txt"

myoutf = "/FILEPATH/lung-prelesion/output/Teixeira_GSE94611_tmoRev_score.txt"

mywt = read.table(myinf2, sep="\t", header=T, quote="", row.names=1)
source("~/WorSpa/system/myRprogram/base5.R")
data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1, stringsAsFactors=F)


##
reg = mywt
xx = base5(data, reg, perm=1000, myoutf, median.norm=T)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.2]LUSC TCGA
rm(list=ls())
myinf1 = "/FILEPATHPubDat/Dataset/Firehose/RNAseqv2_Hiseq/processed/LUSC_RNAseqv2_Tumor_Symbol.rda"
myinf2 = "/FILEPATH/lung-prelesion/output/TCGA-LUNG-GenomicEvent-tempProfile_4base.txt"

myoutf = "/FILEPATH/lung-prelesion/output/TCGA-LUSC_GenomicEvent_tmoRev_score.txt"

mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")

source("base5.R")
load(file= myinf1)
data = mydata	## 262 samples
tag = apply(data>0,1, sum)
data = data[tag>10,]
data = log10(data+1)

##
reg = mywt
xx = base5(data, reg, perm=1000, myoutf, median.norm=T)

EOF

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.3] Mascaux_GSE33479
rm(list=ls())
myinf1 = "/FILEPATHPubDat/Cancer/Lung/PreCancer/Mascaux_GSE33479/Mascaux_GSE33479_Symbol.txt"
myinf2 = "/FILEPATH/lung-prelesion/output/TCGA-LUNG-GenomicEvent-tempProfile_4base.txt"

myoutf = "/FILEPATH/lung-prelesion/output/Mascaux_GSE33479_GenomicEvent_tmoRev_score.txt"
mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")

source("~/WorSpa/system/myRprogram/base5.R")

##
reg = mywt
xx = base5(data, reg, perm=1000, myoutf, median.norm=T)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.4]Bueno_GSE157011
rm(list=ls())
myinf1 = "/FILEPATHPubDat/Cancer/Lung/Bueno_GSE157011/Bueno_GSE157011_Lung_Symbol_expr.txt"
myinf2 = "/FILEPATH/lung-prelesion/output/TCGA-LUNG-GenomicEvent-tempProfile_4base.txt"

myoutf = "/FILEPATH/lung-prelesion/output/Bueno_GSE157011_GenomicEvent_tmoRev_score.txt"

mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")

source("base5.R")
data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)

##
reg = mywt
xx = base5(data, reg, perm=1000, myoutf, median.norm=T)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [3] Statistical Analysis
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [3.1] associatation with survival
rm(list=ls())
myinf1 = "/FILEPATH/lung-prelesion/output/Bueno_GSE157011_GenomicEvent_tmoRev_score.txt"
myinf2 = "/FILEPATHPubDat/Cancer/Lung/Bueno_GSE157011/Clinical_info.txt"


info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = as.numeric(info$os_mo)
e.surv = info$os_event
info = cbind(t.surv, e.surv, info)
tag = !is.na(info$t.surv)
info = info[tag==1,]
dim(info)
info[1,]
se = grep("^CAD", info$Title)
info = info[-se,]
#info = info[-se,]

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp

cor(data, method="s")
##                      uni.adj__CDKN2A__DEL uni.adj__chr9_p21
uni.adj__CDKN2A__DEL            1.0000000         0.9393420
uni.adj__chr9_p21               0.9393420         1.0000000


##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(info)		## 

library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, survreg.pval1,   coxph.pval1, coxph.qval1, hr1)
res

se = grep("uni.adj__", res$name)
xx = res[se,]
xx = xx[order(xx$coxph.pval1), ]


name survreg.pval1 coxph.pval1 coxph.qval1     hr1
3 uni.adj__CDKN2A__DEL   0.006256455 0.007707736  0.01027698 1.00875
4    uni.adj__chr9_p21   0.012146594 0.015426335  0.01542633 1.00946



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.2] corr. with features
rm(list=ls())
myinf1 = "/FILEPATH/lung-prelesion/output/TCGA-LUSC_GenomicEvent_tmoRev_score.txt"
myinf2 = "/FILEPATHPubDat/Dataset/Firehose/done/Thorsson_2018_TCGA_immunelandscape.csv"

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp


uni.adj__CDKN2A__DEL uni.adj__chr9_p21
uni.adj__CDKN2A__DEL            1.0000000         0.9246930
uni.adj__chr9_p21               0.9246930         1.0000000



info = read.table(myinf2, sep=",", header=T, row.names=1, quote="")
se = which(info$TCGA.Study=="LUSC")
info = info[se,]

comxx = intersect(row.names(data), row.names(info))
data = data[comxx, ]
info = info[comxx,]

se = grep("uni.adj__", colnames(data))
data = data[,se]
colnames(data) = c("CDKN2A__DEL", "chr9p21__DEL")

#------------------
subset = c(4:31, 36:63)
info = info[, subset]

xx = cor(info, data, method="s", use="pair")
xx1 = xx[order(xx[,2]), ]
xx2 = xx[order(-xx[,2]), ]

xx1[1:20,]
xx2[1:20,]

CDKN2A__DEL chr9p21__DEL
Macrophage.Regulation                   -0.76263464   -0.8311365
Lymphocyte.Infiltration.Signature.Score -0.62839139   -0.7429141
Stromal.Fraction                        -0.64259040   -0.6832262
TCR.Richness                            -0.58760100   -0.6552021
TCR.Shannon                             -0.56629184   -0.6290427
Leukocyte.Fraction                      -0.53750508   -0.6205011
BCR.Richness                            -0.51741023   -0.5437145
BCR.Shannon                             -0.49964723   -0.5259194
TGF.beta.Response                       -0.59898590   -0.4743222
IFN.gamma.Response                      -0.40640064   -0.4435780

CDKN2A__DEL chr9p21__DEL
Fraction.Altered                  0.58883300   0.59949451
Proliferation                     0.52982765   0.50297982
Number.of.Segments                0.43569278   0.42308327
Homologous.Recombination.Defects  0.43192850   0.40014097


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.3] Predict progression
rm(list=ls())
myinf1 = "/FILEPATH/lung-prelesion/output/Teixeira_GSE94611_tmoRev_score.txt"
myinf2 = "/FILEPATHPubDat/Cancer/Lung/PreCancer/Teixeira_GSE108124/GSE94611_Clinical_info.txt"


data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("\\.up\\.ES", "", colnames(data))

info = read.table(myinf2, sep="\t", header=T, quote="", row.names=1)

se1 = grep("Progressive", info[,1])
se2 = grep("Regressive", info[,1])
sam1 = row.names(info)[se1]
sam2 = row.names(info)[se2]

##
se = which(row.names(data)%in%sam1)
dat1 = data[se,]
se = which(row.names(data)%in%sam2)
dat2 = data[se,]


library(ROCR)
cat = ifelse(row.names(data)%in%sam1, 1, 0)
myauc = rep(0, ncol(data))
names(myauc) = colnames(data)
for(k in 1:ncol(data))
{
  score = data[,k]
  pred <- prediction(score, cat)
  auc.perf <- performance(pred, measure = "auc")
  myauc[k] = as.numeric(auc.perf@y.values)
}
se = grep("uni.adj__", names(myauc))
sort(myauc[se])


uni.adj__CDKN2A__DEL    uni.adj__chr9_p21 
0.8161765            0.8198529 














#@@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#@@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#@@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
The code below are obsoleted code -- DO NOT USE
The code below are obsoleted code -- DO NOT USE


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1.1] Identify chr band deletion more sensitively
## if the bas(Segment_Mean)>0.5 --> amp/del
## a chr band is called amp/del, if it overlap with a segment by >50%
rm(list=ls())
myDir1 = "/FILEPATHPubDat/Dataset/Firehose/CNV/raw/"
myinf2 = "/FILEPATHPubDat/organisms/human/annotation/hg19_cytoBand_Less.txt"

myoutf1 = "/FILEPATH/lung-prelesion/output/TCGA_Firehose_CNV_Chrband_All_V2.txt"
amp.thr = 0.5
del.thr = -0.4
frac.thr = 0.5

myinf1 = dir(myDir1)
mytag = myinf1
for(k in 1:length(mytag))
{
  mytag[k] = unlist(strsplit(mytag[k], "\\.snp__"))[1]
}
myinf1 = paste(myDir1, myinf1, sep="")

for(k in 1:length(myinf1))
{
  cat("\r", k)
  tmp = read.table(myinf1[k], sep="\t",  header=T, stringsAsFactors=F)
  tmp$Sample = paste(mytag[k], tmp$Sample, sep="__")
  xx = as.character(tmp$Chromosome)
  xx[xx=="23"]="X"
  tmp$Chromosome = xx
  if(k==1)
  {
    data = tmp
  }else
  {
    data = rbind(data, tmp)
  }
  
}
mysam = unique(data$Sample)
dim(data)		## 3251138 6
length(mysam)
data$Chromosome = paste("chr", data$Chromosome, sep="")


chrban = read.table(myinf2, sep="\t", header=F,  check.names=F)
row.names(chrban) = paste(chrban[,1], chrban[,4], sep="_")
se = which(chrban$V1=="chrY")
chrban = chrban[-se,]

res = matrix(0, length(mysam), nrow(chrban))
row.names(res)  = mysam
colnames(res) = row.names(chrban)
res = as.data.frame(res)

data = data[data$Segment_Mean>amp.thr | data$Segment_Mean<del.thr, ]
dim(data)		# 1192724  6

for(k in 1:ncol(res))
{
  cat("\r", k)
  tmp.data = subset(data, data$Chromosome==chrban$V1[k])
  L2 = chrban$V2[k]
  R2 = chrban$V3[k]
  L1 = tmp.data$Start
  R1 = tmp.data$End
  L1L2 = L1-L2
  L1R2 = L1-R2
  R1L2 = R1-L2
  R1R2 = R1-R2
  
  ## case 1: overlap
  se = which(L1L2>0 & L1R2<0 & R1L2>0 & R1R2>0)
  if(length(se)>0)
  {
    tmp = tmp.data[se, , drop=F]
    dist = abs(L1R2[se])
    se1 = which((dist/(R2-L2))>frac.thr)
    if(length(se1)>0)
    {
      tmp = tmp[se1, , drop=F]
      xx = tmp$Sample
      res[xx, k] = sign(tmp$Segment_Mean)
    }
  }
  
  # case 2: overlap	
  se = which(L1L2<0 & L1R2<0 & R1L2>0 & R1R2<0)
  if(length(se)>0)
  {
    tmp = tmp.data[se, , drop=F]
    dist = abs(R1L2[se])
    se1 = which((dist/(R2-L2))>frac.thr)
    if(length(se1)>0)
    {
      tmp = tmp[se1, , drop=F]
      xx = tmp$Sample
      res[xx, k] = sign(tmp$Segment_Mean)
    }
  }
  
  ## case 3: segment cover chr band	
  se = which(L1L2<0 & L1R2<0 & R1L2>0 & R1R2>0)
  if(length(se)>0)
  {
    tmp = tmp.data[se, , drop=F]
    xx = tmp$Sample
    res[xx, k] = sign(tmp$Segment_Mean)
  }
  
  ## case 4: chr band covers segment 	
  se = which(L1L2>0 & L1R2<0 & R1L2>0 & R1R2<0)
  if(length(se)>0)
  {
    tmp = tmp.data[se, , drop=F]
    dist = R1[se]-L1[se]
    se1 = which((dist/(R2-L2))>frac.thr)
    if(length(se1)>0)
    {
      tmp = tmp[se1, , drop=F]
      xx = tmp$Sample
      res[xx, k] = sign(tmp$Segment_Mean)
    }
  }	
}

write.table(res, myoutf1, sep="\t", quote=F)
x1 = apply(res==1, 2, sum)
x2 = apply(res==-1, 2, sum)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1.2] Chromsome band amp/del
rm(list=ls())
myinf1 = "/FILEPATH/lung-prelesion/output/TCGA_Firehose_CNV_Chrband_All_V2.txt"
## myinf2 = "/FILEPATHPubDat/Dataset/COSMIC/raw/cancer_gene_census.tsv"

myoutf1 = "/FILEPATH/lung-prelesion/output/TCGA_LUSC_Freq_ChrBand_V2.txt"

fra.thr = 0.1	## in more than 10% samples


data = read.table(myinf1, sep="\t", header=T, row.names=1)
se = grep("LUSC__", row.names(data))
data = data[se,]
dim(data)
row.names(data) = gsub("LUSC__", "", row.names(data))
xx = as.numeric(substr(row.names(data), 14, 15))
se = which(xx==1)
data = data[se,]
dim(data)
row.names(data) = substr(row.names(data), 1, 12)

xx = apply(data==1, 2, sum)
se = which(xx>nrow(data)*fra.thr)
length(se)
dat1 = data[,se]

xx = apply(data==-1, 2, sum)
se = which(xx>nrow(data)*fra.thr)
length(se)
dat2 = data[,se]

##
apply(dat1==1, 2, sum)
# chr2_p15 chr3_q11 chr3_q12 chr3_q13 chr3_q21 chr3_q22 chr3_q23 chr3_q24 
# 55       62       78       83      101      110      140      170 
# chr3_q25 chr3_q26 chr3_q27 chr3_q28 chr3_q29 chr5_p15 chr5_p14 chr5_p13 
# 197      230      271      255      251      132      145      140 
# chr5_p12 chr8_q23 chr8_q24 
# 124       55       58 
# 
# chr2_p15 chr3_q11 chr3_q12 chr3_q13 chr3_q21 chr3_q22 chr3_q23 chr3_q24 
# 55       69       85       85      102      110      141      171 
# chr3_q25 chr3_q26 chr3_q27 chr3_q28 chr3_q29 chr5_p15 chr5_p14 chr5_p13 
# 197      230      271      255      251      135      148      142 
# chr5_p12 chr8_q23 chr8_q24 
# 128       55       58 


se = c("chr2_p15", "chr3_q27", "chr5_p14", "chr8_q24")      
dat1 = dat1[,se]
colnames(dat1) = paste(colnames(dat1), "__AMP", sep="")

apply(dat2==-1, 2, sum)
# chr13_q12 chr13_q13 chr13_q14 chr13_q21 chr13_q22 chr13_q31 chr13_q32 chr13_q33 
# 62        64        69        64        66        63        58        55 
# chr13_q34 chr17_p13 chr17_p12 chr17_p11  chr3_p26  chr3_p25  chr3_p24  chr3_p23 
# 55        68        70        61        83        86        85        88 
# chr3_p22  chr3_p21  chr3_p14  chr3_p13  chr3_p12  chr3_p11  chr4_p16  chr4_p15 
# 87        93       101       115       110        91        72        69 
# chr4_p14  chr4_p13  chr4_q28  chr4_q31  chr4_q32  chr4_q33  chr4_q34  chr4_q35 
# 62        57        54        59        66        73        77        71 
# chr5_q11  chr5_q12  chr5_q13  chr5_q14  chr5_q15  chr5_q21  chr5_q22  chr5_q23 
# 76        88        90        85        82        83        83        87 
# chr5_q31  chr5_q32  chr5_q33  chr5_q34  chr5_q35  chr8_p23  chr8_p22  chr8_p21 
# 81        83        84        87        88       108       108       106 
# chr8_p12  chr9_p24  chr9_p23  chr9_p22  chr9_p21 
# 80        86        84        84        82 



se = c("chr13_q14", "chr3_p13", "chr4_p16", "chr4_q34", "chr5_q13", "chr8_p22", "chr9_p21")


dat2 = dat2[, se, drop=F]
colnames(dat2) = paste(colnames(dat2), "__DEL", sep="")
dat2 = -dat2

res = cbind(dat1, dat2)
dim(res)
write.table(res, myoutf1, sep="\t", quote=F)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1.3] define weight profiles of chrosome bands
rm(list=ls())
myinf1 = "/FILEPATHPubDat/Dataset/Firehose/RNAseqv2_Hiseq/processed/LUSC_RNAseqv2_Tumor_Symbol.rda"
myinf2 = "/FILEPATHPubDat/Dataset/Firehose/Clinical/processed2/LUSC_Clincial_info.txt"
myinf3 = "/FILEPATH/lung-prelesion/output/TCGA_LUSC_Freq_ChrBand_V2.txt"

myoutf1 = "/FILEPATH/lung-prelesion/output/TCGA-LUNG-GenomicEvent-Profile_4base_V2.txt"

#----------------------------
soma = read.table(myinf3, sep="\t", header=T, row.names=1, quote="")
for(k in 1:ncol(soma))
{
  soma[,k] = ifelse(soma[,k]>0, 1, 0)
}

chr2_p15__AMP  chr3_q27__AMP  chr5_p14__AMP  chr8_q24__AMP chr13_q14__DEL 
55            271            148             58             69 
chr3_p13__DEL  chr4_p16__DEL  chr4_q34__DEL  chr5_q13__DEL  chr8_p22__DEL 
115             72             77             90            108 
chr9_p21__DEL 
82 

#---------------------------
load(file= myinf1)
data = mydata	## 262 samples
xx = apply(data>0, 1, sum)
se = which(xx>=10)## >0 in at least 10 samples
data = data[se,]
dim(data)
data = log10(data+1)

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
se = c("age_at_initial_pathologic_diagnosis", "gender", "stage_event.pathologic_stage", "number_pack_years_smoked")
info = info[,se]
colnames(info) = c("age", "gender", "stage", "smoker")
info$gender = as.character(info$gender)
info$age = as.numeric(info$age)
info$smoker = as.numeric(as.character(info$smoker))
xx = as.character(info$stage)
xx = gsub("stage ", "", xx)
se = grep("^iv", xx)
xx[se] = "IV"
se = grep("^iii", xx)
xx[se] = "III"
se = grep("^ii", xx)
xx[se] = "II"
se = grep("^i", xx)
xx[se] = "I"
info$stage = xx

comSam = intersect(colnames(data), row.names(info))
data = data[, comSam]
info = info[comSam,]
raw.info = info
raw.data = data
#~~~~~~~~~~~~~~~~~~~~
tmp = matrix(0, nrow(data), ncol(soma))
row.names(tmp) = row.names(data)
colnames(tmp) = colnames(soma)
p1 = p2 = p3 = p4  =  tmp
## uni.noj --> p1:   univariate using each genomic feature without adjust by clinical variables
## uni.adj --> p2:   univariate using each genomic feature adjusted by clinical variables
## mul.noj --> p3:   mulvariate using all genomic features without adjust by clinical variables
## mul.adj --> p4:   mulvariate using all genomic features adjusted by clinical variables


for(s in 1:ncol(soma))
{
  mut = row.names(soma)[soma[,s]>0]
  wt =  row.names(soma)[soma[,s]==0]	
  mut = mut[mut%in%colnames(data)]
  wt = wt[wt%in%colnames(data)]
  xx = c(mut, wt)
  info = raw.info[xx, ]
  data = raw.data[, xx]
  mut = ifelse(row.names(info)%in%wt, 0, 1)
  info = cbind(mut, info)
  
  pval1 = beta1 = rep(0, nrow(data))
  pval2 = beta2 = rep(0, nrow(data))
  
  for(k in 1:nrow(data))
  {
    cat("\r\r\r", s, "-->", k)
    mytf = as.numeric(data[k,])
    xx = cbind(mytf, info)
    mylog <- lm(mytf ~ mut, data = xx)
    mylog = summary(mylog)
    pval1[k] =  mylog[["coefficients"]][2, 4]
    beta1[k] =  mylog[["coefficients"]][2, 1]
    
    mylog <- lm(mytf ~ mut + age + gender + stage + smoker, data = xx, family = "binomial")
    mylog = summary(mylog)
    pval2[k] =  mylog[["coefficients"]][2, 4]
    beta2[k] =  mylog[["coefficients"]][2, 1]
  }
  xx = -log10(pval1)
  xx = ifelse(beta1>0, xx, -xx)
  p1[,s] = xx
  xx = -log10(pval2)
  xx = ifelse(beta2>0, xx, -xx)
  p2[,s] = xx		
}


#--------------------------
res = p1
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp1 = ifelse(tmp>0, tmp, 0)
  tmp1[tmp1>10] = 10
  res1[,k] = tmp1
  tmp2 = ifelse(tmp<0, -tmp, 0)
  tmp2[tmp2>10] = 10
  res2[,k] = tmp2
}
colnames(res1)= paste(colnames(res1), ".up", sep="")
colnames(res2)= paste(colnames(res2), ".dn", sep="")
myres = cbind(res1, res2)
row.names(myres) = row.names(res)
colnames(myres) = paste("uni.noj", colnames(myres), sep="__")
minv = min(myres)
maxv = max(myres)
myres = (myres-minv)/(maxv-minv)
prof1 = myres


res = p2
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp1 = ifelse(tmp>0, tmp, 0)
  tmp1[tmp1>10] = 10
  res1[,k] = tmp1
  tmp2 = ifelse(tmp<0, -tmp, 0)
  tmp2[tmp2>10] = 10
  res2[,k] = tmp2
}
colnames(res1)= paste(colnames(res1), ".up", sep="")
colnames(res2)= paste(colnames(res2), ".dn", sep="")
myres = cbind(res1, res2)
row.names(myres) = row.names(res)
colnames(myres) = paste("uni.adj", colnames(myres), sep="__")
minv = min(myres)
maxv = max(myres)
myres = (myres-minv)/(maxv-minv)
prof2 = myres


#------------------------------------------------------------
profile = cbind(prof1, prof2)
se1 = grep("\\.up", colnames(profile))
se2 = grep("\\.dn", colnames(profile))
profile = profile[, c(se1, se2)]
dim(profile)	## 20500

write.table(profile, myoutf1, sep="\t", quote=F)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1.4]
rm(list=ls())
myinf1 = "/FILEPATHWorSpa/m1_cancer/Lung/IntGen/data/TCGA_LUSC_Freq_GeneCNV.txt"
myinf2 = "/FILEPATH/lung-prelesion/output/TCGA_LUSC_Freq_ChrBand_V2.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
se = which(colnames(data)=="CDKN2A__DEL")
sam1 = row.names(data)[data[,se]>0]

data = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
se = which(colnames(data)=="chr9_p21__DEL")
sam2 = row.names(data)[data[,se]>0]

sum(sam2%in%sam1)
length(sam1)
length(sam2)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.5] Similarity between CDNA2A-DEL and chr9_p21__DEL
rm(list=ls())
myinf1 = "/FILEPATHWorSpa/m1_cancer/Lung/IntGen/data/TCGA-LUNG-LUSC-GenomicEvent-Profile_4base_Gene.txt"
myinf2 = "/FILEPATH/lung-prelesion/output/TCGA-LUNG-GenomicEvent-Profile_4base_V2.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
se = grep("uni.adj__CDKN2A__DEL", colnames(data))
tmp = data[,se]
sig1 = tmp[,1]- tmp[,2]
names(sig1) = row.names(data)

data = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
se = grep("uni.adj__chr9_p21__DEL", colnames(data))
tmp = data[,se]
sig2 = tmp[,1]- tmp[,2]
names(sig2) = row.names(data)

cor(sig1, sig2, method="s")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [3.1] with CDKN2A-DEL score
rm(list=ls())
myinf1 = "/FILEPATH/lung-prelesion/output/TCGA-LUSC_GenomicEvent_tmoRev_score.txt"
myinf2 = "/FILEPATHPubDat/Dataset/Firehose/CNV/processed/LUSC_CNV_Symbol.rda"
myinf3 = "/FILEPATHPubDat/organisms/human/annotation/Symbol_UCSC_hg19_Gene_Index.txt"
amp.thr = log2(2.8/2)
del.thr = log(1.2/2)
#------------------------
data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp
score = as.numeric(data[, "uni.noj__CDKN2A__DEL"])
names(score) = row.names(data)
#------------------------
info = read.table(myinf3, sep="\t", header=T, row.names=1)
info = info[info$Chr==9, ]
info = info[order(info$Sta), ]
se = which(row.names(info) =="CDKN2A")
info = info[(se-150):(se+150), ]
#------------------------
load(myinf2)
data = mydata
dim(data)
xx = as.numeric(substr(colnames(data), 14, 15))
se = which(xx==1)
data = data[,se]
dim(data)
colnames(data) = substr(colnames(data), 1, 12)
#------------------------
se = which(row.names(info) %in% row.names(data))
info = info[se,]
dim(info)
#------------------------
comxx = intersect(colnames(data), names(score))
data = data[, comxx]
score = score[comxx]
data = data[row.names(info), ]
data = ifelse(data<del.thr, 1, 0)
xx = apply(data, 1, sum)
plot(xx, pch=20)
se = which(row.names(data)=="CDKN2A")
abline(v=se)
myx = info$Sta
plot(myx, xx,  pch=20)
se = which(row.names(data)=="CDKN2A")
abline(v=info$Sta[se])
tag = data[se,]
se = which(tag==1)
dat1 = data[,se]		# samples with CDKN2A del
se = which(tag==0)	
dat2 = data[,se]		# samples without CDKN2A del
#------------------------
xx1 = score[colnames(dat1)]
xx2 = score[colnames(dat2)]
boxplot(list(xx1, xx2))
## -- patient with CDKN2A deletion have higher score
#------------------------
mysiz = apply(dat1, 2, sum)
myy = score[colnames(dat1)]
cor(mysiz, myy, method="s")
## 0.320558 --> the more genes covered by the deletion event in a sample, the higher the CDKN2A-DEL score
##score = as.numeric(data[, "uni.noj__CDKN2A__DEL"]) 	## uni.noj__chr9_p21
score = as.numeric(data[, "uni.noj__chr9_p21"]) 	## uni.noj__chr9_p21
## S= 0.320558 for CDKN2A-del--> the more genes covered by the deletion event in a sample, the higher the CDKN2A-DEL score
## S= 0.4725795 for chr9p21-del --> the more genes covered by the deletion event in a sample, the higher the CDKN2A-DEL score
# -- Indicate the deletion size is positively correlated with CDKN2A-DEL score, indicating it is a surrogate of chr9p21 deletion effect 
# (the reviewer is correct; consistent with previous report)
# -- the chr9p21-del better indicate broad deletion effect of this chromosome region
