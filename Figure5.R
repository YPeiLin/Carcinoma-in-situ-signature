#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [8.4] set variable as ordinal 
rm(list=ls())
myinf1 = "Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf2 = "Mascaux_GSE33479/Clinical_info.txt"


data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("\\.up\\.ES", "", colnames(data))

info = read.table(myinf2, sep="\t", header=T, quote="", row.names=1)
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = grep("uni.adj__", colnames(data))
data = data[,se]
colnames(data) = gsub("uni.adj__", "", colnames(data))

mycat = table(info$histology)
mycat = mycat[c(6, 2, 3, 4, 5, 7, 1, 8)]

res = matrix(0, ncol(data), 7)
row.names(res) = colnames(data)
res1 = res2 = res

for(k in 1:7)
{
  se = which(info$histology==names(mycat)[k])
  dat1 = data[se,]
  se = which(info$histology==names(mycat)[k+1])
  dat2 = data[se,]
  for(i in 1:ncol(data))
  {
    xx1 = dat1[,i]
    xx2 = dat2[,i]
    tmp = wilcox.test(xx1, xx2, alternative="l")$p.value
    res1[i,k] = tmp
    tmp = wilcox.test(xx1, xx2, alternative="g")$p.value
    res2[i,k] = tmp
  }
}

xx1 = apply(res1<0.05, 1, sum)
xx2 = apply(res2<0.05, 1, sum)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [8.5] correlation between signature scores and immune infiltration
rm(list=ls())
dev.off()
myinf1 = "Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf2 = "Mascaux_GSE33479_ImmGen_representative_TILs.txt"
myinf3 = "Mascaux_GSE33479/Clinical_info.txt"
library(reshape)
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
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("_up\\.ES", "", colnames(data))
Imm = data

comxx = intersect(row.names(data), row.names(info))
Sig = Sig[comxx,]
Imm = Imm[comxx,]
info = info[comxx,]

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

# type = "CD8T"
type = "Monocyte"
for(k in 1:ncol(res))
{
  se = which(info$histology==names(mycat)[k]) #stage
  scc = cor(Sig[se,], Imm[se, type], method="s")
  res[,k] = scc
  pval<- sapply(Sig[se,], calculate_correlation_p_value, vector = Imm[se, type])
  res.p[,k]=pval
}

res = round(res, 3)
xx = t(res)
se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL",  "TP53__MUT", "PIK3CA__MUT",  "KRAS__AMP",  "SOX2__AMP")) #CD8T
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL",  "TP53__MUT", "BCL6__AMP",  "KRAS__AMP",  "PIK3CA__AMP")) # MAC
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL", "MLL2__MUT","NF1__MUT", "TP53__MUT", "BCL6__AMP", "EGFR__AMP", "KRAS__AMP", "PIK3CA__MUT", "PIK3CA__AMP"))
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL", "TP53__MUT","MLL2__MUT","PIK3CA__AMP","NF1__MUT"))
xx = xx[,se]
cor= as.data.frame(xx)
cor$stage = row.names(cor)
melted_cor <- melt(cor, id.vars = "stage")
names(melted_cor) =c("stage","sig","cor")

xx = t(round(res.p,3))
se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL",  "TP53__MUT", "PIK3CA__MUT",  "KRAS__AMP",  "SOX2__AMP")) #CD8T
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL",  "TP53__MUT", "BCL6__AMP",  "KRAS__AMP",  "PIK3CA__AMP"))
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL", "MLL2__MUT","NF1__MUT", "TP53__MUT", "BCL6__AMP", "EGFR__AMP", "KRAS__AMP", "PIK3CA__MUT", "PIK3CA__AMP"))
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL", "TP53__MUT","MLL2__MUT","PIK3CA__AMP","NF1__MUT"))
xx = xx[,se]
p=as.data.frame(xx)

p$stage = row.names(p)
melted_p <- melt(p, id.vars = "stage")
names(melted_p) =c("stage","sig","pval")

df <- merge(melted_cor, melted_p,  by = c("stage","sig"), all = FALSE)
custom_order <- c("normal","hyperplasia","metaplasia", "mild dysplasia","moderate dysplasia", "severe dysplasia", "carcinoma in situ","squamous cell carcinoma")

#Remove certain stages
se=which(df$stage %in% c("normal","hyperplasia","moderate dysplasia","severe dysplasia"))
df= df[-se,]
se = which(custom_order %in%c("normal","hyperplasia","moderate dysplasia","severe dysplasia") )
custom_order = custom_order[-se]

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
  geom_point(aes( color = Signifcant,size = abs(cor)),shape=16) +
  scale_size_continuous(range = c(0.1, 13), limits=c(0.1,1), name = 'Correlation')+
  scale_color_manual(values=alpha(c('darkblue',"white",'darkred')))+
  # scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Reds")), name = 'Correlation') +
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


# [7.2] prdictive of CIS progression
rm(list=ls())
myinf1 = "Teixeira_GSE94611_ImmGen_representative_TILs.txt"
myinf2 = "Teixeira_GSE108124/GSE94611_Clinical_info.txt"

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
sort(myauc)

df = as.data.frame(myauc)
df$p= '-'
se =which(df$myauc>0.5)
df[se,]$p = "positive"
df[-se,]$p = "negative"
df[-se,]$myauc = 1- df[-se,]$myauc

pal= c('darkblue','darkred')

rownames(df)
df$name = gsub("\\_up.ES", "",rownames(df))
# bp = barplot(bueno$Bueno_coxph_p,col = col, space =0.5,  ylim=c(0,5) )        
# text(x=bp[,1], y=-0.3, adj=c(1, 1), bueno$name, cex=0.5, srt=65, xpd=TRUE)
# dev.off()
dev.off()
library(ggplot2)

g=ggplot(df, aes(x= reorder(name,-myauc), myauc,fill= p)) +    # ggplot2 plot with modified x-axis labels
  geom_bar(stat = "identity")  +ylim(0,1)+
  scale_fill_manual(values=pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g = g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=18,color="black"))
g=g+theme(axis.text.x = element_text(size=18,color="black",angle = 90,))
g=g+xlab("Immune Cell Type")+ylab("AUC progression")
g=g+labs(title = paste0("GSE157011"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.position = "none")

g
# 
# NavB_up.ES     MemB_up.ES     CD8T_up.ES   NKcell_up.ES Monocyte_up.ES 
# 0.2610294      0.3382353      0.3602941      0.4411765      0.5772059 
# CD4T_up.ES 
# 0.7794118 



#################################################################################
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
myinf1= "Mascaux_GSE33479_Symbol.txt"
data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)

gene_set<-read.table("immune_genesets.txt",sep="\t",header=T,row.names=1)
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
# col_ha = HeatmapAnnotation(df=so[,3:4],col=list(chr3p_DEL=c("0"="white","1"="#F0932B"),VHL__MUT=c("0"="white","1"="#22A6B3")))
pal<-pal_frontiers("default")(8)
names(pal)<-levels(gene_ss$marker)
row_ha = rowAnnotation( df=gene_ss,col=list(marker=pal),show_annotation_name=FALSE)
# ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
# ht<-Heatmap(scaled_mat,bottom_annotation = col_ha, right_annotation=row_ha,
#             show_column_names = FALSE, cluster_rows = FALSE,name="scaled expression",
#             col = col_fun,show_row_dend = FALSE,column_km=2, column_gap = unit(5,"mm"),
#             show_row_names=FALSE,column_title = c("immune hot", "immune cold"),
#             column_title_gp = gpar(fill = c("#E64B35", "#3C5488"),col="white"))
ht<-Heatmap(scaled_mat, right_annotation=row_ha,
            show_column_names = FALSE, cluster_rows = FALSE,name="expression",
            show_row_dend = FALSE,column_km=2, column_gap = unit(5,"mm"),
            show_row_names=FALSE,column_title = c("immune hot", "immune cold"),
            column_title_gp = gpar(fill = c("#E64B35", "#3C5488"),col="white"))
# column_order <- col_order(ht)
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


library(reshape)
myinf1 = "Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf3 = "Mascaux_GSE33479/Clinical_info.txt"
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
df= read.csv("amp-manual-group.csv")
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


#-----------------
# volcano plot
library(ggrepel)
library(stringr)
library(dplyr)
de = res
# de$pval =p.adjust(de$pval, "BH")

de$name = rownames(res)
de$logp=-log10(as.numeric(de$pval))
de$logp[which(de$logp>5)]=5
de$diffexpressed <- "-"
de$diffexpressed[de$Difference>0] <- "UP"
de$diffexpressed[de$Difference<0] <- "DOWN"
de$diffexpressed[de$pval>0.01] <- "Not significant"
# top10 = tail(arrange(de,logp),5)
# top10 =top10$name
# top10
palette  = c("darkred",'#999999',"darkblue")
de$label= "-"

se = which(de$name %in% de$name[which(de$p<0.01)])
de$label[se] =de$name[se] 
de$label[-se]=NA

dev.off()
g=ggplot(data=de, aes(x=as.numeric(Difference), y=as.numeric(logp),col=diffexpressed)) + 
  geom_point(size=2) + 
  scale_color_manual(values=palette)+
  geom_text_repel(aes(label = label),size = 4)+ 
  geom_vline(xintercept=0,  color = "black",size =1)+
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "blue",size =1)
g = g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g = g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=18,color="black"))
g=g+theme(axis.text.x = element_text(size=18,color="black"))
g=g+xlab("Difference")+ylab("P-value")
g=g+labs(title = paste0("GSE33479"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)+ theme(legend.position = "none") 
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())+ ylim(0, 5)+xlim(-40,40)
g
# g=g+theme(legend.text = element_text(size = 15))
print(g)
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
# gene = "NF1-MUT"
gene = "HEY1-AMP"
gene = "TERT-AMP"
# gene="PDE4DIP-MUT"
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
# [M1/M2]
rm(list=ls())
library(reshape)
myinf2 = "Mascaux_GSE33479_ImmGen_representative_TILs.txt"
myinf3 = "Mascaux_GSE33479/Clinical_info.txt"
myinf4 = "M1:M2.csv"


myinf1= "Mascaux_GSE33479_Symbol.txt"
data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)

info = read.table(myinf3, sep="\t", header=T, quote="", row.names=1)

m1 = read.csv(myinf4)
table(m1$M1_36750562) #20
table(m1$M2_36750562) #49
table(m1$M1)#9
table(m1$M2)#18
M1 = m1$X[which(m1$M1_36750562==1)]
M2 = m1$X[which(m1$M2_36750562==1)]
comxx = intersect(row.names(data), toupper(M1)) #19
M1.data= data[comxx,] 


comxx = intersect(row.names(data), toupper(M2)) #41
M2.data= data[comxx,]
rownames(M1.data) = paste0("M1",rownames(M1.data))
rownames(M2.data) = paste0("M2",rownames(M2.data))
#---------------
myx =rbind(M1.data,M2.data) 
hclust_rows <- as.dendrogram(hclust(dist(myx)))
hclust_cols <- as.dendrogram(hclust(dist(t(myx))))

library(RColorBrewer)
se= which(rownames(myx) %in% rownames(M1.data))
rownames(myx)[se]=paste0("M1 ", rownames(myx)[se])
rownames(myx)[-se]=paste0("M2 ", rownames(myx)[-se])

col = colorRampPalette(c("#FDD586","#313695"))(30)
info$phenotypical.stage[which(info$phenotypical.stage=="normal hypofluorescent" | info$phenotypical.stage=="normal normofluorescent")] ="normal"
se = which(colnames(myx)%in%rownames(info[which(info$phenotypical.stage=="normal"),]))
colnames(myx)[se]=paste0("normal",rep(1:length(se)))
se = which(colnames(myx)%in%rownames(info[which(info$phenotypical.stage=="hyperplasia"),]))
colnames(myx)[se]=paste0("hyperplasia",rep(1:length(se)))
se = which(colnames(myx)%in%rownames(info[which(info$phenotypical.stage=="metaplasia"),]))
colnames(myx)[se]=paste0("metaplasia ",rep(1:length(se)))

se = which(colnames(myx)%in%rownames(info[which(info$phenotypical.stage=="mild dysplasia"),]))
colnames(myx)[se]=paste0("metaplasia ",rep(1:length(se)))

heatmap(as.matrix(myx),   
        scale = "column",
        col = col,
        
        Rowv = NA,
        Colv = NA,
        RowSideColors = c(rep("#FEE090",19), rep("#313695", 41)))


mycat = table(info$histology)
mycat = mycat[c(6, 2, 3, 4, 5, 7, 1, 8)]

mycat
res = matrix(0, nrow(myx), length(mycat))
row.names(res) = rownames(myx)
colnames(res) = names(mycat)
res = as.data.frame(res)

res.p = matrix(0, nrow(myx), length(mycat))
row.names(res.p) = colnames(Sig)
colnames(res.p) = names(mycat)
res.p = as.data.frame(res.p)
calculate_correlation_p_value <- function(col, vector) {
  cor_test_result <- cor.test(col, vector, method = "s")
  return(cor_test_result$p.value)
}


type = "Monocyte"
# Monocyte
for(k in 1:ncol(res))
{
  se = which(info$histology==names(mycat)[k]) #stage
  scc = cor(Sig[se,], Imm[se, type], method="p")
  res[,k] = scc
  pval<- sapply(Sig[se,], calculate_correlation_p_value, vector = Imm[se, type])
  res.p[,k]=pval
}

res = round(res, 3)
xx = t(res)
se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL",  "TP53__MUT", "PIK3CA__MUT",  "KRAS__AMP",  "PIK3CA__AMP")) #CD8T
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL",  "TP53__MUT", "BCL6__AMP",  "KRAS__AMP",  "PIK3CA__AMP")) # MAC
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL", "MLL2__MUT","NF1__MUT", "TP53__MUT", "BCL6__AMP", "EGFR__AMP", "KRAS__AMP", "PIK3CA__MUT", "PIK3CA__AMP"))
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL", "TP53__MUT","MLL2__MUT","PIK3CA__AMP","NF1__MUT"))
xx = xx[,se]
cor= as.data.frame(xx)
cor$stage = row.names(cor)
melted_cor <- melt(cor, id.vars = "stage")
names(melted_cor) =c("stage","sig","cor")

xx = t(round(res.p,3))
se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL",  "TP53__MUT", "PIK3CA__MUT",  "KRAS__AMP",  "PIK3CA__AMP")) #CD8T
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL",  "TP53__MUT", "BCL6__AMP",  "KRAS__AMP",  "PIK3CA__AMP"))
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL", "MLL2__MUT","NF1__MUT", "TP53__MUT", "BCL6__AMP", "EGFR__AMP", "KRAS__AMP", "PIK3CA__MUT", "PIK3CA__AMP"))
# se = which(colnames(xx)%in%c("CDKN2A__MUT", "CDKN2A__DEL", "TP53__MUT","MLL2__MUT","PIK3CA__AMP","NF1__MUT"))
xx = xx[,se]
p=as.data.frame(xx)

p$stage = row.names(p)
melted_p <- melt(p, id.vars = "stage")
names(melted_p) =c("stage","sig","pval")

df <- merge(melted_cor, melted_p,  by = c("stage","sig"), all = FALSE)
custom_order <- c("normal","hyperplasia","metaplasia", "mild dysplasia","moderate dysplasia", "severe dysplasia", "carcinoma in situ","squamous cell carcinoma")

#Remove certain stages
se=which(df$stage %in% c("normal","hyperplasia","moderate dysplasia","severe dysplasia"))
df= df[-se,]
se = which(custom_order %in%c("normal","hyperplasia","moderate dysplasia","severe dysplasia") )
custom_order = custom_order[-se]

# Order the "Category" column with the custom order
df$stage <- factor(df$stage, levels = custom_order)




#################################################################################
#CD274
rm(list=ls())
#---------
#Chaoprocess
myinf1 = "Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf2 = "Mascaux_GSE33479/Mascaux_GSE33479_Symbol.txt"
myinf3 = "Mascaux_GSE33479/Clinical_info.txt"
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
se = which(colnames(xx)%in%c("CDKN2A__MUT", "TP53__MUT", "SOX2__AMP", "EGFR__AMP", "KRAS__AMP","EGFR__AMP","CDKN2A__DEL","NFE2L2__MUT" ))
xx = xx[,se]
cor= as.data.frame(xx)
cor$stage = row.names(cor)
melted_cor <- melt(cor, id.vars = "stage")
names(melted_cor) =c("stage","sig","cor")

xx = t(round(res.p,3))
se = which(colnames(xx)%in%c("CDKN2A__MUT", "TP53__MUT", "SOX2__AMP", "EGFR__AMP", "KRAS__AMP","EGFR__AMP","CDKN2A__DEL","NFE2L2__MUT" ))
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

