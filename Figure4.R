#########################
# [8.3] chao heatmap (4a)
rm(list=ls())
myinf1 = "Teixeira_GSE108124_GenomicEvent_iRAS_SQC.txt"
myinf2 = "Teixeira_GSE108124/GSE94611_Clinical_info.txt"

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
#----------------------
library(marray)
library(gplots)
# mycol <- maPalette(low="#313695", mid= "white", high="#FDD586", k=40)
mycol <- maPalette(low="darkblue", mid= "white", high="darkred", k=40)
xx = data
for(k in 1:ncol(xx))
{
  tmp = xx[,k]
  tmp = tmp/mean(abs(tmp))
  tmp[tmp>3]=3
  tmp[tmp<(-3)] = -3
  xx[,k] = tmp
}
str_sub(colnames(xx),start=-5, end=-4) <- "-" 
df= read.csv("amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")
se1=which(colnames(xx)%in%df$rep)
se2=which(!colnames(xx)%in%df$entrezgene_accession)

xx = xx[,c(se1,se2)]

sam.col = rep("#FDD586", nrow(xx))
se = grep("Regressive", info)
sam.col[se] = "#B0C9DE"


myxx= heatmap.2 (as.matrix(t(xx)),
                 
                 # dendrogram control
                 Rowv = TRUE,
                 distfun = dist,
                 hclustfun = hclust,
                 dendrogram = "both",
                 symm = FALSE,
                 
                 # data scaling
                 scale = "none",
                 na.rm=TRUE,
                 col= mycol,
                 
                 # level trace
                 trace= "none",
                 
                 # Row/Column Labeling
                 margins = c(5, 12),
                 
                 # color key + density info
                 key = TRUE,
                 keysize = 1.2,
                 density.info= "none",
                 # labRow = "",
                 labCol = "",
                 
                 
                 # plot labels
                 main = NULL,
                 xlab = NULL,
                 ylab = NULL,
                 
                 ColSideColors = sam.col
                 
)
dev.off()
graphics.off() 
par("mar") 
par(mar=c(1,1,1,1)) 

###########################################################
rm(list=ls())
myinf1 = "Teixeira_GSE94611_GenomicEvent_iRAS_SQC.txt"
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

se = which(row.names(data)%in%sam1)
dat1 = data[se,]

se = which(row.names(data)%in%sam2)
dat2 = data[se,]

#AUC prediction
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
tmp=(myauc[se])
tmp=as.data.frame(tmp)
str_sub(colnames(data),start=-5, end=-4) <- "-" 
rownames(tmp) =gsub("uni.adj__", "",names(data)[se])
tmp$name = gsub("uni.adj__", "",names(data)[se])

df= read.csv("amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")
se1=which(tmp$name%in%df$rep)
se2=which(!tmp$name%in%df$entrezgene_accession)

tmp = tmp[c(se1,se2),]


df = as.data.frame(tmp)
df$p= '-'
df
se =which(df$tmp>0.5)
df[se,]$p = "positive"
df[-se,]$p = "negative"
df[-se,]$tmp = 1- df[-se,]$tmp
rownames(df)

pal= c('darkblue','darkred')
library(ggplot2)
g=ggplot(df, aes(x= reorder(name,-tmp), tmp,fill= p)) +    # ggplot2 plot with modified x-axis labels
  geom_bar(stat = "identity") +geom_hline(yintercept=0.7, linetype="dashed", color = "blue")+ylim(0,1)+
  scale_fill_manual(values=pal)+
  ylab("AUC")+geom_hline(yintercept=0.95, linetype="dashed", color = "red")
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g = g+theme(panel.border = element_rect(size=0.5,fill=NA))
# g=g+theme(axis.text.y = element_text(size=18,color="black"))
# g=g+theme(axis.text.x = element_text(size=18,color="black",angle = 90,))
g=g+xlab("Signature")+ylab("AUC progression")
g=g+labs(title = paste0("GSE157011"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.position = "none")
g=g+theme(axis.text.y = element_text(size=8,color="black"))
g=g+theme(axis.text.x = element_text(size=8,color="black",angle = 90))
g=g+theme(aspect.ratio = 9/16)
g
