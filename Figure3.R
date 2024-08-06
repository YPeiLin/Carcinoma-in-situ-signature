################################################################################
rm(list=ls())
myinf1 = "Teixeira_GSE94611_GenomicEvent_iRAS_SQC.txt"
myinf2 = "Teixeira_GSE108124/GSE94611_Clinical_info.txt"

data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("\\.up\\.ES", "", colnames(data))

#---------
myinf3 = "Sampleinfo.csv"
myinf4 = "Teixeira_GSE108124/Teixeira_mutation.txt"
# mutation compare
sample = read.csv(myinf3,sep=",",header=TRUE)
mut.df = read.table(myinf4,sep="\t",header=TRUE)
se = grep("uni.adj__", names(data))
data= data[,se]
colnames(data)= gsub("uni.adj__", "",colnames(data))

se=grep("Progressive.",sample$Sample.Number..GXN.)
sam.p = sample[se,]
se=grep("Regressive.",sample$Sample.Number..GXN.)
sam.r = sample[se,]
#Recombining only 33 patients
sample = rbind(sam.p,sam.r)

mut="TP53"
se = which(mut.df$Gene==mut)
mut.df = mut.df[se,]
se= intersect( sample$Sample.Number..WGS.,mut.df$Sample) #TP53 36 samples -> 6 overlap; PIK3CA: 3 

mut.sam= sample[which(sample$Sample.Number..WGS.%in%se),]$Sample.Number..GXN.
r1 = paste0("Regressive.",1:16)
r2 = paste0("Progressive.", 1:17)
row.names(data)=append(r1,r2)

#Draw boxplot
tmp=data
tmp$group = rownames(data)%in% mut.sam
dat=tmp[,c(paste0(mut,"__MUT"),"group")]
colnames(dat) = c("score","group")

dat$score = (dat$score - mean(dat$score))/sd(dat$score)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
pal=c( "#CCEBC5","#DECBE4")

se=which(dat$group=="WO (n=27)")
tmp=t.test(dat$score[se],dat$score[-se],alternative = "less")

g_order=names(sort(tapply(dat$score,dat$group,median,na.rm=TRUE), decreasing =TRUE))
dat$group<- factor(dat$group,levels=g_order)
levels(dat$group)= c(paste0(mut," (n=",table(dat$group)[1], ")"), paste0("WO (n=",table(dat$group)[2],")"))
# levels(dat$group)= c(paste0("WO"," (n=",table(dat$group)[1], ")"), paste0(mut," (n=",table(dat$group)[2], ")"))
my_comparisons <- list( c("WO (n=27)","TP53 (n=6)"))
stat.test <- dat %>% 
  wilcox_test(score ~ group, comparisons = my_comparisons,alternative = "less") %>%
  add_xy_position(x = "group") %>%
  mutate(myformatted.p = paste0("p = ", ifelse(p < 1e-5, format(p, scientific = T), signif(p, digits = 1))))
g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA,width=0.5/length(unique(dat$group)))
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Signature score")
g=g+labs(title = paste0("GSE108124_",mut))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 4/3)+ylim(-1.1,1.5)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g+ geom_text(x=2, y=1.4, label=round(tmp$p.value,3),size=4)
print(g)
dev.off()

dat2=dat[which(dat$group=="WO (n=27)"),]

#----------------------

########
# WT p/r boxplot
rm(list=ls())
myinf1 = "Teixeira_GSE94611_GenomicEvent_iRAS_SQC.txt"
myinf2 = "Teixeira_GSE108124/GSE94611_Clinical_info.txt"

data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(cnum*2)]
colnames(data) = gsub("\\.up\\.ES", "", colnames(data))
myinf3 = "Sampleinfo.csv"
myinf4 = "Teixeira_mutation.txt"
# mutation compare
sample = read.csv(myinf3,sep=",",header=TRUE)
mut.df = read.table(myinf4,sep="\t",header=TRUE)
se = grep("uni.adj__", names(data))
data= data[,se]
colnames(data)= gsub("uni.adj__", "",colnames(data))

se=grep("Progressive.",sample$Sample.Number..GXN.)
sam.p = sample[se,]
se=grep("Regressive.",sample$Sample.Number..GXN.)
sam.r = sample[se,]
#Recombining only 33 patients
sample = rbind(sam.p,sam.r)

mut="TP53"
se = which(mut.df$Gene==mut)
mut.df = mut.df[se,]
se= intersect( sample$Sample.Number..WGS.,mut.df$Sample) #TP53 36 samples -> 6 overlap

mut.sam= sample[which(sample$Sample.Number..WGS.%in%se),]$Sample.Number..GXN.
r1 = paste0("Regressive.",1:16)
r2 = paste0("Progressive.", 1:17)
row.names(data)=append(r1,r2)

#Draw boxplot
tmp=data
tmp$group = rownames(data)%in% mut.sam
dat=tmp[,c(paste0(mut,"__MUT"),"group")]
colnames(dat) = c("score","group")

dat$score = (dat$score - mean(dat$score))/sd(dat$score)
dat=dat[which(dat$group==FALSE),]
dat$group=rownames(dat)
dat$group=gsub("\\..*","",dat$group)


library(ggplot2)
library(ggpubr)
library(RColorBrewer)
pal=c( "#FDD586","#B0C9DE")


se=which(dat$group=="Progressive")
tmp=t.test(dat$score[se],dat$score[-se])
tmp

g_order=names(sort(tapply(dat$score,dat$group,median,na.rm=TRUE), decreasing =TRUE))
dat$group<- factor(dat$group,levels=g_order)
g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA,width=0.5/length(unique(dat$group)))
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Progression status")
g=g+labs(title = paste0("GSE108124_TP53WT"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 4/3)
g=g+ylim(-1.1,1.5)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g+ geom_text(x=2, y=1.4, label=round(tmp$p.value,3),size=4)
print(g)
dev.off()



################################################################################
# heatmap by stage
rm(list=ls())

myinf1 = "Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf2 = "Mascaux_GSE33479/Clinical_info.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp

info = read.table(myinf2, sep="\t", header=T, row.names=NULL)
info$sampleid = info$row.names


################################################################################
#Boxplot
name = colnames(data)
se = grep("uni.adj__", name)
uni.data = data[,se]
data=uni.data

data$sampleid = row.names(data)
data=merge(data,info,by.y="sampleid")
data$phenotypical.stage[which(data$phenotypical.stage=="normal hypofluorescent" | data$phenotypical.stage=="normal normofluorescent")] ="normal"
row.names(data)=(data$sampleid)

myx = data[,2:57]

res=NULL

for(i in 1:ncol(myx)){
  score = myx[,i]
  # xx = (score) / sd(abs(score))
  # xx= score
xx=   (score-mean (score) )/ sd(score)
  res = cbind(res,xx)
}

res[res>4] = 4
res[res< -4] = -4
## remove extreme values
thr = 4
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp[tmp<(-thr)] = -thr
  tmp[tmp>thr] = thr
  res[,k] = tmp
}


hist(res)
colnames(res) = colnames(myx)
rownames(res) = rownames(myx)
myx=res
library(RColorBrewer)
str_sub(colnames(myx),start=-5, end=-4) <- "-" 
colnames(myx) =gsub("uni.adj__", "",colnames(myx))
#------------


df= read.csv("amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")

se1=which(colnames(myx)%in%df$rep)
se2=which(!colnames(myx)%in%df$entrezgene_accession)
myx = myx[,c(se1,se2)]


hclust_rows <- as.dendrogram(hclust(dist(myx)))
hclust_cols <- as.dendrogram(hclust(dist(t(myx))))

# heatmap(as.matrix(myx),Rowv = NA,Colv = NA)

library(RColorBrewer)
red=rgb(1,0,0); green=rgb(0,1,0); blue=rgb(0,0,1); white=rgb(1,1,1)
GtoWrange<-colorRampPalette(c(green, white ) )
WtoRrange<-colorRampPalette(c(white, red) )

# col =colorRampPalette(brewer.pal(8, "RdYlBu"))(30)
col =(brewer.pal(8, "RdYlBu"))
library(heatmap)
heatmap(as.matrix(myx),   
        scale = "column",
        col = colorRampPalette(c("darkblue","#D73027"))(10),
        Colv = hclust_cols,
        # Colv = NA,
        # Rowv = hclust_rows,
        Rowv = NA,
        RowSideColors = c(rep(col[8],27), rep(col[7], 15), rep(col[6],15),rep(col[5],13), rep(col[4],13), rep(col[3],12),rep(col[2],13),rep(col[1],14)))

dev.off()

#############################################
# Figure 3
# [8.6] heatmap
rm(list=ls())
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
input1= "Figure_merge.csv"
figure=read.csv(input1, header = T, row.names = 1)

#readfrom above
myx= as.data.frame(myx)
myx[123:124,]=figure[123:124,]
figure= myx

## clinical information
myinf2 = "Mascaux_GSE33479/Clinical_info.txt"
info = read.table(myinf2, sep="\t", header=T, row.names=NULL)
info$sampleid = info$row.names
mycat = table(info$histology)
mycat = mycat[c(6, 2, 3, 4, 5, 7, 1, 8)] # rearrangement of the table names according to the histology normal to advance level
col8= c(brewer.pal(8, "Reds"))
mycol = rep(0, nrow(data))
for(k in 1:length(mycat))
{
  se = which(info$histology==names(mycat)[k])
  mycol[se] = col8[k]
}
sid.col = mycol
## bottom color scaling 
botm.col=figure[123,]
bot_p=ifelse(botm.col<=1.3, "grey","#026c80")


library(gplots)
library(marray)
mycol <- maPalette(low="darkblue", mid= "white", high="darkred", k=40)
data=figure[1:122,]
# class(data)
data <- apply(data, 2, function(x) as.numeric(as.character(x)))

tmp = data
# str_sub(colnames(tmp),start=-5, end=-4) <- "-" 
df= read.csv("/amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")
se1=which(colnames(tmp)%in%df$rep)
se2=which(!colnames(tmp)%in%df$entrezgene_accession)
tmp = tmp[,c(se1,se2)]
data=tmp

# str_sub(colnames(bot_p),start=-5, end=-4) <- "-" 
se1=which(colnames(bot_p)%in%df$rep)
se2=which(!colnames(bot_p)%in%df$entrezgene_accession)
bot_p = bot_p[c(se1,se2)]

#------------ 
# # Convert 'bot_p' to a binary vector (0 and 1)
binary_bot_p = as.integer(bot_p == "#026c80") # Assuming "#8B0000" represents 1
# length(which(bot_p!="grey"))# 23 signatures
# # Create a sorting index based on the binary vector
sort_index = order(binary_bot_p)
# 
# # Reorder your data matrix based on the sorting index
sorted_data = data[, sort_index]
# 
# # Reorder 'bot_p' to match the new column order
sorted_bot_p = bot_p[sort_index]
# 
# # Now plot the heatmap with the reordered data
heatmap.2(as.matrix(sorted_data), trace="none", Colv=NA, margins=c(5,5),
          ColSideColors = sorted_bot_p)

dev.off()
maxx=heatmap.2(as.matrix(data),
                 
                 # deprogram control
                 Rowv = FALSE,
                 Colv= TRUE,
                 distfun = dist,
                 hclustfun = hclust,
                 dendrogram = "column", 						##c("both","row","column","none"),
                 symm = FALSE,
                 
                 # data scaling
                 scale = "none",
                 na.rm=TRUE,
                 col= mycol,
                 
                 # level trace
                 trace= "none",
                 
                 # Row/Column Labeling
                 margins = c(5, 1),
                 ColSideColors = bot_p,
                 RowSideColors = sid.col,
           # cexRow = 1,
           # cexCol = 1,
                 
                 # color key + density info
                 key = TRUE,
                 keysize = 1.5,
                 density.info= "none",
                 labRow = "",
                 labCol = colnames(data),
                 
                 
                 # plot labels
                 main = NULL,
                 xlab = NULL,
                 ylab = NULL
                 
)

#----------------------
# 
# # Convert 'bot_p' to a binary vector (0 and 1)
binary_bot_p = as.integer(bot_p == "#026c80") # Assuming "#8B0000" represents 1
# 

dev.off()
se=which(binary_bot_p==1)
maxx=heatmap.2(as.matrix(data[,se]),
               
               # deprogram control
               Rowv = FALSE,
               Colv= TRUE,
               distfun = dist,
               hclustfun = hclust,
               dendrogram = "column", 						##c("both","row","column","none"),
               symm = FALSE,
               
               # data scaling
               scale = "none",
               na.rm=TRUE,
               col= mycol,
               
               # level trace
               trace= "none",
               
               # Row/Column Labeling
               margins = c(5, 1),
               RowSideColors = sid.col,
               # cexRow = 1,
               # cexCol = 1,
               
               # color key + density info
               key = TRUE,
               keysize = 1.5,
               density.info= "none",
               labRow = "",
               labCol = colnames(data)[se],
               
               
               # plot labels
               main = NULL,
               xlab = NULL,
               ylab = NULL
               
)
dev.off()
########past
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

col8= c(brewer.pal(8, "Reds"))
mycol = rep(0, nrow(data))
for(k in 1:length(mycat))
{
  se = which(info$histology==names(mycat)[k])
  mycol[se] = col8[k]
}
sid.col = mycol

##  scale normalization
myavg = 0
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = (data[,k])/mystd[k]
}

## remove extreme values
thr = 4
for(k in 1:ncol(data))
{
  tmp = data[,k]
  tmp[tmp<(-thr)] = -thr
  tmp[tmp>thr] = thr
  data[,k] = tmp
}
library(gplots)
library(marray)
mycol <- maPalette(low="darkblue", mid= "white", high="darkred", k=40)
myxx= heatmap.2 (as.matrix(data),
                 
                 # dendrogram control
                 Rowv = FALSE,
                 Colv = TRUE,
                 distfun = dist,
                 hclustfun = hclust,
                 dendrogram = "column", 						##c("both","row","column","none"),
                 symm = FALSE,
                 
                 # data scaling
                 scale = "none",
                 na.rm=TRUE,
                 col= mycol,
                 
                 # level trace
                 trace= "none",
                 
                 # Row/Column Labeling
                 margins = c(2, 3),
                 
                 # color key + density info
                 key = TRUE,
                 keysize = 1.5,
                 density.info= "none",
                 labRow = "",
                 labCol = "",
                 
                 
                 # plot labels
                 main = NULL,
                 xlab = NULL,
                 ylab = NULL,
                 RowSideColors = sid.col
                 
)
dev.off()
#################################################################################
# [Trend boxplot]
rm(list=ls())
myinf1 = "Mascaux_GSE33479_GenomicEvent_iRAS.txt"
myinf2 = "Mascaux_GSE33479/Clinical_info.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp

info = read.table(myinf2, sep="\t", header=T, row.names=NULL)
info$sampleid = info$row.names

################################################################################
#Boxplot
name = colnames(data)
se = grep("uni.adj__", name)
uni.data = data[,se]
data=uni.data
##  scale normalization
myavg = 0
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = (data[,k])/mystd[k]
}

## remove extreme values
thr = 4
for(k in 1:ncol(data))
{
  tmp = data[,k]
  tmp[tmp<(-thr)] = -thr
  tmp[tmp>thr] = thr
  data[,k] = tmp
}

data$sampleid = row.names(data)
data=merge(data,info,by.y="sampleid")
data=data[,-1]

normal = data[which(data$phenotypical.stage=="normal hypofluorescent" | data$phenotypical.stage=="normal normofluorescent" ),]
hyperps = data[which(data$phenotypical.stage=="hyperplasia"),]
metaps = data[which(data$phenotypical.stage=="metaplasia"),]
mild.dys = data[which(data$phenotypical.stage=="mild dysplasia"),]
moderate.dys = data[which(data$phenotypical.stage=="moderate dysplasia"),]
severe.dys = data[which(data$phenotypical.stage=="severe dysplasia"),]
cis = data[which(data$phenotypical.stage=="carcinoma in situ"),]
scc= data[which(data$phenotypical.stage=="squamous cell carcinoma"),]


# For example, the CDKN2A-DEL signature had relatively low scores in normal lung and hyperplasia samples 
# but increased significantly at the metaplasia (P < xxx, Wilcox test compared to hyperplasia) and mild dysplasia 
# (P < xxx, Wilcox test compared to metaplasia) stages. 
a= append(normal$uni.adj__CDKN2A__DEL,hyperps$uni.adj__CDKN2A__DEL)

t.test(a, metaps$uni.adj__CDKN2A__DEL) # p-value = 0.0004161
t.test(a, mild.dys$uni.adj__CDKN2A__DEL) #p-value = 0.002748
#default boxplot
dat=as.data.frame(data$uni.adj__CDKN2A__DEL)
# dat=as.data.frame(data$uni.adj__SOX2__AMP)
dat$name = data$histology
names(dat)=c("score","group")

library(Kendall)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
boxplot(score~reorder(group,score, FUN=mean),data=dat,las=2)
summary(MannKendall(dat$score))
tmp=MannKendall(dat$score)


pval = "< 2.22e-16"
pal =(brewer.pal(8, "Reds"))
# pal=c("#D73027", "#F46D43", "#FDAE61", "#FEE090", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")
# g_order=names(sort(tapply(dat$score,dat$group,median,na.rm=TRUE), decreasing =FALSE))
g_order =c("normal",  "hyperplasia", "metaplasia", "mild dysplasia", "moderate dysplasia","severe dysplasia" , "carcinoma in situ" ,"squamous cell carcinoma")
dat$group<- factor(dat$group,levels=g_order)
# levels(dat$group)= c(paste0(mut," (n=",table(dat$group)[1], ")"), paste0("WO (n=",table(dat$group)[2],")"))
# levels(dat$group)= c(paste0("WO"," (n=",table(dat$group)[1], ")"), paste0(mut," (n=",table(dat$group)[2], ")"))
g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=12,color="black"))
g=g+xlab("Patient Group")+ylab("Signature score")
g=g+labs(title = paste0("GSE33479_CDKN2A-DEL"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 9/16)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g+geom_text(x=5.5, y=-2.5, label=paste0('MannKendall Trend \n (P-value ', pval, ')'),size=4)
print(g)
dev.off()

