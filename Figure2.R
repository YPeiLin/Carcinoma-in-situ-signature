# [3.2] survival plot
rm(list=ls())
myinf1 = "GenomicEvent_iRAS_SQC.txt"
myinf2 = "Clinical_info.txt"

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = as.numeric(info$os_mo)
e.surv = info$os_event
info = cbind(t.surv, e.surv, info)
tag = !is.na(info$t.surv)
info = info[tag==1,]
dim(info) #484 10

# info[1,]
se = grep("^CAD", info$Title)
info = info[-se,]

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp
##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(info)		## 

se = grep("uni.adj__", colnames(data))
data = data[,se]
#------------------------------------------------------ 

info=info[,c("t.surv","e.surv","age","Sex","Stage")]


#---------------
score = data$uni.adj__CDKN2A__DEL
scale.score = (score -mean(score))/sd(score)
xx= info
mycox = coxph(Surv(t.surv, e.surv)~scale.score, data = xx) 
summary(mycox)


#Divide risk score into 3 groups

xx = cbind(scale.score, info)

xx$range =xx$scale.score

se = which(xx$scale.score<=mean(scale.score))
xx[se,]$range ="low-value"
se = which(xx$scale.score>mean(scale.score))
xx[se,]$range ="high-value"

model_fit <- survfit(Surv(t.surv, e.surv)~range, data = xx)
g=ggsurvplot(
  model_fit,
  data = xx,
  size = 1,                 # change line size
  pval = TRUE,  
  pval.coord = c(0, 0.03),
  palette  = c("#F8766D","#4d4d4d"),
  title = paste0("CDKN2A-DEL"),
  legend.labs =
    c("High-value","Low-value"),    # Change legend labels
  ggtheme = theme(panel.grid =element_blank())+
    theme(panel.background = element_blank(),
          panel.border = element_rect(size=0.5,fill=NA),
          axis.text.y = element_text(size=20,color="black"),
          axis.text.x = element_text(size=17,color="black"),aspect.ratio = 1,
          legend.text = element_text(size = 15),
          legend.title=element_blank(),
          legend.background = element_rect(color=NA,fill=NA),
          title=element_text(face="bold",size=20),
    ),    
  xlab = "\n Survival Time (days)",
  ylab ="Overall Survival",
)
g


xx = cbind(scale.score,info)
yy=xx

library(gridExtra)
library(ggplot2)

mycox = coxph(Surv(t.surv, e.surv)~scale.score+age+Stage+Sex, data = yy) 
mycox
summary(mycox)

print(ggforest(mycox,data=yy))

###############################################################################
# volcano plot - univariable survreg
rm(list=ls())
myinf1 = "Bueno_GSE157011_GenomicEvent_iRAS_SQC.txt"
myinf2 = "Clinical_info.txt"

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = as.numeric(info$os_mo)
e.surv = info$os_event
info = cbind(t.surv, e.surv, info)
tag = !is.na(info$t.surv)
info = info[tag==1,]
dim(info) #484 10

se = grep("^CAD", info$Title)
info = info[-se,]

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp
##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(info)		## 

se = grep("uni.adj__", colnames(data))
data = data[,se]
#------------------------------------------------------ 
info=info[,c("t.surv","e.surv","age","Sex","Stage")]
xx = cbind(info, data)
xx=xx[which(xx$t.surv>0),]
myx = xx[, 6:ncol(xx)]

#----------------------------------------------------------------
library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))

for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  score=mytf
  # scale.score = (score-min(score))/(max(score)-min(score))
  scale.score = (score - mean(score))/ sd((score))
  # xx = cbind(mytf, info)
  xx = cbind(scale.score, info)
  # xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  # mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = survreg(Surv(t.surv, e.surv)~scale.score, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["scale.score", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~scale.score, xx) 
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


tmp=res$name
str_sub(tmp,start=-5, end=-4) <- "-" 
tmp = gsub("uni.adj__", "",tmp)

df= read.csv("/mount/amos1/home/ylin/lung-prelesion/output/amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")
se1=which(tmp%in%df$rep)
se2=which(!tmp%in%df$entrezgene_accession)
tmp = tmp[c(se1,se2)]

str_sub(res$name,start=-5, end=-4) <- "-" 
res$name = gsub("uni.adj__", "",res$name)
length(which(res[which(res$name%in%tmp),]$survreg.pval1<0.05)) #14 out of 34
length(which(res[which(res$name%in%tmp),]$coxph.pval1<0.05)) 
length(which(res$coxph.qval1<0.05)) #adjusted p-value

################################################################################
# volcano plot
library(ggrepel)
library(stringr)
library(dplyr)
se = grep("uni.adj__", res$name)
de = res[se,]
de$name = gsub("uni.adj__", "",res$name[se])
str_sub(de$name,start=-5, end=-4) <- "-"

#reduce redundancy
df= read.csv("/mount/amos1/home/ylin/lung-prelesion/output/amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")
se1=which(de$name%in%df$rep)
se2=which(!de$name%in%df$entrezgene_accession)
de=de[c(se1,se2),]

de$logp=-log10(as.numeric(de$survreg.pval1))

de$diffexpressed <- "-"
de$diffexpressed[de$hr1>1] <- "UP"
de$diffexpressed[de$hr1<1] <- "DOWN"
de$diffexpressed[de$survreg.pval1>0.05] <- "Not significant"
de$type = de$name
se = grep("-AMP", de$name)
de$type[se] ="Amplification"
se = grep("-DEL", de$name)
de$type[se]="Deletion"
se = grep("-MUT", de$name)
de$type[se]="Mutation"
# palette  = c('#999999',"#e64B35FF","#4dbbd5ff")
shape = c(15,6,4)



top10 = tail(arrange(de,logp),5)
top10 =top10$name
# top10

de$label= "-"
# se = which(de$name %in% top10)
se = which(de$survreg.pval1<0.05)
de$label[se] =de$name[se] 
de$label[-se]=NA

# dev.off()
g=ggplot(data=de, aes(x=as.numeric(hr1), y=logp, label=name,shape=type)) + 
  geom_point(size=2) + 
  scale_shape_manual(values=shape)+
  # scale_color_manual(values=palette)+
  geom_text_repel(aes(label = label),
                  size = 4)+ 
  geom_vline(xintercept=1,  color = "black",size =1)+
  geom_hline(yintercept=1.30102999566, linetype="dashed", color = "blue",size =1)

g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())

g = g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=18,color="black"))
g=g+theme(axis.text.x = element_text(size=18,color="black"))
g=g+xlab("Hazard Ratio  ")+ylab("Log10 p-value")
g=g+labs(title = paste0("GSE157011"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g+theme(legend.text = element_text(size = 15))+ylim(0, 4)+xlim(0.5,1.5)
print(g)



################################################
rm(list=ls())
myinf1 = "Bueno_GSE157011_GenomicEvent_iRAS_SQC.txt"
myinf2 = "Bueno_GSE157011/Clinical_info.txt"

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = as.numeric(info$os_mo)
e.surv = info$os_event
info = cbind(t.surv, e.surv, info)
tag = !is.na(info$t.surv)
info = info[tag==1,]
dim(info) #484 10

se = grep("^CAD", info$Title)
info = info[-se,]

data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp
##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(info)		## 

se = grep("uni.adj__", colnames(data))
data = data[,se]
#------------------------------------------------------ 

info=info[,c("t.surv","e.surv","age","Sex","Stage")]

# Multivariable volcano plot
table(info$Stage)
table(info$Sex)
survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  score=mytf
  scale.score = (score - mean(score))/ sd((score))
  # xx = cbind(mytf, info)
  xx = cbind(scale.score, info)
  # xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  # mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = survreg(Surv(t.surv, e.surv)~scale.score+age+Sex+Stage, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["scale.score", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~scale.score+age+Sex+Stage, xx) 
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
length(which(res$survreg.pval1<0.05))
length(which(res$coxph.pval1<0.05))

################################################################################
# volcano plot
library(ggrepel)
library(stringr)
library(dplyr)
se = grep("uni.adj__", res$name)
de = res[se,]

de$name = gsub("uni.adj__", "",res$name[se])
str_sub(de$name,start=-5, end=-4) <- "-"
# de$logp=-log10(as.numeric(de$coxph.pval1))

#reduce redundancy
df= read.csv("/mount/amos1/home/ylin/lung-prelesion/output/amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")
se1=which(de$name%in%df$rep)
se2=which(!de$name%in%df$entrezgene_accession)
de=de[c(se1,se2),]

de$logp=-log10(as.numeric(de$survreg.pval1))

de$diffexpressed <- "-"
de$diffexpressed[de$hr1>1] <- "UP"
de$diffexpressed[de$hr1<1] <- "DOWN"
de$diffexpressed[de$survreg.pval1>0.05] <- "Not significant"
de$type = de$name
se = grep("-AMP", de$name)
de$type[se] ="Amplification"
se = grep("-DEL", de$name)
de$type[se]="Deletion"
se = grep("-MUT", de$name)
de$type[se]="Mutation"
# palette  = c('#999999',"#e64B35FF","#4dbbd5ff")
shape = c(15,6,4)



top10 = tail(arrange(de,logp),5)
top10 =top10$name
# top10

de$label= "-"
# se = which(de$name %in% top10)
se = which(de$survreg.pval1<0.05)
de$label[se] =de$name[se] 
de$label[-se]=NA
de$diffexpressed
de$type

dev.off()
g=ggplot(data=de, aes(x=as.numeric(hr1), y=logp, label=name,shape=type)) + 
  geom_point(size=2) + 
  scale_shape_manual(values=shape)+
  # scale_color_manual(values=palette)+
  geom_text_repel(aes(label = label),
                  size = 4)+ 
  geom_vline(xintercept=1,  color = "black",size =1)+
  geom_hline(yintercept=1.30102999566, linetype="dashed", color = "blue",size =1)

g = g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g = g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=18,color="black"))
g=g+theme(axis.text.x = element_text(size=18,color="black"))
g=g+xlab("Hazard Ratio  ")+ylab("Log10 p-value")
g=g+labs(title = paste0("GSE157011"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())+ ylim(0, 4)+xlim(0.5,1.5)
# g=g+theme(legend.text = element_text(size = 15))
print(g)


##############################
# [6.2] cross-prediction

rm(list=ls())
myinf1 = "TCGA-LUSC_GenomicEvent_iRAS_SQC.txt"
myinf2 = "TCGA_LUSC_Freq_GenomicEvents_Gene.txt"

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


info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")

#---------------------------
par(mfrow=c(2,2))

se = grep("CDKN2A", colnames(info))
xx = info[,se]
se = which(xx$CDKN2A__MUT>0)
sam1 = row.names(xx)[se]
se = which(xx$CDKN2A__DEL>0)
sam2 = row.names(xx)[se]

tag1 = ifelse(row.names(data)%in%sam1, 1, 0)
tag2 = ifelse(row.names(data)%in%sam2, 1, 0)
xx.mut = data$CDKN2A__MUT
xx.del = data$CDKN2A__DEL

#Use mut signature
myList = list(NULL)
xx = xx.mut
myList[[1]] = xx[tag1==1 & tag2==1]
myList[[2]] = xx[tag1==1 & tag2==0]
myList[[3]] = xx[tag1==0 & tag2==1]
myList[[4]] = xx[tag1==0 & tag2==0]
names(myList) = c("11", "10", "01" ,"00")
boxplot(myList)
wilcox.test(c(myList[[3]]), myList[[4]], alternative="g")
wilcox.test(c(myList[[3]], myList[[1]]), myList[[4]], alternative="g")

library(reshape2)
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rstatix)
pal=c("#FBB4AE","#ABD9E9","#CCEBC5","#FDD586")
g_order=c("11","10","01","00")
dat$group<- factor(dat$group,levels=g_order)

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Mut Signature score")
g=g+labs(title = paste0("CDKN2A-MUT"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g +theme(axis.text.x = element_text(angle = 60,vjust = 1, hjust=1))
print(g)


myList = list(NULL)
xx = xx.del
myList[[1]] = xx[tag1==1 & tag2==1]
myList[[2]] = xx[tag1==1 & tag2==0]
myList[[3]] = xx[tag1==0 & tag2==1]
myList[[4]] = xx[tag1==0 & tag2==0]
names(myList) = c("11", "10", "01" ,"00")
# boxplot(myList)
wilcox.test(c(myList[[2]]), myList[[4]], alternative="g")
wilcox.test(c(myList[[2]], myList[[1]]), myList[[4]], alternative="g")
sapply(myList, length)

library(reshape2)
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rstatix)
pal=c("#FBB4AE","#ABD9E9","#CCEBC5","#FDD586")
g_order=c("11","10","01","00")
dat$group<- factor(dat$group,levels=g_order)
# levels(dat$group)= c(paste0(mut," (n=",table(dat$group)[1], ")"), paste0("WO (n=",table(dat$group)[2],")"))
# levels(dat$group)= c(paste0("WO"," (n=",table(dat$group)[1], ")"), paste0(mut," (n=",table(dat$group)[2], ")"))
my_comparisons <- list(c("10","00"))
stat.test <- dat %>% 
  wilcox_test(score ~ group, comparisons = my_comparisons) %>%
  add_xy_position(x = "group") %>%
  mutate(myformatted.p = paste0("p = ", ifelse(p < 1e-5, format(p, scientific = T), signif(p, digits = 1))))

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Del Signature score")
g=g+labs(title = paste0("CDKN2A-DEL"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g +theme(axis.text.x = element_text(angle = 60,vjust = 1, hjust=1))
print(g)

#---------------------------
se = grep("PIK3CA", colnames(info))
xx = info[,se]
se = which(xx$PIK3CA__MUT>0)
sam1 = row.names(xx)[se]
se = which(xx$PIK3CA__AMP>0)
sam2 = row.names(xx)[se]

tag1 = ifelse(row.names(data)%in%sam1, 1, 0)
tag2 = ifelse(row.names(data)%in%sam2, 1, 0)
xx.mut = data$PIK3CA__MUT
xx.del = data$PIK3CA__AMP


myList = list(NULL)
xx = xx.mut
myList[[1]] = xx[tag1==1 & tag2==1]
myList[[2]] = xx[tag1==1 & tag2==0]
myList[[3]] = xx[tag1==0 & tag2==1]
myList[[4]] = xx[tag1==0 & tag2==0]
names(myList) = c("11", "10", "01" ,"00")
boxplot(myList)
wilcox.test(c(myList[[3]]), myList[[4]], alternative="g")
wilcox.test(c(myList[[3]], myList[[1]]), myList[[4]], alternative="g")

library(reshape2)
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rstatix)
pal=c("#FBB4AE","#ABD9E9","#CCEBC5","#FDD586")
g_order=c("11","10","01","00")
dat$group<- factor(dat$group,levels=g_order)
my_comparisons <- list(c("01","00"))
stat.test <- dat %>% 
  wilcox_test(score ~ group, comparisons = my_comparisons) %>%
  add_xy_position(x = "group") %>%
  mutate(myformatted.p = paste0("p = ", ifelse(p < 1e-5, format(p, scientific = T), signif(p, digits = 1))))

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Mut Signature score")
g=g+labs(title = paste0("PIK3CA-MUT"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g +theme(axis.text.x = element_text(angle = 60,vjust = 1, hjust=1))
# g=g+stat_compare_means(comparisons = my_comparisons,label.y = c(1.5)) # Add pairwise comparisons p-value
print(g)

myList = list(NULL)
xx = xx.del
myList[[1]] = xx[tag1==1 & tag2==1]
myList[[2]] = xx[tag1==1 & tag2==0]
myList[[3]] = xx[tag1==0 & tag2==1]
myList[[4]] = xx[tag1==0 & tag2==0]
names(myList) = c("11", "10", "01" ,"00")
boxplot(myList)
wilcox.test(c(myList[[2]]), myList[[4]], alternative="g")
wilcox.test(c(myList[[2]], myList[[1]]), myList[[4]], alternative="g")
sapply(myList, length)

library(reshape2)
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rstatix)
pal=c("#FBB4AE","#ABD9E9","#CCEBC5","#FDD586")
g_order=c("11","10","01","00")
dat$group<- factor(dat$group,levels=g_order)
my_comparisons <- list(c("10","00"))
stat.test <- dat %>% 
  wilcox_test(score ~ group, comparisons = my_comparisons) %>%
  add_xy_position(x = "group") %>%
  mutate(myformatted.p = paste0("p = ", ifelse(p < 1e-5, format(p, scientific = T), signif(p, digits = 1))))

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Del Signature score")
g=g+labs(title = paste0("PIK3CA-AMP"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g +theme(axis.text.x = element_text(angle = 60,vjust = 1, hjust=1))
print(g)

#-------------------
# se = grep("PIK3CA", colnames(info))
xx = info
se = which(info$PIK3CA__MUT>0)
sam1 = row.names(xx)[se]
se = which(info$SOX2__AMP>0) #change pik3ca to sox2
sam2 = row.names(xx)[se]

tag1 = ifelse(row.names(data)%in%sam1, 1, 0)
tag2 = ifelse(row.names(data)%in%sam2, 1, 0)
xx.mut = data$PIK3CA__MUT
xx.del = data$SOX2__AMP


myList = list(NULL)
xx = xx.mut
myList[[1]] = xx[tag1==1 & tag2==1]
myList[[2]] = xx[tag1==1 & tag2==0]
myList[[3]] = xx[tag1==0 & tag2==1]
myList[[4]] = xx[tag1==0 & tag2==0]
names(myList) = c("11", "10", "01" ,"00")
boxplot(myList)
wilcox.test(c(myList[[3]]), myList[[4]], alternative="g")
wilcox.test(c(myList[[3]], myList[[1]]), myList[[4]], alternative="g")

library(reshape2)
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rstatix)
pal=c("#FBB4AE","#ABD9E9","#CCEBC5","#FDD586")
g_order=c("11","10","01","00")
dat$group<- factor(dat$group,levels=g_order)
my_comparisons <- list(c("01","00"))
stat.test <- dat %>% 
  wilcox_test(score ~ group, comparisons = my_comparisons) %>%
  add_xy_position(x = "group") %>%
  mutate(myformatted.p = paste0("p = ", ifelse(p < 1e-5, format(p, scientific = T), signif(p, digits = 1))))

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Mut Signature score")
g=g+labs(title = paste0("PIK3CA-MUT"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g +theme(axis.text.x = element_text(angle = 60,vjust = 1, hjust=1))
# g=g+stat_compare_means(comparisons = my_comparisons,label.y = c(1.5)) # Add pairwise comparisons p-value
print(g)

myList = list(NULL)
xx = xx.del
myList[[1]] = xx[tag1==1 & tag2==1]
myList[[2]] = xx[tag1==1 & tag2==0]
myList[[3]] = xx[tag1==0 & tag2==1]
myList[[4]] = xx[tag1==0 & tag2==0]
names(myList) = c("11", "10", "01" ,"00")
boxplot(myList)
wilcox.test(c(myList[[2]]), myList[[4]], alternative="g")
wilcox.test(c(myList[[2]], myList[[1]]), myList[[4]], alternative="g")
sapply(myList, length)

library(reshape2)
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rstatix)
pal=c("#FBB4AE","#ABD9E9","#CCEBC5","#FDD586")
g_order=c("11","10","01","00")
dat$group<- factor(dat$group,levels=g_order)
my_comparisons <- list(c("10","00"))
stat.test <- dat %>% 
  wilcox_test(score ~ group, comparisons = my_comparisons) %>%
  add_xy_position(x = "group") %>%
  mutate(myformatted.p = paste0("p = ", ifelse(p < 1e-5, format(p, scientific = T), signif(p, digits = 1))))

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+xlab("Patient Group")+ylab("Del Signature score")
g=g+labs(title = paste0("PIK3CA-AMP"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
g=g +theme(axis.text.x = element_text(angle = 60,vjust = 1, hjust=1))
print(g)

####################
# version2 crosspredict

rm(list=ls())
myinf1 = "TCGA-LUSC_GenomicEvent_iRAS_SQC.txt"
myinf2 = "TCGA_LUSC_Freq_GenomicEvents_Gene.txt"

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


info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")

#---------------------------
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rstatix)



se = grep("CDKN2A", colnames(info))
xx = info[,se]
se = which(xx$CDKN2A__MUT>0)
sam1 = row.names(xx)[se]
se = which(xx$CDKN2A__DEL>0)
sam2 = row.names(xx)[se]

tag1 = ifelse(row.names(data)%in%sam1, 1, 0)
tag2 = ifelse(row.names(data)%in%sam2, 1, 0)
xx.mut = data$CDKN2A__MUT
xx.del = data$CDKN2A__DEL


#Use mut signature
myList = list(NULL)
xx = xx.mut
myList[[1]] = xx[tag1==1]
myList[[2]] = xx[tag1==0]

boxplot(myList)
wilcox.test(c(myList[[1]]), myList[[2]], alternative="g") # p-value = 1.304e-07;#  p-value = 0.0004276
names(myList) = c("MUT+", "MUT-")
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)
pal=c("#FBB4AE","#ABD9E9")

g_order=c("MUT+","MUT-")
dat$group<- factor(dat$group,levels=g_order)

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+ylab("Mut Signature score")
g=g+labs(title = paste0("CDKN2A"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
p=g


myList = list(NULL)
myList[[1]] = xx[tag2==1]
myList[[2]] = xx[tag2==0]
pal =c("#CCEBC5","#FDD586")
names(myList)= c("DEL+","DEL-")
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)

g_order=c("DEL+","DEL-")
dat$group<- factor(dat$group,levels=g_order)

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+ylab("Mut Signature score")
g=g+labs(title = paste0("CDKN2A"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))



#####################
xx = info
se = which(info$PIK3CA__MUT>0)
sam1 = row.names(xx)[se]
se = which(info$PIK3CA__AMP>0) #change pik3ca to sox2
sam2 = row.names(xx)[se]

tag1 = ifelse(row.names(data)%in%sam1, 1, 0)
tag2 = ifelse(row.names(data)%in%sam2, 1, 0)
xx.mut = data$PIK3CA__MUT
# xx.del = data$SOX2__AMP

myList = list(NULL)
xx = xx.mut
myList[[1]] = xx[tag1==1]
myList[[2]] = xx[tag1==0]
wilcox.test(myList[[1]],myList[[2]],alternative = "g") #3.923e-05
names(myList) = c("MUT+", "MUT-")
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)
pal=c("#FBB4AE","#ABD9E9")

g_order=c("MUT+","MUT-")
dat$group<- factor(dat$group,levels=g_order)

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+ylab("Mut Signature score")
g=g+labs(title = paste0("PIK3CA"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))
p=g
print(p)

myList=NULL
myList[[1]] = xx[tag2==1]
myList[[2]] = xx[tag2==0]
wilcox.test(myList[[1]],myList[[2]],alternative = "g")
pal =c("#CCEBC5","#FDD586")
names(myList)= c("AMP+","AMP-")
dat = melt(myList)
dat= as.data.frame(dat)
colnames(dat) = c("score","group")
dat$score = (dat$score - mean(dat$score,na.rm=TRUE))/sd(dat$score,na.rm=TRUE)

g_order=c("AMP+","AMP-")
dat$group<- factor(dat$group,levels=g_order)

g= ggplot(dat, aes(x=group, y=score,fill=group) )+geom_boxplot(color="black",outlier.shape = NA)
g=g+scale_fill_manual(values = pal)
g=g+theme(panel.grid =element_blank())+theme(panel.background = element_blank())
g <- g+theme(panel.border = element_rect(size=0.5,fill=NA))
g=g+theme(axis.text.y = element_text(size=20,color="black"))
g=g+theme(axis.text.x = element_text(size=17,color="black"))
g=g+ylab("Mut Signature score")
g=g+labs(title = paste0("PIK3CA"))
g=g+theme(title=element_text(face="bold",size=20))
g=g+theme(plot.title = element_text(hjust = 0.5))
g=g+theme(aspect.ratio = 1)
g=g+theme(legend.background = element_rect(color=NA,fill=NA))
g=g+theme(legend.title=element_blank())
g=g + theme(legend.position="none")
g=g+theme(legend.text = element_text(size = 15))

print(g)
###################################
# TCGA-volcano plot
myinf1  = "TCGA_12_FireHose_maf_SomaticMutation_count.rda"
myinf3 = "LUSC_Clincial_info.txt"
info = read.table(myinf3, sep="\t", header=T, row.names=1, quote="") 

se = c("vital_status","days_to_death","days_to_last_followup")
info = info[,se]
xx=ifelse(!is.na(info$days_to_death),info$days_to_death, info$days_to_last_followup)
t.surv = as.numeric(xx)
e.surv = ifelse(info[, "vital_status"]=="dead", 1, 0)
info = cbind(t.surv, e.surv, info)
info = info[!is.na(info$t.surv), ]

myinf1 = "TCGA_LUSC_GenomicEvent_iRAS.txt"
data = read.table(myinf1, sep="\t", header=T, row.names=1)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp



##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(info)		## 490

library(survival)
#pathway signatures
survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  # cat("\r", k)
  # k=1
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)$table
  # mycox
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
de = res[se,]
#
# write.csv(de, "/mount/amos1/home/ylin/lung-prelesion/output/lusc/lusc_survreg.txt", quote=F)
de$pvalue=-log10(as.numeric(de$coxph.pval1))
de$loghr= log2(as.numeric(de$hr1))
de$name = gsub("uni.adj__", "",res$name[se])

#--------
str_sub(de$name,start=-5, end=-4) <- "-" 
tmp=de$name
str_sub(tmp,start=-5, end=-4) <- "-" 

df= read.csv("/mount/amos1/home/ylin/lung-prelesion/output/amp-manual-group.csv")
df$entrezgene_accession=paste0(df$entrezgene_accession,"-AMP")
df$rep=paste0(df$rep,"-AMP")
se1=which(tmp%in%df$rep)
se2=which(!tmp%in%df$entrezgene_accession)
tmp = tmp[c(se1,se2)]
de2= de[which(de$name %in%tmp),]
length(which(de2$coxph.pval1<=0.05)) #11 out of 34 significant

de$diffexpressed <- "-"
de$diffexpressed[de$coxph.pval1>0.05] <- "Not significant"
de$diffexpressed[de$loghr >  median(de$loghr) & de$coxph.pval1 < 0.05] <- "UP"
de$diffexpressed[de$loghr < median(de$loghr) & de$coxph.pval1 < 0.05] <- "DOWN"
de$type = de$name
se = grep("_AMP", de$name)
de$type[se] ="Amplification"
se = grep("_DEL", de$name)
de$type[se]="Deletion"
se = grep("_MUT", de$name)
de$type[se]="Mutation"
palette  = c("#4dbbd5ff",'#999999',"#e64B35FF")
shape = c(15,6,4)

library(ggrepel)

ggplot(data=de, aes(x=as.numeric(loghr), y=pvalue, label=name, col=diffexpressed,shape=type)) + 
  geom_point(size=2) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_shape_manual(values=shape)+
  scale_color_manual(values=palette)+
  geom_text_repel() +scale_y_continuous(limits=c(0, 3.5),breaks = c(1,2,3))


#mutation
#TODO change data range
myinf1 = "TCGA_LUSC_Freq_SomaticMutation.txt"
myinf2 =  "TCGA_LUSC_Freq_GeneCNV.txt"
myinf3 = "LUSC_Clincial_info.txt"

info = read.table(myinf3, sep="\t", header=T, row.names=1, quote="") 
se = c("vital_status","days_to_death","days_to_last_followup")
info = info[,se]
xx=ifelse(!is.na(info$days_to_death),info$days_to_death, info$days_to_last_followup)
t.surv = as.numeric(xx)
e.surv = ifelse(info[, "vital_status"]=="dead", 1, 0)
info = cbind(t.surv, e.surv, info)
info = info[!is.na(info$t.surv), ]

mut = read.table(myinf1, sep="\t", header=T, row.names=1)
cnv = read.table(myinf2, sep="\t", header=T, row.names=1)
comxx = intersect(row.names(mut), row.names(cnv))
mut = mut[comxx,]
cnv = cnv[comxx,]

data=cbind(mut,cnv)
# tmp=apply(data,2,range)
##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(info)		## 

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

#MUT
de=res
de$pvalue=-log10(as.numeric(de$coxph.pval1))
de$loghr= log2(as.numeric(de$hr1))

de$diffexpressed <- "-"
de$diffexpressed[de$coxph.pval1>0.05] <- "Not significant"
de$diffexpressed[de$loghr >  median(de$loghr) & de$coxph.pval1 < 0.05] <- "UP"
de$diffexpressed[de$loghr < median(de$loghr) & de$coxph.pval1 < 0.05] <- "DOWN"
de$type = de$name
se = grep("_AMP", de$name)
de$type[se] ="Amplification"
se = grep("_DEL", de$name)
de$type[se]="Deletion"
se = grep("_MUT", de$name)
de$type[se]="Mutation"
palette  = c('#999999',"#e64B35FF","#4dbbd5ff")
shape = c(15,6,4)

library(ggrepel)

p= ggplot(data=de, aes(x=as.numeric(loghr), y=pvalue, label=name, col=diffexpressed,shape=type)) + 
  geom_point(size=2) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_shape_manual(values=shape)+
  scale_color_manual(values=palette)+
  geom_text_repel() +scale_y_continuous(limits=c(0, 3.5),breaks = c(1,2,3))
print(p)
