
# Following code are generated to compute patient-wise signature score, major sections include:
# [1] Identify genomic event highly occured in LUSC
# [2] Define weighted profiles
# [3]Obtain sample-specific enrichment score from BASE algorithm
# [4] Downstream analysis: Association with survival

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [1] Identify genomic event highly occured in LUSC
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.1] Somatic mutations
rm(list=ls())
myinf1  = "TCGA_12_FireHose_maf_SomaticMutation_count.rda" #TCGA
myinf2 = "cancer_gene_census.tsv" #COSMIC
myoutf1 = "TCGA_LUSC_Freq_SomaticMutation.txt"

fra.thr = 0.1	## in more than 10% samples
info = read.table(myinf2, sep="\t", header=T, row.names=1)
cosmic = row.names(info)
length(cosmic)

load(file= myinf1)
data = mydata
se = grep("LUSC__", colnames(data))
data = data[, se]
dim(data)
colnames(data) = gsub("LUSC__", "", colnames(data))


xx = apply(data>0, 1, sum)
se= which(xx>=ncol(data)*fra.thr)
length(se)
data = data[se, ]
dim(data)

se = which(row.names(data)%in%cosmic)
length(se)
data = data[se, ]
dim(data)

res = t(data)
colnames(res) = paste(colnames(res), "__MUT", sep="")
dim(res)
write.table(res, myoutf1, sep="\t", quote=F)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.2] Identify Gene amp/del from TCGA
rm(list=ls())
myinf1 = "LUSC_CNV_Symbol.rda"
myinf2 = "cancer_gene_census.tsv"


myoutf1 = "TCGA_LUSC_Freq_GeneCNV.txt"

fra.thr = 0.1	## in more than 10% samples
amp.thr = log2(2.8/2)
del.thr = log(1.2/2)


info = read.table(myinf2, sep="\t", header=T, row.names=1)
cosmic = row.names(info)
length(cosmic)
cosmic = c(cosmic, "SPOP")

load(myinf1)
data = mydata
se = which(row.names(data)%in%cosmic)
data = data[se,]
dim(data)
xx = as.numeric(substr(colnames(data), 14, 15))
se = which(xx==1)
data = data[,se]
dim(data)
colnames(data) = substr(colnames(data), 1, 12)


xx = apply(data>amp.thr, 1, sum)
se = which(xx>ncol(data)*fra.thr)
length(se)
dat1 = data[se,]

xx = apply(data<del.thr, 1, sum)
se = which(xx>ncol(data)*fra.thr)
length(se)
dat2 = data[se,]

dat1 = t(dat1)
for(k in 1:ncol(dat1))
{
  dat1[,k] = ifelse(dat1[,k]>amp.thr, 1, 0)
}
tmp = apply(dat1, 2, sum)
dat1 = dat1[, order(tmp, decreasing=T)]

apply(dat1, 2, sum)  

xx= cor(dat1)
for(k in 1:ncol(xx))
{
  xx[k,k] = 0
}


## reduce redundancy
## PMS2--> RAC1
## HOXA11 --> HOXA9 HOXA13

mygen = c("SOX2", "PIK3CA", "ETV5", "EIF4A2", "BCL6", "LPP", "TFRC", "MLF1",
          "GMPS", "WWTR1", "TERT", "IL7R", "LIFR", "FOXL2", "GATA2", 
          "WHSC1L1", "FGFR1", "CBLB", "MYC", "TFG", "CCND1", "NDRG1", "RAD21",
          "EXT1", "REL", "XPO1", "CCNE1", "CEBPA", "COX6C", "AKT2",
          "KRAS", "ARNT", "EGFR", "FCGR2B", "SDHC", "CCND2", "MUC1", "PBX1",
          "ASXL1", "HEY1", "HOOK3", "ZNF384", "TPM3")

se = which(colnames(dat1)%in%mygen)
dat1 = dat1[, se]
colnames(dat1) = paste(colnames(dat1), "__AMP", sep="")


dat2 = t(dat2)
for(k in 1:ncol(dat2))
{
  dat2[,k] = ifelse(dat2[,k]<del.thr, 1, 0)
}
tmp = apply(dat2, 2, sum)
dat2 = dat2[, order(tmp, decreasing=T)]

apply(dat2, 2, sum)

xx= cor(dat2)
for(k in 1:ncol(xx))
{
  xx[k,k] = 0
}
apply(xx, 1, max)

mygen = c("CDKN2A", "FOXP1", "PCM1")
se = which(colnames(dat2)%in%mygen)
dat2 = dat2[, se]

colnames(dat2) = paste(colnames(dat2), "__DEL", sep="")

res = cbind(dat1, dat2)
dim(res)
write.table(res, myoutf1, sep="\t", quote=F)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.3] Integate the samples from 1.1 and 1.2
rm(list=ls())
myinf1 = "TCGA_LUSC_Freq_SomaticMutation.txt"
myinf2 = "TCGA_LUSC_Freq_GeneCNV.txt"
myoutf1 = "TCGA_LUSC_Freq_GenomicEvents_Gene.txt"

dat1 = read.table(myinf1, sep="\t", row.names=1, header=T, stringsAsFactors=F)
dat2 = read.table(myinf2, sep="\t", row.names=1, header=T, stringsAsFactors=F)

comxx = union(row.names(dat1), row.names(dat2))
data = matrix(NA, length(comxx), ncol(dat1)+ncol(dat2))
row.names(data) = comxx
colnames(data) = c(colnames(dat1), colnames(dat2))
data = as.data.frame(data)
data[row.names(dat1), 1:ncol(dat1)] = dat1
data[row.names(dat2), (1+ncol(dat1)):ncol(data)] = dat2

dim(data)
write.table(data, myoutf1, sep="\t", quote=F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [2] Define weighted profiles
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list=ls())
myinf1 = "LUSC_RNAseqv2_Tumor_Symbol.rda"
myinf2 = "LUSC_Clincial_info.txt"
myinf3 = "TCGA_LUSC_Freq_GenomicEvents_Gene.txt"
myoutf1 = "TCGA-LUNG-LUSC-GenomicEvent-Profile_4base_Gene.txt"

#----------------------------
soma = read.table(myinf3, sep="\t", header=T, row.names=1, quote="")
for(k in 1:ncol(soma))
{
  soma[,k] = ifelse(soma[,k]>0, 1, 0)
}
apply(soma, 2, sum, na.rm=T)

#---------------------------
load(file= myinf1)
data = mydata	
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
se = which(row.names(soma)%in%colnames(data))
soma = soma[se,]
xx = apply(soma, 2, sum, na.rm=T)
xx


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

#----------------------
raw.soma = soma
comxx = intersect(row.names(raw.soma), colnames(raw.data))
soma = raw.soma[comxx,]
data = raw.data[, comxx]
info = raw.info[comxx,]

nn = ncol(soma)
for(k in 1:nrow(data))
{
  cat("\r",  k)
  mytf = as.numeric(data[k,])
  xx = cbind(mytf, soma)
  mylog <- lm(mytf ~ ., data = xx)
  mylog = summary(mylog)
  pval =  mylog[["coefficients"]][2:(nn+1), 4]
  beta =  mylog[["coefficients"]][2:(nn+1), 1]
  xx = -log10(pval)
  xx = ifelse(beta>0, xx, -xx)
  p3[k,] = xx
  
  xx = cbind(mytf, soma, info)
  mylog <- lm(mytf ~ ., data = xx)
  mylog = summary(mylog)
  pval =  mylog[["coefficients"]][2:(nn+1), 4]
  beta =  mylog[["coefficients"]][2:(nn+1), 1]
  xx = -log10(pval)
  xx = ifelse(beta>0, xx, -xx)
  p4[k,] = xx
}

#--------------------------
res = p1
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp[is.na(tmp)] = 0
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
  tmp[is.na(tmp)] = 0
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


res = p3
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp[is.na(tmp)] = 0
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
colnames(myres) = paste("mul.noj", colnames(myres), sep="__")
minv = min(myres)
maxv = max(myres)
myres = (myres-minv)/(maxv-minv)
prof3 = myres

res = p4
res1 = res2 = res
for(k in 1:ncol(res))
{
  tmp = res[,k]
  tmp[is.na(tmp)] = 0
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
colnames(myres) = paste("mul.adj", colnames(myres), sep="__")
minv = min(myres)
maxv = max(myres)
myres = (myres-minv)/(maxv-minv)
prof4 = myres

#------------------------------------------------------------
profile = cbind(prof1, prof2, prof3, prof4)
se1 = grep("\\.up", colnames(profile))
se2 = grep("\\.dn", colnames(profile))
profile = profile[, c(se1, se2)]
dim(profile)	## 20500
write.table(profile, myoutf1, sep="\t", quote=F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [3]Obtain sample-specific enrichment score from BASE algorithm
rm(list=ls())
myinf1 = "Input-Symbol_expr.txt"
myinf2 = "TCGA-LUNG-LUSC-GenomicEvent-Profile_4base_Gene.txt"
myoutf = "sample_GenomicEvent_iRAS_SQC.txt"

mywt = read.table(myinf2, sep="\t", header=T, quote="", row.names=1)
source("~/WorSpa/system/myRprogram/base5.R")
data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1, stringsAsFactors=F)

##
reg = mywt
xx = base5(data, reg, perm=1000, myoutf, median.norm=T)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# [4] Downstream analysis example: Association with survival
rm(list=ls())
myinf1 = "sample_GenomicEvent_iRAS_SQC.txt" # from section 3
myinf2 = "sample_Clinical_info.txt"


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

