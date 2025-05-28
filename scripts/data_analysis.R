install.packages("qqman")
library(qqman)

Sample QC

imiss<-read.table("sampleQC.imiss",header=T)
head(imiss)
summary(imiss$F_MISS)
hist(imiss$F_MISS)

het<-read.table("sampleQC.het",header=T)
head(het)
summary(het$F)
hist(het$F)

#Identify samples with call rate < 0.98 or inbreeding coefficient < -0.04
samplepb<-rbind(subset(imiss,imiss$F_MISS>=0.02,select=c("FID","IID")),subset(het,het$F<(-0.04),select=c("FID","IID")))
samplepb

write.table(samplepb,file="badsamples.txt",col.names=F,row.names=F,sep="\t",quote=F)

SNP QC

lmiss<-read.table("SNPQC.lmiss",header=T)
head(lmiss)
summary(lmiss$F_MISS)
hist(lmiss$F_MISS)

testmiss<-read.table("SNPQC.missing",header=T)
head(testmiss)
summary(-log10(testmiss$P))
hist(-log10(testmiss$P))

hwe<-read.table("SNPQC.hwe",header=T)
head(hwe)
hwe_ctrl<-subset(hwe,hwe$TEST=="UNAFF")
summary(-log10(hwe_ctrl$P))
hist(-log10(hwe_ctrl$P))

#Identify SNPs with call rate < 0.99 or case/control missing p-val< 0.00001 or hwe p-val < 0.000001
SNPpb<-rbind(subset(lmiss,lmiss$F_MISS>=0.01,select=c("SNP")),subset(testmiss,testmiss$P<0.00001,select=c("SNP")),subset(hwe_ctrl,hwe_ctrl$P<0.000001,select=c("SNP")))
nrow(SNPpb)

write.table(SNPpb,file="badSNPs.txt",col.names=F,row.names=F,quote=F)

###Relatedness######

related<-read.table("Related.genome",header=T)
print(related)
summary(related$PI_HAT)

#Select one sample from each duplicated pair
duplicate<-related[,c("FID1","IID1")]
duplicate

write.table(duplicate,file="dup_samples.txt",col.names=F,row.names=F,quote=F)

###Association testing######
##Post-QC

#Model without covariates
assoc.nocov<-read.table("HCVnocov.assoc.logistic",header=T)
head(assoc.nocov)

qq(assoc.nocov$P)
lambda.nocov<-median(assoc.nocov$STAT^2,na.rm=T)/0.456
lambda.nocov

###Population structure######

pc<-read.table("PCA.eigenvec",header=T)
head(pc)

fam<-read.table("HAfullQC.fam",header=F)
names(fam)<-c("FID","IID","FA","MO","SEX","AFS")

pc<-merge(pc,fam,by=c("FID","IID"))

plot(pc$PC1~pc$PC2,pch=19)
points(pc[pc$AFS==2,]$PC1~pc[pc$AFS==2,]$PC2,pch=19,col="red")

plot(pc$PC1~pc$PC2, lab)
points(pc[pc$AFS==2,]$PC1~pc[pc$AFS==2,]$PC2,col="red")

text(pc$PC1, pc$PC2,
     labels = pc$IID,
     pos = 4,
     cex = 0.8,
     col = "red")

#PC

assoc.PC<-read.table("HAPC.assoc.logistic",header=T)
head(assoc.PC)
assoc.PC<-subset(assoc.PC,assoc.PC$TEST=="ADD")

qq(assoc.PC$P)
lambda.PC<-median(assoc.PC$STAT^2,na.rm=T)/0.456
lambda.PC

plot =manhattan(assoc.PC[!is.na(assoc.PC$P),])
plot

manhattan(assoc.PC, main = "Manhattan Plot", ylim = c(0, 10), cex = 0.6, cex.axis = 1,
col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F)

#top snps
top = subset(assoc.PC,assoc.PC$P<0.00001)
top

#adjusted p-value
subset(assoc.PC,assoc.PC$P<0.00001)$P*373397
subset(assoc.PC,assoc.PC$P<0.00001)$P*92397

#output chr19

chr19<-subset(assoc.PC,assoc.PC$CHR==19)

write.table(chr19,file="chr19.txt",col.names=T,row.names=F,sep="\t",quote=F)


# Model adjusted on 10 PCs
assoc.PC10<-read.table("HAPC10.assoc.logistic",header=T)
assoc.PC10<-subset(assoc.PC10,assoc.PC10$TEST=="ADD")

qq(assoc.PC10$P)
lambda.PC10<-median(assoc.PC10$STAT^2,na.rm=T)/0.456
lambda.PC10
manhattan(assoc.PC10[!is.na(assoc.PC10$P),])

#top snps
subset(assoc.PC10,assoc.PC10$P<0.00001)

#adjusted p-value
subset(assoc.PC10,assoc.PC10$P<0.00001)$P*373397
subset(assoc.PC10,assoc.PC10$P<0.00001)$P*92397

#output chr19

chr19<-subset(assoc.PC,assoc.PC$CHR==19)

write.table(chr19,file="chr19.txt",col.names=T,row.names=F,sep="\t",quote=F)
