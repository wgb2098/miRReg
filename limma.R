rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
library(limma)
getwd()
exp<-read.table("exp_matrix_ges35776.txt", header=T,sep="\t")
names(exp)
d <- exp[, 2:25]
rownames(d) <- exp[, 1]
# Pheno data file
pheno <- read.table("group.txt",header = F)
colnames(pheno) <- c("sample","group")

Group<-factor(pheno$group,levels=levels(pheno$group))

design<-model.matrix(~0+Group)

y <- voom(d,design,plot=TRUE)
colnames(design)
fit <-lmFit(y,design)
##Designing Contrast Matrix for group Differentiation
cont.wt<-makeContrasts(Groupcontrol-Grouptof,levels=design)
fit2 <-contrasts.fit(fit,cont.wt)
fit3<-eBayes(fit2)
DE<-topTable(fit3,number = 'all')
write.csv(DE, "limma_genes.csv")
selected = rownames(DE[DE$P.Value <0.05,])
esetSel = y[selected, ]
heatmap(esetSel@.Data[[1]])
# p<- DE[DE$P.Value <1,]$P.Value
# esetSela <- y[]
# p.adjust(c(0.001,0.02,0.03),"BH")

