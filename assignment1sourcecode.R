# Loading the needed libraries
library(affy)
library(limma)
library(magrittr)
library(hgu133plus2.db)

# Set directory path, read CEL files and normalise with rma
setwd("~/lsm3241/GSE51395/celfiles")
celfiles=ReadAffy()
rmanormalized=rma(celfiles)

# Delineating the experimental conditions
diseasestatus=c(rep("Nontreated",10), 
                rep("Knockdown",10))
modelmatrix=model.matrix(~diseasestatus+0)
colnames(modelmatrix)=c("Knockdown","Nontreated")
contrastmatrix=makeContrasts(Nontreated-Knockdown, levels=modelmatrix) #Since knockdown effect is being studied

# Fitting linear model and identifying significant DGE
fit1=lmFit(rmanormalized, modelmatrix)
contrastfit=contrasts.fit(fit1,contrastmatrix) 
contrastfit=eBayes(contrastfit) 
significanthits=topTable(contrastfit,p.value = 0.05, n=Inf) 

probeMapping=select(hgu133plus2.db,columns=c("SYMBOL",”GENENAME”),keys = rownames(significanthits))
