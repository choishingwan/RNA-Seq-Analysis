#The analysis script without the options. Can turn this into a proper script later on

library(DESeq2)
library(sva)
library(wgcna)
library(affy)
library(affycoretools) 
library(GEOquery)

`%notin%` <- function(x,y) !(x %in% y) 

#####################################
#									#
#									#
#			DESeq2 Analysis			#
#									#
#									#
#####################################
#The fold change is inverted for some reason, invert it so that it is case / control instead

inputFile = "/home/sam/workspace/schizophrenia_mouse/counts/totalCount.csv"
samples = "/home/sam/workspace/schizophrenia_mouse/replication/SampleDE"
clean = TRUE
mgiName = "/home/sam/workspace/schizophrenia_mouse/replication/ensemble2mgi.csv"

if(file.access(inputFile) ==-1){
	stop(sprintf("Count table (%s) does not exist", inputFile));
} else{
	inputTable = read.csv(inputFile, row.names=1, header=T);
}
print("count table obtained")
if(file.access(samples) ==-1){
	stop(sprintf("Sample information file (%s) does not exist", samples));
} else{
	sampleInfo = read.table(samples, row.names=1);
	#We will have to make sure the sample names are correct
	row.names(sampleInfo) = make.names(row.names(sampleInfo))
	inputTable= inputTable[,row.names(sampleInfo)]
}


minSample = min(table(sampleInfo[,1]))
sprintf("sample information obtained with minimum %i samples per group, count table processed", minSample)

if(clean){
	print(paste("Clean table. nrow before cleaning: ", nrow(inputTable), sep=""));
	inputTable = inputTable[row.names(inputTable) %notin% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"),]
	print(paste("Table cleansed. nrow after cleaning: ", nrow(inputTable), sep=""));
}

colData = data.frame(row.names=row.names(sampleInfo), condition=sampleInfo[,1])
dds <- DESeqDataSetFromMatrix(countData = inputTable, colData = colData, design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,levels=unique(sampleInfo[,1]))
colData(dds)$condition <- relevel(colData(dds)$condition,  as.character(sampleInfo[6,1]))

if(minSample <7){
	dds <- DESeq(dds, minReplicatesForReplace=minSample)
}else{
	dds <- DESeq(dds)
}
res=results(dds)
res = res[-grep("ERCC", row.names(res)),]
res$padj = p.adjust(res$pvalue, "BH")
res$bon = p.adjust(res$pvalue, "bonferroni")
use <- res$baseMean >= 10 & !is.na(res$pvalue)
resFilt=as.data.frame(res[use,])

write.csv(resFilt, "DESeq2.result.csv", quote=F)
mgi = read.csv(mgiName, row.names=1) #read this just in case we need it


#####################################
#									#
#									#
#			ComBat Analysis			#
#									#
#									#
#####################################
microArrayPath = "/home/sam/workspace/schizophrenia_mouse/microarray/" #The location of the CEL files

mydata <-  ReadAffy(celfile.path=microArrayPath)
pheno = data.frame(row.names=row.names(pData(mydata)), timepoint= c(rep("GD9.5",5), rep("GD11.5",4),rep("GD13.5",6), "GD9.5"), platform="1", condition="Control")
eset <- rma(mydata) #Normalization

gpl=getGEO("GPL1261", destdir=".")
id2Name=Table(gpl)[,c(1,11)]
id2NameMatrix=matrix(NA, nrow(id2Name), 2)
for(i in 1:nrow(id2Name)){
	id2NameMatrix[i,1] = as.character(id2Name[i,1])
	split = unlist(strsplit(as.character(id2Name[i,2]), split="//")[[1]])
	if(length(split) == 1){
		id2NameMatrix[i,2] =  as.character(split[1]);
	}
	else{
		id2NameMatrix[i,2] = as.character(split[2]);
	}
}
colnames(id2NameMatrix) = c("ID", "Name")
id2NameMatrix = as.data.frame(id2NameMatrix)
micro.info = merge(exprs(eset), id2NameMatrix, by.x="row.names", by.y="ID")

vsd <- varianceStabilizingTransformation(dds)
#Need to perform filtering to filter genes that doesn't pass through the RNA QC
vsd.mgi = merge(assay(vsd), mgi, by="row.names")
vsd.mgi.filt=merge(vsd.mgi, resFilt, by.x="Row.names", by.y="row.names")
vsd.mgi.filt = vsd.mgi.filt[,1:13]
data=merge(micro.info, vsd.mgi.filt, by.x="Name", by.y="MGI.symbol")
data.num=merge(micro.info, vsd.mgi, by.x="Name", by.y="MGI.symbol")
phenoRNA = data.frame(row.names=colnames(assay(vsd)), condition=c(rep("Case", 5), rep("Control", 5)), platform="2", timepoint="GD9.5")
phenotype = rbind(pheno, phenoRNA)
data = data[,-c(2,19,30)]
data.num = data.num[,-c(2,19,30)]
data.unique = data[!duplicated(data$Name),]
data.num.unique = data.num[!duplicated(data.num$Name),]
row.names(data.unique) = data.unique$Name
data.unique = data.unique[,-1]
exclude=row.names(subset(phenotype, phenotype$timepoint!="GD9.5"))
data.input = data.unique[, !(colnames(data.unique) %in% exclude)]
pheno.input = phenotype[!(row.names(phenotype) %in% exclude),]
data.input = data.matrix(data.input)
batch = pheno.input$platform
mod = model.matrix(~as.factor(condition), data=pheno.input)
mod0 = model.matrix(~1,data=pheno.input)
combat_edata = ComBat(dat=data.input, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH",n=nrow(data.num.unique))
bValuesComBat=p.adjust(pValuesComBat, method="bonferroni",n=nrow(data.num.unique))
combatRes = data.frame(row.names=row.names(data.unique), bon=bValuesComBat, pvalue=pValuesComBat)
write.csv(combatRes, "Combat.result.csv", quote=F)


#####################################
#									#
#									#
#		Overlap Analysis			#
#									#
#									#
#####################################

resFilt.mgi = merge(resFilt, mgi, by="row.names");
row.names(resFilt.mgi) = resFilt.mgi$Row.names
resFilt.mgi = resFilt.mgi[,-1]
overlap = merge(resFilt.mgi, combatRes, by.x="MGI.symbol", by.y="row.names")
write.csv(overlap, "Overlap.result.csv", quote=F)
