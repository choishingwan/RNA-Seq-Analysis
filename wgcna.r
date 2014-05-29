#!/usr/bin/env Rscript
#Check if the optparse library is installed, required to parse arguments
`%notin%` <- function(x,y) !(x %in% y) 

if(suppressPackageStartupMessages(require("optparse"))){
	print("optparse is loaded correctly");
} else {
	print ("trying to install optparse")
	install.packages("optparse",repos="http://cran.rstudio.com/");
	if(suppressPackageStartupMessages(require("optparse"))){
		print ("optparse installed and loaded");
	} else{
		stop("could not install optparse");
	}
}



options(stringsAsFactors = FALSE);

#Perform the argument checking before spending time in loading DESeq2
option_list <- list(
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
	make_option(c("-q", "--quietly"), action="store_false",	dest="verbose", help="Print little output"),
	make_option(c("-f", "--file"), type="character",  help="The input count table", dest="inputFile"),
	make_option(c("-p", "--perm"), type="integer", default=1000, help="Number of permutation to perform [default \"%default\"]", metavar="number"),
	make_option(c("-t", "--thread"), type="integer", default=10, help="Number of thread to use [default \"%default\"]", metavar="number"),
	make_option(c("-s", "--softpower"), type="integer", default=0, help="Selected softpower. if 0, will automatically calculate the softpower and select one at the first peak whereas an input of -1 will indicates the use of estimated power from WGCNA [default \"%default\"]", metavar="number"),
	make_option(c("-d", "--de"), type="character", help="csv containing differential expression information. Assuming from DESeq2", metavar="file"),
	make_option(c("-m", "--mgi"), type="character", help="Table to transform count table identifier to MGI id in tab deliminted format with headers start with #. Format should be <Identifier> <MGI>. If it is not a one to one mapping, then the programme will randomly select one. If this file isn't provided, will assumn the identifier to be MGI id", metavar="file"),
	make_option(c("-o", "--out"), type="character", help="Output prefix [default \"count table name\"]", metavar="output")
)

 
opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$inputFile)){
	stop("Require count table input");
}
if(is.null(opt$de)){
	stop("Require differential expression information");
}
if(is.null(opt$out)){
	opt$out = opt$inputFile;
}
if(file.access(opt$inputFile) ==-1){
	stop(sprintf("Count table (%s) does not exist", opt$inputFile));
} else{
	sprintf("reading from %s", opt$inputFile);
	inputTable = read.csv(opt$inputFile, row.names=1, header=T);
}
print("Count table read");
if(file.access(opt$de) ==-1){
	stop(sprintf("DE information file (%s) does not exist", opt$de));
} else{
	sprintf("reading from %s", opt$de);
	deFile = read.csv(opt$de, row.names=1, header=T);
}
print("DE information obtained");



if(suppressPackageStartupMessages(require("WGCNA"))){
	print("WGCNA is loaded correctly");
} else {
	print ("trying to install WGCNA")
	source("http://bioconductor.org/biocLite.R")
	biocLite();
	biocLite("WGCNA");
	if(suppressPackageStartupMessages(require("WGCNA"))){
		print ("WGCNA installed and loaded");
	} else{
		stop("could not install WGCNA");
	}
}
allowWGCNAThreads(opt$thread);

datExpr0 = as.data.frame(t(inputTable));
names(datExpr0) = rownames(inputTable);
rownames(datExpr0) = names(inputTable);

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK){
# Optionally, print the gene and sample names that were removed:

if (sum(!gsg$goodGenes)>0)
	printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
	printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


sampleTrees = hclust(dist(datExpr), method = "average")

#Plot the cluster
pdf(file = paste(opt$out,"sampleClustering.pdf", sep="_"), width = 12, height = 12);
	par(mfrow=c(2,1))
	par(mar = c(0, 4, 2, 0))
	plot(sampleTrees, main = "Sample clustering on all genes",xlab="", sub="", cex = 0.7);
dev.off();

powers = c(c(1:10), seq(from = 12, to=30, by=1))
#If it crash, it is most likely in this part. From my experience, it is usually because it has used too much cpu
sft=NULL;
if(opt$softpower <=0){
	sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
}
prev = NULL;
prevPower=NULL
#should allow users to select their power. three options, user number, estimated power, and our slope calculation model
selectedPower = NULL;
if(opt$softpower == -1){
	selectedPower=sft$powerEstimate
} else if(opt$softpower == 0){
	for(i in 2:length(sft$fitIndices$SFT.R.sq)){
		if(sft$fitIndices$slope[i] <0){
			if(is.null(prev)){
				prev = (sft$fitIndices$SFT.R.sq[i]-sft$fitIndices$SFT.R.sq[i-1])/(sft$fitIndices$Power[i]-sft$fitIndices$Power[i-1])
				prevPower =sft$fitIndices$Power[i-1];
			} else{
				newSlop = (sft$fitIndices$SFT.R.sq[i]-sft$fitIndices$SFT.R.sq[i-1])/(sft$fitIndices$Power[i]-sft$fitIndices$Power[i-1]);
				if(newSlop <=0){
					selectedPower=prevPower
					break;
				}
				prev = newSlop;
				prevPower = sft$fitIndices$Power[i]
			}
		}
	}
} else{
	selectedPower = opt$softpower;
}

sprintf("Selected soft power = %i", selectedPower);
mod= blockwiseModules(datExpr, corType="bicor", power=selectedPower, nThreads=opt$thread, maxBlockSize=nGenes, quickCor=0, impute=TRUE, minModuleSize=30,saveTOMs =TRUE)
ADJ1 = abs(bicor(datExpr))^selectedPower
TOM = TOMsimilarity(ADJ1);
Alldegrees1=intramodularConnectivity(ADJ1, mod$colors)
MEs=mod$MEs
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

includedGene = merge(inputTable, deFile, by="row.names")
geneInfo = data.frame(GeneID=names(datExpr), moduleColor=mod$colors, geneModuleMembership ,Alldegrees1, Bonferroni=includedGene$bon, fc=includedGene$log2FoldChange, pvalue=includedGene$pvalue)
colours = unique(geneInfo$moduleColor)
hyperEnrichment=matrix(NA,length(colours), 3)
row.names(hyperEnrichment)=colours;
colnames(hyperEnrichment)=c("PValue", "Size", "Padj");
for(i in 1:length(colours)){
	current = subset(geneInfo, geneInfo$moduleColor == colours[i])
	msig = nrow(subset(current, current$Bonferroni <=0.05))
	mtotal = nrow(current)
	tsig = nrow(subset(geneInfo, geneInfo$Bonferroni <=0.05))
	hyperEnrichment[i,1]=exp(phyper(msig, tsig, nrow(geneInfo)-tsig, mtotal, lower.tail = FALSE, log.p = TRUE))
	hyperEnrichment[i,2]=mtotal
}
hyperEnrichment = as.data.frame(hyperEnrichment)
hyperEnrichment = hyperEnrichment[order(hyperEnrichment$PValue),]

result = matrix(NA, length(unique(geneInfo$moduleColor)), 2)
row.names(result) = unique(geneInfo$moduleColor)
for(i in 1:nrow(result)){
	current = subset(geneInfo, geneInfo$moduleColor==row.names(result)[i])
	#result[i,1]=pnorm(sum(qnorm(current$pvalue/2, lower.tail=F))/sqrt(nrow(current)), lower.tail=F)*2
	result[i,1]=pnorm(sum(qnorm(current$pvalue/2, lower.tail=F)*sign(current$fc))/sqrt(nrow(current)), lower.tail=F)*2
	if(result[i,1]>1) result[i,1]=1;
	result[i,2] = nrow(current);
}
result = as.data.frame(result)
colnames(result) =c("pvalue", "size")
result$padj = p.adjust(result$pvalue, "BH")
result = result[order(result$pvalue),]
original.sample=data.frame(pvalue=includedGene$pvalue, fc=includedGene$log2FoldChange)
original.sample$value = sign(original.sample$fc)*qnorm(original.sample$pvalue/2, lower.tail=F);
#original.sample$value = qnorm(original.sample$pvalue/2, lower.tail=F);
boostrapSim= matrix(NA, length(unique(result$size)), 3);
boostrapSim[,1] = unique(result$size);
for(i in 1:nrow(boostrapSim)){
	print(paste("Running......Cycle.", i, " out of ", nrow(boostrapSim), sep=""))
	#initialize data.frame
	resample.number=opt$perm
	alpha=0.05
	sample.size=boostrapSim[i,1]
	resample.results<-data.frame("Run.Number"=NULL,"mean"=NULL)
	for(counter in 1:resample.number){
		#temp<-sample(original.sample, size=length(original.sample), replace = TRUE)
		temp<-sample(original.sample$value, size=sample.size, replace = TRUE)
		temp.mean<-sum(temp)/sqrt(sample.size)
		temp.table.row<-data.frame("Run.Number"=counter,"mean"=temp.mean)
		resample.results<-rbind(resample.results,temp.table.row)
	}
	boostrapSim[i,2]= mean(resample.results$mean)
	boostrapSim[i,3]= var(resample.results$mean)
}
colnames(boostrapSim) = c("Size", "Mean", "Var")
boostrapSim= as.data.frame(boostrapSim)
module.sig = matrix(NA,length( unique(geneInfo$moduleColor)),3)
row.names(module.sig) = unique(geneInfo$moduleColor);
colnames(module.sig) = c("size", "score", "final")
for(i in 1:length(unique(geneInfo$moduleColor))){
	current = subset(geneInfo, geneInfo$moduleColor==unique(geneInfo$moduleColor)[i])
	module.sig[i,1] = nrow(current)
	module.sig[i,2] = sum(sign(current$fc)*qnorm(current$pvalue/2, lower.tail=F))/sqrt(nrow(current))
	background = subset(boostrapSim, boostrapSim$Size == nrow(current));
	module.sig[i,3] = (module.sig[i,2]-background$Mean)/sqrt(background$Var)
}
module.sig = as.data.frame(module.sig);
module.sig$pvalue = rep(NA, nrow(module.sig))
for(i in 1:nrow(module.sig)){
	module.sig$pvalue[i] = min(pnorm(module.sig$final[i], lower.tail=T)*2, pnorm(module.sig$final[i], lower.tail=F)*2)
	if(module.sig$pvalue[i] > 1) module.sig$pvalue[i] = 1
}
module.sig = module.sig[order(module.sig$pvalue),]
module.sig$padj = p.adjust(module.sig$pvalue, "BH")
module.sig$bon = p.adjust(module.sig$pvalue, "bonferroni")

sig.modules =merge(module.sig, hyperEnrichment, by="row.names")
subset(sig.modules, sig.modules$bon<=0.05 & sig.modules$Padj<=0.05)
select.modules = subset(sig.modules, sig.modules$bon<=0.05 & sig.modules$Padj<=0.05)
row.names(select.modules) = select.modules$Row.names;
select.modules= select.modules[,-1]
row.names(sig.modules) = sig.modules$Row.names;
sig.modules= sig.modules[,-1]
print("Module significance calculated")
write.csv(sig.modules, paste(opt$out, ".mod.sig.csv", sep=""));
if(nrow(select.modules) == 0){
	stop("There are no significant modules for down-stream cell type analysis");
}
geneInfoBrief = data.frame(Gene.ID=geneInfo$GeneID, moduleColor = geneInfo$moduleColor, bon =geneInfo$Bonferroni, fc=geneInfo$fc)
if(!is.null(opt$mgi) & file.access(opt$mgi) ==-1){
	stop(sprintf("Count table (%s) does not exist", opt$mgi));
} else if(!is.null(opt$mgi)){
	sprintf("reading from %s", opt$mgi);
	mgi = read.csv(opt$mgi, row.names=1, header=T);
	colnames(mgi) = c("GeneID", "MGI")
	geneInfoType = merge(geneInfoBrief, mgi, by="GeneID");
	geneInfoDE = subset(geneInfoType, geneInfoType$bon <= 0.05)
}

 
selectList = row.names(select.modules);
modInfo=subset(geneInfoDE, (moduleColor %in% selectList))
cellTypeEnrichment=userListEnrichment(toupper(modInfo$MGI), modInfo$moduleColor, useBrainLists = TRUE,omitCategories = "grey", outputCorrectedPvalues = TRUE, minGenesInCategory = 30)
cellTypeEnrichment$pValues$CorrectedPvalues = p.adjust(cellTypeEnrichment$pValues$Pvalues,"BH")
cellTypeSig = subset(cellTypeEnrichment$pValues,cellTypeEnrichment$pValues$CorrectedPvalues <=0.05)
write.csv(cellTypeSig, paste(opt$out, ".celltype.csv", sep=""))






