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


#Perform the argument checking before spending time in loading DESeq2
option_list <- list(
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
	make_option(c("-q", "--quietly"), action="store_false",	dest="verbose", help="Print little output"),
	make_option(c("-f", "--file"), type="character",  help="The input count table", dest="inputFile"),
	make_option(c("-i", "--info"), type="character",  help="Sample information file. Should be tab delimited and header should start with #. Should be of the format <Sample Name> <Condition>, and the sample name should be unique", dest="samples", metavar= "sample info"),
	make_option(c("-p", "--plot"), action="store_true", default=FALSE, help="Generate the MA plot and Sample Clustering Graphs"),
	make_option(c("-c", "--clean"), action="store_true", default=FALSE, help="Remove no_feature, ambiguous, too_low_aQual, not_aligned, alignment_not_unique from the count matrix"),
	make_option(c("-w", "--wgcna"), action="store_true", default=FALSE, help="Generate mean centered table for WGCNA"),
	make_option(c("-s", "--spia"), action="store_true", default=FALSE, help="Perform SPIA analysis"),
	make_option(c("-t", "--trans"), type="character", help="Table to transform count table identifier to Entrez id", metavar="file"),
	make_option(c("-e", "--species"), type="character",default="mmu",  help="Animal species used for SPIA, allow hsa or mmu [default \"%default\"]"),
	make_option(c("-d", "--dir"), type="character", help="Directory of the KEGG pathways", metavar="directory"),
	make_option(c("-n", "--perm"), type="integer", default=1000, help="Number of permutation to perform [default \"%default\"]", metavar="number"),
	make_option(c("-o", "--out"), type="character", help="Output prefix [default \"count table name\"]", metavar="output")
)
 
opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$inputFile)){
	stop("Require count table input");
}
if(is.null(opt$samples)){
	stop("Require sample information");
}
if(is.null(opt$out)){
	opt$out = opt$inputFile;
}
if(opt$spia){
	if(is.null(opt$dir) & (opt$species != "mmu" | opt$species != "hsa")){
		stop("Only mmu and hsa are allowed for the default SPIA KEGG dataset")
	}
}
if(opt$perm < 1){
	stop("Number of permutation must be >= 1");
}


#Check if the DESeq2 package is installed
if(suppressPackageStartupMessages(require("DESeq2"))){
	print("DESeq2 is loaded correctly");
} else {
	print ("trying to install DESeq2");
	source("http://bioconductor.org/biocLite.R")
	biocLite();
	biocLite("DESeq2");
	if(suppressPackageStartupMessages(require("DESeq2"))){
		print("DESeq2 installed and loaded");
	} else{
		stop("could not install DESeq2");
	}
}

#Check if the SPIA package is installed
if(opt$spia){
	if(suppressPackageStartupMessages(require("SPIA"))){
		print("SPIA is loaded correctly");
	} else{
		print("trying to install SPIA");
		source("http://bioconductor.org/biocLite.R")
		biocLite();
		biocLite("SPIA");
		if(suppressPackageStartupMessages(require("SPIA"))){
			print("SPIA installed and loaded");
		} else{
			stop("could not install SPIA");
		}
	}
}

if(file.access(opt$inputFile) ==-1){
	stop(sprintf("Count table (%s) does not exist", opt$inputFile));
} else{
	inputTable = read.csv(opt$inputFile, row.names=1, header=T);
}
print("count table obtained")
if(file.access(opt$samples) ==-1){
	stop(sprintf("Sample information file (%s) does not exist", opt$samples));
} else{
	sampleInfo = read.table(opt$samples, row.names=1);
	#We will have to make sure the sample names are correct
	row.names(sampleInfo) = make.names(row.names(sampleInfo))
	inputTable= inputTable[,row.names(sampleInfo)]
}
print("sample information obtained, count table processed")

#print(head(inputTable));
#print(sampleInfo);

if(opt$clean){
	inputTable = inputTable[row.names(inputTable) %notin% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"),]
}

colData = data.frame(row.names=row.names(sampleInfo), condition=sampleInfo[,1])
dds <- DESeqDataSetFromMatrix(countData = inputTable, colData = colData, design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,levels=unique(sampleInfo[,1]))
#We will relevel to the first condition saw in the sample file. As a result of that, assuming control is the first condition observed, a positive logFC = increased expression in Cases w.r.t. Control
colData(dds)$condition <- relevel(colData(dds)$condition,  as.character(sampleInfo[1,1]))
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
#To handle outliers, we use replaceOutliersWithTrimmedMean
dds <- replaceOutliersWithTrimmedMean(dds) 
dds <- DESeq(dds)
res=results(dds)
#We filter genes with read counts less than 10 and that has NA as p-value (outliers that DESeq2 prefer to provide NA)
use <- res$baseMean >= 10 & !is.na(res$pvalue)
resFilt=as.data.frame(res[use,])
#Can filter ERCC transcript here as they are only considered for size factor calculation. However, will need to add one additional parameter. Might do it later.
resFilt$padj = p.adjust(resFilt$pvalue, "BH")
resFilt$bon = p.adjust(resFilt$pvalue, "bonferroni")

#Generate the csv file for DEGs
write.csv(resFilt, paste(opt$out, ".deseq2.csv", sep=""));

#Graph plotting related codes
if(opt$plot){
	#Check the required packages
	if(suppressPackageStartupMessages(require(gplots))){
		print("gplots is loaded correctly");
	} else {
		print ("trying to install gplots")
		install.packages("gplots",repos="http://cran.rstudio.com/");
		if(suppressPackageStartupMessages(require("gplots"))){
			print ("gplots installed and loaded");
		} else{
			print("could not install gplots, will not plot the sample clustering plot");
			opt$plot = FALSE;
		}
	}
	if(suppressPackageStartupMessages(require(RColorBrewer))){
		print("RColorBrewer is loaded correctly");
	} else {
		print ("trying to install RColorBrewer")
		install.packages("RColorBrewer",repos="http://cran.rstudio.com/");
		if(suppressPackageStartupMessages(require("RColorBrewer"))){
			print ("RColorBrewer installed and loaded");
		} else if(opt$plot){
			print("could not install RColorBrewer, will not plot the sample clustering plot");
			opt$plot=FALSE;
		} else{
			print("could not install RColorBrewer, it is also required to plot the sample clustering plot");
		}
	}
	if(opt$plot){
		vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
		distsRL <- dist(t(assay(vsd)))
		mat <- as.matrix(distsRL)
		rownames(mat) <- colnames(mat) <- paste(row.names(colData(dds)),with(colData(dds), condition), sep=":")
		hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
		pdf(paste(opt$out, ".sampleClustering.pdf", sep=""))
			heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
		dev.off()
		pdf(paste(opt$out, ".pca.pdf", sep=""))
			print(plotPCA(vsd, intgroup=c("condition")))
		dev.off()
		pdf(paste(opt$out, ".MAplot.pdf", sep=""))
			DESeq2::plotMA(dds)
		dev.off()
	}	
}

#WGCNA Mean centring
if(opt$wgcna){
	if(! exists("vsd")){
		vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
	}
	count = assay(vsd)
	#only include genes that were passed the filter in DE analysis
	count=count[row.names(resFilt),]
	#We want to mean normalize the data
	count.mean = count;
	for(i in 1:nrow(count.mean)){
		count.mean[i,1:5] = count[i,1:5]-mean(count[i,1:5]);
		count.mean[i,6:10] = count[i,6:10]-mean(count[i,6:10]);
	}
	#Now the count.mean is ready for WGCNA analysis
	write.csv(count.mean, paste(opt$out,"wgcna.csv", sep=""));
}

if(opt$SPIA){
	if(is.na(opt$trans)){
		#Not providing the translation table, will consider the row.name from inputTable to  be entrez id
	}
}