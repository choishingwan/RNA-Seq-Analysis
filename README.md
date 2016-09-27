# NOTE
The code here uses the old RNA Seq pipeline. 
Nowaday, it is better to use tools such as Salmon and tximport as suggested by the big names in RNA Seq. 
For more details, you can read [here](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lect-07-modern-rnaseq/ModernRNAseqAnalysis.pdf)



RNA-Seq-Analysis
================

The programming codes use for RNA Sequencing analysis including the use of STAR, htseq-count, DESeq2 and WGCNA.

perSampleProcess.pl will automatically generate bash scripts for the analysis where user can distribute the jobs manually for maximum efficiency

generateMatrix.pl is used after perSampleProcess.pl so that a count matrix can be generated for down-stream analysis

deAnalysis.r is the Rscript used for the DE analysis. It will try to resolve the dependency wherever possible. 

wgcna.r is the Rscript used for the WGCNA analysis. It will try to resolve the dependency wherever possible

Recommended Reference
================
Reference Fasta   ftp://ftp.ensembl.org/pub/release-74/fasta/mus_musculus/dna/Mus_musculus.GRCm38.74.dna.primary_assembly.fa.gz

Reference GTF     ftp://ftp.ensembl.org/pub/release-74/gtf/mus_musculus/Mus_musculus.GRCm38.74.gtf.gz


Depending Softwares
================
STAR          https://code.google.com/p/rna-star/

htseq-count   http://www-huber.embl.de/users/anders/HTSeq/

R             http://cran.r-project.org/
