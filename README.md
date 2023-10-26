# RNA Sequencing Project for BS6202: Techniques in Biomedical Data Mining Final Project

### Synopsis

Code Repository for our analysis to Understand Biomarkers in Gastric Cancer in the context of overexpression of HOXD9 gene. We downloaded the RNA-seq data from Gene Expression Omnibus database (accession code GSE210016) and found that the sample may contain gDNA / rRNA contamination from a relatively lower uniquely mapping % using STAR aligner. gDNA / rRNA contamination are usually under-reported. We checked for rRNA using bbmap and found that the fastq files contain at least 20% rRNA contamination. Due to time constraint, we did not repeat the analysis after removal of rRNA contamination. (fasta for rRNA is included in the repo). Subsequently, we performed Differential Expression Analysis using DeSeq2 and compared our results to the reference paper. We found that LincRNA level of gene expression increases when using ballgown while protein coding genes level of gene expression increases when we employed STAR and DeSeq2. Ballgown was previously shown to be sensitive to low gene expression count. We also performed Differentially Expressed Alternative Splicing analysis using SplAdder (On a personal level, I wanted to learn how to analyse for Alternative Splicing but spend a lot more time on learning about / identifying gDNA / rRNA contaminants) and found 14 ASEs, mostly retained Introns, which could be False Positives. We also performed GO and KEGG enrichment analysis. We also created machine models to identify potential biomarkers to predict for Gastric Cancer Tumour. 

### Contributions

Irfan contributed to the conception of the project. Irfan contributed to the raw data quality assurance analysis, differences in DEGs results, and Differential Splicing Analysis subsections. Li Ruisi and Hu Yuhan contributed to the DEG analysis, and GO/KEGG analysis. Olivia Loh and Yang Yizhuo contributed to the section on Relative Importance of DEGs as Biomarkers to Predict for Gastric Cancer. A report is also attached.