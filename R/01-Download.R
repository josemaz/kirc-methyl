require(TCGAbiolinks)
# require(SummarizedExperiment)

dir.create("data/RDS",recursive = TRUE)
dir.create("data/TCGA",recursive = TRUE)

#! DOWNLOADING
setwd("data/TCGA")
qry.meth <-  GDCquery(project = "TCGA-KIRC",
	data.category = "DNA Methylation",
	platform="Illumina Human Methylation 450")
GDCdownload(qry.meth)
meths.raw <- GDCprepare(qry.meth, summarizedExperiment=TRUE)
setwd("../..")
saveRDS(meths.raw, file = "data/RDS/meths-raw.rds")
print(dim(meths.raw)) # 485577 x 485
# meths.raw <- readRDS(file = "data/RDS/meths-raw.rds")

\
qry.rna <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")

common.patients <- intersect(substr(getResults(qry.meth, cols = "cases"), 1, 12),
                                 substr(getResults(qry.rna, cols = "cases"), 1, 12))
print(length(common.patients)) #317

setwd("data/TCGA")
qry.rna <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts",
                      barcode = common.patients)
GDCdownload(qry.rna)
rnas.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE)
setwd("../..")
saveRDS(rnas.raw, file = "data/RDS/rnas-raw.rds")
print(dim(rnas.raw)) #56602 x 394




