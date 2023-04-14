r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

cat("Using LiBRARY: ",.libPaths(),"\n")

install.packages("crayon")
require("crayon")

# System Packages
packs <- c("UpSetR", "BiocManager", "DT", "tidyverse", "devtools", 
           "VennDiagram", "ggfortify", "gprofiler2", "cowplot",
           "tictoc", "openxlsx", "ggvenn", "igraph", "doParallel")
noinst <- setdiff(packs, rownames(installed.packages()))
cat(red(noinst) %+% '\n')
install.packages(noinst)
# Bioconductor Packages
packs <- c("TCGAbiolinks", "EDASeq", "edgeR", "EnhancedVolcano",
           "sva", "NOISeq","biomaRt", "DESeq2", "ComplexHeatmap")
noinst <- setdiff(packs, rownames(installed.packages()))
cat(red(noinst) %+% '\n')
BiocManager::install(noinst)
# From GITHUB
packs <- data.frame( p =  c("vqv/ggbiplot") )
packs$p2 <-  sapply(strsplit(packs$p, '/'), "[[", 2)
packs <- merge(packs,installed.packages(),by.x="p2",by.y=0,all.x=TRUE)
packs <- packs[is.na(packs$Package),]
cat(red(packs$p) %+% '\n')
devtools::install_github(packs$p)
