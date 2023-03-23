
require(SummarizedExperiment)

# https://yiweiniu.github.io/blog/2019/06/Cancer-gene-collections/
# X <- read.csv(url("https://raw.githubusercontent.com/biocq/DORGE/master/DORGE_shiny/data/Canyon_info.txt"),sep = '\t')
# table(X$Sample_name)
# https://www.biorxiv.org/content/10.1101/2020.07.21.213702v1.full.pdf

page <- "https://bioinfo.uth.edu/TSGene/Human_TSGs.txt?csrt=18093088261416048799"
TSG <- read.csv(url(page),sep = '\t')
TSG$ftsg <- 1
TSG <- TSG[ , c("GeneSymbol", "ftsg")]

page <- "http://ongene.bioinfo-minzhao.org/ongene_human.txt"
OG <- read.csv(url(page),sep = '\t')
OG$fog <- 1
OG <- OG[ , c("OncogeneName", "fog")]


rnas <- readRDS("data/RDS/rnas-norm.rds")
gdata <- data.frame(gname = rowData(rnas)$gene_name)
rownames(gdata) <- rowData(rnas)$ensembl_gene_id
dim(TSG)
dim(gdata)

df <- merge(gdata,TSG,by.x ="gname",by.y ="GeneSymbol",all.x = TRUE)
df <- merge(df,OG,by.x ="gname",by.y ="OncogeneName",all.x = TRUE)
dim(df)
df$typo <- "None"
df["typo"][df["ftsg"] == 1]<- "tsg"
df["typo"][df["fog"] == 1]<- "og"
qry <- (df["ftsg"] == 1) & (df["fog"] == 1)
df["typo"][qry]<- "both"
df <- df[,c("gname","typo")]
table(df$typo)

write.table(df,file = "data/tables/onco-tsg-db.tsv", sep = '\t', row.names = FALSE)
