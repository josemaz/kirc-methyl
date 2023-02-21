require(TCGAbiolinks)
require(SummarizedExperiment)
source("R/libtools.R")

rnas <- readRDS(file = "data/RDS/rnas-clean.rds")
# meths <- readRDS(file = "data/RDS/meths-clean.rds")


# #! Normalizations
ptm <- proc.time() # Aprox. 1 hr
#! RNAS
rnas.norm <- normal.proc(rnas, "tipo", do.pca = TRUE, expr.tipo="rna")
saveRDS(rnas.norm, file = "data/RDS/rnas-norm.rds")
#! METHS (PCAs por afuera por si truena el X11)
# mypca <- mypca(meths, "tipo",  fout = "PCA-meths-BeforeNorm.png")
# meths.norm <- norm(meths, "tipo", do.pca = FALSE, expr.tipo="meth")
# saveRDS(meths.norm, file = "meths-norm.rds")
# mypca <- mypca(meths, "tipo",  fout = "PCA-meths-AfterNorm.png")
cat("Eleapsed Time", proc.time() - ptm,"\n")

dout="data/tables/expression"
dir.create(dout, recursive = TRUE)
for (tipo in levels(rnas$tipo)){
    df <- assay(rnas[,rnas$tipo == tipo])
    fout <- paste0(dout,"/",tipo,".tsv")
    df <- cbind(gname = rownames(df), df)
    write.table(df, file = fout, sep="\t", quote=FALSE, row.names=FALSE)
}





