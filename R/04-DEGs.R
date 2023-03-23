require(TCGAbiolinks)
require(SummarizedExperiment)
# require(EnhancedVolcano)
require(crayon)

lfc.cut = 2.0
pval = 0.05 #1e-5

pngd <- "data/plots/rna-volcanos/"
suppressWarnings(dir.create(pngd))

#! LOAD RNA DATA NORMALIZED
rnas.norm <- readRDS(file = "data/RDS/rnas-norm.rds")

#! FIRST ANALYSIS
cat(green("DEA Normal Tumor vs Primary Tumor ...\n"))
rnas.norm.nT <- rnas.norm[,rnas.norm$definition == "Solid Tissue Normal"]
rnas.norm.pT <- rnas.norm[,rnas.norm$definition == "Primary solid Tumor"]
dataDEGs <- TCGAanalyze_DEA( assay(rnas.norm.nT), assay(rnas.norm.pT),
                        Cond1type = "Normal",
                        Cond2type = "Tumor",
                        method = "glmLRT")
TCGAVisualize_volcano(
  dataDEGs$logFC,
  dataDEGs$PValue,
  filename = paste0(pngd,"NTvsPT_rna_volcano.png"),
  y.cut = pval,
  x.cut = lfc.cut,
  names = dataDEGs$gene_name,
  legend = "NT vs pT",
  names.fill = FALSE
)

#! SECOND ANALYSIS
cat(green("DEA pairwise Stages ...\n"))
etapas <- c("NT","Stage1", "Stage2", "Stage3","Stage4")
contras <- combn(etapas, 2)
contras <- contras[,1:4] # Hack to make only Ctrl-Stages
degs <- list()
for (i in 1:ncol(contras)){
    grname1 <- contras[1,i]
    grname2 <- contras[2,i]
    cat(grname1,grname2,"\n")
    grupo1 <- rnas.norm[,rnas.norm$tipo == grname1]
    grupo2 <- rnas.norm[,rnas.norm$tipo == grname2]
    dataDEGs <- TCGAanalyze_DEA( assay(grupo1), assay(grupo2),
                        Cond1type = grname1,
                        Cond2type = grname2,
                        method = "glmLRT")
    dataDEGs <- subset(dataDEGs, select = -c(gene_name)) # Fix it merge
    dataDEGs <- merge(dataDEGs,rowData(rnas.norm), by=0)
    TCGAVisualize_volcano(
            dataDEGs$logFC,
            dataDEGs$PValue, #1e-5 o 0.05 
            filename = paste0(pngd, grname1, "vs", grname2,"_rna_volcano.png"),
            y.cut = pval,
            x.cut = lfc.cut,
            names = dataDEGs$gene_name,
            legend = paste0(grname1, " vs ",grname2),
            names.fill = FALSE
    )    
    degs[[paste0(grname1,"_",grname2)]] <- dataDEGs
}


cat(green("Calculating Overexpressed genes ...\n"))
ups <- lapply(degs,function(deg){ 
          return(deg[(deg$logFC > lfc.cut) & (deg$FDR < pval),  ])
        })
lapply(ups,function(up){ print(dim(up)) })

cat(green("Calculating Subexpreseed genes ...\n"))
downs <- lapply(degs,function(deg){ 
          return(deg[ (deg$logFC < (-lfc.cut)) & (deg$FDR < pval),  ]) # pval=0.05
        })
lapply(downs,function(down){ print(dim(down)) })

gdegs <- list( ups = ups, downs = downs) # all gene in all contrasts regulated


#! SAVING
cat(green("Saving DEG data, stages, upregulated, downregulated ...\n"))
saveRDS(degs, file = "data/RDS/degs.rds")
saveRDS(gdegs, file = "data/RDS/gdegs.rds")


