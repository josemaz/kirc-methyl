require(TCGAbiolinks)
require(SummarizedExperiment)
require(tictoc)

meths <- readRDS(file = "data/RDS/meths-clean.rds")

dpng <- "data/plots/meth-stats/"
dir.create(dpng, recursive = TRUE)

# met.sample <- met[sample(nrow(met),100000),]
s.sample <- c(1e2,1e3,1e4,1e5, nrow(meths) )
for (i in s.sample){
    tic(paste0("Heatmap Eleapsed Time: ",i))
    met.sample <- meths[sample(nrow(meths),i),]    
    TCGAvisualize_Heatmap(assay(met.sample),
                col.metadata = colData(met.sample)[,c("barcode","tipo")],
                # row.metadata = row.mdat,
                col.colors = list( grupos = c("NT" = "blue", 
                                                "Stage1"="red",
                                                "Stage2"="green",
                                                "Stage3"="pink",
                                                "Stage4"="grey")
                                ),
                type = "methylation",
                cluster_columns =  TRUE,
                cluster_rows = FALSE,
                filename = paste0(dpng,"heatmap-",i,".pdf")
            )
    toc()
    
}

tic("Mean plot Eleapsed Time: ")
TCGAvisualize_meanMethylation( meths,
            groupCol  = "tipo", sort = "mean.desc", 
            filename=paste0(dpng,"meandesc.png")
        )
toc()