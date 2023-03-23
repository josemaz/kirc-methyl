suppressWarnings(require(foreach))
suppressWarnings(require(doParallel))
suppressWarnings(require(crayon))
suppressWarnings(require(SummarizedExperiment))

#! LOADING
cat(green("LOADING ... \n"))
meths <- readRDS(file = "data/RDS/meths-clean.rds")
cpgs.info <- rowData(meths)
dems <- readRDS(file = "data/RDS/dms.rds")
gdegs <- readRDS(file = "data/RDS/gdegs.rds")

#! create and register cluster
n.cores <- 16
doParallel::registerDoParallel(n.cores)

#! create lists with candidate genes
nominees <- list()
for (nom.dem in names(dems)){ 
    cat(green(nom.dem,"\n"))
    dem <- dems[[nom.dem]]
    hypo <- dem[grep("Hypo", dem$status),]
    hyper <- dem[grep("Hyper", dem$status),]
    m.hypo <- meths[rownames(meths) %in% rownames(hypo) , ]
    m.hyper <- meths[rownames(meths) %in% rownames(hyper) , ]
    gr1name = unlist(strsplit(nom.dem,"_"))[1]
    gr2name = unlist(strsplit(nom.dem,"_"))[2]
    g.up <- gdegs$ups[[nom.dem]]$gene_name
    cat(red("\tOver\n"))
    dfup <- foreach(
        gene = g.up
    ) %dopar% {
        lcpgs <- rownames(cpgs.info[grep(gene, cpgs.info$Gene_Symbol),])
        m <- m.hypo[rownames(m.hypo) %in% lcpgs,]
        if(nrow(m) > 0){
            return(c(gene,nrow(m)*100/length(lcpgs),nrow(m),length(lcpgs)))
        }else{
            return(c())
        }

    }
    cat(red("\tUnder\n"))
    g.down <- gdegs$down[[nom.dem]]$gene_name
    dfdown <- foreach(
        gene = g.down
    ) %dopar% {
        lcpgs <- rownames(cpgs.info[grep(gene, cpgs.info$Gene_Symbol),])
        m <- m.hyper[rownames(m.hyper) %in% lcpgs,]
        if(nrow(m) > 0){
            return(c(gene,nrow(m)*100/length(lcpgs),nrow(m),length(lcpgs)))
        }else{
            return(c())
        }
    }
    dfup[sapply(dfup, is.null)] <- NULL    
    dfup <- as.data.frame(t(as.data.frame(dfup)))
    colnames(dfup) <- c("gname","pct cpgs-meth","ncpgs-meths","total_cpgs")
    rownames(dfup)<-NULL
    dfdown[sapply(dfdown, is.null)] <- NULL
    dfdown <- as.data.frame(t(as.data.frame(dfdown)))
    colnames(dfdown) <- c("gname","pct cpgs-meth","ncpgs-meths","total_cpgs")
    rownames(dfdown)<-NULL
    nominees[[nom.dem]] <-  list(over.hypo = dfup, under.hyper = dfdown, 
                                    mhypo = m.hypo, mhyper = m.hyper)
}

saveRDS(nominees, file = "data/RDS/nominees.rds")

dirtabs = "data/tables/nominees"
dir.create(dirtabs, recursive = TRUE)
for (contra in names(nominees)){
    print(contra)
    fout <- paste0(dirtabs,"/",contra,"-overhypo.tsv")
    write.table(nominees[[contra]]$over.hypo, 
            file = fout, row.names=FALSE, sep="\t", quote=FALSE)
    fout <- paste0(dirtabs,"/",contra,"-underhyper.tsv")
    write.table(nominees[[contra]]$under.hyper, 
            file = fout, row.names=FALSE, sep="\t", quote=FALSE)
}