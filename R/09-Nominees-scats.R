require(foreach)
require(doParallel)
require(crayon)
source("R/libtools.R")

#! create and register cluster
n.cores <- 16
doParallel::registerDoParallel(n.cores)

#! LOADING
rnas <- readRDS(file = "data/RDS/rnas-norm.rds")
grnas <- rnas
rownames(grnas) <- rowData(rnas)$gene_name
meths <- readRDS(file = "data/RDS/meths-clean.rds")
cpgs.info <- rowData(meths)
nominees <- readRDS(file = "data/RDS/nominees.rds")

#! Plotting
stats <- list()
for (contra in names(nominees)){
    print(contra)    

    gr1name <- unlist(strsplit(contra,"_"))[1]
    gr2name <- unlist(strsplit(contra,"_"))[2]

    #! OVEREXPRESSED - HYPOMETHYLATED IN CANCER
    #! plots/scats/NT_Stage1/over
    dout <- paste0('data/plots/scats/',contra,'/over')
    dir.create(dout, recursive =TRUE)
    genes <- nominees[[contra]]$over.hypo
    cat(red("\tOver\n"))
    # for (g in 1:nrow(genes)){
    ltmp <- foreach(
        g = 1:nrow(genes)
    ) %dopar% {
        gene <- genes[g,"gname"]
        # print(gene)
        cpgs <- rownames(cpgs.info[grep(gene, cpgs.info$Gene_Symbol),])
        
        gexpr1 <- as.vector(assay(grnas[rownames(grnas) == gene, grnas$tipo == gr1name]))
        m <- meths[ rownames(meths) %in% cpgs, meths$tipo == gr1name]        
        gcpgs1 <- colMedians(assay(m))
        values1 <- data.frame(gexps = gexpr1, gcpgs = gcpgs1, tipo = gr1name)

        gexpr2 <- as.vector(assay(grnas[rownames(grnas) == gene, grnas$tipo == gr2name]))
        m <- meths[ rownames(meths) %in% cpgs, meths$tipo == gr2name]
        gcpgs2 <- colMedians(assay(m))
        values2 <- data.frame(gexps = gexpr2, gcpgs = gcpgs2, tipo = gr2name)

        vals <- rbind(values1,values2)
        
        fo <- paste0(dout,"/",gene,".png")
        scatplot(vals, fo, median(gcpgs1), median(gcpgs2))

        # dist.filter(values1) # grupo 1
        
        return(c(gene, stats.gaus(vals,gr1name,gr2name)))

    }
    df1 <- as.data.frame(t(as.data.frame(ltmp)))
    rownames(df1) <- NULL
    colnames(df1) <- c("gene", "m",
            "gr1.cpgs.mea","gr1.cpgs.sig", "gr1.cpgs.med",
            "gr2.cpgs.mea","gr2.cpgs.sig", "gr2.cpgs.med",
            "gr1.gexps.mea","gr1.gexps.sig"," gr1.gexps.med",
            "gr2.gexps.mea","gr2.gexps.sig", "gr2.gexps.med" )

    #! UNDEREXPRESSED - HYPERMETHYLATED IN CANCER
    dout <- paste0('data/plots/scats/',contra,'/under')
    dir.create(dout, recursive=TRUE)
    genes <- nominees[[contra]]$under.hyper
    cat(red("\tUnder\n"))
    # for (g in 1:nrow(genes)){
    ltmp <- foreach(
        g = 1:nrow(genes)
    ) %dopar% {
        gene <- genes[g,"gname"]
        cpgs <- rownames(cpgs.info[grep(gene, cpgs.info$Gene_Symbol),])

        gexpr1 <- as.vector(assay(grnas[rownames(grnas) == gene, grnas$tipo == gr1name]))
        m <- meths[ rownames(meths) %in% cpgs, meths$tipo == gr1name]
        gcpgs1 <- colMedians(assay(m))
        values1 <- data.frame(gexps = gexpr1, gcpgs = gcpgs1, tipo = gr1name)

        gexpr2 <- as.vector(assay(grnas[rownames(grnas) == gene, grnas$tipo == gr2name]))
        m <- meths[ rownames(meths) %in% cpgs, meths$tipo == gr2name]
        gcpgs2 <- colMedians(assay(m))
        values2 <- data.frame(gexps = gexpr2, gcpgs = gcpgs2, tipo = gr2name)

        vals <- rbind(values1,values2)

        fo <- paste0(dout,"/",gene,".png")
        scatplot(vals,fo, median(gcpgs1), median(gcpgs2))

        # dist.filter(values1) # grupo 1
        
        return(c(gene, stats.gaus(vals,gr1name,gr2name)))

    }
    df2 <- as.data.frame(t(as.data.frame(ltmp)))
    rownames(df2) <- NULL
    colnames(df2) <- c("gene", "m",
            "gr1.cpgs.mea","gr1.cpgs.sig", "gr1.cpgs.med",
            "gr2.cpgs.mea","gr2.cpgs.sig", "gr2.cpgs.med",
            "gr1.gexps.mea","gr1.gexps.sig"," gr1.gexps.med",
            "gr2.gexps.mea","gr2.gexps.sig", "gr2.gexps.med" )

    stats[[contra]] <-  list(over.hypo = df1, under.hyper = df2)
    # break
}

saveRDS(stats, file = "data/RDS/stats.rds")