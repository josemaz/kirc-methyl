suppressWarnings(require(TCGAbiolinks))
suppressWarnings(require(SummarizedExperiment))
suppressWarnings(require(crayon))
suppressWarnings(require(tictoc))
source("R/libtools.R")


pngd <- "data/plots/meth-volcanos/"
dir.create(pngd, recursive=TRUE)

#! LOAD RNAs and METHs NORMALIZED
rnas <- readRDS(file = "data/RDS/rnas-norm.rds")
# meths <- readRDS(file = "meths-norm.rds")
meths <- readRDS(file = "data/RDS/meths-clean.rds")
cpgs.info <- rowData(meths)

# #! LOAD DEGs
gdegs <- readRDS(file = "data/RDS/gdegs.rds")


# etapas <- c("NT","Stage1", "Stage2", "Stage3","Stage4")
# contras <- data.frame(gr1=character(), gr2=character())
# # for (i in 1:(length(etapas)-1)) contras[i,] <- c(etapas[i],etapas[i+1])
# for (i in 1:(length(etapas)-1)) contras[i,] <- c(etapas[1],etapas[i+1])

grnas <- rnas
rownames(grnas) <- rowData(rnas)$gene_name

#! GENERA LOS RESULTADOS DE METILACION DIFERENCIAL
dems <- list()
for (contra in names(gdegs$up)){
    tic(paste("Time:",contra))
    gr1name = unlist(strsplit(contra,"_"))[1]
    gr2name = unlist(strsplit(contra,"_"))[2]
    plt <- TCGAanalyze_DMC(meths, groupCol = "tipo",
                group1 = gr2name,
                group2 = gr1name,
                diffmean.cut = 0.15,
                p.cut = 0.05,
                save = FALSE,
                plot.filename = paste0(pngd, gr1name, "vs", gr2name,".png"),
                cores = parallel::detectCores())
    print(head(plt))
    dems[[contra]] <- plt
    toc()
}
saveRDS(dems, file = "data/RDS/dms.rds")






# dems <- readRDS(file = "dems.rds")
# #! create and register cluster
# n.cores <- 16
# doParallel::registerDoParallel(n.cores)
# source("../03-1-famous.R")
# source("../03-2-methdriven.R")





















        # # gaussianas en y
        # mixmdl <- normalmixEM(vals$gcpgs, mu = c(mean(gcpgs1), mean(gcpgs2)), k = 3)
        # # mixdensity(mixmdl,paste0(dout,gene,"-k2.png"))
        # mixmdl <- normalmixEM(gcpgs2, k = 3)

        # ggplot() + geom_density(data = vals, aes(x = gcpgs, color = tipo)) +
        # stat_function(geom = "line", fun = plot_mix_comps,
        #   args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
        #   colour = "red", lwd = 1.5) +
        # stat_function(geom = "line", fun = plot_mix_comps,
        #   args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
        #   colour = "blue", lwd = 1.5)

        # ggplot() + geom_density(data = values2, aes(x = gcpgs)) +
        # stat_function(geom = "line", fun = plot_mix_comps,
        #   args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
        #   colour = "red", lwd = 1.5)


# for (nom.dem in names(dems)[1:4]){ # hack: para control vs etapas
#     dout <- paste0(pngd,nom.dem,"/")
#     suppressWarnings(dir.create(dout))
#     cat(green(nom.dem,"\n"))
#     dem <- dems[[nom.dem]]
#     hypo <- dem[grep("Hypo", dem$status),]
#     hyper <- dem[grep("Hyper", dem$status),]
#     m.hypo <- meths[rownames(meths) %in% rownames(hypo) , ]
#     m.hyper <- meths[rownames(meths) %in% rownames(hyper) , ]
#     gr1name = unlist(strsplit(nom.dem,"_"))[1]
#     gr2name = unlist(strsplit(nom.dem,"_"))[2]
#     g.up <- gdegs$ups[[nom.dem]]$gene_name
#     for (gene in g.up){
#         lcpgs <- rownames(cpgs.info[grep(gene, cpgs.info$Gene_Symbol),])
#         m <- m.hypo[rownames(m.hypo) %in% lcpgs,]     
#         if(nrow(m) > 0){
#             print(paste(gene,":", sprintf("%.2f",nrow(m)*100/length(lcpgs)),
#                 "%cpgs hypomethylated", nrow(m), length(lcpgs)))
            
#             gexpr1 <- as.vector(assay(grnas[rownames(grnas) == gene, grnas$tipo == gr1name]))
#             m <- m.hypo[ rownames(m.hypo) %in% lcpgs, m.hypo$tipo == gr1name]
#             gcpgs1 <- colMedians(assay(m))            
#             vals1 <- data.frame(gexps = gexpr1, gcpgs = gcpgs1, tipo = gr1name)
            
#             gexpr2 <- as.vector(assay(grnas[rownames(grnas) == gene, grnas$tipo == gr2name]))
#             m <- m.hypo[ rownames(m.hypo) %in% lcpgs, m.hypo$tipo == gr2name]
#             gcpgs2 <- colMedians(assay(m))
#             vals2 <- data.frame(gexps = gexpr2, gcpgs = gcpgs2, tipo = gr2name)
            
#             vals <- rbind(vals1,vals2)
#             suppressWarnings( scatplot(vals, paste0(dout,gene,".png"), 
#                 median(gcpgs1), median(gcpgs2)) )
#             next
#         }
#         cat(blue(gene," sin CPGS\n"))
#     }
# }









# dems <- list()
# for (i in 1:nrow(contras)){
#     contra <- paste0(contras[i,"gr1"],"_",contras[i,"gr2"])
#     dems[[contra]] <- TCGAanalyze_DMC(meths, groupCol = "tipo",
#                             group1 = contras[i,"gr2"],
#                             group2 = contras[i,"gr1"],
#                             diffmean.cut = 0.15,
#                             p.cut = 0.05,
#                             save = FALSE,
#                             plot.filename = NULL,
#                             cores = parallel::detectCores())
# }
# saveRDS(dems, file = "dems.rds")
# dems <- readRDS(file = "dems.rds")






# dems <- sapply(names(gdegs$up),function(contra){
#     gr1name = unlist(strsplit(contra,"_"))[1]
#     gr2name = unlist(strsplit(contra,"_"))[2]
#     plt <- TCGAanalyze_DMC(meths, groupCol = "tipo",
#                 group1 = gr2name,
#                 group2 = gr1name,
#                 diffmean.cut = 0.15,
#                 p.cut = 0.05,
#                 save = FALSE,
#                 plot.filename = NULL,
#                 cores = parallel::detectCores())
#     print(head(plt))
#     return(plt)
# })


# https://www.roelpeters.be/how-to-add-a-regression-equation-and-r-squared-in-ggplot2/
# https://www.learnbymarketing.com/tutorials/linear-regression-in-r/