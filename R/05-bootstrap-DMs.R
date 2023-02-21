suppressWarnings(require(foreach))
suppressWarnings(require(doParallel))
suppressWarnings(require(crayon))
suppressWarnings(require(SummarizedExperiment))
suppressWarnings(require(TCGAbiolinks))
suppressWarnings(require(tictoc))

#! create and register cluster
n.cores <- 16
doParallel::registerDoParallel(n.cores)

meths <- readRDS(file = "data/RDS/meths-clean.rds")
# met.sample <- met[sample(nrow(meths),10000),]
max.samples <- 31 # stage 2 has 31 samples
max.samples <- 24 # Control has 24 samples
dm.cut <- 0.15
pval <- 0.05
niter <- 50
print(dim(meths))

contras <- combn(levels(meths$tipo), 2)
contras <- contras[,c(1:4,5,8,10)]
# contras <- contras[,c(8,10)]
df <- data.frame(tipo=character(), valor=numeric(), contrast=character())
tic("Eleapsed Time: ")
for (i in 1:ncol(contras)){
    g1 <- contras[1,i]
    g2 <- contras[2,i]
    vs <- paste0(g1,"_",g2)
    cat(green("Contrast: ",vs, "\n"))
    lhyper <- c()
    lhypo <- c()
    for (itera in 1:niter){ # 50 iterations for paper
        cat(green("Iteration: ", itera, "\n"))
        # print(paste(g1,g2))
        muestras1 <- meths[,meths$tipo==g1]$barcode
        muestras1 <- sample(muestras1,max.samples)
        muestras2 <- meths[,meths$tipo==g2]$barcode
        muestras2 <- sample(muestras2,max.samples)
        mets <- meths[,meths$barcode %in% c(muestras1,muestras2)]
        dms <- TCGAanalyze_DMC(mets, groupCol = "tipo",
                                group1 = g1,
                                group2 = g2,
                                diffmean.cut = dm.cut,
                                p.cut = pval,
                                save = FALSE,
                                # plot.filename = paste0(g1,"vs",g2,"_metvolcano.png"),
                                plot.filename = NULL,
                                cores = parallel::detectCores(),
                                )
        df[nrow(df) + 1,] = c("hyper", nrow(dms[dms[3]>dm.cut,]), vs)
        df[nrow(df) + 1,] = c("hypo", nrow(dms[dms[3]<(-dm.cut),]), vs)
    }
}    
# print(df)
df$valor <- as.integer(df$valor)
df$contrast <- as.factor(df$contrast)
saveRDS(df, file = "data/RDS/bootstrap-dms.rds")
toc()








# colnames(dat$up) <- c("NTvsSt1","St1vsSt2","St2vsSt3","St3vsSt4")
# colnames(dat$down) <- c("NTvsSt1","St1vsSt2","St2vsSt3","St3vsSt4")
# print(dat$up)
# print(dat$down)