library(openxlsx)


stats <- readRDS(file = "data/RDS/stats.rds")

for (vs in names(stats)){
    print(vs)
    dout <- paste0("data/tables/stat-meths/",vs,"/")
    suppressWarnings(dir.create(dout,recursive = TRUE))

    oh <- stats[[vs]]$over.hypo
    for (i in 2:ncol(oh)){ class(oh[, i]) = "numeric"}
    oh$diffexpr <- oh$gr2.gexps.mea - oh$gr1.gexps.mea    
    write.csv(x=oh, file=paste0(dout,"over-hypo.csv"), row.names = FALSE)

    uh <- stats[[vs]]$under.hyper
    for (i in 2:ncol(uh)){ class(uh[, i]) = "numeric"}
    uh$diffexpr <- uh$gr2.gexps.mea - uh$gr1.gexps.mea
    write.csv(x=uh, file=paste0(dout,"under-hyper.csv"), row.names = FALSE)
}

#! Filter drivers
drivers <- list()
cuttop <- 20
for (vs in names(stats)){
    oh <- stats[[vs]]$over.hypo
    filtro <- oh[(oh$gr1.cpgs.med>0.6) & (oh$gr2.cpgs.med<0.4),]
    # filtro <- filtro[order(-filtro$diffexpr), ]
    drivers[[vs]]$over.hypo <- filtro$gene[1:cuttop]

    uh <- stats[[vs]]$under.hyper
    filtro <- uh[(uh$gr1.cpgs.med<0.4) & (uh$gr2.cpgs.med>0.6),]
    drivers[[vs]]$under.hyper <- filtro$gene[1:cuttop]
}

#! Write out Excel file
wb <- createWorkbook()
for (vs in names(drivers)){
    print(vs)
    df <- data.frame(Over.Hypo = drivers[[vs]]$over.hypo, 
            Under.Hyper = drivers[[vs]]$under.hyper)
    addWorksheet(wb = wb, sheetName = vs)
    writeData(wb = wb, sheet = vs, x = df, borders = "n")
    print(df)
}
saveWorkbook(wb, "data/tables/drivers.xlsx", overwrite = TRUE)

#! delete NAs
df <- list()
for (vs in names(drivers)){
    v <- drivers[[vs]]$over.hypo
    df[[vs]]$over.hypo <- v[!is.na(v)]
    v <- drivers[[vs]]$under.hyper
    df[[vs]]$under.hyper <- v[!is.na(v)]
}
saveRDS(df, file = "data/RDS/drivers.rds")