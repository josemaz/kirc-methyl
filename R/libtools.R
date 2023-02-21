require(NOISeq)
require(EDASeq)
require(ggplot2)
require(cowplot)
require(MASS)

####################################################################3
#!- Plot a PCA of SummarizedExperiment from a column (namecol)
mypca <- function(dat, namecol, fout = NULL){
    d2 <-  colData(dat)[namecol]
    mydat = NOISeq::readData( assay(dat) , factors = d2)
    # myPCA = dat(mydat, type = "PCA", norm=TRUE)
    myPCA = dat(mydat, type = "PCA", logtransf = F)
    if(!is.null(fout)){
        print(paste0("Writing in: ",fout))
        png(fout)
    }
    explo.plot(myPCA, factor = namecol, plottype = "scores")
    if(!is.null(fout))dev.off()
}

####################################################################3
#!- Normalization
normal.proc <- function(dat, namecol, do.pca = FALSE, expr.tipo="meth",
                                pngd="data/plots/qc/"){
  dir.create(pngd, recursive = TRUE)
  print(paste0("Using normalization method for: ",expr.tipo))
  if(do.pca){
    print("Plotting PCA before Norm")
    mypca( dat, namecol, fout = 
              paste0(pngd,"/PCA-",expr.tipo,"-BeforeNorm.png"))
  }
  fac <- colData(dat)[namecol]
  if(expr.tipo == "rna"){
      print("Starting RNA Normalization")
      # ! Pre Normalization
      ln.data <- withinLaneNormalization(assay(dat), 
                                          rowData(dat)$geneLength, which = "full")
      gcn.data <- withinLaneNormalization(ln.data , rowData(dat)$gcContent,
                                          which = "full")
      # norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
      norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
      noiseqData <- NOISeq::readData( norm.counts, factors = fac)
  } else if(expr.tipo == "meth") {
      print("Starting METHYLATION Normalization")
      noiseqData <- NOISeq::readData( assay(dat), factors = fac)
  }
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  #! Post Normalization
  assay(dat) <- exprs(mydata2corr1)  
  if(do.pca){
    print("Plotting PCA After Norm")
    mypca(dat,namecol,fout = paste0(pngd,
            "PCA-",expr.tipo,"-AfterNorm.png"))
  }
  return(dat)
}

####################################################################3
#!- Calculate all result tables of DEGs multigroup (todos vs. todos)
degs <- function(dats, form = ~ grupo){
    #! DEGs multigroup
    #! Ajuste de datos
    assay(dats) <- round(assay(dats)) # solo enteros
    print("[DEG] Genes duplicados en el gname: ")
    print(rownames(dats[duplicated(rownames(dats)),]))
    dats <- dats[!duplicated(rownames(dats)),]
    #! DEseq2
    ddsSE <- DESeqDataSet(dats, design = ~ grupo)
    dds <- DESeq(ddsSE)
    #! Combiantoria de los contrastes
    combi <- combn(levels(dats$grupo), 2)
    lres <- apply(combi,2, 
      function(x,dds1 = dds){
        print(green(c(x[2],x[1])))
        res <- results(dds1,contrast=c("grupo",x[2],x[1]))
        return(res)
      } 
    )
    names(lres) <- apply(combi,2, function(x) {return(paste0(x[2],"_",x[1]))} )
    # return(list(lres = lres, ups = ups, downs = downs, todos = todos))
    return(lres)
}

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                list(a = format(coef(m)[1], digits = 4),
                b = format(coef(m)[2], digits = 4),
                r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

scatplot <- function(dat, fpng, med1, med2){
  pmain <- ggplot(dat, aes(x=gcpgs, y=gexps)) + 
      # geom_point(size = 2) + 
      geom_point(aes(color=tipo),size = 2) + 
      geom_smooth(formula = y ~ x, method=lm, fullrange=TRUE) +
      geom_vline(xintercept = med1, linetype="dotted", 
                color = "red", size=1.5) +
      geom_vline(xintercept = med2, linetype="dotted", 
                color = "blue", size=1.5) +
      geom_text(x=0.3, y=max(dat$gexps), label = eq(dat$gcpgs,dat$gexps), parse = TRUE) +
      xlim(0, 1) #+ ylim(0, max(vals$gexps)) 
  # xbox <- axis_canvas(pmain, axis = "x", coord_flip = TRUE) + 
  #   geom_boxplot(data = dat, aes(x=gcpgs, y=gexps, color = tipo)) + coord_flip()
  xbox <- axis_canvas(pmain, axis = "x") + 
    geom_density(data = dat, aes(x = gcpgs, color = tipo))
  ybox <- axis_canvas(pmain, axis = "y") + 
    geom_density(data = dat, aes(y = gexps, color = tipo))
  p1 <- insert_xaxis_grob(pmain, xbox, position = "top")
  p2 <- insert_yaxis_grob(p1, ybox, position = "right")
  # ggMarginal(p1, type="density")
  # fpng = paste0(pngdir, "down-NT_Stage1-", gnom, ".png")
  ggsave(fpng, plot=p2, width = 8, height = 8)
}

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

mixdensity <- function(mixm, fpng) {
  # data.frame(x = mixm$x) %>%
  df <- data.frame(x=mixm$x)
  ggplot(df, aes(x=x)) +     
      # geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
      #               fill = "white") +
  geom_histogram(aes(x, ..density..), binwidth = .01, colour = "black", 
          fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
          args = list(mixm$mu[1], mixm$sigma[1], lam = mixm$lambda[1]),
          colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
          args = list(mixm$mu[2], mixm$sigma[2], lam = mixm$lambda[2]),
          colour = "blue", lwd = 1.5) +
  ylab("Density")
}

stats.gaus <- function(valores, gr1n = "NT", gr2n = "Stage1"){
  gr1 <- valores[valores['tipo'] == gr1n,]
  gr2 <- valores[valores['tipo'] == gr2n,]

  #! Ajuste lineal
  mlinear <- lm(gexps~gcpgs, data=valores)
  a <- coef(mlinear)[2] #pendiente del ajuste lineal
  # print(a)

  #! Gauss stats CPGS per group
  fit <- fitdistr(gr1$gcpgs, densfun="normal")
  # gr1.cpgs.mea1 = fit$estimate[1]
  gr1.cpgs.sig = fit$estimate[2]
  gr1.cpgs.mea = mean(gr1$gcpgs)
  gr1.cpgs.med = median(gr1$gcpgs)
  # print(c(gr1.cpgs.mea,gr1.cpgs.sig, gr1.cpgs.med ))
  fit <- fitdistr(gr2$gcpgs, densfun="normal")
  gr2.cpgs.sig = fit$estimate[2]
  gr2.cpgs.mea = mean(gr2$gcpgs)
  gr2.cpgs.med = median(gr2$gcpgs)
  # print(c(gr2.cpgs.mea,gr2.cpgs.sig, gr2.cpgs.med ))

  #! Gauss stats EXPRESION per group
  fit <- fitdistr(gr1$gexps, densfun="normal")
  gr1.gexps.sig = fit$estimate[2]
  gr1.gexps.mea = mean(gr1$gexps)
  gr1.gexps.med = median(gr1$gexps)
  # print(c(gr1.gexps.mea,gr1.gexps.sig, gr1.gexps.med ))
  fit <- fitdistr(gr2$gexps, densfun="normal")
  gr2.gexps.sig = fit$estimate[2]
  gr2.gexps.mea = mean(gr2$gexps)
  gr2.gexps.med = median(gr2$gexps)
  # print(c(gr2.gexps.mea,gr2.gexps.sig, gr2.gexps.med ))

  return(c(a,
            gr1.cpgs.mea,gr1.cpgs.sig, gr1.cpgs.med,
            gr2.cpgs.mea,gr2.cpgs.sig, gr2.cpgs.med,
            gr1.gexps.mea,gr1.gexps.sig, gr1.gexps.med,
            gr2.gexps.mea,gr2.gexps.sig, gr2.gexps.med ))
}


#! TODO: Function to plot individual gausian fit per group
# dist.filter <- function(val){
#   fit <- fitdistr(values1$gexps, densfun="normal")
#   prome <- fit$estimate[1]
#   sigma <- fit$estimate[2]
#   print(paste("mean:",prome,"sigma^2:",sigma2))
#   fit <- fitdistr(values1$gcpgs, densfun="normal")
#   prome <- fit$estimate[1]
#   sigma <- fit$estimate[2]
#   sigma2 <- sigma*sigma
#   print(paste("mean:",prome,"sigma^2:",sigma2))
  
#   # ggplot() + geom_density(data = values1, aes(x = gcpgs)) +
#   #   stat_function(geom = "line", fun = plot_mix_comps,
#   #       args = list(fit$estimate[1], fit$estimate[2], 1),
#   #       colour = "blue", lwd = 1.5) +
#   #   xlim(0, 1)

# }

        
        
        
        