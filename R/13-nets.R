# install.packages("23")
# install.packages("gprofiler2")

require(crayon)
require(igraph)
require(gprofiler2)
require(SummarizedExperiment)
require(openxlsx)
require(RColorBrewer)
coul  <- brewer.pal(3, "Set1")

args<-commandArgs(TRUE)
stopifnot(!is.na(args[1]))
cutmi = args[1]
cat(cyan(cutmi),"\n")
# cutmi = "10k"


cat(cyan("Read RNAs and transform Gene INFO\n"))
rnas <- readRDS("data/RDS/rnas-norm.rds")
df <- as.data.frame(rowData(rnas))
gdata <- as.data.frame(df[,"gene_name"])
colnames(gdata) <- c("gene_name")
rownames(gdata) <- rownames(df)

cat(cyan("Read mi-network files\n"))
fnames <- c("NT","Stage1","Stage2","Stage3","Stage4")
nets.mi <- list()
for (fname in fnames){
  f <- paste0("data/tables/MI/",fname,"-",cutmi,".sif") 
  expr <- read.table(f,sep=' ')
  colnames(expr) <- c("src","mi","dst")
  df2 <- merge(expr,gdata,by.x="src", by.y=0)
  df2 <- merge(df2,gdata,by.x="dst", by.y=0)
  df2 <- df2[,c(4,5,3)]
  colnames(df2) <- c("src","dst","mi")
  nets.mi[[fname]] <- df2
}

#! Read drivers
drivers <- readRDS("data/RDS/drivers.rds")
#! Read db with oncogenes and tumor suppressors
onco.tsg.db <- read.table(file = "data/tables/onco-tsg-db.tsv", 
      sep = '\t',header = TRUE)


#! START MAIN

plot.redes <- function(g,fpath){
  png(fpath)
  cat(cyan("[plot.redes] Writing: ", fpath,"\n"))  
  #! Make a palette of 3 colors
  my_color <- coul[as.numeric(as.factor(V(g)$typo))]
  plot(g, edge.arrow.size = 0.2, vertex.color=my_color,
            vertex.size=8, vertex.label.cex=0.9,
            vertex.label.dist=2.3, main=nom)
  legend("bottomleft", legend=levels(as.factor(V(g)$typo)), 
         col = coul , bty = "n", pch=20 , 
         text.col=coul , horiz = FALSE, inset = c(0.0, 0.0))
  dev.off()
}

#! Plotting and creation of networks by stage and specifications
d <- paste0("data/plots/nets-methdriven/",cutmi)
dir.create(d,recursive = TRUE)
dtablenets <- paste0("data/tables/nets/",cutmi)
dir.create(dtablenets,recursive = TRUE)
nets <- list()
etapas <- c("Stage1","Stage2","Stage3","Stage4")
for (i in etapas){
  net <- nets.mi[[i]]
  genes <- drivers[[paste0("NT_",i)]]$over.hypo
  # print(genes)
  cat(green("=============  Over Hypo: ",i,"\n"))
  df <- rbind( net[net$src %in% genes,], net[net$dst %in% genes,])
  nom <- paste0(dtablenets,"/Over-Hypo-",i,"-net.tsv")
  cat(cyan("[main] Writing: ", nom,"\n"))
  write.table(df,file = nom, sep = '\t', row.names = FALSE, quote = FALSE)
  g1 <- graph_from_data_frame(df, directed=TRUE, vertices=NULL)
  V(g1)$typo <- onco.tsg.db[match(V(g1)$name,onco.tsg.db$gname),]$typo
  nom <- paste0(d,"/Over-Hypo-",i,'.png')
  plot.redes(g1,nom)
  cat(green("=============  Under Hyper: ",i,"\n"))
  genes <- drivers[[paste0("NT_",i)]]$under.hyper
  df <- rbind( net[net$src %in% genes,], net[net$dst %in% genes,])
  nom <- paste0(dtablenets,"/Under-Hyper-",i,"-net.tsv")
  cat(cyan("[main] Writing: ", nom,"\n"))
  write.table(df,file = nom, sep = '\t', row.names = FALSE, quote = FALSE)
  g2 <- graph_from_data_frame(df, directed=FALSE, vertices=NULL)
  V(g2)$typo <- onco.tsg.db[match(V(g2)$name,onco.tsg.db$gname),]$typo
  nom <- paste0(d,"/Under-Hyper-",i,'.png')
  plot.redes(g2,nom)
  nets[[i]] <- list( over.hypo = g1, under.hyper = g2)
}


wb1 <- createWorkbook()
wb2 <- createWorkbook()
for (e in names(nets)){
  cat(green("=============  Over Hypo: ",e,"\n"))
  # cl <- components(g)
  # comps <- lapply(seq_along(cl$csize)[cl$csize > 1], function(x) 
  #                                     V(g)$name[cl$membership %in% x])
  # # print(comps[1])
  # lapply(comps,function(x){
  #   print(x)
  #   gostres <- gost(query = x, organism = "hsapiens", sources="GO")
  #   print(head(gostres$result[,c("p_value","term_name","intersection_size")]))
  # })
  g <- nets[[e]]$over.hypo
  gostres <- gost(query = names(V(g)), organism = "hsapiens", 
                    sources=c("GO","HPA"))
  addWorksheet(wb = wb1, sheetName = e)
  df <- gostres$result[,c("p_value","term_name","intersection_size","source")]
  writeData(wb = wb1, sheet = e, x = df, borders = "n")

  cat(green("=============  Under Hyper: ",e,"\n"))
  g <- nets[[e]]$under.hyper
  gostres <- gost(query = names(V(g)), organism = "hsapiens", 
                sources=c("GO","HPA"))
  addWorksheet(wb = wb2, sheetName = e)
  df <- gostres$result[,c("p_value","term_name","intersection_size","source")]
  writeData(wb = wb2, sheet = e, x = df, borders = "n")
}

d <- "data/tables/enrichment"
dir.create(d, recursive = TRUE)
saveWorkbook(wb1, paste0(d,"/enrich-over-hypo.xlsx"), overwrite = TRUE)
saveWorkbook(wb2, paste0(d,"/enrich-under-hyper.xlsx"), overwrite = TRUE)






# for (i in etapas){
#   print(paste("=============  Over Hypo: ",i))
#   net <- nets.mi[[i]]
#   genes <- drivers[[paste0("NT_",i)]]$over.hypo
#   # print(genes)
#   df <- rbind( net[net$src %in% genes,], net[net$dst %in% genes,])
#   nom <- paste0(dtablenets,"/Over-Hypo-",i,"-net.tsv")
#   print(nom)
#   write.table(df,file = nom, sep = '\t', row.names = FALSE, quote = FALSE)
#   g <- graph_from_data_frame(df, directed=TRUE, vertices=NULL)
#   # png(paste0(d,"/Over-Hypo-",i,'.png'))
#   # V(g)$typo <- onco.tsg.db[match(V(g)$name,onco.tsg.db$gname),]$typo
#   # #! Make a palette of 3 colors
#   # my_color <- coul[as.numeric(as.factor(V(g)$typo))]
#   # plot(g, edge.arrow.size = 0.2, vertex.color=my_color,
#   #           vertex.size=8, vertex.label.cex=0.9,
#   #           vertex.label.dist=2.3, main=nom)
#   # legend("bottomleft", legend=levels(as.factor(V(g)$typo)), 
#   #        col = coul , bty = "n", pch=20 , 
#   #        text.col=coul , horiz = FALSE, inset = c(0.0, 0.0))
#   # dev.off()
#   cat(green("##################################\n"))
#   cat(names(as.list(V(g))),"\n")
#   # cl <- components(g)
#   # comps <- lapply(seq_along(cl$csize)[cl$csize > 1], function(x) 
#   #                                     V(g)$name[cl$membership %in% x])
#   # # print(comps[1])
#   # lapply(comps,function(x){
#   #   print(x)
#   #   gostres <- gost(query = x, organism = "hsapiens", sources="GO")
#   #   print(head(gostres$result[,c("p_value","term_name","intersection_size")]))
#   # })
  
#   print(paste("=============  Under Hyper: ",i))
#   genes <- drivers[[paste0("NT_",i)]]$under.hyper
#   df <- rbind( net[net$src %in% genes,], net[net$dst %in% genes,])
#   nom <- paste0(dtablenets,"/Under-Hyper-",i,"-net.tsv")
#   write.table(df,file = nom, sep = '\t', row.names = FALSE, quote = FALSE)
#   g <- graph_from_data_frame(df, directed=FALSE, vertices=NULL)
#   # nom <- paste0(d,"/Under-Hyper-",i)
#   # png(paste0(nom,'.png'))
#   # V(g)$typo <- onco.tsg.db[match(V(g)$name,onco.tsg.db$gname),]$typo
#   # # Make a palette of 3 colors
#   # my_color <- coul[as.numeric(as.factor(V(g)$typo))]
#   # plot(g, edge.arrow.size = 0.2, vertex.color=my_color,
#   #               vertex.size=8, vertex.label.cex=0.9,
#   #               vertex.label.dist=2.3, main=nom)
#   # legend("bottomleft", legend=levels(as.factor(V(g)$typo)), 
#   #        col = coul , bty = "n", pch=20 , 
#   #        text.col=coul , horiz = FALSE, inset = c(0.0, 0.0))
#   # dev.off()
#   cat(blue("##################################\n"))
#   cat(names(as.list(V(g))),"\n")
#   # cl <- components(g)
#   # comps <- lapply(seq_along(cl$csize)[cl$csize > 1], function(x)
#   #   V(g)$name[cl$membership %in% x])
#   # lapply(comps,function(x){
#   #   print(x)
#   #   gostres <- gost(query = x, organism = "hsapiens", sources="GO")
#   #   print(head(gostres$result[,c("p_value","term_name","intersection_size")]))
#   # })

# }