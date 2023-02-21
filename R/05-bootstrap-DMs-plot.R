library(ggplot2)

df <- readRDS(file = "data/RDS/bootstrap-dms.rds")

dpng <- "data/plots/bootstrap-samples-dms/"
dir.create(dpng, recursive = TRUE)

df1 <- df[df$tipo == "hyper",c(2,3)]
ggplot(df1,aes(contrast,valor)) + 
    geom_boxplot() + coord_flip() + # scale_y_log10() +
    labs( x = "Contrast", y = "Hypermethylated CpGs")
ggsave(paste0(dpng,"hyper-CpGs.png"), width = 10, height = 10, units = "cm")
df1 <- df[df$tipo == "hypo",c(2,3)]
ggplot(df1,aes(contrast,valor)) + 
    geom_boxplot() + coord_flip() + # scale_y_log10() +
    labs( x = "Contrast", y = "Hypomethylated CpGs")
ggsave(paste0(dpng,"hypo-CpGs.png"), width = 10, height = 10, units = "cm")





    # df <- data.frame(hyper=unlist(dat$hyper),Group=rep(names(dat$hyper),each=niter))
    # ro <- rev(levels(as.factor(df$Group)))
    # df$Group <- factor(df$Group,ro)