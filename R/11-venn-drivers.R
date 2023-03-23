library(ggvenn)

drivers <- readRDS("data/RDS/drivers.rds")


# over.hypo <- lapply(drivers, function(x) {
#     return(x$over.hypo)
# })
oh <- list()
uh <- list()
for (vs in names(drivers)){
#     oh[[vs]] <- lapply(drivers[[vs]]$over.hypo, function(x) x[!is.na(x)])
    oh[[vs]] <- drivers[[vs]]$over.hypo
    # under.hyper[[vs]] <- drivers[[vs]]$under.hyper
    uh[[vs]] <- drivers[[vs]]$under.hyper
}
# print(oh)
# print(uh)

dout <- "data/plots/venn/"
dir.create(dout,recursive = TRUE)
p <- ggvenn(oh,  show_percentage = FALSE, 
       show_elements = TRUE, label_sep = "\n", text_size = 3,
       set_name_size = 3
)
ggsave(paste0(dout,"over-hypo.png"),plot = print(p), dpi=300)
p <- ggvenn(uh,  show_percentage = FALSE, 
       show_elements = TRUE, label_sep = "\n", text_size = 3,
       set_name_size = 3
)
ggsave(paste0(dout,"under-hyper.png"),plot = print(p), dpi=300)