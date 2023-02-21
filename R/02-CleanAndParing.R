require(TCGAbiolinks)
require(SummarizedExperiment)
source("R/libtools.R")

rnas.raw <- readRDS(file = "data/RDS/rnas-raw.rds")
print(dim(rnas.raw)) #56602 x 394
meths.raw <- readRDS(file = "data/RDS/meths-raw.rds")
print(dim(meths.raw)) #56602 x 394


#! RNA CLEAN SAMPLES

# clean duplicates in sample_patient_id
rnas <- rnas.raw[,!duplicated(rnas.raw$sample)]
print(dim(rnas)) #391

t <- c("Primary solid Tumor","Solid Tissue Normal")
rnas <- rnas[,rnas$definition %in% t]
t <- c("stage i","stage ii", "stage iii", "stage iv" )
rnas <- rnas[,rnas$tumor_stage %in% t]
rnas$tipo[rnas$tumor_stage=="stage i"] <- "Stage1"
rnas$tipo[rnas$tumor_stage=="stage ii"] <- "Stage2"
rnas$tipo[rnas$tumor_stage=="stage iii"] <- "Stage3"
rnas$tipo[rnas$tumor_stage=="stage iv"] <- "Stage4"
rnas$tipo[rnas$sample_type_id=="11"] <- "NT"
rnas$tipo <- as.factor(rnas$tipo)
dim(rnas) #388
print(table(rnas$tipo))


#! METHYLATION CLEAN SAMPLES
meths <- meths.raw[,!duplicated(meths.raw$sample)]


#! PARING

length(rnas$sample)
length(meths$sample)
common.samples <- intersect(rnas$sample,meths$sample)
length(common.samples)
meths <- meths[,match(common.samples,meths$sample)]
print(dim(meths))
rnas <- rnas[,match(common.samples,rnas$sample)]
meths$tipo <- rnas$tipo
#! print result
print(table(meths$tipo))


#! METHs CpGs CLEAN

print(dim(meths))    #! [1] 485577    342
meths <- subset(meths, subset = (rowSums(is.na(assay(meths))) == 0))
print(dim(meths))     #! [1] 382720    485


#! RNAs gene CLEAN

print(dim(rnas))    #! [1] 56602   342
threshold <- round(dim(rnas)[2]/2)
ridx <- rowSums(assay(rnas) == 0) <= threshold
rnas <- rnas[ridx, ]
print(dim(rnas))    #! [1] 31346   342
ridx <- rowMeans(assay(rnas)) >= 10
rnas <- rnas[ridx, ]
print(dim(rnas))    #! [1] 22199   342

#! Biomart
bm <- read.csv("data/tables/biomart.csv")
genes <- as.data.frame(rowData(rnas))
stopifnot(sum(duplicated(genes))==0)
genes <- merge(genes,bm)
# assay(rnas[rownames(rnas) == "ENSG00000148426",3])
rnas <- rnas[rownames(rnas) %in% genes$ensembl_gene_id, ]
rowData(rnas) <- genes
stopifnot(sum(duplicated(rownames(rnas)))==0)
dim(rnas)
# TODO: validar esto
# borrar <- data.frame(a=genes$ensembl_gene_id, b=rowData(d)$ensembl_gene_id)
# table(ifelse(borrar$a==borrar$b,"Yes","No"))


#! Saving
print("Saving clean-paring data...")
saveRDS(meths, file = "data/RDS/meths-clean.rds")
saveRDS(rnas, file = "data/RDS/rnas-clean.rds")
print(dim(meths))
print(dim(rnas))

