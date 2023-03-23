
ETAPAS = ["Stage1","Stage2","Stage3","Stage4"]
ALLGR = ["NT"] + ETAPAS
CUTSMI = ["10k","50k","100k"]

print(ALLGR)

# GDC deprecated this data
# rule download:
#     input:
#     output:
#         "data/RDS/meths-raw.rds",
#         "data/RDS/rnas-raw.rds"
#     shell:
#         "Rscript R/01-Download.R"

# rule biomart:
#     input:
#     output: "data/tables/biomart.csv"
#     shell: "Rscript R/biomart.R"

rule clean:
    input:
        "data/RDS/rnas-raw.rds",
        "data/RDS/meths-raw.rds",
        "data/tables/biomart.csv"
    output: 
        "data/RDS/rnas-clean.rds",
        "data/RDS/meths-clean.rds"
    shell:
        "Rscript R/02-CleanAndParing.R"

rule normalization:
    input: "data/RDS/rnas-clean.rds"
    output:
        "data/RDS/rnas-norm.rds",
        "data/plots/qc/PCA-rna-AfterNorm.png",
        "data/plots/qc/PCA-rna-BeforeNorm.png",
        expand("data/tables/expression/{etapa}.tsv", etapa=ALLGR)
    shell:
        "Rscript R/03-Normalization.R"

rule bootstrap_dms:
    input: "data/RDS/meths-clean.rds"
    output: "data/RDS/bootstrap-dms.rds"
    log: "logs/salida.log"
    shell: "Rscript R/05-bootstrap-DMs.R &> logs/salida.log"

rule bootstrap_dms_plot:
    input: "data/RDS/bootstrap-dms.rds"
    output: "data/plots/bootstrap-samples-dms/hyper-CpGs.png",
        "data/plots/bootstrap-samples-dms/hypo-CpGs.png"
    log: "logs/salida.log"
    shell: "Rscript R/05-bootstrap-DMs-plot.R"

rule heatmaps_mean:
    input: "data/RDS/meths-clean.rds"
    output: "data/plots/meth-stats/heatmap-100.pdf",
        "data/plots/meth-stats/meandesc.png"
    shell: "Rscript R/06-meth-plots.R"

rule degs:
    input: "data/RDS/rnas-norm.rds"
    output: 
        "data/RDS/degs.rds",
        "data/RDS/gdegs.rds",
        "data/plots/rna-volcanos/NTvsPT_rna_volcano.png"
    shell:
        "Rscript R/04-DEGs.R"

rule dms:
    input: 
        "data/RDS/rnas-norm.rds",
        "data/RDS/meths-clean.rds",
        "data/RDS/gdegs.rds"
    output: 
        "data/RDS/dms.rds",
        # "data/plots/meth-volcanos/NTvsStage1.png"
        expand("data/plots/meth-volcanos/NTvs{etapa}.png", etapa=ETAPAS)
    shell:
        "Rscript R/07-DMs.R"

rule nominees:
    input: 
        "data/RDS/meths-clean.rds",
        "data/RDS/dms.rds", # CpGs differenttialy Methylated
        "data/RDS/gdegs.rds"
    output:
        "data/RDS/nominees.rds",
        #! Only one example
        "data/tables/nominees/NT_Stage1-overhypo.tsv"
    shell:
        "Rscript R/08-Nominees.R"

rule nominees_stats:
    input: 
        "data/RDS/rnas-norm.rds",
        "data/RDS/meths-clean.rds",
        "data/RDS/nominees.rds"
    output:
        #! Also scats were created here
        "data/plots/scats/NT_Stage1/over/ABCA1.png",
        "data/RDS/stats.rds"
    shell:
        "Rscript R/09-Nominees-scats.R"

rule get_drivers:
    input:
        "data/RDS/stats.rds"
    output:
        "data/tables/drivers.xlsx",
        "data/RDS/drivers.rds"
    shell:
        "Rscript R/10-drivers.R"

rule venn_drivers:
    input:
        "data/RDS/drivers.rds"
    output:
        "data/plots/venn/over-hypo.png",
        "data/plots/venn/under-hyper.png"
    shell:
        "Rscript R/11-venn-drivers.R"

rule onco_tsg_db:
    input:
        "data/RDS/rnas-norm.rds"
    output:
        "data/tables/onco-tsg-db.tsv"
    shell:
        "Rscript R/12-onco_tsg_db.R"

rule nets_enrich:
    input:
        expand("data/tables/MI/{etapa}-{cut}.sif", etapa=ALLGR, cut=CUTSMI),
        "data/tables/onco-tsg-db.tsv"
    output:
        # Only a hack
        expand("data/tables/nets/{cut}/Over-Hypo-Stage1-net.tsv", cut=CUTSMI)
    run:
        for c in CUTSMI:
             shell("Rscript R/13-nets.R " + c)

# rule paper:
#     input: "docs/paper.Rmd"
#     output: "docs/paper.pdf"        
#     shell: "cd docs; Rscript -e \"rmarkdown::render('paper.Rmd')\""

rule all:
    input:
        "data/plots/bootstrap-samples-dms/hyper-CpGs.png", 
        "data/plots/meth-stats/meandesc.png", 
        expand("data/plots/meth-volcanos/NTvs{etapa}.png", etapa=ETAPAS),
        expand("data/tables/nets/{cut}/Over-Hypo-Stage1-net.tsv", cut=CUTSMI),
        "data/tables/onco-tsg-db.tsv"