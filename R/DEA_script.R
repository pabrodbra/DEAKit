### ------------------------
### DEA_script.R
### Author: Pablo Rodriguez Brazzarola & Guillermo Lopez Garcia
### Description: Script that executes the DEAKit workflow automatically, using input arguments management.
### Usage (from DEAKit directory): Rscript DEA_script.R -flag1 -flag2 ... (--help for details)
### ------------------------


### Load required optparse package for input arguments management
if(!"optparse" %in% installed.packages()) 
  install.packages("optparse", repos = "https://cran.rstudio.com")
require(optparse)

### Prepare input arguments
option_list <- list(
  make_option(c("-c", "--rcc-directory"),  action="store", type="character", default=file.path("data", "PanCancer Pathways - 57"), 
              dest="rcc_dir", help="RCC data directory path"),
  make_option(c("-r", "--rlf-file"), action="store", type="character", default=file.path("data", "PanCancer Pathways - 57", "NS_CancerPath_C2535.rlf"),
              dest="rlf_file", help="RLF file full path"),
  make_option(c("-k", "--key"), action="store", type="integer", default=1,
              dest="key", help="Key of interest"),
  make_option(c("-l", "--label"), action="store", type="character", default="PCR",
              dest="label", help="Label of key of interest"),
  make_option(c("-i", "--first-value"), action="store", type="character", default="0",
              dest="value_one", help="First value of key of interest"),
  make_option(c("-s", "--second-value"), action="store", type="character", default="1",
              dest="value_two", help="Second value of key of interest"),
  make_option(c("-f", "--fov-threshold"), action="store", type="double", default=80,
              dest="fov", help="Field Of View (FOV) threshold, used in FOV plot"),
  make_option(c("-b", "--bd-min"), action="store", type="double", default=0.05,
              dest="bd_min", help="Binding Density (BD) minimum threshold, used in BD plot"),
  make_option(c("-m", "--bd-max"), action="store", type = "double", default=2.25,
              dest="bd_max", help="Binding Density (BD) maximum threshold, used in BD plot"),
  make_option(c("-d", "--dea-pvalue"), action="store", type="double", default=0.05,
              dest="dea_pvalue", help="Differential Expression Analysis (DEA) p-value threshold"),
  make_option(c("-p", "--path-pvalue"), action="store", type="double", default=0.05,
              dest="path_pvalue", help="Pathway enrichment p-value threshold"),
  make_option(c("-g", "--log-fc"), action="store", type="double", default=1,
              dest="log_fc", help="Log Fold Change (FC) threshold"),
  make_option(c("-q", "--q-value"), action="store", type="double", default=1e-10,
              dest="q_value", help="Q-value threshold used in pathway enrichment"),
  make_option(c("-e", "--seed"), action="store", type="double", default=12345,
              dest="seed", help="Random seed for reproducibility"),
  make_option(c("-o", "--output"), action="store", type="character", default=file.path("output", ""),
              dest="output", help="Output directory path"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              dest="verbose", help="Activate verbose mode")
)

### Parse command line options
opt <- parse_args(OptionParser(option_list=option_list))


### Load required packages and functions for executing the script
source("R/DEA_functions.R")


### Arguments check

# Check path inputs
if(!file.exists(opt$rcc_dir)) {
  stop("An existent directory containing the RCC files must be specified")
}

if(!file.exists(opt$rlf_file)) {
  stop("An existent RLF file must be specified")
}

if(!file.exists(opt$output)) {
  stop("An existent output directory must be specified")
}

# Check positive numeric inputs
if(opt$fov < 0) {
  stop("A positive Field Of View (FOV) threshold must be specified")
}

if(opt$bd_min < 0 || opt$bd_max < 0) {
  stop("Positive Binding Density (BD) minimum and maximum thresholds must be specified")
}

if(opt$dea_pvalue <= 0 || opt$path_pvalue <= 0) {
  stop("Differential Expression Analysis (DEA) and pathway enrichment p-value thresholds must be greater than zero")
}

if(opt$log_fc < 0) {
  stop("A positive log Fold Change (FC) threshold must be specified")
}

if(opt$q_value <= 0) {
  stop("Q-value threshold used in pathway enrichment must be greater than zero")
}


### Prepare functions parameters from received arguments
set.seed(opt$seed)
keys.vector <- unlist(lapply(seq_len(opt$key), function(x) paste("key",x,sep = "")), use.names = FALSE)
key.label <- tail(keys.vector, 1)
LOD.THRESHOLD <- (100-opt$fov)/100 # Dependant of FOV.THRESHOLD
BASE.OUTPUT.NAME <- paste(opt$output, opt$label, sep = "")
OUTPUT.NAME <- paste(BASE.OUTPUT.NAME, opt$value_one, opt$value_two, sep = "-")



### ------------------------
### Load RCC files. Build counts matrix
### ------------------------

# Select files according to the key of interest and its values

if(opt$verbose) {
  message("--------------- Verbose: Loading RCC files ---------------")
}

# Load data + Extract metadata
metadata <- extract.rcc.metadata(opt$rcc_dir, keys.vector, c(opt$value_one, opt$value_two))

# Retrieve Set and Count matrix
rcc.set.and.count.matrix <- get.RCC.set.and.counts(metadata, opt$rcc_dir, opt$rlf_file)
eset <- rcc.set.and.count.matrix$set
counts <- rcc.set.and.count.matrix$count.matrix



### ------------------------
### Quality Control
### ------------------------

if(opt$verbose) {
  message("--------------- Verbose: Performing Quality Control ---------------")
}


# Field of View (FOV) Plot
fov.plot <- plotFOV(eset = eset, metadata = metadata, fov.threshold = opt$fov,
                    comparison.key = key.label, legend.label = opt$label)
pdf(file.path("output/", "FOV.pdf"))
plot(fov.plot)
invisible(dev.off())

# Binding Density (BD) Plots
bd.plot <- plotBD(eset = eset, metadata = metadata, y.thresholds = c(opt$bd_min, opt$bd_max),
                  comparison.key = key.label, legend.label = opt$label)
pdf(file.path("output/", "BD.pdf"))
plot(bd.plot)
invisible(dev.off())

# Positive Controls
pdf(file.path("output/", "Positive-Control.pdf"))
plot(boxplot.expr(eset, is.positive))
invisible(dev.off())

# Negative Controls
pdf(file.path("output/", "Negative-Control.pdf"))
plot(boxplot.expr(eset, is.negative))
invisible(dev.off())


# Noise Threshold
lodcounts <- extract.pred(eset, is.negative)
lod <- mean(lodcounts$count) + 2 * sd(lodcounts$count)

# Housekeeping Genes
housekeeping.boxplot<-boxplot.expr(eset,is.housekeeping)
pdf(file.path("output/", "Housekeeping-Genes.pdf"))
plot(housekeeping.boxplot+geom_hline(yintercept = (lod),colour="red"))
invisible(dev.off())

# Expression of all the housekeeping genes in each sample
housekeeping.condition.boxplot <- boxplot.multiple.conditions(counts, conditions = c("Housekeeping"), 
                                                              title.label = "Housekeeping Genes Expression per Sample")
pdf(file.path("output/", "Housekeeping-Samples.pdf"))
plot(housekeeping.condition.boxplot+geom_hline(yintercept = (lod),colour="red"))
invisible(dev.off())

# Expression of multiple condition genes in each sample
endogenous.housekeeping.condition.boxplot <- boxplot.multiple.conditions(counts, conditions = c("Endogenous", "Housekeeping"), 
                                                                         title.label = "Housekeeping + Endogenous Genes Expression per Sample")
pdf(file.path("output/", "Housekeeping-Endogenous.pdf"))
plot(endogenous.housekeeping.condition.boxplot+geom_hline(yintercept = (lod),colour="red"))
invisible(dev.off())



### ------------------------
### Normalization
### ------------------------

if(opt$verbose) {
  message("--------------- Verbose: Performing normalization ---------------")
}


# Pre-Normalization
metadata$pos_nf = pos.factor(eset) #  normalization
counts <- counts[!grepl("Positive", rownames(counts)),]
counts <- counts[!grepl("Negative", rownames(counts)),]

# Plot pre-normalization
pre.normalization.boxplot <- boxplot.count.matrix(counts, "Pre-Normalization with Positive controls")
pdf(file.path("output/", "Pre-Normalization.pdf"))
plot(pre.normalization.boxplot + geom_hline(yintercept = lod,colour="red"))
invisible(dev.off())


# Normalize using positive controls
ncounts <- counts %*% diag(metadata$pos_nf) 
colnames(ncounts) <- colnames(counts)

# Plot post-normalization
post.normalization.boxplot <- boxplot.count.matrix(ncounts, "Post-Normalization with Positive controls")
pdf(file.path("output/", "Post-Normalization.pdf"))
plot(post.normalization.boxplot + geom_hline(yintercept = lod,colour="red"))
invisible(dev.off())


# Normalize using housekeeping
metadata$hk_nf <- hk.factor(ncounts, lod)
ncounts <- ncounts %*% diag(metadata$hk_nf)
colnames(ncounts) <- colnames(counts)

# Plot normalized housekeeping
housekeeping.post.normalization.boxplot <- boxplot.count.matrix(ncounts, "Post-Normalization with Housekeeping controls")
pdf(file.path("output/", "Housekeeping-Normalization.pdf"))
plot(housekeeping.post.normalization.boxplot + geom_hline(yintercept = lod,colour="red") + geom_smooth(se=T, aes(group=1)))
invisible(dev.off())


# Drop genes below threshold
all.names <- rownames(ncounts)
ncounts <- ncounts[(rowSums(ncounts < lod) < round((LOD.THRESHOLD * ncol(ncounts)),0)),]
filtered.names <- rownames(ncounts)
filter.out.genes <- all.names[all.names %!in% filtered.names]

# Save metadata dataframe and the normalised counts matrix
RDS.METADATA.OUTPUT <- paste(OUTPUT.NAME, "BMI", "metadata.rds", sep = "_")
RDS.NORMALIZED.OUTPUT <- paste(OUTPUT.NAME, "BMI", "HK-normalized-counts.rds", sep = "_")
saveRDS(metadata, RDS.METADATA.OUTPUT) 
saveRDS(ncounts, RDS.NORMALIZED.OUTPUT)


# PCA using variance
pcdata <- scale(ncounts, center = TRUE, scale = TRUE)
pc <- pca.loadings(pcdata, 50)
comps <- data.frame(pc$x)
comps$Name <- rownames(comps)
comps <- comps %>% left_join(metadata, by = c(Name = "Sample.name"))

# Pass labels to legend and/or title
# Select principal components to plot
PRINCIPAL.COMPONENT.X <- 1
PRINCIPAL.COMPONENT.Y <- 2

PCA.PLOT.TITLE <- paste("PCA Plot | Component", PRINCIPAL.COMPONENT.X, "-", PRINCIPAL.COMPONENT.Y, sep = " ")

pdf(file.path("output/", "Samples-PCA.pdf"))
plot(plotPCA(comps, pc, PRINCIPAL.COMPONENT.X, PRINCIPAL.COMPONENT.Y, key.label, 
             legend.label =  opt$label, plot.title = PCA.PLOT.TITLE))
invisible(dev.off())



### ------------------------
### Differential expression analysis
### ------------------------

if(opt$verbose) {
  message("--------------- Verbose: Performing Differential Expression Analysis ---------------")
}

design <- as.formula(paste("~", key.label, sep = ""))
rownames(metadata) <- metadata$Sample.name

dea.o <- getDE(ncounts = round(ncounts), metadata = metadata, design = design)
dea.p <- dea.o %>% dplyr::filter(padj < opt$dea_pvalue) %>% arrange(-log2FoldChange)

# Save DEA Results
P.VALUE.STRING <- paste("PVALUE", opt$dea_pvalue, sep = "-")
DEA.OUTPUT.NAME <- paste(OUTPUT.NAME, "DEA.csv", sep = "_")
DEA.SIGNIFICATIVE.OUTPUT.NAME <- paste(OUTPUT.NAME, P.VALUE.STRING, "DEA.csv", sep = "_")

write_csv(dea.o, DEA.OUTPUT.NAME)
write_csv(dea.p, DEA.SIGNIFICATIVE.OUTPUT.NAME)

# Volcano Plot (first version)
tab = data.frame(logFC = dea.o$log2FoldChange, negLogPval = -log10(dea.o$padj))
tab2 = data.frame(logFC = dea.o$log2FoldChange, negLogPval = -log10(dea.o$padj), Gene=dea.o$gene.name)

pdf(file.path("output/", "Volcano.pdf"))
plot.Volcano(tab, tab2, opt$log_fc, opt$dea_pvalue, plot.title = "Volcano Plot")
invisible(dev.off())

# Volcano Plot (second version)
pdf(file.path("output/", "Volcano_version2.pdf"))
plot.Volcano2(tab2, opt$log_fc, opt$dea_pvalue, plot.title = "Volcano Plot")
invisible(dev.off())



### ------------------------
### Pathway analysis
### ------------------------

if(opt$verbose) {
  message("--------------- Verbose: Performing Pathway Analysis ---------------")
}


# The core script for the pathway analysis expects a dataframe with the DEA named res
res <- getDE.raw(ncounts = round(ncounts), metadata = metadata, design = design) 

# Setup ENTREZ - Takes time
mart <- biomaRt::useMart(biomart = "ensembl", dataset = biomart_dataset)
entrez <- biomaRt::getBM(attributes = c("refseq_mrna", "entrezgene"), mart = mart)
entrez$entrezgene <- as.character(entrez$entrezgene)
entrezsymbol <- biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), mart = mart)
entrezsymbol$entrezgene <- as.character(entrezsymbol$entrezgene)

converted <- unlist(lapply(strsplit(res$rowname, "_", fixed = TRUE), "[", 4))
converted <- unlist(lapply(strsplit(converted, ".", fixed = TRUE), "[", 1))
converted <- paste0("NM_", converted)
res$refseq_mrna <- converted
lost.genes <- length(res$rowname) - length((res %>% data.frame() %>% left_join(entrez, by = "refseq_mrna") %>% 
                                              filter(!is.na(entrezgene)) %>% filter(!is.na(log2FoldChange)) %>% filter(!is.na(lfcSE)))$rowname)

### Enrichments 
ENRICHMENT.ALL.OUTPUT <- paste(OUTPUT.NAME, "Enrichment-ALL.csv", sep = "_")
ENRICHMENT.OVER.OUTPUT <- paste(OUTPUT.NAME, "Enrichment-OVER.csv", sep = "_")
ENRICHMENT.UNDER.OUTPUT <- paste(OUTPUT.NAME, "Enrichment-UNDER.csv", sep = "_")

# Enrich - All
enrich.rs <- enrich.cp(res, entrez, opt$label, type="all", pval.threshold = opt$path_pvalue, lfc.threshold = opt$log_fc)
enrich.rs.summary <- enrich.rs$summary %>% arrange(p.adjust)
enrich.rs.summary <- convert.enriched.ids(enrich.rs.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)

# Enrich - Over
enrich.rs.over <- enrich.cp(res, entrez, opt$label, type="over", pval.threshold = opt$path_pvalue, lfc.threshold = opt$log_fc)
enrich.rs.over.summary <- enrich.rs.over$summary %>% arrange(p.adjust)
enrich.rs.over.summary <- convert.enriched.ids(enrich.rs.over.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)

# Enrich - Under
enrich.rs.under <- enrich.cp(res, entrez, opt$label, type="under", pval.threshold = opt$path_pvalue, lfc.threshold = opt$log_fc)
enrich.rs.under.summary <- enrich.rs.under$summary %>% arrange(p.adjust)
enrich.rs.under.summary <- convert.enriched.ids(enrich.rs.under.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)

write_csv(x = enrich.rs.summary, path = ENRICHMENT.ALL.OUTPUT)
write_csv(x = enrich.rs.over.summary, path = ENRICHMENT.OVER.OUTPUT)
write_csv(x = enrich.rs.under.summary, path = ENRICHMENT.UNDER.OUTPUT)


### Dotplots
DOTPLOT.ENRICHMENT.ALL.TITLE <- "Dotplot of All Enriched Pathways after DEA"
DOTPLOT.ENRICHMENT.OVER.TITLE <- "Dotplot of Overexpressed Enriched, Pathways after DEA"
DOTPLOT.ENRICHMENT.UNDER.TITLE <- "Dotplot of Underexpressed Enriched Pathways after DEA"

# KEGG enrichment in all significantly expressed genes:
pdf(file.path("output/", "Dotplot-All.pdf"))
plot(dotplot(enrich.rs$kg, x="count", showCategory=10, color="qvalue", title = DOTPLOT.ENRICHMENT.ALL.TITLE))
invisible(dev.off())

# KEGG enrichment in over-expressed genes:
pdf(file.path("output/", "Dotplot-Over.pdf"))
plot(dotplot(enrich.rs.over$kg, x="count", showCategory=10, color="qvalue", title = DOTPLOT.ENRICHMENT.OVER.TITLE))
invisible(dev.off())

# KEGG enrichment in under-expressed genes
pdf(file.path("output/", "Dotplot-Under.pdf"))
plot(dotplot(enrich.rs.under$kg, x="count", showCategory=10, color="qvalue", title = DOTPLOT.ENRICHMENT.UNDER.TITLE))
invisible(dev.off())


# Filter by qvalue and by ontology (kegg for pathview)
enrich.kegg.all.summary.filtered <- ontology.enrichment.qval.filter(enrichment.summary = enrich.rs.summary,
                                                                    ontology = "kg",
                                                                    qval = opt$q_value)

# View KEGG Path
# Move to images directory to save .png in that directory
setwd("images")
for(pathid in enrich.kegg.all.summary.filtered$ID){
  viewkeggpath(path=pathid, enrichment = enrich.rs.summary, dea.p = dea.p, output = ".")
}
# Move back to the parent directory
setwd("..")



### ------------------------
### R Markdown Report --- FIGURE OUT TEXT AND PARAMS WANTED FOR THIS REPORT...
### ------------------------

if(opt$verbose) {
  message("--------------- Verbose: Creating R Markdown Report ---------------")
}


# Make all plots as variables, no need to recalculate...
MARKDOWN.PATH <- "reports/DEA_report.Rmd"
MARKDOWN.OUTPUT <- "../docs/DEA_report.html"
markdown.params <- list(
  title = paste("DEAKit Report - Effect of", opt$label, "in gene expression", sep = " "),
  author = "PERSON USING THE TOOLKIT",
  date = Sys.Date(),
  sample_type = "breast cancer",
  sample_associate = opt$label
)
rmarkdown::render("reports/DEA_report.Rmd", output_file = MARKDOWN.OUTPUT)
