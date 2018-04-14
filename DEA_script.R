### ------------------------
### DEA_script.R
### Author: Pablo Rodríguez Brazzarola
### ------------------------

rm(list = ls())
dev.off()

### Load DEA_functions.R
source("DEA_functions.R")

### Argument check (disabled by now)
# PATHS
data.directory <- "data" 
rcc.directory <- file.path(data.directory, "PanCancer Pathways - 57")
rlf.filename <- "NS_CancerPath_C2535.rlf"

# CONSTANTS
KEY.OF.INTEREST <- 1        # Param
LABEL.OF.INTEREST <- "PCR"  # Param
VALUES.OF.KEY.OF.INTEREST <- c("0","1")
SEED <- 12345
FOV.THRESHOLD <- 80 
BD.THRESHOLDS <- c(0.05, 2.25)
P-VALUE.THRESHOLD <- 0.05
LOG.FC.THRESHOLD <- 1
Q.VALUE.THRESHOLD <- 1e-10
# PARAM DEPENDANT
LOD.THRESHOLD <- (100-FOV.THRESHOLD)/100

# SETUP
set.seed(SEED)
keys.vector <- unlist(lapply(seq_len(KEY.OF.INTEREST), function(x) paste("key",x,sep = "")), use.names = FALSE)
key.label <- tail(keys.vector, 1)
OUTPUT.NAME <- paste(LABEL.OF.INTEREST, VALUES.OF.KEY.OF.INTEREST[1], VALUES.OF.KEY.OF.INTEREST[2], sep = "-")
### ------------------------
### Load RCC files. Build counts matrix
### ------------------------
# Seleccion de archivos en base a la llave de interés, y los valores

# Load data + Extract metadata
metadata <- extract.rcc.metadata(rcc.directory, keys.vector, VALUES.OF.KEY.OF.INTEREST)

# Retrieve Set and Count matrix
rcc.set.and.count.matrix <- get.RCC.set.and.counts(metadata, rcc.directory, rlf.filename)
eset <- rcc.set.and.count.matrix$set
counts <- rcc.set.and.count.matrix$count.matrix

### ------------------------
### Quality Control
### ------------------------

# Field of View (FOV) Plots
plotFOV(eset = eset, metadata = metadata, fov_threshold = FOV.THRESHOLD,
        comparison.key = key.label, legend.label = LABEL.OF.INTEREST)

# Binding Density (BD) Plots
plotBD(eset = eset, metadata = metadata,
       comparison.key = key.label, legend.label = LABEL.OF.INTEREST)

# Positive Controls
boxplot.expr(eset,is.positive)

# Negative Controls
boxplot.expr(eset,is.negative)

# Noise Threshold
lodcounts <- extract.pred(eset, is.negative)
lod <- mean(lodcounts$count) + 2 * sd(lodcounts$count)

# Housekeeping Genes
housekeeping.boxplot<-boxplot.expr(eset,is.housekeeping)
housekeeping.boxplot+geom_hline(yintercept = (lod),colour="red")

# Expression of all the housekeeping genes in each sample
housekeeping.condition.boxplot <- boxplot.multiple.conditions(counts, conditions = c("Housekeeping"), 
                                                            title.label = "Housekeeping Genes Expression per Sample")
housekeeping.condition.boxplot+geom_hline(yintercept = (lod),colour="red")

# Expression of multiple condition genes in each sample
endogenous.housekeeping.condition.boxplot <- boxplot.multiple.conditions(counts, conditions = c("Endogenous", "Housekeeping"), 
                                                          title.label = "Housekeeping + Endogenous Genes Expression per Sample")
endogenous.housekeeping.condition.boxplot+geom_hline(yintercept = (lod),colour="red")

### ------------------------
### Normalization
### ------------------------
counts.old <- counts
metadata.old <- metadata
eset.old <- eset

# Pre-Normalization
metadata$pos_nf = pos.factor(eset) # Positive Control normalization
counts <- counts[!grepl("Positive", rownames(counts)),]
counts <- counts[!grepl("Negative", rownames(counts)),]

# Plot pre-normalization
pre.normalization.boxplot <- boxplot.count.matrix(counts, "Pre-Normalization with Positive controls")
pre.normalization.boxplot + geom_hline(yintercept = lod,colour="red")

# Normalize
ncounts <- counts %*% diag(metadata$pos_nf) # Positive Control normalization
colnames(ncounts) <- colnames(counts)

# Plot post-normalization
post.normalization.boxplot <- boxplot.count.matrix(ncounts, "Post-Normalization with Positive controls")
post.normalization.boxplot + geom_hline(yintercept = lod,colour="red")
# ^^ No noticeable change ^^ #

# Housekeeping normalization
# Normalize housekeeping
metadata$hk_nf <- hk.factor(ncounts, lod) # Housekeeping normalization
ncounts <- ncounts %*% diag(metadata$hk_nf) # Housekeeping normalization
colnames(ncounts) <- colnames(counts)

# Plot normalized housekeeping
housekeeping.post.normalization.boxplot <- boxplot.count.matrix(ncounts, "Post-Normalization with Housekeeping controls")
housekeeping.post.normalization.boxplot + geom_hline(yintercept = lod,colour="red") + geom_smooth(se=T, aes(group=1))

# Drop genes below threshold
all.names <- rownames(ncounts)
ncounts <- ncounts[(rowSums(ncounts < lod) < round((LOD.THRESHOLD * ncol(ncounts)),0)),]
filtered.names <- rownames(ncounts)
filter.out.genes <- all.names[all.names %!in% filtered.names]

# Save metadata dataframe and the normalised counts matrix
RDS.METADATA.OUTPUT <- paste("BMI", OUTPUT.NAME, "metadata.rds", sep = "_")
RDS.NORMALIZED.OUTPUT <- paste("BMI", OUTPUT.NAME, "HK-normalized-counts.rds", sep = "_")
saveRDS(metadata, RDS.METADATA.OUTPUT) 
saveRDS(ncounts, RDS.NORMALIZED.OUTPUT)

# PCA
#PCA varianza
pcdata <- scale(ncounts, center = TRUE, scale = TRUE)
pc <- pca.loadings(pcdata, 50)
comps <- data.frame(pc$x)
comps$Name <- rownames(comps)
comps <- comps %>% left_join(metadata, by = c(Name = "Sample.name"))

#- Pass labels to legend and/or title
# Select principal components to plot
PRINCIPAL.COMPONENT.X <- 1
PRINCIPAL.COMPONENT.Y <- 2

PCA.PLOT.TITLE <- paste("PCA Plot | Component", PRINCIPAL.COMPONENT.X, "-", PRINCIPAL.COMPONENT.Y, sep = " ")

plotPCA(comps, PRINCIPAL.COMPONENT.X, PRINCIPAL.COMPONENT.Y, key.label
        , legend.label =  LABEL.OF.INTEREST, plot.title = PCA.PLOT.TITLE)

# Histograma PCA sd???


#counts <- counts.old
#metadata <- metadata.old
#eset <- eset.old
### ------------------------
### Differential expression analysis
### ------------------------
design <- as.formula(paste("~", key.label, sep = ""))
rownames(metadata) <- metadata$Sample.name

dea.o <- getDE(ncounts = round(ncounts), metadata = metadata, design = design)
dea.p <- dea.o %>% dplyr::filter(padj < P-VALUE.THRESHOLD) %>% arrange(-log2FoldChange)

# Save DEA Results
P.VALUE.STRING <- paste("PVALUE", P-VALUE.THRESHOLD, sep = "-")
DEA.OUTPUT.NAME <- paste(OUTPUT.NAME, "DEA.csv", sep = "_")
DEA.SIGNIFICATIVE.OUTPUT.NAME <- paste(OUTPUT.NAME, P.VALUE.STRING, "DEA.csv", sep = "_")

write_csv(dea.o, DEA.OUTPUT.NAME)
write_csv(dea.p, DEA.SIGNIFICATIVE.OUTPUT.NAME)

# Show DEA Results
#dplyr::select(dea.p, type, gene.name, id, baseMean, log2FoldChange, pvalue, padj)

# Volcano Plot
tab = data.frame(logFC = dea.o$log2FoldChange, negLogPval = -log10(dea.o$padj))
tab2 = data.frame(logFC = dea.o$log2FoldChange, negLogPval = -log10(dea.o$padj), Gene=dea.o$gene.name)

plot.Volcano(tab, tab2, LOG.FC.THRESHOLD, P-VALUE.THRESHOLD, plot.title = "Volcano Plot")

### ------------------------
### Pathway analysis
### ------------------------

res<- getDE.raw(ncounts = round(ncounts), metadata = metadata, design = design) # the core script for the pathway analysis expects a dataframe with the DEA named res

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
# (THE COMPARISON ARGUMENT DOES NOTHING???)
# SHOULD GO ANNOTATIONS BE REMOVED?? THEY TAKE THE LONGEST, OR ALLOW TO CHOOSE???
ENRICHMENT.ALL.OUTPUT <- paste(OUTPUT.NAME, "Enrichment-ALL.csv", sep = "_")
ENRICHMENT.OVER.OUTPUT <- paste(OUTPUT.NAME, "Enrichment-OVER.csv", sep = "_")
ENRICHMENT.UNDER.OUTPUT <- paste(OUTPUT.NAME, "Enrichment-UNDER.csv", sep = "_")

# Enrich - All
enrich.rs <- enrich.cp(res, LABEL.OF.INTEREST, type="all", pval.threshold = P-VALUE.THRESHOLD, lfc.threshold = LOG.FC.THRESHOLD)
enrich.rs.summary <- enrich.rs$summary %>% arrange(p.adjust)
enrich.rs.summary <- convert.enriched.ids(enrich.rs.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrich.rs.summary, path = ENRICHMENT.ALL.OUTPUT)

### -- # Faster to get over and under from the all??? # -- ###

# Enrich - Over
enrich.rs.over <- enrich.cp(res, LABEL.OF.INTEREST, type="over", pval.threshold = P-VALUE.THRESHOLD, lfc.threshold = LOG.FC.THRESHOLD)
enrich.rs.over.summary <- enrich.rs.over$summary %>% arrange(p.adjust)
enrich.rs.over.summary <- convert.enriched.ids(enrich.rs.over.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrich.rs.over.summary, path = ENRICHMENT.OVER.OUTPUT)

# Enrich - Under
enrich.rs.under <- enrich.cp(res, LABEL.OF.INTEREST, type="under", pval.threshold = P-VALUE.THRESHOLD, lfc.threshold = LOG.FC.THRESHOLD)
enrich.rs.under.summary <- enrich.rs.under$summary %>% arrange(p.adjust)
enrich.rs.under.summary <- convert.enriched.ids(enrich.rs.under.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrich.rs.under.summary, path = ENRICHMENT.UNDER.OUTPUT)

# Show an enrichment summary
enrich.rs.summary
enrich.rs.over.summary
enrich.rs.under.summary

### Dotplots
DOTPLOT.ENRICHMENT.ALL.TITLE <- "Dotplot of All Enriched Pathways after DEA"
DOTPLOT.ENRICHMENT.OVER.TITLE <- "Dotplot of Overexpressed Enriched, Pathways after DEA"
DOTPLOT.ENRICHMENT.UNDER.TITLE <- "Dotplot of Underexpressed Enriched Pathways after DEA"
# KEGG enrichment in all significantly expressed genes:
dotplot(enrich.rs$kg, x="count", showCategory=10, colorBy="qvalue", title = DOTPLOT.ENRICHMENT.ALL.TITLE)
# KEGG enrichment in over-expressed genes:
dotplot(enrich.rs.over$kg, x="count", showCategory=10, colorBy="qvalue", title = DOTPLOT.ENRICHMENT.OVER.TITLE)
# KEGG enrichment in under-expressed genes
dotplot(enrich.rs.under$kg, x="count", showCategory=10, colorBy="qvalue", title = DOTPLOT.ENRICHMENT.UNDER.TITLE)

# <!-- GO BP Enrichment in all the significantly expressed genes: -->
#dotplot(enrich.rs$bp, x="count", showCategory=15, colorBy="qvalue")
# <!-- GO BP Enrichment in over-expressed genes: -->
#dotplot(enrich.rs.over$bp, x="count", showCategory=10, colorBy="qvalue")
# <!-- GO BP Enrichment in under-expressed genes: -->
#dotplot(enrich.rs.under$bp, x="count", showCategory=10, colorBy="qvalue")

# View KEGG Path
enrich.kegg.all.summary.filtered <- ontology.enrichment.qval.filter(enrichment.summary = enrich.rs.summary,
                                                                    ontology = "kg",
                                                                    qval = Q.VALUE.THRESHOLD)
#
plot(x = enrich.kegg.all.summary.filtered$qvalue)
enrich.kegg.all.summary.filtered

#viewkeggpath(path="hsa04151", enrichment = enrich.rs.summary, dea.p = dea.p)
for(pathid in enrich.kegg.all.summary.filtered$ID){
  viewkeggpath(path=pathid, enrichment = enrich.rs.summary, dea.p = dea.p)
}

