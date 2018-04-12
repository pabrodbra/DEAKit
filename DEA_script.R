### ------------------------
### DEA_script.R
### Author: Pablo Rodríguez Brazzarola
### ------------------------

rm(list = ls())
dev.off()

### Load DEA_functions.R
source("DEA_functions.R")

### Argument check (disabled by now)
# CONSTANTS
number.of.keys <- 2
LABEL.KEY1 <- "pCR"
LABEL.KEY2 <- "NAC"

KEY.OF.INTEREST <- 1
LABEL.OF.INTEREST <- "PCR"


SEED <- 12345
FOV.THRESHOLD <- 80
LOD.THRESHOLD <- (100-FOV.THRESHOLD)/100
BD.THRESHOLDS <- c(0.05, 2.25)
ADJUSTED.P.VALUE <- 0.05
LOG.FC.THRESHOLD <- 1

# PATHS
data.directory <- "data" 
rcc.directory <- file.path(data.directory, "PanCancer Pathways - 57")
rlf.filename <- "NS_CancerPath_C2535.rlf"

# SETUP
set.seed(SEED)
keys.vector <- unlist(lapply(seq_len(number.of.keys), function(x) paste("key",x,sep = "")), use.names = FALSE)

### ------------------------
### Load RCC files. Build counts matrix
### ------------------------
# Seleccion de archivos en base a la llave de interés, y los valores

# Load data + Extract metadata
metadata <- extract.rcc.metadata(rcc.directory, keys.vector)

# Retrieve Set and Count matrix
rcc.set.and.count.matrix <- get.RCC.set.and.counts(metadata, rcc.directory, rlf.filename)
eset <- rcc.set.and.count.matrix$set
counts <- rcc.set.and.count.matrix$count.matrix

### ------------------------
### Quality Control
### ------------------------

# Field of View (FOV) Plots
plotFOV(eset = eset, metadata = metadata, fov_threshold = FOV.THRESHOLD,
        comparison.key = keys.vector[1], legend.label = LABEL.KEY1)

plotFOV(eset = eset, metadata = metadata, fov_threshold = FOV.THRESHOLD,
        comparison.key = keys.vector[2], legend.label = LABEL.KEY2)

# Binding Density (BD) Plots
plotBD(eset = eset, metadata = metadata,
       comparison.key = keys.vector[1], legend.label = LABEL.KEY1)

plotBD(eset = eset, metadata = metadata,
       comparison.key = keys.vector[2], legend.label = LABEL.KEY2)

#/* vv FIX PLOT LABELS vv
# Positive Controls
boxplot.expr(eset,is.positive)

# Negative Controls
boxplot.expr(eset,is.negative)

# Noise Threshold
lodcounts  <-  extract.pred(eset, is.negative)
lod <-  mean(lodcounts$count) + 2 * sd(lodcounts$count)

# Housekeeping Genes
housekeeping.boxplot<-boxplot.expr(eset,is.housekeeping)
housekeeping.boxplot+geom_hline(yintercept = (lod),colour="red")

# ^^ FIX PLOT LABELS ^^ */
############ PORQUE HAY 2 TIPOS DE HOUSEKEEPING???? (V2?)

# Expression of all the housekeeping genes in each sample
condition.boxplot <- boxplot.condition(counts, condition = "Housekeeping",
                                       title.label = "Housekeeping Genes Expression")
condition.boxplot
multiple.condition.boxplot.1 <- boxplot.multiple.conditions(counts, conditions = c("Housekeeping"), 
                                                            title.label = "Housekeeping Genes Expression")
multiple.condition.boxplot.1

# Expression of multiple condition genes in each sample
multiple.condition.boxplot <- boxplot.multiple.conditions(counts, conditions = c("Endogenous", "Housekeeping"), 
                                                          title.label = "Housekeeping + Endogenous Genes Expression")
multiple.condition.boxplot+geom_hline(yintercept = (lod),colour="red")


### ------------------------
### Normalization
### ------------------------
counts.old <- counts
metadata.old <- metadata
eset.old <- eset


# Pre-Normalization
metadata$pos_nf = pos.factor(eset)

# REMOVE CONTROL?? - POSITIVE AND NEGATIVE??
counts <- counts[!grepl("Positive", rownames(counts)),]
counts <- counts[!grepl("Negative", rownames(counts)),]

# Plot pre-normalization
pre.normalization.boxplot <- boxplot.count.matrix(counts, "pre-normalization. Positive controls")
pre.normalization.boxplot + geom_hline(yintercept = lod,colour="red")

# Normalize
ncounts = counts %*% diag(metadata$pos_nf)
colnames(ncounts) = colnames(counts)

# Plot post-normalization
post.normalization.boxplot <- boxplot.count.matrix(ncounts, "post-normalization. Positive controls")
post.normalization.boxplot + geom_hline(yintercept = lod,colour="red")

# Housekeeping normalization

# Normalize housekeeping
metadata$hk_nf <- hk.factor(ncounts, lod)
ncounts <- ncounts %*% diag(metadata$hk_nf)
colnames(ncounts) <- colnames(counts)

# Plot normalized housekeeping
housekeeping.post.normalization.boxplot <- boxplot.count.matrix(ncounts, "post-normalization. Positive controls")
housekeeping.post.normalization.boxplot + geom_hline(yintercept = lod,colour="red") + geom_smooth(se=T, aes(group=1))

# Drop genes

all.names <- rownames(ncounts)
ncounts <- ncounts[(rowSums(ncounts < lod) < round((LOD.THRESHOLD * ncol(ncounts)),0)),]
filtered.names <- rownames(ncounts)
filter.out.genes <- all.names[all.names %!in% filtered.names]

# Save metadata dataframe and the normalised counts matrix
saveRDS(metadata, "BMI-DE_metadata.rds")
saveRDS(ncounts, "BMI-DE_HKnormalised_counts.rds")

# PCA
#PCA varianza
pcdata = scale(ncounts, center = TRUE, scale = TRUE)
pc = pca.loadings(pcdata, 50)
comps = data.frame(pc$x)
comps$Name = rownames(comps)
comps = comps %>% left_join(metadata, by = c(Name = "Sample.name"))

#- Pass labels to legend and/or title
plotPCA(comps, 1, 2, keys.vector[1])
plotPCA(comps, 1, 3, keys.vector[1])

plotPCA(comps, 1, 2, keys.vector[2])
plotPCA(comps, 1, 4, keys.vector[2])

#counts <- counts.old
#metadata <- metadata.old
#eset <- eset.old
### ------------------------
### Differential expression analysis
### ------------------------

# PARAMETRIZABLE ESTA?? -- en funcion dfe lo que se elige al principio
design = ~ key1

#PADJ -- PARAMETRIZABLE TMB
dea.o <- getDE(ncounts = round(ncounts), metadata = metadata, design = design)
dea.p <- dea.o %>% dplyr::filter(padj <ADJUSTED.P.VALUE) %>% arrange(-log2FoldChange)

# Save DEA Results
write_csv(dea.o, "NAC-DEA.csv")
write_csv(dea.p, "NAC-DEA_sign.csv") # PADJ nombre

# Show DEA Results
#dplyr::select(dea.p, type, gene.name, id, baseMean, log2FoldChange, pvalue, padj)

# Volcano Plot

# Include TABs inside volcano
tab = data.frame(logFC = dea.o$log2FoldChange, negLogPval = -log10(dea.o$padj))
tab2 = data.frame(logFC = dea.o$log2FoldChange, negLogPval = -log10(dea.o$padj), Gene=dea.o$gene.name)
# LOGFC THREHSOLD PARAM

lfc <-  LOG.FC.THRESHOLD
pval <- ADJUSTED.P.VALUE

plot.Volcano(tab, tab2, lfc, pval)

### ------------------------
### Pathway analysis --- TODO
### ------------------------

res<- getDE.raw(ncounts = round(ncounts), metadata = metadata, design = design) # the core script for the pathway analysis expects a dataframe with the DEA named res

# Setup ENTREZ - Takes time
mart = biomaRt::useMart(biomart = "ensembl", dataset = biomart_dataset)
entrez = biomaRt::getBM(attributes = c("refseq_mrna", "entrezgene"), mart = mart)
entrez$entrezgene = as.character(entrez$entrezgene)
entrezsymbol = biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), mart = mart)
entrezsymbol$entrezgene = as.character(entrezsymbol$entrezgene)

converted = unlist(lapply(strsplit(res$rowname, "_", fixed = TRUE), "[", 4))
converted = unlist(lapply(strsplit(converted, ".", fixed = TRUE), "[", 1))
converted = paste0("NM_", converted)
res$refseq_mrna = converted
lost.genes = length(res$rowname) - length((res %>% data.frame() %>% left_join(entrez, by = "refseq_mrna") %>% 
                                             filter(!is.na(entrezgene)) %>% filter(!is.na(log2FoldChange)) %>% filter(!is.na(lfcSE)))$rowname)

### Enrichments (INVASIVENESS??)
# Enrich - All
enrich.rs <- enrich.cp(res, "Invasiveness", type="all")
enrich.rs.summary = enrich.rs$summary %>% arrange(p.adjust)
enrich.rs.summary = convert.enriched.ids(enrich.rs.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrich.rs.summary,path = "NAC_enrichment.csv" )

# Enrich - Over
enrich.rs.over = enrich.cp(res, "Invasiveness", type="over")
enrich.rs.over.summary = enrich.rs.over$summary %>% arrange(p.adjust)
enrich.rs.over.summary = convert.enriched.ids(enrich.rs.over.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrich.rs.over.summary,path = "NACover_enrichment.csv" )

# Enrich - Under
enrich.rs.under = enrich.cp(res, "Invasiveness", type="under")
enrich.rs.under.summary = enrich.rs.under$summary %>% arrange(p.adjust)
enrich.rs.under.summary = convert.enriched.ids(enrich.rs.under.summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrich.rs.under.summary,path = "NACunder_enrichment.csv" )

# Show an enrichment summary
enrich.rs.summary
enrich.rs.over.summary
enrich.rs.under.summary

### Dotplots
# <!-- GO BP Enrichment in all the significantly expressed genes: -->
dotplot(enrich.rs$bp, x="count", showCategory=15, colorBy="qvalue")
# KEGG enrichment in all significantly expressed genes:
dotplot(enrich.rs$kg, x="count", showCategory=10, colorBy="qvalue")

# <!-- GO BP Enrichment in over-expressed genes: -->
dotplot(enrich.rs.over$bp, x="count", showCategory=10, colorBy="qvalue")
# KEGG enrichment in over-expressed genes:
dotplot(enrich.rs.over$kg, x="count", showCategory=10, colorBy="qvalue")

# <!-- GO BP Enrichment in under-expressed genes: -->
dotplot(enrich.rs.under$bp, x="count", showCategory=10, colorBy="qvalue")
# KEGG enrichment in under-expressed genes
dotplot(enrich.rs.under$kg, x="count", showCategory=10, colorBy="qvalue")

# View KEGG Path
viewkeggpath(path="hsa04151", enrichment = enrich.rs.summary, dea.p = dea.p)

