### ------------------------
### Install packages ###
### ------------------------
orgdb = "org.Hs.eg.db"
biomart_dataset = "hsapiens_gene_ensembl"
keggname = "hsa"
cran.repo = "https://cran.rstudio.com"

packages <- c("dplyr", "readr", "matrixStats", "tidyr","ggplot2", "cowplot", "rgl", "calibrate", "rmarkdown", "pathfindR")
diff.packages <- setdiff(packages, rownames(installed.packages()))
if (length(diff.packages) > 0) {
  install.packages(diff.packages, dependencies = TRUE, repos = cran.repo)  
}

source("http://bioconductor.org/biocLite.R")
bioconductor.packages <- c("NanoStringQCPro", "DESeq2", "pathview", "biomaRt", "clusterProfiler", "org.Hs.eg.db")
diff.bioconductor <- setdiff( bioconductor.packages, rownames(installed.packages()))
if (length(diff.bioconductor) > 0) {
  biocLite(diff.bioconductor, suppressUpdates = FALSE, ask = FALSE)  
}

package.check <- lapply(c(packages, bioconductor.packages), require, character.only = TRUE)

### COMMENTS
## CoreDEA2
# Helpers really needed? Puede ser parametrizable como con el dotplot
# TODO :: AUTOMATIZE QUALITY CONTROL (REMOVAL OF SAMPLES THAT DO NOT PASS FOV AND BD FILTER)
# WHAT IS BMI??? -> plotFOV.BMI = function(eset, metadata, fov_threshold)
# Multiple Conditions can be applied to 1 condition too, remove for only 1 condition??
# Boxplot Expr - Why does it add 'spot' variable to dataframe?

## CoreDEA4
# What is GSEA
# Some of the results in an enrichment summary are named GO::XXXXXXX -> Should this change?
# How to select which pathway you want to view?


### ------------------------
### Helpers
### ------------------------

my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

is.positive = function(column) {
  return(grepl("Pos", column))
}
is.negative = function(column) {
  return(grepl("Neg", column))
}
is.spikein = function(column) {
  return(grepl("Spike", column))
}
is.ligation = function(column) {
  return(grepl("Ligati", column))
}
is.housekeeping = function(column) {
  return(grepl("Housekee", column))
}

geo.pos = function(eset) {
  counts = exprs(eset)
  hk = counts[grepl("Positive", rownames(counts)), ]
  geo.means = apply(hk, 2, function(col) exp(mean(log(col[col != 0]))))
  return(geo.means)
}

hk.pos = function(counts, lod) {
  hk <- counts[grepl("Housekeeping", rownames(counts)),]
  # Normalize
  abovenoise = rowSums(hk > (lod)) >= (ncol(hk))
  hk.abovenoise = hk[abovenoise,]
  above.mean = (apply(hk.abovenoise,1,mean))>= 200
  hk.sel= hk.abovenoise[above.mean,]
  hk.norm=rownames(hk.sel)
  
  hk = counts[hk.norm,]
  geo.means = apply(hk, 2, function(col) exp(mean(log(col[col != 0]))))
  return(geo.means)
}

pos.factor = function(eset) {
  geo.means = geo.pos(eset)
  nf = mean(geo.means)/geo.means
  return(nf)
}


hk.factor = function(counts, lod) {
  geo.means = hk.pos(counts, lod)
  nf = mean(geo.means) / geo.means
  return(nf)}

### ------------------------
### coreDEA_1.R | Load RCC files
### ------------------------

### Load data + Extract metadata
remove.rows.with.na.and.condition.in.columns <- function(df, columns, condition = c("0","1")) {
  cond1 <- (df[, columns] == condition[1])
  cond2 <- (df[, columns] == condition[2])
  ret <- df[as.logical(rowSums(cbind(cond1, cond2))), ]
  complete.rows <- complete.cases(ret[][columns])
  return(ret[complete.rows, ])
}


extract.rcc.metadata <- function(rcc.directory, keys.vector, condition = c("0","1")){
  key.label <- tail(keys.vector, n=1)
  metadata <- data.frame(fname=list.files(rcc.directory, pattern="*.RCC")) %>%
    tidyr::separate(fname, c("Sample.name", keys.vector), sep="_", remove=F) %>%
    tidyr::separate(key.label, c(key.label, "ext"), sep="\\.") %>%
    dplyr::select(c("fname","Sample.name", key.label))
  
  metadata <- remove.rows.with.na.and.condition.in.columns(metadata, key.label, condition)
  
  return(metadata)
}

### Get the right RCC files into an RCCSet
get.RCC.set.and.counts <- function(sampletype, rcc.directory, rlf.file){
  rccSet <- newRccSet(
    rccFiles               = file.path(rcc.directory, sampletype$fname)
    #,rccCollectorToolExport = file.path(exampleDataDir, "nSolver", "RCC_collector_tool_export.csv")
    ,rlf                    = file.path(rlf.file)
    #,cdrDesignData          = file.path(exampleDataDir, "CDR", "CDR-DesignData.csv")
    #,extraPdata             = file.path(exampleDataDir, "extraPdata", "SampleType.txt")
    #,blankLabel             = "blank"
    # ,experimentData.name    = ""
    # ,experimentData.lab     = " "
    # ,experimentData.contact = ""
    # ,experimentData.title   = ""
    # ,experimentData.abstract= ""
  )
  colnames(rccSet) <- sampletype$Sample.name
  ret <- list(set = rccSet,
              count.matrix = rccSet@assayData$exprs)
  return(ret) # Return the RCC expression set
}

### ------------------------
### coreDEA_2.R | QC and Normalization
### ------------------------

'%!in%' <- function(x,y)!('%in%'(x,y))

### Field of View plots
plotFOV = function(eset, metadata, fov.threshold, comparison.key, legend.label) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    left_join(metadata, by=c("rowname"="Sample.name"))
  pdat$pcounted = pdat$FovCounted/pdat$FovCount * 100
  
  p <- ggplot(pdat, aes(rowname, pcounted, color=pdat[[comparison.key]])) + 
    geom_point() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$pcounted))) +
    ylab("percentage of FOV counted") + xlab("Sample.name") +
    geom_hline(yintercept=fov.threshold, color="red") +
    labs(color=legend.label) +
    scale_colour_brewer(palette="Set1")
  return(p)
}

### Binding Density plot
plotBD <- function(eset, metadata, comparison.key, legend.label, y.thresholds = c(0.05, 2.25)) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    left_join(metadata, by=c("rowname"="Sample.name"))
  p <- ggplot(pdat, aes(rowname, BindingDensity, color=pdat[[comparison.key]])) + 
    geom_point() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$BindingDensity))) +
    ylab("Binding density") + xlab("Sample.name") +
    geom_hline(yintercept=y.thresholds[1], color="red") +
    geom_hline(yintercept=y.thresholds[2], color="red") +
    labs(color=legend.label) +
    scale_colour_brewer(palette="Set1")
  return(p)
}

### Extract predicates
extract.pred = function(eset, predicate, counts=FALSE) {
  if(!counts) {
    counts = data.frame(exprs(eset))
  } else {
    counts = eset
  }
  toplot = counts[predicate(rownames(counts)),] %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("Sample.name", "count", -rowname)
  colnames(toplot) = c("spot", "Sample.name", "count")
  toplot = toplot %>% left_join(metadata, by="Sample.name")
  return(toplot)
}

### BoxPlot
boxplot.count.matrix <- function(counts, title.label = "Genes Expression"){
  current.count <-  counts %>% data.frame() %>% tidyr::gather("sample", "count")
  current.count$sample<-gsub("X","",current.count$sample)
  ggplot(current.count, aes(sample, count)) + geom_boxplot(colour = "black", fill = "#56B4E9",outlier.size = 0.5) +
    scale_y_continuous(trans = "log2") + 
    xlab("Samples") + ylab("Counts") + ggtitle(title.label) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size = 5))
}

boxplot.multiple.conditions <- function(counts, conditions, 
                                        title.label = ""){
  all.conditions <- do.call(rbind, lapply(conditions, function (x) counts[grepl(x, rownames(counts)),]))
  all.conditions.df <- as.data.frame(all.conditions)
  all.conditions.tidy <-tidyr::gather(all.conditions.df)
  colnames(all.conditions.tidy)<-c("sample","count")
  
  ggplot(all.conditions.tidy, aes(sample, count)) + geom_boxplot(colour = "black", fill = "#56B4E9",outlier.size = 0.5) +
    scale_y_continuous(trans = "log2") + 
    xlab("Samples") + ylab("Counts") + ggtitle(title.label) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size = 5))+
    geom_smooth(se=T, aes(group=1))
}

boxplot.expr<-function(ese,predicate) {
  DF<-extract.pred(eset,predicate)%>%
    tidyr::separate(spot,c("a","b","c","d"),sep="_",remove=F)%>%
    tidyr::unite(gene,b,c,sep="_")%>%
    dplyr::select(gene,count)
  ggplot(DF,
         aes(x=gene,y=count)) + geom_boxplot(colour = "black", fill = "#56B4E9")+
    scale_y_continuous(trans = "log2")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
    ylab("Counts") + xlab("Genes")
}

### PCA

pca.loadings = function(object, ntop = 500) {
  object = as.matrix(object)
  rv <- matrixStats::rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  pca <- prcomp(t(object[select, ]))
  percent.var <- pca$sdev^2/sum(pca$sdev^2)
  names(percent.var) = colnames(pca$x)
  pca$percent.var = percent.var
  
  return(pca)
}

plotPCA = function(comps, pc, nc1, nc2, colorby, legend.label = "Key", plot.title = "") {
  c1str = paste0("PC", nc1)
  c2str = paste0("PC", nc2)
  ggplot(comps, aes_string(c1str, c2str, color=colorby)) + geom_point() + theme_bw() + 
    xlab(paste0(c1str, ": ", round(pc$percent.var[nc1] * 100),  "% variance")) +
    ylab(paste0(c2str, ": ", round(pc$percent.var[nc2] * 100), "% variance")) +
    theme(aspect.ratio=1) +
    ggtitle(plot.title) +
    labs(color=legend.label) +
    scale_colour_brewer(palette="Set1")
}

### ------------------------
### coreDEA_3.R | Differential expression analysis
### ------------------------

getDE <- function(ncounts, metadata, design){
  dds <- DESeqDataSetFromMatrix(countData = ncounts, colData = metadata, design = design)
  sizeFactors(dds) <- rep(1, ncol(ncounts))
  dds <- DESeq(dds)
  res <- results(dds) %>% data.frame() %>%
    tibble::rownames_to_column() %>%
    arrange(pvalue) %>%
    mutate(rowname=gsub("NM_", "NM-", rowname)) %>% 
    tidyr::separate(rowname, c("type", "gene.name", "id"), sep="_") %>% 
    mutate(id=gsub("NM-", "NM_", id))
  
  return(res)
}

getDE.raw <- function(ncounts, metadata, design){ ## Gets the same DE data, less processing, for GSEA in coreDEA4.R
  dds = DESeqDataSetFromMatrix(countData = ncounts, colData = metadata, design = design)
  sizeFactors(dds) = rep(1, ncol(ncounts))
  dds = DESeq(dds)
  res = res = results(dds) %>% data.frame() %>% 
    tibble::rownames_to_column() %>% 
    arrange(pvalue)
  return(res)
}

plot.Volcano <- function(tab, tab2, lfc.threshold, pval.threshold, plot.title = ""){
  par(mar = c(5, 4, 4, 4))
  p <- plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
  title(main = plot.title)
  # signGenes <- (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
  points(tab[(abs(tab$logFC) > lfc.threshold), ], pch = 16, cex = 0.8, col = "orange") 
  points(tab[(tab$negLogPval > -log10(pval.threshold)), ], pch = 16, cex = 0.8, col = "green") 
  points(tab[(abs(tab$logFC) > lfc.threshold & tab$negLogPval > -log10(pval.threshold)), ], pch = 16, cex = 0.8, col = "red") 
  abline(h = -log10(pval.threshold), col = "green3", lty = 2) 
  abline(v = c(-lfc.threshold, lfc.threshold), col = "blue", lty = 2) 
  mtext(paste("pval =", pval.threshold), side = 4, at = -log10(pval.threshold), cex = 0.8, line = 0.5, las = 1) 
  mtext(c(paste("-", lfc.threshold, "fold"), paste("+", lfc.threshold, "fold")), side = 3, at = c(-lfc.threshold, lfc.threshold), cex = 0.8, line = 0.5)
  with(subset(tab2, negLogPval > -log10(pval.threshold) & abs(logFC)>1), textxy(logFC, negLogPval, labs=Gene, cex=.4))
}

plot.Volcano2 <- function(tab, lfc.threshold, pval.threshold, plot.title = ""){
  max.abs.log.fc <- max(abs(tab$logFC))
  max.p.value <- max(abs(tab$negLogPval))
  neg.log.p.value <- -log10(pval.threshold)
  subset.tab <- tab[(abs(tab$logFC) > lfc.threshold & tab$negLogPval > -log10(pval.threshold)), ]
  p <- ggplot(data=tab, aes(x=logFC, y=negLogPval))+
    theme(legend.position = "none") +
    ggtitle(plot.title)+
    geom_point(alpha=0.6, size=1.75, color="black")+
    geom_point(data = tab[(tab$negLogPval > -log10(pval.threshold)), ], alpha=0.6, size=1.75, color="green")+
    geom_point(data = subset.tab,
               alpha=0.6, size=1.75, color="red")+
    xlim(c(-max.abs.log.fc,max.abs.log.fc)) + ylim(c(0, max.p.value+1)) +
    xlab(expression(log[2]~fold~change)) + ylab(expression(-log[10]~pvalue)) +
    geom_vline(xintercept=c(-lfc.threshold,lfc.threshold), linetype="dotted", color = "blue", size=0.8) +
    geom_hline(yintercept=c(-neg.log.p.value,neg.log.p.value), linetype="dotted", color = "green3", size=0.8) +
    annotate("text", label= paste("pval =", pval.threshold), x=max.abs.log.fc-0.25, y=neg.log.p.value+0.5, size=4) +
    annotate("text", label= c(paste("-", lfc.threshold, "fold"), paste("+", lfc.threshold, "fold")), 
             x=c(-lfc.threshold-0.5,lfc.threshold+0.5), y=max.p.value-0.25, size=4) +
    annotate("text", label=subset.tab$Gene,
             x = subset.tab$logFC,
             y = subset.tab$negLogPval+0.6, size = 3)
  return(p)
}

### ------------------------
### coreDEA_4.R | Pathway Analysis
### ------------------------

summarize.cp <- function(res, comparison) {
  summaries = data.frame()
  for (ont in names(res)) {
    ontsum = summary(res[[ont]]) # as.data.frame() instead of summary??
    ontsum$ont = ont
    summaries = rbind(summaries, ontsum)
  }
  summaries$comparison = comparison
  return(summaries)
}

enrich.cp <- function(res, entrez, comparison, type="over", pval.threshold = 0.05, lfc.threshold = 1) {
  res = res %>% data.frame()  %>% left_join(entrez, by = "refseq_mrna") %>% filter(!is.na(entrezgene))
  universe = res$entrezgene
  res <- switch(type,
                all = res %>% filter(padj < pval.threshold),
                over = res %>% filter(log2FoldChange > lfc.threshold & padj < pval.threshold),
                under = res %>% filter(log2FoldChange < -lfc.threshold & padj < pval.threshold))
  genes = res$entrezgene
  #mf = enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
  #cc = enrichGO(genes, OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
  #bp = enrichGO(genes, OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
  kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
  # mf = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", 
  # qvalueCutoff = 1, pvalueCutoff = 1)
  # cc = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", 
  # qvalueCutoff = 1, pvalueCutoff = 1)
  # bp = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", 
  # qvalueCutoff = 1, pvalueCutoff = 1)
  # kg = enrichKEGG(gene = genes, universe = universe, organism = keggname, pvalueCutoff = 1, 
  # qvalueCutoff = 1, pAdjustMethod = "BH")
  all = list(kg = kg)#, mf = mf, cc = cc, bp = bp)
  all[["summary"]] = summarize.cp(all, comparison)
  return(all)
}

convert.enriched.ids <- function(res, entrezsymbol) {
  res = res %>% mutate(geneID = strsplit(as.character(geneID), "/")) %>% 
    tidyr::unnest(geneID) %>% 
    left_join(entrezsymbol, by = c(geneID = "entrezgene")) %>% 
    group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, ont, comparison) %>%
    summarise(geneID = paste(geneID, collapse = "/"), symbol = paste(hgnc_symbol, collapse = "/"))
  return(res)
}

# Enrichment Summary Kegg filter
ontology.enrichment.qval.filter <- function(enrichment.summary, ontology = "kg", qval = 0.05){
  return(enrichment.summary[(enrichment.summary$ont == ontology & enrichment.summary$qvalue < qval),])
}

ontology.enrichment.padjust.filter <- function(enrichment.summary, ontology = "kg", p.adj = 0.05){
  return(enrichment.summary[(enrichment.summary$ont == ontology & enrichment.summary$p.adjust < p.adj),])
}

# View KEGG Path
viewkeggpath <- function(path, enrichment, dea.p, output = "."){
  # Get the genes in path from the enrichment summary
  genesymbols <- strsplit(as.character(enrichment[which(enrichment$ID ==path),12]), split = "/")[[1]]
  geneids <- strsplit(as.character(enrichment[which(enrichment$ID ==path),11]), split = "/")[[1]]
  ll <- data.frame(geneid=geneids, gene.name=genesymbols)
  
  # dataframe with the genes and the log2FC
  kk <- filter(dea.p, gene.name %in% genesymbols) %>% 
    dplyr::select(gene.name, log2FoldChange)
  kk <- kk %>% inner_join(ll,by="gene.name")
  gene.vector <- kk$log2FoldChange
  names(gene.vector) <- kk$geneid
  
  pathview(gene.data=gene.vector, pathway.id=path, species="hsa", limit=list(gene=max(abs(gene.vector)), cpd=1), kegg.dir = output)
}

gene.dea.summary <- function(dea.results, adj.pval.threshold = 0.05){
  ret.df <- data.frame(Gene.symbol = dea.results$rowname,
                       logFC = dea.results$log2FoldChange,
                       adj.P.Val = dea.results$pvalue,
                       stringsAsFactors = FALSE
  ) %>% filter(adj.P.Val < adj.pval.threshold) %>%
    tidyr::separate(Gene.symbol, c("type", "Gene.symbol"), sep="\\_")
  
  ret.df$type <- NULL
  return(ret.df)
}

########################################
# GSEA ? | DEPRECATED??


gsea.cp = function(res, comparison) {
  res = res %>% data.frame() %>% left_join(entrez, by = "refseq_mrna") %>% 
    filter(!is.na(entrezgene)) %>% filter(!is.na(log2FoldChange)) %>% filter(!is.na(lfcSE))
  lfc = data.frame(res)[, "log2FoldChange"]
  lfcse = data.frame(res)[, "lfcSE"]
  genes = lfc/lfcse
  names(genes) = res$entrezgene
  genes = genes[order(genes, decreasing = TRUE)]
  cc = gseGO(genes, ont = "CC", OrgDb = orgdb, nPerm = 500, pvalueCutoff = 1, 
             pAdjustMethod = "BH", verbose = TRUE)
  mf = gseGO(genes, ont = "MF", OrgDb = orgdb, nPerm = 500, pvalueCutoff = 1, 
             pAdjustMethod = "BH", verbose = TRUE)
  bp = gseGO(genes, ont = "bp", OrgDb = orgdb, nPerm = 500, pvalueCutoff = 1, 
             pAdjustMethod = "BH", verbose = TRUE)
  # genes = data.frame(res)[, 'log2FoldChange'] names(genes) = res$entrezgene
  # genes = genes[order(genes, decreasing=TRUE)] genes = genes[!is.na(genes)]
  kg = gseKEGG(geneList = genes, organism = keggname, nPerm = 500, pvalueCutoff = 1, 
               verbose = TRUE)
  if (orgdb == "org.Hs.eg.db") {
    do = gseDO(geneList = genes, nPerm = 500, pvalueCutoff = 1, pAdjustMethod = "BH", 
               verbose = TRUE)
    all = list(mf = mf, cc = cc, bp = bp, kg = kg, do = do)
  } else {
    all = list(mf = mf, cc = cc, bp = bp, kg = kg)
  }
  all[["summary"]] = summarize.cp(all, comparison)
  return(all)
}

convert.core.ids = function(row) {
  core_ids = data.frame(entrezgene = unlist(strsplit(row["core_enrichment"], 
                                                     "/")[[1]])) %>% left_join(entrezsymbol, by = "entrezgene")
  core_symbols = unique(core_ids$hgnc_symbol)
  core_symbols = core_symbols[!is.na(core_symbols)]
  names(core_symbols) = NULL
  return(paste(core_symbols, collapse = "/"))
}

convert.gsea.results = function(res) {
  res$symbols = apply(res, 1, convert.core.ids)
  return(res)
}


