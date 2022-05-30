---
title: "Orebro RNA-seq analysis"
author: "Rodrigo Arcoverde Cerveira & Gustav Joas"
date: '`r format(Sys.Date(), "%Y-%m-%d")`'
output: html_document
abstract: Örebro study RNA-seq data analysis.
knit: (function(inputFile, encoding) {
          rmarkdown::render(inputFile,
                            encoding = encoding, 
                            output_file = paste0(
                              xfun::sans_ext(inputFile), '_', Sys.Date(), '.html'),
                                output_dir = "../results/lab_book/")})
editor_options: 
  chunk_output_type: inline
---

```{r global.options, include=FALSE}
#setting up figures and chunks messages

knitr::opts_knit$set(echo = TRUE,
                     root.dir = getwd(),
                     fig.width = 6, fig.height = 5,
                     fig.align = "center",
                     out.width = 768,
                     fig.pos = "H",
                     warning = FALSE, 
                     message = FALSE)
knitr::opts_chunk$set(warning = FALSE,
                      message = FALSE,
                      fig.width = 6, fig.height = 5,
                     fig.align = "center",
                     out.width = 768,
                     fig.pos = "H")

result.dir <- paste0("results/",Sys.Date(),"/")
figures.dir <- paste0("results/",Sys.Date(),"/", "figures/")
## creates result.dir with date in if not existent
ifelse(isFALSE(dir.exists(paste0("../",result.dir))), dir.create(paste0("../",result.dir),recursive = TRUE),"Result directory for today exists already!")
ifelse(isFALSE(dir.exists(paste0("../",figures.dir))), dir.create(paste0("../",figures.dir),recursive = TRUE),"Result directory for today exists already!")

options(stringsAsFactors = FALSE) 
```

## Loading libraries

```{r, message=FALSE}
#if you do not have libraries, they are located in either CRAN or Bioconductor
library(data.table)
library(kableExtra)
library(forcats)
library(caret)
library(C50)
library(mgsub)
library(DESeq2) 
library(pheatmap) 
library(ggplot2) 
library(ggrepel) 
library(RColorBrewer) 
library(limma) 
library(edgeR) 
library(enrichR) 
library(gridExtra) 
library(stringr)
library(ggVennDiagram)
library(dplyr) 
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(factoextra)
library(rafalib)
library(plotly)
library(tibble)
library(devtools)
library(clusterProfiler)
library(enrichplot)
library(pathview)

```


## Download and load data

```{r}
# Reading data
counts_raw <- data.table::fread("../data/proces-data_core-facility/subreadCounts_hg38ens_minus_frag.txt", header = TRUE, sep = "\t") %>%  as.data.frame()

sampleTable <- read.table("../data/metadata/sample_metadata.csv", sep = ",", header = TRUE) %>%
  mutate(dose = as.factor(plyr::mapvalues(visit, from = c(1,2,4,5), to = c("0_dose1","24_dose1","0_dose2","24_dose2"))),
         group = as.factor(group))

# Renaming genes based on emsembl annotation
edb <- EnsDb.Hsapiens.v86
edb_genes <- genes(edb) %>% as.data.frame()

counts_raw <- counts_raw %>%
  dplyr::select(-c(Chr, Start, End, Strand, Length)) %>%
  dplyr::mutate(Geneid = plyr::mapvalues(x = counts_raw$Geneid, from = edb_genes$gene_id, to = edb_genes$symbol, warn_missing = FALSE))

# Aggregate and sum up genes with same gene symbol, which were basically non-coding RNAs
counts_raw <- aggregate(counts_raw[,-1], list(Geneid=counts_raw[,1]), FUN = sum)
rownames(counts_raw) <- counts_raw$Geneid
counts_raw <- dplyr::select(counts_raw, -c(Geneid))

```

## Basic Quality Control

### Inspect raw data

Matching metadata to countTable and plotting the raw counts for first visual inspection.

```{r}
# Synchronize count data with sample table
counts_raw <- counts_raw[, pmatch(sampleTable$sample_ID, colnames(counts_raw))]
colnames(counts_raw) <- sampleTable$sample_ID
all(rownames(sampleTable$sample_ID) == colnames(counts_raw))

# Visualize distribution of raw counts w/ boxplot and density plot
{
  pdf(paste0("../", figures.dir,"raw_counts_QC.pdf"), width = 10,  height = 8, compress = TRUE )
  rafalib::mypar(1,2,mar=c(10,3,3,2))
  boxplot(log2(as.matrix(counts_raw)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Raw data")
  hist(log2(as.matrix(counts_raw)+1),ylab="",las=2,main="Raw data")
  par(mfrow=c(1,1))
  dev.off()
}
```

### Filter data

Plot detection of genes across samples. All samples are more or less close to average so we do not need to discard any samples.

```{r}
{
  pdf(paste0("../", figures.dir,"number_detected_genes.pdf"), width = 10,  height = 8, compress = TRUE )
  par(mar=c(10,4,3,4))
  barplot(colSums(counts_raw>3),ylab="Number of detected genes",las=2)
  abline(h=median(colSums(counts_raw>3)))
  dev.off()
}
```

Removing reads with the log2 of the counts per million (cpm) lower than 1.

```{r}

#Filter low counts
keep_genes <- rowSums( counts_raw > 5 ) >= 8
counts_filtered<- counts_raw[keep_genes,]
sum(keep_genes)


meanLog2CPM <- rowMeans(log2(cpm(counts_raw) + 1)) 
counts_filtered <- counts_raw[meanLog2CPM > 1, ] 


```

Plot detection rate across genes for raw and filtered counts

```{r}
{
  pdf(paste0("../", figures.dir,"detection_rate_raw_filtered.pdf"), width = 10,  height = 8, compress = TRUE )
  par(mar=c(10,4,3,4))
  hist(rowSums(counts_raw>3))
  hist(rowSums(counts_filtered>3))
  par(mfrow=c(1,1))
  dev.off()
}
```

Plot distribution of the filtered counts

```{r}
{
  pdf(paste0("../", figures.dir,"filtered_counts_distribution.pdf"), width = 10,  height = 8, compress = TRUE )
  rafalib::mypar(1,2,mar=c(10,3,3,2))
  boxplot(log2(as.matrix(counts_filtered)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Filtered data")
  hist(log2(as.matrix(counts_filtered)+1),ylab="",las=2,main="Filtered data")
  par(mfrow=c(1,1))
  dev.off()
}
```

### DESeq object creation, Normalization and data Quality Control

Generating the DESeq dataset.

```{r}
# prepare for DESeq
sampleTable$conditions <- (str_c(sampleTable$group, "_", sampleTable$dose ))
sampleTable$conditions <- factor(sampleTable$conditions)
sampleTable$sex <- factor(sampleTable$sex)
sampleTable$dose <- factor(sampleTable$dose)
sampleTable$group <- factor(sampleTable$group)

sampleTable$time <- sampleTable$conditions  
sampleTable$time <- gsub('.{6}$', '', sampleTable$time)
sampleTable$time <- factor(sampleTable$time)
# create 3 DESeq object
dds1 <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = as.data.frame(sampleTable),
                              design = ~  sex + age + conditions)

dds2 <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = as.data.frame(sampleTable),
                              design = ~  sex + age + time)

dds3 <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = as.data.frame(sampleTable),
                              design = ~  sex + age + dose)
# Normalize with rlog
counts_rlog_normalized <- rlog(dds1, blind = TRUE)

# Normalize with variance stabilizing transformation for later PCA and heatmap
counts_vst_normalized <- vst(dds1, blind = TRUE)
```

Plot distribution of data after normalization

```{r}
# haven't decided which plots to use
vst_matrix <- assay(counts_vst_normalized) 
rlog_matrix <- assay(counts_rlog_normalized)
hist(vst_matrix)
hist(rlog_matrix)
boxplot(vst_matrix ,ylab=expression('Log'[2]~'Read counts'),las=2,main="VST")
boxplot(rlog_matrix ,ylab=expression('Log'[2]~'Read counts'),las=2,main="rlog")

```

### Heatmap

```{r}
# Sample heatmap
sampleDist <- cor(vst_matrix, method = "spearman") 

Metadata <- data.frame(sampleTable$group, sampleTable$dose)
names(Metadata) <- c("Group","Dose")
rownames(Metadata) <- sampleTable$sample_ID

# Plot heatmap
colors<-colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(255)

{
  pdf(paste0("../", figures.dir,"heatmap.pdf"), width = 10,  height = 8, compress = TRUE )
 Heatmap <-  pheatmap(sampleDist, 
           color = colors,
           clustering_distance_rows = as.dist(1 - sampleDist),
           clustering_distance_cols = as.dist(1 - sampleDist), 
           show_rownames = F,
           show_colnames = F,
           clustering_method = "ward.D2",
           annotation_col = Metadata)
  par(mfrow=c(1,1))
  dev.off()
}

print(Heatmap)
```

### Principal Component Analyais (PCA)

Dimensionality reduction for evaluating outliers and global sample clusters for both quality control but also for primary exploratory analysis.

```{r}
#Sample PCA
pcaRes <- prcomp(t(assay(counts_vst_normalized)))
varExp <- round(summary(pcaRes)[[1]],digits = 2)

pcaDF <- data.frame(
  PC1 = pcaRes$x[, 1],
  PC2 = pcaRes$x[, 2],
  Group = sampleTable$group,
    Sample = sampleTable$sample_ID,
  Dose = sampleTable$dose)

{ 
  pdf(paste0("../", figures.dir,"PCA.pdf"), width = 10,  height = 8, compress = TRUE )
  pcaPlot <- ggplot( pcaDF, mapping = aes(x = PC1, y = PC2, color = Group, label = Sample)) + 
  geom_point(aes(shape = Dose), size = 3) +
  labs(x = paste0("PC1 (", varExp[1], " %)"), y = paste0("PC2 (", varExp[2], " %)")) + 
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 10)) +
  scale_color_manual(values = brewer.pal(3, "Accent")) +
  cowplot::theme_cowplot()
  par(mfrow=c(1,1))
  dev.off()
}

print(pcaPlot)

## make a scree plot to compute PC variance
pca.var <- pcaRes$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
{
  pdf(paste0("../", figures.dir,"Scree_plot.pdf"), width = 10,  height = 8, compress = TRUE )  
  barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
  par(mfrow=c(1,1))
  dev.off()
}
## exploring which variables contribute the most for each PC
var <- get_pca_var(pcaRes)
var$contrib %>%
  as.data.frame() %>%
  arrange(desc(Dim.1)) %>%
  head()
  
```

## Differential Gene Expression Analysis

```{r}
# Run the DESeq2 analysis
dds1 <-  DESeq(dds1)
dds2 <-  DESeq(dds2)
dds3 <-  DESeq(dds3)

# contrast
resultsNames(dds1)
resultsNames(dds2)
resultsNames(dds3)

# Results dds1

# Baselines DEGs
res_baseline_pos_neg_dose1 <- as.data.frame(results(dds1, contrast = c("conditions", "Pos_0_dose1", "Neg_0_dose1")))
res_baseline_pos_neg_dose2 <- as.data.frame(results(dds1, contrast = c("conditions", "Pos_0_dose2", "Neg_0_dose2")))
res_baseline_pos_doses <- as.data.frame(results(dds1, contrast = c("conditions", "Pos_0_dose2", "Pos_0_dose1")))
res_baseline_neg_doses <- as.data.frame(results(dds1, contrast = c("conditions", "Neg_0_dose2", "Neg_0_dose1")))
#add vaseline for Pos vs neg regardless of dose (theres is a difference here), and then regardles of statu (not here)


# 0-24h Degs
res_neg_dose1 <- as.data.frame(results(dds1, contrast = c("conditions", "Neg_24_dose1", "Neg_0_dose1")))
res_pos_dose1 <- as.data.frame(results(dds1, contrast = c("conditions", "Pos_24_dose1", "Pos_0_dose1")))
res_neg_dose2 <- as.data.frame(results(dds1, contrast = c("conditions", "Neg_24_dose2", "Neg_0_dose2")))
res_pos_dose2 <- as.data.frame(results(dds1, contrast = c("conditions", "Pos_24_dose2", "Pos_0_dose2")))

# Results dds2
res_pos <- as.data.frame(results(dds2, contrast = c("time", "Pos_24", "Pos_0")))
res_neg <- as.data.frame(results(dds2, contrast = c("time", "Neg_24", "Neg_0")))

# Results dds3
res_dose1 <- as.data.frame(results(dds3, contrast = c("dose", "24_dose1", "0_dose1")))
res_dose2 <- as.data.frame(results(dds3, contrast = c("dose", "24_dose2", "0_dose2")))
```


### Volcano plots

Volcano plot for data exploratory analysis, filtering for significant values with False Discovery Rate \< 0.05 and log fold change greater than 1.

```{r}
# check only the negative fold change to see if there's a different pahtway of negative fold changes

results_list <- list(res_baseline_pos_neg_dose1, res_baseline_pos_neg_dose2, res_baseline_pos_doses, res_baseline_neg_doses, res_neg_dose1, res_pos_dose1, res_neg_dose2, res_pos_dose2, res_pos, res_neg, res_dose1, res_dose2)

names(results_list) <- c("Baseline Dose 1", "Baseline Dose 2", "Baseline Conv", "Baseline Naïve", "Naïve Dose 1", "Conv Dose 1", "Naïve Dose 2", "Conv Dose 2", "Conv","Naïve", "Dose 1", "Dose 2") 
    results_list_plot <- lapply(results_list, function(x) {
  x <- x %>% mutate(test_padj = case_when(pvalue < 0.05 & log2FoldChange >= 1 ~ "Up regulated",
                                     pvalue < 0.05 & log2FoldChange <= -1 ~ "Down regulated",
                                     TRUE ~ "Not significant"),
                    test_pvalue = case_when(pvalue < 0.05 & log2FoldChange >= 1 ~ "Up regulated",
                                     pvalue < 0.05 & log2FoldChange <= -1 ~ "Down regulated",
                                     TRUE ~ "Not significant"),
                    )
  })
# color by group dose1 and dose2, or neg and pos. 2 plots colored based on condition
# color by group dose1 and dose2, or neg and pos. 2 plots colored based on condition
# function to plot and save volcano plot
.volcano_plot <-  function(x){
  subsetted <- subset(results_list_plot[[x]] %>% tibble::rownames_to_column("gene_symbol"), abs(log2FoldChange) >= 2 & pvalue < 0.05)
  
  # count up regulated DEGs
  deg_number_up <- subset(results_list_plot[[x]] %>% tibble::rownames_to_column("gene_symbol"), log2FoldChange >= 1 & pvalue < 0.05) %>%
    nrow()
  # count down regulated DEGs
  deg_number_down <- subset(results_list_plot[[x]] %>% tibble::rownames_to_column("gene_symbol"), log2FoldChange <= -1 & pvalue < 0.05) %>%
    nrow()
  
  results_list_plot[[x]] %>% 
  ggplot(aes(x=log2FoldChange,y=-log10(pvalue))) +
  geom_point(aes(colour=test_pvalue), size=3, alpha=0.3) +
  scale_color_manual(values = c("Down regulated"="blue", "Not significant" ="grey80", "Up regulated"="red")) +
  xlab(expression("Fold Change (Log"[2]*")")) +
  ylab(expression("-Log"[10]*"(p-value)")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", size = 1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", size = 1) +
  geom_label_repel(data = subsetted, aes(log2FoldChange,-log10(pvalue),label= gene_symbol), max.overlaps = 10) +
  xlim(-2,6) +
  ylim(0,10) +
 # ggtitle(x) +  
  annotate(
        geom='text',
        colour = "red",
        size = 10,
        x= 4, y=10,
        hjust=0,
        label=paste0(deg_number_up)
    ) +
    annotate(
        geom='text',
        colour = "blue",
        size = 10,
        x= -2, y= 10,
        hjust=0,
        label=paste0(deg_number_down)
    ) +
  cowplot::theme_cowplot(font_size = 20, line_size = 1.5, rel_small = 1.6, rel_tiny = 1,font_family = "Arial",rel_large = 4) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 35))
  }

volcano_plots <- lapply(names(results_list_plot), .volcano_plot)
names(volcano_plots) <- names(results_list_plot)
volcano_plots

lapply(names(volcano_plots), function(x) ggsave(filename = paste0("../", figures.dir, x,"_volcano_plot.png"), width = 7, height = 7,
                                         plot = volcano_plots[[x]]))




```


```{r}

#MA plot
  
 #Counts plot
# plotCounts(d,gene=rownames(res)[1],intgroup="group",normalized=F)
# plotCounts(d,gene=rownames(res)[1],intgroup="group",normalized=T)
```

### Enrichment analysis

Selecting two different databases for search for enrichment analysis. The databases selected were "GO_Biological_Process_2021" and "KEGG_2021_Human". Testing also the Blood Transcriptome modules presented in [Li et al. 2014](10.1038/ni.2789) using some GSEA package.

```{r}
#Enrichment analysis
listEnrichrDbs()

enrichment_list_up <- lapply(results_list_plot, function(x) {
  enrichr(
  genes = rownames(x[x$test_padj == "Up regulated", ]),
  databases = c("GO_Biological_Process_2021", "KEGG_2021_Human","GO_Molecular_Function_2021"))})

enrichment_list_down <- lapply(results_list_plot, function(x) {
  enrichr(
  genes = rownames(x[x$test_padj == "Down regulated", ]),
  databases = c("GO_Biological_Process_2021", "KEGG_2021_Human","GO_Molecular_Function_2021"))})

gmt <- read.gmt("../data/blood_transcription_module/BTM_for_GSEA_20131008.gmt")

```

#### Visualize enriched pathways

```{r echo=TRUE, fig.width=15, fig.height=8}

#Removing items without significant results from list
enrichment_list_up[c("neg_dose1")] <- NULL
enrichment_list_down[c("baseline", "neg_dose1", "pos_dose1", "neg", "dose1")] <- NULL

#Visualize significant terms
enrichment_list <- enrichment_list_up

plots_enrich <- list()
for(pathway in 1:length(enrichment_list_up[[1]])){
  
  rbinded_enrich_list <- lapply(enrichment_list, function(x)x[[pathway]])
  print(names(enrichment_list[[1]])[[pathway]])
  rbinded_enrich_df <- data.table::rbindlist(l = rbinded_enrich_list, idcol = TRUE, fill = TRUE)
  
  plots_enrich[[pathway]] <- rbinded_enrich_df %>%
  group_by(.id) %>%
  arrange(Adjusted.P.value) %>%
  slice(1:10) 
  
}

lapply(plots_enrich, function(x){
  x %>%
  ggplot(aes(x = Term, y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", width = 0.05) + geom_point(size = 3) +
    theme_minimal() +
    theme(text = element_text(size = 10),
          plot.title = element_text(hjust = (5 / nchar(results_list)) * 2), 
          plot.margin = margin(t = 5, r = 50, b = 5, unit = "pt"), 
          axis.text.y = element_text(size = 8)) +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") + 
    theme_bw() +
    labs(y = expression("Adjusted P value, Log"[10]*"")) +
  facet_wrap(~.id)
})


```

### Gene Set Enrichment Analysis with ClusterProfiler

#### Annotaions

```{r, message=F, warning=F}
# SET THE DESIRED ORGANISM HERE
organism <- "org.Hs.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
```

#### Prepare Input

```{r}
process_to_gsea <- function(x){
  original_gene_list <- x$log2FoldChange
  names(original_gene_list) <- rownames(x)
  gene_list<-na.omit(original_gene_list)
  gene_list <- sort(gene_list, decreasing = TRUE)
}

ls_processed <- lapply(results_list, process_to_gsea)

```

#### Gene Set Enrichment

```{r}
list_gse <- lapply(ls_processed, function(x) GSEA(geneList=x, 
            TERM2GENE = gmt,
             minGSSize = 3, 
             maxGSSize = 800, 
            nPermSimple = 10000,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             pAdjustMethod = "fdr"))


```


#### Ridgeplot
Helpful to interpret up/down-regulated pathways.
```{r fig.width=18, fig.height=12}
# change coloring for up and downregulated
# try also to plot as a heatmap for all the group comparisons

rigdgeplots_list <- lapply(list_gse, ridgeplot)
lapply(names(rigdgeplots_list), function(x) ggsave(filename = paste0("../", figures.dir, x,"_BTM_plot.pdf"), 
                                         plot = rigdgeplots_list[[x]], width = 15, height = 10))
rigdgeplots_list

```

#### Dotplot
```{r echo=TRUE, fig.width=15, fig.height=8}
#require(DOSE)
#dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
```

#### GSEA Plot  
Traditional method for visualizing GSEA result.  
  
```{r fig.height=6, fig.width=8}
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
#gseaplot(gse, by = "all", title = gse$Description[55], geneSetID = 55)
#head(gse, 100)

```


### KEGG Gene Set Enrichment Analysis
#### Prepare Input
```{r}
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(ls_processed[[5]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = res_pos[rownames(res_pos) %in% dedup_ids$SYMBOL,]
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID
# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange
# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y
# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
```


#### Create gseKEGG object
 
```{r}
#kegg_organism = "hsa"
#kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "fdr",
               keyType       = "ncbi-geneid")

```
#### Ridgeplot
Helpful to interpret up/down-regulated pathways.
```{r fig.width=18, fig.height=12}
#gseKEGG(kk2) + labs(x = "enrichment distribution")
```

### Pathview
This will create a PNG and *different* PDF of the enriched KEGG pathway.  
```{r, message=F, warning=F, echo = TRUE}
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04620", species = kegg_organism)
# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="map04620", species = kegg_organism, kegg.native = F)
```

```{r pressure, echo=TRUE, fig.cap="KEGG Native Enriched Pathway Plot", out.width = '100%'}
knitr::include_graphics("hsa03015.pathview.png")
```
