# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(edgeR)
library(purrr)
library(NMF)
library(RColorBrewer)
library(ggsci)
library(dplyr)

### HUMANS ###

# Import data
cilia <- readRDS('./Data/Hs.cilia.rds')
t2g <- fread('./Data/Hs79.t2g.csv')
clin <- fread('./Data/Hs.Clinical.csv') %>% filter(Paired == 'Paired') 

# Tissue
files <- file.path('./Data/Counts/Human', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Subject + Tissue)

# Filter, transform counts
dds <- estimateSizeFactors(dds)
keep <- rowSums(cpm(counts(dds)) > 1) >= 3
mat <- counts(dds, normalized = TRUE)[keep, ]
mat <- cpm(mat, log = TRUE, prior.count = 1)
colnames(mat) <- clin$Sample

# Isolate cilia-related genes
cilia_genes <- map_df(seq_along(cilia), function(p) {
  data_frame(Gene = cilia[[p]])
}) %>% unique()

# Identify top 100 cilia-related genes
top <- fread('./Results/Human/Tissue.Tumour_vs_Adrenal.Genes.csv') %>%
  filter(EnsemblID %in% cilia_genes$Gene)
top_cilia_genes <- top$EnsemblID[seq_len(100)]

# Heatmap
rb <- colorRampPalette(rev(brewer.pal(10, 'RdBu')))(n = 256)
top_mat <- mat[top_cilia_genes, ]
aheatmap(top_mat, distfun = 'pearson', scale = 'row', col = rb, hclustfun = 'average',
         main = 'Top 100 Cilia-Related Genes', annCol = list('Tissue' = clin$Tissue),
         annColors = list(pal_npg()(2)[c(2, 1)]), border_color = 'grey60')

# Pathways
eigengenes <- function(dat, pathways) {
  out <- matrix(nrow = length(pathways), ncol = ncol(dat), 
                dimnames = list(names(pathways), colnames(dat)))
  for (m in seq_along(pathways)) {
    mat <- dat[rownames(dat) %in% pathways[[m]], ]
    pca <- prcomp(t(mat))
    out[m, ] <- pca$x[, 1]
  }
  return(out)
}
pathway_mat <- eigengenes(mat, cilia)
aheatmap(pathway_mat, distfun = 'pearson', scale = 'row', col = rb, hclustfun = 'average',
         main = 'Cilia Pathway Eigengenes', annCol = list('Tissue' = clin$Tissue),
         annColors = list(pal_npg()(2)[c(2, 1)]), border_color = 'grey60')



### RATS ###

# Import data
cilia <- readRDS('./Data/Rn.cilia.rds')
kegg <- readRDS('./Data/Rn.kegg.rds')
clin <- fread('./Data/Rn.Clinical.csv') 
t2g <- fread('./Data/Rn83.t2g_Symbol.csv')
files <- file.path('./Data/Counts/Rat', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Condition)

# Filter, transform counts
dds <- estimateSizeFactors(dds)
keep <- rowSums(cpm(counts(dds)) > 1) >= 3
mat <- counts(dds, normalized = TRUE)[keep, ]
mat <- cpm(mat, log = TRUE, prior.count = 1)
colnames(mat) <- clin$Sample

# Isolate cilia-related genes
cilia_genes <- map_df(seq_along(cilia), function(p) {
  data_frame(Gene = cilia[[p]])
}) %>% unique()

# Identify top 100 cilia-related genes, hypoxia vs. normoxia
top <- fread('./Results/Rat/Hypoxia_vs_Normoxia.Genes.csv') %>%
  filter(GeneSymbol %in% cilia_genes$Gene)
top_cilia_genes <- top$GeneSymbol[seq_len(100)]

# Heatmap
top_mat <- mat[top_cilia_genes, clin$Condition %in% c('Hypoxia', 'Normoxia')]
aheatmap(top_mat, distfun = 'pearson', scale = 'row', col = rb, hclustfun = 'average',
         main = 'Top 100 Cilia-Related Genes', cellwidth = 20, cexCol = 1,
         annCol = list('Condition' = clin[Condition %in% c('Hypoxia', 'Normoxia'), Condition]),
         annColors = list(pal_npg()(2)), border_color = 'grey60')

# Identify significantly enriched KEGG pathways
kegg_mat <- eigengenes(mat, kegg[[1]])
anno <- kegg[[2]] %>% filter(Pathway %in% rownames(kegg_mat))
rownames(kegg_mat) <- anno$Description
top <- fread('./Results/Rat/Hypoxia_vs_Normoxia.KEGGpathways.csv') %>%
  filter(FDR <= 0.01) # Alternatively: top 100 by F-test and include all 12 samples?
tmp <- kegg_mat[top$Description, clin$Condition %in% c('Hypoxia', 'Normoxia')]
aheatmap(tmp, distfun = 'pearson', scale = 'row', col = rb, hclustfun = 'average',
         main = 'Enriched KEGG Pathway Eigengenes', cellwidth = 20, cexCol = 1,
         annCol = list('Condition' = clin[Condition %in% c('Hypoxia', 'Normoxia'), Condition]),
         annColors = list(pal_npg()(2)), border_color = 'grey60')


