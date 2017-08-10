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
files <- file.path('./Data/Counts/Human', clin$Sample, 'abundance.h5')
names(files) <- clin$Sample
txi <- tximport(files, type = 'kallisto', tx2gene = t2g)
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Subject + Tissue)

# Filter, transform counts
dds <- estimateSizeFactors(dds)
keep <- rowSums(cpm(counts(dds)) > 1) >= 3
mat <- counts(dds, normalized = TRUE)[keep, ]
mat <- cpm(mat, log = TRUE, prior.count = 1)

# Isolate cilia-related genes
cilia_genes <- seq_along(cilia) %>%
  map_df(~ data_frame(Gene = cilia[[.x]])) %>% 
  unique(.)

# Identify top 100 cilia-related genes
top <- fread('./Results/Human/Tissue.Tumour_vs_Adrenal.Genes.csv') %>%
  filter(EnsemblID %in% cilia_genes$Gene)
top_cilia_genes <- top$EnsemblID[seq_len(100)]

# Heatmap
rb <- colorRampPalette(rev(brewer.pal(10, 'RdBu')))(256)
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
aheatmap(pathway_mat, distfun = 'pearson', scale = 'row', col = rb, 
         border_color = 'grey60', hclustfun = 'average', 
         main = 'Cilia Pathway Eigengenes', 
         annCol = list('Tissue' = clin$Tissue), 
         annColors = list(pal_d3()(2)))

# For particular pathways
anno <- fread('./Data/Hs.anno.csv')
df <- seq_along(cilia) %>%
  map_df(~ data_frame(Pathway = names(cilia)[.x],
                    EnsemblID = cilia[[.x]])) %>%
  inner_join(anno, by = 'EnsemblID') %>%
  select(Pathway, EnsemblID, GeneSymbol)

# For example
tmp <- df %>% filter(Pathway == 'HALLMARK_HEDGEHOG_SIGNALING') 
tmp_mat <- mat[match(tmp$EnsemblID, rownames(mat)), ] %>% na.omit(.)
tmp <- tmp %>% filter(EnsemblID %in% rownames(mat))  # Some don't pass the filter
rownames(tmp_mat) <- tmp$GeneSymbol                  # Makes reading the heatmap easier
aheatmap(tmp_mat, distfun = 'pearson', scale = 'row', col = rb, hclustfun = 'average',
         main = 'HALLMARK_HEDGEHOG_SIGNALING', annCol = list('Tissue' = clin$Tissue),
         annColors = list(pal_npg()(2)[c(2, 1)]), border_color = 'grey60')

# As a function:
pathway_heatmap <- function(pathway) {
  
  tmp <- df %>% filter(Pathway == pathway)
  tmp_mat <- mat[match(tmp$EnsemblID, rownames(mat)), ] %>% na.omit(.)
  tmp <- tmp %>% filter(EnsemblID %in% rownames(mat))
  rownames(tmp_mat) <- tmp$GeneSymbol
  aheatmap(tmp_mat, distfun = 'pearson', scale = 'row', col = rb, hclustfun = 'average',
           annCol = list('Tissue' = clin$Tissue), 
           annColors = list(pal_d3()(2)), border_color = 'grey60')
  
}

# Add columns to pathway CSV
genes <- fread('./Results/Human/Tissue.Tumour_vs_Adrenal.Genes.csv')
pways <- fread('./Results/Human/Tissue.Tumour_vs_Adrenal.Pathways.csv')
up <- down <- sig_up <- sig_down <- double(length = nrow(pways))
for (i in seq_len(nrow(pways))) {
  tmp <- df %>% filter(Pathway == pways$Pathway[i])
  genes_tmp <- genes %>% filter(EnsemblID %in% tmp$EnsemblID)
  up[i] <- sum(genes_tmp$logFC > 0L)
  down[i] <- sum(genes_tmp$logFC < 0L)
  sig_up[i] <- sum(genes_tmp$q.value <= 0.01 & genes_tmp$logFC >= 2L)
  sig_down[i] <- sum(genes_tmp$q.value <= 0.01 & genes_tmp$logFC <= -2L)
}
pways <- pways %>%
  mutate(Up = up, Down = down, Sig_Up = sig_up, Sig_Down = sig_down)
fwrite(pways, './Results/Human/Tissue.Tumour_vs_Adrenal.Pathways.csv')


### RATS ###

# Import data
cilia <- readRDS('./Data/Rn.cilia.rds')
kegg <- readRDS('./Data/Rn.kegg.rds')
clin <- fread('./Data/Rn.Clinical.csv') %>% 
  filter(Condition %in% c('IFT88_KD', 'Scrambled'))
t2g <- fread('./Data/Rn83.t2g_Symbol.csv')
files <- file.path('./Data/Counts/Rat', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Condition)

# Filter, transform counts
dds <- estimateSizeFactors(dds)
keep <- rowSums(cpm(counts(dds)) > 1) >= 3
mat <- counts(dds, normalized = TRUE)[keep, ]
mat <- cpm(mat, log = TRUE, prior.count = 1)
colnames(mat) <- clin$Sample

# For particular pathways
df <- fread('IPA.csv')
logs <- fread('./Data/Hs2Rn.csv')
pathway_heatmap <- function(pathway) {
  tmp <- df %>% filter(Function == pathway)
  x <- tmp$Molecules %>%
    strsplit(',') %>%
    unlist(.)
  tmp <- data.frame(Hs_Gene = x)
  y <- tmp %>% inner_join(logs, by = 'Hs_Gene')
  tmp_mat <- mat[match(y$Rn_Gene, rownames(mat)), ] %>% na.omit(.)
  aheatmap(tmp_mat, distfun = 'pearson', scale = 'row', col = rb, 
           hclustfun = 'average', annCol = list('Condition' = clin$Condition), 
           annColors = list(pal_d3()(2)), border_color = 'grey60',
           cellwidth = 15)
}







# Isolate cilia-related genes
cilia_genes <- seq_along(cilia) %>% 
  map_df(~ data_frame(Gene = cilia[[.x]])) %>% 
  unique(.)

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

# For particular pathways
df <- seq_along(kegg$p2g) %>% 
  map_df(~ data_frame(Pathway = kegg$anno[.x, 1],
                  Description = kegg$anno[.x, 2],
                         Gene = kegg$p2g[[.x]]))
clin <- clin %>% filter(Condition %in% c('IFT88_KD', 'Scrambled'))
mat <- mat[, match(clin$Sample, colnames(mat))]
pathway_heatmap <- function(pathway) {
  tmp <- df %>% filter(Pathway == pathway)
  tmp_mat <- mat[match(tmp$Gene, rownames(mat)), ] %>% na.omit(.)
  aheatmap(tmp_mat, distfun = 'pearson', hclustfun = 'average',
           scale = 'row', col = rb, cellwidth = 20,
           annCol = list('Condition' = clin$Condition), 
           annColors = list(pal_d3()(2)), border_color = 'grey60')
}


  
  
  
  
  
  
  
  
  
  
  






