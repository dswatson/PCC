# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(DESeq2)
library(edgeR)
library(qusage)
library(dplyr)
library(doMC)
registerDoMC(4)
set.seed(123)

# Import data
anno <- fread('./Data/Rn.anno.csv')
cilia <- readRDS('./Data/Rn.cilia.rds')
kegg <- readRDS('./Data/Rn.kegg.rds')
clin <- fread('./Data/Rn.Clinical.csv') 
t2g <- fread('./Data/Rn83.t2g_Symbol.csv')
files <- file.path('./Data/Counts/Rat', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)

# Results function
res <- function(contrast) {
  
  # Genes
  dds <- dds[rowSums(counts(dds)) > 1L, ]
  dds <- DESeq(dds)
  as.data.frame(results(dds, filterfun = ihw, alpha = 0.05, contrast = contrast)) %>%
    mutate(GeneSymbol = rownames(dds),
              AvgExpr = log2(baseMean / 1e6L)) %>%
    na.omit() %>%
    rename(logFC = log2FoldChange,
         p.value = pvalue,
             FDR = padj) %>%
    arrange(p.value) %>%
    inner_join(anno, by = 'GeneSymbol') %>%
    select(GeneSymbol, Description, AvgExpr,
           logFC, p.value, FDR) %>%
    fwrite(paste0('./Results/Rat/', contrast[2], '_vs_', contrast[3], '.Genes.csv'))
  
  # Pathways
  n <- ncol(dds)
  p <- ncol(model.matrix(design(dds), colData(dds)))
  keep <- rowSums(cpm(counts(dds)) > 1L) >= 3L
  dds <- dds[keep, ]
  dds_res <- results(dds, filterfun = ihw, alpha = 0.05, 
                     contrast = contrast, independentFiltering = FALSE)
  mean <- dds_res$log2FoldChange
  SD <- dds_res$lfcSE
  sd.alpha <- rep(1L, times = nrow(dds))
  dof <- rep((ncol(dds) - p), times = nrow(dds)) 
  names(mean) <- names(SD) <- names(sd.alpha) <- names(dof) <- rownames(dds)
  cnts <- counts(dds, normalized = TRUE)
  cnts <- cpm(cnts, log = TRUE, prior.count = 1L)
  signal_mat <- assays(dds)[['mu']] / normalizationFactors(dds)
  signal_mat <- cpm(signal_mat, log = TRUE, prior.count = 1L)
  resid_mat <- cnts - signal_mat
  
  # Cilia
  overlap <- sapply(cilia, function(g) sum(g %in% rownames(dds)))
  cilia <- cilia[overlap > 1L]
  res <- newQSarray(mean = mean, SD = SD, sd.alpha = sd.alpha, dof = dof,
                    labels = rep('resid', n))
  res <- aggregateGeneSet(res, cilia, 2L^14L)     
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE) 
  qsTable(res, number = Inf, sort.by = 'p') %>%
    rename(Pathway = pathway.name,
             logFC = log.fold.change,
           p.value = p.Value) %>%
    fwrite(paste0('./Results/Rat/', contrast[2], '_vs_', contrast[3], 
                  '.CiliaPathways.csv'))
  
  # KEGG
  overlap <- sapply(kegg$p2g, function(g) sum(g %in% rownames(dds)))
  kegg$p2g <- kegg$p2g[overlap > 1L]
  res <- newQSarray(mean = mean, SD = SD, sd.alpha = sd.alpha, dof = dof,
                    labels = rep('resid', n))
  res <- aggregateGeneSet(res, kegg$p2g, 2L^14L)     
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE) 
  qsTable(res, number = Inf, sort.by = 'p') %>%
    rename(Pathway = pathway.name,
             logFC = log.fold.change,
           p.value = p.Value) %>%
    inner_join(kegg$anno, by = 'Pathway') %>%
    select(Pathway, Description, logFC, p.value, FDR) %>%
    fwrite(paste0('./Results/Rat/', contrast[2], '_vs_', contrast[3], 
                  '.KEGGPathways.csv'))
  
}

# Execute in parallel
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Condition)
cont_list <- list(c('Condition', 'Hypoxia', 'Normoxia'),
                  c('Condition', 'IFT88_KD', 'Scrambled'),
                  c('Condition', 'SDHB_KD', 'Scrambled'),
                  c('Condition', 'VHL_KD', 'Scrambled'))
foreach(contrast = cont_list) %dopar% res(contrast)


