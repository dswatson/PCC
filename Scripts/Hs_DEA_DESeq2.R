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
anno <- fread('./Data/Hs.anno.csv')
cilia <- readRDS('./Data/cilia.rds')
t2g <- fread('./Data/Hs79.t2g.csv')
clin <- fread('./Data/Hs.Clinical.csv', stringsAsFactors = TRUE) 

# Results funcion
res <- function(dds, cont) {
  
  # Genes
  dds <- dds[rowSums(counts(dds)) > 1L, ]
  dds <- DESeq(dds)
  as.data.frame(results(dds, filterfun = ihw, alpha = 0.05, contrast = cont)) %>%
    mutate(EnsemblID = rownames(dds),
             AvgExpr = log2(baseMean)) %>%
    na.omit() %>%
    rename(logFC = log2FoldChange,
         p.value = pvalue,
             FDR = padj) %>%
    arrange(p.value) %>%
    mutate(Rank = row_number()) %>%
    inner_join(anno, by = 'EnsemblID') %>%
    select(Rank, EnsemblID, GeneSymbol, Description, 
           AvgExpr, logFC, p.value, FDR) %>%
    fwrite(paste0('./Results/Human/DESeq2/', cont[1], '.', 
                  cont[2], '_vs_', cont[3], '.Genes.csv'))
  
  # Pathways
  n <- ncol(dds)
  p <- ncol(model.matrix(design(dds), colData(dds)))
  max_rps <- max(table(pheno[[cont[1]]]))
  keep <- rowSums(cpm(counts(dds)) > 1L) >= max_rps
  dds <- dds[keep, ]
  dds_res <- results(dds, independentFiltering = FALSE)
  mean <- dds_res$log2FoldChange
  SD <- dds_res$lfcSE
  sd.alpha <- rep(1L, times = nrow(dds))
  dof <- rep((ncol(dds) - p), times = nrow(dds)) 
  names(mean) <- names(SD) <- names(sd.alpha) <- names(dof) <- rownames(dds)
  nf <- normalizationFactors(dds) + 1L
  cnts <- log2((counts(dds) + 0.5) / nf * 1e6L)
  signal_mat <- log2((assays(dds)[['mu']] + 0.5) / nf * 1e6L)
  resid_mat <- cnts - signal_mat
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
    mutate(Rank = row_number()) %>%
    select(Rank, Pathway, logFC, p.value, FDR) %>%
    fwrite(paste0('./Results/Human/DESeq2/', cont[1], '.',
                  cont[2], '_vs_', cont[3], '.Pathways.csv'))
  
}

# Tissue
pheno <- clin %>% filter(Paired == 'Paired')
files <- file.path('./Data/Counts/Human', pheno$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)
dds <- DESeqDataSetFromTximport(txi, colData = pheno, design = ~ Subject + Tissue)
contrast <- c('Tissue', 'Tumour', 'Adrenal')
res(dds, contrast)

# Mutation: vs. Neg
pheno <- clin %>% filter(Tissue == 'Tumour')
files <- file.path('./Data/Counts/Human', pheno$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)
dds <- DESeqDataSetFromTximport(txi, colData = pheno, design = ~ Mutation)
cont_list <- list(c('Mutation', 'Pseudohypoxia', 'Negative'),
                  c('Mutation', 'Ret', 'Negative'),
                  c('Mutation', 'Unknown', 'Negative'))
foreach(contrast = cont_list) %dopar% res(dds, contrast)

# Mutation: pseudohypoxia
design(dds) <- formula(~ Pseudohypoxia)
contrast <- c('Pseudohypoxia', 'VHL', 'SDHB')
res(dds, contrast)

# Location
design(dds) <- formula(~ Location)
dds$Location <- droplevels(dds$Location)
cont_list <- list(c('Location', 'Abdominal', 'Adrenal'),
                  c('Location', 'HNPGL', 'Adrenal'))
foreach(contrast = cont_list) %dopar% res(dds, contrast)


