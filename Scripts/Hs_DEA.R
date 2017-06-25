# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(DESeq2)
library(edgeR)
library(qusage)
library(qvalue)
library(dplyr)
library(doMC)
registerDoMC(4)
set.seed(123)

# Import data
anno <- fread('./Data/Hs.anno.csv')
cilia <- readRDS('./Data/Hs.cilia.rds')
t2g <- fread('./Data/Hs79.t2g.csv')
clin <- fread('./Data/Hs.Clinical.csv', stringsAsFactors = TRUE) 

# Results funcion
res <- function(dds, cont) {
  
  # Genes
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  dds %>%
    results(contrast = cont, filterfun = ihw, alpha = 0.01, tidy = TRUE) %>%
    na.omit(.) %>%
    mutate(AvgExpr = log2(baseMean / 1e6L)) %>%
    rename(EnsemblID = row,
               logFC = log2FoldChange,
             p.value = pvalue,
             q.value = padj) %>%
    arrange(p.value) %>%
    inner_join(anno, by = 'EnsemblID') %>%
    select(EnsemblID, GeneSymbol, Description, 
           AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Human/', cont[1], '.', 
                  cont[2], '_vs_', cont[3], '.Genes.csv'))
  
  # Pathways
  n <- ncol(dds)
  p <- ncol(model.matrix(design(dds), data = colData(dds)))
  max_rps <- max(table(pheno[[cont[1]]]))
  keep <- rowSums(cpm(counts(dds)) > 1) >= max_rps
  dds <- dds[keep, ]
  dds_res <- results(dds, contrast = cont, independentFiltering = FALSE)
  mean <- dds_res$log2FoldChange
  SD <- dds_res$lfcSE
  sd.alpha <- rep(1, times = nrow(dds))
  dof <- rep((ncol(dds) - p), times = nrow(dds)) 
  names(mean) <- names(SD) <- names(sd.alpha) <- names(dof) <- rownames(dds)
  cnts <- dds %>% 
    counts(normalized = TRUE) %>% 
    cpm(log = TRUE, prior.count = 1)
  signal_mat <- (assays(dds)[['mu']] / normalizationFactors(dds)) %>%
    cpm(log = TRUE, prior.count = 1)
  resid_mat <- cnts - signal_mat
  overlap <- sapply(cilia, function(i) sum(i %in% rownames(dds)))
  cilia <- cilia[overlap > 1]
  res <- newQSarray(mean = mean, SD = SD, sd.alpha = sd.alpha, dof = dof,
                    labels = rep('resid', n))
  res <- aggregateGeneSet(res, cilia, 2^14)     
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE) 
  qsTable(res, number = Inf, sort.by = 'p') %>%
    rename(Pathway = pathway.name,
             logFC = log.fold.change,
           p.value = p.Value) %>%
    mutate(q.value = qvalue(p.value)$qvalues) %>%
    select(-FDR) %>%
    fwrite(paste0('./Results/Human/', cont[1], '.',
                  cont[2], '_vs_', cont[3], '.Pathways.csv'))
  
}

# Tissue
pheno <- clin %>% filter(Paired == 'Paired')
files <- file.path('./Data/Counts/Human', pheno$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)
dds <- DESeqDataSetFromTximport(txi, colData = pheno, design = ~ Subject + Tissue)
contrast <- c('Tissue', 'Tumour', 'Adrenal')
res(dds, contrast)

# Mutation: vs. Neg
pheno <- clin %>% filter(Tissue == 'Tumour')
files <- file.path('./Data/Counts/Human', pheno$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)
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


