# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(dplyr)

# Import data
clin <- fread('./Data/Clinical_Rn.csv') %>%
  mutate(Condition = relevel(as.factor(Condition), ref = 'Normoxia'))
t2g <- fread('./Data/Rn83.t2g.csv')
files <- file.path('./Data/Counts/Rat', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread)

# Filter, design
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Condition)
dds <- estimateSizeFactors(dds)
mat <- counts(dds, normalized = TRUE)
dds <- dds[rowSums(mat >= 1) >= 3, , drop = FALSE]
des <- model.matrix(~ Condition, data = clin)
colnames(des) <- gsub('Condition', '', colnames(des))
coefs <- colnames(des)[-1]

# DESeq
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, modelMatrix = des)
dds <- nbinomWaldTest(dds, modelMatrix = des, betaPrior = FALSE)

# Export results
for (coef in coefs) {
  as.data.frame(results(dds, name = coef, filterfun = ihw)) %>%
    mutate(Gene = rownames(dds),
        AvgExpr = log2(baseMean)) %>%
    rename(logFC = log2FoldChange,
         p.value = pvalue,
         q.value = padj) %>%
    na.omit() %>%
    arrange(p.value) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, Gene, AvgExpr, logFC, p.value, q.value) %>%
    fwrite(paste0('./Results/Rat/DESeq2/', paste0(coef, '.Genes.txt')), sep = '\t')
}