# Load libraries, set seed
library(data.table)
library(tximport)
library(limma)
library(edgeR)
library(qusage)
library(dplyr)
set.seed(123)

# Import data
clin <- fread('./Data/Clinical_Hs.csv') %>% filter(Paired == 'Paired')
t2g <- fread('./Data/Hs79.t2g.csv')
# t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
files <- file.path('./Data/Counts/Human', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')

# Filter, transform counts
keep <- rowSums(cpm(txi$counts) >= 1L) >= 2L
y <- DGEList(txi$counts[keep, , drop = FALSE])
y <- calcNormFactors(y)

# Fit model
des <- model.matrix(~ Subject + Tissue, data = clin)
colnames(des)[grepl('Tumour', colnames(des))] <- 'Tumour'
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))

# Export results
topTable(fit, coef = 'Tumour', number = Inf, sort.by = 'none') %>%
  rename(AvgExpr = AveExpr,
         p.value = P.Value,
             FDR = adj.P.Val) %>%
  mutate(Gene = rownames(v)) %>%
  arrange(p.value) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, Gene, AvgExpr, logFC, p.value, FDR) %>%
  fwrite('./Results/Human/Tumour.Genes.csv')

# KEGG pathways
kegg <- read.gmt('./Data/c2.cp.kegg.v5.2.symbols.gmt')
resid_mat <- residuals(fit, v$E)
mean <- fit$coefficients[, 'Tumour'] 
SD <- sqrt(fit$s2.post) * fit$stdev.unscaled[, 'Tumour']
sd.alpha <- SD / (fit$sigma * fit$stdev.unscaled[, 'Tumour'])
sd.alpha[is.infinite(sd.alpha)] <- 1L
dof <- fit$df.total

# Run QuSAGE functions
res <- newQSarray(mean = mean,                       # Create QSarray obj
                    SD = SD,
              sd.alpha = sd.alpha,
                   dof = dof,
                labels = rep('resid', ncol(fit)))
res <- aggregateGeneSet(res, kegg, 2L^14L)           # PDF per gene set
res <- calcVIF(resid_mat, res, useCAMERA = FALSE)    # VIF on resid_mat

# Export
qsTable(res, number = Inf, sort.by = 'p') %>%
  rename(Pathway = pathway.name,
           logFC = log.fold.change,
         p.value = p.Value) %>%
  fwrite('./Results/Human/Tumour.Pathways.csv')


