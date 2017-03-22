# Load libraries, register cores, set seed
library(data.table)
library(tximport)
library(limma)
library(edgeR)
library(qusage)
library(dplyr)
library(doMC)
registerDoMC(4)
set.seed(123)

# Import data
anno <- fread('./Data/Rn.anno.csv')
kegg <- readRDS('./Data/Rn.kegg.rds')
clin <- fread('./Data/Rn.Clinical.csv') 
t2g <- fread('./Data/Rn83.t2g.csv')
files <- file.path('./Data/Counts/Rat', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')

# Filter, transform counts
keep <- rowSums(cpm(txi$counts) >= 1L) >= 3L
y <- DGEList(txi$counts[keep, ])
y <- calcNormFactors(y)
overlap <- sapply(kegg$p2g, function(g) sum(g %in% rownames(y)))
kegg$p2g <- kegg$p2g[overlap > 1L]

# Results function
res <- function(coef) {

  # Genes
  topTable(fit, coef = coef, number = Inf, sort.by = 'none') %>%
    rename(AvgExpr = AveExpr,
           p.value = P.Value,
               FDR = adj.P.Val) %>%
    mutate(EnsemblID = rownames(v)) %>%
    arrange(p.value) %>%
    mutate(Rank = row_number()) %>%
    inner_join(anno, by = 'EnsemblID') %>%
    select(Rank, EnsemblID, GeneSymbol, Description, 
           AvgExpr, logFC, p.value, FDR) %>%
    fwrite(paste0('./Results/Rat/', coef, '.Genes.csv'))

  # Pathways
  resid_mat <- residuals(fit, v$E)
  mean <- fit$coefficients[, coef] 
  SD <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  sd.alpha <- SD / (fit$sigma * fit$stdev.unscaled[, coef])
  sd.alpha[is.infinite(sd.alpha)] <- 1L
  dof <- fit$df.total
  res <- newQSarray(mean = mean, SD = SD, sd.alpha = sd.alpha, dof = dof,
                    labels = rep('resid', ncol(resid_mat)))
  res <- aggregateGeneSet(res, kegg$p2g, 2L^14L)           
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE)    
  qsTable(res, number = Inf, sort.by = 'p') %>% 
    mutate(Rank = row_number()) %>%
    rename(Pathway = pathway.name,
             logFC = log.fold.change,
           p.value = p.Value) %>%
    inner_join(kegg$anno, by = 'Pathway') %>%
    select(Rank, Pathway, Description, logFC, p.value, FDR) %>%
    fwrite(paste0('./Results/Rat/', coef, '.Pathways.csv'))

}

# Normoxia vs. hypoxia
clin <- clin %>%
  mutate(Condition = relevel(as.factor(Condition), ref = 'Normoxia'))
des <- model.matrix(~ Condition, data = clin)
colnames(des) <- gsub('Condition', '', colnames(des))
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))
res('Hypoxia')

# Knock downs
clin <- clin %>%
  mutate(Condition = relevel(Condition, ref = 'Scr_KD'))
des <- model.matrix(~ Condition, data = clin)
colnames(des) <- gsub('Condition', '', colnames(des))
coefs <- colnames(des)[3:6]
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))
foreach (coef = coefs) %dopar% res(coef)


