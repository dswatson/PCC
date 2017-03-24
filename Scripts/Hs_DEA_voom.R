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
anno <- fread('./Data/Hs.anno.csv')
cilia <- readRDS('./Data/Hs.cilia.rds')
t2g <- fread('./Data/Hs79.t2g.csv')
clin <- fread('./Data/Hs.Clinical.csv')
files <- file.path('./Data/Counts/Human', clin$Sample, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, reader = fread, 
                countsFromAbundance = 'lengthScaledTPM')
colnames(txi$counts) <- clin$Sample

# Results function
res <- function(comp, coef) {
  
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
    fwrite(paste0('./Results/Human/voom/', comp, '.', coef, '.Genes.csv'))
  
  # Pathways
  mean <- fit$coefficients[, coef] 
  SD <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  sd.alpha <- SD / (fit$sigma * fit$stdev.unscaled[, coef])
  dof <- fit$df.total
  resid_mat <- residuals(fit, v$E)
  overlap <- sapply(cilia, function(g) sum(g %in% rownames(fit)))
  cilia <- cilia[overlap > 1L]
  res <- newQSarray(mean = mean, SD = SD, sd.alpha = sd.alpha, dof = dof,
                    labels = rep('resid', ncol(resid_mat)))
  res <- aggregateGeneSet(res, cilia, 2L^14L)         
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE)
  qsTable(res, number = Inf, sort.by = 'p') %>%
    mutate(Rank = row_number()) %>%
    rename(Pathway = pathway.name,
           logFC = log.fold.change,
           p.value = p.Value) %>%
    select(Rank, Pathway, logFC, p.value, FDR) %>%
    fwrite(paste0('./Results/Human/voom/', comp, '.', coef, '.Pathways.csv'))
  
}

# Tissue
pheno <- clin %>% filter(Paired == 'Paired')
cnts <- txi$counts[, match(pheno$Sample, colnames(txi$counts))]
max_rps <- max(table(pheno$Tissue))
keep <- rowSums(cpm(cnts) > 1L) >= max_rps
y <- DGEList(cnts[keep, ])
y <- calcNormFactors(y)
des <- model.matrix(~ Subject + Tissue, data = pheno)
colnames(des)[grepl('Tumour', colnames(des))] <- 'Tumour_vs_Adrenal'
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))
res(comp = 'Tissue', coef = 'Tumour_vs_Adrenal')

# Mutation: vs. Neg
pheno <- clin %>% 
  filter(Tissue == 'Tumour') %>%
  mutate(Mutation = relevel(as.factor(Mutation), ref = 'Negative'),
         Pseudohypoxia = relevel(as.factor(Pseudohypoxia), ref = 'SDHB'),
         Location = relevel(as.factor(Location), ref = 'Adrenal'))
cnts <- txi$counts[, match(pheno$Sample, colnames(txi$counts))]
max_rps <- max(table(pheno$Mutation))
keep <- rowSums(cpm(cnts) > 1L) >= max_rps
y <- DGEList(cnts[keep, ])
y <- calcNormFactors(y)
des <- model.matrix(~ Mutation, data = pheno)
colnames(des) <- gsub('Mutation', '', colnames(des))
colnames(des)[2:4] <- paste0(colnames(des)[2:4], '_vs_Neg')
coefs <- colnames(des)[2:4]
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))
foreach (coef = coefs) %dopar% res(comp = 'Mutation', coef = coef)

# Mutation: pseudohypoxia
des <- model.matrix(~ Pseudohypoxia, data = pheno)
colnames(des)[3] <- 'VHL_vs_SDHB'
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))
res(comp = 'Mutation', coef = 'VHL_vs_SDHB')

# Location
des <- model.matrix(~ Location, data = pheno)
colnames(des) <- gsub('Location', '', colnames(des))
colnames(des)[2:3] <- paste0(colnames(des)[2:3], '_vs_Adrenal')
coefs <- colnames(des)[2:3]
v <- voom(y, des)
fit <- eBayes(lmFit(v, des))
foreach(coef = coefs) %dopar% res(comp = 'Location', coef = coef)


