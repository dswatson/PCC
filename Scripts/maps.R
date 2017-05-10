# Load libraries
library(data.table)
library(AnnotationHub)
library(ensembldb)
library(biomaRt)
library(qusage)
library(limma)
library(purrr)
library(dplyr)

# Map transcripts to gene symbols
ah <- AnnotationHub()
rn83 <- query(ah, c('ensembl', 'gtf', '83', 'rattus norvegicus'))  # Load genome build
db <- ensDbFromAH(ah[names(rn83)])		 						                 # Create ensembldb
edb <- EnsDb(db)
transcripts(edb, columns = c('tx_id', 'gene_id'), return.type = 'data.frame') %>%
  rename(ensembl_id = gene_id) %>%
  fwrite('./Data/Rn83.t2g.csv')							   

# Map Ensembl to Entrez IDs
transcripts(edb, columns = c('tx_id', 'gene_id'), return.type = 'data.frame') %>%
  rename(ensembl_id = gene_id) %>%
  fwrite('./Data/Ensembl.Rn83.t2g.csv')


### Pathway Lists ###

# Human
t2g <- fread('./Data/Hs79.t2g.csv')
e2e <- getBM(attributes = c('ensembl_gene_id', 'entrezgene'),
                filters = 'ensembl_gene_id', 
                 values = t2g$ensembl_id, 
                   mart = useDataset(dataset = 'hsapiens_gene_ensembl', 
                                        mart = useMart('ensembl'))) %>%
  rename(ensembl_id = ensembl_gene_id, entrez_id = entrezgene) %>%
  na.omit() 
kegg <- getGeneKEGGLinks() %>%
  rename(Pathway = PathwayID) %>%
  mutate(entrez_id = as.numeric(GeneID)) %>%
  inner_join(e2e, by = 'entrez_id') %>%
  data.table()
anno <- getKEGGPathwayNames(species.KEGG = 'hsa', remove.qualifier = TRUE) %>%
  rename(Pathway = PathwayID) %>%
  data.table()
kegg_list <- lapply(unique(kegg$Pathway), function(p) kegg[Pathway == p, ensembl_id])
names(kegg_list) <- unique(kegg$Pathway)
kegg_list <- list('p2g' = kegg_list, 'anno' = anno)
saveRDS(kegg_list, './Data/Hs.kegg.rds')

# Custom cilia-related pathways
msig <- read.gmt('./Data/msigdb.v5.2.entrez.gmt')
cilia <- map(peep, function(x) msig[[x]])  # peep is vector of pathway names
df <- map_df(seq_along(cilia), function(p) {
  data_frame(Pathway = names(cilia)[p],
             entrez_id = cilia[[p]])
}) %>%
  mutate(entrez_id = as.numeric(entrez_id)) %>%
  inner_join(e2e, by = 'entrez_id') %>%
  data.table()


# Rat
t2g <- fread('./Data/Rn83.t2g_Ensembl.csv')
mart <- useDataset('rnorvegicus_gene_ensembl', mart = useMart('ensembl'))
e2e <- getBM(attributes = c('ensembl_gene_id', 'entrezgene'), 
                filters = 'ensembl_gene_id', 
                 values = unique(t2g$ensembl_id), 
                   mart = useDataset(dataset = 'rnorvegicus_gene_ensembl', 
                                        mart = useMart('ensembl'))) %>%
  rename(ensembl_id = ensembl_gene_id, entrez_id = entrezgene) %>%
  na.omit()
kegg <- getGeneKEGGLinks(species.KEGG = 'rno') %>%
  rename(Pathway = PathwayID) %>%
  mutate(entrez_id = as.numeric(GeneID)) %>%
  inner_join(e2e, by = 'entrez_id') %>%
  data.table()
anno <- getKEGGPathwayNames(species.KEGG = 'rno', remove.qualifier = TRUE) %>%
  rename(Pathway = PathwayID) 
kegg_list <- lapply(unique(kegg$Pathway), function(p) kegg[Pathway == p, ensembl_id])
names(kegg_list) <- unique(kegg$Pathway)
kegg_list <- list('p2g' = kegg_list, 'anno' = anno)
saveRDS(kegg_list, './Data/Rn.kegg.rds')

# Turn pathway .tab files into well formatted .rds lists
dir <- './Data/Rn.pathways/'
files <- list.files(dir)
pways <- data.table(Pathway = NA, Gene = NA)
for (file in files) {
  pways <- read.delim(paste0(dir, file)) %>%
    select(Term, Symbol) %>%
    rename(Pathway = Term, Gene = Symbol) %>%
    rbind(pways)
}
pways <- pways %>% 
  na.omit() %>% 
  unique() %>%
  data.table()
pways <- pways[, Count := .N, by = Pathway][Count > 1L]
pway_list <- lapply(unique(pways$Pathway), function(p) pways[Pathway == p, Gene])
names(pway_list) <- unique(pways$Pathway)
saveRDS(pway_list, './Data/Rn.cilia.rds')


### Annotation Files ###

# Human
t2g <- fread('./Data/Ensembl.Hs79.Tx.csv')
mart <- useDataset(dataset = 'hsapiens_gene_ensembl', mart = useMart('ensembl'))
attr <- listAttributes(mart) # SO MUCH CRAP
getBM(filters = 'ensembl_gene_id', values = t2g$ensembl_id, mart = mart,
	  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description')) %>%
  rename(EnsemblID = ensembl_gene_id,
  	    GeneSymbol = hgnc_symbol,
  	   Description = description) %>%
  fwrite('./Data/Hs.anno.csv')

# Rat
t2g <- fread('./Data/Ensembl.Rn83.t2g.csv')
mart <- useDataset('rnorvegicus_gene_ensembl', mart = useMart('ensembl'))
getBM(filters = 'ensembl_gene_id', values = t2g$ensembl_id, mart = mart,
	  attributes = c('ensembl_gene_id', 'rgd_symbol', 'description')) %>%
  rename(EnsemblID = ensembl_gene_id,
  	    GeneSymbol = rgd_symbol,
  	   Description = description) %>%
  fwrite('./Data/Rn.anno.csv')




