wd <- getwd()
setwd(wd)

# Read files
rhino <- read.csv(file = "rhino_top.csv", sep = ",",row.name = 1, header = TRUE)
RSV <- read.csv(file = "RSV_top.csv", sep = ",",row.name = 1, header = TRUE)
infA <- read.csv(file = "infA_top.csv", sep = ",",row.name = 1, header = TRUE)
all <- read.csv(file = "all_top.csv", sep = ",",row.name = 1, header = TRUE)

library(biomaRt)

# ---------------------------------------------------------------------------------------------
# Create TABLE (GENE SYMBOL, P-VALUE) FOR EACH DATASET
# ---------------------------------------------------------------------------------------------

# All BATCHES

#convert probeID to gene symbols
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup_all <- getBM(
  mart = mart,
  attributes = c(
    "affy_hg_u133a_2",
    "external_gene_name"),
  filter = "affy_hg_u133a_2",
  values = all$probeID
)

head(annotLookup_all, 20)

indicesLookup_all <- match(all$probeID, annotLookup_all$affy_hg_u133a_2)

genesymb_all <- annotLookup_all[indicesLookup_all, "external_gene_name"]
all$gene.symb <- genesymb_all

all_table <- data.frame(all$gene.symb, all$pvalue)
#remove rows with NA
all_table <- all_table[complete.cases(all_table),]



# RHINO

#convert probeID to gene symbols
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup_rhino <- getBM(
  mart = mart,
  attributes = c(
    "affy_hg_u133a_2",
    "external_gene_name"),
  filter = "affy_hg_u133a_2",
  values = rhino$probeID
)

head(annotLookup_rhino, 20)

indicesLookup_rhino <- match(rhino$probeID, annotLookup_rhino$affy_hg_u133a_2)

genesymb_rhino <- annotLookup_rhino[indicesLookup_rhino, "external_gene_name"]
rhino$gene.symb <- genesymb_rhino

rhino_table <- data.frame(rhino$gene.symb, rhino$pvalue)
#remove rows with NA
rhino_table <- rhino_table[complete.cases(rhino_table),]


# RSV

#convert probeID to gene symbols
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup_RSV <- getBM(
  mart = mart,
  attributes = c(
    "affy_hg_u133a_2",
    "external_gene_name"),
  filter = "affy_hg_u133a_2",
  values = RSV$probeID
)

head(annotLookup_RSV, 20)

indicesLookup_RSV <- match(RSV$probeID, annotLookup_RSV$affy_hg_u133a_2)

genesymb_RSV <- annotLookup_RSV[indicesLookup_RSV, "external_gene_name"]
RSV$gene.symb <- genesymb_RSV

RSV_table <- data.frame(RSV$gene.symb, RSV$pvalue)
#remove rows with NA
RSV_table <- RSV_table[complete.cases(RSV_table),]


# infA

#convert probeID to gene symbols
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup_infA <- getBM(
  mart = mart,
  attributes = c(
    "affy_hg_u133a_2",
    "external_gene_name"),
  filter = "affy_hg_u133a_2",
  values = infA$probeID
)

head(annotLookup_infA, 20)

indicesLookup_infA <- match(infA$probeID, annotLookup_infA$affy_hg_u133a_2)

genesymb_infA <- annotLookup_infA[indicesLookup_infA, "external_gene_name"]
infA$gene.symb <- genesymb_infA

infA_table <- data.frame(infA$gene.symb, infA$pvalue)
#remove rows with NA
infA_table <- infA_table[complete.cases(infA_table),]

# ---------------------------------------------------------------------------------------------
# NETWORK ANALYSIS PATHFINDR FOR EACH DATASET
# ---------------------------------------------------------------------------------------------

library("KEGGREST")
library("KEGGgraph")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("pathfindR")

# ALL

#pathway enrichment
all_demo <- run_pathfindR(all_table, visualize_enriched_terms = FALSE)
#cluster togheter similar enriched terms
all_demo_cluster <- cluster_enriched_terms(all_demo)
term_gene_graph(all_demo)


# RHINO

#pathway enrichment
rhino_demo <- run_pathfindR(rhino_table, visualize_enriched_terms = FALSE)
#cluster togheter similar enriched terms
rhino_demo_cluster <- cluster_enriched_terms(rhino_demo)
term_gene_graph(rhino_demo)

# RSV

#pathway enrichment
RSV_demo <- run_pathfindR(RSV_table, visualize_enriched_terms = FALSE)
#cluster togheter similar enriched terms
RSV_demo_cluster <- cluster_enriched_terms(RSV_demo)
term_gene_graph(RSV_demo)

# INFA

#pathway enrichment
infA_demo <- run_pathfindR(infA_table, visualize_enriched_terms = FALSE)
#cluster togheter similar enriched terms
infA_demo_cluster <- cluster_enriched_terms(infA_demo)
term_gene_graph(infA_demo)


