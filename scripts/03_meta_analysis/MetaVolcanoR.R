# Meta-Analysis of RNA-seq data using MetaVolcanoR
library(tidyverse)
library(MetaVolcanoR)
library(biomaRt)
library(ensembldb)
library(AnnotationHub)


# Load one study
GSE262612_diffexp <- prepare_deseq2(read.csv("outputs/DESeq2/GSE262612.csv"), gene_column = "Gene_ID")
GSE253274_diffexp <- prepare_deseq2(read.csv("outputs/DESeq2/GSE253274.csv"), gene_column = "Gene_ID")
GSE239806_diffexp <- prepare_deseq2(read.csv("outputs/DESeq2/GSE239806.csv"), gene_column = "Gene_ID")
GSE135668_diffexp <- prepare_deseq2(read.csv("outputs/DESeq2/GSE135668.csv"), gene_column = "Gene_ID")

# Load all diffexp datasets
studies <- list(
    GSE262612 = GSE262612_diffexp,
    GSE253274 = GSE253274_diffexp,
    GSE239806 = GSE239806_diffexp,
    GSE135668 = GSE135668_diffexp
)

# Meta-Analysis
# Run Random Effects Model meta-analysis
rem_results <- rem_mv(diffexp = studies, metathr = 0.01)
meta_results <- rem_results@metaresult

# Offline annotation
ah <- AnnotationHub()
# Query for Oncorhynchus mykiss EnsDb
query_res <- query(ah, c("EnsDb", "Oncorhynchus mykiss"))
edb <- query_res[[length(query_res)]] # Select the latest release

# Query local database without network calls to Ensembl
annotations <- select(
    edb,
    keys = meta_results$Symbol,
    keytype = "GENEID",
    columns = c("GENEID", "GENENAME", "DESCRIPTION")
)

# Deduplicate 1-to-many mappings if present in annotations
annotations_clean <- annotations  |>
    as_tibble() |>
    group_by(GENEID) |>
    summarize(
        GENENAME = paste(unique(na.omit(GENENAME)), collapse = "; "),
        DESCRIPTION = paste(unique(na.omit(DESCRIPTION)), collapse = "; "),
        .groups = "drop"
    )

# Merge annotations into meta_results
meta_annotated <- meta_results |>
    as_tibble() |>
    left_join(annotations_clean, by = c("Symbol" = "GENEID")) |>
    relocate(GENENAME, DESCRIPTION, .after = Symbol)

meta_annotated <- meta_annotated |>
    dplyr::rename(
        Gene_ID = Symbol,
        Symbol = GENENAME
)

write.csv(meta_annotated, "outputs/Meta_Analysis/REM_Results.csv", row.names = FALSE)

# Filtering based on statistical significance
significant_genes <- meta_annotated |>
    dplyr::filter(randomP < 0.05)

# Filtering based reliable effect size estimates
key_genes <- significant_genes |>
    dplyr::filter(abs(randomSummary) >= 1)
write.csv(meta_annotated, "outputs/Meta_Analysis/Filtered_Meta_Results.csv", row.names = FALSE)

# Load Venn
intersects <- read.csv("outputs/Venn/intersect_genes.csv")
meta_annotated |>
    dplyr::filter(Gene_ID %in% intersects$Gene_ID)
