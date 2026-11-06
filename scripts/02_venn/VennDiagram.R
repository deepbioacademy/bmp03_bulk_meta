# Load packages
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)

# Load DESeq2 datasets
GSE262612_diffexp <- read.csv("outputs/DESeq2/GSE262612.csv")
GSE253274_diffexp <- read.csv("outputs/DESeq2/GSE253274.csv")
GSE239806_diffexp <- read.csv("outputs/DESeq2/GSE239806.csv")
GSE135668_diffexp <- read.csv("outputs/DESeq2/GSE135668.csv")

# Filter significantly expressed genes
GSE262612_diffexp <- GSE262612_diffexp |>
    filter(padj < 0.05 & abs(log2FoldChange) > 1)

GSE253274_diffexp <- GSE253274_diffexp |>
    filter(padj < 0.05 & abs(log2FoldChange) > 1)

GSE239806_diffexp <- GSE239806_diffexp |>
    filter(padj < 0.05 & abs(log2FoldChange) > 1)

GSE135668_diffexp <- GSE135668_diffexp |>
    filter(padj < 0.05 & abs(log2FoldChange) > 1)


# Create a list of gene set
gene_sets <- list(
    GSE262612 = GSE262612_diffexp$Gene_ID,
    GSE253274 = GSE253274_diffexp$Gene_ID,
    GSE239806 = GSE239806_diffexp$Gene_ID,
    GSE135668 = GSE135668_diffexp$Gene_ID
)

# Create VennDiagram
png("figures/Venn/VennDiagram.png", units = "in", width = 6, height = 6, res = 2000)
venn_plot <- venn.diagram(
    x = gene_sets,
    filename = NULL,
    fill = c("#99d8c9", "#addd8e", "#bcbddc", "#fec44f"),
    alpha = 0.6,
    cex = 0.8,
    cat.cex = 1,
    cat.col = c("#2ca25f", "#31a354", "#756bb1", "#e6550d"),
    cat.pos = c(-20, 20, -20, 20),
    cat.dist = c(0.22, 0.22, 0.1, 0.09),
    margin = 0.1,
    lwd = 2,
    disable.logging = TRUE
)
grid.draw(venn_plot)
dev.off()

# Retrieve intersection of three genes
intersect <- process_region_data(Venn(gene_sets))

# Extract genes that are common across all studies
intersect_genes <- intersect |>
    filter(id == "1/2/3/4") |>
    pull(item) |>
    as.data.frame()
colnames(intersect_genes) <- "Gene_ID"

write.csv(intersect_genes, "outputs/Venn/intersect_genes.csv", row.names = FALSE)
