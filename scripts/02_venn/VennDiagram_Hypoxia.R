# Load packages
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)

# Load DESeq2 datasets
GSE253274_diffexp <- read.csv("outputs/DESeq2/GSE253274.csv")
GSE239806_diffexp <- read.csv("outputs/DESeq2/GSE239806.csv")

# Filter significantly expressed genes
GSE253274_diffexp <- GSE253274_diffexp |>
    filter(padj < 0.05 & abs(log2FoldChange) > 1)

GSE239806_diffexp <- GSE239806_diffexp |>
    filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Create a list of gene set
gene_sets <- list(
    GSE253274 = GSE253274_diffexp$Gene_ID,
    GSE239806 = GSE239806_diffexp$Gene_ID
)
png("figures/Venn/VennDiagram_Hypoxia.png", units = "in", width = 6, height = 6, res = 2000)
venn_plot <- venn.diagram(
    x = gene_sets,
    filename = NULL,
    fill = c("#bcbddc", "#addd8e"),
    alpha = 0.6,
    cex = 0.8,
    cat.cex = 1,
    cat.col = c("#756bb1", "#2ca25f"),
    cat.pos = c(0, 0),
    cat.dist = c(0.02 , 0.02),
    margin = 0.1,
    lwd = 2,
    scaled = FALSE,
    disable.logging = TRUE
)
grid.draw(venn_plot)
dev.off()


# Retrieve intersection of three genes
intersect <- process_region_data(Venn(gene_sets))

# Extract genes that are common across all studies
intersect_genes <- intersect |>
    filter(id == "1/2") |>
    pull(item) |>
    as.data.frame()
colnames(intersect_genes) <- "Gene_ID"

write.csv(intersect_genes, "outputs/Venn/intersect_genes_hypoxia.csv", row.names = FALSE)

