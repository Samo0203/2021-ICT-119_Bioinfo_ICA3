# Differential Gene Expression Analysis
# Leukemia: ALL vs AML

# Load required packages
library(limma)
library(ggplot2)
library(ggrepel)

# Load saved expression data
expression_data <- readRDS(
  "data/expression_data.rds"
)

# Load patient group information
patient_groups <- read.csv(
  "data/patient_groups.csv"
)

# Fix sample ID data type
patient_groups$Sample <- as.character(
  patient_groups$Sample
)

# Create group factor
group <- factor(
  patient_groups$Group,
  levels = c("ALL", "AML")
)

# Check sample alignment
stopifnot(
  identical(
    colnames(expression_data),
    patient_groups$Sample
  )
)

# Check group information
print(table(group))

# Create design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit linear model
fit <- lmFit(
  expression_data,
  design
)

# Create contrast
contrast.matrix <- makeContrasts(
  AML_vs_ALL = AML - ALL,
  levels = design
)

# Apply contrast
fit2 <- contrasts.fit(
  fit,
  contrast.matrix
)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract differential expression results
all_results <- topTable(
  fit2,
  coef = "AML_vs_ALL",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Add probe IDs
all_results$Probe <- rownames(all_results)

all_results <- all_results[
  ,
  c("Probe", setdiff(names(all_results), "Probe"))
]

# Identify significant probes
significant_results <- subset(
  all_results,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

# Create output folders
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save differential expression results
write.csv(
  all_results,
  "results/tables/differential_expression_all.csv",
  row.names = FALSE
)

write.csv(
  significant_results,
  "results/tables/significant_probes.csv",
  row.names = FALSE
)

# Display top 10 results
print(head(all_results, 10))

# Display analysis summary
cat("\nDifferential Expression Analysis Summary\n")

cat("Total probes analysed:", nrow(all_results), "\n")

cat("Significant probes:", nrow(significant_results), "\n")

cat(
  "Higher expression in AML:",
  sum(significant_results$logFC > 1),
  "\n"
)

cat(
  "Higher expression in ALL:",
  sum(significant_results$logFC < -1),
  "\n"
)

# Volcano Plot

# Classify probe expression
all_results$Regulation <- "Not Significant"

all_results$Regulation[
  all_results$adj.P.Val < 0.05 &
    all_results$logFC > 1
] <- "Higher in AML"

all_results$Regulation[
  all_results$adj.P.Val < 0.05 &
    all_results$logFC < -1
] <- "Higher in ALL"

# Calculate negative log10 adjusted p-value
all_results$negLog10AdjP <- -log10(
  pmax(
    all_results$adj.P.Val,
    .Machine$double.xmin
  )
)

# Select top 10 significant probes for labels
top_labels <- head(
  significant_results[
    order(significant_results$adj.P.Val),
  ],
  10
)

# Create volcano plot
volcano_plot <- ggplot(
  all_results,
  aes(
    x = logFC,
    y = negLog10AdjP,
    color = Regulation
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 1.5
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_text_repel(
    data = top_labels,
    inherit.aes = FALSE,
    aes(
      x = logFC,
      y = -log10(
        pmax(adj.P.Val, .Machine$double.xmin)
      ),
      label = Probe
    ),
    color = "black",
    size = 3,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "Higher in AML" = "red",
      "Higher in ALL" = "blue",
      "Not Significant" = "grey"
    )
  ) +
  
  labs(
    title = "Volcano Plot: AML vs ALL",
    subtitle = "Differential Gene Expression Analysis",
    x = "Log2 Fold Change (AML vs ALL)",
    y = "-Log10 Adjusted P-value",
    color = "Expression Group"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(size = 11),
    legend.position = "right"
  )

# Display volcano plot
print(volcano_plot)

# Save volcano plot
ggsave(
  filename = "results/figures/volcano_plot.png",
  plot = volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# Create patient groups
group <- factor(patient_groups$Group, levels = c("ALL", "AML"))

# Check patient group counts
table(group)

# Create design matrix
design <- model.matrix(~0 + group)

# Rename design matrix columns
colnames(design) <- levels(group)

# Display design matrix
design

# Fit linear model
fit <- lmFit(expression_data, design)

# Define comparison: AML vs ALL
contrast.matrix <- makeContrasts(
  AML_vs_ALL = AML - ALL,
  levels = design
)

# Apply contrast
fit2 <- contrasts.fit(fit, contrast.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results for all probes
all_results <- topTable(
  fit2,
  coef = "AML_vs_ALL",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Add probe IDs as a column
all_results$Probe <- rownames(all_results)

# Move Probe column to the first position
all_results <- all_results[
  , c("Probe", setdiff(colnames(all_results), "Probe"))
]

# Identify significantly differentially expressed probes
significant_results <- subset(
  all_results,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

# Save all results
write.csv(
  all_results,
  "results/tables/differential_expression_all.csv",
  row.names = FALSE
)

# Save significant results
write.csv(
  significant_results,
  "results/tables/significant_probes.csv",
  row.names = FALSE
)

# Display top 10 results
head(all_results, 10)

# Count significant probes
nrow(significant_results)

# Count probes with higher expression in AML
upregulated_probes <- sum(significant_results$logFC > 1)

# Count probes with higher expression in ALL
downregulated_probes <- sum(significant_results$logFC < -1)

# Display counts
cat("Total significant probes:", nrow(significant_results), "\n")
cat("Higher expression in AML:", upregulated_probes, "\n")
cat("Higher expression in ALL:", downregulated_probes, "\n")

# Load required packages
library(ggplot2)
library(ggrepel)

# Classify probes
all_results$Significance <- "Not Significant"

all_results$Significance[
  all_results$adj.P.Val < 0.05 & all_results$logFC > 1
] <- "Higher in AML"

all_results$Significance[
  all_results$adj.P.Val < 0.05 & all_results$logFC < -1
] <- "Higher in ALL"

# Calculate negative log10 adjusted p-value
all_results$neg_log10_adj_p <- -log10(
  pmax(all_results$adj.P.Val, .Machine$double.xmin)
)

# Select top 10 significant probes for labels
# Select top 10 significant probes for labels
top_labels <- head(
  significant_results[order(significant_results$adj.P.Val), ],
  10
)

# Calculate negative log10 adjusted p-value for labels
top_labels$neg_log10_adj_p <- -log10(
  pmax(top_labels$adj.P.Val, .Machine$double.xmin)
)

# Create volcano plot
volcano_plot <- ggplot(
  all_results,
  aes(x = logFC, y = neg_log10_adj_p, color = Significance)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_labels,
    aes(
      x = logFC,
      y = neg_log10_adj_p,
      label = Probe
    ),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    title = "Volcano Plot: ALL vs AML",
    x = "Log Fold Change (AML - ALL)",
    y = "-Log10 Adjusted P-value",
    color = "Expression"
  ) +
  theme_minimal()

# Display plot
print(volcano_plot)

# Save plot
ggsave(
  "results/figures/volcano_plot.png",
  plot = volcano_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# Load heatmap package
library(pheatmap)

# Select the top 50 most significant probes
top_50_probes <- head(significant_results$Probe, 50)

# Extract expression values
heatmap_data <- expression_data[top_50_probes, ]

# Prepare patient group annotations
annotation_col <- data.frame(
  Group = patient_groups$Group
)

# Match annotation row names with expression sample names
rownames(annotation_col) <- patient_groups$Sample

# Create heatmap
pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Top 50 Differentially Expressed Probes"
)

# Save heatmap
png(
  "results/figures/differential_expression_heatmap.png",
  width = 1200,
  height = 1000,
  res = 150
)

pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Top 50 Differentially Expressed Probes"
)

dev.off()

# Select the top 50 significant probes
top_50_probes <- head(significant_results$Probe, 50)

# Extract expression values
heatmap_data <- expression_data[top_50_probes, ]

# Prepare patient group annotations
annotation_col <- data.frame(
  Group = patient_groups$Group
)

# Match sample IDs
rownames(annotation_col) <- patient_groups$Sample

# Generate and save heatmap
pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Top 50 Differentially Expressed Probes",
  filename = "results/figures/differential_expression_heatmap.png",
  width = 10,
  height = 8
)
