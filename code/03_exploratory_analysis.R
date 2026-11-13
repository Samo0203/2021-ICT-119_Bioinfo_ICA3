# Exploratory Data Analysis (EDA)
# Load saved data


library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

# Load saved datasets
expression_data <- readRDS("data/expression_data.rds")

sample_information <- read.csv(
  "data/sample_information.csv"
)

patient_groups <- read.csv(
  "data/patient_groups.csv"
)

# Check the data
dim(expression_data)
head(sample_information)
table(patient_groups$Group)


# Visualisation : Patient Distribution

patient_distribution <- data.frame(
  Group = names(table(patient_groups$Group)),
  Count = as.integer(table(patient_groups$Group))
)

ggplot(
  patient_distribution,
  aes(x = Group, y = Count, fill = Group)
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = Count),
    vjust = -0.5,
    size = 5
  ) +
  labs(
    title = "Distribution of Leukemia Patient Groups",
    x = "Leukemia Type",
    y = "Number of Patients"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Save the chart
ggsave(
  "results/figures/patient_distribution.png",
  width = 7,
  height = 5,
  dpi = 300
)


# Visualisation : Gene Expression Distribution

# Convert expression matrix into long format

# Convert expression matrix into long format
expression_long <- as.data.frame(
  as.table(expression_data)
)

colnames(expression_long) <- c(
  "Probe",
  "Sample",
  "Expression"
)

# Add patient group information
expression_long$Group <- patient_groups$Group[
  match(
    expression_long$Sample,
    patient_groups$Sample
  )
]

# Calculate the 1st and 99th percentiles
lower_limit <- quantile(
  expression_long$Expression,
  0.01,
  na.rm = TRUE
)

upper_limit <- quantile(
  expression_long$Expression,
  0.99,
  na.rm = TRUE
)

# Create improved boxplot
expression_plot <- ggplot(
  expression_long,
  aes(x = Sample, y = Expression, fill = Group)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  coord_cartesian(
    ylim = c(lower_limit, upper_limit)
  ) +
  labs(
    title = "Gene Expression Distribution Across Samples",
    subtitle = "Displaying the central 98% of expression values",
    x = "Patient Samples",
    y = "Expression Level",
    fill = "Leukemia Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )

# Display chart
print(expression_plot)

# Save chart
ggsave(
  "results/figures/expression_distribution.png",
  plot = expression_plot,
  width = 12,
  height = 6,
  dpi = 300
)

# Visualisation : Principal Component Analysis


# Perform PCA
pca_result <- prcomp(
  t(expression_data),
  center = TRUE,
  scale. = FALSE
)

# Calculate percentage of variance explained
variance_explained <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

# Create PCA data frame
pca_data <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# Add patient groups
pca_data$Group <- patient_groups$Group[
  match(
    pca_data$Sample,
    patient_groups$Sample
  )
]

# Create PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = Group)
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of Leukemia Gene Expression Data",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 2),
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 2),
      "% variance)"
    ),
    color = "Leukemia Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )

# Display PCA plot
print(pca_plot)

# Save PCA plot
ggsave(
  "results/figures/pca_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Display variance explained
variance_explained[1:5]


# Visualisation 4: Gene Expression Heatmap

# Calculate variance for each probe
probe_variance <- apply(
  expression_data,
  1,
  var,
  na.rm = TRUE
)

# Select the 50 most variable probes
top_variable_probes <- names(
  sort(
    probe_variance,
    decreasing = TRUE
  )
)[1:50]

# Extract expression values
heatmap_data <- expression_data[
  top_variable_probes,
]

# Create patient group annotations
annotation_col <- data.frame(
  Group = patient_groups$Group
)

rownames(annotation_col) <- patient_groups$Sample

# Create heatmap
pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_method = "complete",
  main = "Top 50 Variable Gene Expression Probes",
  filename = "results/figures/gene_expression_heatmap.png",
  width = 10,
  height = 8
)

# Confirm selected probes
length(top_variable_probes)
