# Load required libraries
library(golubEsets)
library(Biobase)

# Load the Golub leukemia dataset
data("Golub_Merge")

# Extract gene expression measurements
expression_data <- exprs(Golub_Merge)

# Extract patient information
sample_information <- pData(Golub_Merge)

# Display dataset dimensions
dim(expression_data)

# Display first few rows of patient information
head(sample_information)

# Display available patient information columns
colnames(sample_information)

# Count ALL and AML samples
table(sample_information$ALL.AML)
