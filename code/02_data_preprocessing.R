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


#  Data Preprocessing

# Check missing values in expression data
sum(is.na(expression_data))

# Check infinite values
sum(is.infinite(expression_data))

# Check missing values in patient information
colSums(is.na(sample_information))

# Create patient group information
patient_group <- factor(
  sample_information$ALL.AML,
  levels = c("ALL", "AML")
)

# Check patient groups
table(patient_group)

# Check expression data structure
str(expression_data)

# Save expression data
saveRDS(
  expression_data,
  "data/expression_data.rds"
)

# Save patient information
write.csv(
  sample_information,
  "data/sample_information.csv",
  row.names = FALSE
)

# Save patient groups
write.csv(
  data.frame(
    Sample = colnames(expression_data),
    Group = patient_group
  ),
  "data/patient_groups.csv",
  row.names = FALSE
)

# Confirm saved files
list.files("data")
