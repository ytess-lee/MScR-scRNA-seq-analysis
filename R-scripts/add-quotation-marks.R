# Load required library
library(readr)

# Specify the path to the CSV file
file_path <- '/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/Book1.csv'


# Function to format the gene names
format_gene_names <- function(file_path) {
  # Read the CSV file
  data <- read_csv(file_path, col_names = FALSE)
  
  # Extract the first column (assuming gene names are in the first column)
  gene_names <- data[[1]]
  
  # Trim any leading/trailing whitespace and remove empty lines
  gene_names <- trimws(gene_names)
  gene_names <- gene_names[gene_names != ""]
  
  # Add quotes around each gene name
  formatted_names <- paste0("\"", gene_names, "\"")
  
  # Concatenate the names with commas, except for the last one
  result <- paste(c(formatted_names[-length(formatted_names)], formatted_names[length(formatted_names)]), collapse = ", ")
  
  return(result)
}


# Format the gene names
formatted_gene_names <- format_gene_names(file_path)

# Print the result
cat("Formatted gene names:\n")
cat(formatted_gene_names)

