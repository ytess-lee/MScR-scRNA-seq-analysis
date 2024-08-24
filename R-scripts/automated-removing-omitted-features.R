# this script is to automatically remove omitted features by provide gene list
# it will remove the rows containing the exact gene name and save it to a new csv
# this script REMOVE EXTRACT GENE FOR NO REASON
# UPDATE: this script works, it will spot duplicated genes in the findallarker gene list, that's why it's removing 'more than it should'

library(readr)
library(dplyr)
library(stringr)

# Function to remove rows containing specified gene names
remove_genes_from_csv <- function(input_file, output_file, gene_list) {
  # Read the CSV file into a dataframe
  data <- read_csv(input_file)
  
  # Examine unique gene names in the data
  unique_genes <- unique(data$gene)
  cat("Unique genes in the data:\n")
  print(unique_genes)
  
  # Clean the gene column in the dataframe and the gene list
  data <- data %>%
    mutate(gene = str_trim(gene)) # Remove leading/trailing whitespace
  
  # Clean gene list
  clean_gene_list <- str_trim(gene_list)
  
  # Examine the cleaned unique gene names
  cleaned_unique_genes <- unique(data$gene)
  cat("Cleaned unique genes in the data:\n")
  print(cleaned_unique_genes)
  
  # Remove rows that contain any of the gene names in the cleaned gene list
  filtered_data <- data %>%
    filter(!gene %in% clean_gene_list)
  
  # Write the filtered data to a new CSV file
  write_csv(filtered_data, output_file)
}

# User-defined list of gene names
gene_list <- c(
  "AL450332.1", "VEGFD", "AL136084.3", "MIAT", "GGT5", "AC110611.1", 
  "IL15", "KLF4", "AL161629.1", "PCDHB15", "AC091563.1", "ACSS1", 
  "CDS1", "CAPN8", "RTN4RL1", "Z99289.3", "SERPINB4", "ADAMTSL2", 
  "CCL8", "KLHDC8A", "S1PR5", "AC018647.1", "DCN", "C5orf46", 
  "ACTBL2", "AC245041.2", "AC097480.1", "CRB2", "LRRN4", "LINC01164", 
  "CCN4", "BMP10", "SERPINB7", "EGFL6", "ALPK2", "MAB21L2", "WNT2", 
  "BCL2A1", "AC112721.2", "DLX1", "C1orf115", "SAMD9", "OTOA", 
  "AC068987.3", "AC137932.2", "SLC14A1", "LINC02202", "ABCG2", 
  "CCR5AS", "AC007563.2", "CMKLR1", "SLFN11", "HSPA12B", "DIRC3-AS1", 
  "AC007406.3", "LYPD5", "LAMB3", "STEAP2-AS1", "SPOCK3", "HIST1H3B", 
  "PCYT1B", "AL359513.1", "HOXA-AS2", "TNR", "METTL7A", "AL356010.2", 
  "GLP1R", "AL022328.4", "SLC18A1", "AC091178.2", "LINC02055", 
  "TMOD1", "GDF3", "SNAP91", "ABCG1", "AC012078.2", "AL391261.2", 
  "IRF6", "C2CD4A", "AC004381.1", "MPP4", "CDC25A", "AL021155.5", 
  "TMEM156", "MYEOV", "AL355075.4", "PANTR1", "CBLN1", "EPN3", 
  "VSTM1", "AC239799.2", "Z93241.1", "LRRC38", "TNFSF11"
  )
# File paths
input_file <- "240805-finallmarker-ECs.csv"
output_file <- "o240805-finallmarker-ECs-AUTOMATED.csv"

# Remove specified genes and save the result to a new CSV file
remove_genes_from_csv(input_file, output_file, gene_list)
