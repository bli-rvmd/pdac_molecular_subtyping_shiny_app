### FIXME!
## 1. add actual scores of classifiers (done)
## 2. add missingness of genes as a table to each sample of each classifier (done)
## 3. add "Download" button (full results including predictions of each classifier and missingness of genes) (done)
## 4. add feature of allowing gene missingness to calculate score | or having to have all genes available to calculate score 
## 4.1 if allowing to calculate score despite missingness (either deem missingness as NA, by default, such that for Moffitt/Purist a 0.5 coefficient would be assigned and for other methods cor.test(...) would remove genes of missingness, or convert genes of missingness to 0)
## 4.2 on frontend, add two variables - one controlling max_percentage_missingness_allowed, and the other if_or_not_converting_NA_to_0

## 5. fix app 'plot' function when there's gene missingness (Error: Some genes in order_genes_by are not found in the expression_data.)


# Load necessary libraries
library(shiny)
library(DT) # for data table display
library(dplyr)
library(openxlsx)
library(readxl)
library(tools)
library(bslib)
library(jsonlite)

list_top_genes <- fromJSON("/Users/bli/Documents/Projects/WebApps/pdac_molecular_subtyping_shiny_app/list_top_genes.json")

source("/Users/bli/Documents/Projects/WebApps/pdac_molecular_subtyping_shiny_app/utility_functions.R")

## test predict_pdac_subtype function 
df_exp <- read_excel("/Users/bli/Documents/Projects/PDAC_classification/data/bailey_dataset_normalized_expression_grch37_annotated_genes.xlsx")
colnames(df_exp)[1] <- "gene_name"


dat_tmp <- predict_pdac_subtype(df_exp, list_top_genes[["Bailey"]])


############## TEST BELOW #############

dat_tmp_ss <- predict_pdac_subtype_ss(df_exp, "Moffitt")


#### FIXME! show list of gene names and bifurcate 

tmp <- generate_pdac_subtype_genedata(df_exp, list_top_genes$Collisson)

#######


generate_pdac_subtype_genedata_ss <- function(df_exp, params) {
  
  gene_pair_df <- do.call(rbind, lapply(seq_along(params$gene_pairs), function(i) {
    pair <- params$gene_pairs[[i]]
    coef <- params$coefs[i]
    
    # Dynamically identify sample columns
    sample_cols <- setdiff(names(df_exp), "gene_name")
    
    # Prepare data for the two genes
    pair_data <- df_exp %>%
      filter(gene_name %in% pair) %>%
      complete(gene_name = factor(pair, levels = pair), fill = setNames(as.list(rep(NA, length(sample_cols))), sample_cols))
    
    # Combine data into one row
    if (nrow(pair_data) == 2) {
      combined_samples <- pair_data %>%
        select(-gene_name) %>%
        summarise(across(everything(), ~ paste(., collapse = "|"))) %>%
        as.list()
      
      data.frame(
        gene_pair = paste(pair, collapse = "|"),
        coef = coef,
        combined_samples
      )
    } else {
      # Handle missing data for all sample columns
      missing_samples <- setNames(as.list(rep("NA|NA", length(sample_cols))), sample_cols)
      
      data.frame(
        gene_pair = paste(pair, collapse = "|"),
        coef = coef,
        missing_samples
      )
    }
  })) %>%
    bind_rows()
  
  return (gene_pair_df)
  
}
