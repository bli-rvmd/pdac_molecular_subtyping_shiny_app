library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)


softmax <- function(x) {
  exp_x <- exp(x)               # Exponentiate each value
  probabilities <- exp_x / sum(exp_x) # Normalize by the sum of exponentials
  return(probabilities)
}


#### Script for centroid-based classifications made from clustering-based molecular subtyping
predict_pdac_subtype <- function(df_exp, sig_gene_list, max_perc_missing = 0, convert_na_to_zero = FALSE) {
  
  # samples 
  samples <- colnames(df_exp)
  samples <- samples[samples != "gene_name"]
  
  classes <- names(sig_gene_list) 
  
  df_coefs <- data.frame()
  
  df_coefs <- do.call(cbind, lapply(classes, function(cls) {
    
    df_tmp <- sig_gene_list[[cls]] %>%
      left_join(df_exp, by = "gene_name")
    
    coefs <- sapply(samples, function(s) {
      cor.test(df_tmp[[s]], df_tmp[["value"]], method = "spearman")$estimate
    })
    
    as.data.frame(coefs)
    
  }))
  colnames(df_coefs) <- classes
  
  # convert coefs to probabilities by softmax
  # df_probs <- as.data.frame(t(apply(df_coefs, 1, softmax)))
  # df_probs$Prediction <- apply(df_probs, 1, function(x) names(x)[which.max(x)])
  # df_probs$sample_id <- samples
  
  # add predictive class
  df_coefs$Prediction <- apply(df_coefs, 1, function(x) names(x)[which.max(x)])
  df_coefs$sample_id <- samples
  
  # reorder columns 
  df_coefs <- df_coefs %>%
    select(sample_id, Prediction, dplyr::everything())
  
  return(df_coefs)
  # return(df_probs)
  
}


# function to generate gene data 
generate_pdac_subtype_genedata <- function(df_exp, sig_gene_list) {
  
  df_gene <- as.data.frame(sig_gene_list) %>%
    select(1, ends_with("value")) 
  
  colnames(df_gene)[1] <- "gene_name"
  
  df_gene <- df_gene %>%
    left_join(df_exp, by = "gene_name")
  
  return (df_gene)
  
}


generate_heatmap_with_annotations <- function(expression_data, classification_df, gene_name_column_string, order_by = NULL, gene_limit = NULL, order_genes_by = NULL) {
  # expression_data: data frame of normalized gene expression with genes as rows and samples as columns, including a gene name column
  # classification_df: data frame of multiple classifiers with each column representing a different classification, plus a "sample_id" column to match samples
  # gene_name_column_string: string used to identify the column in expression_data that contains gene names
  # order_by: the name of the classifier column by which the samples (columns) will be ordered
  # gene_limit: upper limit of genes above which gene names will be omitted
  # order_genes_by: a vector of gene names to subset and order the genes to plot
  
  # Ensure classification_df contains a 'sample_id' column
  if (!"sample_id" %in% colnames(classification_df)) {
    stop("The classification_df must contain a 'sample_id' column.")
  }
  
  # Find the column in expression_data that matches the gene_name_column_string
  gene_name_col_index <- grep(gene_name_column_string, colnames(expression_data), ignore.case = TRUE)
  
  # Check if we found exactly one gene name column
  if (length(gene_name_col_index) != 1) {
    stop("Could not find exactly one column matching the gene_name_column_string.")
  }
  
  
  # Extract gene names and remove the gene name column from expression_data
  gene_names <- expression_data[, gene_name_col_index]
  
  # 
  # print(gene_name_col_index)
  # print(length(gene_names))
  # print(length(unique(gene_names)))
  # 
  # print(gene_names[duplicated(gene_names)])
  # 
   
  expression_data <- expression_data[, -gene_name_col_index]
  rownames(expression_data) <- gene_names  # Set gene names as rownames
  
  # Match samples in expression_data to classification_df by the sample_id column
  if (!all(colnames(expression_data) %in% classification_df$sample_id)) {
    stop("Not all sample IDs in expression_data are present in classification_df$sample_id")
  }
  
  
  
  # Reorder classification_df to match the order of sample IDs in expression_data
  classification_df <- classification_df[match(colnames(expression_data), classification_df$sample_id), ]
  
  # If an order_by classifier is provided, reorder the columns by this classifier
  if (!is.null(order_by)) {
    if (!(order_by %in% colnames(classification_df))) {
      stop("The specified classifier name does not exist in classification_df")
    }
    ordered_samples <- order(classification_df[[order_by]])
    expression_data <- expression_data[, ordered_samples]
    classification_df <- classification_df[ordered_samples, ]
  }
  
  # Subset and reorder the genes based on order_genes_by if provided
  if (!is.null(order_genes_by)) {
    if (!all(order_genes_by %in% rownames(expression_data))) {
      stop("Some genes in order_genes_by are not found in the expression_data.")
    }
    expression_data <- expression_data[order_genes_by, , drop = FALSE]  # Subset and order the genes
  }
  
  # Determine whether to show gene names based on gene_limit
  show_gene_names <- ifelse(!is.null(gene_limit) && nrow(expression_data) > gene_limit, FALSE, TRUE)
  
  # Create a top annotation for each classifier using ComplexHeatmap
  top_annotation <- HeatmapAnnotation(df = classification_df[, -which(colnames(classification_df) == "sample_id")])  # Exclude sample_id from annotations
  
  # Set the color palette for the heatmap
  colors <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Generate the heatmap with gene names annotated based on the gene_limit and multiple classification bars on top
  p <- Heatmap(
    expression_data,
    name = "Expression",             # Name of the heatmap legend
    col = colors,                    # Define color scale
    top_annotation = top_annotation, # Add multiple annotations (classifiers) on top
    show_row_names = show_gene_names, # Show or omit gene names based on gene_limit
    cluster_columns = is.null(order_by),  # Cluster samples unless ordered by classifier
    cluster_rows = is.null(order_genes_by)  # Cluster genes unless order_genes_by is provided
  )
  
  return(p)
}


######
# R script for Single Sample kTSP Binary Classifier

Single_Sample_kTSP_Binary_Classifier <- function(params) {
  
  # Initialize the classifier with parameters
  gene_pairs <- params$gene_pairs
  coefs <- params$coefs
  intercept <- params$intercept
  classifier_name <- params$classifier_name
  
  # Function for logistic regression classification
  logreg <- function(X, output_type = "class") {
    
    # X should be a named vector (gene expression values) where the names are gene symbols (Hugo Symbol format)
    xs <- numeric(length(gene_pairs))
    
    for (i in seq_along(gene_pairs)) {
      g <- gene_pairs[[i]]
      gene_1 <- g[1]
      gene_2 <- g[2]
      
      # Check if both genes exist in the expression data
      if (!(gene_1 %in% names(X)) || !(gene_2 %in% names(X))) {
        xs[i] <- 0.5  # If either gene is missing, assign 0.5
      } else if (X[gene_1] > X[gene_2]) {
        xs[i] <- 1  # Gene 1 has higher expression
      } else if (X[gene_1] < X[gene_2]) {
        xs[i] <- 0  # Gene 2 has higher expression
      } else {
        xs[i] <- 0.5  # If expressions are equal, assign 0.5
      }
    }
    
    # Perform classification by summing weighted values and adding intercept
    classification <- sum(xs * coefs) + intercept
    
    if (output_type == "class") {
      if (classification > 0) {
        return("Basal")
      } else {
        return("Classical")
      }
    } else if (output_type == "score") {
      return(classification)
    }
  }
  
  # Return a list of functions to simulate class methods
  list(
    logreg = logreg,
    classifier_name = classifier_name
  )
}


# Moffitt parameters
moffitt_params <- list(
  classifier_name = "Moffitt",
  gene_pairs = list(
    c("CD109", "GPR160"),
    c("SLC2A1", "AGR2"),
    c("KRT16", "SLC44A4"),
    c("CTSV", "TMEM45B"),
    c("KRT6A", "BCAS1"),
    c("B3GNT5", "VSIG2"),
    c("MET", "TFF3"),
    c("CHST6", "PLA2G10"),
    c("SERPINB5", "HPGD"),
    c("DCBLD2", "PLS1"),
    c("IL20RB", "FAM3D"),
    c("PPP1R14C", "SYTL2"),
    c("NAB1", "PLEKHA6"),
    c("MSLN", "CAPN9")
  ),
  coefs = c(0.87, 1.22, 0.52, 1.43, 0.70, 0.41, 0.72, 0.80, 0.76, 1.40, 1.33, 1.58, 0.41, 1.58),
  intercept = -7.16
)

# PurIST parameters
purist_params <- list(
  classifier_name = "PurIST",
  gene_pairs = list(
    c("GPR87", "REG4"),
    c("KRT6A", "ANXA10"),
    c("BCAR3", "GATA6"),
    c("PTGES", "CLDN18"),
    c("ITGA3", "LGALS4"),
    c("C16orf74", "DDC"),
    c("S100A2", "SLC40A1"),
    c("KRT5", "CLRN3")
  ),
  coefs = c(1.994, 2.031, 1.618, 0.922, 1.059, 0.929, 2.505, 0.485),
  intercept = -6.815
)


# Initialize the classifier
moffitt_classifier <- Single_Sample_kTSP_Binary_Classifier(moffitt_params)

purist_classifier <- Single_Sample_kTSP_Binary_Classifier(purist_params)


predict_pdac_subtype_ss <- function(df_exp, method = "Moffitt", max_perc_missing = 0, convert_na_to_zero = FALSE) {
  
  # samples 
  samples <- colnames(df_exp)
  samples <- samples[samples != "gene_name"]
  
  idx_gene_name <- which(colnames(df_exp) == "gene_name")
  
  dat <- t(df_exp[, -idx_gene_name])
  colnames(dat) <- df_exp$gene_name
  
  if (method == "Moffitt") {
    
    return (data.frame(
      sample_id = samples, 
      Prediction = sapply(1:nrow(dat), function(idx) {moffitt_classifier$logreg(dat[idx, ], output_type = "class")}), 
      Score = sapply(1:nrow(dat), function(idx) {moffitt_classifier$logreg(dat[idx, ], output_type = "score")})
    ))
  }
  
  if (method == "PurIST") {
    
    return (data.frame(
      sample_id = samples, 
      Prediction = sapply(1:nrow(dat), function(idx) {purist_classifier$logreg(dat[idx, ], output_type = "class")}), 
      Score = sapply(1:nrow(dat), function(idx) {purist_classifier$logreg(dat[idx, ], output_type = "score")})
    ))
    
  }
  
}


generate_pdac_subtype_genedata_ss <- function(df_exp, method = "Moffitt") {
  
  if (method == "Moffitt") {
    params <- moffitt_params
  } else if (method == "PurIST") {
    params <- purist_params  
  } else {
    stop ("Check input option of 'method' of function 'generate_pdac_subtype_genedata_ss!")
  }
  
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


## function to calculate median expression from input data by each subtype
calculate_median_exp <- function(df, df_class, gene_id_col_name = "symbol") {
  
  df_t <- df %>%
    column_to_rownames(var = gene_id_col_name) %>%
    t() %>%
    as.data.frame()
  
  df_res <- df_t %>%
    rownames_to_column(var = "sample_id") %>%
    left_join(df_class %>% select(sample_id, "class_name"), by = "sample_id") %>%
    column_to_rownames(var = "sample_id") %>%
    group_by(class_name) %>%
    summarize(dplyr::across(dplyr::everything(), median, na.rm = TRUE))
  
  df_output <- df_res %>%
    column_to_rownames(var = "class_name") %>%
    t() %>%
    as.data.frame()
  
  return (df_output)
  
}


## function to calculate average expression of each subgroup
calculate_mean_exp <- function(df, df_class, gene_id_col_name = "symbol") {
  
  df_t <- df %>%
    column_to_rownames(var = gene_id_col_name) %>%
    t() %>%
    as.data.frame()
  
  df_res <- df_t %>%
    rownames_to_column(var = "sample_id") %>%
    left_join(df_class %>% select(sample_id, "class_name"), by = "sample_id") %>%
    column_to_rownames(var = "sample_id") %>%
    group_by(class_name) %>%
    summarize(dplyr::across(dplyr::everything(), mean, na.rm = TRUE))
  
  df_output <- df_res %>%
    column_to_rownames(var = "class_name") %>%
    t() %>%
    as.data.frame()
  
  return (df_output)
  
}