#Imports
library(tidyverse)
library(DESeq2)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x 1+m) tibble with a 'gene' column followed by
#' sample names as column names.
#'
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#'
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){
  #reading in all columns as numeric excluding gene (is this neccessary?)
  tibble <- read_delim(filename, delim = "\t", col_types = cols(.default = "n", gene = "c"))
  return(tibble)
}


#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns with sample names as column names.
#'
#' @return tibble: a (n x 1+m) tibble with a 'gene' column followed by m columns
#' of raw counts with genes that have zero variance across samples removed
#'
#' @note (g >= n)
#'
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
  #calculate the variance of each row (excluding the 'gene' column bc it's non-numeric)
  variances <- verse_counts %>%
    select(-gene) %>%
    apply(1, var)
  
  #filter out rows where variance is 0
  filtered_counts <- as_tibble(verse_counts[variances != 0, ])
  
  return(filtered_counts)
}


#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

timepoint_from_sample <- function(x) {
  #strip of "v", any "_", and numbers after the _
  stripped_items <- gsub("^[v]|_\\d+$", "", x) #regex
  return(stripped_items)
}

#' Grab sample replicate number from sample name
#'
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

sample_replicate <- function(x) {
  #regex for including only numbers after the underscore
  numbers_after_underscore <- gsub(".*_(\\d+)$", "\\1", x)
  return(numbers_after_underscore)
}

#' Generate sample-level metadata from sample names.
#'
#' Will include columns named "sample", "timepoint", and "replicate" that store
#' sample names, sample time points, and sample replicate, respectively.
#'
#'
#' @param sample_names vector: character vector of length (_S_) consisting of sample
#' names from count data.
#'
#' @return tibble: a (_S_ x 3) tibble with column names "sample",
#' "timepoint", and "replicate". "sample"holds sample_names; "timepoint"
#' stores sample time points; and "replicate" stores sample replicate
#'
#' @note _S_ < m
#'
#' @example `meta <- meta_info_from_labels(colnames(count_data)[colnames(count_data)!='gene'])`

meta_info_from_labels <- function(sample_names) {
  # Generate the timepoint and replicate information
  timepoints <- sapply(sample_names, FUN=timepoint_from_sample)
  replicates <- sapply(sample_names, FUN=sample_replicate)
  
  # Create a tibble with the required columns
  meta_data <- tibble(
    sample = sample_names,
    timepoint = timepoints,
    replicate = replicates)
  
  return(meta_data)
}


#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @return tibble of read totals from each sample. A tibble can be `(1 x _S_)` 
#' with sample names as columns names OR `(_S_ x 2)` with columns ("sample", "value")
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
  #find the total read counts for each sample by summing across rows
  total_counts <- colSums(select(count_data, -gene))
  
  # Create a tibble with two columns ("sample", "value") 
  library_size <- tibble(
    sample = names(total_counts),
    value = total_counts
  )
  
  return(as_tibble(library_size))
}


#' Normalize raw count data to counts per million WITH pseudocounts using the
#' following formula:
#'     count / (sample_library_size/10^6)
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m columns of cpm normalized read counts
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(count_data) {
  # find the library size for each sample
  library_size <- get_library_size(count_data)
  
  # find counts per million (CPM) with pseudocounts
  cpm_data <- count_data %>%
    mutate(across(-gene, ~ .x / (sum(.x) / 1e6))) %>%
    mutate(across(where(is.numeric), ~ round(.x + 1e-6, 3)))
  
  return(cpm_data)
}

#' Normalize raw count matrix using DESeq2
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @param meta_data tibble: sample-level information tibble corresponding to the
#' count matrix columns
#'
#' @return tibble: DESeq2 normalized count matrix
#' @export
#'
#' @examples
deseq_normalize <- function(count_data, meta_data) {
  # create count matrix without gene column
  count_matrix <- as.matrix(count_data[, -1])
  
  # creating the deseq object
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                colData = meta_data, 
                                design = ~ 1)
  
  # normalize counts using DESeq2
  dds <- DESeq(dds)
  normalized_counts <- counts(dds, normalized = TRUE)
  
  #reate normalized tibble with the gene column back
  normalized_tibble <- bind_cols(gene = count_data$gene, normalized_counts) %>% 
  as_tibble()
  
  #return normalized tibble
  return(normalized_tibble)
}

#' Perform and plot PCA using processed data.
#'
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs.
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(data, meta, title="") {
  # perform PCA on everything except the gene column
  pca_results <- prcomp(t(data %>% select(-gene)), center = TRUE, scale. = TRUE)
  
  # only grab PC1 and PC2 from the PCA results
  pc1and2 <- pca_results$x[, 1:2]
  
  # convert to dataframe and include sample names
  pcs_df <- as.data.frame(pc1and2) %>% rownames_to_column(var = "sample")
  
  # combine with metadata
  all_data <- left_join(pcs_df, meta, by = 'sample')
  
  # calculate explained variance
  e_v <- pca_results$sdev^2 / sum(pca_results$sdev^2)
  pc1_e_v <- round(e_v[1] * 100, 1) 
  pc2_e_v <- round(e_v[2] * 100, 1)
  
  # plot !
  pca_plot <- ggplot(all_data, aes(x = PC1, y = PC2, color = timepoint)) + 
    geom_point() +
    labs(
      title = title,
      x = paste0("PC 1: ", pc1_e_v, "% variance"),
      y = paste0("PC 2: ", pc2_e_v, "% variance")
    )
  
  return(pca_plot)
}


#' Plot gene count distributions for each sample using boxplots.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#'
#' @example `plot_sample_distributions(data, scale_y_axis=TRUE, title='Raw Count Distributions')`

plot_sample_distributions <- function(data, scale_y_axis=FALSE, title="") {
  # Pivot longer to have the samples and their values
  data_long <- data %>% 
    select(-gene) %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "value")
  
  # Depending on the parameter setting, return the different plots
  if (scale_y_axis == FALSE) {
    plot <- ggplot(data_long, aes(x = sample, y = value, color = sample)) + 
      geom_boxplot() +
      labs(title = title, y = "Counts", x = "Sample")
  } else {
    plot <- ggplot(data_long, aes(x = sample, y = value, color = sample)) + 
      geom_boxplot() +
      labs(title = title, y = "Counts", x = "Sample") + 
      scale_y_log10()
  }
  return(plot)
}


#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#'
#' @example `plot_variance_vs_mean(data, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`

plot_variance_vs_mean <- function(data, scale_y_axis=FALSE, title=""){
  # find the mean and variance for each gene across samples
  gene_stats <- data %>%
  rowwise() %>% #going row by row; genes are in rows
  mutate(
    mean_count = mean(c_across(where(is.numeric)), na.rm = TRUE),  # mean of all numeric columns
    variance = var(c_across(where(is.numeric)), na.rm = TRUE)      # variance of all numeric columns
  ) %>%
  ungroup() #ungroup after performing calculations

# rank genes by mean count; not sure why this is important but vastly changes the shape of the graph
gene_stats <- gene_stats %>%
  mutate(rank = rank(mean_count))

# add a line representing the average mean and variance values
mean_var <- gene_stats %>%
  summarise(avg_mean = mean(mean_count, na.rm = TRUE), avg_var = mean(variance, na.rm = TRUE))

# plot scatter plot
plot <- ggplot(gene_stats, aes(x = rank, y = variance)) +
  geom_point(alpha = 0.5) +  # alpha > transparency
  labs(title = title, x = "Gene Rank by Mean Count", y = "Variance") +
  theme_minimal() + 
  geom_hline(yintercept = mean_var$avg_var, color = "darkred", linetype = "dashed") +
  geom_vline(xintercept = mean(gene_stats$rank), color = "lightblue", linetype = "dashed")

# scale the y-axis to log10 if specified
if (scale_y_axis) { 
  plot <- plot + scale_y_log10()
}

return(plot)
}

