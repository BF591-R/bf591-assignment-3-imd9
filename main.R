library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
  data <- read.csv(intensity_data, sep = delimiter)
  return(data)
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  v <- pca_results$sdev^2
  ve <- v / sum(v)
  return(ve)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PCA names, variance explained by each PC, and the
#' cumulative sum of variance explained
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#'
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained and the cumulative variance explained
#' @export
#'
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  var_tib <- tibble(variance_explained = pca_ve) %>% 
  mutate(principal_components = factor(colnames(pca_results$x), levels=colnames(pca_results$x)), cumulative = cumsum(variance_explained))
  return(var_tib)
}

#' Define a function that creates a bar plot of the variance explained by each
#' PC along with a scatter plot showing the cumulative sum of variance explained
#' using ggplot2
#' @param variance_tibble (tibble): the tibble gnerated in the previous function
#' that contains each PC label, the variance explained by each PC, and the 
#' cumulative sum of variance explained
#'
#' @return A ggplot with a barchart representing individual variance
#'   explained and a scatterplot (connected with a line) that represents the
#'   cumulative sum of PCs
#' @export
#'
#' @examples
plot_pca_variance <- function(variance_tibble) {
  ls_cols <- c('Cumulative'='black')
  bar_cols <-c('Variance Explained'='light blue')
  variance_tibble %>%
    ggplot() + 
    geom_bar(aes(x=principal_components, y=variance_explained, fill='Variance Explained'), stat='identity', color='black') + 
    geom_line(aes(x=principal_components, y=cumulative, group=1, color='Cumulative')) +
    geom_point(aes(x=principal_components, y=cumulative, group=1, color='Cumulative')) +
    scale_colour_manual(name="Cumulative",values=ls_cols) +
    scale_fill_manual(name="Variance Explained",values=bar_cols) +
    labs(x='PC', y='% variance') +
    theme_classic(base_size=8) + 
    theme(axis.text.x=element_text(angle=90,hjust=1))
}

#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
  meta <- readr::read_csv(metadata) %>% 
    select(geo_accession, SixSubtypesClassification)

  labeled <- pca_results$x %>% 
    as_tibble(rownames='geo_accession') %>% 
    left_join(meta, by='geo_accession')
  
  biplot <- labeled %>% 
    ggplot() + 
    geom_point(aes(x=PC1, y=PC2, color=SixSubtypesClassification)) +
    theme_classic()
  return(biplot)
}


#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_csv (str): The path to the differential expression results
#'   file we have provided
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the CSV.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_csv, fdr_threshold) {
sig_ids <- read_data(diff_exp_csv, ',') %>% 
  as_tibble(rownames='probeid') %>% 
  filter(padj < fdr_threshold) %>% 
  pull(probeid)
  return(sig_ids)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
  de_intensity <- intensity %>% 
    as_tibble(rownames='probeid') %>% 
    filter(probeid %in% sig_ids_list) %>% 
    column_to_rownames(var='probeid') %>% 
    as.matrix()
  return(de_intensity)
}

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
  col.pal <- RColorBrewer::brewer.pal(num_colors, palette)
  return(heatmap(de_intensity, col=col.pal))
}