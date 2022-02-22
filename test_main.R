#!/usr/bin/Rscript
source("main.R")
library(testthat)

#generate consistent sample data for all tests using it
fx_intensity_df <- data.frame(GSM971958 = runif(1000, 0, 15),
                              GSM971959 = runif(1000, 0, 15),
                              GSM971967 = runif(1000, 0, 15),
                              GSM971968 = runif(1000, 0, 15))

row.names(fx_intensity_df) <- paste0(1:1000, "_s_at")
fx_intensity_mat <- as.matrix(fx_intensity_df)


test_that("read_data loads correctly", {
  #checks dimensions and is.data.frame()
  
  res <- read_data("example_intensity_data.csv", " ")
  
  expect_equal(dim(res), c(54675, 35))
  expect_true(is.data.frame(res))
  
})

test_that("calculate_variance squares each element in a vector and divides by sum", {
  #mimic the object by using a named list so input arg is the same
  
  test_obj <- list(sdev = c(sort(runif(20, min=0, max=100), decreasing=TRUE)))
  test_res <- test_obj$sdev^2 / sum(test_obj$sdev^2)
  
  expect_true(all(dplyr::near(calculate_variance_explained(test_obj), test_res, tol=.1)))
})

test_that("make_variance_tibble", {
  #just checks if length of tibble is the same as number of PCs, and is_tibble()
  
  test_pca_res <- prcomp(scale(t(fx_intensity_mat)), scale=FALSE, center=FALSE)
  test_ve <- test_pca_res$sdev^2 / sum(test_pca_res$sdev^2)
  
  function_tib <- make_variance_tibble(test_ve, test_pca_res)
  
  expect_equal(dim(function_tib), c(length(test_ve), 3))
  expect_true(is_tibble(function_tib))
})

test_that("plot_pca_variance creates a ggplot object with the correct geoms", {
  pcs <- factor(paste0("PC", 1:35), levels = paste0("PC", 1:35), ordered = T)
  tibble(variance_explained = runif(35, 1.0e-30, 1.5e-01),
         principal_components = pcs) %>%
    mutate(cumulative = cumsum(variance_explained)) -> variance_tibble
  plot <- plot_pca_variance(variance_tibble)
  
  # testing all geoms in the ggplot object
  geoms <- c()
  for (geom in plot$layers) {
    geoms <- c(geoms, class(geom$geom)[1])
  }
  expect_true("GeomBar" %in% geoms)
  expect_true("GeomLine" %in% geoms)
  expect_true("GeomPoint" %in% geoms)
})

test_that("list_significant_probes returns the probeids with padj < threshold", {
  #test.csv has two probes < .01, one == .01, and one > .01
  #saved as data file to mimic function args
  
  expect_equal(list_significant_probes('test.csv', .01), c('1_s_at', '2_s_at'))
})

test_that("return_de_intensity returns a matrix and the selected probes", {
  #tests is.matrix() and since function subsets the matrix should be identical
  #no need to worry about rounding, dimensions, etc. ?
  
  test_mat <- fx_intensity_mat[c('1_s_at', '2_s_at'),]
  function_mat <- return_de_intensity(fx_intensity_df, c('1_s_at', '2_s_at'))
  
  expect_true(is.matrix(function_mat))
  expect_true(identical(function_mat, test_mat))
})

test_that("plot_heatmap has correct rows and cols", {
  #this one is an interesting case, since its object (what the function returns)
  #doesn't have as many descriptive elements as ggplot's. We can at least 
  #confirm that their function creates a heatmap and that it has a predictable 
  #number of rows/columns. Since we don't know if their earlier functions will 
  #work, we will create a fake de_intensity matrix to test with.
  
  heatmap <- plot_heatmap(fx_intensity_mat, 11, 'RdBu')
  expect_equal(length(heatmap$colInd), 4)
  expect_equal(length(heatmap$rowInd), 1000)
  expect_equal(class(heatmap), "list")
})
