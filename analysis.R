#!/usr/bin/env Rscript --vanilla

# General data handling
library(tidyverse)
# Rarefy + handle microbiome data
library(phyloseq)
# CLR transformation
library(compositions)
# Parallel analysis
library(furrr)
# 3D matrix helper function
library(abind)
# SPRING
#install.packages("devtools")
#devtools::install_github("GraceYoon/SPRING")
#install.packages("pak")
pak::pak("https://github.com/GraceYoon/SPRING")
library(SPRING)

read_data <- function(countfile, taxfile, metadatafile){
  counts <- read_csv(countfile, show_col_types = FALSE) |>
    phyloseq::otu_table(taxa_are_rows = TRUE)
  taxa <- readr::read_delim(taxfile, delim = ";", show_col_types = FALSE) |>
    as.matrix() |>
    phyloseq::tax_table()
  metadata <- metadatafile |>
    read_csv(show_col_types = FALSE, col_types = cols(DAYOFYEAR = col_character())) |>
    tibble::column_to_rownames(var = "SEQUENCEFILE") |>
    phyloseq::sample_data()
  return(phyloseq(counts, taxa, metadata))
}

# Read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 4){
  stop("Usage: analysis.R <countfile> <taxfile> <metadatafile> <outfile>")
}

data <- read_data(
  countfile = args[[1]],
  taxfile = args[[2]],
  metadatafile =  args[[3]]
) |>
  # Remove site 6
  subset_samples(SITE != "6") |>
  # Agglomerate taxa at the genus level
  tax_glom(taxrank = "genus")

# Remove low abundance OTUs
data <- prune_taxa(taxa_sums(data) / sum(taxa_sums(data)) > 0.01, data)

# What I propose is to construct n networks
# by using n different random subsamples of the data
# using rarefaction
sample_size <- min(sample_sums(data))
# SPRING: this is the approach I think it makes more sense
rarefied_network_spring <- function(seed){
  rarefied <- data |>
    # This function is crazy slow, I wonder why...
    # In my head, it should not be that expensive
    rarefy_even_depth(
      sample_size, rngseed = seed,
      replace = TRUE, trimOTUs = FALSE, verbose = FALSE
    )
  fit.spring <- otu_table(rarefied) |>
    t() |>
    as.matrix() |>
    SPRING(
      Rmethod = "approx", quantitative = TRUE,
      verbose = FALSE, seed = seed,ncores = 5
    )
  opt.index <- fit.spring$output$stars$opt.index
  list(
    beta = fit.spring$fit$est$beta[[opt.index]],
    adj = fit.spring$fit$refit$stars
  )
}

plan(multicore, workers=5)
seeds <- 1:100
nets <- future_map(
  seeds, rarefied_network_spring,.progress = TRUE, .options = furrr_options(seed = TRUE)
)

into3D <- \(mat_list){
  mat_list |> map(as.matrix) |> abind(along = 0)
}

adj_matrix <- nets |> map("adj") |> into3D()
beta_matrix <- nets |> map("beta") |> into3D()

list(
  adj_matrix = adj_matrix,
  beta_matrix = beta_matrix,
  physeq = data
) |>
  write_rds(args[[4]])
