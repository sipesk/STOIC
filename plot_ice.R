library(tidyverse)
library(ggsci)
library(phyloseq)
data <- read_rds("ice_spring.Rds")
physeq_ice_genus <- data$physeq
mat_spring <- data$adj_matrix
nreps <- dim(mat_spring)[1]
# Count how many times each edge appears in each network across
# rarefied datasets
counts <- apply(mat_spring, c(2, 3), sum)
rownames(counts) <- taxa_names(physeq_ice_genus)
colnames(counts) <- taxa_names(physeq_ice_genus)

hist(
  counts[upper.tri(counts, diag = FALSE)] / nreps,
  main = "Presence probability of edges",
  xlab = "Fraction of networks",
  ylab = "Number of edges"
  )

# We select those that appear in at least 50% of the networks
tozero_spring <- which(counts / nreps < 0.5)
above_spring <- ifelse(counts/nreps < 0.5, 0, 1)

# With this approach, 10% of the edges are kept
# i.e. we have, hopefully, make the graph sparse and
# enriched in the most important edges
mean(above_spring[upper.tri(above_spring)])

# Estimated partial correlation coefficient, same as negative precision matrix.
partial_correlations <- 1:nreps |>
  # Get the i-matrix
  map(\(i) data$beta_matrix[i,,]) |>
  # Symmetrize the i-matrix
  map(\(mat) SpiecEasi::symBeta(mat, mode = 'maxabs')) |>
  map(as.matrix) |>
  abind::abind(along = 0)

# Compute median of the partial correlations
median_pcor <- apply(partial_correlations, c(2, 3), median)
# Get only the edges above 0.5 prob
weights <- above_spring*median_pcor

# Load necessary libraries
library(GGally)
library(network)

# Assuming 'weights' is an adjacency matrix
net <- network(weights, directed = FALSE)

# Extract the upper triangular part of the adjacency matrix, excluding the diagonal
edge_weights <- weights[upper.tri(weights)]
edge_weights <- edge_weights[edge_weights != 0]

# Set the edge weights as an attribute of the network
set.edge.attribute(net, "weight", edge_weights)

# Define edge colors based on weight (blue for negative, red for positive)
edge_colors <- ifelse(edge_weights < 0, "blue", "red")
set.edge.attribute(net, "color", edge_colors)

# Set the edge width proportional to the absolute value of the weights
edge_widths <- abs(edge_weights) * 3  # Scaled by a factor for better visualization
set.edge.attribute(net, "width", edge_widths)

# Add edge labels (rounded to 3 decimal places)
set.edge.attribute(net, "label", round(edge_weights, 3))

taxa <- tax_table(data$physeq)
phyla <- gsub(" p__", "", taxa[rownames(weights), "phylum"])
set.vertex.attribute(net, "phylum", as.character(phyla))
# Plot the network with custom node and edge attributes
p1 <- ggnet2(net, 
       edge.color = "color", 
       edge.size = "width", 
       edge.label = "label",
       edge.alpha = 0.8,
       edge.label.size = 4,
       edge.node.size = 1,
       node.color = "phylum",  # Apply the color palette based on phylum
       node.size = 5,                # Adjust node size for better visibility
       label = TRUE,                 # Add node labels
       label.size = 4,
) +
  scale_color_uchicago(name = "Phylum") +
  theme(legend.position = "bottom") +
  ggtitle("Ice data network plot", subtitle = "Edges show median partial correlation coefficient. We show edges above 0.5 probability of absence (i.e. they were present in more than 50% of 100 repeated rarefactions).")

ggsave("network_plot.pdf", p1, width = 11.69, height = 8.27, units = "in", dpi = "retina")
