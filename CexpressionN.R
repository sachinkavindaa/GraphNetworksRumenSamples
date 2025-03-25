library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(readxl)
library(ggplot2)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(igraph)

data <- read_excel("ORF_matrix_all_heritable-2.xlsx")
data <- as.data.frame(data)
rownames(data) <- data[, 1]  # Replace "1" with the correct column index
data <- data[, -1] 

correlation_matrix <- cor(data, method = "spearman")
head(correlation_matrix[,1:30])
threshold <- 0.5

adjacency_matrix <- ifelse(abs(correlation_matrix) > threshold, 1, 0)
head(adjacency_matrix[,1:20])

diag(adjacency_matrix) <- 0
head(adjacency_matrix[,1:20])

graph <- graph_from_adjacency_matrix(adjacency_matrix , 
                                     mode = "undirected")
graph


#ceb <- cluster_edge_betweenness(graph)

ceb <- cluster_louvain(graph)

community_membership <- membership(ceb)
head(community_membership)
length(ceb)

cluster_sizes <- table(community_membership)
cluster_sizes

igraph::V(graph)$cluster_membership <- community_membership

cluster_colors <- rainbow(max(community_membership))
head(cluster_colors)

V(graph)$color_variable <- cluster_colors[community_membership]

vertex_degrees <- igraph::degree(graph)
head(vertex_degrees)

igraph::V(graph)$degree <- vertex_degrees

vertex.attr <- list(
  cluster_membership = V(graph)$cluster_membership,
  name = V(graph)$name,
  degrees = V(graph)$degree
  #color = V(graph)$color_variable)

graphml_file <- "Test_unweighted_undirected.graphml"

write_graph(
  graph,
  graphml_file,
  format  = "graphml")

vertex_data <- data.frame(vertex.attr)

write.csv(vertex_data, "vertex_data19_unweighted_undirected.csv", row.names = FALSE)

degree_centrality <- degree(graph)
mean(degree_centrality)

average.path.length(graph)

betweenness_centrality <- betweenness(graph)
mean(betweenness_centrality)


save(betweenness_centrality, file = "betweenness_centrality19.RData")

V(graph)$betweenness = betweenness_centrality

closeness_centrality <- igraph::closeness(graph)
V(graph)$closeness = closeness_centrality

save(closeness_centrality, file = "closeness_centrality19.RData")


Eigenvectors  = eigen_centrality(graph)$vector
V(graph)$eigen = Eigenvectors

vertex.attr.ext <- list(
  cluster_membership = V(graph)$cluster_membership,
  name = V(graph)$name,
  degrees = V(graph)$degree,
  color = V(graph)$color_variable, 
  betweeness = V(graph)$betweenness,
  closeness = V(graph)$closeness,
  eigenvectors = V(graph)$eigen)

vertex_data_ext <- data.frame(vertex.attr.ext)
head(vertex_data_ext)

write.csv(vertex_data_ext, "vertex_data_unweighted_undirected19.csv", row.names = FALSE)

clustering_coefficient <- transitivity(graph, type = "global")
clustering_coefficient

ddist <- igraph::degree.distribution(graph)

df <- data.frame(Degree = as.factor((seq_along(ddist)) - 1),
                 Fraction = ddist)

ggplot(data = df, aes(x = Degree, y = Fraction, group = 1)) +
  geom_line() +
  geom_point() +
  theme_bw()

cwt <- cluster_walktrap(graph)
membership2 <- membership(cwt)
table(membership2)

save(membership2, cwt, file = "community_results.RData")

write.csv(table(membership2), file = "community_sizes.csv")

adjacency_matrix_1 <- ifelse(abs(correlation_matrix) > threshold, as.numeric(correlation_matrix), 0)
head(adjacency_matrix_1[,1:40])
diag(adjacency_matrix_1) <- 0

graph1 <- graph_from_adjacency_matrix(adjacency_matrix_1, 
                                      mode = "undirected", weighted = TRUE)
graph1

adjacency_matrix <- ifelse(abs(correlation_matrix) > threshold, as.numeric(abs(correlation_matrix)), 0)
head(adjacency_matrix[,1:40])
diag(adjacency_matrix) <- 0

# Now the weighted option is going to be set to TRUE
graph2 <- graph_from_adjacency_matrix(adjacency_matrix , 
                                      mode = "undirected", weighted = TRUE)
graph2

head(E(graph2)$weight)

vertex_degrees <- igraph::degree(graph2)

cfg <- cluster_fast_greedy(graph2)
membership2 <- membership(cfg)
table(membership2)

community_membership_1 <- membership(cfg)

igraph::V(graph2)$cluster_membership <- community_membership_1

vertex_degrees_1 <- igraph::degree(graph2)

igraph::V(graph2)$degree <- vertex_degrees_1

cluster_colors_1 <- viridis(max(community_membership_1))

V(graph2)$color_variable <- cluster_colors_1[community_membership_1]

E(graph2)$corr_scores = E(graph1)$weight

vertex.attr = list(
  cluster_membership = V(graph2)$cluster_membership,
  name = V(graph2)$name,
  degrees = V(graph2)$degree,
  color = V(graph2)$color_variable,
  size  = V(graph2)$degree)

edge.attr = list(
  corr_scores = round(E(graph2)$corr_scores, 3 ))

graphml_file_1 <- "graph_mRNA_with_metadata_undirected_WEIGHTED.graphml"

write_graph(
  graph2,
  graphml_file_1,
  format  = "graphml")


library(RCy3)

# Start Cytoscape (must be open)
cytoscapePing() 

# Create network from igraph object
createNetworkFromIgraph(
  graph2,
  title = "mRNA Network",
  collection = "My Networks"
)






