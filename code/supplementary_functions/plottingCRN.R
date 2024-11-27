
library(igraph)

set.seed(1)

wd = "/home/raucello/EpiCell_CDifficile"

# setting GEMM tags
model.type = "E_coli"
model.name = "e_coli_core"

# setting results directory
result.dir = paste0(wd, "/input/", model.type, "/", model.name, "/global_sens_analysis")

model.file = paste0(wd, "/input/", model.type, "/", model.name, "/", model.name, ".RData")
all.react.file = paste0(wd, "/input/", model.type, "/", model.name, "/all.react.rds")

load(model.file)
all.react = readRDS(all.react.file)

S = model@S

colnames(S) = model@react_id
rownames(S) = model@met_id

rownames(S) = gsub("\\(", replacement = "_", rownames(S))
rownames(S) = gsub("\\[", replacement = "_", rownames(S))
rownames(S) = gsub("\\)", replacement = "", rownames(S))
rownames(S) = gsub("\\]", replacement = "", rownames(S))

# create graph from incidence matrix
G = graph_from_incidence_matrix(S)

# Which vertices are not isolates?
which_ids <- which(degree(G, mode = "total") > 0)

# Getting the subgraph
G_sub = induced_subgraph(G, which_ids)

plot(G_sub,
     vertex.size = degree(G_sub)/50 + 1,
     edge.arrow.size = 0.25,
     vertex.label.color = "#2a475e",
     vertex.label.family = "Helvetica",
     vertex.label.cex = 0.5,
     vertex.label.dist = 0, # Distance between the label and the vertex
     vertex.label.degree = 0,
     edge.arrow.size =.4,
     edge.label.color = "darkred",
     layout= 2)
