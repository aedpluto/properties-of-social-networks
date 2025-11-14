library(igraph)
library(blockmodeling)

# import data
data <- data.frame(read.csv("edgelist.csv", header=TRUE))
redactions <- read.csv("individualInfo.csv", header=TRUE)

# redact data
for (n in redactions$name) {
  data$nameIn[data$nameIn == n] <- redactions$num[redactions$name == n]
  data$nameOut[data$nameOut == n] <- redactions$num[redactions$name == n]
}

# convert data to a matrix
M <- as.matrix(data)
g <- graph_from_edgelist(M)

# plot data
set.seed(2)
png("my-plot.png", width=1200, height=1200)
par(mar=c(1,1,3,1))
plot(g, vertex.size=4, edge.arrowsize=0.1, edge.arrow.width=0.5, vertex.label.cex=3, edge.curved=0.1, asp=0, vertex.label.dist=0.5) # layout=layout.circle
title("Social group", cex.main=4)
dev.off()

# degree and diameter
cat("Average degree:", mean(degree(g, mode="in"))) # be careful since directed
cat("Edge density:", edge_density(g, loops=FALSE))
cat("Diameter:", max(shortest.paths(g, mode="out")))

# calculate and print the centrality scores

dC <- sort(degree(g), decreasing = TRUE)[1:5] # mode="in"
cC <- sort(closeness(g), decreasing = TRUE)[1:5]
bC <- sort(betweenness(g), decreasing = TRUE)[1:5]
prC <- sort(page.rank(g)$vector, decreasing = TRUE)[1:5]

degreeC <- V(g)$name[degree(g) %in% dC]
closenessC <- V(g)$name[closeness(g) %in% cC]
betweennessC <- V(g)$name[betweenness(g) %in% bC]
pageRankC <- V(g)$name[page.rank(g)$vector %in% prC]

cat("Degree central nodes:", degreeC)
cat("Closeness central nodes:", closenessC)
cat("Betweenness central nodes:", betweennessC)
cat("PageRank central nodes:", pageRankC)

# select the most central nodes
toselect <- c()

for (m in c(degreeC, closenessC, betweennessC, pageRankC))
  for (v in m)
    toselect <- c(toselect, v)

print(toselect)

# graph the most central nodes
centralityScores <- cbind(degreeC, closenessC, betweennessC, pageRankC)
rownames(centralityScores) <- c("1", "2", "3", "4", "5")
colnames(centralityScores) <- c("Degree", "Closeness", "Betweenness", "Page Rank")
print(centralityScores)

set.seed(2)
png("centralityScores.png", width=1200, height=1200)
par(mfrow=c(2,2), mar=c(5,10,12,5))
barplot(height=dC, names.arg=degreeC, xlab="Individual", ylab="Score", main="Degree centrality", col="red", cex.main=2, cex.names=2, cex.axis=2, cex.lab=2)
barplot(height=cC, names.arg=closenessC, xlab="Individual", ylab="Score", main="Closeness centrality", col="orange", cex.main=2, cex.names=2, cex.axis=2, cex.lab=2)
barplot(height=bC, names.arg=betweennessC, xlab="Individual", ylab="Score", main="Betweenness centrality", col="darkgreen", cex.main=2, cex.names=2, cex.axis=2, cex.lab=2)
barplot(height=prC, names.arg=pageRankC, xlab="Individual", ylab="Score", main="PageRank centrality", col="navy", cex.main=2, cex.names=2, cex.axis=2, cex.lab=2)
dev.off()

# cliques
largeClq <- largest.cliques(g) # ignoring edge direction
otherClq <- cliques(g, min=8, max=9)

# calculate cliques where edges are reciprocated
reciprocatedEdges <- which(is.mutual(g))
g1 <- subgraph.edges(g, reciprocatedEdges, delete.vertices = TRUE)
g1_undirected <- as.undirected(g1, mode = "mutual")
clqs <- cliques(g1_undirected, min=6)
print(largest.cliques(g1_undirected))
c1 <- c("1", "2", "3", "4", "12", "14", "22")
c2 <- c("1", "2", "3", "4", "12", "15", "22")

# plot the cliques
 png("largest-clique.png", width=1200, height=1200)
 par(mar=c(1,1,3,1), mfrow=c(1,2))
 set.seed(2)
 V(g)$color <- ifelse(V(g)$name %in% c1, "blue", "orange")
 plot(g, asp=0, vertex.size=7, edge.arrowsize=0.1, edge.arrow.width=0.5, vertex.label.cex=3, edge.curved=0.1, asp=0, vertex.label.dist=0.5)
 title("Clique 1", cex.main=4)
 set.seed(2)
 V(g)$color <- ifelse(V(g)$name %in% c2, "darkgreen", "orange")
 plot(g, asp=0, vertex.size=7, edge.arrowsize=0.1, edge.arrow.width=0.5, vertex.label.cex=3, edge.curved=0.1, asp=0, vertex.label.dist=0.5)
 title("Clique 2", cex.main=4)
 dev.off()

# k-cores
numCores <- coreness(g)
table(numCores)

set.seed(2)
png("k-cores.png", width=1200, height=1200)
par(mfrow=c(1,1), mar=c(1,1,3,1))
plot(g, vertex.color=rainbow(6)[unique(as.factor(coreness(g)))], asp=0, vertex.size=4, edge.arrowsize=0.1, edge.arrow.width=0.5, vertex.label.cex=3, edge.curved=0.1, vertex.label.dist=0.5)
legend("bottomleft",bty = "n",legend=unique(as.factor(coreness(g))), cex = 3, fill = rainbow(6)[unique(as.factor(coreness(g)))])
title("Coreness of the social network", cex.main=4)
dev.off()

# block modelling
Gclass2 <- optRandomParC(M=as.matrix(get.adjacency(g)), k=2, rep=10, approach="ss", blocks = "com")
Gclass3 <- optRandomParC(M=as.matrix(get.adjacency(g)), k=3, rep=10, approach="ss", blocks = "com")
Gclass4 <- optRandomParC(M=as.matrix(get.adjacency(g)), k=4, rep=10, approach="ss", blocks = "com")
Gclass5 <- optRandomParC(M=as.matrix(get.adjacency(g)), k=5, rep=10, approach="ss", blocks = "com")

set.seed(2)
png("block-modelling-1.png", width=2500, height=500)
par(mfrow=c(1,4))
plot(Gclass2, main="")
title("Two Block Partition", cex.main=4)
plot(Gclass3, main="") 
title("Three Block Partition", cex.main=4)
plot(Gclass4, main="") 
title("Four Block Partition", cex.main=4)
plot(Gclass5, main="") 
title("Five Block Partition", cex.main=4)
dev.off()

# plot structural equivalence
set.seed(2)
png("block-modelling-2.png", width=2500, height=800)
par(mfrow=c(1,4))

set.seed(2)
plot(g, vertex.color=Gclass2$best$best1$clu, asp=0, vertex.size=8, edge.arrowsize=0.1, edge.arrow.width=0.5, vertex.label.cex=3, vertex.label.dist=0.5) # edge.curved=0.1,
title("Two Block Partition", cex.main=4)

set.seed(2)
plot(g, vertex.color=Gclass3$best$best1$clu, asp=0, vertex.size=8, edge.arrowsize=0.1, edge.arrow.width=0.5, vertex.label.cex=3, vertex.label.dist=0.5)
title("Three Block Partition", cex.main=4)

set.seed(2)
plot(g, vertex.color=Gclass4$best$best1$clu, asp=0, vertex.size=8, edge.arrowsize=0.1, edge.arrow.width=0.5, vertex.label.cex=3, vertex.label.dist=0.5)
title("Four Block Partition", cex.main=4)

set.seed(2)
plot(g, vertex.color=Gclass5$best$best1$clu, asp=0, vertex.size=8, edge.arrowsize=0.1, edge.arrow.width=0.5, vertex.label.cex=3, vertex.label.dist=0.5)
title("Five Block Partition", cex.main=4)
dev.off()

# transitivity
cat("Global transitivity:", transitivity(g, type="global"))
cat("Local transitivity:", transitivity(g, type="local", vids=V(g)))

# reciprocity
cat("Reciprocity:", reciprocity(g))

