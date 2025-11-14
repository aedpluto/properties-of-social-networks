### Assignment 2 ###

library(igraph)
library(graphlayouts)
library(ergm)
library(intergraph)
library(network)

# importing data
data <- data.frame(read.csv("edgelist.csv", header=TRUE))
redactions <- read.csv("individualInfo.csv", header=TRUE)

# redacting data
for (n in redactions$name) {
  data$nameIn[data$nameIn == n] <- redactions$num[redactions$name == n]
  data$nameOut[data$nameOut == n] <- redactions$num[redactions$name == n]
}

# converting data to a matrix
M <- as.matrix(data)
g <- graph_from_edgelist(M)

# plotting data
set.seed(2)
png("networkPlotted.png", width=1200, height=1200)
par(mar=c(1,1,3,1))
plot(g, vertex.size=6, edge.arrowsize=1, edge.arrow.width=2, vertex.label.cex=3, edge.curved=0.15, asp=0, vertex.label.dist=1.1)
title("Social group", cex.main=4)
dev.off()


### Assortativity and communities ###

# choosing subgroups: uniA = 1, uniB = 2, other = 3
universities <- c(1,1,1,1,3,3,1,2,1,2,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,3,3,1,3,3)

uniA <- ifelse(universities == 1, 1, 2)
uniB <- ifelse(universities == 2, 1, 2)
other <- ifelse(universities == 3, 1, 2)

# calculate modularity using chosen subgroups
modTotal <- modularity(g, membership = universities)
cat("Modularity of the network:", modTotal)

# assortativity
yearGroup <- c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 5, 4, 3, 5, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2, 4)
assortU <- assortativity_nominal(g, as.factor(universities), directed = TRUE)
assortYG <- assortativity_nominal(g, as.factor(yearGroup), directed = TRUE)
cat("Assortativity according to university:", assortU)
cat("Assortativity according to year group:", assortYG)

# community detection algorithm
communityDetection1 <- cluster_infomap(g)
print(communityDetection1)

communityDetection2 <- cluster_walktrap(g)
print(communityDetection2)

set.seed(2)
png("communityDetectionInfoMap.png", width=1200, height=1200)
par(mar=c(1,1,3,1))
plot(g, vertex.size=6, edge.arrowsize=1, edge.arrow.width=2, vertex.label.cex=3, edge.curved=0.15, asp=0, vertex.label.dist=1.1, vertex.color=membership(communityDetection1))
title("Community detection (InfoMap)", cex.main=4)
dev.off()

set.seed(2)
png("communityDetectionWalkTrap.png", width=1200, height=1200)
par(mar=c(1,1,3,1))
plot(g, vertex.size=6, edge.arrowsize=1, edge.arrow.width=2, vertex.label.cex=3, edge.curved=0.15, asp=0, vertex.label.dist=1.1, vertex.color=membership(communityDetection2))
title("Community detection (WalkTrap)", cex.main=4)
dev.off()


### Small-world and scale-free networks ###

# calculate the degree distribution
deg_table <- table(degree(g, mode="all"))
deg_dist <- as.numeric(deg_table) / sum(deg_table)
degree_distribution(g)

png("degreeDistribution.png", width=700, height=700)
par(mfrow=c(2,2))
title("Measures of degree distribution", cex.main=4)

hist(degree(g), main="Degree Distribution", xlab="Degree", ylab="Frequency")

plot(c(0:36), degree_distribution(g, cumulative = FALSE), 
     type = "b", main ="Degree distribution", 
     xlab="Degree", ylab ="p(degree)")

plot(c(0:36), degree_distribution(g, cumulative = TRUE), 
     type = "b", main = "Cumulative degree distribution", 
     xlab="Degree", ylab ="p(degree) > x")

# calculate cumulative degree distribution against log scale
plot(c(0:36), degree_distribution(g, cumulative = TRUE), 
     log="xy", type = "b", main = "Cumulative degree distribution (Log-Log)", 
     xlab="Degree (log)", ylab ="p(degree) > x (log)")

dev.off()
par(mfrow=c(1,1))

# average diameter between two nodes
cat("Average degree: ", mean(degree(g)))
cat("Average path length:", average.path.length(g, directed = TRUE))
cat("Diameter:",diameter(g, directed = TRUE))
# clustering of the network
cat("Global clustering coefficient:", transitivity(g, type="global"))

# compare to random network
n <- 33
p <- mean(degree(g)) / (n-1)
gRandom <- sample_gnp(n, p, directed = TRUE, loops = FALSE)
cat("\nMean degree of random net:", mean(degree(gRandom)))
cat("Average path length random:",mean_distance(gRandom, directed = TRUE, unconnected = FALSE))
cat("Average path length random 2:", average.path.length(gRandom, directed = TRUE))
cat("Global clustering coefficient random:",transitivity(gRandom))


### ERGMs ###

# change format of network for ergm package
D1 <- asNetwork(g)
set.vertex.attribute(D1, "university", as.character(universities))
set.vertex.attribute(D1, "year", as.character(yearGroup))

# set up model
model1 <- ergm(D1 ~ edges + mutual + nodefactor("university") + nodematch("university") + nodematch("year"))
summary(model1)

# Hypothesis 1: Reciprocal edges
plogis(coef(model1)[['edges']]) # null hypothesis
plogis(coef(model1)[['edges']] + coef(model1)[['mutual']]) # calculate the estimated probability of a reciprocal tie between two nodes

# Hypothesis 2: More friends if from university B/C
plogis(coef(model1)[['edges']] + coef(model1)[['nodefactor.university.2']]) # calculate the estimated probability of a tie between a node from university B and another university
plogis(coef(model1)[['edges']] + coef(model1)[['nodefactor.university.3']])

# Hypothesis 3: More friends if at same university than similar age
plogis(coef(model1)[['edges']] + coef(model1)[['nodematch.university']]) # calculate the estimated probability of a tie between a node of the same university and of the same age category
plogis(coef(model1)[['edges']] + coef(model1)[['nodematch.year']])

# Calculate goodness of fit
gofm1 <- gof(model1)
plot(gofm1)
