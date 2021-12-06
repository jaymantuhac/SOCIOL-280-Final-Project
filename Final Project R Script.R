# SOCIOL 280: Final Project Script

#Load packages
library(pacman)
pacman::p_load(ggplot2, 
               sna,
               statnet,
               RSiena,
               psych,
               car,
               stringr,
               gdata,
               cluster,
               reshape2,
               dplyr,
               openxlsx,
               gridExtra,
               influenceR,
               igraph,
               cowplot,
               ggpubr,
               magick,
               ergMargins,
               tnet,
               tidyverse,
               coda,
               statnet)

nodes <- read.csv('nodeList_covidDummy.csv')
nodes <- nodes[!(nodes$Id == 55 | nodes$Id == 126),] #Nodes 55 and 126 don't appear on sociogram

#Edgelist to Sociomatrix Conversion

# specify the edgelist as an object called 'playNet'
edgelist <- read.csv('edgeList_covidDummy.csv')
head(edgelist)

# split the edgelist into edges (ego, alter) and edge attributes
if(dim(edgelist)[2] > 2) { 
  edges <- edgelist[,c(1,2)]
  edgeAttrs <- edgelist[,-c(1,2)] 
}

# specify which IDs to include in the network
#  - the code below uses all IDs in the edgelist
#  - an alternative would be to only specify nodes that meet some criteria
rawIDs <- sort(unique(unlist(edges)))    # the IDs to include in the network
( N <- length(rawIDs) )          # check that this is the number of nodes you expected
newIDs <- seq(1:length(rawIDs))  # new IDs
lookUp <- data.frame(rawIDs, newIDs)    # link the old IDs to the new IDs
elNew <- matrix(NA, nrow=dim(edges)[1], ncol=2)     # create a new edgelist object filled with NA
elNew[] <- lookUp$newIDs[match(unlist(edges), lookUp$rawIDs)]   # recode old IDs into new IDs

# bring any attributes back into the edgelist
if(dim(edgelist)[2] > 2) { 
  elNew <- data.frame(elNew, edgeAttrs) 
}

# remove edges that weren't recoded
#  - this occurs if you specify rawIDs that are a subset of nodes in the edgelist
elNew <- elNew[apply(elNew, 1, function(x) sum(is.na(x)))==0,]  

# examine the new edgelist 	
head(elNew)
tail(elNew)

# turn the edgelist into a matrix (called 'mat')
mat <- matrix(0, nrow=N, ncol=N)  # create a NxN sociomatrix filled with 0's

# populate the matrix
for (i in 1:dim(elNew)[1]) {
  mat[elNew[i,1],elNew[i,2]] <- 1           # each edge is coded as 1
  #	mat[elNew[i,1],elNew[i,2]] <- elNew[i,3]  # or, to use a third column with edgeweights
}

# examine the sociomatrix
table(mat, useNA='always')  # are the number of ties what you expect?	
table(diag(mat), useNA='always')  # are there any self-ties - there are 2 in this
diag(mat) <- 0    # set self-ties to 0
gplot(mat, label=rawIDs)

# Note: Matrix has 263 nodes and 297 ties

sociomatrix <- mat

#Initial Visualizations
N <- dim(sociomatrix)[1]

race <- rep("NA", length.out = N)
race[nodes[, 'race.ethnicity'] == 'White'] <- 'white'
race[nodes[, 'race.ethnicity'] == 'Black'] <- 'gray50'
race[nodes[, 'race.ethnicity'] == 'NA'] <- 'gray90'

raceNum <- rep("NA", length.out = N)
raceNum[nodes[, 'race.ethnicity'] == 'White'] <- 1
raceNum[nodes[, 'race.ethnicity'] == 'Black'] <- 2

gender <- rep("NA", length.out = N)
gender[nodes[, 'gender'] == 'male'] <- 'green'
gender[nodes[, 'gender'] == 'female'] <- 'yellow'
gender[nodes[, 'gender'] == 'NA'] <- 'red'

genderref <- rep("NA", length.out = N)
genderref[nodes[, 'gender'] == 'male'] <- 1
genderref[nodes[, 'gender'] == 'female'] <- 2

zip <- rep("NA", length.out = N)
zip[nodes[, 'Zip.Code'] == 30306] <- '#9c482f'
zip[nodes[, 'Zip.Code'] == 30324] <- '#e63192'
zip[nodes[, 'Zip.Code'] == 30329] <- '#af0efa'
zip[nodes[, 'Zip.Code'] == 30338] <- '#771358'
zip[nodes[, 'Zip.Code'] == 30340] <- '#649832'
zip[nodes[, 'Zip.Code'] == 30341] <- '#e6d2b6'
zip[nodes[, 'Zip.Code'] == 30602] <- '#efcf4f'
zip[nodes[, 'Zip.Code'] == 30605] <- '#69e32b'
zip[nodes[, 'Zip.Code'] == 32606] <- '#490362'
zip[nodes[, 'Zip.Code'] == 35064] <- '#8b3292'
zip[nodes[, 'Zip.Code'] == 37116] <- '#4f06a8'

gplot(sociomatrix, vertex.col = race)
gplot(sociomatrix, vertex.col = gender)
gplot(sociomatrix, vertex.col = zip) #Visaulization reveals a large scale outbreak at zip code 30338

# Measures of centrality

outdegree <- sna::degree(sociomatrix, cmode = 'outdegree')
mean(outdegree) # Potential proxy for R0 value of the COVID-19 alpha variant?
hist(outdegree) 
as.data.frame(outdegree) %>% 
  ggplot(aes(x = outdegree)) +
  geom_histogram(color = "black",
                 fill = "white") +
  geom_vline(aes(xintercept = mean(outdegree)),
             color = "blue",
             size = 1) +
  labs(title = "Outdegree Distribution",
       x = "Outdegree",
       y = "Frequency")

indegree <- sna::degree(sociomatrix, cmode = 'indegree')
mean(indegree)
hist(indegree) 

# No histogram for indegree, may be too similar to the outdegree

sna::betweenness(sociomatrix, cmode = "directed")

sna::dyad.census(sociomatrix) #295 asymmetric dyads - 295 pairs of infections

## Subgroups

clique_census <- clique.census(sociomatrix, clique.comembership = "bysize")
clique_census$clique.count
clique_census$cliques

# Notes from meeting with Prof Schaefer 11/18
# Dive into how this data came to be
# Calculate homophily on each of the attributes
# Use ERGMs to put in multiple terms (e.g. gender, race, etc.)
# Do outdegree distribution of the zip code of interest

#Deep dive into zip code 30338
zip_30338 <- nodes %>% 
  filter(nodes$Zip.Code == 30338) %>% 
  select(Id) #Selects nodes with zip code 30338

outdegree_30338 <- sna::degree(sociomatrix, nodes = unlist(zip_30338), cmode = 'outdegree')
mean(outdegree_30338) #R0 of 0.93 (the infection should be dying down???)
as.data.frame(outdegree_30338) %>% 
  ggplot(aes(x = outdegree_30338)) +
  geom_histogram(color = "black",
                 fill = "white") +
  geom_vline(aes(xintercept = mean(outdegree_30338)),
             color = "blue",
             size = 1) +
  labs(title = "Outdegree Distribution of Zip Code 30338",
       x = "Outdegree",
       y = "Frequency") 

indegree_30338 <- sna::degree(sociomatrix, nodes = unlist(zip_30338), cmode = 'indegree')
mean(indegree) #Indegree of 1.13 
as.data.frame(indegree_30338) %>% 
  ggplot(aes(x = indegree_30338)) +
  geom_histogram(color = "black",
                 fill = "white") +
  geom_vline(aes(xintercept = mean(indegree_30338)),
             color = "blue",
             size = 1) +
  labs(title = "Outdegree Distribution of Zip Code 30338",
       x = "Outdegree",
       y = "Frequency") 

# Individual level homophily measures for zip code

sociomatrixNA <- sociomatrix
sociomatrixNA[sociomatrixNA == 0] <- NA

#Make zip code attribute numeric
zipref <- rep(0, length.out = N)
zipref[nodes[, 'Zip.Code'] == 30306] <- 30306
zipref[nodes[, 'Zip.Code'] == 30324] <- 30324
zipref[nodes[, 'Zip.Code'] == 30329] <- 30329
zipref[nodes[, 'Zip.Code'] == 30338] <- 30338
zipref[nodes[, 'Zip.Code'] == 30340] <- 30340
zipref[nodes[, 'Zip.Code'] == 30341] <- 30341
zipref[nodes[, 'Zip.Code'] == 30602] <- 30602
zipref[nodes[, 'Zip.Code'] == 30605] <- 30605
zipref[nodes[, 'Zip.Code'] == 32606] <- 32606
zipref[nodes[, 'Zip.Code'] == 35064] <- 35064
zipref[nodes[, 'Zip.Code'] == 37116] <- 37116
zipref[is.na(nodes[, 'Zip.Code'])] <- 0

alter_zip <- t(t(sociomatrixNA) * zipref)
ego_zip <- sociomatrix * zipref 

# now compare ego's value to alter's value to see if they match exactly 
ego_zip==alter_zip  # trues and falses when ties exist; NAs if no tie
table(ego_zip==alter_zip)  # Note: each tie counted twice
rowSums(ego_zip==alter_zip, na.rm=T)  # this gives us the number of alters who are similar
rowSums(ego_zip==alter_zip, na.rm=T) / rowSums(sociomatrix)  # and here is the proportion similar
mean_rowsums <- mean(rowSums(ego_zip==alter_zip, na.rm=T) / rowSums(sociomatrix), na.rm=T)  # mean proportion similar
hist(rowSums(ego_zip==alter_zip, na.rm=T) / rowSums(sociomatrix))
abline(v = mean_rowsums, col = 'blue', lwd = 2)

# what would we get by chance?
# permute attribute and recalculate, several times
# - percent similar on categorical attribute (race)
nPerms <- 10000
sims <- rep(0,nPerms)
for (i in 1:nPerms) {
  altersZipPerm <- t(t(sociomatrixNA) * sample(zipref, N))   # the sample function effectively permutes the zip code vector
  sims[i] <- mean(rowSums(ego_zip==alter_zip, na.rm=T) / rowSums(sociomatrix), na.rm=T)  # save mean similarity
}
mean(sims)   # proportion of alters similar expected by chance
(obsSim <- mean(rowSums(ego_zip==alter_zip, na.rm=T) / rowSums(sociomatrix), na.rm=T))  # observed proportion similar
hist(sims, xlim=c(min(sims, obsSim), max(sims,obsSim)))
abline(v=obsSim, col='red', lwd=2)
mean(obsSim < sims)   #  prop. of tests where observed similarity fell below expected (rare)
mean(obsSim > sims)   #  prop. of tests where observed similarity was greater than expected (most of the time)

#In short, not sure if the MCMC gives any insight. Move forward with calculating the EI Index

# - Spatial segregation index (see Bojanowski & Corten 2014)
#install.packages("netseg")
library(netseg)
library(igraph)
library(network)
sociomatrix.ig <- graph.adjacency(sociomatrix, mode=c("directed"), diag=FALSE)
fr <- layout.fruchterman.reingold(sociomatrix.ig)  
vertex_attr(sociomatrix.ig) <- list(zip = zipref,
                                    color = zip)

# Calculate homophily at the network level & adjust for baseline homophily
#  For a categorical attribute
#  - Calculate the E-I index 
net2 <- network::network(sociomatrix)  # convert into a "network" network object (compatible w/ sna package)
net9 <- network::network(1-sociomatrix)  # flip all ties
(oT <- mixingmatrix(net2,"zip") ) # observed ties - shows how many ties are between individuals of each race
nT <- mixingmatrix(net9,"zip")  # null ties
ei(sociomatrix.ig,"zip") # range: -1 (all ties within) to 1 (all ties between)
# Baseline EI of 0.26 - implies a good number of homophilous ties

#  - Calculate an odds ratio
#    Odds of observing homophilous tie relative to odds of observing heterophilous tie
network::set.vertex.attribute(net2,"zip", zipref)
network::set.vertex.attribute(net9,"zip", zipref)
obsOR <- (sum(diag(oT)) * (sum(nT)-sum(diag(nT)))) / 
  (sum(diag(nT)) * (sum(oT)-sum(diag(oT))))
print(obsOR)  # odds of homophilous tie 1.71 times greater than odds of heterophilous tie

orwg(sociomatrix.ig,"zip") #We get the same OR as the code above

# ERGM Development

#Implement other vertex attributes
network::set.vertex.attribute(net2, "race", race)
network::set.vertex.attribute(net2, "raceNum", raceNum)
network::set.vertex.attribute(net2, "gender", gender)
network::set.vertex.attribute(net2, "genderref", genderref)

model1 <- ergm(net2 ~ edges)
summary(model1)

model2 <- ergm(net2 ~ edges + nodefactor('zip'))
summary(model2)

model3 <- ergm(net2 ~ edges + nodefactor('zip') + nodematch('zip')) 
summary(model3) #Best fit model (lowest AIC) - gender/race don't really contribute to odds of a tie

model4 <- ergm(net2 ~ edges + nodefactor('zip') + nodematch('zip') +
                 nodematch('gender')) 
summary(model4) 

model5 <- ergm(net2 ~ edges + nodefactor('zip') + nodematch('zip') +
                 nodematch('gender') + nodematch('race')) 
summary(model5) #Most comprehensive model

model6 <- ergm(net2 ~ edges + nodefactor('zip') +
                 nodematch('gender') + nodematch('race'))
summary(model6)

#Generate odds/probabilities for each covariate
inv.logit <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

coef(model5) #odds of tie observation
inv.logit(coef(model5)) #Probability of tie observation

#Comparison with simulated model1 (model with just edges)
net_densities_2 <- unlist(lapply(hundred_simulations_2, network.density))

hist(net_densities_2, xlab = "Density", main = "", col = "lightgray")
abline(v = network.density(net2), col = "red", lwd = 3, lty = 2)
abline(v = mean(net_densities_2), col = "blue", lwd = 3, lty = 1)
# Model 5 significantly deviates from what would be simulated, suggests there truly is a better fit

# Test ERGM godness of fit of most comprehensive model compared to model with just edges
# gof() compares the fit of our model to simulated version of our model

gof_stats <- gof(model5)
plot(gof_stats)

# Goodness of fit reveals that model 5 has more deviation than simulation compared to the null (model1) model
# Model 5 is much more comprehensive and is a better fit with our data!

