library(caper)

setwd("/Users/julievanderhoop/Documents/R/BreathCounts")
data <- read.csv("MortolaData.csv")

tree <- read.nexus("nature05634-s2-revised_Lesley.nex")
tree <- tree[[1]] #Take only the tree with best dates

#Plot the tree
#par(mfrow = c(1,1), mar = c(0,0,0,0) + 0.5)
#plot(tree, show.tip.label = F)

#Plot the data
par(mfrow = c(1,1), mar = c(4,4,0,0) + 0.5)
plot(log(f) ~ log(M), data = data)

#Run the PGLS
breathf.cdat <- comparative.data(phy = tree,
                              data = data,
                              names.col = "Latin_Name")
m1 <- pgls(log(f) ~ log(M), data = breathf.cdat, lambda = "ML")
summary(m1) # can see here how branch length kappa = 1, how lambda = 0.31, and model equation
# log(Ms) = -1.21 +1.04 log(M) + 0.02[log*M)]^2 (slightly different from manuscript result)

# to check which tips are not matched:
breathf.cdat$dropped$unmatched.rows
# search within tips:
# tree$tip.label[grep("Pelecanus", tree$tip.label)]

x.vals <- seq(from = min(data$logM), to = max(data$logM), length.out = 100)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["logM"] * x.vals)

## do this with RESP data -- FAKE DATA but want to test tree 
TESTbreath <- read.csv("BreathCounts_PGLSdata")
colnames(TESTbreath) <- c("M","IBI","spp")
TESTbreath$Latin_Name <- c("Phocoena_phocoena","Tursiops_truncatus",
                           "Mesoplodon_densirostris","Globicephala_macrorhynchus",
                           "Ziphius_cavirostris","Globicephala_melas",
                           "Physeter_catodon")
plot(log(60/IBI) ~ log(M), data = TESTbreath)

TESTbreath.cdat <- comparative.data(phy = tree,
                                 data = TESTbreath,
                                 names.col = "Latin_Name")
m1 <- pgls(log(60/IBI) ~ log(M), data = TESTbreath.cdat, lambda = "ML")
summary(m1) 



# to check which tips are not matched:
TESTbreath.cdat$dropped$unmatched.rows

x.vals <- seq(from = min(log(TESTbreath$M)), to = max(log(TESTbreath$M)), length.out = 100)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["log(M)"] * x.vals)

