library(caper)

setwd("/Users/julievanderhoop/Documents/R/BreathCounts")
data <- read.csv("MortolaData.csv")

tree <- read.nexus("nature05634-s2-revised_Lesley.nex")
tree <- tree[[1]] #Take only the tree with best dates

#Plot the tree
par(mfrow = c(1,1), mar = c(0,0,0,0) + 0.5)
plot(tree, show.tip.label = F)

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
tree$tip.label[grep("Pelecanus", tree$tip.label)]


x.vals <- seq(from = min(data$logM), to = max(data$logM), length.out = 100)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["logM"] * x.vals)
