library(caper)

data <- read.csv("Skel mass data.csv")

tree <- read.nexus("nature05634-s2-revised_Lesley.nex")
tree <- tree[[1]] #Take only the tree with best dates

#Plot the tree
par(mfrow = c(1,1), mar = c(0,0,0,0) + 0.5)
plot(tree, show.tip.label = F)

#Plot the data
par(mfrow = c(1,1), mar = c(4,4,0,0) + 0.5)
plot(logS ~ logM, data = data)

#Run the PGLS
skel.cdat <- comparative.data(phy = tree,
                              data = data,
                              names.col = "Mesquite_Name")
m1 <- pgls(logS ~ logM + logM2, data = skel.cdat, lambda = "ML")
summary(m1) # can see here how branch length kappa = 1, how lambda = 0.31, and model equation
# log(Ms) = -1.21 +1.04 log(M) + 0.02[log*M)]^2 (slightly different from manuscript result)

x.vals <- seq(from = min(data$logM), to = max(data$logM), length.out = 100)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["logM"] * x.vals +
        coef(m1)["logM2"] * (x.vals ^ 2))
