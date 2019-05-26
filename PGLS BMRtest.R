library(caper)

# load FMR data
data2 <- read.csv("FMRdata.csv")
colnames(data2) <- c("Family","Genus","Species","Common","Common1",
                     "M","FMR","Method","Source1","Source2")

png('FMRscaling_1.png', width = 8, height = 6, units = "in", res = 300)
plot(log10(FMR) ~ log10(M), data = data2[which(data2$Method == "measurement"),],
     xlim = c(-2,5.5), ylim = c(0.5,8.5), cex = 1, cex.axis = 1.5)
m1 <- lm(log10(FMR) ~log10(M), data = data2[which(data2$Method == "measurement"),])
summary(m1)

# add best fit line
x.vals <- seq(from = min(log(data2$M)), to = max(log(data$M)), length.out = 100)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["log10(M)"] * x.vals, lty = "dashed")
# lines(x = x.vals, y = coef(m1)["(Intercept)"] + 1 * x.vals, lty = "dashed")
dev.off()

# if want to make on log axes 
# lines(x = x.vals,
# y = log(coef(m1)["(Intercept)"]) * x.vals ^(coef(m1)["log10(M)"]), lty = "dashed")


# add marine mammals (tried to do automatically like [which(which(data2$Method != "estimate")<43),],pch = 19)
# but just couldn't get it)
png('FMRscaling_2.png', width = 8, height = 6, units = "in", res = 300)
plot(log10(FMR) ~ log10(M), data = data2[which(data2$Method == "measurement"),],
     xlim = c(-2,5.5), ylim = c(0.5,8.5), cex = 1, cex.axis = 1.5)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["log10(M)"] * x.vals, lty = "dashed")
points(log10(FMR) ~ log10(M), data = data2[c(1,9,14:21,24:34,36,38:43),],pch = 19, cex = 1.5)

# add Maresh proposal 
x.vals <- x.vals <- seq(from = min(1), to = max(log10(data$M[1:43])), length.out = 100)
lines(x = x.vals,
      y = log10(732) + 0.45 * x.vals)
dev.off()



#### load BMR data
dataBMR <- read.csv("BMRdata.csv")
colnames(dataBMR) <- c("Family","Genus","Species","Common","Common1","Common2",
                     "M","BMR","Method","Source1","Source2")
dataBMR$M <- as.numeric(dataBMR$M)

plot(log10(BMR) ~ log10(M), data = dataBMR)



#tree <- read.nexus("nature05634-s2-revised_Lesley.nex")
#tree <- tree[[1]] #Take only the tree with best dates

#Plot the tree
#par(mfrow = c(1,1), mar = c(0,0,0,0) + 0.5)
#plot(tree, show.tip.label = F)

#Plot the data
#par(mfrow = c(1,1), mar = c(4,4,0,0) + 0.5)
plot(logBMR ~ logM, data = data)

#Run the PGLS
BMR.cdat <- comparative.data(phy = tree,
                              data = data,
                              names.col = "Mesquite_Name")
m1 <- pgls(logBMR ~ logM, data = BMR.cdat, lambda = "ML")
summary(m1) # can see here how branch length kappa = 1, how lambda = 0.31, and model equation
# log(Ms) = -1.21 +1.04 log(M) + 0.02[log*M)]^2 (slightly different from manuscript result)

x.vals <- seq(from = min(data$logM), to = max(data$logM), length.out = 100)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["logM"] * x.vals)
