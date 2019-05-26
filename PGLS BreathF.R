library(caper)

setwd("/Users/julievanderhoop/Documents/R/BreathCounts")
data <- read.csv("MortolaData.csv")

tree <- read.nexus("nature05634-s2-revised_Lesley.nex")
tree <- tree[[1]] #Take only the tree with best dates

#Plot the tree
#par(mfrow = c(1,1), mar = c(0,0,0,0) + 0.5)
#plot(tree, show.tip.label = F)

#Plot the data
# par(mfrow = c(1,1), mar = c(4,4,0,0) + 0.5)
nonrum <- data[which(data$Ruminant == 0 & data$Marine == 0 & data$Aquatic.Bird == 0 & data$Terrestrial.Bird == 0),]
rum <- data[which(data$Ruminant == 1 & data$Marine == 0 & data$Aquatic.Bird == 0 & data$Terrestrial.Bird == 0),]
png('Mortola_NonRuminants.png', width = 8, height = 6, units = "in", res = 300)
plot((f) ~ (M), data = nonrum, log = "xy", cex = 1.5, cex.axis = 1.5)
m <- lm(log(f) ~ log(M), data = nonrum)
x.vals <- seq(from = min(nonrum$M), to = max(nonrum$M), length.out = 100)
lines(x = x.vals,
      y = exp(coef(m)["(Intercept)"]) * x.vals^coef(m)["log(M)"] )
dev.off()
# or plot in ggplot
library(ggplot2)
library(scales)     # Need the scales package
ggplot(nonrum, aes(x = M, y = f)) +
  geom_point() +
  scale_y_continuous(trans=log10_trans()) + 
  scale_x_continuous(trans=log10_trans()) + 
  stat_smooth(method = "lm", col = "gray")


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
# plot(breathf.cdat$phy) # plot tree

x.vals <- seq(from = min(data$logM), to = max(data$logM), length.out = 100)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["logM"] * x.vals)


## do this with RESP data -- FAKE DATA but want to test tree 
TESTbreath <- read.csv("/Users/julievanderhoop/Documents/MATLAB/BreathCounts/BreathCounts_PGLSdata")
colnames(TESTbreath) <- c("Mb","mnIBI","mdIBI","spp","mnFreq")
TESTbreath$Latin_Name <- c("Phocoena_phocoena",
                           "Tursiops_truncatus",
                           "Globicephala_macrorhynchus",
                           "Mesoplodon_densirostris",
                           "Globicephala_melas",
                           "Ziphius_cavirostris",
                           "Physeter_catodon",
                           "Hyperoodon_ampullatus",
                           "Delphinapterus_leucas",
                           "Balaenoptera_acutorostrata",
                           "Grampus_griseus",
                           "Megaptera_novaeangliae",
                           "Eschrichtius_robustus",
                           "Eubalaena_glacialis",
                           "Balaenoptera_physalus",
                           "Balaenoptera_musculus")

TESTbreath.cdat <- comparative.data(phy = tree,
                                 data = TESTbreath,
                                 names.col = "Latin_Name")

# do some plotting 
branches <- TESTbreath.cdat$phy$edge
species <- TESTbreath.cdat$phy$tip.label
brlength <- TESTbreath.cdat$phy$edge.length
nodes <- TESTbreath.cdat$phy$Nnode

mycol<-c("blue", "blue", "blue", "red", "red", "red","blue","red") # dummy
p1 <- plot(TESTbreath.cdat$phy, adj=0, label.offset=1.75, lwd=2)
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=2)
#nodelabels()
#tiplabels()
#rot.phy <- rotate(TESTbreath.cdat$phy, node=18)
#ladder.phy <- ladderize(TESTbreath.cdat$phy)

# fit PGLS to breath data
m1 <- pgls(log(60/mnIBI) ~ log(Mb), data = TESTbreath.cdat, lambda = "ML")
summary(m1) 

# fit alternative models
m05 <- pgls(log(60/mnIBI) ~ offset(-0.5*log(Mb)), data = TESTbreath.cdat, lambda = 1)
m25 <- pgls(log(60/mnIBI) ~ offset(-0.25*log(Mb)), data = TESTbreath.cdat, lambda = 1)

anova(m1, m05, m25)
AIC(m1, m05, m25)

# to check which tips are not matched:
TESTbreath.cdat$dropped$unmatched.rows

x.vals <- seq(from = min(log(TESTbreath$Mb)), to = max(log(TESTbreath$Mb)), length.out = 100)

plot(log(60/mnIBI) ~ log(Mb), data = TESTbreath)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["log(M)"] * x.vals)
# add lines for 0.5 and 0.25 
lines(x = x.vals,
      y = 20*coef(m05)["(Intercept)"] +
        -0.5 * x.vals, lty = 2)
lines(x = x.vals,
      y = 10*coef(m25)["(Intercept)"] +
        -0.25 * x.vals, lty = 2)

# plot mortola data and ours
png('Mortola_NonRuminantsCetaceans.png', width = 8, height = 6, units = "in", res = 300)
plot((f) ~ (M), data = nonrum, log = "xy", 
     cex = 1, cex.axis = 1.5, 
     xlim = c(0.01,80000), ylim = c(0.1, 200))
x.vals <- seq(from = min(nonrum$M), to = max(nonrum$M), length.out = 100)
lines(x = x.vals,
      y = exp(coef(m)["(Intercept)"]) * x.vals^coef(m)["log(Mb)"] )
points(60/mnIBI ~ Mb, data = TESTbreath, pch = 19, cex = 1.5)
m1 <- lm(log(60/mnIBI) ~ log(Mb), data = TESTbreath)
x.vals <- seq(from = min(TESTbreath$Mb), to = max(TESTbreath$Mb), length.out = 100)
lines(x = x.vals,
      y = exp(coef(m1)["(Intercept)"]) * x.vals^coef(m1)["log(Mb)"] )

dev.off()

