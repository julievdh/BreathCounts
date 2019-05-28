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
summary(m1) # can see here how branch length kappa = 1, how lambda = 0.898, and model equation
# log(IBI) = 3.8 + 0.33*log(M)

# to check which tips are not matched:
breathf.cdat$dropped$unmatched.rows
# search within tips:
# tree$tip.label[grep("Pelecanus", tree$tip.label)]
# plot(breathf.cdat$phy) # plot tree

x.vals <- seq(from = min(data$logM), to = max(data$logM), length.out = 100)
lines(x = x.vals,
      y = coef(m1)["(Intercept)"] +
        coef(m1)["logM"] * x.vals)


## do this with resp data output from matlab roughCountData.m 
breathdat <- read.csv("/Users/julievanderhoop/Documents/MATLAB/BreathCounts/BreathCounts_PGLSdata")
colnames(breathdat) <- c("Mb","mnIBI","mdIBI","spp","mnFreq","n")
breathdat$Latin_Name <- c("Phocoena_phocoena",
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

# subset into original n = 7
original <- subset(breathdat, spp <= 7)
# subset into odontocetes
odont <- subset(breathdat, spp <= 9 | spp == 11)
# subset into mysticetes
mysti <- subset(breathdat, spp == 10 | spp >= 12)
  
breath.cdat <- comparative.data(phy = tree,
                                 data = breathdat,
                                 names.col = "Latin_Name")
species <- breath.cdat$phy$tip.label

# mycol<-c("blue", "blue", "blue", "red", "red", "red","blue","red") # dummy
plot(breath.cdat$phy, adj=0, label.offset=5, lwd=2)
# tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=2)
tiplabels(breath.cdat$data$n,frame = "none", adj = -0.2)
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

# get residuals from model
d1 <- data.frame(order = c(1:16),id=breath.cdat$phy$tip.label) # these are the tip order
# but the data go in with a different order, so have to match data to tips 
mtch <- match(breathdat$Latin_Name,d1$id) # find matches between data frames
for (i in 1:16){
  breathdat$residuals[mtch == i] <- m1$residuals[i] # assign weight based on match ID
breathdat$order[mtch == i] <- d1$order[i]
}

# plotting
# library(phytools)
#phytools::plotTree.barplot(TESTbreath.cdat$phy,TESTbreath$residuals,
#        args.barplot=list(col=sapply(TESTbreath$residuals,
#        function(x) if(x>=0) "blue" else "red"),
#        xlim=c(-1.2,1.2)))
# phytools::plotTree.barplot(TESTbreath.cdat$phy,TESTbreath$residuals, xlim=c(-1.2,1.2))

# to check which tips are not matched:
breath.cdat$dropped$unmatched.rows

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

