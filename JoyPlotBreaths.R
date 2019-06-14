# try ggjoy 
library(ggplot2)
library(ggjoy)
library(plyr)

# 
# ggplot(diamonds, aes(x = price, y = cut, fill = cut)) + 
#   geom_joy(scale = 4) + theme_joy() +
#   scale_fill_cyclical(values = c("blue", "green")) +
#   scale_x_continuous(expand = c(0, 0))  +
#   scale_y_discrete(expand = c(0.01, 0))

# load in data
IBI <- read.csv('JoyPlotBreathDatTest_export.csv') 
colnames(IBI) <- c("weight","sppweight","IBI","spp","file","weightfile")

# RBG to Hex for Colormap --  7 species
x <- c("68 1 84",
       "68 58 131",
       "49 104 142",
       "33 145 140",
       "53 183 121",
       "144 215 67",
       "253 231 37")
cmap <- sapply(strsplit(x, " "), function(x)
  rgb(x[1], x[2], x[3], 150, maxColorValue=255))

# IBI versus WEIGHT 
r2<-ddply(IBI, .(weightfile), summarize, mean=mean(IBI))
colnames(r2)<- c("weightfile","IBI")
r1<-ddply(IBI, .(weightfile), summarize, median=median(IBI))
colnames(r1)<- c("weightfile","IBI")

png('JoyPlot_test.png', width = 8, height = 6, units = "in", res = 300)
ggplot(IBI, aes(x = log10(IBI), y = as.factor(weightfile), fill = as.factor(spp))) +
 geom_joy(scale = 10) + theme_joy() +
  scale_fill_cyclical(values = cmap) +  # get values for VERIDIS
 scale_y_discrete(expand = c(0.01, 0)) +   # will generally have to set the `expand` option
 scale_x_continuous(limits = c(0, 4)) +      # for both axes to remove unneeded padding
  labs(x="log(Inter Breath Interval) (sec)",y="Individual Weights")   
dev.off()

# IBI$dive = c(2*60,5*60,45*60,20*60,3600,15*60,3600)
r2<-ddply(IBI, .(spp), summarize, mean=mean(IBI))
colnames(r2)<- c("spp","IBI")
r1<-ddply(IBI, .(spp), summarize, median=median(IBI))
colnames(r1)<- c("spp","IBI")
library(prettyR)
r3<-ddply(IBI, .(spp), summarize, mode=Mode(IBI))
colnames(r3)<- c("spp","IBI")
r3$IBI <- as.numeric(r3$IBI)

IBI$spp <- as.factor(IBI$spp)

ggplot(IBI, aes(x = log10(IBI), y = spp, fill = spp, height = ..density..))+
  geom_density_ridges(scale = 4, stat = "density") +
  scale_fill_cyclical(values = cmap) +  # get values for VERIDIS
  scale_y_discrete(expand = c(0.01, 0))+    # will generally have to set the `expand` option
  scale_x_continuous(limits = c(0, 4)) + # for both axes to remove unneeded padding
  #geom_point(data=r2, cex = 2) +
  #geom_point(data=r1, shape = 21, cex = 2) +
  #geom_point(data=r3,shape = 22, cex = 2) +
  labs(x="log(Inter Breath Interval) (sec)",y="Species") 

ggplot(IBI, aes(x = IBI, y = as.factor(spp), fill = as.factor(spp), height = ..density..)) +
  geom_density_ridges(scale = 4, stat = "density") +
  scale_fill_cyclical(values = cmap) +  # get values for VERIDIS
  scale_y_discrete(expand = c(0.01, 0)) +   # will generally have to set the `expand` option
  scale_x_continuous(limits = c(0,120), breaks = c(0,10,60,120)) + # for both axes to remove unneeded padding
  labs(x="Inter Breath Interval (sec)",y="Species")                       

## plot only gm11_248c
gm11_148c <- IBI[which(IBI$file == "59"),]
r2<-ddply(gm11_148c, .(spp), summarize, mean=mean(IBI))
colnames(r2)<- c("spp","IBI")
r1<-ddply(gm11_148c, .(spp), summarize, median=median(IBI))
colnames(r1)<- c("spp","IBI")
library(prettyR)
r3<-ddply(gm11_148c, .(spp), summarize, mode=Mode(IBI))
colnames(r3)<- c("spp","IBI")
r3$IBI <- as.numeric(r3$IBI)

ggplot(gm11_148c, aes(x = (IBI), y = as.factor(weight), fill = as.factor(spp))) +
  geom_joy(scale = 10) + theme_joy() +
  scale_fill_cyclical(values = cmap[5]) +  # get values for VERIDIS
  # geom_point(data=r2, cex = 2) +
  scale_y_discrete(expand = c(0.01, 0)) +   # will generally have to set the `expand` option
  scale_x_continuous(limits = c(0, 300),breaks = c(0,16,34,120)) +      # for both axes to remove unneeded padding
  # geom_point(data=r1, shape = 21, cex = 2) +
  # geom_point(data=r3,shape = 22, cex = 2) + 
   labs(x="Inter Breath Interval (sec)",y="Individual Weights")   

ggplot(gm11_148c, aes(x = log10(IBI), y = as.factor(weight), fill = as.factor(spp))) +
  geom_joy(scale = 10) + theme_joy() +
  scale_fill_cyclical(values = cmap[5]) +  # get values for VERIDIS
  scale_y_discrete(expand = c(0.01, 0)) +   # will generally have to set the `expand` option
  scale_x_continuous(limits = c(0, 5)) +      # for both axes to remove unneeded padding
  labs(x="log(Inter Breath Interval) (sec)",y="Individual Weights")   


# load in other lit values
Lit <- read.csv("LitRespData.csv")

RoosEstdist <- data.frame()
for (n in 11:20){
N <- round(Lit$Tag.duration[n]*Lit$Mean[n]) # should be total n 
mean <- Lit$Mean[n]
sd <- Lit$SD[n]
lowerBound <- Lit$Min[n]
upperBound <- Lit$Max[n] 

RoosEstdist$n <- distEst(N,mean,sd,lowerBound,upperBound)
}

