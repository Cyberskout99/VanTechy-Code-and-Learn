library(dplyr)
library(sqldf)
library(psych)

#### Import Data Set ####
ALL_gaits2a <- read.csv("ALL_dataset2.csv",na.strings=c("","NA"), fileEncoding="UTF-8-BOM")

####  Principal Components Example ####
## drop out character variable "Grouping" (column #3) ##

ALL_gaits2 <- select(ALL_gaits2a, c(-3))

cor(ALL_gaits2)  ## quick output of correlations in console ('stats' package)
corr.test(ALL_gaits2)  ## reports all correlations & p-values ('psych' package)

View(cor(ALL_gaits2))  ## Usually works to display data set in a "spreadsheet" format
View(corr.test(ALL_gaits2))  ## Doesn't work for all data objects, unfortunately

GaitsPCA <- princomp(ALL_gaits2, cor=T)
GaitsPCA ##  Eigenvalues 

summary(GaitsPCA)  ##  Proportion of variance explained
loadings(GaitsPCA) ##  Loading values for each Eigenvector

graphics.off()
par("mar")
par(mar=c(1,1,1,1))

plot(GaitsPCA,type="lines") ## scree plot
biplot(GaitsPCA)  ## dots & arrows plot

#### Gaits K-means Clustering ####

## Step #1 - Create a matrix of the clustering data ##
## Can specify either individual variables or column ranges ##
## If specifying column range, no need to rename columns later ##

#### Residualized Gait Data (factor out Age), PCA ####

lmVeloc <- lm(VelocNorm~Age, data=ALL_gaits2)
resVeloc <- residuals(lmVeloc)

lmCad <- lm(CadenceNorm~Age, data=ALL_gaits2)
resCad <- as.numeric(residuals(lmCad))

lmStride <- lm(StrideRNorm~Age, data=ALL_gaits2)
resStride <- residuals(lmStride)

lmWalkR <- lm(WalkRatio~Age, data=ALL_gaits2)
resWalkR <- residuals(lmWalkR)

lmCycTime<- lm(CycleTimeR~Age, data=ALL_gaits2)
resCycTime <- residuals(lmCycTime)

lmSSupP <- lm(SingleSuppR~Age, data=ALL_gaits2)
resSSupP <- residuals(lmSSupP)

lmSSupT <- lm(S_Supp_Time_R~Age, data=ALL_gaits2)
resSSupT <- residuals(lmSSupT)

lmDSupP <- lm(DoubleSuppR~Age, data=ALL_gaits2)
resDSupP <- residuals(lmDSupP)

lmHeelOOT <- lm(HeelOffOnTimeR~Age, data=ALL_gaits2)
resHeelOOT <- residuals(lmHeelOOT)

lmHeelOOP <- lm(HeelOffOnPercR~Age, data=ALL_gaits2)
resHeelOOP <- residuals(lmHeelOOP)

resGaits <- cbind(resVeloc,resCad, resStride, resWalkR,
                  resCycTime, resSSupP, resSSupT, resDSupP,
                  resHeelOOT, resHeelOOP) ## Create matrix of residualized data ##


resGaitsIDs <- select(ALL_gaits2a, c(1:3)) ## Extract PTID, Age, & Grouping from original data

## Step #2 - Standardize variables in matrix, update row names ##

resGaits2 <- scale(resGaits) # Apply scaling to residualized data
View(resGaits2)  ##  Data has been scaled in terms of Standard Deviation
rownames(resGaits2) <- resGaitsIDs$PTID  ##  Add handles to residualized data

## Step #3 - Run the K-Means algorithm, specify 3 clusters (a priori) ##
model <- kmeans(resGaits2, 3)
print(model)

## breakdown ##
model$cluster  ## clustering membership ##
model$size  ## number of units within clusters ##
model$centers  ## reports value centers ##
model$totss ## reports total sum of squares ##
model$withinss ## reports within sum of squares ##
model$tot.withinss ## total within-cluster sum of squares ##

## plot cluters & scenters via scatterplot ##
centers <- data.frame(model$centers)

## Add Timepoint Data to matrix of Scaled cluster data ##
resGaits4 <- cbind(resGaitsIDs,resGaits2)
View(resGaits4)
names(resGaits4)[3] <- c("Grouping")
resGaits4$Grouping <- as.factor(resGaits4$Grouping)

## Scatterplot of specific data, also display centers ##
require(ggplot2)
ggplot()+geom_point(data=resGaits4, aes(x=resCycTime, y=resVeloc, color = Grouping))+
  geom_point(data=centers, aes(x=resCycTime[1], y=resVeloc[1]), size=10, color="Red", shape=8)+
  geom_point(data=centers, aes(x=resCycTime[2], y=resVeloc[2]), size=10, color="Blue", shape=8)+
  geom_point(data=centers, aes(x=resCycTime[3], y=resVeloc[3]), size=10, color="Yellow", shape=8)+
  theme(legend.position="right")+
  labs(color="Grouping")

## Can add center labels to the scatterplot using "geomtext" function ##
ggplot()+geom_point(data=resGaits4, aes(x=resCycTime, y=resVeloc, color = Grouping))+
  geom_point(data=centers, aes(x=resCycTime[1], y=resVeloc[1]), size=8,color="Red", shape=8)+
  geom_text(data=centers, aes(x=resCycTime[1], y=resVeloc[1]-0.09), label="Cluster 1")+
  geom_point(data=centers, aes(x=resCycTime[2], y=resVeloc[2]), size=8,color="Blue", shape=8)+
  geom_text(data=centers, aes(x=resCycTime[2], y=resVeloc[2]-0.09), label="Cluster 2")+
  geom_point(data=centers, aes(x=resCycTime[3], y=resVeloc[3]), size=8, color="Yellow", shape=8)+
  geom_text(data=centers, aes(x=resCycTime[3], y=resVeloc[3]-0.09), label="Cluster 3")+
  labs(color="Grouping")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey"),panel.grid.minor = element_line(colour = "lightgrey"))+
  ylab("Residual Velocity")+
  xlab("Residual Cycle Time")

ggplot()+geom_point(data=resGaits4, aes(x=resCycTime, y=resVeloc, color = Grouping))+
  geom_point(data=centers, aes(x=resCycTime[1], y=resVeloc[1], size=10), color="Red", shape=8)+
  geom_point(data=centers, aes(x=resCycTime[2], y=resVeloc[2], size=10), color="Blue", shape=8)+
  geom_point(data=centers, aes(x=resCycTime[3], y=resVeloc[3], size=10), color="Yellow", shape=8)+
  theme(legend.position="right")

## 'ggfortify' library, add lassoes to the scatterplot ##
library(ggfortify)

## PCA plot via 'ggplot2' package, base R 'stats' package ##
## Loadings option displays eigenvectors ##
autoplot(prcomp(resGaits2), data=resGaits4, colour='Grouping');  ## datapoints 
autoplot(prcomp(resGaits2), data=resGaits4, colour='Grouping', loadings = T, label=T) ## with vectors

## K-means clusters via 'ggplot2' package, base R 'stats' package ##
## Generalizes full dataset instead of specific attributes ##
autoplot(kmeans(resGaits2,3), data=resGaits3) ## colored clusters 
autoplot(kmeans(resGaits2,3), data=resGaits3, label=T, label.size=3) ## with labels

## Lassoes around the clusters, 'cluster' package ##
library(cluster)
autoplot(pam(resGaits4[-5],3), data=resGaits2, frame=TRUE, frame.type='norm', label=T)
