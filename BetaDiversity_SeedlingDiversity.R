library(reshape2)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(MuMIn)
library(lsmeans)
library(effects)
library(betapart)

###############################################################################
################################ Beta diversity ###############################
###############################################################################


betamatrix <- datat.cast[datat.cast$brn == "NO",9:60]
rownames(betamatrix) <- datat.cast[datat.cast$brn == "NO",]$ID

for  (i in 1:length(betamatrix$NONE))
{
  betamatrix[i,which(colnames(betamatrix) == "NONE")] <- ifelse(sum(betamatrix[i,]) == 0, 1, 0)
}

betamatrix <- betamatrix[betamatrix$NONE == 0,]
betamatrix <- as.data.frame(ifelse(betamatrix == 0, 0, 1))
betamatrix <- subset(betamatrix, select = -NONE)

beta1 <- beta.pair(betamatrix)

beta1.sne <- as.matrix(beta1$beta.sne)
beta1.sne.melt <- melt(data = beta1.sne)

vals <- strsplit(as.character(beta1.sne.melt$Var1), " ") 
x <- as.data.frame(do.call("rbind", vals)) 
beta1.sne.melt3 <- cbind(x, beta1.sne.melt)

vals <- strsplit(as.character(beta1.sne.melt$Var2), " ") 
x <- as.data.frame(do.call("rbind", vals)) 
nes.melt <- cbind(x, beta1.sne.melt3)
colnames(nes.melt) <- c("Site.x", "LUH.x", "CT.x", "SP.x", "ET.x", "EP.x", "Site.y", "LUH.y", "CT.y", "SP.y", "ET.y", "EP.y", "ID.x", "ID.y", "nes")

nes.melt <- nes.melt[nes.melt$Site.x == nes.melt$Site.y,]
nes.melt <- nes.melt[nes.melt$LUH.x == nes.melt$LUH.y,]
nes.melt <- nes.melt[nes.melt$CT.x == nes.melt$CT.y,]

nes.melt$pair <- ifelse(nes.melt$ET.x == "E" & nes.melt$ET.y == "P" , "A", "P") 
nes.melt$pair <- ifelse(nes.melt$ET.x == "P" & nes.melt$ET.y == "E" , "A", nes.melt$pair) 

hist(sqrt(nes.melt$nes))

mnes1 <- lmer(sqrt(nes) ~ LUH.x*CT.x*pair + (1|Site.x), data = nes.melt)
sumnes <- summary(mnes1)
plot(mnes1)
hist(resid(mnes1))
qqnorm(resid(mnes1))
lsnes <- lsmeansLT(mnes1)
plot(allEffects(mnes1))

