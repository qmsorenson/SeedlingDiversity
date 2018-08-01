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


#############################################################################################
###############################  richness without burned  ###################################
#############################################################################################

data.cast <- read.csv("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/data.cast.csv")


rich.a <- lmer(logrich.a ~ ET*LUH*CT + (1|Site/LUH:CT/SP), data = data.cast[data.cast$brn == "NO",])
sra <- step(rich.a)
lsmra<- lsmeansLT(rich.a)
#plot(rich.a)
#hist(resid(rich.a))

rich.s <- lmer(logrich.s ~ ET*LUH*CT + (1|Site/LUH:CT/SP), data = data.cast[data.cast$brn == "NO",])
srs <- step(rich.s)
lsmrs<- lsmeansLT(rich.s)
#plot(rich.s)
#hist(resid(rich.s))

rich.v <- lmer(logrich.v ~ ET*LUH*CT + (1|Site/LUH:CT/SP), data = data.cast[data.cast$brn == "NO",])
srv <- step(rich.v)
lsmrv<- lsmeansLT(rich.v)
#plot(rich.v)
#hist(resid(rich.v))

save(lsmra, lsmrs, lsmrv, file = "G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/lsmr.RData")



######### Environmental variables ############

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/Enviro.data.RData")

###### all seedlings ######
rich.a.whc <- lmer(logrich.a ~ ET*LUH*CT*whc +  (1|Site/SP2), data = data.cast.whc)
sraw <- step(rich.a.whc) # ET*whc + LUH*CT
rich.a.whc <- lmer(logrich.a ~ ET*whc + LUH*CT +  (1|Site/SP2), data = data.cast.whc)
plot(rich.a.whc)
hist(resid(rich.a.whc))

rich.a.lig <- lmer(logrich.a ~ ET*LUH*CT*X5cm + (1|Site/LUH:CT/SP), data = ear.light)
sral <- step(rich.a.lig) #ET + LUH*CT + X5cm:CT + X5cm
rich.a.lig <- lmer(logrich.a ~ ET + LUH*CT + X5cm:CT + X5cm + (1|Site/LUH:CT/SP), data = ear.light)
plot(rich.a.lig)
hist(resid(rich.a.lig))

rich.a.comp<- lmer(logrich.a ~ ET*LUH*comp*CT + (1|Site/SP), data = data.castcomp)
srac <- step(rich.a.comp)
plot(rich.a.comp)
hist(resid(rich.a.comp))

rich.a.dbh <- lmer(logrich.a ~ ET*LUH*BA +(1|Site/LUH:CT/SP), data = data.cast.dbh[data.cast.dbh$CT == "U",])
anova(rich.a.dbh)
srad <- step(rich.a.dbh) # ET*LUH*dbh
plot(rich.a.dbh)
hist(resid(rich.a.dbh))

rich.a.gre <- lmer(logrich.a ~ ET*LUH*CT*green + (1|Site/LUH:CT/SP), data = ear.gre)
srag <- step(rich.a.gre) #ET + CT*LUH + gre
rich.a.gre <- lmer(logrich.a ~ ET + CT*LUH + green + (1|Site/LUH:CT/SP), data = ear.gre)
plot(rich.a.gre)
hist(resid(rich.a.gre))




save(rich.a.whc, rich.a.lig, rich.a.dbh, rich.a.gre, file = "G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/AllseededEnviroOutput.RData")


########## seeded species ##########

rich.s.whc <- lmer(logrich.s ~ ET*LUH*CT*whc + (1|Site/LUH:CT/SP), data = data.cast.whc)
srsw <- step(rich.s.whc)
plot(rich.s.whc)
hist(resid(rich.s.whc))

rich.s.lig <- lmer(logrich.s ~ ET*LUH*CT*X5cm + (1|Site/LUH:CT/SP), data = ear.light)
srsl <- step(rich.s.lig)
plot(rich.s.lig)
hist(resid(rich.s.lig))

rich.s.comp<- lmer(logrich.s ~ ET*LUH*comp*CT + (1|Site/LUH:CT/SP), data = data.castcomp)
srsc <- step(rich.s.comp)
plot(rich.s.comp)
hist(resid(rich.s.comp))

rich.s.dbh <- lmer(logrich.s ~ ET*LUH*BA+ (1|Site/LUH:CT/SP), data = data.cast.dbh[data.cast.dbh$CT == "U",])
anova(rich.s.dbh)
srsd <- step(rich.s.dbh)
plot(rich.s.dbh)
hist(resid(rich.s.dbh))

rich.s.gre <- lmer(logrich.s ~ ET*LUH*CT*green + (1|Site/LUH:CT/SP), data = ear.gre)
srsg <- step(rich.s.gre)
plot(rich.s.gre)
hist(resid(rich.s.gre))


########### volunteer species ############

rich.v.whc <- lmer(logrich.v ~ CT*LUH*ET*whc + (1|Site/LUH:CT/SP), data = data.cast.whc)
srvw <- step(rich.v.whc)
plot(rich.v.whc)
hist(resid(rich.v.whc))

rich.v.lig <- lmer(logrich.v ~ ET*LUH*CT*X5cm + (1|Site/LUH:CT/SP), data = ear.light)
srvl <- step(rich.v.lig)
anova(rich.v.lig)
plot(rich.v.lig)
hist(resid(rich.v.lig))

rich.v.comp<- lmer(logrich.v ~ ET*LUH*comp*CT + (1|Site/LUH:CT/SP), data = data.castcomp)
srvc <- step(rich.v.comp)
plot(rich.v.comp)
hist(resid(rich.v.comp))

rich.v.dbh <- lmer(logrich.v ~ ET*LUH*BA+ (1|Site/LUH:CT/SP), data = data.cast.dbh[data.cast.dbh$CT == "U",])
anova(rich.v.dbh)
srvd <- step(rich.v.dbh)
plot(rich.v.dbh)
hist(resid(rich.v.dbh))

rich.v.gre <- lmer(logrich.v ~ ET*LUH*CT*green + (1|Site/LUH:CT/SP), data = ear.gre)
srvg <- step(rich.v.gre)
plot(rich.v.gre)
hist(resid(rich.v.gre))




############################################################################################
#################################   MLM TRAITS #############################################
#############################################################################################

library("lme4")
load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/trait.dataframes.RData")

#binoomial models

m1 <- glmer(pa ~ LUH*CT*lnsm + ET + (1|Site/LUH:CT/SP2) + (1+ LUH*CT|Spp) + (1|obs), nAGQ = 0, family = "binomial", data=sm.melt)
summary(m1)
hist(resid(m1))
qqnorm(residuals(m1))
plot(allEffects(m1try))


m2 <- glmer(pa ~ LUH + CT:LUH + CT*ET*lnsla + (1|Site/LUH:CT/SP2) + (1+ ET*CT|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=sla.melt)

summary(m2)
hist(resid(m2))
qqnorm(residuals(m2))
plot(allEffects(m2))


m3 <- glmer(pa ~ LUH + LUH:CT + CT*ET*lnht + (1|Site/LUH:CT/SP2) + (1+ CT*ET|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=ht.melt)

summary(m3)
hist(resid(m3))
qqnorm(residuals(m3))
plot(allEffects(m3))

## moedls for plotting significant interactions
m1sm <- glmer(pa ~ LUH*CT*lnsm + (1|Site/LUH:CT/SP2) + (1+ LUH*CT|Spp) + (1|obs), nAGQ = 0, family = "binomial", data=sm.melt)

m2sla <- glmer(pa ~ CT*ET*lnsla + (1|Site/LUH:CT/SP2) + (1+ ET*CT|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=sla.melt)

m3ht <- glmer(pa ~ CT*ET*lnht + (1|Site/LUH:CT/SP2) + (1+ CT*ET|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=ht.melt)

save(m1sm, m2sla, m3ht, file = "G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/trait.models.RData")



#poisson models

sm_p1 <- glmer(value ~ LUH*CT*ET*lnsm - LUH:CT:ET:lnsm - LUH:ET:lnsm - LUH:CT:ET - LUH:ET - LUH:CT:lnsm -LUH:lnsm+ (1|Site/LUH:CT/SP2) + (1+ CT*ET|Spp) + (1|obs), nAGQ = 0, family = "poisson", data=sm.melt)

summary(sm_p1)

datat.melt$pa <- ifelse(datat.melt$value == 0, 0, 1)
datat.melt$obs <- 1:nrow(datat.melt)
datat.melt$SP2 <- paste(datat.melt$LUH, datat.melt$CT, datat.melt$SP)

#trait.cor <- left_join(sm, sla2, by = c("Spp" = "SppCode"))
#trait.cor <- left_join(trait.cor, ht, by = c("Spp" = "SppCode"))
ggplot(data.frame(x1=datat.melt$lnsm,pearson=residuals(m1,type="pearson")),
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw()

m8 <- glmer(pa ~ LUH*CT + ET + (0+CT|Spp) + (0+LUH|Spp) + (0+LUH:CT|Spp) + (1|Site/LUH:CT/SP2) + (1|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data = datat.melt)
summary(m8)
plot(allEffects(m7))

datat.melt2 <- aggregate(pa ~ Site + LUH)

m4 <- glmer(pa ~ LUH*CT +  ET + (1|Site/LUH:CT/SP2) + (1|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data = datat.melt)
summary(m4)
plot(allEffects(m4))


##############################################################################################################
#############################          Species Accumulation Curves           #################################
##############################################################################################################

datat.cast <- read.csv("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/data.cast.csv")

library(vegan)
betamatrixRTE <- datat.cast[datat.cast$LUH == "R" & datat.cast$CT == "T" & datat.cast$ET == "E" ,9:60]
specRTE <- specaccum(betamatrixRTE, permutations = 10000)

betamatrixRTP <- datat.cast[datat.cast$LUH == "R" & datat.cast$CT == "T" & datat.cast$ET == "P" ,9:60]
specRTP <- specaccum(betamatrixRTP)

betamatrixMTE <- datat.cast[datat.cast$LUH == "M" & datat.cast$CT == "T" & datat.cast$ET == "E" ,9:60]
specMTE <- specaccum(betamatrixMTE)

betamatrixMTP <- datat.cast[datat.cast$LUH == "M" & datat.cast$CT == "T" & datat.cast$ET == "P" ,9:60]
specMTP <- specaccum(betamatrixMTP)

betamatrixRUE <- datat.cast[datat.cast$LUH == "R" & datat.cast$CT == "U" & datat.cast$ET == "E" ,9:60]
specRUE <- specaccum(betamatrixRUE)

betamatrixRUP <- datat.cast[datat.cast$LUH == "R" & datat.cast$CT == "U" & datat.cast$ET == "P" ,9:60]
specRUP <- specaccum(betamatrixRUP)

betamatrixMUE <- datat.cast[datat.cast$LUH == "M" & datat.cast$CT == "U" & datat.cast$ET == "E" ,9:60]
specMUE <- specaccum(betamatrixMUE)

betamatrixMUP <- datat.cast[datat.cast$LUH == "M" & datat.cast$CT == "U" & datat.cast$ET == "P" ,9:60]
specMUP <- specaccum(betamatrixMUP)

plot(specRTE, ci = 0, col = "darkolivegreen3", lwd = 3, xlab = "Samples", ylab = "Number of species") #gold2, darkolivegreen
plot(specRTP, add = TRUE, ci = 0, col = "darkolivegreen3", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specMTE, add = TRUE, ci = 0, col = "gold2", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specMTP, add = TRUE, ci = 0, col = "gold2", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specRUE, add = TRUE, ci = 0, col = "darkolivegreen", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specRUP, add = TRUE, ci = 0, col = "darkolivegreen", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specMUE, add = TRUE, ci = 0, col = "gold4", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specMUP, add = TRUE, ci = 0, col = "gold4", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")

#Volunteer only


VbetamatrixRTE <- datat.cast[datat.cast$LUH == "R" & datat.cast$CT == "T" & datat.cast$ET == "E" ,names(datat.cast) %in% volunteerlist]
VspecRTE <- specaccum(VbetamatrixRTE, permutations = 10000)

VbetamatrixRTP <- datat.cast[datat.cast$LUH == "R" & datat.cast$CT == "T" & datat.cast$ET == "P" ,names(datat.cast) %in% volunteerlist]
VspecRTP <- specaccum(VbetamatrixRTP, permutations = 10000)

VbetamatrixMTE <- datat.cast[datat.cast$LUH == "M" & datat.cast$CT == "T" & datat.cast$ET == "E" ,names(datat.cast) %in% volunteerlist]
VspecMTE <- specaccum(VbetamatrixMTE, permutations = 10000)

VbetamatrixMTP <- datat.cast[datat.cast$LUH == "M" & datat.cast$CT == "T" & datat.cast$ET == "P" ,names(datat.cast) %in% volunteerlist]
VspecMTP <- specaccum(VbetamatrixMTP, permutations = 10000)

VbetamatrixRUE <- datat.cast[datat.cast$LUH == "R" & datat.cast$CT == "U" & datat.cast$ET == "E" ,names(datat.cast) %in% volunteerlist]
VspecRUE <- specaccum(VbetamatrixRUE, permutations = 10000)

VbetamatrixRUP <- datat.cast[datat.cast$LUH == "R" & datat.cast$CT == "U" & datat.cast$ET == "P" ,names(datat.cast) %in% volunteerlist]
VspecRUP <- specaccum(VbetamatrixRUP, permutations = 10000)

VbetamatrixMUE <- datat.cast[datat.cast$LUH == "M" & datat.cast$CT == "U" & datat.cast$ET == "E" ,names(datat.cast) %in% volunteerlist]
VspecMUE <- specaccum(VbetamatrixMUE, permutations = 10000)

VbetamatrixMUP <- datat.cast[datat.cast$LUH == "M" & datat.cast$CT == "U" & datat.cast$ET == "P" ,names(datat.cast) %in% volunteerlist]
VspecMUP <- specaccum(VbetamatrixMUP, permutations = 10000)

plot(VspecRTE, ci = 0, col = "darkolivegreen3", lwd = 3, xlab = "Samples", ylab = "Number of species") #gold2, darkolivegreen
plot(VspecRTP, add = TRUE, ci = 0, col = "darkolivegreen3", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(VspecMTE, add = TRUE, ci = 0, col = "gold1", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(VspecMTP, add = TRUE, ci = 0, col = "gold1", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(VspecRUE, add = TRUE, ci = 0, col = "darkolivegreen", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(VspecRUP, add = TRUE, ci = 0, col = "darkolivegreen", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(VspecMUE, add = TRUE, ci = 0, col = "gold3", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(VspecMUP, add = TRUE, ci = 0, col = "gold3", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")




data.cast.whc$SP2 <- paste(data.cast.whc$LUH, data.cast.whc$CT, data.cast.whc$SP)
ear.light$SP2 <- paste(ear.light$LUH, ear.light$CT, ear.light$SP)
data.castcomp$SP2 <- paste(data.castcomp$LUH, data.castcomp$CT, data.castcomp$SP)
data.cast.dbh$SP2 <- paste(data.cast.dbh$LUH, data.cast.dbh$CT, data.cast.dbh$SP)



##############################################################################################################
######################################          Trash bin           ##########################################
##############################################################################################################


#data.cast.whc <- data.cast.whc[!is.na(data.cast.whc$whc),]
#data.cast.whc$obs <- 1:nrow(data.cast.whc)
#data.cast.whc$obs <- as.character(data.cast.whc$obs)
#glmtry <- glmer(rich.a ~ ET*whc + CT*LUH + (1|Site/LUH:CT/SP) + (1|obs), data = data.cast.whc, family = "poisson")

