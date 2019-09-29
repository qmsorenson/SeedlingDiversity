library(reshape2)
library(lme4)
library(lmerTest)
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(MuMIn)
library(lsmeans)
library(effects)
library(phia)
library(gganimate)
library(piecewiseSEM)


#############################################################################################
#######################################   Richness   ########################################
#############################################################################################

data.cast <- read.csv("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/data.cast.csv")


rich.a <- lmer(logrich.a ~ ET*LUH*CT + (1|Site/LUH:CT) + (1|Site:LUH:CT:SP), data = data.cast[data.cast$brn == "NO",])
anova(rich.a)
summary(rich.a)
sra <- step(rich.a)
lsmra <- lsmeansLT(rich.a)
tab_model(rich.a)
#plot(rich.a)
#hist(resid(rich.a))

rich.s <- lmer(logrich.s ~ ET*LUH*CT + (1|Site/LUH:CT/SP), data = data.cast[data.cast$brn == "NO",])
anova(rich.s)
srs <- step(rich.s)
lsmrs<- lsmeansLT(rich.s)
#plot(rich.s)
#hist(resid(rich.s))

rich.v <- lmer(logrich.v ~ ET*LUH*CT + (1|Site/LUH:CT/SP), data = data.cast[data.cast$brn == "NO",])
anova(rich.v)
srv <- step(rich.v)
lsmrv<- lsmeansLT(rich.v)
#plot(rich.v)
#hist(resid(rich.v))

save(lsmra, lsmrs, lsmrv, file = "G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/lsmr.RData")


########################################################################################%%%%%
###############################  Environmental variables  ###################################
########################################################################################%%%%%


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
testInteractions(rich.a.dbh, custom = c(list(LUH = c(1,-1)), list(ET = c(0,1))), slope = "BA", adjustment="none") 
# R:E-P  0.003541, M:E-P 0.238, E:R-M 0.3057, P:R-M 0.1537
testInteractions(rich.a.dbh, custom = c(list(LUH = c(1,-1)), list(ET = c(0,1))), adjustment="none") 
# R:E-P 0.01406, M:E-P 0.1166, E:R-M 0.2309, P:R-M 0.1461
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

m1 <- glmer(pa ~ LUH*CT*scale(log(SeedMass_mg)) +ET +  (1|Site/LUH:CT) + (1|Site:LUH:CT:SP2) + (1 + LUH*CT|Spp) , nAGQ = 0, family = "binomial", data=sm.melt)# %>% select(-lnsm) %>% filter(Spp %in% as.list(sm.spplist %>% filter(Spp_n >= 8) %>% select(Spp))$Spp))
testInteractions(m1, custom = c(list(LUH = c(1,-1)), list(CT = c(0,1))), slope = "lnsm", adjustment="none") # R:T-U 0.0001138, M:T-U 0.09371, T:R-M  0.4, U:R-M 0.06195
testInteractions(m1, custom = c(list(LUH = c(1,-1)), list(CT = c(0,1))), adjustment="none") # R:T-U 4.454e-08, M:T-U 0.003106, T:R-M 0.6534, U:R-M 0.1393

summary(m1)
tab_model(m1)
sem.model.fits(m1)
hist(resid(m1))
qqnorm(residuals(m1))
plot(allEffects(m1try))


m2 <- glmer(pa ~ LUH + CT:LUH + CT*ET*scale(sla) + (1|Site/LUH:CT/SP2) + (1+ ET*CT|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=sla.melt)# %>% select(-lnsla) %>% filter(Spp %in% as.list(sla.spplist %>% filter(Spp_n >= 5) %>% select(Spp))$Spp))
Anova(m2)
rsquared(m2)
summary(m2)
testInteractions(m2, custom = c(list(CT = c(1, -1)), list(ET = c(0, 1))), slope = "lnsla", adjustment="none") # T:E-P 0.827, U:E-P 0.01952, E:T-U 0.1435, P:T-U 0.004004
testInteractions(m2, custom = c(list(CT = c(1, -1)), list(ET = c(0, 1))), adjustment="none") # T:E-P 0.01, U:E-P 3.548e-10, E:T-U 0.000239, P:T-U 1.919e-09
hist(resid(m2))
qqnorm(residuals(m2))
plot(allEffects(m2))


m3 <- glmer(pa ~ LUH + LUH:CT + CT*ET*scale(log(MaxPhotosynthetic_cm)) + (1|Site/LUH:CT/SP2) + (1+ CT*ET|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=ht.melt)# %>% select(-lnht) %>% filter(Spp %in% as.list(ht.spplist %>% filter(Spp_n >= 5) %>% select(Spp))$Spp))
testInteractions(m3, custom = c(list(CT = c(1, -1)), list(ET = c(0,1))), slope = "lnht", adjustment="none") # T:E-P 0.09361, U:E-P  0.2215, E:T-U 0.6927, P:T-U 0.06173
testInteractions(m3, custom = c(list(CT = c(1, -1)), list(ET = c(0,1))), adjustment="none") # T:E-P 0.003427, U:E-P 8.203e-08, E:T-U 0.0003088, P:T-U 2.437e-07
Anova(m3)
summary(m3)
sem.model.fits(m3)
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

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/data.cast.RData")


library(vegan)
betamatrixRTE <- data.cast[data.cast$LUH == "R" & data.cast$CT == "T" & data.cast$ET == "E" ,9:60]
specRTE <- specaccum(betamatrixRTE, method = "coleman", permutations = 10000)

betamatrixRTP <- data.cast[data.cast$LUH == "R" & data.cast$CT == "T" & data.cast$ET == "P" ,9:60]
specRTP <- specaccum(betamatrixRTP, method = "coleman", permutations = 10000)

betamatrixMTE <- data.cast[data.cast$LUH == "M" & data.cast$CT == "T" & data.cast$ET == "E" ,9:60]
specMTE <- specaccum(betamatrixMTE, method = "coleman", permutations = 10000)

betamatrixMTP <- data.cast[data.cast$LUH == "M" & data.cast$CT == "T" & data.cast$ET == "P" ,9:60]
specMTP <- specaccum(betamatrixMTP, method = "coleman", permutations = 10000)

betamatrixRUE <- data.cast[data.cast$LUH == "R" & data.cast$CT == "U" & data.cast$ET == "E" ,9:60]
specRUE <- specaccum(betamatrixRUE, method = "coleman", permutations = 10000)

betamatrixRUP <- data.cast[data.cast$LUH == "R" & data.cast$CT == "U" & data.cast$ET == "P" ,9:60]
specRUP <- specaccum(betamatrixRUP, method = "coleman", permutations = 10000)

betamatrixMUE <- data.cast[data.cast$LUH == "M" & data.cast$CT == "U" & data.cast$ET == "E" ,9:60]
specMUE <- specaccum(betamatrixMUE, method = "coleman", permutations = 10000)

betamatrixMUP <- data.cast[data.cast$LUH == "M" & data.cast$CT == "U" & data.cast$ET == "P" ,9:60]
specMUP <- specaccum(betamatrixMUP, method = "coleman", permutations = 10000)

plot(specRTE, ci = 0, col = "darkolivegreen3", lwd = 3, xlab = "Samples", ylab = "Number of species") #gold2, darkolivegreen
plot(specRTP, add = TRUE, ci = 0, col = "darkolivegreen3", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specMTE, add = TRUE, ci = 0, col = "gold2", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specMTP, add = TRUE, ci = 0, col = "gold2", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specRUE, add = TRUE, ci = 0, col = "darkolivegreen", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specRUP, add = TRUE, ci = 0, col = "darkolivegreen", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specMUE, add = TRUE, ci = 0, col = "gold4", lwd = 3, xlab = "Samples", ylab = "Number of species")
plot(specMUP, add = TRUE, ci = 0, col = "gold4", lty = "dashed", lwd = 3, xlab = "Samples", ylab = "Number of species")

data.cast.nms <- data.cast 
rownames(data.cast.nms) <- data.cast$ID
data.cast.nms <- data.cast.nms %>%
  select(ACAGRA:VITROT) %>%
  mutate(FAKE <- 1)

seedling_nms <- metaMDS(data.cast.nms, k = 2)

data.scores <- as.data.frame(scores(seedling_nms))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
 # create a column of site names, from the rownames of data.scores
data.scores <- cbind(data.scores,  data.cast %>% select(Site, LUH, CT, SP, ET, EP))


species.scores <- as.data.frame(scores(seedling_nms, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=interaction(CT,ET),colour=interaction(CT,ET)),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal()

  
#Volunteer only


VbetamatrixRTE <- data.cast[data.cast$LUH == "R" & data.cast$CT == "T" & data.cast$ET == "E" ,names(data.cast) %in% volunteerlist]
VspecRTE <- specaccum(VbetamatrixRTE, permutations = 10000)

VbetamatrixRTP <- data.cast[data.cast$LUH == "R" & data.cast$CT == "T" & data.cast$ET == "P" ,names(data.cast) %in% volunteerlist]
VspecRTP <- specaccum(VbetamatrixRTP, permutations = 10000)

VbetamatrixMTE <- data.cast[data.cast$LUH == "M" & data.cast$CT == "T" & data.cast$ET == "E" ,names(data.cast) %in% volunteerlist]
VspecMTE <- specaccum(VbetamatrixMTE, permutations = 10000)

VbetamatrixMTP <- data.cast[data.cast$LUH == "M" & data.cast$CT == "T" & data.cast$ET == "P" ,names(data.cast) %in% volunteerlist]
VspecMTP <- specaccum(VbetamatrixMTP, permutations = 10000)

VbetamatrixRUE <- data.cast[data.cast$LUH == "R" & data.cast$CT == "U" & data.cast$ET == "E" ,names(data.cast) %in% volunteerlist]
VspecRUE <- specaccum(VbetamatrixRUE, permutations = 10000)

VbetamatrixRUP <- data.cast[data.cast$LUH == "R" & data.cast$CT == "U" & data.cast$ET == "P" ,names(data.cast) %in% volunteerlist]
VspecRUP <- specaccum(VbetamatrixRUP, permutations = 10000)

VbetamatrixMUE <- data.cast[data.cast$LUH == "M" & data.cast$CT == "U" & data.cast$ET == "E" ,names(data.cast) %in% volunteerlist]
VspecMUE <- specaccum(VbetamatrixMUE, permutations = 10000)

VbetamatrixMUP <- data.cast[data.cast$LUH == "M" & data.cast$CT == "U" & data.cast$ET == "P" ,names(data.cast) %in% volunteerlist]
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

