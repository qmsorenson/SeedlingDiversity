setwd("C:/Users/Quinn Sorenson/Documents/R/Git_R/Root_exclosures/data")

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


data <- read.csv("C:/Users/Quinn Sorenson/Documents/R/Git_R/Root_exclosures/data/Initial survival and height.csv")
data$Sdling <- ifelse(data$SppCode == "NONE", 0, data$Sdling)
data <- data %>% 
  rename("Site" = "S_ID", "Spp" = "SppCode", "brn" = "Burned.")
data$ID <- paste(data$Site, data$LUH, data$CT, data$SP, data$ET, data$EP)


#### Seedling richness ####

data.cast <- dcast(data, Site + LUH + CT + SP + EP + ET + ID + brn ~ Spp, value.var = "Sdling", fun.aggregate = mean)
data.cast[is.na(data.cast)] <- 0
data.cast.check <- data.cast
data.cast$rich.a <- 0
for  (i in 1:length(data.cast$rich.a))
{
  data.cast[i, which(colnames(data.cast) == "rich.a")] <- sum(ifelse(data.cast[i,c(which(colnames(data.cast) == "ACAGRA"):which(colnames(data.cast) == "VITROT"))] == 0, 0, 1))
}
data.cast$rich.s <- 0
for  (i in 1:length(data.cast$rich.s))
{
  data.cast[i,which(colnames(data.cast) == "rich.s")] <- sum(ifelse(data.cast[i, names(data.cast) %in% c("EUPCUN", "VERANG", "CORMAJ", "SILCOM") ] == 0, 0, 1))
}
volunteerlist <- colnames(data.cast[i,c(which(colnames(data.cast) == "ACAGRA"):which(colnames(data.cast) == "VITROT"))])
volunteerlist <- volunteerlist[volunteerlist != "CORMAJ" & volunteerlist != "EUPCUN" & volunteerlist != "SILCOM" & volunteerlist != "VERANG"]
data.cast$rich.v <- 0
for  (i in 1:length(data.cast$rich.v))
{
  data.cast[i,which(colnames(data.cast) == "rich.v")] <- sum(ifelse(data.cast[i,names(data.cast) %in% volunteerlist] == 0, 0, 1))
}

data.cast$logrich.a <- log(data.cast$rich.a + 1)
data.cast$logrich.s <- log(data.cast$rich.s + 1)
data.cast$logrich.v <- log(data.cast$rich.v + 1)

#save(data.cast, file = "G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/data.cast.RData")

#### Environmental variables ####


# Water holding capacity

whc.id <- read.csv("WHC_ID.csv")
whc.ww <- read.csv("WHC_WW.csv")
whc.dw <- read.csv("WHC_DW.csv")
whc.tin <- read.csv("WHC_tin.csv")

whc.id$tin.id <- paste(whc.id$tin,whc.id$round)
whc.ww$tin.id <- paste(whc.ww$tin,whc.ww$round)
whc.dw$tin.id <- paste(whc.dw$tin,whc.dw$round)

whc <- left_join(whc.id, whc.ww[,c(2,4)], by = "tin.id")
whc <- left_join(whc, whc.dw[,c(2,6)], by = "tin.id")
whc <- left_join(whc, whc.tin, by = "tin")
tin.ave <- mean(whc$mass, na.rm = TRUE)
whc[c("mass")][is.na(whc[c("mass")])] <- tin.ave
whc <- whc[!is.na(whc$wet),]
whc <- whc[!is.na(whc$dry),]
whc$ID <- paste(whc$Site, whc$LUH, whc$CT, whc$SP, whc$ET, whc$EP)

whc$whc <- ((whc$wet-whc$mass)-(whc$dry-whc$mass))/(whc$dry-whc$mass)

data.cast.whc <- left_join(data.cast[data.cast$brn == "NO",], whc[,c(14,15)], by = "ID")
data.cast.whc$whc2 <- data.cast.whc$whc^2

# Light intensity

light<- read.csv("Light.csv")
names <- c("date", "time", "Site", "LUH", "CT", "SP", "EP", "ET", "angle", "X5cm", "X15cm", "X25cm", "notes")
colnames(light) <- names
light$ID <- paste(light$Site, light$LUH, light$CT, light$SP, light$ET, light$EP)
light2 <- aggregate(cbind(X5cm, X15cm, X25cm) ~ Site + LUH + CT + SP + EP + ET +ID, data = light, FUN = mean)
ear.light <- left_join(data.cast[data.cast$brn == "NO",], light2[,7:10], by = "ID")

# Compaction

comp <- read.csv("compaction.csv")
names <- c("Site", "LUH", "CT", "SP", "SS", "1", "2", "3", "4", "Notes")
colnames(comp) <- names
comp$Notes <- NULL
comp <- melt(data = comp, id = c("Site", "LUH", "CT", "SP", "SS"), variable.name = "EP")
comp$ID2 <- paste(comp$Site, comp$LUH, comp$CT, comp$SP, comp$EP)
library(dplyr)
names(comp)[names(comp)=="value"] <- "comp"
comp2 <- aggregate(comp ~ Site + LUH + CT + SP + EP + ID2, data = comp, FUN = mean)
data.cast$ID2 <- paste(data.cast$Site, data.cast$LUH, data.cast$CT, data.cast$SP, data.cast$EP)
data.castcomp <- left_join(data.cast[data.cast$brn == "NO",], comp2[,names(comp2) %in% c("comp", "ID2")], by = "ID2")

# Surrounding tree basal area

DBH <- read.csv("TreeDBH.csv")
DBH1 <- DBH[DBH$Spp != "",]
DBH1 <- aggregate(BA ~ Site + LUH + CT + SP, data = DBH1, FUN = sum)
DBH1$ID2 <- paste(DBH1$Site, DBH1$LUH, DBH1$CT, DBH1$SP)
data.cast.dbh <- data.cast[data.cast$brn == "NO",]
data.cast.dbh$ID2 <- paste(data.cast.dbh$Site, data.cast.dbh$LUH, data.cast.dbh$CT, data.cast.dbh$SP)
data.cast.dbh <- left_join(data.cast.dbh, DBH1[,5:6], by = "ID2")
data.cast.dbh$basal <- data.cast.dbh$BA
data.cast.dbh$BA <- scale(data.cast.dbh$BA)
data.cast.dbh$BA <- as.numeric(data.cast.dbh$BA)

# Surrounding vegetation cover

veg <- read.csv("Vegetation Survey.csv")
names <- c("order", "date", "time", "obs", "plot", "Site", "LUH", "CT", "SP", "DT", "ET", "EP", "Spp", "cover", "green", "leaf", "wood", "bare", "nonvasc")
colnames(veg) <- names
gre <- veg[veg$DT == "A",]
gre$ID <- paste(gre$Site, gre$LUH, gre$CT, gre$SP, gre$ET, gre$EP)
gre$ltb <- gre$leaf/(gre$leaf + gre$bare)
ear.gre <- left_join(data.cast[data.cast$brn == "NO",], gre[,names(gre) %in% c("ID", "green", "ltb")], by = "ID")

# combine together for cor table

data.cast.all <- data.cast.whc %>%
  select(Site, LUH, CT, ET, EP, SP, ID, whc) %>%
  left_join(ear.light %>% select(X5cm, X15cm, X25cm, ID), by = "ID") %>%
  left_join(data.cast.dbh %>% select(basal, ID), by = "ID") %>%
  left_join(ear.gre %>% select(green, ID), by = "ID") %>%
  left_join(data.castcomp %>% select(comp, ID2), by = c("ID" = "ID2"))

#save(data.cast.whc, ear.light, data.castcomp, data.cast.dbh, ear.gre, data.cast.all, file = "G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/Enviro.data.RData")


#### MLM functional traits ####

setwd("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data")
sd <- read.csv("BatchedDispSeedTraits.csv") #load seed traits
sd$SppCode <- substr(sd$SppCode, 1, 6)
#sd$SppCode <- ifelse(substr(sd$SppCode, 1,3) == "AND", "ANDSPP", sd$SppCode)
sd$SppCode <- ifelse(sd$SppCode == "LIASEC" | sd$SppCode == "LIAGRA", "LIASPP", sd$SppCode)
#sd$SppCode <- ifelse(sd$SppCode == "DICLEU", "DICSPP", sd$SppCode)
lw <- read.csv("QuantitativeTraits.csv") #load leaf weight
la <- read.csv("TotalLeafArea14-16.csv") #load leaf area
ht <- read.csv("CorridorProject_FunTraits_Height_forQuinn_04May2018.csv")
sl <- as.data.frame(colnames(data.cast[i,c(which(colnames(data.cast) == "ACAGRA"):which(colnames(data.cast) == "VITROT"))]))
colnames(sl)[1] <- "Spp"

sm <- left_join(sl, sd, by = c("Spp"="SppCode"))
sm$SeedMass_mg <- ifelse(sm$DispersuleSameAsSeed == TRUE, sm$DispersuleMass_mg, sm$SeedMass_mg) 
sm <- sm[!is.na(sm$SeedMass_mg),]
sm <- aggregate(SeedMass_mg ~ Spp, data = sm, FUN = mean)

sla <- left_join(la, lw, by= "CollectionName")
sla$sla <- sla$Total.Area/sla$DryWgt_g
sla <- sla[!(sla$PhotoName1 == "EREMIE1a*"),]
sla2 <- aggregate(sla~ SppCode, data = sla, FUN = mean)
LIA <- as.data.frame(cbind("LIASPP", mean(c(sla2[sla2$SppCode == "LIASEC", colnames(sla2) %in% c("sla")], sla2[sla2$SppCode == "LIAGRA", colnames(sla2) %in% c("sla")]))))
colnames(LIA) <- c("SppCode", "sla")
LIA$sla <- as.numeric(as.character(LIA$sla))
sla2 <- rbind(sla2, LIA)
sla2$SppCode <- as.character(sla2$SppCode)
sla2$SppCode <- ifelse(substr(sla2$SppCode, 1,6) == "DICLEU", "DICLEU", sla2$SppCode)
sla2 <- aggregate(sla~ SppCode, data = sla2, FUN = mean)

ht$SppCode <- as.character(ht$SppCode)
ht <- aggregate(cbind(MaxPhotosynthetic_cm, MaxReproductive_cm, MinReproductive_cm) ~ SppCode + Population, data = ht, FUN = mean)
ht <- aggregate(cbind(MaxPhotosynthetic_cm, MaxReproductive_cm, MinReproductive_cm) ~ SppCode, data = ht, FUN = mean)
ht$SppCode <- ifelse(substr(ht$SppCode, 1,6) == "DICLEU", "DICLEU", ht$SppCode)
ht <- aggregate(cbind(MaxPhotosynthetic_cm, MaxReproductive_cm, MinReproductive_cm) ~ SppCode, data = ht, FUN = mean)
LIAht <- as.data.frame(cbind("LIASPP", mean(c(ht[ht$SppCode == "LIASEC", colnames(ht) %in% c("MaxPhotosynthetic_cm")], ht[ht$SppCode == "LIAGRA", colnames(ht) %in% c("MaxPhotosynthetic_cm")])), mean(c(ht[ht$SppCode == "LIASEC", colnames(ht) %in% c("MaxReproductive_cm")], ht[ht$SppCode == "LIAGRA", colnames(ht) %in% c("MaxReproductive_cm")])), mean(c(ht[ht$SppCode == "LIASEC", colnames(ht) %in% c("MinReproductive_cm")], ht[ht$SppCode == "LIAGRA", colnames(ht) %in% c("MinReproductive_cm")]))))
colnames(LIAht) <- c("SppCode", "MaxPhotosynthetic_cm", "MaxReproductive_cm", "MinReproductive_cm")
LIAht$MaxPhotosynthetic_cm <- as.numeric(as.character(LIAht$MaxPhotosynthetic_cm))
LIAht$MaxReproductive_cm <- as.numeric(as.character(LIAht$MaxReproductive_cm))
LIAht$MinReproductive_cm <- as.numeric(as.character(LIAht$MinReproductive_cm))
ht <- rbind(ht, LIAht)
ht$SppCode <- ifelse(ht$SppCode == "PINPAL", "PINTAE", ht$SppCode)


datat.cast <- dcast(data, Site + LUH + CT + SP + EP + ET + ID + brn ~ Spp, value.var = "Sdling", fun.aggregate = mean)
datat.cast[is.na(datat.cast)] <- 0
datat.melt <- melt(data = datat.cast, id = c("Site" , "LUH" , "CT" , "SP" , "ET" , "EP" , "ID", "brn") , variable.name = "Spp")
datat.melt$brn <- as.character(datat.melt$brn)
datat.melt <- datat.melt[datat.melt$brn == "NO",]

sm.melt <- left_join(datat.melt, sm, by = "Spp")
sm.melt <- sm.melt[!is.na(sm.melt$SeedMass_mg),]

sm.melt$pa <- ifelse(sm.melt$value == 0, 0, 1)
sm.melt$obs <- 1:nrow(sm.melt)
sm.melt$lnsm <- scale(log(sm.melt$SeedMass_mg))
sm.melt$SP2 <- paste(sm.melt$LUH, sm.melt$CT, sm.melt$SP)


sla.melt <- left_join(datat.melt, sla2, by = c("Spp" = "SppCode"))
sla.melt <- sla.melt[!is.na(sla.melt$sla),]

sla.melt$pa <- ifelse(sla.melt$value == 0, 0, 1)
sla.melt$obs <- 1:nrow(sla.melt)

#sla.melt <- sla.melt[!sla.melt$sla < 50000,]
sla.melt$lnsla <- scale(log(sla.melt$sla))
sla.melt$SP2 <- paste(sla.melt$LUH, sla.melt$CT, sla.melt$SP)
sla.whc.melt <- left_join(sla.melt, data.cast.whc %>% select("ID", "whc"), by = "ID")
sla.whc.melt$lnsla <- scale(log(sla.whc.melt$sla))
sla.whc.melt$lnwhc <- scale(sla.whc.melt$whc)

ht.melt <- left_join(datat.melt, ht, by = c("Spp" = "SppCode"))
ht.melt <- ht.melt[!is.na(ht.melt$MaxPhotosynthetic_cm),]
ht.melt$pa <- ifelse(ht.melt$value == 0, 0, 1)
ht.melt$obs <- 1:nrow(ht.melt)
ht.melt$lnht <- scale(log(ht.melt$MaxPhotosynthetic_cm))
ht.melt$SP2 <- paste(ht.melt$LUH, ht.melt$CT, ht.melt$SP)

save(sm.melt, sla.melt, ht.melt, file = "R_dataframes/trait.dataframes.RData")


### check for unused species ###

fulll <- colnames(data.cast[,c(which(colnames(data.cast) == "ACAGRA"):which(colnames(data.cast) == "VITROT"))])
#redul <- unique(sm.melt$Spp)
#redul <- unique(sla.melt$Spp)
redul <- unique(ht.melt$Spp)
fulll[!(fulll %in% redul)]

### make histogram of number of seedlings in root exclosures

data.cast.check$tot.seed <- 0
for  (i in 1:nrow(data.cast.check))
{
  data.cast.check[i, which(colnames(data.cast.check) == "tot.seed")] <- sum(data.cast[i,c(which(colnames(data.cast) == "ACAGRA"):which(colnames(data.cast) == "VITROT"))])
}

nrow(data.cast.check[data.cast.check$tot.seed <= 20, ])
hist(data.cast.check$tot.seed, breaks = 188)

### write files for chtc ###

#write.csv(sm.melt, "C:/Users/Quinn Sorenson/Google Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/CHTC files/sm_BootMer0.csv")
#write.csv(sla.melt, "C:/Users/Quinn Sorenson/Google Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/CHTC files/sla_BootMer0.csv")
#write.csv(ht.melt, "C:/Users/Quinn Sorenson/Google Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/CHTC files/ht_BootMer0.csv")