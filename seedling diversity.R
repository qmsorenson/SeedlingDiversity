setwd("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/data") #Set working directory

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


data <- read.csv("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/data/Initial survival and height.csv")
data$Sdling <- ifelse(data$SppCode == "NONE", 0, data$Sdling)
library(plyr)
data <- plyr::rename(data, replace = c("S_ID" = "Site", "SppCode" = "Spp", "Burned." = "brn"))
data$ID <- paste(data$Site, data$LUH, data$CT, data$SP, data$ET, data$EP)
datat <- data


library(reshape2)
data.cast <- dcast(data, Site + LUH + CT + SP + EP + ET + ID + brn ~ Spp, value.var = "Sdling", fun.aggregate = mean)
data.cast[is.na(data.cast)] <- 0
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


data.cast$logCORMAJ <- log(data.cast$CORMAJ + 1)
data.cast$logEUPCUN <- log(data.cast$EUPCUN + 1)
data.cast$logSILCOM <- log(data.cast$SILCOM + 1)
data.cast$logVERANG <- log(data.cast$VERANG + 1)
#data.cast2 <- data.cast
#data.cast2[1:380, 7:93] <- ifelse(data.cast2[1:380, 7:93] == 0, 0, 1)
#data.sp.num <- summarise_each(data.cast2, funs(sum))
#data.sp.num[1,7:93]<-data.sp.num[1,7:93]*386
#data.sp.num <- t(data.sp.num)
data.cast2 <- dcast(data, Site + LUH + CT + SP + EP + ET + ID + brn ~ Spp, value.var = "remain", fun.aggregate = mean)

datat.cast <- dcast(data, Site + LUH + CT + SP + EP + ET + ID + brn ~ Spp, value.var = "Sdling", fun.aggregate = mean)
datat.cast[is.na(datat.cast)] <- 0
#datat.melt <- melt(data = datat.cast, id = c("Site" , "LUH" , "CT" , "SP" , "ET" , "EP" , "ID", "brn") , variable.name = "Spp")
#datat.melt <- left_join(datat.melt, sl, by = "Spp")
#datat.melt <- datat.melt[!is.na(datat.melt$SeedMass_mg),]

#datat.melt2 <- left_join(datat.melt, sla2, by = c("Spp" = "Spp.x"))
#datat.melt2 <- datat.melt2[!is.na(datat.melt2$sla),]


library(grid)
library(gridExtra)

fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

rotatedAxisElementText = function(angle,position='x',Size,Color){
  angle     = angle[1]; 
  position  = position[1]
  positions = list(x=0,y=90,top=180,right=270)
  if(!position %in% names(positions))
    stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
  if(!is.numeric(angle))
    stop("'angle' must be numeric",call.=FALSE)
  rads  = (angle - positions[[ position ]])*pi/180
  hjust = 0.5*(1 - sin(rads))
  vjust = 0.5*(1 + cos(rads))
  element_text(angle=angle,vjust=vjust,hjust=hjust, size = Size, color = Color)
}

options(na.action = "na.fail")
options(na.action = "na.omit")

library(dplyr)
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

ht.melt <- left_join(datat.melt, ht, by = c("Spp" = "SppCode"))
ht.melt <- ht.melt[!is.na(ht.melt$MaxPhotosynthetic_cm),]
ht.melt$pa <- ifelse(ht.melt$value == 0, 0, 1)
ht.melt$obs <- 1:nrow(ht.melt)
ht.melt$lnht <- scale(log(ht.melt$MaxPhotosynthetic_cm))
ht.melt$SP2 <- paste(ht.melt$LUH, ht.melt$CT, ht.melt$SP)


fulll <- colnames(data.cast[,c(which(colnames(data.cast) == "ACAGRA"):which(colnames(data.cast) == "VITROT"))])
#redul <- unique(sm.melt$Spp)
#redul <- unique(sla.melt$Spp)
redul <- unique(ht.melt$Spp)
fulll[!(fulll %in% redul)]

#write.csv(sm.melt, "C:/Users/Quinn Sorenson/Google Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/CHTC files/sm_BootMer0.csv")
#write.csv(sla.melt, "C:/Users/Quinn Sorenson/Google Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/CHTC files/sla_BootMer0.csv")
#write.csv(ht.melt, "C:/Users/Quinn Sorenson/Google Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/CHTC files/ht_BootMer0.csv")
#############################################################################################
###############################  richness without burned  ###################################
#############################################################################################

hist(data.cast$rich.a)
hist(data.cast$rich.s)
hist(data.cast$rich.v)
data.cast$logrich.a <- log(data.cast$rich.a + 1)
data.cast$logrich.s <- log(data.cast$rich.s + 1)
data.cast$logrich.v <- log(data.cast$rich.v + 1)
hist(data.cast$logrich.a)
hist(data.cast$logrich.s)
hist(data.cast$logrich.v)


rich.a <- lmer(logrich.a ~ ET*LUH*CT + (1|Site/LUH:CT/SP), data = data.cast[data.cast$brn == "NO",])
sra <- step(rich.a)
plot(rich.a)
hist(resid(rich.a))
lsmra<- lsmeansLT(rich.a)

rich.s <- lmer(logrich.s ~ ET*LUH*CT + (1|Site/LUH:CT/SP), data = data.cast[data.cast$brn == "NO",])
srs <- step(rich.s)
plot(rich.s)
hist(resid(rich.s))
lsmrs<- lsmeansLT(rich.s)

rich.v <- lmer(logrich.v ~ ET*LUH*CT + (1|Site/LUH:CT/SP), data = data.cast[data.cast$brn == "NO",])
srv <- step(rich.v)
plot(rich.v)
hist(resid(rich.v))
lsmrv<- lsmeansLT(rich.v)

### Plots generated 4.25.18 ####

### bar plots ###

pcas<- ggplot(data = lsmra$lsmeans.table[is.na(lsmra$lsmeans.table[,3]) & is.na(lsmra$lsmeans.table[,2]),], 
              aes(x = factor(ET, labels = c("E" = "Excluded", "P" = "Present")), y = exp(Estimate)-1, fill = factor(ET))) + 
  
  scale_fill_manual(values=c("Grey50", "white")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1), position=position_dodge(.9),  width = 0.1) + 
  labs(title = expression(paste("a) All species"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3 ), labels = fmt_dcimals(1)) +
  #ylab("") + 
  ylab(NULL) + 
  xlab(NULL) +
  #guides(fill=guide_legend(title="Belowground \ncompetition", keywidth=0.4,
  #                        keyheight=0.6,
  #                       default.unit="inch")) +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(),
        text = element_text(size = 30),
        #legend.spacing = unit(2,"cm"),
        #legend.key = element_rect(size = 10),
        #legend.key.size = unit(2, 'lines'),
        axis.text.x = element_text(size = 34, color = "black"),#, hjust = 1, vjust =.97),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")

pcss <- ggplot(data = lsmrs$lsmeans.table[is.na(lsmrs$lsmeans.table[,3]) & is.na(lsmrs$lsmeans.table[,2]),], 
               aes(x = factor(ET, labels = c("E" = "Excluded", "P" = "Present")), y = exp(Estimate)-1, fill = factor(ET))) + 
  scale_fill_manual(values=c("Grey50", "white")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1), position=position_dodge(.9),  width = 0.1) + 
  labs(title = expression(paste("b) Seeded species"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3 ), labels = fmt_dcimals(1)) +
  #ylab("") + 
  ylab(NULL) + 
  xlab(NULL) +
  #guides(fill=guide_legend(title="Belowground \ncompetition", keywidth=0.4,
  #                        keyheight=0.6,
  #                       default.unit="inch")) +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(),
        text = element_text(size = 30),
        #legend.spacing = unit(2,"cm"),
        #legend.key = element_rect(size = 10),
        #legend.key.size = unit(2, 'lines'),
        axis.text.x = element_text(size = 34, color = "black"),#, hjust = 1, vjust =.97),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")


pcvs <- ggplot(data = lsmrv$lsmeans.table[is.na(lsmrv$lsmeans.table[,3]) & is.na(lsmrv$lsmeans.table[,2]),], 
               aes(x = factor(ET, labels = c("E" = "Excluded", "P" = "Present")), y = exp(Estimate)-1, fill = factor(ET))) + 
  scale_fill_manual(values=c("Grey50", "white")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1), position=position_dodge(.9),  width = 0.1) + 
  labs(title = expression(paste("c) Volunteer species"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3 ), labels = fmt_dcimals(1)) +
  #ylab("") + 
  ylab(NULL) + 
 xlab(label = NULL) +#"Belowground competition") +
  #guides(fill=guide_legend(title="Belowground \ncompetition", keywidth=0.4,
  #                        keyheight=0.6,
  #                       default.unit="inch")) +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(),
        text = element_text(size = 30),
        #legend.spacing = unit(2,"cm"),
        #legend.key = element_rect(size = 10),
        #legend.key.size = unit(2, 'lines'),
        axis.text.x = element_text(size = 34, color = "black"),#, hjust = 1, vjust =.97),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")




### interaction plots

pias <- ggplot(data = lsmra$lsmeans.table[complete.cases(lsmra$lsmeans.table[,2:3]) & is.na(lsmra$lsmeans.table[,1]),], aes(x = factor(LUH), y = exp(Estimate)-1, colour = CT, shape = CT)) + 
  geom_point(data = lsmra$lsmeans.table[complete.cases(lsmra$lsmeans.table[,2:3]) & is.na(lsmra$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, fill = factor(LUH)), size = 10, show.legend = FALSE) +
  geom_line(data = lsmra$lsmeans.table[complete.cases(lsmra$lsmeans.table[,2:3]) & is.na(lsmra$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, group = CT, linetype = CT), show.legend = FALSE, size = 1.7 ) + 
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1, width = 0.1)) + 
  labs(title = "",#expression(paste("a) All species")),
       x = NULL, y = NULL, shape = "Tree density", 
       linetype = "Tree density") + 
  scale_x_discrete(labels = c("M" = "Post-ag", "R" = "Remnant")) +
  scale_y_continuous(labels = fmt_dcimals(1)) +
  scale_colour_manual(values=c("Grey60", "Grey20")) +
  scale_fill_manual(values=c("gold2", "darkolivegreen")) +
  scale_shape_manual(values=c(21, 21)) +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(), legend.title.align = .5,
        title =element_text(size=34),
        text = element_text(size = 30), axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")

piss<- ggplot(data = lsmrs$lsmeans.table[complete.cases(lsmrs$lsmeans.table[,2:3]) & is.na(lsmrs$lsmeans.table[,1]),], aes(x = factor(LUH), y = exp(Estimate)-1, colour = CT, shape = CT)) + 
  geom_point(data = lsmrs$lsmeans.table[complete.cases(lsmrs$lsmeans.table[,2:3]) & is.na(lsmrs$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, fill = factor(LUH)), size = 10, show.legend = FALSE) +
  geom_line(data = lsmrs$lsmeans.table[complete.cases(lsmrs$lsmeans.table[,2:3]) & is.na(lsmrs$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, group = CT, linetype = CT), show.legend = FALSE, size = 1.7 ) + 
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1, width = 0.1)) + 
  labs(title = " ",#expression(paste("b) Seeded species")),
       x = NULL, y = NULL, shape = "Tree density", 
       linetype = "Tree density") + 
  scale_x_discrete(labels = c("M" = "Post-ag", "R" = "Remnant")) +
  scale_y_continuous(labels = fmt_dcimals(0)) +
  scale_colour_manual(values=c("Grey60", "Grey20")) +
  scale_fill_manual(values=c("gold2", "darkolivegreen")) +
  scale_shape_manual(values=c(21, 21)) +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(), legend.title.align = .5,
        title =element_text(size=34),
        text = element_text(size = 30), axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")

pivs<- ggplot(data = lsmrv$lsmeans.table[complete.cases(lsmrv$lsmeans.table[,2:3]) & is.na(lsmrv$lsmeans.table[,1]),], aes(x = factor(LUH), y = exp(Estimate)-1, colour = CT, shape = CT)) + 
  geom_point(data = lsmrv$lsmeans.table[complete.cases(lsmrv$lsmeans.table[,2:3]) & is.na(lsmrv$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, fill = factor(LUH)), size = 10, show.legend = FALSE) +
  geom_line(data = lsmrv$lsmeans.table[complete.cases(lsmrv$lsmeans.table[,2:3]) & is.na(lsmrv$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, group = CT, linetype = CT), show.legend = FALSE, size = 1.7 ) + 
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1, width = 0.1)) + 
  labs(title = "",#expression(paste("c) Volunteer species")),
       x = NULL, y = NULL, shape = "Tree density", 
       linetype = "Tree density") + 
  scale_x_discrete(labels = c("M" = "Post-ag", "R" = "Remnant")) +
  #xlab(label = "Land-use History") +
  scale_y_continuous(labels = fmt_dcimals(0)) +
  scale_colour_manual(values=c("Grey60", "Grey20")) +
  scale_fill_manual(values=c("gold2", "darkolivegreen")) +
  scale_shape_manual(values=c(21, 21)) +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(), legend.title.align = .5,
        title =element_text(size=34),
        text = element_text(size = 30), axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")

grid.arrange(pcas, pias, pcss, piss, pcvs, pivs, ncol = 2, left = textGrob("Species richness", gp=gpar(fontsize=34), rot = 90), bottom = textGrob("Belowground Competition     Land-use history  ", gp=gpar(fontsize =34))) # pdf dim 10.75 16.5





means1 <- aggregate(cbind(logrich.a, logrich.s, logrich.v) ~ Site + ET + LUH + CT, data = data.cast[data.cast$brn == "NO",], FUN = mean)
means <- aggregate(cbind(logrich.a, logrich.s, logrich.v) ~ ET + CT + LUH, data = means1, FUN = mean)
means$mean.a <- exp(means$logrich.a) - 1
means$mean.s <- exp(means$logrich.s) - 1
means$mean.v <- exp(means$logrich.v) - 1
err.a <- aggregate(logrich.a ~ ET + CT + LUH, data = means1, FUN = st.err)
err.s <- aggregate(logrich.s ~ ET + CT + LUH, data = means1, FUN = st.err)
err.v <- aggregate(logrich.v ~ ET + CT + LUH, data = means1, FUN = st.err)
means$logerr.a <- err.a$logrich.a
means$logerr.s <- err.s$logrich.s
means$logerr.v <- err.v$logrich.v
means$up.a <- exp(means$logrich.a + means$logerr.a) - 1
means$lo.a <- exp(means$logrich.a - means$logerr.a) - 1
means$up.s <- exp(means$logrich.s + means$logerr.s) - 1
means$lo.s <- exp(means$logrich.s - means$logerr.s) - 1
means$up.v <- exp(means$logrich.v + means$logerr.v) - 1
means$lo.v <- exp(means$logrich.v - means$logerr.v) - 1


pa1 <- ggplot(data=lsmra$lsmeans.table[complete.cases(lsmra$lsmeans.table[,1:3]),], aes(x=LUH:CT, y=exp(Estimate), fill=factor(ET, labels = c("E" = "Without \nbelowground \ncompetition",
                                                                                                                                              "P" = "With \nbelowground \ncompetition")))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c("#999999", "#FFFFFF")) +
  geom_errorbar(aes(ymax = exp(Estimate+`Standard Error`), ymin=exp(Estimate-`Standard Error`)), 
                position=position_dodge(.9), width = 0.3) +
  scale_x_discrete(labels = c("M:T" = "Post-ag \nThinned", "M:U" = "Post-ag \nUnthinned", "R:T" = "Remnant \nThinned", "R:U" = "Remnant \nUnthinned")) +
  scale_y_continuous(expand = c(0,0), limit = c(0,5.5)) +
  ylab("Species richness \n(# species / exclosure)") +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none") +
  labs(title ="", fill ="")

ps2 <- ggplot(data=lsmrs$lsmeans.table[complete.cases(lsmrs$lsmeans.table[,1:3]),], aes(x=LUH:CT, y=exp(Estimate), fill=factor(ET, labels = c("E" = "Without \nbelowground \ncompetition",
                                                                                                                                              "P" = "With \nbelowground \ncompetition")))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c("#999999", "#FFFFFF")) +
  geom_errorbar(aes(ymax = exp(Estimate+`Standard Error`), ymin=exp(Estimate-`Standard Error`)), 
                position=position_dodge(.9), width = 0.3) +
  scale_x_discrete(labels = c("M:T" = "Post-ag \nThinned", "M:U" = "Post-ag \nUnthinned", "R:T" = "Remnant \nThinned", "R:U" = "Remnant \nUnthinned")) +
  scale_y_continuous(expand = c(0,0), limit = c(0,3.75)) +
  ylab("Seeded species richness \n(# species / exclosure)") +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none") +
  labs(title ="", fill ="")

pv3 <- ggplot(data=lsmrv$lsmeans.table[complete.cases(lsmrv$lsmeans.table[,1:3]),], aes(x=LUH:CT, y=exp(Estimate), fill=factor(ET, labels = c("E" = "Without \nbelowground \ncompetition",
                                                                                                                                              "P" = "With \nbelowground \ncompetition")))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c("#999999", "#FFFFFF")) +
  geom_errorbar(aes(ymax = exp(Estimate+`Standard Error`), ymin=exp(Estimate-`Standard Error`)), 
                position=position_dodge(.9), width = 0.3) +
  scale_x_discrete(labels = c("M:T" = "Post-ag \nThinned", "M:U" = "Post-ag \nUnthinned", "R:T" = "Remnant \nThinned", "R:U" = "Remnant \nUnthinned")) +
  scale_y_continuous(expand = c(0,0), limit = c(0,2.75)) +
  ylab("Volunteer species richness \n(# species / exclosure)") +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none") +
  labs(title ="", fill ="")

grid.arrange(pa1, ps2, pv3, ncol = 1, bottom = textGrob("Land-use history and Canopy Thinning", gp=gpar(fontsize=12)))

### VOLUNTEER V SEEDED ####

pda <- data.cast[data.cast$brn == "NO",]

ggplot(pda, aes(x = logrich.s, y = logrich.v)) +
  geom_count(mapping = NULL, data = NULL, stat = "sum",
             position = "identity", na.rm = FALSE, show.legend = NA,
             inherit.aes = TRUE) +
  geom_smooth( method = "lm", se = FALSE, color = "black") + 
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))


##### attempt at by species bar ######

data.cast$logCORMAJ <- ifelse(data.cast$CORMAJ == 0, 0, 1)
data.cast$logEUPCUN <- ifelse(data.cast$EUPCUN == 0, 0, 1)
data.cast$logSILCOM <- ifelse(data.cast$SILCOM == 0, 0, 1)
data.cast$logVERANG <- ifelse(data.cast$VERANG == 0, 0, 1)

means1 <- aggregate(cbind(logCORMAJ, logVERANG, logEUPCUN, logSILCOM) ~ Site + ET + LUH + CT, data = data.cast[data.cast$brn == "NO",], FUN = mean)
means <- aggregate(cbind(logCORMAJ, logVERANG, logEUPCUN, logSILCOM) ~ ET + CT + LUH, data = means1, FUN = mean)
means$total <- means$logCORMAJ + means$logVERANG + means$logEUPCUN + means$logSILCOM
means$pC <- means$logCORMAJ/means$total
means$pV <- means$logVERANG/means$total
means$pE <- means$logEUPCUN/means$total
means$pS <- means$logSILCOM/means$total

######### Environmental variables ############

# must make these first (below): data.cast.whc ear.light data.castcomp data.cast.dbh ear.gre

###### all seedlings ######
data.cast.whc$SP2 <- paste(data.cast.whc$LUH, data.cast.whc$CT, data.cast.whc$SP)
ear.light$SP2 <- paste(ear.light$LUH, ear.light$CT, ear.light$SP)
data.castcomp$SP2 <- paste(data.castcomp$LUH, data.castcomp$CT, data.castcomp$SP)
data.cast.dbh$SP2 <- paste(data.cast.dbh$LUH, data.cast.dbh$CT, data.cast.dbh$SP)


data.cast.whc <- data.cast.whc[!is.na(data.cast.whc$whc),]
#data.cast.whc$obs <- 1:nrow(data.cast.whc)
#data.cast.whc$obs <- as.character(data.cast.whc$obs)
#glmtry <- glmer(rich.a ~ ET*whc*CT*LUH + (1|Site/LUH:CT/SP) + (1|obs), data = data.cast.whc, family = "poisson")
rich.a.whc <- lmer(logrich.a ~ ET*LUH*CT*whc +  (1|Site/LUH:CT/SP), data = data.cast.whc)
pract <- lmer(logrich.a ~ ET*whc + (1|Site/LUH:CT/SP), data = data.cast.whc)
sraw <- step(rich.a.whc)
plot(rich.a.whc)
hist(resid(rich.a.whc))
data.cast.whc$predicted <- predict(pract)

rich.a.lig <- lmer(logrich.a ~ ET*LUH*CT*X5cm + (1|Site/LUH:CT/SP2), data = ear.light)
sral <- step(rich.a.lig)
plot(rich.a.lig)
hist(resid(rich.a.lig))

rich.a.comp<- lmer(logrich.a ~ ET*LUH*comp*CT + (1|Site/SP2), data = data.castcomp)
srac <- step(rich.a.comp)
plot(rich.a.comp)
hist(resid(rich.a.comp))

rich.a.dbh <- lmer(logrich.a ~ ET*LUH*BA +(1|Site/LUH:CT/SP), data = data.cast.dbh)
anova(rich.a.dbh)
srad <- step(rich.a.dbh)
plot(rich.a.dbh)
hist(resid(rich.a.dbh))

rich.a.dbh <- lmer(logrich.a ~ ET*CT*LUH*BA + (1|Site/LUH:CT/SP), data = data.cast.dbh)
anova(rich.a.dbh)
srad <- step(rich.a.dbh)
plot(rich.a.dbh)
hist(resid(rich.a.dbh))

rich.a.gre <- lmer(logrich.a ~ ET*LUH*CT*green + (1|Site/LUH:CT/SP), data = ear.gre)
srag <- step(rich.a.gre)
plot(rich.a.gre)
hist(resid(rich.a.gre))

############################################################################################
#########################################   SEM   ##########################################
############################################################################################

library(piecewiseSEM)
semrichall <- read.csv("data.cast.all.csv")
library(nlme)
semrichall$SP2 <- paste(semrichall$LUH, semrichall$CT, semrichall$SP)

#mod.1 <- list( 
 # lme(logrich.a ~  LUH*CT*ET + comp + whc + BA + BA:CT + green +X5cm , random = ~1|Site/SP2, na.action = na.omit, data = semrichall), #logrich.a ~  LUH*CT + ET + LUH:ET + BA:LUH:ET +LUH:BA + BA:ET + comp + whc + ET:whc + green + CT:X5cm + X5cm 
  #lme(comp ~ LUH*CT, random = ~1|Site/SP2, na.action = na.omit, data = semrichall),
#  lme(whc ~ LUH*CT + BA + BA:CT , random = ~1|Site/SP2, na.action = na.omit, data = semrichall),
 # lme(BA ~ LUH*CT*whc, random = ~1|Site, na.action = na.omit, data = semrichall),
  #lme(green ~ whc +comp + BA, random = ~1|Site/SP2, na.action = na.omit, data = semrichall),
  #lme(X5cm ~ green + BA +LUH*CT, random = ~1|Site/SP2, na.action = na.omit, data = semrichall))







mod.1 <- list( 
  lme(logrich.a ~ LUH*CT + comp + ET  + whc+ green  + X5cm +BA , random = ~1|Site/SP2, na.action = na.omit, data = semrichall), #logrich.a ~  LUH*CT + ET + LUH:ET + comp + whc + ET:whc + BA:LUH:ET +LUH:BA + BA:ET + green + CT:X5cm + X5cm 
  lme(comp ~ LUH*CT, random = ~1|Site/SP2, na.action = na.omit, data = semrichall),
  lme(whc ~ LUH*CT, random = ~1|Site/SP2, na.action = na.omit, data = semrichall),
  lme(BA ~ LUH*CT + whc, random = ~1|Site, na.action = na.omit, data = semrichall),
  lme(green ~ LUH*CT, random = ~1|Site/SP2, na.action = na.omit, data = semrichall),
  lme(X5cm ~ green + BA + CT+LUH, random = ~1|Site/SP2, na.action = na.omit, data = semrichall))

sem.fit(mod.1, semrichall)
sem.coefs(mod.1, coreop)

######## Plots ###############



#wch
pw1 <- ggplot(data.cast.whc, 
              aes(x= whc, y = logrich.a,linetype = factor(ET), shape = ET)) + 
  geom_point() + labs(colour = "Legend", y = NULL,x = "Soil water holding capacity (%)") +
  geom_smooth( method = "lm", se = FALSE, colour = "black") + 
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  scale_shape_manual(values=c(19, 4), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"), legend.position="none")

#light
pl2 <- ggplot(ear.light, 
              aes(x= X5cm, y = logrich.a, colour =  CT)) + 
  geom_point() + labs(colour = "Legend", y = NULL,x = expression(paste("Light intensity (", mu, "mol/s)"))) +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_colour_manual(values=c("Grey60","royalblue3")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"), legend.position="none")

#dbh
pd3 <- ggplot(data.cast.dbh, 
              aes(x= basal, y = logrich.a, colour =  factor(LUH),linetype = factor(ET), shape = ET)) + 
  geom_point() + labs(colour = "Legend", y = NULL,x = expression(paste("Basal area (cm"^"2", ")"))) +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_colour_manual(values=c("Grey20","coral3")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  scale_shape_manual(values=c(19, 4), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"), legend.position="none")

#gre
pg4 <- ggplot(ear.gre, 
              aes(x= green, y = logrich.a)) + 
  geom_point() + labs(colour = "Legend", y = NULL,x = "Vegetation cover (%)") +
  geom_smooth(method = "lm", se = FALSE, colour="black") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"), legend.position="none")

grid.arrange(pw1, pl2, pd3, pg4, ncol = 2, left = textGrob("Species richness \nln(# species / exclosure)", gp=gpar(fontsize=14), rot = 90)) 


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

rich.s.dbh <- lmer(logrich.s ~ ET*LUH*CT*BA - ET:LUH:CT:BA - ET:CT:BA - LUH:CT:BA - ET:LUH:CT - ET:LUH:BA - CT:BA- LUH:BA - ET:LUH-ET:BA- ET:CT- BA+ (1|Site/LUH:CT/SP), data = data.cast.dbh)
anova(rich.s.dbh)
srsd <- step(rich.s.dbh)
plot(rich.s.dbh)
hist(resid(rich.s.dbh))

rich.s.gre <- lmer(logrich.s ~ ET*LUH*CT*green + (1|Site/LUH:CT/SP), data = ear.gre)
srsg <- step(rich.s.gre)
plot(rich.s.gre)
hist(resid(rich.s.gre))


######## Plots ###############

#wch
ggplot(data.cast.whc, 
       aes(x= whc, y = logrich.s,linetype = factor(ET))) + 
  geom_point() + labs(title = "All richness", colour = "Legend", y = "ln # of spp", x = "Soil water holding capacity (%)") +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#light
ggplot(ear.light, 
       aes(x= X5cm, y = logrich.s, colour =  factor(CT))) + 
  geom_point() + labs(title = "All richness", colour = "Legend", y = "ln # of spp", x = "PPFD") +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_colour_manual(values=c("Grey60","blue")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#dbh
ggplot(data.cast.dbh, 
       aes(x= BA, y = logrich.s, colour =  factor(CT),linetype = factor(ET))) + 
  geom_point() + labs(title = "All richness", colour = "Legend", y = "ln # of spp", x = "DBH") +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_colour_manual(values=c("black","red")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#gre
ggplot(ear.gre, 
       aes(x= green, y = logrich.s, color = LUH:CT)) + 
  geom_point() + labs(title = "All richness", colour = "Legend", y = "ln # of spp", x = "% vegetation surrounding") +
  geom_smooth( method = "lm", se = FALSE) + 
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

########### volunteer species ############

###### all seedlings ######

rich.v.whc <- lmer(logrich.v ~ CT +ET*whc + (1|Site/LUH:CT/SP), data = data.cast.whc)
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

rich.v.dbh <- lmer(logrich.v ~ ET*LUH*CT*BA - ET:LUH:CT:BA - ET:CT:BA - ET:LUH:CT - LUH:CT:BA - CT:BA - ET:CT - LUH:CT + (1|Site/LUH:CT/SP), data = data.cast.dbh)
anova(rich.v.dbh)
srvd <- step(rich.v.dbh)
plot(rich.v.dbh)
hist(resid(rich.v.dbh))

rich.v.gre <- lmer(logrich.v ~ ET*LUH*CT*green + (1|Site/LUH:CT/SP), data = ear.gre)
srvg <- step(rich.v.gre)
plot(rich.v.gre)
hist(resid(rich.v.gre))


######## Plots ###############

#wch
ggplot(data.cast.whc, 
       aes(x= whc, y = logrich.v,linetype = factor(ET))) + 
  geom_point() + labs(title = "All richness", colour = "Legend", y = "ln # of spp", x = "Soil water holding capacity (%)") +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#light
ggplot(ear.light, 
       aes(x= X5cm, y = logrich.v, colour =  paste(factor(CT), factor(LUH)), linetype = factor(ET))) + 
  geom_point() + labs(title = "All richness", colour = "Legend", y = "ln # of spp", x = "PPFD") +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_colour_manual(values=c("black", "red","Grey60","blue")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#dbh
ggplot(data.cast.dbh, 
       aes(x= BA, y = logrich.v, colour =  factor(CT),linetype = factor(ET))) + 
  geom_point() + labs(title = "All richness", colour = "Legend", y = "ln # of spp", x = "DBH") +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_colour_manual(values=c("black","red")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#gre
ggplot(ear.gre, 
       aes(x= green, y = logrich.v)) + 
  geom_point() + labs(title = "All richness", colour = "Legend", y = "ln # of spp", x = "% vegetation surrounding") +
  geom_smooth( method = "lm", se = FALSE) + 
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#############################################################################################

############################################################################################
#################################   MLM TRAITS #############################################
#############################################################################################

library("lme4")

#binoomial models

m1 <- glmer(pa ~ LUH*CT*lnsm + ET + (1|Site/LUH:CT/SP2) + (1+ LUH*CT|Spp) + (1|obs), nAGQ = 0, family = "binomial", data=sm.melt)

summary(m1)
hist(resid(m1))
qqnorm(residuals(m1))
plot(allEffects(m1))




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

####Plots
m1sm <- glmer(pa ~ LUH*CT*lnsm + (1|Site/LUH:CT/SP2) + (1+ LUH*CT|Spp) + (1|obs), nAGQ = 0, family = "binomial", data=sm.melt)


m2sla <- glmer(pa ~ CT*ET*lnsla + (1|Site/LUH:CT/SP2) + (1+ ET*CT|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=sla.melt)

m3ht <- glmer(pa ~ CT*ET*lnht + (1|Site/LUH:CT/SP2) + (1+ CT*ET|Spp) + (1|obs), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=ht.melt)

#### Seed mass ####



predict_sm.melt <- sm.melt %>%
  expand(nesting(LUH,CT),lnsm)

predict_sm.melt$pa <- predict(m1sm, newdata = predict_sm.melt, re.form = NA, type="response")
#predict_sm.melt_no_et <- aggregate(pa ~ LUH + CT + lnsm, data = predict_sm.melt, FUN = mean)

psm <- ggplot(sm.melt, aes(x = lnsm, y = pa, color = interaction(LUH,CT))) +
  labs(#title = expression(italic("Coreopsis major")), 
    #colour = "Legend", 
    y = NULL,#"Log total biomass (mg)", 
    x = "Seed mass") +
    geom_line(data = predict_sm.melt, size = 2) +
  scale_colour_manual(values=c("gold2","darkolivegreen3", "gold4", "darkolivegreen")) +
  #scale_linetype_manual(values = c("solid", "twodash"), name = "Legend") +
  scale_y_continuous(breaks = c(.01, .03)) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 24), axis.text.x = element_text(size = 24, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 24,  Color = "black"), 
        legend.position="none")



#### Specific leaf area ####

predict_sla.melt <- sla.melt %>%
  expand(nesting(CT,ET),lnsla)

predict_sla.melt$pa <- predict(m2sla, newdata = predict_sla.melt, re.form = NA, type="response")
#predict_sla.melt_no_luh <- aggregate(pa ~ ET + CT + lnsla, data = predict_sla.melt, FUN = mean)

psla <- ggplot(ht.melt, aes(x = lnsla, y = pa, color = CT, linetype = ET)) +
  labs(#title = expression(italic("Coreopsis major")), 
    #colour = "Legend", 
    y = NULL,#"Log total biomass (mg)", 
    x = "Specific leaf area") +
  geom_line(data = predict_sla.melt, size = 2) +
  scale_colour_manual(values=c("khaki3","khaki4")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  scale_y_continuous(breaks = c(.01, .03)) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 24), axis.text.x = element_text(size = 24, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 24,  Color = "black"), 
        legend.position="none")


#### Max ht ####

predict_ht.melt <- ht.melt %>%
  expand(nesting(CT,ET),lnht)

predict_ht.melt$pa <- predict(m3ht, newdata = predict_ht.melt, re.form = NA, type="response")
#predict_ht.melt_no_luh <- aggregate(pa ~ ET + CT + lnht, data = predict_ht.melt, FUN = mean)

pht <- ggplot(ht.melt, aes(x = lnht, y = pa, color = CT, linetype = ET)) +
  labs(#title = expression(italic("Coreopsis major")), 
    #colour = "Legend", 
    y = NULL,#"Log total biomass (mg)", 
    x = "Plant height") +
  geom_line(data = predict_ht.melt, size = 2) +
  scale_colour_manual(values=c("khaki3","khaki4")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  scale_y_continuous(breaks = c(.01, .03)) +
  scale_x_continuous(breaks = c(-2,0,2,4)) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 24), axis.text.x = element_text(size = 24, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 24,  Color = "black"), 
        legend.position="none")


grid.arrange(psm, psla, pht, ncol = 1, left = textGrob("Probability of occurrence", gp=gpar(fontsize=34), rot = 90)) 

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

#### Plots ####




###############################################################################
################################ Beta diversity ###############################
###############################################################################


library(betapart)
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



############################## SAC ####################################
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

