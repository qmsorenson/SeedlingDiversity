library(reshape2)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(Hmisc)
library(MuMIn)
library(lsmeans)
library(effects)
library(ggplot2)
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


##############################################################################################################
######################################           Figure 1           ##########################################
##############################################################################################################

### Belowground competition bar plots --- all, seeded, volunteer ###

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/lsmr.RData")

# all species
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

# seeded species
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

# volunteer species
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




### Land-use history x Canopy thinning interaction plots --- all, seeded, volunteer ###

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



##############################################################################################################
######################################           Figure 2           ##########################################
##############################################################################################################

### Environmental variable relationships

###### all seedlings ######

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/Enviro.data.RData")
load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/AllseededEnviroOutput.RData")

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
  scale_colour_manual(values=c("khaki3","khaki4")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"), legend.position="none")

#dbh
pd3 <- ggplot(data.cast.dbh[data.cast.dbh$CT == "U",], 
              aes(x= basal, y = logrich.a, colour =  factor(LUH),linetype = factor(ET), shape = ET)) + 
  geom_point() + labs(colour = "Legend", y = NULL,x = expression(paste("Basal area (cm"^"2", ")"))) +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_colour_manual(values=c("gold3","darkolivegreen3")) +
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




##############################################################################################################
######################################           Figure 3           ##########################################
##############################################################################################################


### MLM figures for seed mass, SLA, Max Height

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/trait.models.RData")

# Seed mass 

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



# Specific leaf area

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


# Max ht 

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


##############################################################################################################
#####################################           Appendix 1           #########################################
##############################################################################################################


### Bar graphs of richness

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


##############################################################################################################
#####################################           Appendix 2           #########################################
##############################################################################################################


### VOLUNTEER V SEEDED ####

pda <- data.cast[data.cast$brn == "NO",]

ggplot(pda, aes(x = logrich.s, y = logrich.v)) +
  geom_point(position = "jitter") +
  #geom_count(mapping = NULL, data = NULL, stat = "sum",
  # position = "identity", na.rm = FALSE, show.legend = NA,
  #inherit.aes = TRUE) +
  geom_smooth( method = "lm", se = FALSE, color = "black") + 
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))


##############################################################################################################
#####################################           Appendix 3           #########################################
##############################################################################################################

### Environmental variables seeded and unseeded species


########## seeded species ##########

#wch
ggplot(data.cast.whc, 
       aes(x= whc, y = logrich.s,linetype = factor(ET))) + 
  geom_point() + labs(title = "seeded richness", colour = "Legend", y = "ln # of spp", x = "Soil water holding capacity (%)") +
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
  geom_point() + labs(title = "Seeded richness", colour = "Legend", y = "ln # of spp", x = "PPFD") +
  geom_smooth( method = "lm", se = FALSE) + 
  scale_colour_manual(values=c("Grey60","blue")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#dbh
ggplot(data.cast.dbh[data.cast.dbh$CT == "U",], 
       aes(x= BA, y = logrich.s, colour =  factor(LUH),linetype = factor(ET))) + 
  geom_point() + labs(title = "Seeded richness", colour = "Legend", y = "ln # of spp", x = "DBH") +
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
  geom_point() + labs(title = "Seeded richness", colour = "Legend", y = "ln # of spp", x = "% vegetation surrounding") +
  geom_smooth( method = "lm", se = FALSE) + 
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))


########### volunteer species ############

#wch
ggplot(data.cast.whc, 
       aes(x= whc, y = logrich.v,linetype = factor(ET))) + 
  geom_point() + labs(title = "Volunteer richness", colour = "Legend", y = "ln # of spp", x = "Soil water holding capacity (%)") +
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
  geom_point() + labs(title = "Volunteer richness", colour = "Legend", y = "ln # of spp", x = "PPFD") +
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
  geom_point() + labs(title = "Volunteer richness", colour = "Legend", y = "ln # of spp", x = "DBH") +
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
  geom_point() + labs(title = "Volunteer richness", colour = "Legend", y = "ln # of spp", x = "% vegetation surrounding") +
  geom_smooth( method = "lm", se = FALSE) + 
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