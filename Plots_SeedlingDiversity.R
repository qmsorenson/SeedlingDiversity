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
library(tidyverse)

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


#########################################################################################################%%%%%
######################################           Figure 1           ##########################################
#########################################################################################################%%%%%

### Belowground competition bar plots --- all, seeded, volunteer ###

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/lsmr.RData")

# all species
pcas<- ggplot(data = lsmra$lsmeans.table[is.na(lsmra$lsmeans.table[,3]) & is.na(lsmra$lsmeans.table[,2]),], 
              aes(x = factor(ET, labels = c("E" = "Excluded", "P" = "Present")), y = exp(Estimate)-1, fill = factor(ET))) + 
  
  scale_fill_manual(values=c(rgb(105, 91, 78, maxColorValue = 255), "white")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1), position=position_dodge(.9),  width = 0.1) + 
  #labs(title = expression(paste("a) All species"))) +
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
        axis.text.x = element_text(size = 30, color = "black"),#, hjust = 1, vjust =.97),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")


### Land-use history x Canopy thinning interaction plots --- all, seeded, volunteer ###

lstab <- lsmra$lsmeans.table[complete.cases(lsmra$lsmeans.table[,2:3]) & is.na(lsmra$lsmeans.table[,1]),]
lstab$LUH <- factor(lstab$LUH, levels = c("R","M"))

pias <- ggplot(data = lstab, aes(x = factor(LUH), y = exp(Estimate)-1, colour = CT, shape = CT)) + 
  geom_line(data = lsmra$lsmeans.table[complete.cases(lsmra$lsmeans.table[,2:3]) & is.na(lsmra$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, group = CT, linetype = CT), show.legend = FALSE, size = 1.7 ) + 
  geom_point(data = lstab, aes(y = exp(Estimate)-1, fill = interaction(LUH:CT)), size = 10, stroke = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1, width = 0.1)) + 
  labs(title = NULL,#expression(paste("a) All species")),
       x = NULL, y = NULL, shape = "Tree density", 
       linetype = "Tree density") + 
  scale_x_discrete(labels = c("M" = "Post-ag.", "R" = "Remnant")) +
  scale_y_continuous(labels = fmt_dcimals(1)) +
  scale_colour_manual(values=c("#988F76", "#695B4E")) +
  scale_fill_manual(values=c("#9DBE78", "#66883F", "#EEDB74", "#E6C936")) +
  scale_shape_manual(values=c(21, 21)) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(), legend.title.align = .5,
        title =element_text(size=34),
        text = element_text(size = 30), axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")


grid.arrange(pcas, pias, ncol = 2, left = textGrob("        Species density", gp=gpar(fontsize=30), rot = 90), bottom = textGrob("Belowground Competition     Land-use history  \n ", gp=gpar(fontsize =34))) # pdf dim 10.75 5.9



#########################################################################################################%%%%%
######################################           Figure 2           ##########################################
#########################################################################################################%%%%%

### Environmental variable relationships

###### all seedlings ###

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/Enviro.data.RData")
load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/AllseededEnviroOutput.RData")


#scale_colour_manual(values=c("#988F76", "#695B4E")) +
# scale_fill_manual(values=c("#9DBE78", "#66883F", "#EEDB74", "#E6C936")) +

#wch
pw1 <- ggplot(data.cast.whc, 
              aes(x= whc, y = logrich.a, linetype = factor(ET, labels = c("BC Permitted", "BC Excluded")), shape = ET)) + 
  geom_point(alpha = .15, colour = "#35291F", show.legend=F) + labs(colour = "Legend", y = NULL,x = "Soil water holding capacity (%)") +
  geom_smooth( method = "lm", se = FALSE, colour = "#35291F", show.legend = F) + 
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  scale_shape_manual(values=c(19, 1), name = "Legend") +
  guides(linetype = guide_legend(direction = "horizontal", title = NULL), color = NULL ) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_rect(colour = 'Grey20', fill = 'white', linetype='solid'), 
        legend.key = element_blank(), legend.title.align = .5, aspect.ratio = .9,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 10,  Color = "black"), legend.position=c(.5,.9))
ggsave(pw1, file="Figure 4/Fig4_whc.pdf", width=2.94, height=2.67)

#light
pl2 <- ggplot(ear.light, 
              aes(x= X5cm, y = logrich.a, colour =  factor(CT, labels = c("T" = "Thinned",
                                                                          "U" = "Unthinned")))) + 
  geom_point(alpha = .15, show.legend = F) + labs(colour = "Legend", y = NULL,x = expression(paste("Light intensity (", mu, "mol/s)"))) +
  geom_smooth( method = "lm", se = FALSE, show.legend = F) + 
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_colour_manual(values=c("#CBC469", "#8B864E"))+ #, guide = guide_legend) +
  guides(colour = guide_legend(direction = "horizontal", title = NULL)) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_rect(colour = 'Grey20', fill = 'white', linetype='solid'), 
        legend.key = element_blank(), legend.title.align = .5, aspect.ratio = .9,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 10, Color = "black"), 
        legend.position=c(.5,.9), legend.box = "horizontal") # 3.32 x 2.67 in
ggsave(pl2, file="Figure 4/Fig4_lig.pdf", width=2.94, height=2.67)


#dbh
pd3 <- ggplot(data.cast.dbh[data.cast.dbh$CT == "U",], 
              aes(x= basal, y = logrich.a, colour =  factor(LUH, labels = c("Post.ag", "Remnant")),linetype = factor(ET), shape = ET)) + 
  geom_point(alpha = .15, show.legend = F) + labs(colour = "Legend", y = NULL,x = expression(paste("Basal area (cm"^"2", ")"))) +
  geom_smooth( method = "lm", se = F, show.legend = F) + 
  scale_colour_manual(values=c("#D7C769", "#74924F", "#D7C769", "74924F")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend", guide = FALSE) +
  scale_shape_manual(values=c(19, 1), name = "Legend") +
  guides(colour = guide_legend(direction = "horizontal", title = NULL)) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_rect(colour = 'Grey20', fill = 'white', linetype='solid'),
        legend.key = element_blank(), legend.title.align = .5, aspect.ratio = .9,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y",Size = 10, Color = "black"),
        legend.position=c(.5,.9), legend.box = "horizontal") # 3.32 x 2.67 in
ggsave(pd3, file="Figure 4/Fig4_dbh.pdf", width=2.94, height=2.67)

#gre
pg4 <- ggplot(ear.gre, 
              aes(x= green, y = logrich.a)) + 
  geom_point(alpha = .15, colour="#35291F") + labs(colour = "Legend", y = NULL,x = "Vegetation cover (%)") +
  geom_smooth(method = "lm", se = FALSE, colour="#35291F", show.legend = F) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5, aspect.ratio = .9,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y",Size = 10, Color = "black"), legend.position="none")
ggsave(pg4, file="Figure 4/Fig4_veg.pdf", width=2.94, height=2.67)


grid.arrange(pw1, pl2, pd3, pg4, ncol = 2, left = textGrob("Species density \nln(# species / exclosure)", gp=gpar(fontsize=14), rot = 90))




#########################################################################################################%%%%%
######################################           Figure 3           ##########################################
#########################################################################################################%%%%%


### MLM figures for seed mass, SLA, Max Height

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/trait.models.RData")
load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/trait.dataframes.RData")


# Seed mass 

predict_sm.melt <- sm.melt %>%
  expand(nesting(LUH,CT),lnsm)

predict_sm.melt$pa <- predict(m1sm, newdata = predict_sm.melt, re.form = NA, type="response")
#predict_sm.melt_no_et <- aggregate(pa ~ LUH + CT + lnsm, data = predict_sm.melt, FUN = mean)

psm <- ggplot(sm.melt, aes(x = lnsm, y = pa, color = factor(interaction(LUH,CT), labels = c("Post-ag. thinned", "Remnant thinned", "Post-ag. unthinned", "Remnant unthinned")))) +
  labs(#title = expression(italic("Coreopsis major")), 
    #colour = "Legend", 
    y = NULL,#"Log total biomass (mg)", 
    x = "Seed mass") +
  geom_line(data = predict_sm.melt, size = 2) +
  scale_colour_manual(values=c("gold2","darkolivegreen3", "gold4", "darkolivegreen")) +
  #scale_linetype_manual(values = c("solid", "twodash"), name = "Legend") +
  scale_y_log10(breaks = c(0.001, .020)) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.title=element_blank(),
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 30), axis.text.x = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 30,  Color = "black"), 
        legend.position=c(.74,.12), legend.box = "horizontal", legend.key.width = unit(3,"line"))

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
  scale_y_log10(breaks = c(0.001, .020)) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 30), axis.text.x = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 30,  Color = "black"), 
        legend.position="none")


# Max ht 

predict_ht.melt <- ht.melt %>%
  expand(nesting(CT,ET),lnht)

predict_ht.melt_pred_se <- as.data.frame(predict(m3ht, newdata = predict_ht.melt, re.form = NA, type="response", se.fit = TRUE))#,interval = "confidence")
#predict_ht.melt_no_luh <- aggregate(pa ~ ET + CT + lnht, data = predict_ht.melt, FUN = mean)
predict_ht.melt <- cbind(predict_ht.melt, predict_ht.melt_pred_se)
predict_ht.melt <- predict_ht.melt %>% rename(pa = fit)
predict_ht.melt$seu <- predict_ht.melt$pa + predict_ht.melt$se.fit
predict_ht.melt$sel <- predict_ht.melt$pa - predict_ht.melt$se.fit

pht <- ggplot(ht.melt, aes(x = lnht, y = pa, color = factor(CT, labels = c("Thinned", "Unthinned")), linetype = factor(ET, labels = c("BC Excluded", "BC Permitted")))) +
  labs(#title = expression(italic("Coreopsis major")), 
    #colour = "Legend", 
    y = NULL,#"Log total biomass (mg)", 
    x = "Plant height") +
  geom_ribbon(data = predict_ht.melt, aes(ymin = pa - se.fit, ymax = pa + se.fit), fill = "grey") +
  geom_line(data = predict_ht.melt, size = 2) +
  scale_colour_manual(values=c("khaki3","khaki4")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Legend") +
  scale_y_log10(breaks = c(0.001, .020)) +
  scale_x_continuous(breaks = c(-2,0,2,4)) +
  guides( color = guide_legend(title=NULL, size =1, title.position= "bottom"), 
          linetype = guide_legend(title = "Belowground \ncompetition")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_rect(fill = NA), legend.title.align = .5,
        text = element_text(size = 30), axis.text.x = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 30,  Color = "black"), 
        legend.position= "none" )#c(.82,.20), legend.box = "vertical",
       # legend.title = element_text(size = 24), legend.key.width = unit(3,"line"))  



grid.arrange(psm, psla, pht, ncol = 1, left = textGrob("Probability of occurrence (log)", gp=gpar(fontsize=36), rot = 90)) #2.75in x 6.95 -> 8.25, 20.85


#########################################################################################################%%%%%
#####################################           Appendix 1           #########################################
#########################################################################################################%%%%%


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


#########################################################################################################%%%%%
#####################################           Appendix 2           #########################################
#########################################################################################################%%%%%


### VOLUNTEER V SEEDED ####

pda <- data.cast[data.cast$brn == "NO",]

ggplot(pda, aes(x = logrich.s, y = logrich.v)) +
  #geom_point()#position = "jitter") +
  geom_count(mapping = NULL, data = NULL, stat = "sum",
   position = "identity", na.rm = FALSE, show.legend = NA,
  inherit.aes = TRUE) +
  geom_smooth( method = "lm", se = FALSE, color = "black") + 
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.title.align = .5,
        text = element_text(size = 12), axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))


#########################################################################################################%%%%%
#####################################           Appendix 3           #########################################
#########################################################################################################%%%%%

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



#########################################################################################################%%%%%
######################################           Figure S1           ########################################
#########################################################################################################%%%%%

### Seeded species ###

pcss<- ggplot(data = lsmrs$lsmeans.table[is.na(lsmrs$lsmeans.table[,3]) & is.na(lsmrs$lsmeans.table[,2]),], 
              aes(x = factor(ET, labels = c("E" = "Excluded", "P" = "Present")), y = exp(Estimate)-1, fill = factor(ET))) + 
  
  scale_fill_manual(values=c(rgb(105, 91, 78, maxColorValue = 255), "white")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1), position=position_dodge(.9),  width = 0.1) + 
  #labs(title = expression(paste("a) All species"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2), labels = fmt_dcimals(1)) +
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
        text = element_text(size = 24),
        #legend.spacing = unit(2,"cm"),
        #legend.key = element_rect(size = 10),
        #legend.key.size = unit(2, 'lines'),
        axis.text.x = element_text(size = 24, color = "black"),#, hjust = 1, vjust =.97),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 24, Color = "black"), legend.position="none")


### Land-use history x Canopy thinning interaction plots --- all, seeded, volunteer ###

lstab <- lsmrs$lsmeans.table[complete.cases(lsmrs$lsmeans.table[,2:3]) & is.na(lsmrs$lsmeans.table[,1]),]
lstab$LUH <- factor(lstab$LUH, levels = c("R","M"))
piss <- ggplot(data = lstab, aes(x = factor(LUH), y = exp(Estimate)-1, colour = CT, shape = CT)) + 
  geom_line(data = lsmrs$lsmeans.table[complete.cases(lsmrs$lsmeans.table[,2:3]) & is.na(lsmrs$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, group = CT, linetype = CT), show.legend = FALSE, size = 1.7 ) + 
  geom_point(data = lstab, aes(y = exp(Estimate)-1, fill = interaction(LUH:CT)), size = 10, stroke = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1, width = 0.1)) + 
  labs(title = NULL,#expression(paste("a) All species")),
       x = NULL, y = NULL, shape = "Tree density", 
       linetype = "Tree density") + 
  scale_x_discrete(labels = c("M" = "Post-ag.", "R" = "Remnant")) +
  scale_y_continuous(labels = fmt_dcimals(1)) +
  scale_colour_manual(values=c("#988F76", "#695B4E")) +
  scale_fill_manual(values=c("#9DBE78", "#66883F", "#EEDB74", "#E6C936")) +
  scale_shape_manual(values=c(21, 21)) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(), legend.title.align = .5,
        title =element_text(size=24),
        text = element_text(size = 24), axis.text.x = element_text(size = 24, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 24, Color = "black"), legend.position="none")


### Non-seeded species

pcvs<- ggplot(data = lsmrv$lsmeans.table[is.na(lsmrv$lsmeans.table[,3]) & is.na(lsmrv$lsmeans.table[,2]),], 
              aes(x = factor(ET, labels = c("E" = "Excluded", "P" = "Present")), y = exp(Estimate)-1, fill = factor(ET))) + 
  
  scale_fill_manual(values=c(rgb(105, 91, 78, maxColorValue = 255), "white")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1), position=position_dodge(.9),  width = 0.1) + 
  #labs(title = expression(paste("a) All species"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2), labels = fmt_dcimals(1)) +
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
        text = element_text(size = 24),
        #legend.spacing = unit(2,"cm"),
        #legend.key = element_rect(size = 10),
        #legend.key.size = unit(2, 'lines'),
        axis.text.x = element_text(size = 24, color = "black"),#, hjust = 1, vjust =.97),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 24, Color = "black"), legend.position="none")


### Land-use history x Canopy thinning interaction plots --- all, seeded, volunteer ###

lstab <- lsmrv$lsmeans.table[complete.cases(lsmrv$lsmeans.table[,2:3]) & is.na(lsmrv$lsmeans.table[,1]),]
lstab$LUH <- factor(lstab$LUH, levels = c("R","M"))
pivs <- ggplot(data = lstab, aes(x = factor(LUH), y = exp(Estimate)-1, colour = CT, shape = CT)) + 
  geom_line(data = lsmrv$lsmeans.table[complete.cases(lsmrv$lsmeans.table[,2:3]) & is.na(lsmrv$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, group = CT, linetype = CT), show.legend = FALSE, size = 1.7 ) + 
  geom_point(data = lstab, aes(y = exp(Estimate)-1, fill = interaction(LUH:CT)), size = 10, stroke = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1, width = 0.1)) + 
  labs(title = NULL,#expression(paste("a) All species")),
       x = NULL, y = NULL, shape = "Tree density", 
       linetype = "Tree density") + 
  scale_x_discrete(labels = c("M" = "Post-ag.", "R" = "Remnant")) +
  scale_y_continuous(labels = fmt_dcimals(1)) +
  scale_colour_manual(values=c("#988F76", "#695B4E")) +
  scale_fill_manual(values=c("#9DBE78", "#66883F", "#EEDB74", "#E6C936")) +
  scale_shape_manual(values=c(21, 21)) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(), legend.title.align = .5,
        title =element_text(size=24),
        text = element_text(size = 24), axis.text.x = element_text(size = 24, color = "black"),
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 24, Color = "black"), legend.position="none")



figs1 <- grid.arrange(pcss, piss, pcvs, pivs, ncol = 2, left = textGrob("Species density", gp=gpar(fontsize=26), rot = 90), bottom = textGrob("Belowground Competition     Land-use history  \n ", gp=gpar(fontsize =26))) # pdf dim 10.75 11.8
ggsave(figs1, file = "Figure S1/FigS1.pdf", width = 10.75, height = 11.8)





#########################################################################################################%%%%%
######################################           Figure S2           ########################################
#########################################################################################################%%%%%


load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/Enviro.data.RData")
load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/AllseededEnviroOutput.RData")

vio_whc <- ggplot(data=data.cast.all, aes(x = interaction(LUH, CT), y = whc, color = interaction(LUH, CT))) + #, fill = interaction(LUH, CT)
  geom_violin(size = 1) +
  #geom_boxplot(width=0.1) +
  geom_jitter(width = .02, alpha = .25) +
  #geom_dotplot(binaxis='y', stackdir='center',binwidth = .01, dotsize=.7, aes(fill = interaction(LUH, CT))) +
  scale_color_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_fill_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_x_discrete(labels = c("M.T" = "Post-ag \nThinned", "M.U" = "Post-ag \nUnthinned", "R.T" = "Remnant \nThinned", "R.U" = "Remnant \nUnthinned")) +
  ylab("Soil water holding capacity (%)") +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none", 
        text = element_text(size = 14), axis.text.x = element_text(size = 12, color = "black"), aspect.ratio = .9,
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 16, Color = "black")) +
  labs(title ="", fill ="")
ggsave(vio_whc, file="Figure S2/FigS2_whc.pdf", width=4.5, height=4.15)


vio_lig <- ggplot(data=data.cast.all, aes(x = interaction(LUH, CT), y = X5cm, color = interaction(LUH, CT))) + #, fill = interaction(LUH, CT)
  geom_violin(size = 1) +
  #geom_boxplot(width=0.1) +
  geom_jitter(width = .02, alpha = .25) +
  #geom_dotplot(binaxis='y', stackdir='center',binwidth = 15, dotsize=1.5, aes(fill = interaction(LUH, CT))) +
  scale_color_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_fill_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_x_discrete(labels = c("M.T" = "Post-ag \nThinned", "M.U" = "Post-ag \nUnthinned", "R.T" = "Remnant \nThinned", "R.U" = "Remnant \nUnthinned")) +
  ylab(expression(paste("Light intensity (", mu, "mol/s)"))) +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none", 
        text = element_text(size = 14), axis.text.x = element_text(size = 12, color = "black"), aspect.ratio = .9,
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 16, Color = "black")) +
  labs(title ="", fill ="")
ggsave(vio_lig, file="Figure S2/FigS2_lig.pdf", width=4.5, height=4.15)


vio_dbh_t <- ggplot(data=data.cast.all %>% filter(CT == "T", basal< 750), aes(x = LUH, y = basal, color = LUH)) + #, fill = interaction(LUH, CT)
  geom_violin(size = 1) +
  #geom_boxplot(width=0.1) +
  geom_jitter(width = .02, alpha = .25) +
  #geom_dotplot(binaxis='y', stackdir='center',binwidth = .01, dotsize=.7, aes(fill = interaction(LUH, CT))) +
  scale_color_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_fill_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_x_discrete(labels = c("M" = "Post-ag \nThinned", "M.U" = "Post-ag \nUnthinned", "R" = "Remnant \nThinned", "R.U" = "Remnant \nUnthinned")) +
  ylab(expression(paste("Basal area (cm"^"2", ")"))) +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none", 
        text = element_text(size = 14), axis.text.x = element_text(size = 12, color = "black"), aspect.ratio = 1.8,
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 16, Color = "black")) +
  labs(title ="", fill ="")
ggsave(vio_dbh_t, file="Figure S2/FigS2_dbh_t.pdf", width=4.5, height=4.15)


vio_dbh_u <- ggplot(data=data.cast.all %>% filter(CT == "U"), aes(x = LUH, y = basal, color = LUH)) + #, fill = interaction(LUH, CT)
  geom_violin(size = 1) +
  #geom_boxplot(width=0.1) +
  geom_jitter(width = .02, alpha = .25) +
  #geom_dotplot(binaxis='y', stackdir='center',binwidth = .01, dotsize=.7, aes(fill = interaction(LUH, CT))) +
  scale_color_manual(values=c("#AD8516", "#66883F")) +
  scale_fill_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_x_discrete(labels = c("M." = "Post-ag \nThinned", "M" = "Post-ag \nUnthinned", "R." = "Remnant \nThinned", "R" = "Remnant \nUnthinned")) +
  ylab(expression(paste("Basal area (cm"^"2", ")"))) +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none", 
        text = element_text(size = 14), axis.text.x = element_text(size = 12, color = "black"), aspect.ratio = 1.8,
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 16, Color = "black")) +
  labs(title ="", fill ="")
ggsave(vio_dbh_u, file="Figure S2/FigS2_dbh_u.pdf", width=4.5, height=4.15)


vio_veg <- ggplot(data=data.cast.all, aes(x = interaction(LUH, CT), y = green, color = interaction(LUH, CT))) + #, fill = interaction(LUH, CT)
  geom_violin(size = 1) +
  #geom_boxplot(width=0.1) +
  geom_jitter(width = .02, alpha = .25) +
  #geom_dotplot(binaxis='y', stackdir='center',binwidth = .01, dotsize=.7, aes(fill = interaction(LUH, CT))) +
  scale_color_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_fill_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_x_discrete(labels = c("M.T" = "Post-ag \nThinned", "M.U" = "Post-ag \nUnthinned", "R.T" = "Remnant \nThinned", "R.U" = "Remnant \nUnthinned")) +
  ylab("Vegentation cover (%)") +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none", 
        text = element_text(size = 14), axis.text.x = element_text(size = 12, color = "black"), aspect.ratio = .9,
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 16, Color = "black")) +
  labs(title ="", fill ="")
ggsave(vio_veg, file="Figure S2/FigS2_veg.pdf", width=4.5, height=4.15)

vio_comp <- ggplot(data=data.castcomp, aes(x = interaction(LUH, CT), y = comp, color = interaction(LUH, CT))) + #, fill = interaction(LUH, CT)
  geom_violin(size = 1) +
  #geom_boxplot(width=0.1) +
  #geom_jitter(width = .02, alpha = .25) +
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 6, dotsize=.7, aes(fill = interaction(LUH, CT))) +
  scale_color_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_fill_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_x_discrete(labels = c("M.T" = "Post-ag \nThinned", "M.U" = "Post-ag \nUnthinned", "R.T" = "Remnant \nThinned", "R.U" = "Remnant \nUnthinned")) +
  ylab("compaction") +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none", 
        text = element_text(size = 14), axis.text.x = element_text(size = 12, color = "black"), aspect.ratio = .9,
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 16, Color = "black")) +
  labs(title ="", fill ="")
ggsave(vio_comp, file="Figure S2/FigS2_comp.pdf", width=4.5, height=4.15)






gif <- ggplot(data=data.cast.all %>% gather(Lig_level, ppfd, X5cm, X15cm, X25cm) , aes(x = interaction(LUH, CT), y = ppfd, color = interaction(LUH, CT))) + #, fill = interaction(LUH, CT)
  geom_violin(size = 2) +
  #geom_boxplot(width=0.1) +
  #geom_jitter(width = .02, alpha = .25) +
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 15, dotsize=1.5, aes(fill = interaction(LUH, CT))) +
  labs(title = "{closest_state}", x = NULL, y = expression(paste("Light intensity (", mu, "mol/s)"))) +
  transition_states(Lig_level, state_length = 15, transition_length = 8) +
  scale_color_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_fill_manual(values=c("#EDDB73", "#9DBE78","#AD8516", "#66883F")) +
  scale_x_discrete(labels = c("M.T" = "Post-ag \nThinned", "M.U" = "Post-ag \nUnthinned", "R.T" = "Remnant \nThinned", "R.U" = "Remnant \nUnthinned")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(), legend.position="none", 
        text = element_text(size = 40), axis.text.x = element_text(size = 36, color = "black"), aspect.ratio = .9,
        axis.text.y = rotatedAxisElementText(angle = 90, position = "y", Size = 36, Color = "black")) +
  ease_aes('sine-in-out')
  
gifbig <- animate(gif, width = 880, height = 740)
anim_save(file = "Light_gif.gif", animation = gifbig) 

##################################################################################################################
##############################################   TRASH   #########################################################
##################################################################################################################



##############################################################################################################
######################################           Figure 1           ##########################################
##############################################################################################################

### Belowground competition bar plots --- all, seeded, volunteer ###

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/lsmr.RData")

pcas<- ggplot(data = lsmra$lsmeans.table[is.na(lsmra$lsmeans.table[,3]) & is.na(lsmra$lsmeans.table[,2]),], 
              aes(x = factor(ET, labels = c("E" = "Excluded", "P" = "Present")), y = exp(Estimate)-1, fill = factor(ET))) + 
  
  scale_fill_manual(values=c(rgb(105, 91, 78, maxColorValue = 255), "white")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1), position=position_dodge(.9),  width = 0.1) + 
  #labs(title = expression(paste("a) All species"))) +
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
        axis.text.x = element_text(size = 30, color = "black"),#, hjust = 1, vjust =.97),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")


### Land-use history x Canopy thinning interaction plots --- all, seeded, volunteer ###

lstab <- lsmra$lsmeans.table[complete.cases(lsmra$lsmeans.table[,2:3]) & is.na(lsmra$lsmeans.table[,1]),]
lstab$LUH <- factor(lstab$LUH, levels = c("R","M"))
pias <- ggplot(data = lstab, aes(x = factor(LUH), y = exp(Estimate)-1, colour = CT, shape = CT)) + 
  geom_line(data = lsmra$lsmeans.table[complete.cases(lsmra$lsmeans.table[,2:3]) & is.na(lsmra$lsmeans.table[,1]),], aes(y = exp(Estimate)-1, group = CT, linetype = CT), show.legend = FALSE, size = 1.7 ) + 
  geom_point(data = lstab, aes(y = exp(Estimate)-1, fill = interaction(LUH:CT)), size = 10, stroke = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymax = exp(Estimate +`Standard Error`)-1, ymin=exp(Estimate -`Standard Error`)-1, width = 0.1)) + 
  labs(title = NULL,#expression(paste("a) All species")),
       x = NULL, y = NULL, shape = "Tree density", 
       linetype = "Tree density") + 
  scale_x_discrete(labels = c("M" = "Post-ag.", "R" = "Remnant")) +
  scale_y_continuous(labels = fmt_dcimals(1)) +
  scale_colour_manual(values=c("#988F76", "#695B4E")) +
  scale_fill_manual(values=c("#9DBE78", "#66883F", "#EEDB74", "#E6C936")) +
  scale_shape_manual(values=c(21, 21)) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Legend") +
  theme(panel.background = element_rect(fill = NA, color = "black"), panel.grid = element_blank(),
        legend.background = element_blank(), aspect.ratio = .9,
        legend.key = element_blank(), legend.title.align = .5,
        title =element_text(size=34),
        text = element_text(size = 30), axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"), legend.position="none")


grid.arrange(pcas, pias, ncol = 2, left = textGrob("        Species density", gp=gpar(fontsize=30), rot = 90), bottom = textGrob("Belowground Competition     Land-use history  \n ", gp=gpar(fontsize =34))) # pdf dim 10.75 5.9

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


