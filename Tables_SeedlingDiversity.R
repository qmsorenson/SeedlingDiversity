library(dplyr)
library(sjPlot)


##### Table S1 ---------------------------------------------------------------

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/Enviro.data.RData")

cor(data.cast.all %>% select(whc, X5cm, basal, green), use = "complete.obs", method = "pearson")


##### Table S2 ---------------------------------------------------------------

load("G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/R_dataframes/trait.dataframes.RData")

sm.spplist <- sm.melt %>%
  group_by(Spp) %>%
  summarise(Spp_n = sum(pa), SeedMass_mg = mean(SeedMass_mg)) %>%
  arrange(Spp_n)

sla.spplist <- sla.melt %>%
  group_by(Spp) %>%
  summarise(Spp_n = sum(pa),sla = mean(sla)) %>%
  arrange(Spp_n)

ht.spplist <- ht.melt %>%
  group_by(Spp) %>%
  summarise(Spp_n = sum(pa), MaxPhotosynthetic_cm = mean(MaxPhotosynthetic_cm)) %>%
  arrange(Spp_n)

spplist <- sm.spplist %>% left_join(sla.spplist %>% select(Spp, sla), by = "Spp")
spplist <- spplist %>% left_join(ht.spplist %>% select(Spp, MaxPhotosynthetic_cm), by = "Spp")
#write.csv(spplist, "G:/My Drive/Graduate School/Research/Remnant/Root exclosure/Seedling diversity/data/SpeciesList.csv")
cor(spplist %>% select(SeedMass_mg, sla, MaxPhotosynthetic_cm),  use = "complete.obs", method = "pearson")


##### Table S3 -------------------------------------------------------------

tab_model(rich.a, dv.labels = " ", digits.p = 4, file = "Tables/Table_S3.xls",
          pred.labels = c("Intercept",
                                                "Belowground competition (BGC)",
                                                "Land-use history (LUH)",
                                                "Canopy thinning (CT)",
                                                "LUH:CT"))


##### Table3 S4-7 -------------------------------------------------------------

tab_model(rich.a.whc, dv.labels = " ", digits.p = 4,
          pred.labels = c("Intercept",
                                                                                        "Belowground competition (BGC)",
                                                                                        "Land-use history (LUH)",
                                                                                        "Canopy thinning (CT)",
                                                                                        "Soil water holding capacity(WHC)",
                                                                                        "LUH:CT",
                                                                                        "BGC:WHC"
                                                                                       ))

tab_model(rich.a.lig, dv.labels = " ", digits.p = 4, file = "Tables/Table_S5.xls",
          pred.labels = c("Intercept",
                                                                                        "Belowground competition (BGC)",
                                                                                        "Land-use history (LUH)",
                                                                                        "Canopy thinning (CT)",
                                                                                        "Light intensity (LI)",
                                                                                        "LUH:CT",
                                                                                        "CT:LI"
                                                                                        ))

tab_model(rich.a.dbh, dv.labels = " ", digits.p = 4, file = "Tables/Table_S6.xls",
          pred.labels = c("Intercept",
                                                                                        "Belowground competition (BGC)",
                                                                                        "Land-use history (LUH)",
                                                                                        "Basal area (BA)",
                                                                                        "BGC:LUH",
                                                                                        "BGC:BA",
                                                                                        "LUH:BA",
                                                                                        "BGC:LUH:BA"))

tab_model(rich.a.gre, dv.labels = " ", digits.p = 4, file = "Tables/Table_S7.xls",
          pred.labels = c("Intercept",
                                                                                        "Belowground competition (BGC)",
                                                                                        "Land-use history (LUH)",
                                                                                        "Canopy thinning (CT)",
                                                                                        "vegetation cover",
                                                                                        "LUH:CT"))


##### Tables S8-10 -------------------------------------------------------------

tab_model(m1, digits.p = 4, dv.labels = " ", 
          file = "Tables/Table_S8.xls",
          pred.labels = c("Intercept",
                                            "Belowground competition (BGC)",
                                            "Land-use history (LUH)",
                                            "Canopy thinning (CT)",
                                            "Seed mass (SM)",
                                            "LUH:CT",
                                            "LUH:SM",
                                            "CT:SM",
                                            "LUH:CT:SM"
                                            ))

tab_model(m2, digits.p = 4, dv.labels = " ", file = "Tables/Table_S9.xls", pred.labels = c("Intercept",
                                            "Belowground competition (BGC)",
                                            "Land-use history (LUH)",
                                            "Canopy thinning (CT)",
                                            "Specific leaf area (SLA)",
                                            "LUH:CT",
                                            "BGC:CT",
                                            "CT:SLA",
                                            "BGC:SLA",
                                            "BGC:CT:SLA"))

tab_model(m3, digits.p = 4, dv.labels = " ", file = "Tables/Table_S10.xls",
          pred.labels = c("Intercept",
                                            "Belowground competition (BGC)",
                                            "Land-use history (LUH)",
                                            "Canopy thinning (CT)",
                                            "Max. Height (MH)",
                                            "LUH:CT",
                                            "BGC:CT",
                                            "CT:MH",
                                            "BGC:MH",
                                            "BGC:CT:MH"))


##### Table S11 -------------------------------------------------------------

tab_model(rich.s, rich.v, digits.p = 4, dv.labels = " ", file = "Tables/Table_S11.xls",
          pred.labels = c("Intercept",
                                                       "Belowground competition (BGC)",
                                                       "Land-use history (LUH)",
                                                       "Canopy thinning (CT)",
                                                       "LUH:CT"))


##### mods for Table S3-11 -------------------------------------------------------------

#rich.a <- lmer(logrich.a ~ ET+LUH*CT + (1|Site/LUH:CT) + (1|Site:LUH:CT:SP), data = data.cast[data.cast$brn == "NO",])

rich.a.whc <- lmer(logrich.a ~ ET + LUH*CT + ET*whc + (1|Site/LUH:CT) + (1|Site:LUH:CT:SP), data = data.cast.whc)
anova(rich.a.whc)
summary(rich.a.whc)

#rich.a.lig <- lmer(logrich.a ~ ET + LUH*CT + CT*X5cm + (1|Site/LUH:CT) + (1|Site:LUH:CT:SP), data = ear.light)

#rich.a.dbh <- lmer(logrich.a ~ ET*LUH*BA +(1|Site/LUH), data = data.cast.dbh[data.cast.dbh$CT == "U",])

#rich.a.gre <- lmer(logrich.a ~ ET+LUH*CT+green + (1|Site/LUH:CT) + (1|Site:LUH:CT:SP), data = ear.gre)

#m1 <- glmer(pa ~ ET + LUH*CT*scale(log(SeedMass_mg)) + (1|Site) + (1|Site:LUH:CT) + (1|Site:LUH:CT:SP2) + (1 + LUH*CT|Spp),
#            nAGQ = 0, family = "binomial", data=sm.melt)

#m2 <- glmer(pa ~ ET + LUH + CT:LUH + CT*ET*scale(sla) + (1|Site) + (1|Site:LUH:CT) + (1|Site:LUH:CT:SP2) + (1 + ET*CT|Spp),
#            family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=sla.melt)

#m3 <- glmer(pa ~ LUH + LUH:CT + CT*ET*scale(log(MaxPhotosynthetic_cm)) + (1|Site) + (1|Site:LUH:CT) + (1|Site:LUH:CT:SP2) +
#             (1+ CT*ET|Spp), family = "binomial", control=glmerControl(calc.derivs=F), nAGQ = 0, data=ht.melt)

#rich.s <- lmer(logrich.s ~ ET+LUH*CT + (1|Site/LUH:CT) + (1|Site:LUH:CT:SP), data = data.cast[data.cast$brn == "NO",])

#rich.v <- lmer(logrich.v ~ ET+CT + (1|Site/LUH:CT) + (1|Site:LUH:CT:SP), data = data.cast[data.cast$brn == "NO",])

