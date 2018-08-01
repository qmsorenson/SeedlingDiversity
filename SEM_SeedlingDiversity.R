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