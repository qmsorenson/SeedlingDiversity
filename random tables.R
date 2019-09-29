library(sjPlot)
library(sjmisc)
data(iris)
head(iris)

tab_df(iris,
       file="sjt_des.doc")

tab_model(m1)



###############################################################

#plot_grid() 
#Arrange list of plots as grid
#Description
#Plot multiple ggplot-objects as a grid-arranged single plot.

#sjp.pca Plot PCA results
#Description
#Performs a principle component analysis on a data frame or matrix (with varimax or oblimin rotation)
#and plots the factor solution as ellipses or tiles.
#In case a data frame is used as argument, the cronbach’s alpha value for each factor scale will
#be calculated, i.e. all variables with the highest loading for a factor are taken for the reliability test.
#The result is an alpha value for each factor dimension.

