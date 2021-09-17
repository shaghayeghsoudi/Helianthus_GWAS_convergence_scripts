#https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
rm(list=ls())
library(factoextra)
setwd("~/Desktop")
env_data<-read.csv(file = "env_data_all_species.csv",sep = ",", header = TRUE)
env_data<-env_data[!(env_data$species=="H.petiolaris.canescens"),]

#env_data_a<-na.omit(env_data_a)
env_data$Elevation<-as.integer(gsub(",","",env_data$Elevation))
env_data_good<-na.omit(env_data)
rownames(env_data_good)<-env_data_good[,1]

env_data_a<-env_data_good[,c(-1,-2)]

data_climate_good_scaled<-scale(env_data_a,center = TRUE, scale = TRUE)
#summary(data_climate_good_scaled)
#colnames(data_climate_good_scaled)


pca_res <- prcomp(data_climate_good_scaled)
#autoplot(pca_res, data = env_data_good, colour = 'species',alpha = 0.5)
#autoplot(pca_res, data = env_data_good, colour = 'species', loadings = TRUE)


autoplot(pca_res, data = env_data_good, colour = 'species', loadings = TRUE,
        loadings.colour = "blue",size=3,
        loadings.label = TRUE, loadings.label.size = 3,loadings.label.colour = "blue",alpha = 0.7)



####################################
##### selected variable (AHM, MAT)
env_data
ahm<-env_data[,c("species","Latitude")]
p <- ggplot(ahm, aes(x=species, y=Latitude,fill=species)) + 
  geom_violin(trim=FALSE)

#Add median and quartile
p+geom_boxplot(width=0.1)
p + geom_jitter(shape=16, position=position_jitter(0.2))



data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

p + stat_summary(fun.data=data_summary)









