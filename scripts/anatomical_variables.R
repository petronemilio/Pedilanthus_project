#####Load anatomical variables #####
anatomical_variables<- read.csv("Data/Raw_data/Anatomical_variables.csv")
anatomical_variables_cacho <- read.csv("Data/Raw_data/anatomical_variables_cacho.csv") 

boxplot(log10(anatomical_variables_cacho$vessel_wall_thickness) ~
          log10(anatomical_variables_cacho$vessel_diameter))

plot(log10(anatomical_variables$Grosor_pared_vaso) ~ log10(anatomical_variables$Diametro_de_vaso))
mean()

vd.mean.by.sample <-aggregate(anatomical_variables$Diametro_de_vaso,
                               by=list(anatomical_variables$Muestra), mean, na.rm=TRUE)

summary(vd.mean.by.sample)
