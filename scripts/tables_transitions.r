##
#Prueba de bondad de ajuste
1-pchisq(0.47,3) ###para el ejercicio a mano
chicharo <- chisq.test(c(315,108,101,32), p=c(9/16, 3/16, 3/16, 1/16))
chicharo

tabla<-read.csv("../Data/Pedilanthus/tabla_transiciones.csv", sep=",")



#PLot longitudes
read.delim("../Data/Pedilanthus/P_bracteatus/P_bracteatus_long.txt")
