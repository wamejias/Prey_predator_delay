#Descripción:
#Este script puede ser visto como una plantilla 
#que puede adaptarse para realizar simulaciones numéricas
#de modelos que involucran sistemas 
#de Ecuaciones Diferenciales con Retardo (EDR).

#Para más información acerca de la resolución de EDR con R
#se recomienda ir a los capítulas 6 y 7 de:
#Soetaert, K., Cash, J., & Mazzia, F. (2012). Solving Differential Equations in R. Springer, Berlin, Heidelberg.

#Begin script

#Librerias
library(ggplot2)
library(reshape2)
library(gridExtra)
library(deSolve)
library(latex2exp)

#Se define el lado derecho del sistema de EDR
LVdede <- function(t, y, p) {
  if (t > tau1) Lag1 <- lagvalue(t - tau1) else Lag1 <- yini
  dy1=r*((K-y[1])/K)*y[1]-((g*y[1]^2)*y[2])/(a^2+y[1]^2)
  dy2=-b*y[2]+((g*Lag1[1]^2)/(a^2+Lag1[1]^2))*Lag1[2]*exp(-Lag1[2]/N)
  list(c(dy1, dy2))
}

#Listas de parametros
lista_r=c(0.05,0.1,0.2) #Tasa nacimiento presas mensual
lista_K=c(200, 500) #Capacidad de presas
lista_a=c(0.1,0.5,0.8)
lista_b=c(0.05,0.1) #Mortaliadad mensual puma
lista_g=c(0.1,0.5,0.8)
lista_N=c(1,2) #relacionado con capacidad de carga dos situaciones
tau1 <- 27 #Meses
yini <- c(y1 = 200, y2 = 1)

#Ciclo de graficos tau fijo

contador=1  
for (i in lista_r){
  for(j in lista_K){
    for(k in lista_a){
      for(l in lista_N){
        for (m in lista_g){
          for (n in lista_b){
            times <- seq(from = 0, to = 1800, by = 0.1)
            r <- i; K <- j; a <- k; b <- n; g <-m; N<-l 
            leyenda<-paste("r =",r,"K =",K,"a =", a,"\U1D6FD =",b,"\U1D6FE=",g, "N =",N)
            yout1<- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0))
            yout1 <- cbind(yout1 , c1 = "Density puma and prey")
            yout1 <- cbind(yout1 , c2 = "Numeric simulation")
            p1<-ggplot(data=yout1, aes(x=time, y=log(y1),group=1,colour = "Log(Prey)",))+
              geom_line()+
              geom_line(data=as.data.frame(yout1), aes(x=time, y=y2,group=1,colour = "Predator (Puma)"))+
              theme_bw()+
              labs(y = expression("Density N°ind/100km"^2) ,x="Time (Months)", color = "Species", linetype="Delay")+
              facet_wrap(~ c2)
            p2<-ggplot(data=yout1, aes(x=y1, y=y2 ))+
              geom_path(colour='red')+
              theme_bw()+
              labs(y = expression("Predator N°ind/100km"^2) ,x=expression("Prey N°ind/100km"^2), color = "Delay")+
              expand_limits(x = 0,y = 0)+
              geom_label(aes(x = max(y1)/1.9 , y = max(y2)+0.25, label =leyenda), stat = "unique",color="black",size = 2.8)+
              facet_wrap(~ c1)
            png(paste0("Grafico",contador,".png"), width = 4800, height = 2400, res=600)
            grid.arrange(p2,p1,nrow=1)
            dev.off() 
            print(contador)
            contador=contador+1
          }
        }
      }
    }
  }
}


#Ciclo de graficos tau variable

contador=1  
for (i in lista_r){
  for(j in lista_K){
    for(k in lista_a){
      for(l in lista_N){
        for (m in lista_g){
          for (n in lista_b){
            times <- seq(from = 0, to = 1800, by = 0.1)
            r <- i; K <- j; a <- k; b <- n; g <-m; N<-l 
            tau1 <- 27 #Meses
            yout1<- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0))
            yout1 <- cbind(yout1 , c1 = "Density puma and prey")
            yout1 <- cbind(yout1 , c2 = "Numeric simulation")
            tau1 <-13 #Meses
            yout2<- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0))
            yout2 <- cbind(yout2 , c1 = "Density puma and prey")
            yout2 <- cbind(yout2 , c2 = "Numeric simulation")
            tau1 <-5  #Meses
            yout3 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0))
            yout3 <- cbind(yout3 , c1 = "Density puma and prey")
            yout3 <- cbind(yout3 , c2 = "Numeric simulation")
            leyenda<-paste("r =",r,"K =",K,"a =", a,"\U1D6FD =",b,"\U1D6FE=",g, "N =",N)
            p1<-ggplot(data=yout1, aes(x=time, y=log(y1),group=1,colour = "Log(Prey)",linetype ="\U1D70F=27"))+
              geom_line()+
              geom_line(data=as.data.frame(yout2), aes(x=time, y=log(y1),group=1,colour = "Log(Prey)",linetype ="\U1D70F=13"))+
              geom_line(data=as.data.frame(yout3), aes(x=time, y=log(y1),group=1,colour = "Log(Prey)",linetype ="\U1D70F=5"))+  
              geom_line(data=as.data.frame(yout1), aes(x=time, y=y2,group=1,colour = "Predator (Puma)",linetype ="\U1D70F=27"))+
              geom_line(data=as.data.frame(yout2), aes(x=time, y=y2,group=1,colour = "Predator (Puma)",linetype ="\U1D70F=13"))+
              geom_line(data=as.data.frame(yout3), aes(x=time, y=y2,group=1,colour = "Predator (Puma)",linetype ="\U1D70F=5"))+
              theme_bw()+
              labs(y = expression("Density N°ind/100km"^2) ,x="Time (Months)", color = "Species")+
              facet_wrap(~ c2)
            p2<-ggplot(data=yout1, aes(x=y1, y=y2,colour = "\U1D70F=27" ))+
              geom_path()+
              geom_path(data=as.data.frame(yout2), aes(x=y1, y=y2, group=1, colour = "\U1D70F=13" ))+
              geom_path(data=as.data.frame(yout3), aes(x=y1, y=y2,group=1, colour = "\U1D70F=5" ))+
              theme_bw()+
              labs(y = expression("Predator N°ind/100km"^2) ,x=expression("Prey N°ind/100km"^2), color = "Delay")+
              expand_limits(x = 0,y = 0)+
              geom_label(aes(x = max(y1)/1.9 , y = max(y2)+0.25, label =leyenda), stat = "unique",color="black",size = 2.8)+
              facet_wrap(~ c1)
            png(paste0("Grafico_tau_var",contador,".png"), width = 4800, height = 2400, res=600)
            grid.arrange(p2,p1,nrow=1)
            dev.off() 
            print(contador)
            contador=contador+1
          }
        }
      }
    }
  }
}

########Simulación  con retiro de Species ###############

#Funciones retiro de Species
#Retiro de presas

rootfun <- function(t, y, p)
  return(y[1] - Ycrit)
eventfun <- function(t, y, p)
  return (c(y[1] * 0.5, y[2]))

Ycrit <- exp(4)
#Listas de parametros
r=0.05 #Tasa nacimiento presas mensual
k=200 #Capacidad de presas
a=0.8
b=0.1 #Mortaliadad mensual puma
g=0.8
N=2 #relacionado con capacidad de carga dos situaciones
tau1 <- 27 #Meses
yini <- c(y1 = 200, y2 = 1)
times <- seq(from = 0, to = 1200, by = 0.1)

LVdede <- function(t, y, p) {
  if (t > tau1) Lag1 <- lagvalue(t - tau1) else Lag1 <- yini
  dy1=r*((k-y[1])/k)*y[1]-((g*y[1]^2)*y[2])/(a^2+y[1]^2)
  dy2=-b*y[2]+((g*Lag1[1]^2)/(a^2+Lag1[1]^2))*Lag1[2]*exp(-Lag1[2]/N)
  list(c(dy1, dy2))
}

yout1 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, rootfun = rootfun,
                           events = list(func = eventfun, root = TRUE)))
yout1 <- cbind(yout1 , c1 = "Density puma and prey")
yout1 <- cbind(yout1 , c2 = "Numeric simulation")

leyenda<-paste("r =",r,"K =",k,"a =", a,"\U1D6FD=",b,"\U1D6FE =",g, "N =",N)
p1<-ggplot(data=yout1, aes(x=time, y=log(y1),group=1,colour = "Log(Prey)",linetype ="\U1D70F=27"))+
  geom_line()+
  geom_line(data=as.data.frame(yout1), aes(x=time, y=y2,group=1,colour = "Predator (Puma)",linetype ="\U1D70F=27"))+
  theme_bw()+
  labs(y = expression("Density N°ind/100km"^2) ,x="Time (Months)", color = "Species", linetype="Delay")+
  facet_wrap(~ c2)
p2<-ggplot(data=yout1, aes(x=y1, y=y2,colour = "\U1D70F=27" ))+
  geom_path()+
  theme_bw()+
  labs(y = expression("Predator N°ind/100km"^2) ,x=expression("Prey N°ind/100km"^2), color = "Delay")+
  expand_limits(x = 0,y = 0)+
  geom_label(aes(x = max(y1)/1.9 , y = max(y2)+0.25, label =leyenda), stat = "unique",color="black",size = 2.8)+
  facet_wrap(~ c1)
png("Retiro_50_presas.png", width = 4800, height = 2400, res=600)
grid.arrange(p2,p1,nrow=1)
dev.off() 

#Retiro de predador

rootfun <- function(t, y, p)
  return(y[2] - Ycrit)
eventfun <- function(t, y, p)
  return (c(y[1], y[2]*0.5))

Ycrit <- 4
#Listas de parametros
r=0.05 #Tasa nacimiento presas mensual
k=200 #Capacidad de presas
a=0.8
b=0.1 #Mortaliadad mensual puma
g=0.8
N=2 #relacionado con capacidad de carga dos situaciones
tau1 <- 27 #Meses
yini <- c(y1 = 200, y2 = 1)
times <- seq(from = 0, to = 600, by = 0.1)

LVdede <- function(t, y, p) {
  if (t > tau1) Lag1 <- lagvalue(t - tau1) else Lag1 <- yini
  dy1=r*((k-y[1])/k)*y[1]-((g*y[1]^2)*y[2])/(a^2+y[1]^2)
  dy2=-b*y[2]+((g*Lag1[1]^2)/(a^2+Lag1[1]^2))*Lag1[2]*exp(-Lag1[2]/N)
  list(c(dy1, dy2))
}

th<-0.1
yout1 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, rootfun = rootfun,
                           events = list(func = eventfun, root = TRUE)))
yout1 <- cbind(yout1 , c1 = "Density puma and prey")
yout1 <- cbind(yout1 , c2 = "Numeric simulation")
leyenda<-paste("r =",r,"K =",k,"a =", a,"\U1D6FD=",b,"\U1D6FE =",g, "N =",N)
p1<-ggplot(data=yout1, aes(x=time, y=log(y1),group=1,colour = "Log(Prey)",linetype ="\U1D70F=27"))+
  geom_line()+
  geom_line(data=as.data.frame(yout1), aes(x=time, y=y2,group=1,colour = "Predator (Puma)",linetype ="\U1D70F=27"))+
  theme_bw()+
  labs(y = expression("Density N°ind/100km"^2) ,x="Time (Months)", color = "Species", linetype="Delay")+
  facet_wrap(~ c2)
p2<-ggplot(data=yout1, aes(x=y1, y=y2,colour = "\U1D70F=27" ))+
  geom_path()+
  theme_bw()+
  labs(y = expression("Predator N°ind/100km"^2) ,x=expression("Prey N°ind/100km"^2), color = "Delay")+
  expand_limits(x = 0,y = 0)+
  geom_label(aes(x = max(y1)/1.9 , y = max(y2)+0.25, label =leyenda), stat = "unique",color="black",size = 2.8)+
  facet_wrap(~ c1)
png("Retiro_50_predador.png", width = 4800, height = 2400, res=600)
grid.arrange(p2,p1,nrow=1)
dev.off() 


##### Equilibrio positivo #########

#Funciones retiro de Species
#Retiro de presas

rootfun <- function(t, y, p)
  return(y[1] - Ycrit)
eventfun <- function(t, y, p)
  return (c(y[1] * 0.5, y[2]))

Ycrit <- exp(5)
#Listas de parametros
r=0.05 #Tasa nacimiento presas mensual
k=200 #Capacidad de presas
a=0.8
b=0.05 #Mortaliadad mensual puma
g=0.8
N=1 #relacionado con capacidad de carga dos situaciones
tau1 <- 27 #Meses
yini <- c(y1 = 200, y2 = 1)
times <- seq(from = 0, to = 1200, by = 0.1)

LVdede <- function(t, y, p) {
  if (t > tau1) Lag1 <- lagvalue(t - tau1) else Lag1 <- yini
  dy1=r*((k-y[1])/k)*y[1]-((g*y[1]^2)*y[2])/(a^2+y[1]^2)
  dy2=-b*y[2]+((g*Lag1[1]^2)/(a^2+Lag1[1]^2))*Lag1[2]*exp(-Lag1[2]/N)
  list(c(dy1, dy2))
}

yout1 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, rootfun = rootfun,
                           events = list(func = eventfun, root = TRUE)))
yout1 <- cbind(yout1 , c1 = "Density puma and prey")
yout1 <- cbind(yout1 , c2 = "Numeric simulation")

leyenda<-paste("r =",r,"K =",k,"a =", a,"\U1D6FD=",b,"\U1D6FE =",g, "N =",N)
p1<-ggplot(data=yout1, aes(x=time, y=log(y1),group=1,colour = "Log(Prey)",linetype ="\U1D70F=27"))+
  geom_line()+
  geom_line(data=as.data.frame(yout1), aes(x=time, y=y2,group=1,colour = "Predator (Puma)",linetype ="\U1D70F=27"))+
  theme_bw()+
  labs(y = expression("Density N°ind/100km"^2) ,x="Time (Months)", color = "Species", linetype="Delay")+
  facet_wrap(~ c2)
p2<-ggplot(data=yout1, aes(x=y1, y=y2,colour = "\U1D70F=27" ))+
  geom_path()+
  theme_bw()+
  labs(y = expression("Predator N°ind/100km"^2) ,x=expression("Prey N°ind/100km"^2), color = "Delay")+
  expand_limits(x = 0,y = 0)+
  geom_label(aes(x = max(y1)/1.9 , y = max(y2)+0.25, label =leyenda), stat = "unique",color="black",size = 2.8)+
  facet_wrap(~ c1)
png("Retiro_50_presas_b.png", width = 4800, height = 2400, res=600)
grid.arrange(p2,p1,nrow=1)
dev.off() 

#Retiro de predador

rootfun <- function(t, y, p)
  return(y[2] - Ycrit)
eventfun <- function(t, y, p)
  return (c(y[1], y[2]*0.5))

Ycrit <- 3
#Listas de parametros
r=0.05 #Tasa nacimiento presas mensual
k=200 #Capacidad de presas
a=0.8
b=0.05 #Mortaliadad mensual puma
g=0.8
N=1 #relacionado con capacidad de carga dos situaciones
tau1 <- 27 #Meses
yini <- c(y1 = 200, y2 = 1)
times <- seq(from = 0, to = 1200, by = 0.1)

LVdede <- function(t, y, p) {
  if (t > tau1) Lag1 <- lagvalue(t - tau1) else Lag1 <- yini
  dy1=r*((k-y[1])/k)*y[1]-((g*y[1]^2)*y[2])/(a^2+y[1]^2)
  dy2=-b*y[2]+((g*Lag1[1]^2)/(a^2+Lag1[1]^2))*Lag1[2]*exp(-Lag1[2]/N)
  list(c(dy1, dy2))
}

th<-0.1
yout1 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, rootfun = rootfun,
                           events = list(func = eventfun, root = TRUE)))
yout1 <- cbind(yout1 , c1 = "Density puma and prey")
yout1 <- cbind(yout1 , c2 = "Numeric simulation")
leyenda<-paste("r =",r,"K =",k,"a =", a,"\U1D6FD=",b,"\U1D6FE =",g, "N =",N)
p1<-ggplot(data=yout1, aes(x=time, y=log(y1),group=1,colour = "Log(Prey)",linetype ="\U1D70F=27"))+
  geom_line()+
  geom_line(data=as.data.frame(yout1), aes(x=time, y=y2,group=1,colour = "Predator (Puma)",linetype ="\U1D70F=27"))+
  theme_bw()+
  labs(y = expression("Density N°ind/100km"^2) ,x="Time (Months)", color = "Species", linetype="Delay")+
  facet_wrap(~ c2)
p2<-ggplot(data=yout1, aes(x=y1, y=y2,colour = "\U1D70F=27" ))+
  geom_path()+
  theme_bw()+
  labs(y = expression("Predator N°ind/100km"^2) ,x=expression("Prey N°ind/100km"^2), color = "Delay")+
  expand_limits(x = 0,y = 0)+
  geom_label(aes(x = max(y1)/1.9 , y = max(y2)+0.25, label =leyenda), stat = "unique",color="black",size = 2.8)+
  facet_wrap(~ c1)
png("Retiro_50_predador_b.png", width = 4800, height = 2400, res=600)
grid.arrange(p2,p1,nrow=1)
dev.off() 

