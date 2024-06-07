##Function Stability analysis (SFA)##

#Librerias requeridas
library(rootSolve)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(deSolve)
options(warn=-1)

FSA <- function(r,K,a,b,g,N,x_0,y_0,tau1,simu_time) {
  #Calcular Puntos de Equilibrio
  equations <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    f1 <- r*((K-x)/K)*x-((g*x^2)*y)/(a^2+x^2)
    f2 <- -b*y+((g*x^2)/(a^2+x^2))*y*exp(-y/N)
    return(c(f1, f2))
  }
  initial_guess <- c(x = K, y = 0)  # Punto de inicio
  result <- multiroot(equations, start = initial_guess)
  eq_x<-round(result$root[1],2)
  eq_y<-round(result$root[2],2)
  
  initial_guess <- c(x = 0, y = 0)  # Punto de inicio
  result <- multiroot(equations, start = initial_guess,maxiter=999,positive=TRUE,useFortran=FALSE)  ####
  eq_x_0<-round(result$root[1],2)
  eq_y_0<-round(result$root[2],2)
  
  for (i in 1:999) {
    initial_guess <- c(x = K, y = i)  # Punto de inicio
    result <- multiroot(equations, start = initial_guess,maxiter=999,positive=TRUE,useFortran=FALSE)
    eq_x_temp<-round(result$root[1],2)
    eq_y_temp<-round(result$root[2],2)
    if (eq_x_temp > 0 & eq_y_temp >0) {
      break
    }
  }
  
  eq_x_1<-eq_x_temp
  eq_y_1<-eq_y_temp
  
  for (j in 1:999) {
    initial_guess <- c(x =j+10, y = j)  # Punto de inicio
    result <- multiroot(equations, start = initial_guess,maxiter=99999,positive=TRUE,useFortran=FALSE)
    eq_x_temp1<-round(result$root[1],2)
    eq_y_temp1<-round(result$root[2],2)
    if (eq_x_temp1 > 0 & eq_y_temp1 >0) {
      break
    }
  }
  eq_x_2<-eq_x_temp1
  eq_y_2<-eq_y_temp1
  
  for (j in 1:K/2) {
    initial_guess <- c(x =K/2-j, y = K/2-j)  # Punto de inicio
    result <- multiroot(equations, start = initial_guess,maxiter=9999,positive=TRUE,useFortran=FALSE)
    eq_x_temp1<-round(result$root[1],2)
    eq_y_temp1<-round(result$root[2],2)
    if (eq_x_temp1 > 0 & eq_y_temp1 >0) {
      break
    }
  }
  eq_x_3<-eq_x_temp1
  eq_y_3<-eq_y_temp1
  
  cond<-round((g*K^2)/(b*(a^2+K^2)),2)
  
  lista_tau<-c()
  c_tau_1<-0
  
  #Condición 3.8
  print("--------------------- Reporte de Puntos de equilibrio ------------------------")
  if(cond>1){
    print("Por la condición 3.8 al menos hay una solución positiva del sistema")
  }else {
    print("Por la condición 3.8 no hay soluciones positivas del sistema")
  }
  print("-------------------------------------------------------------------------------")
  
  #Analisis de estabilidad
  if(eq_y_0==0 & eq_x_0==0){
    print("Existe el punto de equilibrio es del tipo (0,0)")
    print("El punto de equilibrio es inestable")
    print("-------------------------------------------------------------------------------")
  }
  
  if(eq_y==0 & eq_x>0){
    m<- r+b
    n<- -((g*K^2)/(a^2+K^2))
    p<- r*b  
    q<- -((g*r*K^2)/(a^2+K^2))
    print("Existe el punto de equilibrio es del tipo (K,0)")
    print(paste("El valor de equilibrio las presas es",eq_x))
    print(paste("El valor de equilibrio del depredador es",eq_y))
    if (b>g){
      print("(K, 0) del sistema es asintóticament globalmente estable")}
    else if ( b>(g*K^2)/(a^2+K^2)){
      print("Se cumple con condición H1")
      print("(K, 0) del sistema es asintóticamente estable para todo τ ≥ 0")}
    else{print('El criterio no entrega información')}
    print("-------------------------------------------------------------------------------")
  }
  
  if(eq_y_1>0 & eq_x_1>0){
    print("Existe el punto de equilibrio es del tipo (X*,Y*)")
    print(paste("El valor de equilibrio las presas es",eq_x_1))
    print(paste("El valor de equilibrio del depredador es",eq_y_1))
    if( ((r*eq_x_1/K)+b*(eq_y_1/N))/(r*(1-(eq_x_1/K)))>(eq_x_1^2-a^2)/(eq_x_1^2+a^2) ){
      print("Se cumple con condición H2")
      if(1+(1+(eq_y_1/N)*(r*(1-(eq_x_1/K))*(((2*a^2)/(a^2+eq_x_1^2))-1)+((r*eq_x_1)/(K))))>
         ((2*g*a^2*eq_x_1^2)/((a^2+eq_x_1^2)^2))*r*(1-(eq_x_1/K))){
        print("Se cumple con condición H3")
        
        m<- r*(1-(eq_x_1/K))*((a^2-eq_x_1^2)/(a^2+eq_x_1^2)) + b + (r*eq_x_1/K)
        n<- b*((eq_y_1/N)-1)
        p<- b*(m-b)
        q<- b*((1-(eq_y_1/N))*(m-b)-exp(eq_y_1/N)*(2*r*a^2/(a^2+eq_x_1^2))*(1-(eq_x_1/K)))
        
        if (2*p+n^2-m^2<0 & p^2-q^2>0){
          print("(x*, y*) del sistema es asintóticamente estable para todo τ ≥ 0")
        }
        else if (p^2-q^2<0){
          print("(x*, y*) es asintóticamente estable para τ ∈ [0, τ_0+), e inestable para τ ∈(τ_0+ , +∞)")
          w_1<-sqrt(0.5*(2*p+n^2-m^2)+0.5*((2*p+n^2-m^2)^2-4*(p^2-q^2))^(1/2)) 
          c_tau_1<-round((1/w_1)*acos((q*(w_1^2-p)-m*n*w_1^2)/(q^2+n^2*w_1^2))+(2*3.1416/(w_1)),2)
          print( paste("Con τ_0+ =", c_tau_1))
          lista_tau<-append(lista_tau, c_tau_1 )
        }
        else if (p^2-q^2>0 & n^2-m^2+2*p>0 & (n^2-m^2+2*p)^2>4*(p^2-q^2)) {
          print("Existe un entero positivo k tal que se presentan k cambios de estabilidad")
          c_tau_1<-0.2*tau1
          c_tau_2<-0.8*tau1
          lista_tau<-append(lista_tau, c_tau_1 )
          lista_tau<-append(lista_tau, c_tau_2 )
        }
        else{print('El criterio no entrega información')}
      }else {print("No se cumple con condición H3")}
    }else{print("No se cumple con condición H2")}}
  
  if((eq_y_2>0 & eq_x_2>0)&(eq_y_2!=eq_y_1 | eq_x_2!=eq_x_1)){
    print("-------------------------------------------------------------------------------")
    print("Existe un segundo punto de equilibrio es del tipo (X*,Y*)")
    print(paste("El valor de equilibrio las presas es",eq_x_2))
    print(paste("El valor de equilibrio del depredador es",eq_y_2))
    if( ((r*eq_x_2/K)+b*(eq_y_2/N))/(r*(1-(eq_x_2/K)))>(eq_x_2^2-a^2)/(eq_x_2^2+a^2) ){
      print("Se cumple con condición H2")
      if(1+(1+(eq_y_2/N)*(r*(1-(eq_x_2/K))*(((2*a^2)/(a^2+eq_x_2^2))-1)+((r*eq_x_2)/(K))))>
         ((2*g*a^2*eq_x_2^2)/((a^2+eq_x_2^2)^2))*r*(1-(eq_x_2/K))){
        print("Se cumple con condición H3")
        
        m<- r*(1-(eq_x_2/K))*((a^2-eq_x_2^2)/(a^2+eq_x_2^2)) + b + (r*eq_x_2/K)
        n<- b*((eq_y_2/N)-1)
        p<- b*(m-b)
        q<- b*((1-(eq_y_2/N))*(m-b)-exp(eq_y_2/N)*(2*r*a^2/(a^2+eq_x_2^2))*(1-(eq_x_2/K)))
        
        if(2*p+n^2-m^2<0 & p^2-q^2>0){
          print("(x*, y*) del sistema es asintóticamente estable para todo τ ≥ 0")
        }
        else if(p^2-q^2<0){
          print("(x*, y*) es asintóticamente estable para τ ∈ [0, τ_0+), e inestable para τ ∈(τ_0+ , +∞)")
          w_1<-sqrt(0.5*(2*p+n^2-m^2)+0.5*((2*p+n^2-m^2)^2-4*(p^2-q^2))^(1/2)) 
          c_tau_1<-round((1/w_1)*acos((q*(w_1^2-p)-m*n*w_1^2)/(q^2+n^2*w_1^2))+(2*3.1416/(w_1)),2)
          print( paste("Con τ_0+ =", c_tau_1))
          lista_tau<-append(lista_tau, c_tau_1 )
        }
        else if(p^2-q^2>0 & n^2-m^2+2*p>0 & (n^2-m^2+2*p)^2>4*(p^2-q^2)) {
          print("Existe un entero positivo k tal que se presentan k cambios de estabilidad")
          c_tau_1<-1.2*tau1
          c_tau_2<-1.5*tau1
          lista_tau<-append(lista_tau, c_tau_1 )
          lista_tau<-append(lista_tau, c_tau_2 )}
          
        else{print('El criterio no entrega información')}
     }else {print("No se cumple con condición H3")}
    }else {print("No se cumple con condición H2")}}
  
  if((eq_y_3>0 & eq_x_3>0)&((eq_y_3!=eq_y_2 & eq_x_3!=eq_x_1) | (eq_x_3!=eq_x_2 &eq_y_3!=eq_y_1) )){
    print("-------------------------------------------------------------------------------")
    print("Existe un tercer punto de equilibrio es del tipo (X*,Y*)")
    print(paste("El valor de equilibrio las presas es",eq_x_3))
    print(paste("El valor de equilibrio del depredador es",eq_y_3))
    if( ((r*eq_x_3/K)+b*(eq_y_3/N))/(r*(1-(eq_x_3/K)))>(eq_x_3^2-a^2)/(eq_x_3^2+a^2) ){
      print("Se cumple con condición H2")
      if(1+(1+(eq_y_3/N)*(r*(1-(eq_x_3/K))*(((2*a^2)/(a^2+eq_x_3^2))-1)+((r*eq_x_3)/(K))))>
         ((2*g*a^2*eq_x_3^2)/((a^2+eq_x_3^2)^2))*r*(1-(eq_x_3/K))){
        print("Se cumple con condición H3")
        
        m<- r*(1-(eq_x_3/K))*((a^2-eq_x_3^2)/(a^2+eq_x_3^2)) + b + (r*eq_x_3/K)
        n<- b*((eq_y_3/N)-1)
        p<- b*(m-b)
        q<- b*((1-(eq_y_3/N))*(m-b)-exp(eq_y_3/N)*(2*r*a^2/(a^2+eq_x_3^2))*(1-(eq_x_3/K)))
        
        if(2*p+n^2-m^2<0 & p^2-q^2>0){
          print("(x*, y*) del sistema es asintóticamente estable para todo τ ≥ 0")
        }
        else if(p^2-q^2<0){
          print("(x*, y*) es asintóticamente estable para τ ∈ [0, τ_0+), e inestable para τ ∈(τ_0+ , +∞)")
          w_1<-sqrt(0.5*(2*p+n^2-m^2)+0.5*((2*p+n^2-m^2)^2-4*(p^2-q^2))^(1/2)) 
          c_tau_1<-round((1/w_1)*acos((q*(w_1^2-p)-m*n*w_1^2)/(q^2+n^2*w_1^2))+(2*3.1416/(w_1)),2)
          print(paste("Con τ_0+ =", c_tau_1))
          lista_tau<-append(lista_tau, c_tau_1)
          
        }
        else if(p^2-q^2>0 & n^2-m^2+2*p>0 & (n^2-m^2+2*p)^2>4*(p^2-q^2)) {
          print("Existe un entero positivo k tal que se presentan k cambios de estabilidad")
          c_tau_1<-1.5*tau1
          c_tau_2<-1.8*tau1
          lista_tau<-append(lista_tau, c_tau_1 )
          lista_tau<-append(lista_tau, c_tau_2 )}
          
        else{print('El criterio no entrega información')}
      }else {print("No se cumple con condición H3")}
    }else {print("No se cumple con condición H2")}}
  
  #Simulación y gráficas
  LVdede <- function(t, y, p, tau) {
    if (t > tau) Lag1 <- lagvalue(t - tau) else Lag1 <- yini
    dy1=r*((K-y[1])/K)*y[1]-((g*y[1]^2)*y[2])/(a^2+y[1]^2)
    dy2=-b*y[2]+((g*Lag1[1]^2)/(a^2+Lag1[1]^2))*Lag1[2]*exp(-Lag1[2]/N)
    list(c(dy1, dy2))
  }
  yini <- c(y1 = x_0, y2 = y_0)
  times <- seq(from = 0, to = simu_time, by = 0.1)
  
  Funtion_plot<-function(tau_list){
    
    tau_list<-round(tau_list,1)
    
    tau<-tau_list[1]
    yout1<- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, tau=tau))
    yout1 <- cbind(yout1 , c1 = "Density puma and prey")
    yout1 <- cbind(yout1 , c2 = "Numeric simulation")
    
    tau <-tau_list[2]
    yout2<- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, tau=tau))
    yout2 <- cbind(yout2 , c1 = "Density puma and prey")
    yout2 <- cbind(yout2 , c2 = "Numeric simulation")
    
    tau <-tau_list[3]
    yout3 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, tau=tau))
    yout3 <- cbind(yout3 , c1 = "Density puma and prey")
    yout3 <- cbind(yout3 , c2 = "Numeric simulation")
    
    tau <-tau_list[4]
    yout4 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, tau=tau))
    yout4 <- cbind(yout4 , c1 = "Density puma and prey")
    yout4 <- cbind(yout4 , c2 = "Numeric simulation")
    
    tau <-tau_list[5]
    yout5 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, tau=tau))
    yout5 <- cbind(yout5 , c1 = "Density puma and prey")
    yout5 <- cbind(yout5 , c2 = "Numeric simulation")
    
    tau <-tau_list[6]
    yout6 <- as.data.frame(dede(func = LVdede, y = yini, times = times, parms = 0, tau=tau))
    yout6 <- cbind(yout6 , c1 = "Density puma and prey")
    yout6 <- cbind(yout6 , c2 = "Numeric simulation")
    
    leyenda<-paste("Parameters","\n", "r =",r,"K =",K,"\n", "a =", a,"β =",b,"\n"," γ=",g, "N =",N,"\n",
                   'Initial =(', x_0,',',y_0,')')
    
    yout1$tau <- paste("\U1D70F=",tau_list[1])
    yout2$tau <- paste("\U1D70F=",tau_list[2])
    yout3$tau <- paste("\U1D70F=",tau_list[3])
    yout4$tau <- paste("\U1D70F=",tau_list[4])
    yout5$tau <- paste("\U1D70F=",tau_list[5])
    yout6$tau <- paste("\U1D70F=",tau_list[6])
    
    yout1$tau <- factor(yout1$tau, levels = paste("\U1D70F=", tau_list))
    
    # Combinar los dataframes hacia abajo (row-binding)
    merged_yout <- rbind(yout1, yout2, yout3, yout4, yout5, yout6)
    
    p4<-ggplot(data=merged_yout, aes(x=time, y=log(y1),group=1, color = "Log(Prey)"))+
      geom_line()+
      geom_line(data=merged_yout, aes(x=time, y=y2,group=1, color= "Predator (Puma)"))+
      labs(y = expression("Density N°ind/100km"^2) ,x="Time (Months)",color  = "Species")+
      theme_bw()+
      scale_x_continuous(limits = c(0, simu_time),
                         breaks = c(0, simu_time/2, simu_time))+
      facet_wrap(~ tau)
    
    p3 <- ggplot(data = merged_yout, aes(x = y1, y = y2)) +
      geom_path(color = "red") +
      theme_bw() +
      labs(y = expression("Predator N°ind/100km"^2),x = expression("Prey N°ind/100km"^2),
           shape = "Equilibrium") +
      geom_point(aes(x = eq_x, y = eq_y, shape = "(K,0)")) +
      geom_point(aes(x = eq_x_0, y = eq_y_0, shape = "(0,0)")) +
      expand_limits(x = 0, y = 0) +
      facet_wrap(~tau)
    
    if(eq_y_1>0 & eq_x_1>0){
      p3<-p3+geom_point(aes(x=eq_x_1, y=eq_y_1, shape="(x*,y*)"))}
    if((eq_y_2>0 & eq_x_2>0)&(eq_y_2!=eq_y_1 | eq_x_2!=eq_x_1)){
      p3<-p3+geom_point(aes(x=eq_x_2, y=eq_y_2, shape="(x*,y*)"))}
    if((eq_y_3>0 & eq_x_3>0)&((eq_y_3!=eq_y_2 & eq_x_3!=eq_x_1) | (eq_x_3!=eq_x_2 &eq_y_3!=eq_y_1) )){
      p3<-p3+geom_point(aes(x=eq_x_3, y=eq_y_3, shape="(x*,y*)"))}
    
    p3<-p3+labs(tag=leyenda)+
      coord_cartesian(clip = "off") +
      theme(plot.tag.position = c(0.9, 0.8), plot.tag = element_text(size=9))
    grid.arrange(p3,p4,nrow=1)
  }
  
  if(c_tau_1==0){
    lista<-seq(min(tau1)*0.2, max(tau1)*1.2, length.out = 6)
    Funtion_plot(tau_list= lista)}
  if(length(lista_tau)==1){
    lista<-seq(min(c_tau_1)*0.2, max(c_tau_1)*1.2, length.out = 6)
    Funtion_plot(tau_list= lista)}
  if(length(lista_tau)>=2){
    lista<-seq(min(lista_tau)*0.2, max(lista_tau)*1.2, length.out = 6)
    Funtion_plot(tau_list= lista)
  }
}

#NOTA: Toma el tau inicial (tau1) como referencia cuando 
#      hay estabilidad absoluta para realizar el grafico

#Ejemplo de como aplicar la función 

#Example 1 es asintóticamente estable para τ ∈ [0, τ_0+)
png("Example1.png", width = 6000, height = 2400, res=600)
FSA(r=0.1,K=300,a=0.8,b=0.05,
    g=0.8,N=1,x_0=1,y_0=1,tau1=35, 1200)
dev.off()

#Example 2 (x*, y*) asintóticamente estable para todo τ ≥ 0"
png("Example2.png", width = 6000, height = 2400, res=600)
FSA(r=0.08,K=200,a=0.05,b=0.5,
    g=0.8,N=2,x_0=250,y_0=3,tau1=30, 2400)
dev.off()

#Example 3 (K,0) asintoticamente globalmente estable
png("Example3.png", width = 6000, height = 2400, res=600)
FSA(r=0.05,K=50,a=0.1,b=0.1,
    g=0.08,N=2,x_0=5,y_0=1,tau1=50, 1200)
dev.off()

#Example 4 (x*, y*) asintóticamente estable para τ ∈ [0, τ_0+), e inestable para τ ∈(τ_0+ , +∞)
#(x*, y*) asintóticamente estable para todo τ ≥ 0
png("Example4.png", width = 6000, height = 2400, res=600)
FSA(r=0.2,K=100,a=0.1,b=0.2,
    g=0.5,N=2,x_0=0.5,y_0=0.5,tau1=5,1200)
dev.off()

#Example 5
png("Example5.png", width = 6000, height = 2400, res=600)
FSA(r=0.2,K=300,a=0.8,b=0.05,
    g=0.9,N=5,x_0=150,y_0=12,tau1=1, 1200)
dev.off()






