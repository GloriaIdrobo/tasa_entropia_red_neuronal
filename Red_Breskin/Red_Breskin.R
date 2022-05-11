#_______________________________Creacion,analisis topologico y Calculo de h para la Red Breskin_______________________________

#Librerias 
#install.packages("igraph") 
#install.packages("repr")
#install.packages("SciViews")
library("igraph")
library("repr")
library("SciViews")

#Red neuronal
Nodos=2000

#Creamos un vector donde se guardaran los números aleatorios de las posiciones
Posicionx<-seq(1,Nodos)
Posiciony<-seq(1,Nodos)

for (i in 1:Nodos) {
  x<-runif(1, min=0, max=2) # Crea los numeros aleatorios entre 0 y 2 que reprensenta la posicion de las neuronas 
  y<-runif(1, min=0, max=2)
  Posicionx[i]<- x # Guarda el número aleatorio en el vector
  Posiciony[i]<- y
}
#Vectores de posicion de las neuronas
Posicionx
Posiciony

#Graficia de ubicacion aleatoria de las neuronas en la rejilla
plot(Posicionx, Posiciony, xlab="X(mm)", ylab="Y(mm)",pch=20, col="red")
#_________________________________________________________________________________________________________________

#Distancia entre cada una de las neuronas hacia las demas
Lij <- matrix(data=0, nrow=Nodos, ncol=Nodos)
for(i in 1:Nodos){
  for (j in 1:Nodos) {
    DISTANCIAiAj <-sqrt((Posicionx[j]-Posicionx[i])^2 + (Posiciony[j]-Posiciony[i])^2)
    Lij[i,j]<-DISTANCIAiAj
  }
}
Lij
#_________________________________________________________________________________________________________________

#funcion de probabildad que une las neuronas
ProConexion<-matrix(data=0, nrow=Nodos, ncol=Nodos)
pijneuronas<-0
for (i in 1:Nodos) {
  for (j in 1:Nodos) {
    Lsubij<-Lij[i,j]
    Lsubij
    if ((000.1<Lsubij) & (Lsubij<= 1)) pijneuronas=(1/(exp((Lsubij/0.3)-(1.2/Lsubij))+1)) else  pijneuronas <- 0
    ProConexion[i,j]<- pijneuronas
  }
}
ProConexion[]

#Creacion de neuronas con data.frame
Neuronasposiciones <- data.frame( Neuronas=1:Nodos,
                                  PosicionX=Posicionx,
                                  PosicionY=Posiciony )
Neuronasposiciones 
#_________________________________________________________________________________________________________________

#Creación de la matriz de adyacencia segun la probabilidad 
MatrizAdyacencia<-matrix(data=0, nrow=Nodos, ncol=Nodos)

for (i in 1:Nodos) {
  for (j in 1:Nodos) {
    if (ProConexion[i,j]==0) MatrizAdyacencia[i,j]=0 else  MatrizAdyacencia[i,j]=1
  }
}
MatrizAdyacencia[]
write.table(MatrizAdyacencia, "Matriz_A_Breskin.txt") 

red_B <- graph_from_adjacency_matrix(MatrizAdyacencia)

#________________________________PROPIEDADES TOPOLÓGICAS DE LA RED______________________________________________

k<-degree(red_B, mode=c("all")) #Grado de cada nodo
hist(k, breaks = 30) # Histograma del Grado
Grado_medio <- mean(k)
Grado_medio
gsize(red_B) # Cantidad de enlaces que tiene la red
gd <- edge_density(red_B) #Densidad de la red= a proporción de enlaces presentes de todos los enlaces posibles en la red
diameter(red_B, directed=FALSE) # Diametro de la red
dis<-mean_distance(red_B, directed = F) # Longitud promedio de camino más corto
dis
Camino_corto<-get_diameter(red_B) #Nodos de la longitud promedio de camino más corto. Figura 3.6 libro
Camino_corto
#Codigo que muestra solo los nodos que creran la lngitud de camino más corto en la red
longitud_camino<-length(Camino_corto)
Posicionx_camino_corto<-c()
Posiciony_camino_corto<-c()
for (i in 1:longitud_camino) {
  s<-Camino_corto[i]
  posix<-Posicionx[s]
  posiy<-Posiciony[s]
  Posicionx_camino_corto[i]<-posix
  Posiciony_camino_corto[i]<-posiy
}
Posicionx_camino_corto

#Graficia los nodos del camino corto
plot(Posicionx_camino_corto, Posiciony_camino_corto,xlab = "X (mm)",ylab = "Y (mm)", pch = 20,tcl=0.2, color = "Pink",las=1)

centralidad<-centr_clo(red_B, mode="all", normalized=T)# Centralidad de cercania
nodo_maxima_centralidad<-which.max(closeness)
nodo_maxima_centralidad

bett<-betweenness(red_B, directed=T)
nodo_maxima_centralidad_intermedia<-which.max(bett) #Centralidad de intermediación
nodo_maxima_centralidad_intermedia

Wc <- cluster_walktrap(red_B) 
modularity(Wc) #Modularidad en la red

#________________________________________________Análisis tasa de entropía________________________________________
MATRIZADY<-MatrizAdyacencia # Matriz 
k<-degree(red_B) #Grado de cada nodo
Num_nodos<-Nodos #Cantidad de nodos en la red

#Se crea la matriz de transición de probabilidad PI vacia
PI <- matrix(data=0, nrow=Num_nodos, ncol=Num_nodos)

#Se crea un vector que almacenará los valores de la tasa de entropia para cada alpha de -3 a 3
H<-seq(1,61) 

#Valores del sesgo en los que se va a evaluar la tasa de entropía
ValoresAlpha<- seq(1,61) #se va a evaluar de alpha -3 a alpha 3 con intervalos de 0.1 es decir 60 valores de alpha
Alphapasos<--3.1
for (i in 1:61) {
  Alphapasos<-Alphapasos+0.1
  ValoresAlpha[i]<-Alphapasos
}
ValoresAlpha[] 

#Calculo de la ecuación de la tasa de entropía para 61 valores de alpha
for (alph in 1:61)
{
  alpha <- ValoresAlpha[alph]  #En cada paso alpha cambia de valor para calcular h
  Ci <-seq(1,Num_nodos) # Se crea el vector Ci que guardara el valor de ci para cada i 
  
  #Calculo de la matriz de transición de probabilidad PI
  # ______________________________________________________________________________________
  for (i in 1:Num_nodos)
  { 
    denominadortotal <-0 
    #Bucle para calcular la parte del denominador de la ecuación de PI para cada i
    for (s in 1:Num_nodos)
    {
      denominador<-MATRIZADY[i,s]*(k[s])^alpha
      denominadortotal <- denominador+ denominadortotal
    }
    denominadortotal # valor del denominador para un i
    
    #Bucle para calcular la parte del numerador de la ecuación PI
    for (j in 1:Num_nodos)
    {
      Ady<-MATRIZADY[i,j] # Valor a_{i,j} de la matriz de adyacencia
      gr<-k[j]      # Valor del grado que tiene el nodo j, es decir el nodo destino
      numerador <- Ady*(gr)^alpha
      pij<- numerador/denominadortotal # calculo de la probabilidad de transición de i a j
      PI[i,j] <- pij # Asigna la probabilidad de transición de i a j en la matriz PI
    }
    
    #parte de w*
    ci<-0
    #Bucle para calcular el ci
    for (j in 1:Num_nodos)
    {
      c<-MATRIZADY[i,j]*(k[j])^alpha
      ci <- ci+c
    }
    Ci[i]<-ci
    
  }
  # PI # Mostrar matriz de transición 
  #Ci # C_{i} 
  #_______________________________________________________________________________
  # Distribución estacionaria W*
  #Ahora vamos a calcular el numerador de la ecuación de distribución estacionaria W*
  Numeradori <-seq(1,Num_nodos) # Guardará el valor del numerador de cada i
  
  #Bucle para calcular el numerador
  for (i in 1:Num_nodos)
  {
    nu<-Ci[i]*(k[i])^alpha
    Numeradori[i] <-nu
  }
  Numeradori
  Denominadorsuma<-sum(Numeradori)
  
  #RESULTADO DE W*
  #Bucle para calcular el W*
  Wstar<-seq(1,Num_nodos) # Creación del vector que guardara la información de cada W*_{i}
  for (i in 1:Num_nodos)
  {
    wi<-Numeradori[i]/Denominadorsuma
    Wstar[i]<-wi
  }
  # Wstar # Mostrar la DISTRIBUCIÓN ESTACIONARIA W*
  
  #__________________________________________________________________________________________
  #Calculo de la tasa de entropía
  # H = PI_{i,j} W* ln PI{i,j}
  Ensayoln<- ln(25)
  h<- c(0,0,0,0,0) # Vector que guardara el valor de la tasa de entropia para cada i
  # Creamos el ciclo sobre i para encontrar cada h_{i} y luego sumarlos. Ver libro
  for (i in 1:Num_nodos)
  {
    mul<-0
    for (j in 1:Num_nodos)
    {
      mmmm<-PI[i,j]
      # cuando la matriz tiene PI_{i,j}=0 el valor de h se hace cero. ver ecuación. si no calcula el logaritmo natural de PI
      if (mmmm <= 0) PIlnPIW <-0 else logaritmoN <- ln(PI[i,j])
      
      PIlnPIW<--PI[i,j]*Wstar[i]*logaritmoN
      mul<-mul+PIlnPIW
    }
    h[i]<-mul
  }
  H[alph]<-sum(h) # Valor total de la tasa de entropía para alpha
  #_____________________________________________________________________________________________
  #Creación de un txt para guardar los valores de la tasa de entropía en funcion de Alpha
  
  ALPHAvsH <- data.frame(ValoresAlpha,H )
  write.table(ALPHAvsH, "h(alpha)_red_Breskin.txt")  
}
ALPHAvsH

plot(ValoresAlpha,H)


