#_______________________________Creacion, análisis topológico y Calculo de h en funcion de alpha para la red Hellwig_______________________________
#Para instalar los paquetes elimine el numeral(#) y corra el codigo
#install.packages("igraph") 
#install.packages("repr")
#install.packages("SciViews")
#install.packages("scatterplot3d")

library("readxl")   
library("SciViews")   
library("igraph")   
library("repr")
library(scatterplot3d)

#Constantes
#Red neuronal
Nodos=2000
Nodoscapa<-Nodos/2 # Cantidad de nodos que va haber en cada capa

#Caso 2->2
A_proba<-459.3
wproba<-489.9
y_0<- 0.0714
x_proba<- -35.33
Aproba<-A_proba/(wproba*sqrt(pi/2))

#Caso 2->3
A_proba23<-987.1
wproba23<-620.8
y_023<- -0.0231
x_proba23<- -259.4
Aproba23<-A_proba/(wproba23*sqrt(pi/2))

#Caso 3->2
A_proba32<-325.7
wproba32<-451.0
y_032<- -0.0121
x_proba32<- -72.86
Aproba32<-A_proba/(wproba32*sqrt(pi/2))

#Caso 3->3
A_proba33<-270.0
wproba33<-296.7
y_033<- 0.0720
x_proba33<- 5.705
Aproba33<-A_proba/(wproba33*sqrt(pi/2))

#________________________________CREACION DE NODOS ALEATORIOS EN LA RED________________________________________

#Creamos un vector donde se guardaran los números aleatorios de las posiciones
Posicionx<-seq(1,Nodos)
Posiciony<-seq(1,Nodos)
Posicionz<-seq(1,Nodos)

#Creamos el un bucle para generar las posiciones aleatorias

for (i in 1:Nodoscapa) {
  xcapa1<-runif(1, min=0, max=2) # Crea las posiciones aleatorios entre 0 y 2 donde estaran las neuronas de la capa 1
  ycapa1<-runif(1, min=0, max=2)
  zcapa1<-runif(1, min=0, max=1) # Crea las posiciones aleatorios entre 0 y 1 donde estaran las neuronas de la capa 1
  xcapa2<-runif(1, min=0, max=2) 
  ycapa2<-runif(1, min=0, max=2)
  zcapa2<-runif(1, min=1, max=2) # Crea las posiciones aleatorios entre 0 y 2  donde estaran las neuronas de la capa2
  Posicionx[i]<- xcapa1 # Guarda el número aleatorio en el vector
  Posiciony[i]<- ycapa1
  Posicionz[i]<- zcapa1
  #Un solo vector para crear la caja de 2000 nodos
  Posicioncapa2<-Nodoscapa+i
  Posicionx[Posicioncapa2]<- xcapa2 # Guarda el número aleatorio en el vector
  Posiciony[Posicioncapa2]<- ycapa2
  Posicionz[Posicioncapa2]<- zcapa2
}

#Graficia de ubicacion aleatoria de las neuronas en el cubo
scatterplot3d(Posicionx, Posiciony, Posicionz,xlab = "X (mm)",ylab = "Y (mm)",zlab = "Z (mm)", pch = 20,tcl=0.2, color = c("Pink", "skyblue")[1+(Posicionz>1)],las=1)

#Creacion data.frame donde muestra la posicion de cada una de las neuronas
#definir capa de cada neurona
CAPA1 <-rep(1,Nodoscapa)
Caja_con_capas<-seq(1:Nodos)
Caja_con_capas[1:Nodoscapa] <- 1
Caja_con_capas[Nodoscapa:Nodos] <- 2

Neuronasposiciones <- data.frame( Neuronas=1:Nodos,
                                  Capa=Caja_con_capas,
                                  PosicionX=Posicionx,
                                  PosicionY=Posiciony)
Neuronasposiciones

#______________PROCEDIMIENTO PARA CONECTAR LAS NEURONAS SEGUN PARAMETROS EXPERIMENTALES__________________________________

#Distancia entre cada una de las neuronas hacia las demas
Lij <- matrix(data=0, nrow=Nodos, ncol=Nodos)

for(i in 1:Nodos){
  for (j in 1:Nodos) {
    DISTANCIAiAj <-sqrt((Posicionx[j]-Posicionx[i])^2 + (Posiciony[j]-Posiciony[i])^2 + (Posicionz[j]-Posicionz[i])^2 )
    Lij[i,j]<-DISTANCIAiAj
  }
}

Lij
#________________________Probabilidad de conexion entre el par de neuronas_______________________

ProConexion<-matrix(data=0, nrow=Nodos, ncol=Nodos)
pijneuronas<-0

for (i in 1:Nodos) {
  for (j in 1:Nodos) {
    Lsubij<-Lij[i,j]*1000 # pasa las longitudes en micrometros
    Lsubij
    if ((Posicionz[i]<1) & (Posicionz[j]<1))#Caso 2->2
      if ((150<Lsubij) & (Lsubij<=500)) pijneuronas=( y_0 + Aproba*(exp((-2*((Lsubij - x_proba)^2))/( wproba)^2)))  else  pijneuronas <- 0
    else  (if ((Posicionz[i]<1) & (Posicionz[j]>1)) #Caso2->3
      if ((150<Lsubij) & (Lsubij<=500)) pijneuronas=( y_023 + Aproba23*(exp((-2*((Lsubij - x_proba23)^2))/( wproba23)^2))) else  pijneuronas <- 0
      else  (if ((Posicionz[i]>1) & (Posicionz[j]>1)) #Caso3->3
        if ((150<Lsubij) & (Lsubij<=500)) pijneuronas=( y_033 + Aproba33*(exp((-2*((Lsubij - x_proba33)^2))/( wproba33)^2))) else  pijneuronas <- 0
        else  (if ((150<Lsubij) & (Lsubij<=500)) pijneuronas=( y_032 + Aproba32*(exp((-2*((Lsubij - x_proba32)^2))/( wproba32)^2))) else  pijneuronas <- 0
               
        )
      )
    )
    ProConexion[i,j]<- pijneuronas
  }
}
ProConexion
write.table(ProConexion, "ProConexion_pares_de_neuronas.txt")

#___________________________CREACION DE LA MATRIZ DE ADYACENCIA SEGÚN LA PROBABILIDAD DE CONEXION __________________

MAdyaW<-matrix(data=0, nrow=Nodos, ncol=Nodos)

for (i in 1:Nodos) {
  for (j in 1:Nodos) {
    pcone<-ProConexion[i,j]
    if (ProConexion[i,j]<=0) MAdyaW[i,j]=0 else  MAdyaW[i,j]=pcone
  }
}

#Creación de la matriz de adyacencia de 1 y ceros
MAya_NO_W<-matrix(data=0, nrow=Nodos, ncol=Nodos)

for (i in 1:Nodos) {
  for (j in 1:Nodos) {
    if (ProConexion[i,j]==0) MAya_NO_W[i,j]=0 else  MAya_NO_W[i,j]=1
  }
}
MAya_NO_W[]

red_H <- graph_from_adjacency_matrix(MAya_NO_W)

#________________________________PROPIEDADES TOPOLÓGICAS DE LA RED______________________________________________

k<-degree(red_H, mode=c("all")) #Grado de cada nodo
hist(k, breaks = 30) # Histograma del Grado
Grado_medio <- mean(k)
Grado_medio
gsize(red_H) # Cantidad de enlaces que tiene la red
gd <- edge_density(red_H) #Densidad de la red= a proporción de enlaces presentes de todos los enlaces posibles en la red
diameter(red_H, directed=TRUE) # Diametro de la red
dis<-mean_distance(red_H, directed = F) # Longitud promedio de camino más corto
dis
Camino_corto<-get_diameter(red_H) #Nodos de la longitud promedio de camino más corto. Figura 3.6 libro
Camino_corto
       #Codigo que muestra solo los nodos que creran la lngitud de camino más corto en la red
        longitud_camino<-length(Camino_corto)
        Posicionx_camino_corto<-c()
        Posiciony_camino_corto<-c()
        Posicionz_camino_corto<-c()
        for (i in 1:longitud_camino) {
        s<-Camino_corto[i]
        posix<-Posicionx[s]
        posiy<-Posiciony[s]
        posiz<-Posicionz[s]
        Posicionx_camino_corto[i]<-posix
        Posiciony_camino_corto[i]<-posiy
        Posicionz_camino_corto[i]<-posiz
        }
        Posicionx_camino_corto

        #Graficia los nodos del camino corto
        scatterplot3d(Posicionx_camino_corto, Posiciony_camino_corto, Posicionz_camino_corto,xlab = "X (mm)",ylab = "Y (mm)",zlab = "Z (mm)", pch = 20,tcl=0.2, color = "Pink",las=1)

centralidad<-centr_clo(red_H, mode="all", normalized=T)# Centralidad de cercania
nodo_maxima_centralidad<-which.max(closeness)
nodo_maxima_centralidad

bett<-betweenness(red_H, directed=T)
nodo_maxima_centralidad_intermedia<-which.max(bett) #Centralidad de intermediación
nodo_maxima_centralidad_intermedia

Wc <- cluster_walktrap(red_H) 
modularity(Wc) #Modularidad en la red

#________________________________________________Análisis tasa de entropía________________________________________
MATRIZADY<-MAdyaW # Matriz pesada
k<-degree(red_H) #Grado de cada nodo
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
  write.table(ALPHAvsH, "h(alpha)_red_Hellwin.txt")  
}
ALPHAvsH

plot(ValoresAlpha,H)


