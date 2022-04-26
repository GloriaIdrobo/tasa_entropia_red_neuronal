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

red <- graph_from_adjacency_matrix(MatrizAdyacencia)
plot(red, edge.arrow.size=.000001, edge.curved=0.000001, vertex.size=6,
     vertex.color="red" , vertex.frame.color="black", vertex.label.color="black",
     vertex.label.cex=0.00001,edge.color="Gray")


#________________________________________Análisis Topológico __________________________________________

# Cantidad de enlaces que tiene la red
gsize(red)
# Diametro de la red
diameter(red, directed=F)
# Camino de diametro de la red
get_diameter(red)

# Grado de cada nodo
k<-degree(red, mode=c("all"))
#Tabla que muestra nodos tiene grado k
table(k)
# Histograma de grado
hist(k, breaks = 30)
# Muestra el nodo que tiene mayor grado
which.max(k)
#Grafica de distribución de grado
deg.dist <- degree_distribution(red, cumulative=T, mode="all")
plot( x=0:max(k), y=1-deg.dist, pch=19, cex=0.5, col="red", 
      xlab="k", ylab="P(k)")

write.table(deg, "Gradonodos.csv") 
write.table(deg.dist, "Probabilidad_Grado.csv") 

# Centralidad intermedia de cada uno de los nodos
betweenness_red<-betweenness(red, directed = FALSE, normalized = T)
# Muestra el histograma de betweenness
hist(betweenness_red, breaks = 80)
which.max(betweenness_red)

# Centralidad de cercania de cada uno de los nodos
closeness<-closeness(red_H, mode="all")
which.max(closeness)
hist(closeness)

# Detecta las comunidades que tiene la red
kc = fastgreedy.community(red)
# Determina el tamaño de cada comunidad
sizes(kc)
# muestra los nodos que pertenecen a cada comunidad
membership(kc)
# Perform edge-betweenness community detection on network graph
gc = edge.betweenness.community(red)
# Determine sizes of each community
sizes(gc)

fc <- cluster_fast_greedy(red)
membership(fc)
sizes(fc)
# Densidad de la red
gd <- edge_density(red)

# Create numerical vector of vertex eigenvector centralities 
ec <- as.numeric(eigen_centrality(red) )

#__________________________CALCULO DE LA TASA DE ENTROPIA h __________________________________________________

Num_nodos<-length(V(red)) #Cantidad de nodos de la red
Num_nodos

MATRIZADY<-MAdyaW
k<-degree(red)

#Se crea la matriz de transición de probabilidad PI vacia
PI <- matrix(data=0, nrow=Num_nodos, ncol=Num_nodos)

#Se crea un vector que almacenará los valores de la tasa de entropia para cada alpha
H<-seq(1,61) 

#Valores del sesgo en los que se va a evaluar la tasa de entropía
ValoresAlpha<- seq(1,61) #se va a evaluar de alpha -3 a alpha 3 con intervalos de 0.1 es decir 60 valores de alpha
Alphapasos<--3.1
for (i in 1:61) {
  Alphapasos<-Alphapasos+0.1
  ValoresAlpha[i]<-Alphapasos
}
ValoresAlpha[] 

#Calculo de la ecuación de la tasa de entropía para 60 valores de alpha
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
      # cuando la matriz tiene PI_{i,j}=0 el valor de h se hace cero. ver ecuación. si no calcula el logaritmo natural de PI
      if (PI[i,j] <= 0) PIlnPIW <-0 else logaritmoN <- ln(PI[i,j])
      
      PIlnPIW<--PI[i,j]*Wstar[i]*logaritmoN
      mul<-mul+PIlnPIW
    }
    h[i]<-mul
  }
  H[alph]<-sum(h) # Valor total de la tasa de entropía para alpha
  #_____________________________________________________________________________________________
  #Creación de un txt para guardar los valores de la tasa de entropía en funcion de Alpha
  
  H(Alpha) <- data.frame(ValoresAlpha,H )
  write.table(ALPHAvsH, "h(alpha)_red_Bresking.txt")  
}
H(Alpha) # Muestra los resultados de h en funcio de alpha

