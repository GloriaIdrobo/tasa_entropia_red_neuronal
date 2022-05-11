#_______________________________ Percolación y Calculo de la tasa de entropía sin sesgo (alpha = 0) para la Red percolada de Hellwig_______________________________

#install.packages("igraph") 
#install.packages("repr")
#install.packages("SciViews")
library("igraph")
library("repr")
library("SciViews")

#Importamos la matriz de adyacencia de la red de 2000 nodos creada previamente
MatrizA <- read.table("Matriz_A_Hellwig.txt")
MATRIZADY <- as.matrix(MatrizA)       # Convertimos en una matriz el txt importado
View(MatrizA) # Permite ver la matriz de adyacencia que hemos importado

#Creación de la matriz de adyacencia de 1 y ceros
MAya_NO_W<-matrix(data=0, nrow=Nodos, ncol=Nodos)

for (i in 1:Nodos) {
  for (j in 1:Nodos) {
    if (ProConexion[i,j]==0) MAya_NO_W[i,j]=0 else  MAya_NO_W[i,j]=1
  }
}
MAya_NO_W[]

#Con la matriz de adyacencia se crea la red
red_NO_W_H<- graph_from_adjacency_matrix(MAya_NO_W)
Num_nodos<-length(V(red_NO_W_H)) #Cantidad de nodos de la red
Num_nodos

#____________________________Proceso percolativo

CNE<-100 #Cantidad de Nodos a Eliminar
#Neurona aleatoria que será eliminada
Num_nodos<-length(V(red_NO_W_H)) #Cantidad de nodos de la red
EliminaNodos<- seq(1,CNE) # Vector que guardará la información de nodos que se eliminaran
Elinodos<-sample(1:Num_nodos,CNE,replace=F) #Genera numeros aleatorios entre 1 y el número total de nodos
EliminaNodos<-Elinodos # Nodos que se eliminan de la red
red_NO_W_H<-delete.vertices(graph = red_NO_W_H,EliminaNodos) #Elimina los nodos selecionados aleatoriamente
Num_nodos<-length(V(red_NO_W_H)) #Cantidad de nodos de la red
Num_nodos
MAdyaW<-MAdyaW[,-EliminaNodos] #Elimina las filas y la columnas de los nodos a eliminar en la matriz pesada
MAdyaW<-MAdyaW[-EliminaNodos,]
MAdyaW
dim(MAdyaW)
write.table(MAdyaW, "MAdyaW.txt")  

MATRIZADY<-MAdyaW
k<-degree(red_NO_W_H) # Calcula el grado de cada uno de los nodos 
red_NO_W_H[]

Num_nodos<-length(V(red_NO_W_H)) #Cantidad de nodos de la red
Num_nodos

#Se crea la matriz de transición de probabilidad PI vacia
PI <- matrix(data=0, nrow=Num_nodos, ncol=Num_nodos)

#Se crea un vector que almacenará los valores de la tasa de entropia para cada alpha
#H<-seq(1,35) #Solo en la primera vuelta se corre que es para crear el vector, luego no se corre
#Ntotal<-seq(1,35)

alpha <- 0  #En cada paso alpha cambia de valor para calcular h
Ci <-seq(1,Num_nodos) # Se crea el vector Ci que guardara el valor de ci para cada i 

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
    gr<-k[i]      # Valor del grado que tiene el nodo j, es decir el nodo destino
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
#Calculo de la matriz de transición de probabilidad PI
for (i in 1:Num_nodos) {
  for (j in 1:Num_nodos) {
    pruebaa<-is.na(PI[i,j])
    sino<-PI[i,j]
    if (pruebaa==TRUE) PI[i,j]=0 else  PI[i,j]=sino
  }
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

Ntotal[1]<-Num_nodos # Guarda la cantidad de nodos en la red en la posición [i]
H[1]<-sum(h) # Guarda el valor de la tasa de entropía en la red para el número de nodos
H_percolada<- data.frame(Ntotal,H ) 
H_percolada # 
write.table(H_percolada, "Resultados_finales_Hellwig_percolada.csv")  
