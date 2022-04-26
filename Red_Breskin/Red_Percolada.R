#_______________________________ Percolación y Calculo de h para la Red percolada de Breskin_______________________________

#install.packages("igraph") 
#install.packages("repr")
#install.packages("SciViews")
library("igraph")
library("repr")
library("SciViews")

#Importamos la matriz de adyacencia de la red de 2000 nodos creada previamente
MatrizA <- read.table("Matriz_A_Breskin.txt")
MATRIZADY <- as.matrix(MatrizA)       # Convertimos en una matriz el txt importado
View(MatrizA) # Permite ver la matriz de adyacencia que hemos importado

#Con la matriz de adyacencia se crea la red
red <- graph_from_adjacency_matrix(MATRIZADY)
#______________________________________________________________________________________________
#Proceso percolativo- Eliminación de nodos(neuronas)

CNE<-100 #Cantidad de Nodos a Eliminar
#Neurona aleatoria que será eliminada
Num_nodos<-length(V(red)) #Cantidad de nodos de la red
EliminaNodos<- seq(1,CNE) # Vector que guardará la información de nodos que se eliminaran
Elinodos<-sample(1:Num_nodos,CNE,replace=F) #Genera numeros aleatorios entre 1 y el número total de nodos
EliminaNodos<-Elinodos # Nodos que se eliminan de la red
red<-delete.vertices(graph = red,EliminaNodos) #Elimina los nodos selecionados aleatoriamente

#Propiedades de la red percolada
Num_nodos<-length(V(red)) #Cantidad de nodos de la red
Num_nodos
MATRIZADY<-red[] 
MATRIZADYP <- as.matrix(MATRIZADY)

#Definimos las propiedades para el calculo de la tasa de entropia
Num_nodos<-length(V(red)) #Cantidad de nodos de la red
Num_nodos
k <- degree(red) # Se guarda en la variable k el valor deL grado que tiene cada uno de los nodos

#Se crea la matriz de transición de probabilidad PI vacia
PI <- matrix(data=0, nrow=Num_nodos, ncol=Num_nodos)
#Se crea un vector que almacenará los valores de la tasa de entropia para cada alpha
#H<-seq(1,20) 
#Ntotal<-seq(1,20)
#___________________________________________________________________________________________

#Calculo de la ecuación de la tasa de entropía para 60 valores de alpha

alpha <-0 #En cada paso alpha cambia de valor para calcular h
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
#  Wstar # Mostrar la DISTRIBUCIÓN ESTACIONARIA W*

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
    if (PI[i,j] == 0) PIlnPIW <-0 else logaritmoN <- ln(PI[i,j])
    
    PIlnPIW<--PI[i,j]*Wstar[i]*logaritmoN
    mul<-mul+PIlnPIW
  }
  h[i]<-mul
}

Ntotal[20]<-Num_nodos
H[20]<-sum(h) # Valor total de la tasa de entropía para alpha
#_____________________________________________________________________________________________
#Creación de un txt para guardar los valores de la tasa de entropía en funcion de Alpha

H_percolada<- data.frame(Ntotal,H )
write.table(H_percolada, "h(alpha)_red_percolada_Breskin.txt")  
H_percolada

