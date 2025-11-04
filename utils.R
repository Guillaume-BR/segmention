# Fonctions utilitaires pour le parsing de fichiers FASTA

# Fonction pour parser un fichier FASTA et le transformer en liste R
parse_fasta <- function(filepath) {
  # Lire le contenu du fichier
  lines <- readLines(filepath, warn = FALSE)
  
  # Initialiser les structures de données
  sequences <- list()
  current_header <- NULL
  current_sequence <- ""
  
  for (line in lines) {
    # Supprimer les espaces en début et fin de ligne
    line <- trimws(line)
    
    # Ignorer les lignes vides
    if (nchar(line) == 0) {
      next
    }
    
    # Vérifier si la ligne est un header (commence par >)
    if (startsWith(line, ">")) {
      # Si on avait déjà une séquence en cours, la sauvegarder
      if (!is.null(current_header)) {
        sequences[[current_header]] <- current_sequence
      }
      
      # Commencer une nouvelle séquence
      current_header <- substr(line, 2, nchar(line))  # Enlever le >
      current_sequence <- ""
    } else {
      # Ajouter la ligne à la séquence courante
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  # Sauvegarder la dernière séquence
  if (!is.null(current_header)) {
    sequences[[current_header]] <- current_sequence
  }
  
  return(sequences)
}

LetterToNumber <- function(letter){
  if (letter=="A") return(1)
  if (letter=="C") return(2)
  if (letter=="G") return(3)
  if (letter=="T") return(4)
}

mat_trans <- function(liste){
  P=table(c(NA,liste),c(liste,NA))
  return(P/rowSums(P))
}

Loi_Stationnaire<-function(P){
  E <- eigen(t(P))
  return(Re(E$vectors[,1])/sum(Re(E$vectors[,1])))
}

ForBack<-function(ADN,pi,f0,f1) 
  # ADN la séquence d'ADN observée, 
  # pi la matrice de transition, 
  # f0 et f1 les deux lois d'emission
{
  NumADN<-as.numeric(sapply(ADN,LetterToNumber))
  l=length(ADN)
  
  #matrices pour stocker
  P=matrix(ncol = 2,nrow = l) #prédiction
  Fil=matrix(ncol = 2,nrow =l ) #filtrage
  L=matrix(ncol = 2,nrow =l ) #lissage
  
  
  #i=1
  i=1
  #initialisation
  P[i,]=Loi_Stationnaire(pi)
  Fil[i,]= P[i,]*c(f0[NumADN[i]],f1[NumADN[i]])/sum(P[i,]*c(f0[NumADN[i]],f1[NumADN[i]]))
  
  #i=2...l
  for(i in 2:l)
  {
    #prediction
    P[i,]=Fil[i-1,]%*%pi
    
    #filtration
    Fil[i,]=P[i,]*c(f0[NumADN[i]],f1[NumADN[i]])/sum(P[i,]*c(f0[NumADN[i]],f1[NumADN[i]]))
  }
  
  #lissage
  i=l
  L[i,]=Fil[i,]
  for(i in l:2)
  { 
    L[i-1,]=Fil[i-1,]*((L[i,]/P[i,])%*%t(pi))
  }
  return(list(L=L,P=P,Fil=Fil))
}

NormTheta<-function(pi,f1,f2){
  d1 <- sqrt(sum(diag(t(pi)%*%pi)))#norme matricielle sur pi (remarque : calcul de la trace par sum(diag())
  d2 <- 1/2*(sum(abs(f1))+sum(abs(f2)))#somme des normes sur les lois f1 et f2
  return(d1+d2)
}

Baum_Welch<-function(ADN, pi0,f10,f20,epsilon)
{
  n=length(ADN)
  delta=2*epsilon
  NumADN<-as.numeric(sapply(ADN,LetterToNumber))
  count=0
  while (delta>epsilon)
  { 
    pi1=matrix(0,ncol=2,nrow=2)
    f11=vector(length=length(f10))
    f21=vector(length=length(f20))
    
    FB <- ForBack(ADN,pi0,f10,f20)
    
    #pi1
    for (q in 1:2){
      for (r in 1:2){
        pi1[q,r] <- sum(pi0[q,r]*FB$Fil[1:n-1,q]*FB$L[2:n,r]/FB$P[2:n,r])
      }
    }
    #normalisation
    pi1 <- pi1/rowSums(pi1)
    
    #f11 et f21
    for (x in 1:4){
      f11[x] <- sum(FB$L[1:n,1]*(NumADN==x))
      f21[x] <- sum(FB$L[1:n,2]*(NumADN==x))
    }
    #normalisation
    f11 <- f11/sum(f11)
    f21 <- f21/sum(f21)
    
    delta=NormTheta(pi0-pi1,f10-f11,f20-f21)/NormTheta(pi0,f10,f20)
    pi0<-pi1
    f10<-f11
    f20<-f21
    count=count+1
    
  }
  Theta<-list(pi=pi0,f1=f10,f2=f20,Numit=count)
  return(Theta)
}


LogVraisemblance<-function(ADN,pi,f1,f2)
{
  n=length(ADN)
  NumADN<-as.numeric(sapply(ADN,LetterToNumber))
  mu <- Loi_Stationnaire(pi)
  
  #initialisation
  alpha <- c(mu[1]*f1[NumADN[1]],mu[2]*f2[NumADN[1]])
  f <- diag(c(f1[NumADN[1]],f2[NumADN[1]]))
  for (i in 2:n){
    alpha <- alpha%*%pi%*%f
    f <- diag(c(f1[NumADN[i]],f2[NumADN[i]]))
  }
  return(log(sum(alpha)))
}

