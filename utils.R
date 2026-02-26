# Fonctions utilitaires pour le parsing de fichiers FASTA et algorithmes HMM

# =========================================================================
# PARSING FASTA
# =========================================================================
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

# =========================================================================
# CONVERSION LETTRE → NOMBRE
# =========================================================================
LetterToNumber <- function(letter) {
  letter <- toupper(letter)
  if (letter == "A") return(1)
  if (letter == "C") return(2)
  if (letter == "G") return(3)
  if (letter == "T") return(4)
  # Pour les caractères invalides, retourner 1 par défaut
  return(1)
}

# ✅ Convertir une séquence ADN (chaîne) en vecteur de nombres
ADNToNumeric <- function(adn_string) {
  # Découper la chaîne en caractères individuels
  chars <- strsplit(adn_string, "")[[1]]
  # Convertir chaque caractère en nombre
  numbers <- sapply(chars, LetterToNumber, USE.NAMES = FALSE)
  return(as.numeric(numbers))
}

# =========================================================================
# MATRICE DE TRANSITION
# =========================================================================
mat_trans <- function(liste) {
  P <- table(c(NA, liste), c(liste, NA))
  P <- P / rowSums(P)
  return(as.matrix(P))
}

# =========================================================================
# LOI STATIONNAIRE (distribution d'équilibre)
# =========================================================================
Loi_Stationnaire <- function(P) {
  P <- as.matrix(P)  # ✅ Conversion explicite en matrice
  E <- eigen(t(P))
  result <- Re(E$vectors[, 1])
  result <- result / sum(result)
  return(as.numeric(result))  # ✅ Conversion en vecteur numérique
}

# =========================================================================
# ALGORITHME FORWARD-BACKWARD (inférence)
# =========================================================================
ForBack <- function(ADN, pi, f0, f1) {
  # Convertir les entrées en types corrects
  ADN <- as.character(ADN)
  pi <- as.matrix(pi)
  f0 <- as.numeric(f0)
  f1 <- as.numeric(f1)
  
  # ✅ Conversion numérique correcte de la séquence ADN
  NumADN <- ADNToNumeric(ADN)
  l <- length(NumADN)
  
  if (l < 1) {
    stop("Séquence ADN vide ou invalide")
  }
  
  # Initialiser les matrices de résultats
  P <- matrix(0, ncol = 2, nrow = l)      # Prédiction
  Fil <- matrix(0, ncol = 2, nrow = l)    # Filtrage
  L <- matrix(0, ncol = 2, nrow = l)      # Lissage
  
  # ───────────────────────────────────────────────────────────────
  # FORWARD (Prédiction + Filtrage)
  # ───────────────────────────────────────────────────────────────
  
  # i = 1 (initialisation)
  i <- 1
  P[i, ] <- Loi_Stationnaire(pi)
  
  # Émission pour la première position
  emission <- c(f0[NumADN[i]], f1[NumADN[i]])
  emission <- as.numeric(emission)
  
  # Filtrage initial
  pred_emiss <- as.numeric(P[i, ] * emission)
  sum_pred_emiss <- sum(pred_emiss)
  
  if (sum_pred_emiss > 0) {
    Fil[i, ] <- pred_emiss / sum_pred_emiss
  } else {
    Fil[i, ] <- c(0.5, 0.5)  # Valeur par défaut
  }
  
  # i = 2 à l (boucle forward)
  if (l > 1) {
    for (i in 2:l) {
      # Prédiction: P[i] = Fil[i-1] * pi (produit matrice-vecteur)
      P[i, ] <- as.numeric(Fil[i - 1, ] %*% pi)
      
      # Émission pour position i
      emission <- c(f0[NumADN[i]], f1[NumADN[i]])
      emission <- as.numeric(emission)
      
      # Filtrage
      pred_emiss <- as.numeric(P[i, ] * emission)
      sum_pred_emiss <- sum(pred_emiss)
      
      if (sum_pred_emiss > 0) {
        Fil[i, ] <- pred_emiss / sum_pred_emiss
      } else {
        Fil[i, ] <- Fil[i - 1, ]  # Garder valeur précédente
      }
    }
  }
  
  # ───────────────────────────────────────────────────────────────
  # BACKWARD (Lissage)
  # ───────────────────────────────────────────────────────────────
  
  i <- l
  L[i, ] <- Fil[i, ]
  
  if (l > 1) {
    for (i in (l - 1):1) {
      # Division élément par élément
      division <- L[i + 1, ] / P[i + 1, ]
      
      # Éviter les divisions par zéro
      division <- ifelse(is.nan(division) | is.infinite(division), 0, division)
      
      # Vecteur de transition
      transition_factor <- as.numeric(division %*% t(pi))
      transition_factor <- ifelse(is.nan(transition_factor) | is.infinite(transition_factor), 1, transition_factor)
      
      L[i, ] <- as.numeric(Fil[i, ] * transition_factor)
      
      # Normaliser si nécessaire
      sum_L <- sum(L[i, ])
      if (sum_L > 0) {
        L[i, ] <- L[i, ] / sum_L
      }
    }
  }
  
  return(list(L = L, P = P, Fil = Fil))
}

# =========================================================================
# NORME DE THETA (pour convergence)
# =========================================================================
NormTheta <- function(pi, f1, f2) {
  pi <- as.matrix(pi)
  f1 <- as.numeric(f1)
  f2 <- as.numeric(f2)
  
  # Norme matricielle sur pi (trace)
  d1 <- sqrt(sum(diag(t(pi) %*% pi)))
  
  # Somme des normes sur f1 et f2
  d2 <- 0.5 * (sum(abs(f1)) + sum(abs(f2)))
  
  return(d1 + d2)
}

# =========================================================================
# ALGORITHME BAUM-WELCH (optimisation du modèle)
# =========================================================================
Baum_Welch <- function(ADN, pi0, f10, f20, epsilon) {
  # Conversion des entrées
  ADN <- as.character(ADN)
  pi0 <- as.matrix(pi0)
  f10 <- as.numeric(f10)
  f20 <- as.numeric(f20)
  epsilon <- as.numeric(epsilon)
  
  # ✅ Conversion numérique correcte
  NumADN <- ADNToNumeric(ADN)
  n <- length(NumADN)
  
  if (n < 2) {
    warning("Séquence trop courte pour Baum-Welch")
    return(list(pi = pi0, f1 = f10, f2 = f20, Numit = 0))
  }
  
  delta <- 2 * epsilon
  count <- 0
  
  # Itérations Baum-Welch
  while (delta > epsilon && count < 100) {  # Limite d'itérations
    pi1 <- matrix(0, ncol = 2, nrow = 2)
    f11 <- numeric(4)
    f21 <- numeric(4)
    
    # Forward-Backward
    FB <- ForBack(ADN, pi0, f10, f20)
    
    # ─────────────────────────────────────────────────────────────
    # Mise à jour de pi (matrice de transition)
    # ─────────────────────────────────────────────────────────────
    for (q in 1:2) {
      for (r in 1:2) {
        # Accès aux lignes correctes (1 à n-1 et 2 à n)
        indices_up <- 1:(n - 1)
        indices_down <- 2:n
        
        numerator <- sum(
          pi0[q, r] * FB$Fil[indices_up, q] * 
            (FB$L[indices_down, r] / FB$P[indices_down, r])
        )
        
        # Gestion des NaN
        if (is.na(numerator) || is.nan(numerator)) {
          numerator <- 0
        }
        
        pi1[q, r] <- as.numeric(numerator)
      }
    }
    
    # Normalisation de pi1
    row_sums <- rowSums(pi1)
    row_sums <- ifelse(row_sums == 0, 1, row_sums)
    pi1 <- pi1 / row_sums
    
    # ─────────────────────────────────────────────────────────────
    # Mise à jour de f1 et f2 (matrices d'émission)
    # ─────────────────────────────────────────────────────────────
    for (x in 1:4) {
      match_x <- as.numeric(NumADN == x)
      f11[x] <- sum(FB$L[1:n, 1] * match_x)
      f21[x] <- sum(FB$L[1:n, 2] * match_x)
    }
    
    # Normalisation de f11 et f21
    f11 <- as.numeric(f11)
    f21 <- as.numeric(f21)
    
    sum_f11 <- sum(f11)
    sum_f21 <- sum(f21)
    
    f11 <- if (sum_f11 > 0) f11 / sum_f11 else rep(0.25, 4)
    f21 <- if (sum_f21 > 0) f21 / sum_f21 else rep(0.25, 4)
    
    # ─────────────────────────────────────────────────────────────
    # Calcul de la convergence
    # ─────────────────────────────────────────────────────────────
    norm_new <- NormTheta(pi1, f11, f21)
    norm_old <- NormTheta(pi0, f10, f20)
    
    if (norm_old > 0) {
      delta <- abs(norm_new - norm_old) / norm_old
    } else {
      delta <- 0
    }
    
    # Mise à jour pour itération suivante
    pi0 <- pi1
    f10 <- f11
    f20 <- f21
    count <- count + 1
  }
  
  return(list(pi = pi0, f1 = f10, f2 = f20, Numit = count))
}

# =========================================================================
# VRAISEMBLANCE LOGARITHMIQUE
# =========================================================================
LogVraisemblance <- function(ADN, pi, f1, f2) {
  ADN <- as.character(ADN)
  pi <- as.matrix(pi)
  f1 <- as.numeric(f1)
  f2 <- as.numeric(f2)
  
  # ✅ Conversion numérique correcte
  NumADN <- ADNToNumeric(ADN)
  n <- length(NumADN)
  
  mu <- Loi_Stationnaire(pi)
  
  # Initialisation
  alpha <- c(mu[1] * f1[NumADN[1]], mu[2] * f2[NumADN[1]])
  
  # Forward
  if (n > 1) {
    for (i in 2:n) {
      f <- diag(c(f1[NumADN[i]], f2[NumADN[i]]))
      alpha <- as.numeric(alpha %*% pi %*% f)
    }
  }
  
  result <- log(sum(alpha))
  
  return(as.numeric(result))
}