#!/usr/bin/env Rscript

# Script de test pour la fonction parse_fasta

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
    if (substr(line, 1, 1) == ">") {
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

# Test avec le fichier example.fasta
cat("Test de la fonction parse_fasta\n")
cat("================================\n\n")

sequences <- parse_fasta("example.fasta")

cat("Nombre de séquences:", length(sequences), "\n\n")

for (i in seq_along(sequences)) {
  name <- names(sequences)[i]
  sequence <- sequences[[i]]
  cat(sprintf("Séquence %d: %s\n", i, name))
  cat(sprintf("  Longueur: %d nucléotides\n", nchar(sequence)))
  cat(sprintf("  Aperçu: %s...\n\n", substr(sequence, 1, 50)))
}

cat("Structure de la liste R:\n")
str(sequences)
