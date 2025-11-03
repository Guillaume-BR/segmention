#!/usr/bin/env Rscript

# Script de test pour la fonction parse_fasta

# Charger les fonctions utilitaires
source("utils.R")

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
