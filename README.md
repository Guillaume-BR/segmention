# segmention
Application Shiny permettant de déterminer la probabilité qu'une région d'une séquence d'ADN soit codante ou non codante

## Description

Cette application Shiny permet de :
- Charger une séquence ADN au format FASTA
- Transformer la séquence en une liste R pour l'analyse
- Visualiser les informations sur les séquences chargées

## Installation

Pour utiliser cette application, vous devez avoir R et le package Shiny installés.

```r
# Installer Shiny si nécessaire
install.packages("shiny")
```

## Utilisation

Pour lancer l'application :

```r
library(shiny)
runApp()
```

Ou depuis la ligne de commande :

```bash
R -e "shiny::runApp()"
```

## Format FASTA

L'application accepte les fichiers au format FASTA standard :
- Les en-têtes commencent par `>`
- Les séquences peuvent être sur une ou plusieurs lignes
- Extensions acceptées : `.fasta`, `.fa`, `.txt`

### Exemple de fichier FASTA

```
>sequence1 Example DNA sequence 1
ATGCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGA
>sequence2 Example DNA sequence 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAG
```

Un fichier exemple (`example.fasta`) est fourni dans le dépôt pour tester l'application.

## Structure du projet

- `app.R` : Application Shiny principale
- `example.fasta` : Fichier FASTA d'exemple pour les tests
- `test_parse.R` : Script de test pour la fonction de parsing FASTA

## Fonctionnalités

### Chargement de fichier FASTA
L'application permet de charger un fichier FASTA via une interface utilisateur intuitive.

### Transformation en liste R
Le fichier FASTA est automatiquement parsé et transformé en une liste R où :
- Les noms des éléments de la liste correspondent aux en-têtes des séquences
- Les valeurs sont les séquences ADN sous forme de chaînes de caractères

### Affichage des résultats
L'application affiche :
- Le nombre de séquences chargées
- Les informations détaillées pour chaque séquence (nom, longueur)
- La structure de la liste R créée
- Un aperçu des séquences

## Code de segmentation

Une fois les séquences chargées et transformées en liste R, vous pouvez appliquer vos algorithmes de segmentation et d'analyse sur ces données. 
