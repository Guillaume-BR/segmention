# Application Shiny pour l'analyse de séquences ADN
# Permet de charger une séquence ADN au format FASTA et la transformer en liste R

library(shiny)

# Charger les fonctions utilitaires
source("utils.R")

# Lancer l'application
shinyApp(ui = ui, server = server)