library(shiny)

ui <- fluidPage(
  titlePanel("Application de Segmentation de Séquences ADN"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Charger une séquence FASTA"),
      fileInput("fasta_file", 
                "Choisir un fichier FASTA",
                accept = c(".fasta", ".fa", ".txt", ".fna")),
      hr(),
      helpText("Chargez un fichier au format FASTA pour analyser la séquence ADN."),
      helpText("Le fichier sera transformé en liste R pour l'analyse."),
      h4("Sélection de la séquence à analyser:"),
      
      # UI dynamique pour sélectionner la séquence
      uiOutput("sequence_selector"),
      
      actionButton("analyze_btn", "Analyser la séquence sélectionnée", 
                   class = "btn-primary", width = "100%"),
      hr()
    ),
    
    mainPanel(
      # Message si aucun fichier n'est chargé
      conditionalPanel(
        condition = "output.fileUploaded == false",
        p("Veuillez charger un fichier FASTA pour commencer.")
      ),
      
      # Affichage des onglets uniquement si fichier chargé
      conditionalPanel(
        condition = "output.fileUploaded == true",
        tabsetPanel(
          id = "main_tabs",
          tabPanel(
            title = "Résultats globaux",
            h4("Informations sur les séquences chargées:"),
            verbatimTextOutput("sequence_info"),
            hr(),
            h4("Structure de la liste R:"),
            verbatimTextOutput("list_structure"),
            hr(),
            h4("Aperçu des séquences:"),
            verbatimTextOutput("sequence_preview")
          ),
          
          tabPanel(
            title = "Séquence sélectionnée",
            h4("Analyse de la séquence sélectionnée:"),
            verbatimTextOutput("selected_analysis"),
            hr()
          ),
          
          tabPanel(
            title = "Graphique HMM",
            h4("Modèle de Markov Caché - Analyse de la séquence:"),
            plotOutput("sequence_plot", width = "100%", height = "600px"),
            hr(),
            p("Le graphique affiche :"),
            p("🔴 Courbe rouge = Prédiction"),
            p("🔵 Courbe bleue = Filtrage"),
            p("🟠 Courbe orange = Lissage")
          )
        )
      )
    )
  )
)