# Application Shiny pour l'analyse de séquences ADN
# Permet de charger une séquence ADN au format FASTA et la transformer en liste R

library(shiny)

# Charger les fonctions utilitaires
source("utils.R")

# Interface utilisateur
ui <- fluidPage(
  titlePanel("Application de Segmentation de Séquences ADN"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Charger une séquence FASTA"),
      fileInput("fasta_file", 
                "Choisir un fichier FASTA",
                accept = c(".fasta", ".fa", ".txt")),
      hr(),
      helpText("Chargez un fichier au format FASTA pour analyser la séquence ADN."),
      helpText("Le fichier sera transformé en liste R pour l'analyse.")
    ),
    
    mainPanel(
      h3("Résultats"),
      conditionalPanel(
        condition = "output.fileUploaded",
        h4("Informations sur les séquences chargées:"),
        verbatimTextOutput("sequence_info"),
        hr(),
        h4("Structure de la liste R:"),
        verbatimTextOutput("list_structure"),
        hr(),
        h4("Aperçu des séquences:"),
        verbatimTextOutput("sequence_preview")
      ),
      conditionalPanel(
        condition = "!output.fileUploaded",
        p("Veuillez charger un fichier FASTA pour commencer.")
      )
    )
  )
)

# Serveur
server <- function(input, output, session) {
  
  # Variable réactive pour stocker la liste des séquences
  sequences_list <- reactiveVal(NULL)
  
  # Observer le chargement du fichier
  observeEvent(input$fasta_file, {
    req(input$fasta_file)
    
    # Valider l'extension du fichier
    file_ext <- tools::file_ext(input$fasta_file$name)
    valid_extensions <- c("fasta", "fa", "txt", "FASTA", "FA", "TXT")
    
    if (!file_ext %in% valid_extensions) {
      showNotification(
        paste("Extension de fichier non valide:", file_ext, 
              ". Extensions acceptées: .fasta, .fa, .txt"),
        type = "error",
        duration = 5
      )
      return()
    }
    
    # Parser le fichier FASTA
    tryCatch({
      sequences <- parse_fasta(input$fasta_file$datapath)
      
      # Vérifier que le fichier contient au moins une séquence
      if (length(sequences) == 0) {
        showNotification(
          "Le fichier ne contient aucune séquence FASTA valide.",
          type = "warning",
          duration = 5
        )
        return()
      }
      
      sequences_list(sequences)
      showNotification(
        paste("Chargement réussi:", length(sequences), "séquence(s) trouvée(s)"),
        type = "message",
        duration = 3
      )
    }, error = function(e) {
      showNotification(
        paste("Erreur lors de la lecture du fichier:", e$message),
        type = "error",
        duration = 5
      )
    })
  })
  
  # Vérifier si un fichier est chargé
  output$fileUploaded <- reactive({
    !is.null(sequences_list())
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  # Afficher les informations sur les séquences
  output$sequence_info <- renderPrint({
    req(sequences_list())
    sequences <- sequences_list()
    
    cat("Nombre de séquences:", length(sequences), "\n\n")
    
    for (i in seq_along(sequences)) {
      name <- names(sequences)[i]
      seq_length <- nchar(sequences[[i]])
      cat(sprintf("Séquence %d: %s\n", i, name))
      cat(sprintf("  Longueur: %d nucléotides\n\n", seq_length))
    }
  })
  
  # Afficher la structure de la liste R
  output$list_structure <- renderPrint({
    req(sequences_list())
    str(sequences_list())
  })
  
  # Afficher un aperçu des séquences
  output$sequence_preview <- renderPrint({
    req(sequences_list())
    sequences <- sequences_list()
    
    for (i in seq_along(sequences)) {
      name <- names(sequences)[i]
      sequence <- sequences[[i]]
      
      cat(sprintf("=== %s ===\n", name))
      
      # Afficher les 100 premiers caractères
      if (nchar(sequence) > 100) {
        cat(substr(sequence, 1, 100), "...\n\n")
      } else {
        cat(sequence, "\n\n")
      }
    }
  })
}

# Lancer l'application
shinyApp(ui = ui, server = server)
