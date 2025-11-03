# Application Shiny pour l'analyse de séquences ADN
# Permet de charger une séquence ADN au format FASTA et la transformer en liste R

library(shiny)

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
    
    # Parser le fichier FASTA
    tryCatch({
      sequences <- parse_fasta(input$fasta_file$datapath)
      sequences_list(sequences)
    }, error = function(e) {
      showNotification(
        paste("Erreur lors de la lecture du fichier:", e$message),
        type = "error"
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
