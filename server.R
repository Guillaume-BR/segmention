library(shiny)

# Charger les fonctions utilitaires
source("utils.R")

server <- function(input, output, session) {
  
  # Variable réactive pour stocker la liste des séquences
  sequences_list <- reactiveVal(NULL)
  # Stocke le nom (clé) de la séquence actuellement sélectionnée/analysee
  selected_sequence_name <- reactiveVal(NULL)
  
  # Observer le chargement du fichier
  observeEvent(input$fasta_file, {
    req(input$fasta_file)
    
    # Valider l'extension du fichier
    file_ext <- tools::file_ext(input$fasta_file$name)
    valid_extensions <- c("fasta", "fa", "txt", "FASTA", "FA", "TXT", "fna")
    
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
      
      # S'assurer que chaque séquence a un nom (utile si parse_fasta ne les fournit pas)
      nms <- names(sequences)
      if (is.null(nms)) nms <- rep("", length(sequences))
      nms <- ifelse(is.na(nms) | nms == "", paste0("sequence_", seq_along(sequences)), nms)
      names(sequences) <- nms
      
      sequences_list(sequences)
      # réinitialiser la sélection précédente
      selected_sequence_name(NULL)
      
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
  
  # UI dynamique: menu de sélection des séquences une fois le fichier chargé
  output$sequence_selector <- renderUI({
    req(sequences_list())
    seqs <- sequences_list()
    choices <- names(seqs)
    selectInput("sequence_choice", "Choisir la séquence :", choices = choices, selected = choices[1])
  })
  
  # Quand l'utilisateur clique sur le bouton "Analyser", on fixe la séquence sélectionnée
  observeEvent(input$analyze_btn, {
    req(sequences_list())
    chosen <- input$sequence_choice
    if (is.null(chosen) || chosen == "") {
      showNotification("Veuillez choisir une séquence avant d'analyser.", type = "warning")
      return()
    }
    selected_sequence_name(chosen)
    showNotification(paste("Séquence sélectionnée pour analyse :", chosen), type = "message", duration = 2)
  })
  
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
  
  # Afficher les détails de la séquence sélectionnée (après clic sur le bouton)
  output$selected_analysis <- renderPrint({
    req(selected_sequence_name(), sequences_list())
    seqs <- sequences_list()
    sel_name <- selected_sequence_name()
    seq <- seqs[[sel_name]]
    if (is.null(seq)) {
      cat("Séquence non trouvée.\n")
      return()
    }
    
    cat("Nom de la séquence sélectionnée:", sel_name, "\n")
    cat("Longueur :", nchar(seq), "nucléotides\n\n")
    
    # Affichage d'un aperçu plus long (jusqu'à 1000 caractères)
    max_chars <- 1000
    if (nchar(seq) > max_chars) {
      cat("Aperçu (premiers", max_chars, "caractères):\n")
      cat(substr(seq, 1, max_chars), "...\n")
    } else {
      cat("Séquence complète :\n")
      cat(seq, "\n")
    }
    
    # Ici, vous pouvez ajouter des analyses supplémentaires (composition, GC%, segmentation, ...)
    # Exemple simple: composition des nucléotides
    chars <- strsplit(toupper(seq), "")[[1]]
    nucs <- table(chars)
    cat("\nComposition des nucléotides:\n")
    print(nucs)
  })
}