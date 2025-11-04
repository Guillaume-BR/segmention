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
  observeEvent(input$analyze_btn, {
    
    # -------------------------------
    # 1️⃣ Récupérer la séquence sélectionnée
    # -------------------------------
    req(input$selected_sequence)
    seq_selected <- sequences()[[input$selected_sequence]]
    if (is.null(seq_selected)) {
      output$selected_analysis <- renderPrint({ cat("Séquence non trouvée.\n") })
      return()
    }
    
    # -------------------------------
    # 2️⃣ Analyse texte
    # -------------------------------
    output$selected_analysis <- renderPrint({
      cat("Nom de la séquence sélectionnée:", input$selected_sequence, "\n")
      cat("Longueur :", nchar(seq_selected), "nucléotides\n\n")
      
      # Aperçu jusqu'à 1000 caractères
      max_chars <- 1000
      if (nchar(seq_selected) > max_chars) {
        cat("Aperçu (premiers", max_chars, "caractères):\n")
        cat(substr(seq_selected, 1, max_chars), "...\n")
      } else {
        cat("Séquence complète :\n")
        cat(seq_selected, "\n")
      }
      
      # Composition des nucléotides
      chars <- strsplit(toupper(seq_selected), "")[[1]]
      nucs <- table(chars)
      cat("\nComposition des nucléotides:\n")
      print(nucs)
    })
    
    # -------------------------------
    # 3️⃣ Calcul HMM / Baum-Welch
    # -------------------------------
    pi <- matrix(c(0.65,0.55,0.35,0.45), nrow=2, ncol=2)
    f0 <- matrix(c(0.4,0.1,0.1,0.4), nrow=4)
    f1 <- matrix(c(0.25,0.25,0.25,0.25), nrow=4)
    epsilon <- 0.001
    
    BW <- Baum_Welch(seq_selected, pi, f0, f1, epsilon)
    best_model <- ForBack(seq_selected, BW$pi, BW$f1, BW$f2)
    
    # -------------------------------
    # 4️⃣ Graphique HMM
    # -------------------------------
    output$sequence_plot <- renderPlot({
      plot(best_model$P[,1], type="l", col="red", ylab="Probabilité", xlab="Position")
      lines(best_model$F[,1], col="blue")
      lines(best_model$L[,1], col="orange", lwd=2)
      legend("bottomright", legend=c("Prediction","Filtrage","Lissage"),
             col=c("red","blue","orange"), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n",
             lty=1)
    })
  })
}