library(shiny)

# Charger les fonctions utilitaires
source("utils.R")

server <- function(input, output, session) {
  
  # =========================================================================
  # VARIABLES RÉACTIVES
  # =========================================================================
  sequences_list <- reactiveVal(NULL)
  selected_sequence_name <- reactiveVal(NULL)
  hmm_results <- reactiveVal(NULL)
  analysis_data <- reactiveVal(NULL)
  
  # =========================================================================
  # 1️⃣ CHARGEMENT DU FICHIER FASTA
  # =========================================================================
  observeEvent(input$fasta_file, {
    req(input$fasta_file)
    
    file_ext <- tolower(tools::file_ext(input$fasta_file$name))
    valid_extensions <- c("fasta", "fa", "txt", "fna")
    
    if (!file_ext %in% valid_extensions) {
      showNotification(
        paste("❌ Extension invalide:", file_ext),
        type = "error",
        duration = 5
      )
      return()
    }
    
    tryCatch({
      sequences <- parse_fasta(input$fasta_file$datapath)
      
      if (length(sequences) == 0) {
        showNotification(
          "⚠️ Aucune séquence FASTA trouvée.",
          type = "warning",
          duration = 5
        )
        return()
      }
      
      nms <- names(sequences)
      if (is.null(nms) || all(nms == "")) {
        nms <- paste0("sequence_", seq_along(sequences))
      } else {
        nms <- ifelse(is.na(nms) | nms == "", 
                      paste0("sequence_", seq_along(sequences)), 
                      nms)
      }
      names(sequences) <- nms
      
      sequences_list(sequences)
      selected_sequence_name(NULL)
      hmm_results(NULL)
      analysis_data(NULL)
      
      showNotification(
        paste("✅ Chargé:", length(sequences), "séquence(s)"),
        type = "message",
        duration = 3
      )
    }, error = function(e) {
      showNotification(
        paste("❌ Erreur:", e$message),
        type = "error",
        duration = 5
      )
    })
  })
  
  # =========================================================================
  # 2️⃣ INDICATEUR DE FICHIER CHARGÉ
  # =========================================================================
  output$fileUploaded <- reactive({
    !is.null(sequences_list())
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  # =========================================================================
  # 3️⃣ SÉLECTEUR DYNAMIQUE
  # =========================================================================
  output$sequence_selector <- renderUI({
    req(sequences_list())
    seqs <- sequences_list()
    choices <- names(seqs)
    
    selectInput(
      "sequence_choice", 
      "Choisir la séquence:", 
      choices = choices, 
      selected = choices[1]
    )
  })
  
  # =========================================================================
  # 4️⃣ BOUTON "ANALYSER"
  # =========================================================================
  observeEvent(input$analyze_btn, {
    req(sequences_list())
    req(input$sequence_choice)
    
    seq_name <- input$sequence_choice
    seq_selected <- sequences_list()[[seq_name]]
    
    if (is.null(seq_selected) || seq_selected == "") {
      showNotification("❌ Séquence non trouvée.", type = "error")
      return()
    }
    
    selected_sequence_name(seq_name)
    
    showNotification(
      paste("🔍 Analyse en cours..."),
      type = "message",
      duration = 1
    )
    
    # =====================================================================
    # AFFICHAGE TEXTE
    # =====================================================================
    output$selected_analysis <- renderPrint({
      cat("╔════════════════════════════════════════════════════════════╗\n")
      cat("║          ANALYSE DE LA SÉQUENCE SÉLECTIONNÉE              ║\n")
      cat("╚════════════════════════════════════════════════════════════╝\n\n")
      
      cat("📌 Nom:", seq_name, "\n")
      cat("📏 Longueur:", nchar(seq_selected), "nucléotides\n\n")
      
      max_chars <- 1000
      cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
      cat("SÉQUENCE:\n")
      cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
      
      if (nchar(seq_selected) > max_chars) {
        cat(substr(seq_selected, 1, max_chars), "\n")
        cat("... [+", nchar(seq_selected) - max_chars, "caractères]\n\n")
      } else {
        cat(seq_selected, "\n\n")
      }
      
      chars <- strsplit(toupper(seq_selected), "")[[1]]
      nucs <- table(chars)
      
      cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
      cat("COMPOSITION:\n")
      cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
      
      for (nuc in c("A", "C", "G", "T")) {
        count <- if (nuc %in% names(nucs)) nucs[nuc] else 0
        pct <- (count / length(chars)) * 100
        cat(sprintf("%s: %4d (%.2f%%)\n", nuc, count, pct))
      }
      cat("\n")
    })
    
    # =====================================================================
    # CALCUL HMM
    # =====================================================================
    tryCatch({
      print("DEBUG: Début calcul HMM")
      
      pi <- matrix(c(0.65, 0.55, 0.35, 0.45), nrow = 2, ncol = 2)
      f0 <- matrix(c(0.4, 0.1, 0.1, 0.4), nrow = 4)
      f1 <- matrix(c(0.25, 0.25, 0.25, 0.25), nrow = 4)
      epsilon <- 0.001
      
      BW <- Baum_Welch(seq_selected, pi, f0, f1, epsilon)
      print(paste("DEBUG: Baum-Welch done, iterations:", BW$Numit))
      
      best_model <- ForBack(seq_selected, BW$pi, BW$f1, BW$f2)
      print(paste("DEBUG: ForBack done, dim P:", nrow(best_model$P), "x", ncol(best_model$P)))
      
      # Stocker dans une liste réactive
      analysis_data(list(
        model = best_model,
        seq_name = seq_name,
        iterations = BW$Numit
      ))
      
      print("DEBUG: analysis_data mise à jour")
      
      showNotification(
        paste("✅ Analyse réussie! (", BW$Numit, "itérations )"),
        type = "message",
        duration = 3
      )
      
    }, error = function(e) {
      print(paste("DEBUG: Erreur HMM -", e$message))
      showNotification(
        paste("❌ Erreur:", e$message),
        type = "error",
        duration = 5
      )
    })
  })
  
  # =====================================================================
  # ✅ GRAPHIQUE HMM - AVEC observe() POUR FORCER LE RENDU
  # =====================================================================
  observe({
    # Cette dépendance réactive force le rendu quand analysis_data() change
    data <- analysis_data()
    
    print(paste("DEBUG: observe() appelé, data is.null =", is.null(data)))
    
    output$sequence_plot <- renderPlot({
      print("DEBUG: renderPlot appelé DEPUIS observe()")
      
      # Récupérer les données
      data <- analysis_data()
      
      if (is.null(data)) {
        print("DEBUG: data is NULL")
        plot(1, 1, type = "n", axes = FALSE)
        text(1, 1, "Veuillez analyser une séquence d'abord", col = "gray")
        return()
      }
      
      model <- data$model
      seq_name <- data$seq_name
      
      print(paste("DEBUG: model non-null, dimensions:", nrow(model$P), "x", ncol(model$P)))
      
      if (is.null(model) || is.null(model$P)) {
        print("DEBUG: model$P is NULL")
        plot(1, 1, type = "n", axes = FALSE)
        text(1, 1, "Erreur: modèle invalide", col = "red")
        return()
      }
      
      n <- nrow(model$P)
      if (n == 0) {
        print("DEBUG: n = 0")
        plot(1, 1, type = "n", axes = FALSE)
        text(1, 1, "Erreur: séquence trop courte", col = "red")
        return()
      }
      
      print(paste("DEBUG: n =", n, ", création du graphique"))
      
      positions <- 1:n
      pred <- as.numeric(model$P[, 1])
      filt <- as.numeric(model$Fil[, 1])
      liss <- as.numeric(model$L[, 1])
      
      print(paste("DEBUG: pred length =", length(pred)))
      print(paste("DEBUG: pred first 5:", paste(head(pred, 5), collapse=", ")))
      
      y_max <- max(c(pred, filt, liss), na.rm = TRUE)
      if (is.na(y_max) || is.infinite(y_max) || y_max == 0) {
        y_max <- 1
      }
      
      print(paste("DEBUG: y_max =", y_max))
      print("DEBUG: plot() appelé")
      
      # ✅ CRÉER LE GRAPHIQUE
      plot(positions, pred, 
           type = "l", 
           col = "red", 
           lwd = 2.5,
           ylab = "Probabilité", 
           xlab = "Position",
           main = paste("Modèle de Markov Caché (HMM) -", seq_name),
           ylim = c(0, y_max * 1.2))
      
      lines(positions, filt, col = "blue", lwd = 2.5)
      lines(positions, liss, col = "orange", lwd = 2.5)
      
      grid(col = "gray80", lty = "dotted")
      
      legend("topright", 
             legend = c("Prédiction", "Filtrage", "Lissage"),
             col = c("red", "blue", "orange"), 
             lty = 1,
             lwd = 2.5,
             bty = "n")
      
      print("DEBUG: graphique terminé ✅")
    })
  })
  
  # =========================================================================
  # AUTRES AFFICHAGES
  # =========================================================================
  output$sequence_info <- renderPrint({
    req(sequences_list())
    sequences <- sequences_list()
    
    cat("SÉQUENCES CHARGÉES\n")
    cat("═════════════════════════════════════════\n\n")
    cat("Nombre:", length(sequences), "\n\n")
    
    for (i in seq_along(sequences)) {
      name <- names(sequences)[i]
      seq_length <- nchar(sequences[[i]])
      cat(sprintf("%d. %s (%d nt)\n", i, name, seq_length))
    }
  })
  
  output$list_structure <- renderPrint({
    req(sequences_list())
    str(sequences_list(), max.level = 2)
  })
  
  output$sequence_preview <- renderPrint({
    req(sequences_list())
    sequences <- sequences_list()
    
    cat("APERÇU DES SÉQUENCES\n")
    cat("═════════════════════════════════════════\n\n")
    
    for (i in seq_along(sequences)) {
      name <- names(sequences)[i]
      sequence <- sequences[[i]]
      
      cat(sprintf("%d. %s\n", i, name))
      cat(sprintf("   Longueur: %d nt\n", nchar(sequence)))
      cat("   Aperçu: ")
      
      if (nchar(sequence) > 100) {
        cat(substr(sequence, 1, 100), "...\n\n")
      } else {
        cat(sequence, "\n\n")
      }
    }
  })
}