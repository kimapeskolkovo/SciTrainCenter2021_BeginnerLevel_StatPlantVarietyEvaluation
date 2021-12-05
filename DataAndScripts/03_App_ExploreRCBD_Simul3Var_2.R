library(shiny) 

ui <- fluidPage(tabsetPanel(
  tabPanel("Main",
  titlePanel(h2("Explore the key elements of RCBD by simulating datasets using six varieties")),
  
  sidebarPanel(width = 2, h4("Parameters to explore"),br(),
   
   sliderInput(inputId = "EffSizeVariety", 
               label = "Effect size for varieties", 
               value = 0.6, min = 0, max = 6, step = 0.1),
   
   sliderInput(inputId = "EffSizeBlock", 
               label = "Effect size for blocks", 
               value = 1, min = 0, max = 6, step = 0.1),
   
   sliderInput(inputId = "Sigma2E", 
               label = "Residual variance", 
               value = 4, min = 0.5, max = 20, step = 0.5),
   
   numericInput(inputId = 'NumBlocks', 
                label = 'Number of blocks (2 to 15)', 
                value = 3, min = 2, max = 15),
   
   numericInput(inputId = 'Mu', label = "Grand mean", 
                value = 100),
   
   sliderInput(inputId = "Iter", 
               label = "num. of simulated datasets", 
               value = 30, min = 2, max = 200, step = 2),
   
   actionButton("simule", "Simulate"), 
   br(), 
   br(), 
  ),
  
  mainPanel(h3("Datasets and Analysis"),
            fluidRow(
              textOutput("NumVarieties"), 
              splitLayout(style = "border: 1px solid silver:", cellWidths = c("45%", "45%"), 
                          plotOutput("Graph2"), 
                          plotOutput("Graph3")
                          )
            ), 

            plotOutput("Graph1"),

            h3("Efficiency in detecting differences among varieties"),
            "P-values of analyses",
            verbatimTextOutput("StatsPvalues"),

    # tableOutput("Combin"),
    # verbatimTextOutput("Xmatrix"), 
    
    # verbatimTextOutput("BlockEffectSize"), 
    # verbatimTextOutput("VarEffectSize"),     
    # tableOutput("Betas"),
    
    # "Predicted mean values",
    # tableOutput("Predicted0"),
    h3("Expected mean values of the randomized design:"),
    h4("Yield_{Genotype in Block} = Mu + Genotype effect + Block effect"),
    fluidRow(
      column(4, tableOutput("Predicted")),
      column(4, verbatimTextOutput("VarEffects")), 
      column(4, verbatimTextOutput("BlockEffects"))
    ), 
    
    # verbatimTextOutput("BlockEffectsRecalc"),
    # textOutput("BlockEffectSizeRecalc"),
    # verbatimTextOutput("VarEffectsRecalc"),
    # textOutput("VarEffectSizeRecalc"),
    
    br(),
    h3("Simulated Values: Expected mean values PLUS Residual Variation (non controlled environmental effect)"),
    verbatimTextOutput("Donnees"),


    # "Simulated effect sizes for Varieties per samples",
    # verbatimTextOutput("VarietyEffectSizeObs"),
    # "Simulated effect sizes for Blocks per samples",
    # verbatimTextOutput("BlockEffectSizeObs"),

    h3("P-values of ANOVAs on simulated data"),
    tableOutput("Pvalues")   

  )
), # fin tabPanel Main

tabPanel("Means",
         h3("Grand mean per simulated dataset, as differences from Grand Mean"),
         verbatimTextOutput("MuSimulated"),
         
         h3("Simulated means per Variety"),
         verbatimTextOutput("MoyVarietyObs"),
         
         h3("Simulated means per Block"),
         verbatimTextOutput("MoyBlockObs")
         ),

tabPanel("Anova details",
         # h3(paste("ANOVA tables for 'Full model' of the simulated datasets")),
         # verbatimTextOutput("AnovasFull"),
         # 
         # h3(paste("ANOVA tables for 'Reduced model' of the simulated datasets")),
         # verbatimTextOutput("AnovasReduced"),
         
         fluidRow(
           column(6, verbatimTextOutput("AnovasFull")),
           column(6, verbatimTextOutput("AnovasReduced"))
         )
)

) # fin tabsetPanel()
) ## fin ui()


# Dans cette application, après avoir appuyé sur un bouton, un nombre aléatoire et des données pour une table et un tracé sont générés.
# Ensuite, le nombre, les données pour la table et un tracé sont renvoyés sous forme de liste et rendus avec les fonctions render* appropriées.

server <- function(input, output) { 
  
  Analyses <- eventReactive(input$simule, {   ## a priori, cela evite d'utiliser reactive()
    
    ##################   general definitions  ################## 
    
    # To erase all graphs
    graphics.off()
    # To erase objects from the working space - Clean up of the memory
    rm(list = ls())
    # use of the constraint 'set-to-zero' for ANOVAs 
    options(contrasts=c('contr.treatment','contr.poly'))
    
    ## Loading of the R packages needed for the analysis.
    library(ggplot2) # Needed for some graphs (e.g. bwplots)
    library(xtable)
    library(gridExtra)
    library(dplyr)
    library(desplot)
    library(RColorBrewer)
    
    # # for debugging
    # input <- list(NumBlocks = 3,
    #               EffSizeBlock = 0.6,
    #               EffSizeVariety = 1,
    #               Sigma2E = 4,
    #               Mu = 100)
    
    # Case with 6 varieties
    NumVars <- 6
    N <- input$NumBlocks * NumVars
    
    ## permute a set of NumVars varieties. 
    library(gtools)
    Lines1 <- c(paste("Variety", 1:NumVars, sep=""))
    # Lines2 <- c(paste("Variety", 7:12, sep=""))
    
    permut1 <- permutations(n = NumVars, r = NumVars, Lines1)
    # permut2 <- permutations(n=6, r=6, Lines2)
    
    choix1 <- sample(1:dim(permut1)[1], 1) #; choix2 <- sample(1:dim(permut2)[1], 1)
 #   Genotypes <- c(permut1[choix1, ], permut2[choix2, ])
    Genotypes <- c(permut1[choix1, ])
    
    ## the design
    Combinations <- expand.grid(genotype = Genotypes, 
                                           block = paste('Block', seq(1:input$NumBlocks), sep = ''))
    
    ## the design matrix  
    X <- model.matrix( ~ Combinations$genotype + Combinations$block)

    ## from Cohens' f computations, we get:
    SumBiCarre <- (N/NumVars) * input$EffSizeBlock^2 * input$Sigma2E 
    SumGiCarre <- (N/input$NumBlocks) * input$EffSizeVariety^2 * input$Sigma2E 
    #
    # on cree trois classes effects: moins, zero, plus
    # astuce: (titi %% 3 ) -1
    titi <- seq(1:input$NumBlocks)
    Blocki <- ((titi %% 3) - 1 ) * sqrt(SumBiCarre / sum(abs((titi %% 3) - 1)) )
    if (input$NumBlocks %in% c(5,8,11,14)){Blocki[input$NumBlocks] <- 0}
    # to check
    # sum(Blocki)  # sum-to-zero non compulsory
    # t(Blocki) %*% Blocki  # to check  with SumBiCarre
    
    toto <- seq(1:NumVars)
    Varietyi<- ((toto %% 3) - 1 )* sqrt(SumGiCarre/ (sum(abs((toto %% 3) - 1))) )
    if (NumVars %in% c(5,8,11,14)){Varietyi[NumVars] <- 0}
    # to check
    # sum(Varietyi)   # sum-to-zero non compulsory
    # t(Varietyi) %*% Varietyi  # to check  with SumGiCarre

    # Vector of effects for the two fixed factors, from Effects size.  Set_to_zero constraint for the first level
    #  Mu <- 100   # to vary, cf ui()
    beta <- c(input$Mu,
             Varietyi[-1],
             Blocki[-1])

    ## predicted values using  P = G + E ; E being the controlled environmental variation
    predits <- X %*% beta
    PrAdits <- as.data.frame(cbind( predits, Combinations ))
    colnames(PrAdits) <- c('predits','Variety','Block')

    ## we now need to permute the positions of the Varieties within each block
    truc <- c()
    for (k in 1:input$NumBlocks){
      permutRowsWithinBlocks <- permutations(n = NumVars, r = NumVars, seq(1:NumVars))
      truc <- c( truc, permutRowsWithinBlocks[ sample(1:dim(permutRowsWithinBlocks)[1], 1), ] + (k-1)*NumVars )  
    }
    # to check permutation trick
    # truc
    # print(cbind( truc, c(1:N), PrAdits))
    
    Tempo <- PrAdits[truc, 1:2]
    Predits <- cbind(Tempo, PrAdits[3])

    ## (re)compute effects sizes using Cohen's f formulas to confirm fixed effects values are OK
    ## these computations of effect sizes can be also verified @ https://webpower.psychstat.org/models/means04/effectsize.php
    (MoyBlock <- PrAdits %>% group_by(Block) %>% summarize(moy = mean(predits)) %>% select(moy) %>%  - input$Mu %>% pull)
    (MoyBlock <- Predits %>% group_by(Block) %>% summarize(moy = mean(predits)) %>% select(moy) %>%  - input$Mu %>% pull)
    (BlockEffectSizeRecalc <- sqrt( NumVars/N * t(MoyBlock) %*% MoyBlock / input$Sigma2E ))

    (MoyVar <- PrAdits %>% group_by(Variety) %>% summarize(moy = mean(predits)) %>% select(moy) %>% - input$Mu %>% pull)
    (MoyVar <- Predits %>% group_by(Variety) %>% summarize(moy = mean(predits)) %>% select(moy) %>% - input$Mu %>% pull)
    (VarEffectSizeRecalc <- sqrt( input$NumBlocks/N * t(MoyVar) %*% MoyVar / input$Sigma2E ))
    
    
    ############################  
    ############################  a partir de la, les ennuis commencaient  ############################  
    ############################  
    ############################   simulations of datasets
    
    # inititalize vectors to record P-values of ANOVA
    donnees <- list()
    PvalsBlockFull <- c()
    PvalsVarietyFull <- c()
    PvalsVarietyRed <- c()
    
    # grBlockPlots <- list()
    # grVarietyPlots <- list()
    
    Moy <- list()  ## temoin positif
    
    for (i in 1:input$Iter) {
      #############  partie test positif
      # generate data for a table and plot 
      data <- rnorm(n = 100, mean = input$Mu, sd = 1.414) 
      Moy[[i]] <- mean(data)
      #############  fin partie test positif
      
      # generates UNcontrolled environmental variation. This is the RESIDUAL variation.
      e <- rnorm(N, 0, sqrt(input$Sigma2E))   # pas tout a fait exact. devrait etre des tirages independants. code qq part
           
      #     # Generate Phenotypic data, such as P = G + E + e    
      y <- predits + e
           
      # the definitive dataframes
      donnees[[i]] <- as.data.frame(cbind( y, Combinations )) 

    }  ## fin des input$Iter iterations

    ############################     
    ## Starting from here, all computations will use the list donnees that contains le simulated data
    ############################ 
    
    ## Empirical effects and empirical effect sizes
    MuObs <- lapply(donnees, function(x){mean(x$y)})  ## attention, petite astuce : x$VariableATraiter
    
    MoyVarObs <- lapply(donnees, function(x){ 
      MuObs <- mean(x$y)
      x %>% group_by(genotype) %>% summarize(moy = mean(y)) %>% select(moy) %>% - MuObs %>% pull
      })
    
    VarEffectSizeObs <- lapply(MoyVarObs, function(x){
      sqrt( input$NumBlocks/N * t(x) %*% x / input$Sigma2E )
    })
 
    MoyBlockObs <- lapply(donnees, function(x){ 
      MuObs <- mean(x$y)
      x %>% group_by(block) %>% summarize(moy = mean(y)) %>% select(moy) %>% - MuObs %>% pull
    })
    
    BlockEffectSizeObs <- lapply(MoyBlockObs, function(x){
      sqrt( NumVars/N * t(x) %*% x / input$Sigma2E )
    })
    
    
    ############################## 
    ############################## enfin la partie utile !

    # ANOVA model that account for Genotype AND blocks
    Full <- lapply(donnees, function(x){ 
      tableAnova <- with(x, summary(aov( y ~ genotype + block)))  ## syntaxe oubliee de with() ... du temps perdu
      list(
      tableAnova = tableAnova,
      PvalsVarietyFull = unlist(tableAnova)[[13]],
      PvalsBlockFull = unlist(tableAnova)[[14]] 
    )
    })
    # print("############### Full starts")
    # print(Full)
    # print("############### Full ends")
    
    
    # ANOVA model that account for only for Genotype
    Reduced <- lapply(donnees, function(x){
      tableAnova <- summary(aov( y ~ genotype, data = x))
      list(
        tableAnova = tableAnova,
        PvalsVarietyRed = unlist(tableAnova)[[9]]
      )
    })    
    # print("############### Reduced starts")
    # print(Reduced)
    # print("############### Reduced ends")

    # # extract all nth elements of list lst:
    # #    sapply(lst, "[", n)
    Pvalues <- data.frame(VarietyFull = as.vector(unlist( sapply(Full, "[", 2) )),
                          BlockFull = as.vector(unlist( sapply(Full, "[", 3) )),
                          VarietyRed = as.vector(unlist( sapply(Reduced, "[", 2) ))
                          )
    # print(Pvalues)
    
  ## Some summary stats on pvalues
    StatsPvalues <- list(
  StatsPvalVarietyFull = Pvalues %>% select(VarietyFull) %>%
    summarize(n_pval = n(),
              n_sm5 = sum(. <= 0.05),
              p_sm5 = round(n_sm5 / n_pval, 3)) ,

  StatsPvalBlockFull = Pvalues %>% select(BlockFull) %>%
    summarize(n_pval = n(),
              n_sm5 = sum(. <= 0.05),
              p_sm5 = round(n_sm5 / n_pval,3)) ,

  StatsPvalVarietyRed = Pvalues %>% select(VarietyRed) %>%
    summarize(n_pval = n(),
              n_sm5 = sum(. <= 0.05),
              p_sm5 = round(n_sm5 / n_pval, 3))
    )

    ######################################
    ######################################  lastly : the graphics !!!!!!
    ######################################
    ######################################  il semble qu'il faut generer les graphiques a **l'exterieur du eventReactive**  - donc cf plus bas
    ######################################  le dataframe Pvalues est passe en list()
    ######################################  

    
    ############################ 
    ############################  impeccable !
    # return all object as a list
    list(NumVars = NumVars, Combinations = Combinations, Xmatrix = X,
         BlockEffectSize = Blocki, VarEffectSize = Varietyi, Betas = beta,
         Predicted0 = PrAdits, Predicted = Predits, MoyBlock = MoyBlock, BlockEffectSizeRecalc = BlockEffectSizeRecalc, MoyVar = MoyVar, VarEffectSizeRecalc = VarEffectSizeRecalc,
         Donnees = donnees, MuSimulated = MuObs, 
         MoyVarietyObs = MoyVarObs, MoyBlockObs = MoyBlockObs, VarietyEffectSizeObs = VarEffectSizeObs, BlockEffectSizeObs = BlockEffectSizeObs,
         Pvalues = Pvalues, StatsPvalues = StatsPvalues,
         # Grafik1 = Grafik1,  ## je laisse comme un trace; les graphiques sont crees dans les renderPlot correspondants
         Full = Full, Reduced = Reduced
    )
    
    ##################  excellent example qui marche !!
    # # draw a random number between 1 and 100 and print it 
    # random <- sample(1:100, 1) 
    # print(paste0("The number is: ", random))   ## output in console
    # 
    # # generate data for a table and plot 
    # data <- rnorm(10, mean = 100) 
    # table <- matrix(data, ncol = 2) 
    # 
    # # create a plot 
    # Figure <- plot(1:length(data), data, pch = 16, xlab ="", ylab = "") 
    # 
    # # return all object as a list 
    # list(random = random, Figure = Figure, table = table) 
  }) ## fin eventReactive Analyses

  
  ############################   
  ### Les fonctions d'affichage des objets generes dans Analyses
  ############################   
  
  
  output$NumVarieties <- renderText({
    paste('There is a fixed number of', Analyses()$NumVars, 'varieties in the dataset')
  })
  
  output$Combin <-  renderTable({ Analyses()$Combinations })
  output$Xmatrix <- renderPrint({ Analyses()$Xmatrix })
  
  output$BlockEffectSize <- renderPrint({ 
    print(sum(Analyses()$BlockEffectSize))  # sum-to-zero non compulsory
    print(t(Analyses()$BlockEffectSize) %*% Analyses()$BlockEffectSize)  # to check  with SumBiCarre
    paste("Fixed block effects", Analyses()$BlockEffectSize) 
  })

  output$VarEffectSize <- renderPrint({ 
    print(sum(Analyses()$VarEffectSize))  # sum-to-zero non compulsory
    print(t(Analyses()$VarEffectSize) %*% Analyses()$VarEffectSize)  # to check  with SumBiCarre
    paste("Fixed variety effects", Analyses()$VarEffectSize) 
    })
  
  output$Betas <- renderTable({ paste('betas', Analyses()$Betas) })
  
  output$VarEffects <- renderPrint({
    bidule <- function(x){mean(x) - input$Mu}
    print("Mu")
    print(input$Mu)
    print("Effects of Varieties")
    aggregate(Analyses()$Predicted$predits, by = list(Analyses()$Predicted$Variety), FUN = bidule) 
    })
  
  output$BlockEffects <- renderPrint({
    bidule <- function(x){mean(x) - input$Mu}
    # print("Blocks effects")
    # print(Analyses()$BlockEffectSize)
    print("Effects of Blocks")
    aggregate(Analyses()$Predicted$predits, by = list(Analyses()$Predicted$Block), FUN = bidule) 
  })
 
  output$Predicted0 <- renderTable({Analyses()$Predicted0})
  output$Predicted <- renderTable({Analyses()$Predicted})
  
  output$BlockEffectsRecalc <- renderText({ paste('Block *effects* computed from predicted data', Analyses()$MoyBlock) })
  output$BlockEffectSizeRecalc <- renderText({ paste('Block *effect size* computed from predicted data:', Analyses()$BlockEffectSizeRecalc) })
  output$VarEffectsRecalc <- renderText({ paste('Variety *effects* computed from predicted data',  Analyses()$MoyVar) })
  output$VarEffectSizeRecalc <- renderText({ paste('Variety *effect size* computed from predicted data:', Analyses()$VarEffectSizeRecalc) })
  
  output$Donnees <- renderPrint({
 #   print(Analyses()$Donnees)  ## output on console
    Analyses()$Donnees
  })
  
  output$MuSimulated <- renderPrint({ Analyses()$MuSimulated })
  
  output$MoyVarietyObs <- renderPrint({ Analyses()$MoyVarietyObs })
  output$MoyBlockObs <- renderPrint({ Analyses()$MoyBlockObs })
  
  output$VarietyEffectSizeObs <- renderPrint({ 
    paste0("Simulated effect sizes For Varities per samples. Expected value:", input$EffSizeVariety)
    Analyses()$VarietyEffectSizeObs })
  
  output$BlockEffectSizeObs <- renderPrint({ 
    paste0("Simulated effect sizes For Blocks per samples. Expected value:", input$EffSizeBlock)
    Analyses()$BlockEffectSizeObs })
  
  output$Pvalues <- renderTable({
    # print(Analyses()$Pvalues)  ## output on console
    Analyses()$Pvalues
  }, digits = 3)

  output$StatsPvalues <- renderPrint({
    Analyses()$StatsPvalues
  })
  
  
  #####################################
  #####################################  il semble qu'il faut generer les graphiques a **l'exterieur du eventReactive**
  #####################################
  
  output$Graph1 <- renderPlot({ 
    Pvalues <- Analyses()$Pvalues
    
    ######################################
    ######################################  lastly : the graphics !!!!!!
    ######################################
    grVarietyFull <- ggplot(Pvalues) +
      aes(x = VarietyFull, y = ..count../sum(..count..), fill = as.factor(VarietyFull <= 0.05) ) + ylab('Percentage of P-values') + xlab('P-values') +
      ggtitle(paste('Variety - Full model; effect size = ',input$EffSizeVariety )) +
      scale_fill_manual("P<0.05", values = c("gray70", "green")) +
      geom_histogram( breaks = seq(0, 1, by = 0.05)) + xlim(0, 1)
    
    grBlockFull <- ggplot(Pvalues) +
      aes(x = BlockFull, y = ..count../sum(..count..), fill = as.factor(BlockFull <= 0.05) ) + ylab('Percentage of P-values') + xlab('P-values') +
      ggtitle(paste('Block - Full model; effect size = ',input$EffSizeBlock,'; # of blocks:', input$NumBlocks )) +
      scale_fill_manual("P<0.05", values = c("gray70", "green")) +
      geom_histogram( breaks = seq(0, 1, by = 0.05)) + xlim(0, 1)
    
    grVarietyRed <- ggplot(Pvalues) +
      aes(x = VarietyRed, y = ..count../sum(..count..), fill = as.factor(VarietyRed <= 0.05)  ) + ylab('Percentage of P-values') + xlab('P-values') +
      ggtitle('Variety -  model w/o block') +
      scale_fill_manual("P<0.05", values = c("gray70", "green")) +
      geom_histogram( breaks = seq(0, 1, by = 0.05)) + xlim(0, 1)
    
    grid.arrange(grVarietyFull, 
                            grBlockFull, 
                            grVarietyRed, 
                            ncol = 2, nrow = 2)
     })
  
  
  output$Graph2 <- renderPlot({ 
    RCBD <- Analyses()$Predicted
    RCBD$row <- rep( c(1:input$NumBlocks), each = Analyses()$NumVars)
    RCBD$col <- rep( c(1:Analyses()$NumVars), input$NumBlocks)
    
    desplot(RCBD, Block ~ col+row,
            out1=Block, out2 = Variety, text=Variety, cex=1, aspect=200/400)
  })
  
  output$Graph3 <- renderPlot({ 
    RCBD <- Analyses()$Predicted
    RCBD$row <- rep( c(1:input$NumBlocks), each = Analyses()$NumVars)
    RCBD$col <- rep( c(1:Analyses()$NumVars), input$NumBlocks)

    Palette <- brewer.pal(n = 9, name = "BuGn")
    desplot(RCBD, predits ~ col+row,
            out1=Block, out2 = Variety, text = Variety, cex=1, aspect=250/400, col.regions = Palette)
  })
  
 output$AnovasFull <- renderPrint({
   AnovasFull <- sapply(Analyses()$Full, "[", 1)
   print("ANOVA tables of 'Full Model'")
   AnovasFull
 })

 output$AnovasReduced <- renderPrint({
   AnovasReduced <- sapply(Analyses()$Reduced, "[", 1)
   print("ANOVA tables of 'Reduced Model'")
   AnovasReduced
 })  

} ## fin server()

####################################################
shinyApp(ui = ui, server = server) 