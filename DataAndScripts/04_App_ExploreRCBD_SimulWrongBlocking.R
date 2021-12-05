library(shiny) 

ui <- fluidPage(titlePanel(h3("Explore the effects of a wrong blocking strategy in RCBD designs by simulating datasets using 12 varieties and 3 blocks")),
  
  sidebarPanel(width = 2, h4("Parameters to explore"),br(),
   
   sliderInput(inputId = "EffSizeVariety", 
               label = "Effect size for varieties", 
               value = 0.7, min = 0, max = 6, step = 0.1),
   
   sliderInput(inputId = "EffSizeBlock", 
               label = "Effect size for blocks", 
               value = 0.7, min = 0, max = 6, step = 0.1),
   
   sliderInput(inputId = "Sigma2E", 
               label = "Residual variance", 
               value = 4, min = 0.5, max = 20, step = 0.5),
   
   numericInput(inputId = 'Mu', label = "Grand mean", 
                value = 100),
   
   sliderInput(inputId = "Iter", 
               label = "num. of simulated datasets", 
               value = 6, min = 2, max = 200, step = 2),
   
   actionButton("simule", "Simulate"), 
   br(), 
   br(), 
  ),
  mainPanel(h3("Datasets"),
            fluidRow(
              splitLayout(style = "border: 1px solid silver:", cellWidths = c("47%", "47%"), 
                          plotOutput("Graph2"), 
                          plotOutput("Graph3")
              ),
              splitLayout(style = "border: 1px solid silver:", cellWidths = c("47%", "47%"), 
                          plotOutput("Graph4"), 
                          plotOutput("Graph5")
              ),
              h3("Analyses"),
              plotOutput("Graph1")
            ),
    

    h3("Efficiency in detecting differences among varieties"),
    verbatimTextOutput("StatsPvalues"),
    
    # tableOutput("Combin"),
    # verbatimTextOutput("Xmatrix"), 
    
    # verbatimTextOutput("BlockEffectSize"), 
    # verbatimTextOutput("VarEffectSize"),     
    # tableOutput("Betas"),

    # verbatimTextOutput("BlockEffectsRecalc"),
    # textOutput("BlockEffectSizeRecalc"),
    # verbatimTextOutput("VarEffectsRecalc"),
    # textOutput("VarEffectSizeRecalc"),
    
    br(),
    h3("P-Values of analyses"),    
    tableOutput("Pvalues"),   
    
    br(),
    "Predicted mean values",
    tableOutput("Predicted"),
    h3("Simulated Values"),
    verbatimTextOutput("Donnees"),
    "Grand mean per simulated dataset",
    verbatimTextOutput("MuSimulated"),
    "Simulated means per Variety, as differences from grand mean of dataset",
    verbatimTextOutput("MoyVarietyObs"),
    # "Simulated effect sizes for Varieties per samples",
    # verbatimTextOutput("VarietyEffectSizeObs"),
    # "Simulated effect sizes for Blocks per samples",
    # verbatimTextOutput("BlockEffectSizeObs"),
    

    
    # br(), 
    # h3("temoin positif") ,   
    # tableOutput("Moyennes")       
    # affichage des trucs generes dans Bidules
    # verbatimTextOutput("text"), 
    # plotOutput("plot"), 
    # tableOutput("table")      
  )


) 

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
    library(gridExtra)
    library(dplyr)
    library(desplot)
    library(RColorBrewer)
    
    # Case with 3 varieties
    NumVars <- 12
    NumBlocks <- 3
    N <- NumBlocks * NumVars
    
    # ## for debugging
    # input <- list(EffSizeBlock = 0.7, EffSizeVariety = 0.7,
    #            Sigma2E = 4, Mu = 100,
    #            Iter = 4)
    
    # the controlled environmental (ie non genetic) variation
    #  EffSizeBlock <- 0.7      # to vary  cf ui()
    # the variation of the trait due to genetic factors
    #  EffSizeVariety <- 0.7    # to vary  cf ui() - relates to heritability ( Regular level )
    # the UNcontrolled environmental (ie non genetic) variation
    #  Sigma2E <- 4             # to vary cf ui()
    
    ## permute a set of 12 varieties. Too large so permute two subsets of varieties
    library(gtools)
    Lines1 <- c(paste("Variety", 1:6, sep=""))
    Lines2 <- c(paste("Variety", 7:12, sep=""))
    
    permut1 <- permutations(n=6, r=6, Lines1)
    permut2 <- permutations(n=6, r=6, Lines2)
    
    choix1 <- sample(1:dim(permut1)[1], 1) ; choix2 <- sample(1:dim(permut2)[1], 1)
    Genotypes <- c(permut1[choix1, ], permut2[choix2, ])
    
    ## the design
    Combinations <- expand.grid(genotype = Genotypes, ### permutations to create different expected datasets
                                           block = paste('Block', seq(1:NumBlocks), sep = ''))
    ## wrong block design 
    Combinations$blockbad <- c("Block1", "Block1", "Block1", "Block1",
                         "Block2", "Block2", "Block2", "Block2",
                         "Block3", "Block3", "Block3", "Block3",
                         "Block2", "Block2", "Block2", "Block2",
                         "Block3", "Block3", "Block3", "Block3",
                         "Block1", "Block1", "Block1", "Block1",
                         "Block3", "Block3", "Block3", "Block3",
                         "Block1", "Block1", "Block1", "Block1",
                         "Block2", "Block2", "Block2", "Block2"
    )
    
    ## the design matrix  
    X <- model.matrix( ~ Combinations$genotype + Combinations$block)

    ## from Cohens' f computations, we get:
    SumBiCarre <- (N/NumVars) * input$EffSizeBlock^2 * input$Sigma2E 
    SumGiCarre <- (N/NumBlocks) * input$EffSizeVariety^2 * input$Sigma2E 
    #
    # on cree trois classes effects: moins, zero, plus
    # astuce: (titi %% 3 ) -1
    titi <- seq(1:NumBlocks)
    Blocki <- ((titi %% 3) - 1 ) * sqrt(SumBiCarre / sum(abs((titi %% 3) - 1)) )
    # to check
    # sum(Blocki)  # sum-to-zero non compulsory
    # t(Blocki) %*% Blocki  # to check  with SumBiCarre
    
    toto <- seq(1:NumVars)
    Varietyi<- ((toto %% 3) - 1 )* sqrt(SumGiCarre/ (sum(abs((toto %% 3) - 1))) )
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
    Predits <- as.data.frame(cbind( predits, Combinations ))
    colnames(Predits) <- c('predits','Variety','Block', 'BlockWrong')

    ## (re)compute effects sizes using Cohen's f formulas to confirm fixed effects values are OK
    ## these computations of effect sizes can be also verified @ https://webpower.psychstat.org/models/means04/effectsize.php
    MoyBlock <- Predits %>% group_by(Block) %>% summarize(moy = mean(predits)) %>% select(moy) %>% - input$Mu %>% pull
    BlockEffectSizeRecalc <- sqrt( NumVars/N * t(MoyBlock) %*% MoyBlock / input$Sigma2E )

    MoyVar <- Predits %>% group_by(Variety) %>% summarize(moy = mean(predits)) %>% select(moy) %>% - input$Mu %>% pull
    VarEffectSizeRecalc <- sqrt( NumBlocks/N * t(MoyVar) %*% MoyVar / input$Sigma2E )
    
    
    ############################  
    ############################   simulations of datasets
    ############################  
    
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

    }  ## fin de la boucle pour generer input$Iter jeux de donnees

    print(donnees)
    
    
    ############################     
    ## Starting from here, all computations will use the list donnees that contains the simulated data
    ############################ 
    
    ## Empirical effects and empirical effect sizes
    MuObs <- lapply(donnees, function(x){mean(x$y)})  ## attention, petite astuce : x$VariableATraiter
    
    MoyVarObs <- lapply(donnees, function(x){ 
      MuObs <- mean(x$y)
      x %>% group_by(genotype) %>% summarize(moy = mean(y)) %>% select(moy) %>% - MuObs %>% pull
      })
    
    VarEffectSizeObs <- lapply(MoyVarObs, function(x){
      sqrt( NumBlocks/N * t(x) %*% x / input$Sigma2E )
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
    ############################## 
    
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
    
    FullWrong <- lapply(donnees, function(x){ 
      tableAnova <- with(x, summary(aov( y ~ genotype + blockbad)))  ## syntaxe oubliee de with() ... du temps perdu
      list(
        tableAnova = tableAnova,
        PvalsVarietyFullWrong = unlist(tableAnova)[[13]],
        PvalsBlockFullWrong = unlist(tableAnova)[[14]] 
      )
    })
    
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
                          VarietyFullWrong = as.vector(unlist( sapply(FullWrong, "[", 2) )),
                          BlockFullWrong = as.vector(unlist( sapply(FullWrong, "[", 3) )),
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
  
  StatsPvalVarietyFullWrong = Pvalues %>% select(VarietyFullWrong) %>%
    summarize(n_pval = n(),
              n_sm5 = sum(. <= 0.05),
              p_sm5 = round(n_sm5 / n_pval, 3)) ,
  
  StatsPvalBlockFullWrong = Pvalues %>% select(BlockFullWrong) %>%
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
         Predicted = Predits, MoyBlock = MoyBlock, BlockEffectSizeRecalc = BlockEffectSizeRecalc, MoyVar = MoyVar, VarEffectSizeRecalc = VarEffectSizeRecalc,
         Donnees = donnees, MuSimulated = MuObs, 
         MoyVarietyObs = MoyVarObs, VarietyEffectSizeObs = VarEffectSizeObs, BlockEffectSizeObs = BlockEffectSizeObs,
         Pvalues = Pvalues, StatsPvalues = StatsPvalues,
         # Grafik1 = Grafik1,  ## je laisse comme un trace; les graphiques sont crees dans les renderPlot correspondants
         Moyennes = Moy  ## temoin positif
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
      ggtitle(paste('Block - Full model; effect size = ',input$EffSizeBlock)) +
      scale_fill_manual("P<0.05", values = c("gray70", "green")) +
      geom_histogram( breaks = seq(0, 1, by = 0.05)) + xlim(0, 1)
    
    grVarietyFullWrong <- ggplot(Pvalues) +
      aes(x = VarietyFullWrong, y = ..count../sum(..count..), fill = as.factor(VarietyFullWrong <= 0.05) ) + ylab('Percentage of P-values') + xlab('P-values') +
      ggtitle(paste('Variety - WRONG blocking; effect size = ',input$EffSizeVariety )) +
      scale_fill_manual("P<0.05", values = c("gray70", "darkgreen")) +
      geom_histogram( breaks = seq(0, 1, by = 0.05)) + xlim(0, 1)
    
    grBlockFullWrong <- ggplot(Pvalues) +
      aes(x = BlockFullWrong, y = ..count../sum(..count..), fill = as.factor(BlockFullWrong <= 0.05) ) + ylab('Percentage of P-values') + xlab('P-values') +
      ggtitle(paste('Block - WRONG blocking; effect size = ',input$EffSizeBlock)) +
      scale_fill_manual("P<0.05", values = c("gray70", "darkgreen")) +
      geom_histogram( breaks = seq(0, 1, by = 0.05)) + xlim(0, 1)
    
    grVarietyRed <- ggplot(Pvalues) +
      aes(x = VarietyRed, y = ..count../sum(..count..), fill = as.factor(VarietyRed <= 0.05)  ) + ylab('Percentage of P-values') + xlab('P-values') +
      ggtitle('Variety -  model w/o block') +
      scale_fill_manual("P<0.05", values = c("gray70", "green")) +
      geom_histogram( breaks = seq(0, 1, by = 0.05)) + xlim(0, 1)
    
    grid.arrange(grVarietyFull, 
                 grBlockFull, 
                 grVarietyFullWrong,
                 grBlockFullWrong,
                 grVarietyRed, 
                 ncol = 2, nrow = 3)
     })
  
  
  ################ plot of RCBD
  output$Graph2 <- renderPlot({ 
    RCBD <- Analyses()$Predicted
    RCBD$row <- c(2,2,1,1,2,2,1,1,1,2,2,1,
                     4,3,3,4,4,3,3,4,4,3,4,3,
                     5,6,6,5,6,5,6,5,6,5,5,6)
    RCBD$col <- c(1,2,1,2,3,4,4,3,5,5,6,6,
                     3,4,3,4,5,5,6,6,1,1,2,2,
                     5,6,5,6,1,1,2,2,3,4,3,4)

    desplot(RCBD, Block ~ col+row,
            out1=Block, out2 = Variety, text=Variety, cex=1, aspect=200/400)
    
  })
  
  ################ plot of RCBD
  output$Graph3 <- renderPlot({ 
    RCBD <- Analyses()$Predicted
    RCBD$row <- c(2,2,1,1,2,2,1,1,1,2,2,1,
                      4,3,3,4,4,3,3,4,4,3,4,3,
                      5,6,6,5,6,5,6,5,6,5,5,6)
    RCBD$col <- c(1,2,1,2,3,4,4,3,5,5,6,6,
                      3,4,3,4,5,5,6,6,1,1,2,2,
                      5,6,5,6,1,1,2,2,3,4,3,4)

    desplot(RCBD, BlockWrong ~ col+row,
            out1=BlockWrong, out2 = Variety, text=Variety, cex=1, aspect=200/400)

  })
  
  ################ plot of RCBD
  output$Graph4 <- renderPlot({ 
    RCBD <- Analyses()$Predicted
    RCBD$row <- c(2,2,1,1,2,2,1,1,1,2,2,1,
                  4,3,3,4,4,3,3,4,4,3,4,3,
                  5,6,6,5,6,5,6,5,6,5,5,6)
    RCBD$col <- c(1,2,1,2,3,4,4,3,5,5,6,6,
                  3,4,3,4,5,5,6,6,1,1,2,2,
                  5,6,5,6,1,1,2,2,3,4,3,4)
    
    Palette <- brewer.pal(n = 6, name = "BuGn")
    desplot(RCBD, predits ~ col+row,
            out1=Block, out2 = Variety, text=Variety, cex=1, aspect=200/400, col.regions = Palette)
  })
  
  
  ################ plot of RCBD
  output$Graph5 <- renderPlot({ 
    RCBD <- Analyses()$Predicted
    RCBD$row <- c(2,2,1,1,2,2,1,1,1,2,2,1,
                  4,3,3,4,4,3,3,4,4,3,4,3,
                  5,6,6,5,6,5,6,5,6,5,5,6)
    RCBD$col <- c(1,2,1,2,3,4,4,3,5,5,6,6,
                  3,4,3,4,5,5,6,6,1,1,2,2,
                  5,6,5,6,1,1,2,2,3,4,3,4)
    
    Palette <- brewer.pal(n = 6, name = "BuGn")
    desplot(RCBD, predits ~ col+row,
            out1=BlockWrong, out2 = Variety, text=Variety, cex=1, aspect=200/400, col.regions = Palette)
    
  })
  
  
  ##########" temoin positif
  output$Moyennes <- renderTable({
    # print(Analyses()$Moyennes)  ## output on console
    unlist(Analyses()$Moyennes)
  })
  

  
  
  # ### fonctions d'affichage des trucs generes dans Bidules
  # output$text <- renderText({ 
  #   # print the random number after accessing "model" with brackets. 
  #   # It doesn't re-run the function. 
  #   youget <- paste0("After using Bidules()$random you get: ", Bidules()$random, 
  #                    ". Compare it with a value in the console") 
  #   print(youget)    ## output in console 
  #   youget ## pour ui()
  # }) 
  # 
  # output$plot <- renderPlot({ 
  #   # render saved plot 
  #   Bidules()$Figure 
  # }) 
  # 
  # output$table <- renderTable({
  #   Bidules()$table
  # })
} ## fin server()

####################################################
shinyApp(ui = ui, server = server) 
####################################################
