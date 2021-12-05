#' ---
#' title:  |
#'   | Case Study  : A RCBD trial for wheat - a simple analysis using fixed models
#' author:
#'   - Prof L. Gentzbittel Skoltech, Digital Agriculture Laboratory ^[l.gentzbittel@skoltech.ru]
#'   - Prof C. Ben, Skoltech, Digital Agriculture Laboratory ^[c.ben@skoltech.ru]
#' date: "April, 2nd 2021 - Skoltech"
#' output: 
#'   pdf_document:
#'     keep_tex: true
#' use_bookdown: TRUE
#' latex_engine: xelatex
#' header-includes:
#'   - \usepackage{bbold}
#'   - \def\+#1{\mathbf{#1}}
#' geometry: left = 2cm, right = 1.5cm, top = 1.5cm, bottom = 1.5cm
#' ---
#' 
#' 
#+ echo = FALSE, message = FALSE, warning = FALSE
# just forgot these lines. They are used to generate the printed version
knitr::opts_chunk$set(fig.width = 4.5, fig.height = 4.5, fig.align="center",
                      warning = FALSE, message = FALSE,
                      tidy.opts=list(width.cutoff = 80), tidy = FALSE
                      )
# Now we resume to nominal situation
#' 
#' 
#' # CASE STUDY PRESENTATION 
#' 
#' Seven winter wheat cultivars were assessed for yield in a RCBD with four blocks
#' 
#' # PREPARATION OF THE WORKING INTERFACE IN R 
#' 
### I. Set working directory
# On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
# Choose the directory containing the datafile and the associated R script.

### II. Installation R packages needed for the analysis on RStudio:
# Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
# Comment #1: R package installation requires a connection to internet
# Comment #2: Once packages have been installed, no need to re-install 
# them again when you close-open again RStudio.

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# use of the constraint 'set-to-zero' for ANOVAs ## will see later in this script
options(contrasts=c('contr.treatment','contr.poly'))
# we can also use 'contr.sum' for a 'sum-to-zero' constraint
#'
#' # LOADING REQUIRED METHODS FOR ANALYSIS 
#' 
## Loading of the R packages needed for the analysis.
library(ggplot2) # Needed for some graphs (e.g. bwplots)
library(gridExtra)
library(agricolae) # For multiple mean comparisons
library(multcomp)  # for alternative multiple mean comparisons
library(multcompView)
library(emmeans)   # for alternative multiple mean comparisons
library(car)       # for Levene's test
library(openxlsx)   ## to import Excel files
#'
#'
#' # STARTING THE ANALYSIS 
#' 
###################
# Data import in R
###################

WheatYield  <- read.xlsx("02_Wheat7Var4Blocks.xlsx", sheet = 1)

WheatYield
str(WheatYield )  ## check if all columns are of the expected type: numeric, or factors, ...

## need to convert characters data into factor data -- a weakness of read.xlsx()
WheatYield$Genotype  <- as.factor(WheatYield$Genotype)
WheatYield$Block  <- as.factor(WheatYield$Block)
WheatYield$Yield <- as.numeric(WheatYield$Yield)   ## may not be required on your computer. 
## also a weakness of read.xlsx()
str(WheatYield )

attach(WheatYield )# It avoids having to specify the name of the dataframe in R commands 
## i.e. it is no more useful to write Dataframe$factor or Dataframe$variable


##Check for balanced dataset :
table(Genotype,Block)

## of course the phenotypic value is expected to be the mean of a microplot 
## or of randomly chosen plants ; not one plant !
## if several plants per plot ->  better work with the mean per plot to reduce variance
## variance of mean = variance of raw data / nbr of plants

########################
# CHECK POINT !
# Identify the factors,
# propose a practical set up of this design, in the field or in greenhouse
# which practical data (or informations) are lacking ?
#
########################

########################
# Graphic visualizations
########################

## Boxplots to reveal the distribution and variance of the measured traits 
## depending on the different factors of interest

#Individual graphs
x11()
ggplot(WheatYield) +
    aes(x = Genotype, y = Yield) +
    geom_boxplot() + geom_jitter(width = 0.1) + theme(axis.text.x = element_text(angle = 90))

x11()
ggplot(WheatYield) +
  aes(x = Block, y = Yield) +
  geom_boxplot() + geom_jitter(width = 0.10)+
  theme(axis.text.x = element_text(angle = 90))

#2 graphs on the same window. Put each graphic in an object
graf1 <- ggplot(WheatYield) +
    aes(x = Genotype, y = Yield) +
    geom_boxplot() +
    geom_jitter(width = 0.15) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_linedraw()  ## to evidence quick customisation of figures

graf2 <- ggplot(WheatYield) +
     aes(x = Block, y = Yield) +
    geom_boxplot() +
    geom_jitter(width = 0.15)+
    theme(axis.text.x = element_text(angle = 90)) +
    theme_minimal()  ## to evidence quick customisation of figures

x11()  
grid.arrange(graf1, graf2, ncol = 2, nrow = 1)  ## display the graphical objects


#xyplots if you need to see any individual performances
graf3 <- ggplot(WheatYield) +
    aes(x = Block, y = Yield) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45)) +
    facet_grid( . ~ Genotype)


graf4 <- ggplot(WheatYield) +
    aes(x = Genotype, y = Yield) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45)) +
    facet_grid( . ~ Block)

x11()
grid.arrange(graf3, graf4, ncol = 1, nrow = 2)



###############  ANOVA 
model1 <- aov( Yield ~ Block + Genotype )  
summary(model1)
## conclusions ?


#Test for ANOVA pre-requisites
#Normality of ANOVA residuals
shapiro.test(residuals(model1))
#Variance homogeneity of ANOVA residuals
leveneTest(Yield, Genotype)
leveneTest(Yield, Block)
## because only one residual per combination Block x Genotype


#### Multiple mean comparisons - Tukey HSD
print(HSD.test(model1,"Genotype"))
x11()
plot(HSD.test(model1,"Genotype"), las = 2)
## this graph is done with "base graphics" but not with "ggplot graphics". Syntax and options are different - see 'agricolae' tutorial

## Multiple mean comparisons - Newman-Keuls
print(SNK.test(model1,"Genotype"))
x11()
plot(SNK.test(model1,"Genotype"), las = 2)
## this graph is done with "base graphics" but not with "ggplot graphics". Syntax and options are different - see 'agricolae' tutorial


## Adjusted means (because the design is balanced, adjusted means = means from data)
AdjustMoys1 <- emmeans(model1,
                     pairwise ~ Genotype, 
                     adjust = "tukey")
## Options are "tukey", "scheffe", "sidak", "bonferroni", "dunnettx", "mvt", and "none"
AdjustMoys1

## Multiple comparisons. another way to compute them
## Note that the adjust= option should also be applied to the cld function if a compact letter display is requested.
(CompMoys1 <- multcomp::cld(AdjustMoys1[[1]], 
    alpha = 0.05 ,
    Letters = letters ,
    adjust = "tukey"
    )
)


########################
# CHECK POINT !
########################
#
# what if we did not define blocks ?
#
########################

ReducedModel <- aov( Yield ~ Genotype)

summary(ReducedModel)

## Multiple mean comparisons - Newman-Keuls
print(SNK.test(ReducedModel, "Genotype"))

## alternative computation
AdjustMoys4 <- emmeans(ReducedModel,  
                     pairwise ~ Genotype, 
                     adjust = "tukey")
AdjustMoys4

(CompMoys4 <- cld(AdjustMoys4[[1]],
    alpha = 0.05 ,
    Letters = letters ,
    adjust = "tukey"
    )
)
# compare to CompMoys1
## we can follow up by doing formal tests on the power of the experiment.


#' 
#' # A POINT OF THEORY : HOW ARE *really* COMPUTED ANOVA TABLES 
#' 
## aov() is fitting a GLM to test factor effects
## To explore the **reality** of computing ANOVA tables
## We will explore this point in "regular course", for unbalanced data

(SST <- t(Yield) %*% Yield)  ## Total sum-of-squares

M0 <- aov( Yield ~ 1)  ## adjusting a mean
summary(M0)

M1 <- aov( Yield ~ 1 + Genotype)  
## eq to aov( Yield ~ Genotype) allows to show we are fitting the mean AND Genotype
summary(M1)

M2 <- aov( Yield ~ 1 + Block)  
## eq to aov( Yield ~ Block) allows to shom we are fitting the mean AND Block
summary(M2)

Mfinal <- aov( Yield ~ 1 + Genotype + Block)
summary(Mfinal)
