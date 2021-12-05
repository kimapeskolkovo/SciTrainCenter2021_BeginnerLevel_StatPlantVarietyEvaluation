#' ---
#' title:  |
#'   | Case Study  : a multi-environment trial  for yield  using RCBD
#' author:
#'   - Prof L. Gentzbittel Skoltech, Digital Agriculture Laboratory ^[l.gentzbittel@skoltech.ru]
#'   - Prof C. Ben, Skoltech, Digital Agriculture Laboratory ^[c.ben@skoltech.ru]
#' date: "April, 5th 2021 - Skoltech"
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
#' Ten genotypes are assessed at five different locations.  
#' Within-site variability is controlled (assessed) using RCBD in each site, with 4 blocks per site.
#'
#' This is a first analysis to familiarize yourself with the method. YET it is not the state-of-the-art
#' method that can be requested by reviewers or shareholders.
#' A 'modern' and detailed analysis of this type of data will be carried out in the "Advanced Course"
#' 
#' # PREPARATION OF THE WORKING INTERFACE IN R 
#' 
### I. Set working directory
#On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
#Choose the directory containing the datafile and the associated R script.

### II. Installation R packages needed for the analysis on RStudio:
#Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
#Comment #1: R package installation requires a connection to internet
#Comment #2: Once packages have been installed, no need to re-install them again when you close-open again RStudio.

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# use of the constraint 'set-to-zero' for ANOVAs ## will see later in this script
options(contrasts=c('contr.treatment','contr.poly'))
#can also use 'contr.sum' for a 'sum-to-zero' constraint

#'
#' # LOADING REQUIRED METHODS FOR ANALYSIS 
#' 
## Loading of the R packages needed for the analysis.
library(car)         # Levene's test
library(agricolae)   # Newman-Keuls & Tukeys tests
library(ggplot2)
library(dplyr)
library(openxlsx)    # to load excel files
#'
#'
#' # STARTING THE ANALYSIS 
#'
#'  
## loading data file
Produc <- read.xlsx("07_MET_beginnerLevel.xlsx", sheet = 1)
str(Produc)

## The person who typed the data was lazy and has just indicated blocks and environments with numbers.
## need to transform numbers into factors, and modify names of levels in the same operation


Produc$Gen <- factor(paste("Gen", Produc$Gen, sep = ""))
Produc$Env <- factor(paste("Env", Produc$Env, sep = ""))
Produc$Rep <- factor(paste("Block", Produc$Rep, sep = ""))
Produc$Repb <- factor(paste("Block",Produc$Repb, sep = "/"))  ## To indicate that blocks are nested in environments
Produc$Prod <- as.numeric(Produc$Prod)     ## may not be required on your computer. also a weakness of read.xlsx()
str(Produc)


############ 1. visualisations

## Production per genotype, in each environment, using raw data from blocks
x11()
(graf1 <- ggplot(Produc, aes(x = Env, y = Prod, fill = Env)) +
    geom_boxplot() +
    facet_wrap( ~ Gen)
)


## A summary : average production per site and per genotype, by creasing order (add sd and n )
Summaries <- Produc %>%
    group_by(Env, Gen) %>%
     summarise(avgProduc = mean(Prod, na.rm = TRUE),
               sdProduc = sd(Prod, na.rm = TRUE),
               nbData = n()
               ) %>%
     arrange(desc(avgProduc)) %>%
    print(n = Inf)  # to seee all data


## The standard - and UNinformative -- barplot. 
## Evidence for information loss when compared to boxplots.
x11()
(graf2 <- ggplot( Summaries, aes( x = Env, y = avgProduc, fill = Env)) +
    geom_bar(stat = "identity") +
    facet_wrap( ~ Gen))
## the ratio information/ink is very low for this figure

## A classical interaction plot:
x11()
interaction.plot(x.factor = Produc$Env, trace.factor = Produc$Gen, response = Produc$Prod,
                 type = "b",  # lines and points
                 xlab = "Sites", ylab = "Prod", trace.label = "Genotypes",
                 pch=c(1:10), cex=2, lty=1, lwd=2, col = as.numeric(Produc$Gen))
## your conclusions ?
## "Egular" and "Advanced" training sessions will provide you with methodes to go further in the analysis
## and understand the pattern of variation


###############
###############  old school analysis : use of aov()
###############


## model w/accounting for blocks in each site.
## this goes down to a CRD in a each site 

model1 <- aov( Prod ~ Gen * Env , data = Produc)
summary(model1)


## Model with blocks nested in sites 
model2 <- aov( Prod ~ Rep %in% Env + Gen * Env , data = Produc)
summary(model2)
## Check and understand Df. Why 15Df for blocks within sites ?


## INCORRECT MODEL : the breeder forget that blocks are nested within sites
## Model with blocks not nested in sites 
model2BAD <- aov( Prod ~ Rep + Gen * Env , data = Produc)
summary(model2BAD)


### HORRIBLE MODEL - TOTALLY WRONG - why ?
model2HORROR <- aov( Prod ~ Rep * Gen * Env , data = Produc)
summary(model2HORROR)


## This syntax using alternative coding of blocks is also acceptable.
## you will understand why the computation of Df is correct during the "Regular Course" 
modele3 <- aov( Prod ~ Gen * Env + Repb, data = Produc)
summary(modele3)


# test Gaussian distribution of residuals - use model2.
shapiro.test(model2$residuals)  # explore using graphics to decide if it is a concern or not really
x11()
qqnorm(model2$residuals)

# test for variance homogeneity. Hand-made Levene's test !
## Caution: there is only one numerical value par combination Genotype * site * block thus NO variance ;-)
## eg:
summary(aov( abs(model2$residuals) ~ Produc$Gen:Produc$Env:Produc$Repb) )  # no test !

## if me make the asumption the the variability within blocks is homogeneous, we can test such as :
summary( aov(abs(model2$residuals) ~ Produc$Gen:Produc$Env) )
## NOT very good. We need to explore why. 

## Let's have a look at residuals variances per site:
Produc$Residuals <- model2$residuals
head(Produc) ## to see beginining of dataframe
tail(Produc) ## to see the end

## A summary of residuals:
Produc %>%
    group_by(Env) %>%
     summarise(avgResiduals = mean(Residuals, na.rm = TRUE),  ## to check one property of residuals
               varResiduals = sd(Residuals, na.rm = TRUE)^2,
               nbData = n()
               )
### Does our hypothesis make sense ?
## Sure we would need to use a ANOVA method that authorize unequal variances
## please attend the "Regular" and "Advanced" training sessions :-)






