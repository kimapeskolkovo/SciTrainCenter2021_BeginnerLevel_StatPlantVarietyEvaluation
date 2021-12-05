
#' ---
#' title:  |
#'   | Case Study  : Analysis of a seriously unbalanced design
#' author:
#'   - Prof L. Gentzbittel ^[l.gentzbittel@skoltech.ru] Digital Agriculture Laboratory - Skoltech
#'   - Prof C. Ben ^[c.ben@skoltech.ru] Digital Agriculture Laboratory - Skoltech 
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
#' Ten genotypes are assessed at two different locations.  
#' For different reasons, different numbers of measurements or replicates are performed at each site
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
library(car)         # Levene's test and types of SS
library(agricolae)   # Newman-Keuls & Tukeys tests
library(ggplot2)
library(dplyr)
library(openxlsx)    # to load excel files

#'
#' # STARTING THE ANALYSIS 
#'
## loading data file
StemLength <- read.xlsx("08_UnbalancedDesign.xlsx", sheet = 1)
str(StemLength)

## check the design
with(StemLength,
     table(Genotype, Site))  ## seriously unbalanced


## ANOVA tables of unbalanced design depends on order of factor fitted in the model !
## fit ANOVA2 mmodels, with interaction:
ana1 <- aov( StemLength ~ 1 + Genotype * Site, data = StemLength)
summary(ana1)

ana2 <- aov( StemLength ~ 1 + Site * Genotype, data = StemLength)
summary(ana2)
## clear case !

##########
## Taking unbalanceness into account using different types of sum-of-squares
##########

# create the 'most complete' model
model <- lm( StemLength ~ 1 + Site * Genotype, data = StemLength )

anova(model)  ## type I analysis, the default method for R ## OK

Anova(model, type = "II")  ## interesting

Anova(model, type = "III")  ## usual in agronomy


