#' ---
#' title:  |
#'   | Case Study  : a three-generation-means analysis in barley
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
#+ echo = FALSE, message = FALSE, warning = FALSE
# just forgot these lines. They are used to generate the printed version
knitr::opts_chunk$set(fig.width = 5, fig.height=6, warning=FALSE, message=FALSE,
                      tidy.opts=list(width.cutoff=80), tidy=TRUE)
# Now we resume to nominal situation
#' 
#' # CASE STUDY PRESENTATION 
#' 
#' The objective of this script is to create a short and simple script to explore a barley trial. Two parental lines, 
#' their F1 and F2 offsprings were sown in a same place. Three phenotypes are recorded on each plant:  
#' 
#' 1. number of grains per ear; 
#' 2. the presence/absence of long awns on the ear
#' 3. resistance to *Puccinia hordei*  
#' 
#' The goals are:
#' 
#'  1. to test for putative differences between the different generations
#'  2. to test for simple genetic model for both qualitative traits -- awn and resistance ; and to evaluate if they are linked
#'  3. to test if one of the morphologic trait is linked to Grains per ear.
#' 
#' # PREPARATION OF THE WORKING INTERFACE IN R 
#' 
### I. Set working directory  ####
# On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
# Choose the directory containing the datafile and the associated R script.

### II. Possibly, installation of new R packages needed for the analysis on RStudio:
# Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
# Comment #1: R package installation requires a connection to internet
# Comment #2: Once packages have been installed, 
# no need to re-install them again when you close-open again RStudio.

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())

#'
#' # LOADING REQUIRED METHODS FOR ANALYSIS 
#' 
## In this example, we will use  R-base graphics.
## We will use the newer 'ggplot2' graphic package in other examples

library(Hmisc)      ## for describe()
library(openxlsx)   ## to import Excel files
library(agricolae)  ## for Newman-Keuls

#'
#' # STARTING THE ANALYSIS
#' 
####################################
# Import of data
####################################

## before loading data, open the excel file. Inspect organisation. 
## Understand what are the factors, what are the variables.

## import from Excel to R - We will see other methods later in Regular Training
Barley <- read.xlsx("01_UsingR_BarleyPreBreeding_YieldRustEarAwns.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

# The data
Barley

# Structure of dataset -- important, to check if data import is OK
str(Barley)     ## important.

# A quick description of all columns of the dataset
describe(Barley)

# A shortcut to avoid typing dataframe name in subsequent analyses
attach(Barley)


#################
# 1.  Graphic Analysis of data
#################

# boxplot per population
x11()  ## opens a graphic window to display the figure
boxplot(GrainPerEar ~ Population,
        las = 2,  # variety names written vertical
        ylab='number of grain per ear',
        col = "grey75", border = "darkred")


## The order of population is not convenient (alphabetic order)
## We can reorganise the 'population factor in a suitable order
Population <- factor(Population, levels = c("MU302", "Thibault", "F1","F2"))
 
## Back to the boxplot to view re-ordering of names
x11()
boxplot(GrainPerEar ~ Population,
        ylab='number of grain per ear',
        col = "lightblue", border = "darkgreen")

## please note that the variability of F2 values is greater than that of parental lines and F1: 
## Any ideas why ? 

#############
# let's test if it exists a difference among populations
#############
##
################  IMPORTANT : WHY do we have the right to do this comparison ?
##
## aov() is the fonction to adjust an ANOVA model to the data
## if the present case, we will fit a model with ONE factor
##  ResultOfAnova <- aov( VariableToTest ~ Factor )  ##
##
## the model is :
##
## GrainPerEar = mu + Population Effect + residual variability

ana1 <- aov(GrainPerEar ~ Population)  
summary(ana1)
## Your conclusions ?  

## Post-hoc analysis/ what are the populations which significantly differ 
## regarding the number of grains per ear
MultCompTest  <- SNK.test( ana1, trt = "Population", console = TRUE ) 

x11()
plot(MultCompTest)  ## simple but useful figure


## Just to see if we would get the same results for the difference among the parental lines 
## using a t-test for means
GrainP1 <- GrainPerEar[ Population == 'MU302' ]  ## [ ] is the operator to subset among the data 
GrainP2 <- GrainPerEar[ Population == 'Thibault']

t.test(GrainP1, GrainP2, var.equal = FALSE)  ## in case the variances among P1 and P2 are different
## if the variances among P1 and P2 are different (or among the other populations) ... 
## this is an issue for ANOVA. We will discuss it in a few minutes)



# We can also use an ANOVA for a factor with only two levels
## subset() subset a dataframe. It is an alternative to []
toto <- subset( Barley,  subset = Population %in% c('Thibault','MU302') )  
toto
## the t-test can be related to an ANOVA with only two levels of a factor
anaParents <- aov( GrainPerEar ~ Population, data = toto)
summary(anaParents)


## ------------------  STOP HERE FOR THAT STEP !
## ------------------  At a next session, the following will be performed:
## ------------------ Statistical analysis of qualitative traits (type of segregation, possible linkage)
## ------------------ Statistical analysis of possible linkage between yield and qualitative traits


## ------------------  Qualitative characters

## awness :2 phenotypic categories in the F2 population in segregation -> 
## one locus with recessive/dominance relationships
## resistance / 2 phenotypic categories in the F2 population in segregation -> 
## one locus  with recessive/dominance relationships

###################################
#  to test AT THE SAME TIME if our genetic hypothesis is true AND if the locus are linked or not, 
## we will test the expected segregation in a F2
#
#  we observe 4 phenotypic classes in F2 ->  our basis hypothesis is :  ?
## what are the expected genetic and phenotypic formulas in F2 given our basis hypothesis ?

TwoTraits <- table( EarAwn[Population == 'F2'], Resistance[Population == 'F2'] )
TwoTraits  ## is a table

segreg <- as.vector( TwoTraits )  ## table as vector for next computations
segreg

# Test to fit the theoretical distribution, using a ChiSquare test :
chisq.test(segreg,  # observed distribution to test
           p = c( 9/16, 3/16, 3/16, 1/16))   ## expected distribution of segregation of two 
                                             ## unlinked loci with Recessive/Dominance

# this test is approximate because one case is less than 5
## conclusions ?  Is Awness a possible marker for resistance:susceptibility to brown rust ?



###################################
# Linkage between Awness and Grain per Ear ?

## Use only F2 data, so we create a new dataframe with only F2 data

F2Data <- Barley[ Population == 'F2', ] ## select all lines where population equals F2; and all columns 
F2Data


# Ear awness and Grain per Ear
# A graphic
x11()
boxplot(GrainPerEar ~ EarAwn, data = F2Data )

### test for difference using ANOVA - note there are only two levels of EarAwn 
## and a student t-test may have been sufficient

GrainAwn  <- aov(GrainPerEar ~ EarAwn, data = F2Data)
summary(GrainAwn)
## conclusions ? Is awness a good marker for high yield ?


###################################
# Linkage between Resistance and Grain per Ear ?

# A graphic
x11()
boxplot(GrainPerEar ~ Resistance, data = F2Data )

### test for difference using ANOVA - note there are only two levels of Resistance 
## and a student t-test may have been sufficient

GrainResistance  <- aov(GrainPerEar ~ Resistance, data = F2Data)
summary(GrainResistance)

## conclusions ? Is resistance a good marker for high yield  ?  
## Will it be easy to breed for resistance AND high yielding varieties ?


################  in fact :
################  you've done your first QTL detection !!



