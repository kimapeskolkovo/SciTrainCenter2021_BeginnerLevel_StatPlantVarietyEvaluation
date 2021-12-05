
###############################################################################################
###                  Scientific Training Center of Plant Biotechnologies                    ###
###                  "Estimation of the parental lines genetic value                        ###
###                            for quantitative traits"                                     ###
###                               Beginner course level                                     ###
###                    25 January 2021, Skoltech, Moscow, Russian Federation                ###
###                                                                                         ###
###                                       ----------                                        ###
###                                                                                         ###
###     Case Study  : tests for cytoplasmic inheritance and heterosis in bread wheat        ###
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Pr Cecile Ben - Skoltech          ###
###############################################################################################

######### CASE STUDY PRESENTATION #############################################################

## A company aims to obtain information on the genetic control of the quantitative level of grain protein, in particular about the possible existence
## of an heterosis effect (and more particularly of superdominance) and about possible cytoplasmic inheritance.  A comparative test in the field of the
## protein content in the grain is carried out under cultivation under low nitrogen inputs.  The performance of the 2 parental lines (Opale and E508),
## reciprocal F1 hybrids F1OxE (using the Opale variety as female parent and E508 as male parent) and F1ExO (using variety E508 as female parent and
## Opale as male parent), and segregating F2 populations from the two reciprocal F1 crosses (F2OxE and F2ExO) are evaluated under the same environmental
## conditions
## data is observed yield from 10 m2 microplots

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
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

#############################################################################################################

## Loading of the R packages needed for the analysis.
library(Hmisc)
library(agricolae)
library(car)
library(openxlsx)
library(dplyr)
library(ggplot2)

#################
# Import data
#################
 
WheatProt <- read.xlsx("09_WheatProteinCytoplasmicHeredity.xlsx",sheet = 1)

WheatProt

## reorder factor level per generation (useless alphabetic order by default)
WheatProt$Variety <- factor(WheatProt$Variety, levels = c("Opale", "E508",
                                                             "F1ExO", "F1OxE",
                                                          "F2ExO", "F2OxE" ))
WheatProt$Population <- factor(WheatProt$Population, levels = c("parent", "F1", "F2" ))

WheatProt$ProteinContent <- as.numeric(WheatProt$ProteinContent)     ## may not be required on your computer. also a weakness of read.xlsx()
str(WheatProt)   ## OK !

describe(WheatProt)

# To avoid typing the dataframe name
attach(WheatProt)

#################
# Graphic analysis
#################
x11()
boxplot(ProteinContent ~ Variety, las = 2, ylab = 'Protein content in grain')

## discuss the quality of the data set (eg existence of outliers, etc.),
## discuss the potential existence of an effect of variety on the protein level 
## Some plants are outliers. Would you select them ?  transgressive segregation or environmental effect ?

## more sophisticated plot to show differences in variances between the F2 segregating populations and the non-segregation populations
x11()
ggplot(WheatProt) + aes( x = ProteinContent, col = Variety, fill = Variety) +
    geom_density(alpha = 0.1, position = "identity") +
    facet_grid(Population ~ .) + 
    xlab("Protein content in grain (as percentage of grain mass)") + ylab("Percent of population") +
    ggtitle("Assessment for improvment of protein content in winter wheat")


#################
# First Analysis : does protein content significantly vary among populations ?
#################

ana1 <- aov(ProteinContent ~ Variety)

summary(ana1)
## your first conclusions ?

## shapiro.test(residuals(ana1)) # residuals are OK
## x11()
## qqnorm(residuals(ana1))

# test de Levene's test for homogeneity of variances
leveneTest(ana1, center = "mean")
summary(aov(abs(residuals(ana1)) ~ Variety))  ## to show how the Levene's test is computed

# Your conclusions ? Is the heterogeneity of variances expected ? What to do?


#################
# Let's carry out some analyses on populations with no VG : parents and F1's 
#################

NonSegregPop <- subset( WheatProt,
                       subset = Variety %in% c('Opale','E508','F1OxE', 'F1ExO') )
## visualisation
x11()
bwplot( ProteinContent ~ Variety , data = NonSegregPop )


### differences among populations ?
anaNonSegreg <- aov( ProteinContent ~ Variety , data = NonSegregPop)
summary(anaNonSegreg)

## Normality of residuals
shapiro.test(residuals(anaNonSegreg)) #
# Levene's test
leveneTest(anaNonSegreg)
### conclusions ? 


#################
# Multiple means comparisons
##  this is the KEY result to decide about heterosis and cytoplasmic inheritance and differences among parents
#################
print( HSD.test(anaNonSegreg, "Variety"))
print( SNK.test(anaNonSegreg, "Variety"))

## Post-hoc analysis : simple but useful figure
MultCompTest  <- SNK.test( anaNonSegreg, trt = "Variety", console = TRUE ) 
x11()
plot(MultCompTest)  

######## questions :
## • Is the difference observed in the protein content between the 2 parental varieties due to hereditary factors?
## • Is there a heterosis effect for protein content ?
## • Is there cytoplasmic inheritance controlling the agronomic trait of interest ?

    
###################################
# Computations of heritabilities 

## Computing observed variances per population
( ObsVariances <- WheatProt %>%
    group_by(Variety) %>%
     summarise(ObsVar = sd(ProteinContent, na.rm = TRUE)^2,  ## to check one property of residuals
               nbData = n()
               )
)

## VP in non-segregating population is only due to VE. Take an average estimate
VarE <- ObsVariances %>% filter( Variety %in% c('Opale','E508','F1OxE', 'F1ExO')) %>% summarize( VE = mean(ObsVar) )
VarE
VarE <- ObsVariances %>% filter( !Variety %in% c("F2ExO","F2OxE")) %>% summarize( VE = mean(ObsVar) )
VarE

## VP = VG + VE thus
## VG = VP - VE
# Let's consider the VP of the F2 created from the most efficient F1 : F1OxE
VP  <- ObsVariances$ObsVar[6]  ## 6th line of results of the column ObsVar
VP
### so 
VarG <- VP - VarE
VarG

### Broad sense of heritability (not the one you can use for breeding)
### Caution ! you should NOT use the VP of all data, or the VP of all F2
## BUT the VP used to compute the VG

H2 <- ( VarG / VP ) * 100  
H2


