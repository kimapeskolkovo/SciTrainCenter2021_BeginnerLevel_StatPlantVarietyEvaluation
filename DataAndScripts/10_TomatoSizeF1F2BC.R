
###############################################################################################
###                  Scientific Training Center of Plant Biotechnologies                    ###
###                  "Estimation of the parental lines genetic value                        ###
###                            for quantitative traits"                                     ###
###                               Beginner course level                                     ###
###                    25 January 2021, Skoltech, Moscow, Russian Federation                ###
###                                                                                         ###
###                                       ----------                                        ###
###                                                                                         ###
###               Case Study  : - Genotype value for fruit diameter in tomato -             ###
###                      Analysis of heterosis, additivity and dominance                    ### 
###                        and heritabilities in the narrow and broad sense                 ###
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Pr Cecile Ben - Skoltech          ###
###############################################################################################

######### CASE STUDY PRESENTATION #############################################################

## As part of a tomato breeding program, two varieties were crossed to improve fruit diameter.
## In order to assess the interest of crossing these varieties, 
## we compare the performance of the parental varieties P1 and P2, of one F1 (P1xP2) and F2 offsprings 
## and of the two back-crosses (P1xP2) x P1 and (P1xP2) x P2.
## To this aim, the fruit diameter was evaluated under controlled environmental conditions and identical 
## for all the individuals assessed in the greenhouse.

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

 
Tomato <- read.xlsx("10_TomatoSize.xlsx", sheet = 1)

Tomato$Variety <- factor(Tomato$Variety, levels = c("P1","P2",
                                                    "F1P1xP2",
                                                    "F1xP1","F1xP2",
                                                    "F2"))
Tomato$Population <- factor(Tomato$Population, levels = c("parental","F1","BC","F2"))
Tomato$FruitSize <- as.numeric(Tomato$FruitSize) ## may not be need with your computer.

Tomato
str(Tomato)

## What are the data ?
describe(Tomato)

# To avoid typing the dataframe name
attach(Tomato)

#################
# Graphics
#################
# boxplot 
x11()
boxplot(FruitSize ~ Variety, las = 2, ylab = 'Tomato fruit diameter (mm)')

## more sophisticated plot to show differences in variances between the F2 & BC segregating populations and the non-segregation populations
x11()
ggplot(Tomato) + aes( x = FruitSize, col = Variety, fill = Variety) +
    geom_density(alpha = 0.1, position = "identity") +
    facet_grid(Population ~ .) + 
    xlab("Tomato fruit diameter (mm)") + ylab("Percent of population") +
    ggtitle("Assessment for improvment of fruit diameter in tomato")


#################
# First analysis : does fruit diameter significantly differ among populations ?
#################

ana1 <- aov(FruitSize ~ Variety)

summary(ana1)
# Conclusions ?

# Checking required conditions: residuals follow a Gaussian distribution ?
shapiro.test(residuals(ana1)) # 

# Checking required conditions: variances of the residuals are homogeneous 

## introducing Bartlett test for general comparisons of variances
#Bartlett's test (Snedecor and Cochran, 1983) is used to test if k samples have equal variances. Equal variances across samples is called homogeneity of variances. Some statistical tests, for example the analysis of variance, assume that variances are equal across groups or samples. The Bartlett test can be used to verify that assumption.
#Bartlett's test is sensitive to departures from normality. That is, if your samples come from non-normal distributions, then Bartlett's test may simply be testing for non-normality. The Levene test is an alternative to the Bartlett test that is less sensitive to departures from normality. 
bartlett.test(FruitSize ~ Variety)              # Test on raw data of fruit diameter
bartlett.test(abs(residuals(ana1)) ~ Variety)   # Test on residuals of ANOVA

# test of Brown-Forsythe, a variation of Levene's test
leveneTest(ana1, 'mean')
summary( aov(abs(residuals(ana1)) ~ Variety) )

# Genuine Levene's test : using SS of residuals to the median (not to the mean). More robust
leveneTest(ana1, 'median')  ## differences are negligeable in that example
## Your conclusions ?



###############################################################
# Let's carry out some analyses on populations with no VG : parents and F1's 
#################

NonSegregPop <- subset( Tomato,
                       subset = Population %in% c('parental', 'F1') )
## visualisation
x11()
bwplot( FruitSize ~ Variety , data = NonSegregPop,  ylab = 'Tomato fruit diameter (mm)' )

### differences among populations ?
anaNonSegreg <- aov( FruitSize ~ Variety , data = NonSegregPop)
summary(anaNonSegreg)

## Normality of residuals
shapiro.test(residuals(anaNonSegreg)) #
# Levene's test
leveneTest(anaNonSegreg)
### conclusions ? 

### A bonus: write table from R into excel format ###
write.xlsx( NonSegregPop,'TomatoHeterosis.xlsx')
write.xlsx( anaNonSegreg,'AnovaTomatoHeterosis.xlsx')


#################
# Multiple means comparisons
##  this is the KEY result to decide about heterosis and differences among parents
#################
comparisons <- TukeyHSD(anaNonSegreg) ## from base R
comparisons

print(SNK.test(anaNonSegreg, "Variety"))  ## from agricolae
print(HSD.test(anaNonSegreg, "Variety"))

## Post-hoc analysis : simple but useful figure
MultCompTest  <- SNK.test( anaNonSegreg, trt = "Variety", console = TRUE ) 
x11()
plot(MultCompTest)  

# The fruit diameters of the 3 varieties (the 2 parents and the F1 hybrid) are significantly different
# We work in controlled conditions so we consider that the effect of environmental factors is the same for all genotypes.
# Warning! The environment-genotype interaction is not taken into account
# In this hypothesis, if we find that the means of the 2 parents and of F1 for the diameter of the fruit are statistically different then the difference is due to hereditary factors.
# There is therefore genetic variability for the genes involved in controlling our trait of interest. PRE-REQUIREMENTS FOR ANY SELECTION
# The studied trait is subject to a superdominance effect.

###################################
# Computation of dominance ratio  -- need to compute the values of 'a' and 'd'
###################################

## Computing observed means per population
( ObsMeans <- Tomato %>%
      group_by(Variety) %>%
      summarise(ObsMean = mean(FruitSize, na.rm = TRUE),  
                nbData = n()
      )
)

## the method is exactly the same as for the computation of variances
## we use a system of simultaneaou equations based on formulas of parametric means

(Xmeans <- matrix(c(1,1,0,       ## P1
                    1,-1,0,      ## P2
                    1,0,1,       ## F1
                    1,0.5,0.5,   ## BC F1xP1
                    1,-0.5,0.5,  ## BC F1xP2
                    1,0,0.5),    ## F2
                  byrow = TRUE, nrow = 6, ncol = 3)
)

( PopMeans  <- ObsMeans %>% select(ObsMean) %>% pull() )

## Least Square estimates of mu, 'a' and 'd'
BetaMoy <- solve( t(Xmeans) %*% Xmeans ) %*% t(Xmeans) %*% PopMeans
BetaMoy
## Alternatively (see above)
BetaMoy1 <- qr.solve(Xmeans, PopMeans)
names(BetaMoy1) <- c("mu","a","d")


#Degree of dominance: d/a

D2 <- BetaMoy1[3] / BetaMoy1[2]
D2
### Do we reveal super-dominance ? Is it consistent with heterosis hypothesis ?

##########################
### in the "Regular" level training, we will  pursue this analysis.
##
## we will test if the additive-dominance model is sufficient or if we would need a more detailed model
## we will compute the confidence intervals for 'a' and 'd' and other exciting things
##########################


###################################
# Computations of heritabilities   - Broad sense heritability
###################################

## Computing observed variances per population
( ObsVariances <- Tomato %>%
    group_by(Variety) %>%
     summarise(ObsVar = sd(FruitSize, na.rm = TRUE)^2, ## var = sd^2
               nbData = n()
               )
)

## VP in non-segragating population is only due to VE. Take a average estimate of VE
VarE <- ObsVariances %>% filter( !Variety %in% c("F1xP1","F1xP2", "F2")) %>% summarize( VarE = mean(ObsVar) ) %>% pull()
VarE

## Use VP of the F2 to compute H2  -- see below for reason to work with F2 at this stage
( VP <- ObsVariances$ObsVar[6] )
( VarG <- VP - VarE )

H2 <- (VarG / VP )*100
H2


###################################
# Computations of heritabilities   - Narrow sense heritability
###################################
##
## This requires to estimate VA.
## To get the best estimation of VA and VD, we will use the information of ALL populations.
## We need to resolve a system of simultaneous equations
##
## the best approach is to use some tricks of matrix algebra to get the values

###################################
# Variance estimations:  VA and VD
###################################

ObsVariances
## PopVars contains the different observed VPs of the populations
( PopVars <- ObsVariances %>% select(ObsVar) %>% pull() )

## Xvars is a representation of the system of simultaneous equations
Xvars <- matrix(c(1, 0, 0, 0, 	    #P1
                 1, 0, 0, 0,		#P2
                 1, 0, 0, 0,		#F1 P1xP2
                 1, 0.5, 1, -1,	    #BC1 F1xP1
                 1, 0.5, 1, 1,		#BC2 F1xP2
                 1, 1, 1, 0		    #F2
                 ),
               byrow = TRUE, nrow = 6, ncol = 4)

## Least Square estimation of the variances -- Good old formula to know by heart !
BetaVar <- solve( t(Xvars) %*% Xvars ) %*% t(Xvars) %*% PopVars
row.names(BetaVar)  <-  c("VE","VA","VD","VAD")  ## add names for more legibility
BetaVar

## It exists a R function to compute the above Least Square estimation  
( BetaVar1 <- qr.solve(Xvars, PopVars) )
names(BetaVar1) <-  c("VE","VA","VD","VAD")

###################################
###  Alternate computation: to use the pooled estimate of VE from the non-segregating populations

## Xvar2 is a representation of the system of simultaneous equations
Xvars2 <- matrix(c(1,0,0,0, 	    # mean estimate from P1, P2, F1 P1xP2
                  1,0.5,1,-1,	#BC1 F1xP1
                  1,0.5,1,1,    #BC1 F1xP2
                  1,1,1,0		#F2
                  ),
                byrow = TRUE, nrow = 4, ncol = 4)

## Values of observed VE and VP
( PopVars2 <- c(VarE, ObsVariances %>% filter( !Variety %in% c("P1","P2","F1P1xP2") ) %>% select(ObsVar) %>% pull() ) )

BetaVar2 <- qr.solve(Xvars2, PopVars2)
names(BetaVar2) <-  c("VE","VA","VD","VAD")

### compare both computations
BetaVar1
BetaVar2  ## your comments ?

#######
####### Narrow-sense heritability
h2 <- BetaVar1[2] / ( BetaVar1[1] + BetaVar1[2] + BetaVar1[3] )  ## neglecting VAD
h2
## can be written also
BetaVar1[2] / sum( BetaVar1[-4] )  ## neglecting VAD


####### Broad-sense heritability
H2new <- ( BetaVar1[2] + BetaVar1[3] ) / sum( BetaVar1[-4] )
H2new
## previously computed value for H2:
H2




##########################
### in the "Regular" level training, we will  pursue this analysis.
##
## we will test if the additive-dominance model is sufficient or if we would need a more detailed model
## we will compute the confidence intervals for 'a' and 'd' and other exciting things
##########################
