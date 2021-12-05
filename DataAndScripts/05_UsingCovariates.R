
###############################################################################################
###                  Scientific Training Center of Plant Biotechnologies                    ### 
###                  "Experimental design and biostatistics applied                         ###
###                      to plant breeding and biotechnologies"                             ###
###                               Beginner course level                                     ### 
###                    20-22 January 2021, Skoltech, Moscow, Russian Federation             ### 
###                                                                                         ### 
###                                       ----------                                        ###
###                                                                                         ###
###     Case Study  : on the use of ratios or covariates for 'standardising' measurements   ###
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Pr Cecile Ben - Skoltech          ###
###############################################################################################

######### CASE STUDY PRESENTATION #############################################################

## Data for evaluation of the activity of one enzyme in diseased and control cases.
## The totalN per sample was also obtained, to account for putative different
## quantities of enzymes per sample.


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
library(ggplot2)  ## graphic library
library(agricolae)
library(openxlsx)

# Loading data
Biochem <- read.xlsx("05_BiochemistryData.xlsx", sheet = 1)

Biochem$Status  <- as.factor(Biochem$Status)
str(Biochem)


## correlation between total N / activity
with(Biochem,
     cor(EnzymeActivity, totalN) )  ## not outstanding ...


## POST strategy
Post <- lm(EnzymeActivity ~ Status, data = Biochem)
anova(Post)
summary(Post)    ## see Adjusted R-squared


## FRACTION strategy 1
Ratio1 <- lm( EnzymeActivity / totalN ~ Status, data = Biochem)
anova(Ratio1)
summary(Ratio1)  ## see Adjusted R-squared


## FRACTION strategy 2  -- just to show it is identical as previous one
Ratio2 <- lm( (totalN - EnzymeActivity) / totalN ~ Status, data = Biochem)
anova(Ratio2)
summary(Ratio2)  ## see Adjusted R-squared


## CHANGE strategy is NOT meaningful in this case. Why ?


## ANCOVA/GLM strategy
Ancova1 <- lm( EnzymeActivity ~ Status +  totalN, data = Biochem )
anova(Ancova1)
summary(Ancova1)  ## see Adjusted R-squared

Biochem$Ancova1 <- predict(Ancova1)  ## predicted values ie w/o residuals


## ANCOVA/GLM strategy - different slopes for control and infected
Ancova2 <- lm( EnzymeActivity ~ Status * totalN, data = Biochem )
anova(Ancova2)
summary(Ancova2)  ## see Adjusted R-squared

Biochem$Ancova2 <- predict(Ancova2)  ## predicted values 


## --------------------------- A figure is worth than 1000 words   ---

x11()
ggplot(Biochem, aes(x = totalN, y = EnzymeActivity, color = Status)) +
    geom_point( cex = 3 ) + 
    geom_line(aes( y = Ancova1), linewidth = 1.5) +
    theme_bw()


x11()
ggplot(Biochem, aes(x = totalN, y = EnzymeActivity, color = Status)) +
    geom_point( cex = 3 ) + 
    geom_line(aes( y = Ancova2), linewidth = 1.5) +
    theme_bw()

## -------------------------------------------------------------

## How these models describe the data ? 
summary(Post)$adj.r.squared
summary(Ratio1)$adj.r.squared
summary(Ancova1)$adj.r.squared
summary(Ancova2)$adj.r.squared

## is Ancova2 signif better than Ancova1 ?
anova(Ancova1, Ancova2)
## When you use anova(lm.1, lm.2), it performs the F-test to compare lm.1 and lm.2
## (i.e. it tests whether reduction in the residual sum of squares are statistically significant or not). Note that this makes sense only if lm.1 and lm.2 are nested models.

## -------------------------------------------------------------

Comparison1 <- SNK.test(Post, "Status")    # POST as variable
Comparison2 <- SNK.test(Ratio1, "Status")  # FRACTION  as variable
Comparison4 <- SNK.test(Ancova1, "Status") # Ancova1 as variable
Comparison5 <- SNK.test(Ancova2, "Status") # Ancova2 as variable


x11()
par(mfrow=c(2,3))
plot(Comparison1, main = "POST method")
plot(Comparison2, main = "FRACTION method")
plot(Comparison4, main = "ANCOVA simple method")
plot(Comparison5, main = "ANCOVA interaction method")



## ##### if needed, to get values of SD etc ....
## ## multiple means comparison with POST as variable
## CompMoy1 <- lsmeans(Post,
##                     pairwise ~ Status,
##                     adjust = "tukey")
## cld(CompMoy1[[1]],
##      alpha = 0.05,
##      Letters = letters,      ### Use lower-case letters for .group
##      adjust = "tukey")       ### Tukey-adjusted comparisons


## ## multiple means comparison with FRACTION  as variable
## CompMoy2 <- lsmeans(Ratio1,
##                     pairwise ~ Status,
##                     adjust = "tukey")
## cld(CompMoy2[[1]],
##      alpha = 0.05,
##      Letters = letters,      ### Use lower-case letters for .group
##      adjust = "tukey")       ### Tukey-adjusted comparisons


## ## multiple means comparison with ANCOVA1
## CompMoy4 <- lsmeans(Ancova1,
##                     pairwise ~ Status,
##                     adjust = "tukey")
## cld(CompMoy4[[1]],
##      alpha = 0.05,
##      Letters = letters,      ### Use lower-case letters for .group
##      adjust = "tukey")       ### Tukey-adjusted comparisons


## ## multiple means comparison with ANCOVA2
## CompMoy5 <- lsmeans(Ancova2,
##                     pairwise ~ Status,
##                     adjust = "tukey")
## cld(CompMoy5[[1]],
##      alpha = 0.05,
##      Letters = letters,      ### Use lower-case letters for .group
##      adjust = "tukey")       ### Tukey-adjusted comparisons
## ## Conclusions ?
