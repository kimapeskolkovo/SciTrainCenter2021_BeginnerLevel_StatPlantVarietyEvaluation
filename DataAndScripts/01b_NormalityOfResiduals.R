#' ---
#' title:  |
#'   | Case Study  : Normality of data? No ! Normality of the residuals of ANOVA.
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
knitr::opts_chunk$set(fig.width = 5, fig.height=6, warning=FALSE, message=FALSE,
                      tidy.opts=list(width.cutoff=80), tidy=TRUE)
# Now we resume to nominal situation
#' 
#' 
#' # CASE STUDY PRESENTATION 
#' 
#' The objective of this script is to exemplify the requirement of normality of the residuals of
#' the ANOVA model -- NOT of the raw data
#'
#' PREPARATION OF THE WORKING INTERFACE IN R 
#' 
### I. Set working directory
# On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
# Choose the directory containing the datafile and the associated R script.

### II. Possibly, installation of new R packages needed for the analysis on RStudio:
# Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
# Comment #1: R package installation requires a connection to internet
# Comment #2: Once packages have been installed, no need to re-install them 
# again when you close-open again RStudio.

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())

#'
#' # LOADING REQUIRED METHODS FOR ANALYSIS 
#' 

library(ggplot2)    # a new graphic library - Will be presented in details in 'Regular' course
library(car)        # Levene's test for homogeneity of variances 
library(agricolae)  # the Newman-Keuls test for multiple mean comparisons

#'
#' # LOADING/CREATING THE DATA and STARTING ANALYSIS
#' 
set.seed(145)  # set random seed generator of computer. 
               # All participants will get the 'same' random data set

#' We  generate yield data for 100 plants of three hybrid maize varieties (G1, G2 and G3). 
#' The environmental variance $\sigma^2_E$ equals to 1.5. 
#' The phenotypic means of G1 equals 120kg/ha, of G2 equals 115kg/ha and of G3 equals 113kg/ha.
## generate data :
Yield <- data.frame( yield = c(rnorm(100, 120, 1.5),
                               rnorm(100, 115, 1.5),
                               rnorm(100, 113, 1.5)),
                    genotype = rep(c("G1", "G2", "G3"), each = 100)
                    )
( head(Yield, 12) )  ## print 12 first lines of Yield on screen
( tail(Yield, 12) )  ## print 12 last lines of Yield on screen

## histograms of yield
x11()  ## open a graphic window
## syntax for ggplot will be explained later - do not spend time to understand/recall for now
ggplot(Yield) + aes( yield) +
    geom_histogram(binwidth = 0.8, col = "black", fill = "pink", alpha = 0.2) +
    ggtitle("All yield data")

## histograms per variety
x11()
ggplot( Yield) + aes( yield) +
        geom_histogram(binwidth = 0.8, col = "black", fill = "pink", alpha = 0.2) +
    facet_grid( genotype ~ .) +
    ggtitle("Yield per genotype")

# Anova of yield data
model1 <- aov( yield ~ genotype, data = Yield)
summary(model1)  ## your conclusions ?

x11()
## quantile - quantile plot are the standard tool to visualise
## observed distribution vs expected distribution
qqnorm(residuals(model1))  ## qqnorm() draws a quantile-quantile plot for normal (gaussian) distribution
qqline(residuals(model1), col = 2)

# normality test for yield variable of the Yield dataframe
shapiro.test(Yield$yield)

# normality test for the residudals of the model
shapiro.test(residuals(model1))

# homogeneity of the variances of residuals in the different varieties
## Q: is Sigma2E the same for all varieties ?
leveneTest( aov( yield ~ genotype, data = Yield) )  ## leveneTest automagically test the residuals
## your conclusions ?

## comparing the varieties using the Neuman-Keuls post-hoc test.
SNK.test(model1, "genotype", group = TRUE, console = TRUE)
