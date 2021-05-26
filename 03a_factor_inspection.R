############### STEP 3a: Visual Inspection of Time Courses ###############
# This script enables visual inspection of the factors.
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

## Check if necessary packages are installed and if not
# install them
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(RColorBrewer)) install.packages("RColorBrewer")
if(!require(reshape)) install.packages("reshape")
if(!require(xlsx)) install.packages("xlsx")
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(xlsx)

## Load data
load("results/01_data_import/erpdata.Rdata")


######### Plot the Time Courses (i.e., factor loadings) #########

for (iFile in c("results/02bc_rotation_score/rotfit_ad23_geomin0.01.Rdata",
                "results/02bc_rotation_score/rotfit_ch21_geomin0.01.Rdata")){
  
  ## Load the results file from the previous script.
  # Please remember that this script needs to be run once for each group.
  load(iFile)
  
  ## Create a separate data object containing the latency (for the x-axis)
  # and the unstandardized loadings.
  loadings = data.frame(lat, Factor = unclass(rotFit$loadings) * rotFit$varSD)
  
  ## Reshape the data set into wide format (for the plotting function)
  loadings = melt(loadings, id.vars = "lat", variable_name = "Factor")
  
  ## Mark the factors which shall be highlighted in the loading figure
  # Based on all results we labeled the factors as follows:
  ## Adults
  # "Factor 2" P2
  # "Factor 5" early P3a
  # "Factor 3" late P3a
  # "Factor 1" LDN
  ## Children
  # "Factor 2" P2
  # "Factor 6" earlyP3a
  # "Factor 5" late P3a
  # "Factor 1" LDN
  
  
  ### This marks the factors that should be highlighted in the data.
  if (rotFit$group == "ad"){
    loadings$highlight <- factor(ifelse(loadings$Factor %in% paste0("Factor.", c(1,2,3,5)),
                                        yes = 1,
                                        no = 0))
  } else if (rotFit$group == "ch") {
    loadings$highlight <- factor(ifelse(loadings$Factor %in% paste0("Factor.", c(1,2,5,6)),
                                        yes = 1,
                                        no = 0))
  } else {
    stop("Unknown group.")
  }
  
  factorColors = colorRampPalette(brewer.pal(12, "Paired"))(efaFit$factors)
  
  ggplot(data = loadings, aes(x = lat, y = value, color = Factor, size = highlight)) +
    geom_line() + 
    scale_size_manual(values = c(0.7,1.5)) +
    guides(size = FALSE, color = guide_legend(ncol = 2)) +
    ylim(-1,5) + 
    xlab("Time [ms]") +
    ylab("Unstandardized Loadings") +
    labs(title = ifelse(efaFit$group == "ad", "Adult EFA Geomin (0.01)", "Child EFA Geomin (0.01)")) + 
    theme_classic() +
    scale_color_manual(values = factorColors) 
  
  ggsave(paste0("results/03a_factor_inspection/rotfit_", efaFit$group, efaFit$factors, "_geomin0.01.eps"), device = "eps",
         width = 6, height = 3.5) 
}

################  Match factors #################

## load both EFA results to match factors 
load("results/02bc_rotation_score/rotfit_ch21_geomin0.01.Rdata")
rotFit -> rotFit_childEFA
scores -> scores_childEFA
load("results/02bc_rotation_score/rotfit_ad23_geomin0.01.Rdata")
scores -> scores_adultEFA
rotFit -> rotFit_adultEFA

## compute average factor scores across participants
avr_scr_adEFA  <- aggregate(. ~ cond + chan, data = scores_adultEFA, FUN = mean)
avr_scr_chEFA  <- aggregate(. ~ cond + chan, data = scores_childEFA, FUN = mean)

####### Match time courses and topographies #######
## Correlation of factor loadings
# The rows are factors from the adultEFA, the columns from the childEFA.
# Higher correlations indicate higher similarity.
Rloadings <- cor(rotFit_adultEFA$loadings, rotFit_childEFA$loadings)

## Correlation of factor scores 
# Note: The excluded columns contain index variables for electrode sites etc.
# and need to be removed for the correlations to be meaningful.
Rtopo_sta <- cor(avr_scr_adEFA[avr_scr_adEFA$cond == "sta",-c(1:4)],
                 avr_scr_chEFA[avr_scr_chEFA$cond == "sta",-c(1:4)])

Rtopo_nov <- cor(avr_scr_adEFA[avr_scr_adEFA$cond == "nov",-c(1:4)],
                 avr_scr_chEFA[avr_scr_chEFA$cond == "nov",-c(1:4)])


rownames(Rloadings) <- rownames(Rtopo_sta) <- rownames(Rtopo_nov) <- paste0("adFA", 1:nrow(Rloadings))
colnames(Rloadings) <- colnames(Rtopo_sta) <- colnames(Rtopo_nov) <- paste0("chFA_", 1:ncol(Rloadings))

#### CONVENIENCE FUNCTIONS
source("99_inspection_tools.R")

# In these matrices, the rows still represent Factors from the adult EFA
# but the similarities in the columns were sorted
# and the respective factor was named.
# Example:
#        [,1]           [,2]           [,3]                     
# adFA1  "chF1_r:0.95"  "chF11_r:0.42" "chF13_r:0.39"
# This means that adult Factor 1 is most similar to child Factor 1 (followed by 11 and 13),
# with correlations of 0.95, 0.42 and 0.39

sort_similarities(similarities = Rloadings)
sort_similarities(similarities = Rtopo_sta)
sort_similarities(similarities = Rtopo_nov)


#### EXPORT TO EXCEL AS TABLES
## Format matrices for export
Rloadings <- apply(Rloadings, MARGIN = 1:2, function(x) sprintf("%.2f",x))
Rtopo_sta <- apply(Rtopo_sta, MARGIN = 1:2, function(x) sprintf("%.2f",x))
Rtopo_nov<- apply(Rtopo_nov, MARGIN = 1:2, function(x) sprintf("%.2f",x))


# Export to Excel-Files
write.xlsx(x = Rloadings[1:6,1:6], file = "results/03a_factor_inspection/Tables.xlsx", sheetName = "Table1")
write.xlsx(x = Rtopo_sta[1:6,1:6], file = "results/03a_factor_inspection/Tables.xlsx", sheetName = "Table2", append = TRUE)
write.xlsx(x = Rtopo_nov[1:6,1:6], file = "results/03a_factor_inspection/Tables.xlsx", sheetName = "Table3", append = TRUE)



