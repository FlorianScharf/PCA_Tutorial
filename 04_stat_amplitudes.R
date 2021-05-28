############### STEP 4: Statistical analysis of amplitude effects ###############
# This script runs inferential statistics on the extracted factors.
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

# empty workspace
rm(list=ls())

## Load necessary packages 
if(!require(afex)) install.packages("afex")
if(!require(BayesFactor)) install.packages("BayesFactor")
if(!require(bayestestR)) install.packages("bayestestR")

library(afex)
library(psych)
library(BayesFactor)
library(bayestestR)
library(ggplot2)

## Load the data from both EFAs and put the results into separate objects
load("results/02bc_rotation_score/rotfit_ad23_geomin0.01.Rdata",  temp_env <- new.env())
adEFA <- as.list.environment(temp_env)
load("results/02bc_rotation_score/rotfit_ch21_geomin0.01.Rdata",  temp_env <- new.env())
chEFA <- as.list.environment(temp_env)

## Compute the dependent variables for statistical analyses 

## Adult EFA
# Compute unstandardized factor loadings
adEFA$rotFit$loadings_unstd <- adEFA$rotFit$loadings * adEFA$rotFit$varSD
# Compute the peak loading for each factor
adEFA$peak_loadings <- apply(adEFA$rotFit$loadings_unstd, MARGIN = 2, FUN = max)
# Multiply peak loadings with factor scores
adEFA$scores[,-c(1:4)] = as.matrix(adEFA$scores[,-c(1:4)]) %*% diag(adEFA$peak_loadings)

## Child EFA
# Compute unstandardized factor loadings
chEFA$rotFit$loadings_unstd <- chEFA$rotFit$loadings * chEFA$rotFit$varSD
# Compute the peak loading for each factor
chEFA$peak_loadings <- apply(chEFA$rotFit$loadings_unstd, MARGIN = 2, FUN = max)
# Multiply peak loadings with factor scores
chEFA$scores[,-c(1:4)] = as.matrix(chEFA$scores[,-c(1:4)]) %*% diag(chEFA$peak_loadings)

## Relable the factors which we labeled in the text.
# Based on all results we labeled the factors as follows:
## Adults
# "Factor 1" LDN
# "Factor 2" P2
# "Factor 3" late P3a
# "Factor 5" early P3a

## Children
# "Factor 1" LDN
# "Factor 2" P2
# "Factor 5" late P3a
# "Factor 6" earlyP3a


# Which factors need to be relabeled?
highlighted_FA_ad <- names(adEFA$scores) %in% paste0("X", c(1,2,3,5))
highlighted_FA_ch <- names(chEFA$scores) %in% paste0("X", c(1,2,5,6))

# Rename the respective columns
names(adEFA$scores)[highlighted_FA_ad] <- c("LDN", "P2", "lP3a", "eP3a")
names(chEFA$scores)[highlighted_FA_ch] <- c("LDN", "P2", "lP3a", "eP3a")

# Create a new data set containing both groups and only the relevant data
scores = rbind(chEFA$scores[, c("group", "cond", "subj", "chan", "P2", "eP3a", "lP3a", "LDN")], adEFA$scores[, c("group", "cond", "subj", "chan", "P2", "eP3a", "lP3a", "LDN")])

# Compute difference scores for plotting purposes
scores_diff = scores[scores$cond == "nov",]
scores_diff[, 5:8] = scores[scores$cond == "nov", 5:8] - scores[scores$cond == "sta", 5:8]
scores_diff[, "cond"] = "nov-sta"

# Export final file, if you prefer, you can analyze the scores in a different
# software such as JASP or SPSS.
write.csv(scores, file = "results/04_stat_amplitudes/final_scores.csv")

########### Factorwise-statistical analysis ###########
##### Factor P2 #####

# Descriptive Statistics
describeBy(P2 ~ cond + group, data = scores[scores$chan == "Cz",], mat = TRUE, digits = 3)

# Frequentist ANOVA
fit <- aov_ez(data = scores[scores$chan == "Cz",], dv = "P2", id = "subj", between = "group", within = "cond", type = 3)
fit

# Bayes ANOVA
res = anovaBF(P2 ~ (cond * group) + subj, data = scores, 
              whichRandom = "subj", 
              whichModels = "withmain", 
              rscaleRandom = "nuisance", 
              rscaleFixed = "medium", iterations = 10000)
sort(res)
bayesfactor_inclusion(res, match_models = TRUE)

# Frequentist Pairwise comparisons:
#  nov vs. sta 
#  for Adults
t.test(scores[scores$group == "ad" & scores$chan == "Cz" & scores$cond == "nov",]$P2, 
scores[scores$group == "ad" & scores$chan == "Cz" & scores$cond == "sta",]$P2, 
paired = TRUE)

#  for Children
t.test(scores[scores$group == "ch" & scores$chan == "Cz" & scores$cond == "nov",]$P2, 
       scores[scores$group == "ch" & scores$chan == "Cz" & scores$cond == "sta",]$P2, 
       paired = TRUE)


# Bayesian Pairwise comparisons:
#  nov vs. sta 
#  for Adults
ttestBF(scores[scores$group == "ad" & scores$chan == "Cz" & scores$cond == "nov",]$P2, 
        scores[scores$group == "ad" & scores$chan == "Cz" & scores$cond == "sta",]$P2, 
        paired = TRUE)

#  for Children
ttestBF(scores[scores$group == "ch" & scores$chan == "Cz" & scores$cond == "nov",]$P2, 
        scores[scores$group == "ch" & scores$chan == "Cz" & scores$cond == "sta",]$P2, 
        paired = TRUE)

## Plot interaction
ggplot(scores[scores$chan == "Cz",], aes(x = group, y = P2, color = cond)) +
        geom_violin(position = position_dodge(0.8)) +
        geom_boxplot(width = 0.1, position = position_dodge(0.8)) +
        scale_y_reverse() +
        labs(x = "Age group", y = "Amplitude [µV]", color = "Sound type")

ggplot(scores_diff[scores_diff$chan == "Cz",], aes(color = group, y = P2, x = group)) +
        geom_violin(position = position_dodge(0.8)) +
        geom_boxplot(width = 0.1, position = position_dodge(0.8)) +
        scale_y_reverse() +
        labs(x = "Age group", y = "Nov - Sta difference amplitude [µV]", color = "Age group")

##### Factor eP3a #####

# Descriptive Statistics
fit <- describeBy(eP3a ~ cond + group, data = scores[scores$chan == "Cz",], mat = TRUE, digits = 3)
fit 

# Frequentist ANOVA
aov_ez(data = scores[scores$chan == "Cz",], dv = "eP3a", id = "subj", between = "group", within = "cond", type = 3)

# Bayes ANOVA
res = anovaBF(eP3a ~ (cond * group) + subj, data = scores, 
              whichRandom = "subj", 
              whichModels = "withmain", 
              rscaleRandom = "nuisance", 
              rscaleFixed = "medium", iterations = 10000)
sort(res)
bayesfactor_inclusion(res, match_models = TRUE)

# Frequentist Pairwise comparisons:
#  nov vs. sta 
#  for Adults
t.test(scores[scores$group == "ad" & scores$chan == "Cz" & scores$cond == "nov",]$eP3a, 
       scores[scores$group == "ad" & scores$chan == "Cz" & scores$cond == "sta",]$eP3a, 
       paired = TRUE)

#  for Children
t.test(scores[scores$group == "ch" & scores$chan == "Cz" & scores$cond == "nov",]$eP3a, 
       scores[scores$group == "ch" & scores$chan == "Cz" & scores$cond == "sta",]$eP3a, 
       paired = TRUE)


# Bayesian   Pairwise comparisons:
#  nov vs. sta 
#  for Adults
ttestBF(scores[scores$group == "ad" & scores$chan == "Cz" & scores$cond == "nov",]$eP3a, 
        scores[scores$group == "ad" & scores$chan == "Cz" & scores$cond == "sta",]$eP3a, 
        paired = TRUE)

#  for Children
ttestBF(scores[scores$group == "ch" & scores$chan == "Cz" & scores$cond == "nov",]$eP3a, 
        scores[scores$group == "ch" & scores$chan == "Cz" & scores$cond == "sta",]$eP3a, 
        paired = TRUE)


## Plot interaction
ggplot(scores[scores$chan == "Cz",], aes(x = group, y = eP3a, color = cond)) +
        geom_violin(position = position_dodge(0.8)) +
        geom_boxplot(width = 0.1, position = position_dodge(0.8)) +
        scale_y_reverse() +
        labs(x = "Age group", y = "Amplitude [µV]", color = "Sound type")


ggplot(scores_diff[scores_diff$chan == "Cz",], aes(color = group, y = eP3a, x = group)) +
        geom_violin(position = position_dodge(0.8)) +
        geom_boxplot(width = 0.1, position = position_dodge(0.8)) +
        scale_y_reverse() +
        labs(x = "Age group", y = "Nov - Sta difference amplitude [µV]", color = "Age group")


##### Factor lP3a #####

# Descriptive Statistics
describeBy(lP3a ~ cond + group, data = scores[scores$chan == "Fz",], mat = TRUE, digits = 3)

# Frequentist ANOVA
fit <- aov_ez(data = scores[scores$chan == "Fz",], dv = "lP3a", id = "subj", between = "group", within = "cond", type = 3)
fit

# Bayes ANOVA
res = anovaBF(lP3a ~ (cond * group) + subj, data = scores, 
              whichRandom = "subj", 
              whichModels = "withmain", 
              rscaleRandom = "nuisance", 
              rscaleFixed = "medium", iterations = 10000)
sort(res)
bayesfactor_inclusion(res, match_models = TRUE)

# Frequentist Pairwise comparisons:
#  nov vs. sta 
#  for Adults
t.test(scores[scores$group == "ad" & scores$chan == "Fz" & scores$cond == "nov",]$lP3a, 
       scores[scores$group == "ad" & scores$chan == "Fz" & scores$cond == "sta",]$lP3a, 
       paired = TRUE)

#  for Children
t.test(scores[scores$group == "ch" & scores$chan == "Fz" & scores$cond == "nov",]$lP3a, 
       scores[scores$group == "ch" & scores$chan == "Fz" & scores$cond == "sta",]$lP3a, 
       paired = TRUE)


# Bayesian   Pairwise comparisons:
#  nov vs. sta 
#  for Adults
ttestBF(scores[scores$group == "ad" & scores$chan == "Fz" & scores$cond == "nov",]$lP3a, 
        scores[scores$group == "ad" & scores$chan == "Fz" & scores$cond == "sta",]$lP3a, 
        paired = TRUE)

#  for Children
ttestBF(scores[scores$group == "ch" & scores$chan == "Fz" & scores$cond == "nov",]$lP3a, 
        scores[scores$group == "ch" & scores$chan == "Fz" & scores$cond == "sta",]$lP3a, 
        paired = TRUE)

## Plot interaction
ggplot(scores[scores$chan == "Fz",], aes(x = group, y = lP3a, color = cond)) +
        geom_violin(position = position_dodge(0.8)) +
        geom_boxplot(width = 0.1, position = position_dodge(0.8)) +
        scale_y_reverse() +
        labs(x = "Age group", y = "Amplitude [µV]", color = "Sound type")

ggplot(scores_diff[scores_diff$chan == "Fz",], aes(color = group, y = lP3a, x = group)) +
        geom_violin(position = position_dodge(0.8)) +
        geom_boxplot(width = 0.1, position = position_dodge(0.8)) +
        scale_y_reverse() +
        labs(x = "Age group", y = "Nov - Sta difference amplitude [µV]", color = "Age group")


##### Factor LDN #####

# Descriptive Statistics
describeBy(LDN ~ cond + group, data = scores[scores$chan == "F4",], mat = TRUE, digits = 3)

# Frequentist ANOVA
fit <- aov_ez(data = scores[scores$chan == "F4",], dv = "LDN", id = "subj", between = "group", within = "cond", type = 3)
fit

# Bayes ANOVA
res = anovaBF(LDN ~ (cond * group) + subj, data = scores, 
              whichRandom = "subj", 
              whichModels = "withmain", 
              rscaleRandom = "nuisance", 
              rscaleFixed = "medium", iterations = 10000)
sort(res)
bayesfactor_inclusion(res, match_models = TRUE)

# Frequentist Pairwise comparisons:
#  nov vs. sta 
#  for Adults
t.test(scores[scores$group == "ad" & scores$chan == "F4" & scores$cond == "nov",]$LDN, 
       scores[scores$group == "ad" & scores$chan == "F4" & scores$cond == "sta",]$LDN, 
       paired = TRUE)

#  for Children
t.test(scores[scores$group == "ch" & scores$chan == "F4" & scores$cond == "nov",]$LDN, 
       scores[scores$group == "ch" & scores$chan == "F4" & scores$cond == "sta",]$LDN, 
       paired = TRUE)


# Bayesian Pairwise comparisons:
#  nov vs. sta 
#  for Adults
ttestBF(scores[scores$group == "ad" & scores$chan == "F4" & scores$cond == "nov",]$LDN, 
        scores[scores$group == "ad" & scores$chan == "F4" & scores$cond == "sta",]$LDN, 
        paired = TRUE)

#  for Children
ttestBF(scores[scores$group == "ch" & scores$chan == "F4" & scores$cond == "nov",]$LDN, 
        scores[scores$group == "ch" & scores$chan == "F4" & scores$cond == "sta",]$LDN, 
        paired = TRUE)

## Plot interaction
ggplot(scores[scores$chan == "F4",], aes(x = group, y = LDN, color = cond)) +
        geom_violin(position = position_dodge(0.8)) +
        geom_boxplot(width = 0.1, position = position_dodge(0.8)) +
        scale_y_reverse() +
        labs(x = "Age group", y = "Amplitude [µV]", color = "Sound type")

ggplot(scores_diff[scores_diff$chan == "F4",], aes(color = group, y = LDN, x = group)) +
        geom_violin(position = position_dodge(0.8)) +
        geom_boxplot(width = 0.1, position = position_dodge(0.8)) +
        scale_y_reverse() +
        labs(x = "Age group", y = "Nov - Sta difference amplitude [µV]", color = "Age group")



