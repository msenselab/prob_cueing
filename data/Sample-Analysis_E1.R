#################################################################################
### Sample Analysis for MS-ASP-LT Experiment 1                      #############
#################################################################################
#################################################################################
### written by Marian Sauter, February 2019 

# VARIABLE                  VALUE     MEANING
# --------------------------------------------------
# subject_nr                numeric   subject identifier (not unique!)
# Session                   numeric   session (1 or 2)
# accuracy                  numeric   accuracy, not meaningful
# response_time             numeric   response time in ms
# correct                   boolean   0 = incorrect, 1 = correct
# correct_response          char      what the correct response was (m or x/y)
# count_trial_sequence      numeric   trial counter
# block_count               numeric   block counter
# condition                 char      frequent distractor half (top or bottom)
# distractor                char      present/absent
# distractor_color          char      color as hex value
# distractor_location       numeric   location in degrees, (0 = 3 o'clock, 270 = 12 o'clock)
# distractor_orient         numeric   distractor orientation
# distractor_pos            numeric   distractor in frequent or rare half (9 = freq, 1 = rare, 0 = absent)
# distractor_type           char      orient or color distractor
# target_color              char      color as hex value
# target_identity           char      i or ! (bang)
# target_location           char      location in degrees, (0 = 3 o'clock, 270 = 12 o'clock)
# target_orient             numeric   distractor orientation
# target_pos                numeric   target in frequent or rare distractor half (9 = freq, 1 = rare)
# pid                       char      subject identifier made from distractor_type + subject_nr
# epoch                     numeric   epochs from two consecutive blocks
# experiment                char      experiment identifier

#### Set working directory and load libraries
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ez) # ezANOVA-Funktions
library(plyr) # ddply!!!
library(lsr) # cohensD
library(ggplot2) # plots
library(BayesFactor) # Bayes

#### Load helper functions
if(!exists("normDataWithin", mode="function")) source("http://www.mariansauter.de/files/R/r_functions.R")
if(!exists("summarySE", mode="function")) source("http://www.mariansauter.de/files/R/r_functions.R")
if(!exists("summarySEwithin", mode="function")) source("http://www.mariansauter.de/files/R/r_functions.R")
if(!exists("multiplot", mode="function")) source("http://www.mariansauter.de/files/R/r_functions.R")
if(!exists("article.ttest", mode="function")) source("http://www.mariansauter.de/files/R/r_functions_article.R")
if(!exists("article.ttest_full", mode="function")) source("http://www.mariansauter.de/files/R/r_functions_article.R")
if(!exists("article.ezANOVA", mode="function")) source("http://www.mariansauter.de/files/R/r_functions_article.R")

##########################################
####### Read in data
# This is the data used in Sauter et al 2019;
# ITE and td-distance effects are already removed
data <- read.table(file="data_LT_E1.csv", sep=",", dec=".", head=TRUE)

##########################################
####### Some more preprocessing

### Remove error trials 
n_filtered = nrow(data)
data <- subset(data, correct == 1)
nrow(data)/n_filtered


##########################################
##########################################
### ANALYSIS
##########################################
##########################################

#############################
### OVERALL ANOVA
#############################
# Aggregate data by distractor_type + Session + distractor_pos for each subject
# calculates median response times
agg <- ddply(data, ~ pid + distractor_type + Session + distractor_pos, .fun=summarise, 
                N = length(response_time),
                RTm = median(response_time))
 
agg$distractor_pos <- factor(agg$distractor_pos) # Set num as factor for aov
agg$Session <- factor(agg$Session) # Set num as factor for aov

# Aggregate over all participants and output table for means/sd/se/ci
(agg_td2 <- summarySEwithin(agg, measurevar="RTm", withinvars=c("Session","distractor_pos"), idvar="pid"))

# Calculate repeated-measures ANOVA
aov <- ezANOVA(agg, 
                 dv=.(RTm), 
                 within=.(distractor_pos, Session), 
                 between=.(distractor_type),
                 wid=.(pid),
                 type=2,
                 detailed = TRUE)

# Output in article-format
article.ezANOVA(aov)

#########################################
### Specific distractor intereference tests
#########################################
same <- subset(data, distractor_type == "orient") # Create subset for same-dim distractors
diff <- subset(data, distractor_type == "color") # Create subset for diff-dim distractors

# Exemplatory subsetting for different circumstances
SAME_SESSION_1 <- subset(same, Session == 1) # Define 'constants'
DIFF_SESSION_1 <- subset(diff, Session == 1) # Define 'constants'

# Set actual test object 'test' to one of the constants
# Reason: this part of the code is used multiple times for comparisons
test <- SAME_SESSION_1 # THIS IS CURRENTLY TESTED

same_freq <- subset(test, distractor_pos == 9) # frequent distractors
same_rare <- subset(test, distractor_pos == 1) # rare distractors
same_absent <- subset(test, distractor_pos == 0) # absent distractors
same_present <- subset(test, distractor_pos != 0) # present distractors

# calc median RTs for all of the above
agg_freq <- ddply(same_freq, ~ pid, summarise, RTm = median(response_time))
agg_rare <- ddply(same_rare, ~ pid, summarise, RTm = median(response_time))
agg_absent <- ddply(same_absent, ~ pid, summarise, RTm = median(response_time))
agg_present <- ddply(same_present, ~ pid, summarise, RTm = median(response_time))

# First comparison: present RTs vs. absent RTs
(ttest = t.test(agg_present$RTm, agg_absent$RTm, paired=TRUE, alternative="greater"))
d <- cohensD(agg_present$RTm, agg_absent$RTm, method = "paired")
(BF = ttestBF(agg_present$RTm, agg_absent$RTm, paired=TRUE, nullInterval=c(0,Inf)))
BFsamples <- posterior(BF, 1,iterations = 10000)
HPDinterval(BFsamples)

mean(agg_present$RTm)
mean(agg_absent$RTm)
article.ttest_full(ttest, d, BF[1]) # Output article-style

# Second comparison: frequent RTs vs. rare RTs
(ttest = t.test(agg_freq$RTm, agg_rare$RTm, paired=TRUE, alternative="less"))
d <- cohensD(agg_freq$RTm, agg_rare$RTm, method = "paired")
(BF = ttestBF(agg_freq$RTm, agg_rare$RTm, paired=TRUE, nullInterval=c(-Inf,0)))
BFsamples <- posterior(BF, 1,iterations = 10000)
HPDinterval(BFsamples)
mean(agg_freq$RTm)
mean(agg_rare$RTm)
article.ttest_full(ttest, d, BF[1]) # Output article-style

# Interference on frequent positions
inter_freq <- agg_freq$RTm - agg_absent$RTm 
mean(inter_freq)
(ttest = t.test(inter_freq))
d <- cohensD(inter_freq)
(BF = ttestBF(inter_freq,mu = 0))
BFsamples <- posterior(BF, iterations = 1000)
HPDinterval(BFsamples)
article.ttest_full(ttest, d, BF) # Output article-style

# Interference on rare positions
inter_rare <- agg_rare$RTm - agg_absent$RTm 
mean(inter_rare)
(ttest = t.test(inter_rare))
d <- cohensD(inter_rare)
(BF = ttestBF(inter_rare,mu = 0))
BFsamples <- posterior(BF, iterations = 1000)
HPDinterval(BFsamples)
article.ttest_full(ttest, d, BF) # Output article-style
