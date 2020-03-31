library(ez) # ezANOVA-Funktions
library(tidyverse)

source('fun.R')

data_e1 <- read.table(file="./data/data_LT_E1.csv", sep=",", dec=".", head=TRUE)
data_e2 <- read.table(file="./data/data_LT_E2.csv", sep=",", dec=".", head=TRUE)

data_e1$distractor_pos <- factor(data_e1$distractor_pos, levels=c(0, 1, 9),
                                 labels=c("Absent", "Rare region", "Freq. region"))
data_e2$distractor_pos <- factor(data_e2$distractor_pos, levels=c(0, 1, 9), 
                                 labels=c("Absent", "Rare region", "Freq. region"))

# Distractor location effect
data_e1 %>% 
  filter(correct==1, Session==1, distractor_type=="orient") %>% 
  group_by(distractor_pos, block_count) %>% summarize(mRT=mean(response_time)) %>% 
  ggplot(aes(x=block_count, color=distractor_pos, group = distractor_pos, y=mRT)) + 
  geom_point() + geom_line() + theme_bw() + 
  labs(x="Block count", y = "Reaction time (ms)", color = "Distractor position")

# Distractor position repetition effect
data_e1 %>% mutate(dprep=ifelse(distractor == "present",
                                ifelse(distractor_location == lag(distractor_location), 
                                       "Repeat", "No repeat"), "Absent")) %>% 
  filter(!is.na(dprep), correct==1, Session==1, distractor_type=="orient") %>% 
  group_by(dprep, block_count) %>% summarize(mRT=mean(response_time)) %>% 
  ggplot(aes(x=block_count, color=dprep, group = dprep, y=mRT)) + geom_point() + geom_line() +
  theme_bw() + labs(x="Block count", y = "Reaction time (ms)", color = "Distractor position")

# Target position repetition effect
data_e1 %>% mutate(tprep = ifelse(target_location == lag(target_location), 
                                  "Repeat", "No repeat")) %>% 
  filter(!is.na(tprep), correct==1, Session==1, distractor_type=="orient") %>% 
  group_by(tprep, block_count) %>% summarize(mRT=mean(response_time)) %>% 
  ggplot(aes(x=block_count, color=tprep, group = tprep, y=mRT)) + geom_point() + geom_line() +
  theme_bw() + labs(x="Block count", y = "Reaction time (ms)", color = "Target position")

# Distractor - target coincidence effect
data_e1 %>% mutate(dtcoinc = ifelse(lag(distractor)=="present", ifelse(target_location == lag(distractor_location), 
                                    "Coincident", "Non-Coincident"),"Dist. absent")) %>% 
  filter(!is.na(dtcoinc), correct==1, Session==1, distractor_type=="orient") %>% 
  group_by(dtcoinc, block_count, distractor) %>% summarize(mRT=mean(response_time)) %>% 
  ggplot(aes(x=block_count, color=dtcoinc, group = dtcoinc, y=mRT)) + 
  geom_point() + geom_line() + theme_bw() +
  labs(x="Block count", y = "Reaction time (ms)", color = "DT coincidence") + 
  facet_wrap(~distractor)

# Target - distractor coincidence effect (not sure there is one?)
data_e1 %>% mutate(tdcoinc = ifelse(distractor_location == lag(target_location), 
                                    "Coincident", "Non-Coincident")) %>% 
  filter(!is.na(tdcoinc), correct==1, Session==1, 
         distractor_type=="orient", distractor=="present") %>% 
  group_by(tdcoinc, block_count) %>% summarize(mRT=mean(response_time)) %>% 
  ggplot(aes(x=block_count, color=tdcoinc, group = tdcoinc, y=mRT)) + 
  geom_point() + geom_line() + theme_bw() + 
  labs(x="Block count", y = "Reaction time (ms)", color = "TD coincidence")
