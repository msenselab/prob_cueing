library(ez) # ezANOVA-Funktions
library(tidyverse)

source('fun.R')

data_e1 <- read.table(file="./data/data_LT_E1.csv", sep=",", dec=".", head=TRUE)
data_e2 <- read.table(file="./data/data_LT_E2.csv", sep=",", dec=".", head=TRUE)

subs <- unique(filter(data_e1, distractor_type=="orient")$subject_nr)

data_pred <- data.frame() 
fit_param <- data.frame()
# STAN fitting is really slow, so lets do ML (least squares) for now
#smod <- stan_model('model.stan') 
for(s in subs) {
  sub <- filter(data_e1, subject_nr==s, Session==1, distractor_type=="orient")
  
  # STAN fit
#  fit <- sampling(smod, list(N=dim(sub)[1], rt=sub$response_time, dpos=replacePos(sub$distractor_location), 
#                            tpos=replacePos(sub$target_location)))

  # ML fit
  fitRT <- function(p) {
    predictRT(sub, p)
  }

  param <- optim(c(10, 0.85, 630, 430, 220, 840), fitRT, method="SANN", control= list(maxit = 30000))
  print(param)
  param <- optim(param$par, fitRT, control= list(maxit = 30000))
  print(param)
  fit_param <- rbind(fit_param, data.frame(subject=sub$subject_nr[1], a0 = param$par[1], mem = param$par[2], k = param$par[3],
                                           t_a = param$par[4], k_d = param$par[5], d_a = param$par[6]))
  sub <- addPredictRT(sub, param$par)
  data_pred <- rbind(data_pred, sub)
}