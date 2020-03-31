# Transforms the angular positions used in the experiments to numbers from 1 to 10 for use in the updatePrior function
replacePos <- function(pos) {
  
  angular_positions <- c(seq(30,150,30), seq(210,330,30)) # The positions in degrees on circle used in exp
  for(i in 1:10) { 
    p <- angular_positions[i]
    pos <- replace(pos, which(pos==p), i)
  }
  return(pos)
}

# Updates a Dirichlet distribution based on the list of distactor positions on each trial (pos). Uses
# forgetting similar to the Dynamic Belief Model of Yu & Cohen, but applying the forgetting to distribution parameters 
# rather than the distribution itself. Returns a matrix with the parameters of the 
# Dirichlet distribution on each trial.
updatePriorA <- function(pos, a0=1, m=0.5) { 
  
  pos <- replacePos(pos)
  a <- rep(a0,10)
  mu <- matrix(rep(0,10*length(pos)),length(pos),10)
  mu[1,] <- a/sum(a)
  for(i in 2:length(pos)) {
    prev_pos <- pos[i-1]
    if (prev_pos > 0) {  # distractor present 
     a[prev_pos] <-  a[prev_pos]+1
    }
    a <- a*m+(1-m)*rep(a0,10)
    mu[i,] <- a/sum(a)
  }
  return(mu)
}

# Predicts an RT on each trial based on learning the distractor distributions and assuming that target RTs are slowed
# when appearing in a probable distractor location
predictRT <- function(sub, param) {
  a0 <- param[1]
  mem <- param[2]
  k <- param[3]
  t_a <- param[4]
  k_d <- param[5]
  d_a <- param[6]
  
  if(mem<0) {
    return(10^10)
    }

  if(a0<0) {
    return(10^10)
    }
  
  mu <- updatePriorA(sub$distractor_location, a0=a0, m=mem)
  tpos <- replacePos(sub$target_location)
  tprob <- mu[cbind(seq_along(tpos), tpos)]
  dpos <- replacePos(sub$distractor_location)
  dprob <- rep(0, length(dpos))
  dpres <- which(dpos>0)
  dp_mu <- mu[dpres,]
  dprob[dpres] <- dp_mu[cbind(seq_along(dpos[dpres]), dpos[dpres])]
  pred_rt <-  k + t_a*tprob + as.numeric(dpos>0) * k_d - d_a*dprob
  obs_rt <- sub$response_time
  correct <- which(sub$correct==1)
  obs_rt <- obs_rt[correct]
  pred_rt <- pred_rt[correct]
  
  return(sum((pred_rt - obs_rt)^2) + as.numeric(a0<0)*10^10 + as.numeric(mem<0)*10^10 + as.numeric(mem>1)*10^10)
}

addPredictRT <- function(sub, param) {
  a0 <- param[1]
  mem <- param[2]
  k <- param[3]
  t_a <- param[4]
  k_d <- param[5]
  d_a <- param[6]
  
  mu <- updatePriorA(sub$distractor_location, a0=a0, m=mem)
  tpos <- replacePos(sub$target_location)
  sub$tprob <- mu[cbind(seq_along(tpos), tpos)]
  dpos <- replacePos(sub$distractor_location)
  dprob <- rep(0, length(dpos))
  dpres <- which(dpos>0)
  dp_mu <- mu[dpres,]
  dprob[dpres] <- dp_mu[cbind(seq_along(dpos[dpres]), dpos[dpres])]
  sub$dprob <- dprob
  sub <- mutate(sub, pred_rt =  k + t_a*tprob + as.numeric(distractor_location>0) * k_d - d_a*dprob)

  return(sub)
}
