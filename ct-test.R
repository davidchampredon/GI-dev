library(GI, lib.loc = './lib')
library(plyr)

# Calendar times when the GI were observed:
t.obs <- c(6,7,10,21,25,40,40,45,47,47,47,47,55,59,59,59)
# GI at the associated calendar times:
gi.obs <- c(5,6,8,9,10,7,9,9,11,9,10,11,11,10,11,12)

err.fct <- function(x, t.obs, gi.obs) {
    R0 <- x['R0']
    m  <- x['gimean']
    tu <- unique(t.obs)
    
    G <- GI.resude(cal.times.fwdbck = tu,
                   R0,
                   alpha = 0,
                   kappa = 0,
                   GI_span = 20,
                   GI_mean = m,
                   GI_var = NA,
                   GI_type = 'pois',
                   horizon = 200)
    
    gbckmean <- G$bck.mean 
    
    dat.obs <- ddply(data.frame(t.obs,gi.obs),
                     "t.obs", summarize, m=mean(gi.obs),n=length(gi.obs))
    # error is weighted by the number of observations:
    err <- sqrt(sum( dat.obs$n * (gbckmean - dat.obs$m)^2 )/sum(dat.obs$n))
}

fit <- optim(par = c(R0=2, gimean=9),  # starting values
             fn = err.fct,
             t.obs=t.obs, gi.obs=gi.obs)

# Fitted parameters:
fitted.prm <- fit$par
print(fitted.prm)

