
# ---- Libraries ----
library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(devtools) 
library(profvis)


# Load (and download if needed) my libraries
use.local.GI <- TRUE

if(use.local.GI){
    library(GI, lib.loc = './lib')
}

if(!use.local.GI){
    lib.GI <- try(library(GI))
    if(class(lib.GI)=='try-error'){
        install_github("davidchampredon/GI", 
                       build_vignettes = FALSE, force=TRUE)
        library(GI)
    }
}

lib.seminribm <- try(library(seminribm))
if(class(lib.seminribm)=='try-error'){
    install_github("davidchampredon/seminribm", 
                   build_vignettes = FALSE, force=TRUE)
    library(seminribm)
}

set.seed(1234)


# ---- Generate data from an individual-based model ----

horizon <- 100
popSize <- 2e4  # warning: above 5e3 takes long time if nE and nI large!

initInfectious  <- 2
R0              <- 3.0
latent_mean     <- 2
infectious_mean <- 4
nE              <- 1
nI              <- 1
calc_WIW_Re     <- FALSE
doExact         <- FALSE
timeStepTauLeap <- 0.1
rnd_seed        <- 1234

gi.mean.true <- latent_mean + (nI+1)/2/nI * infectious_mean

target.val <- c(R0, gi.mean.true)

# See: ?seminribm_run
sim <- seminribm_run(horizon,
                     popSize ,
                     R0 ,
                     latent_mean ,
                     infectious_mean,
                     nE ,
                     nI ,
                     initInfectious ,
                     doExact ,
                     timeStepTauLeap,
                     rnd_seed ,
                     calc_WIW_Re)

# Retrieve backward generation intervals from simulation:
gi.true  <- sim$GI_bck
at <- sim$acq_times

df <- data.frame(at=at, gi.true=gi.true, rt = round(at))
df2 <- df %>%
    group_by(rt) %>%
    summarise(bb = mean(gi.true))



# ---- Sampled imperfectly observed GIs ----

# We assume that :
# - not all GIs are observed
# - there is an observation error

# Sample the GIs observed:
prop.observed <- 0.999  # proportion of GIs observed
n.obs <- min(length(gi.true), round(prop.observed*popSize) ) # number of bckwd GIs observed
idx.obs <- sample(x = 1:length(gi.true), size = n.obs, replace = FALSE)
gi.obs.true <- gi.true[idx.obs]
at.obs <- at[idx.obs]

# Add observation error:
sd.err <- 0.001
gi.obs <- rnorm(n = n.obs, mean = gi.obs.true, sd = sd.err)
gi.obs[gi.obs<1] <- 1
gi.obs <- round(gi.obs)


df.gi.obs <- data.frame(t = at.obs, 
                        gi.obs.true = gi.obs.true, 
                        gi.obs = gi.obs)
df.gi.obs <- df.gi.obs[order(df.gi.obs$t),]

# Visualize 'true' vs observed:
plot(gi.obs.true, gi.obs, las=1,
     main = 'observation error for backward GI')
grid()
abline(a = 0,b=1, lty=2)

# plot observed GIs as function of infectee's acquisition time:
df.b <- data.frame(at.obs, gi.obs) %>%
    mutate(rt = round(at.obs))

df.b2 <- df.b %>%
    group_by(rt) %>%
    summarise(gi.obs_mean = mean(gi.obs)) 

df.b2%>%
    ggplot(aes(x=rt,y=gi.obs_mean)) + 
    geom_point(data = df.b, 
               aes(x=at.obs,y=gi.obs), 
               alpha=0.15, colour="orange",
               pch=16,size=4) +
    geom_abline(slope = 1, intercept = 0, linetype=2, colour = "grey")+
    geom_line(size=1.5) + 
    geom_point(size=2) +
    ggtitle('Backward Generation Intervals (line: daily mean)')+
    xlab('calendar time')+ylab('days')


# ---- Fit model from GIs ----

fxd.prm.resude <- list(horizon=horizon, 
                       alpha=0, 
                       kappa=0, 
                       GI_span = 20, 
                       GI_var = 5, 
                       GI_type = 'pois', 
                       dt = 1.0)


fxd.prm.seminr <- list(horizon=horizon, 
                       nE=nE, 
                       nI=nI, 
                       latent_mean=latent_mean, 
                       dt = 0.5)


R0.rng     <- seq(1.5, 6, by=0.25)
gimean.rng <- seq(2, 10, by=0.5)
CI <- 0.90
do.plot <- TRUE

# See: ?gi_ct_fit
if(FALSE){
    fit.resude <- gi_ct_fit(t.obs = at.obs, 
                            gi.obs = gi.obs, 
                            model.epi = 'resude', 
                            fxd.prm = fxd.prm.resude,
                            R0.rng = R0.rng, 
                            gimean.rng = gimean.rng,
                            CI = CI,
                            do.plot = do.plot,
                            R0.true = R0,
                            gimean.true = gi.mean.true)
}


if(TRUE){ 
    
    # STOPPED HERE
    # Tidy everything.
    
    # Fit with SEmInR only, as the data were generated with this model. 
    # (although there is a link with ReSuDe)
    
    # Make a function of all this, with the goal of increasing 
    # 'fit.first.n.obs' to show how the fit is better as it increases.
    # 
    # Also, try to speed up further... (although running on HPC should be quick).
    
    cal.t.bck <- 1:horizon
    
    z.true.det <- GI.seminr(latent_mean = latent_mean, 
                            infectious_mean = infectious_mean,
                            R0 = R0, 
                            nE = nE, 
                            nI = nI,
                            cal.times.fwdbck = cal.t.bck, 
                            horizon = horizon, 
                            calc.fwd = FALSE)
    
    # Fit on only the n first observations:
    fit.first.n.obs <- 15  # (30: 8min on 4 cpus)
    gi.1st.obs <- df.gi.obs[df.gi.obs$t< fit.first.n.obs,]
    gi.1st.obs$t.obs <- round(gi.1st.obs$t)
    
    # Fitting to observed bckwd GIs:
    fit.seminr <- gi_ct_fit(t.obs  = gi.1st.obs$t.obs, # at.obs,cal.t.bck, #
                            gi.obs = gi.1st.obs$gi.obs, # gi.obs,round(z.true.det$bck.mean), #
                            model.epi = 'seminr',
                            fxd.prm = fxd.prm.seminr,
                            R0.rng = R0.rng,
                            gimean.rng = gimean.rng,
                            CI = CI,
                            R0.true = R0,
                            gimean.true = gi.mean.true,
                            do.plot = TRUE)
    
    z.fit <- GI.seminr(latent_mean = latent_mean, 
                       infectious_mean = infectious_mean,
                       R0 = fit.seminr$R0.best, 
                       nE = nE, 
                       nI = nI,
                       cal.times.fwdbck = cal.t.bck, 
                       horizon = horizon, 
                       calc.fwd = FALSE)
    
    plot(x=df$at, 
         y=df$gi.true, 
         col=rgb(0,0.3,0,0.1), 
         las=1,
         xlab='calendar time',
         ylab='Backward GI',
         log='y')
    lines(df2$rt, df2$bb, lwd=3)
    lines(cal.t.bck, z.true.det$bck.mean, col='black', lwd=2, lty=2)
    lines(cal.t.bck, z.fit$bck.mean, col=rgb(0,0,1,0.5), lwd=6, lty=1)
    abline(v=fit.first.n.obs, lty=4)
}


if(FALSE){ # TESTING.... 
    
    profvis(expr = 
    {
        z <- nllk(R0 = 3, 
                  gimean = 3,
                  t.obs = cal.t.bck, #at.obs,
                  gi.obs = round(z.true.det$bck.mean), #gi.obs,
                  model.epi = 'seminr', 
                  fxd.prm = fxd.prm.seminr)
    })
}


if(0){
    library(bbmle)
    
    fr2 <- gi_ct_fit_mle2(t.obs  = gi.1st.obs$t.obs, #at.obs, 
                          gi.obs = gi.1st.obs$gi.obs, #gi.obs, 
                          model.epi = 'seminr', 
                          fxd.prm = fxd.prm.seminr,
                          start.optim = c(R0=2,  gimean=5), 
                          CI = CI,
                          do.plot = FALSE)
    
    
    fr2 <- gi_ct_fit_mle2(t.obs = at.obs, 
                          gi.obs = gi.obs, 
                          model.epi = 'resude', 
                          fxd.prm = fxd.prm.resude,
                          start.optim = c(R0=2,  gimean=5), 
                          CI = CI,
                          do.plot = FALSE)
    
}