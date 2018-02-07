
# ---- Libraries ----
library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(devtools) 



# Load (and download if needed) my libraries
lib.GI <- try(library(GI))
if(class(lib.GI)=='try-error'){
    install_github("davidchampredon/GI", 
                   build_vignettes = FALSE, force=TRUE)
    library(GI)
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
popSize <- 2e3  # warning: above 5e3 takes long time!

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

target.val <- c(R0, latent_mean+infectious_mean/2)

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
b  <- sim$GI_bck
at <- sim$acq_times

# ---- Sampled imperfectly observed GIs ----

# We assume that :
# - not all GIs are observed
# - there is an observation error

# Sample the GIs observed:
prop.observed <- 0.99  # proportion of GIs observed
n.obs <- min(length(b), round(prop.observed*popSize) ) # number of bckwd GIs observed
idx.obs <- sample(x = 1:length(b), size = n.obs, replace = FALSE)
gi.obs.true <- b[idx.obs]
at.obs <- at[idx.obs]

# Add observation error:
sd.err <- 0.01
gi.obs <- rnorm(n = n.obs, mean = gi.obs.true, sd = sd.err)
gi.obs[gi.obs<1] <- 1
gi.obs <- round(gi.obs)

# Visualize 'true' vs observed:
plot(gi.obs.true, gi.obs, las=1,
     main = 'observation error for backward GI')
grid()
abline(a = 0,b=1, lty=2)

# plot observed GIs as function of infectee's acquisition time:
df.b <- data.frame(at.obs, gi.obs) %>%
    mutate(rt = round(at.obs))

df.b %>%
    group_by(rt) %>%
    summarise(gi.obs_mean = mean(gi.obs)) %>%
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


# ---- Fit ReSuDe model from GIs ----

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
                       dt = 1.0)


R0.rng     <- seq(1, 6, by=1)
gimean.rng <- seq(2, 8, by=1)
CI <- 0.95
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
                        do.plot = do.plot)
}


if(FALSE){ # Takes too much time... 
    fit.seminr <- gi_ct_fit(t.obs = at.obs,
                            gi.obs = gi.obs,
                            model.epi = 'seminr',
                            fxd.prm = fxd.prm.seminr,
                            R0.rng = R0.rng,
                            gimean.rng = gimean.rng,
                            CI = CI,
                            do.plot = do.plot)
}


if(0){
    library(bbmle)
    fr2 <- gi_ct_fit_mle2(t.obs = at.obs, 
                          gi.obs = gi.obs, 
                          model.epi = 'resude', 
                          fxd.prm = fxd.prm.resude,
                          start.optim = c(R0=2,  gimean=5), 
                          CI = CI,
                          do.plot = FALSE)
    
}