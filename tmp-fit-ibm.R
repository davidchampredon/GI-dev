library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())

library(devtools) 
install_github("davidchampredon/seminribm", build_vignettes = TRUE, force=TRUE)
install_github("davidchampredon/GI", build_vignettes = FALSE, force=TRUE)
library(seminribm)
library(GI)

set.seed(1234)

# ---- Generate data from an individual-based model ----

horizon <- 100
popSize <- 5e3

initInfectious  <- 2
R0              <- 3.0
latent_mean     <- 2
infectious_mean <- 4
nE              <- 6
nI              <- 6
calc_WIW_Re     <- FALSE
doExact         <- FALSE
timeStepTauLeap <- 0.1
rnd_seed        <- 1234

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

# Retrieve backward generation intervals:
b  <- sim$GI_bck
at <- sim$acq_times


# ---- Sampled imperfectly observed GIs ----

# We assume that :
# - not all GIs are observed
# - there is an observation error

prop.observed <- 0.1  # proportion of GIs observed
n.obs <- round(prop.observed*popSize)  # number of bckwd GIs observed
idx.obs <- sample(x = 1:length(b), size = n.obs, replace = FALSE)

gi.obs.true <- b[idx.obs]
at.obs <- at[idx.obs]

# Add observation error:
gi.obs <- rnorm(n = n.obs, mean = gi.obs.true, sd = 0.5)
gi.obs[gi.obs<1] <- 1
gi.obs <- round(gi.obs)
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


# ---- Fit model from GIs ----

fxd.prm.seminr <- list(horizon=horizon, 
                       nE=nE, 
                       nI=nI, 
                       latent_mean=latent_mean, 
                       dt = 1.0)

fxd.prm.resude <- list(horizon=horizon, alpha=0, kappa=0, 
                       GI_span=20, 
                       GI_var=NULL, 
                       GI_type='pois', 
                       dt=1.0)

R0.rng     <- seq(1, 6, by=0.25)
gimean.rng <- seq(2, 8, by=0.5)
CI <- 0.80
do.plot <- TRUE

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

fit.resude <- gi_ct_fit(t.obs = at.obs, 
                        gi.obs = gi.obs, 
                        model.epi = 'resude', 
                        fxd.prm = fxd.prm.resude,
                        R0.rng = R0.rng, 
                        gimean.rng = gimean.rng,
                        CI = CI,
                        do.plot = do.plot)






