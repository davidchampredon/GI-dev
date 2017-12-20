library(GI, lib.loc = './lib')
library(dplyr)

set.seed(12345)

# Simulate 'true' data:
R0.true <- 3
dt <- 0.25
hz <- 100
nE <- nI <- 4
latent_mean <- 2
infectious_mean <- 4

# Calendar times when the GI were observed:
t.obs <- sort(sample(x = 1:(hz/2), 
                     size = 35, 
                     replace = TRUE))

G.true <- GI.seminr(latent_mean = latent_mean,
                    infectious_mean = infectious_mean, 
                    R0 = R0.true, 
                    nE = nE, 
                    nI = nI,
                    cal.times.fwdbck = t.obs,
                    horizon = hz, 
                    dt = dt)

gibck.true <-  G.true$bck.mean

# Imperfectly observed GI at the associated calendar times:
gi.obs <- rpois(n = length(gibck.true),lambda = gibck.true)
gi.obs[gi.obs==0] <- 1

par(mfrow=c(1,2))
plot(x=1:length(G.true$incidence) * dt, 
     y=G.true$incidence, 
     typ='l', xlab='day',ylab='incidence') ; grid()
plot(t.obs, gibck.true, typ='b',ylim=range(c(0,gibck.true,gi.obs)))
points(x=t.obs, y=gi.obs, pch=16)

# Model fixed parameters (not fitted):
fxd.prm.seminr <- list(horizon=hz, nE=nE, nI=nI, latent_mean=latent_mean, dt=dt)
fxd.prm.resude <- list(horizon=hz, alpha=0, kappa=0, GI_span=20, 
                       GI_var=NULL, GI_type='pois', dt=dt)

# Fit (one for each model):

R0.rng     <- seq(1, 6, by=1) #0.25)
gimean.rng <- seq(2, 8, by=1) #0.5)
CI <- 0.95
do.plot <- TRUE

if(0){ # Model SEmInR takes a long time
gi_ct_fit(t.obs = t.obs,
          gi.obs = gi.obs,
          model.epi = 'seminr',
          fxd.prm = fxd.prm.seminr,
          R0.rng = R0.rng,
          gimean.rng = gimean.rng,
          CI = CI,
          do.plot = do.plot)
}

gi_ct_fit(t.obs = t.obs, 
          gi.obs = gi.obs, 
          model.epi = 'resude', 
          fxd.prm = fxd.prm.resude,
          R0.rng = R0.rng, 
          gimean.rng = gimean.rng,
          CI = CI,
          do.plot = do.plot)


gi_ct_fit_mle2(t.obs = t.obs, 
          gi.obs = gi.obs, 
          model.epi = 'resude', 
          fxd.prm = fxd.prm.resude,
          start.optim = c(R0=2, gimean=2),
          CI = CI,
          do.plot = do.plot)


message('|=== CT TEST DONE ===|')
