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
t.obs <- sort(sample(x = 1:(hz/2), size = 15, replace = TRUE))

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


fxd.prm <- list(horizon=hz, nE=nE, nI=nI, latent_mean=latent_mean, dt=dt)

nllk(R0 = 2.3, gimean = 3, t.obs = t.obs, gi.obs = gi.obs, model.epi = 'seminr', fxd.prm = fxd.prm)

gi_ct_fit(t.obs = t.obs, 
          gi.obs = gi.obs, 
          model.epi = 'seminr', 
          fxd.prm = fxd.prm,
          R0.rng = seq(2,6,by=1), 
          gimean.rng = seq(1,6,by=1),
          CI = 0.95)





if(FALSE){
    llk <- function(R0, infectmean, t.obs, gi.obs) {
        
        # Calculate the backward GI for given parameter values
        # at the observed dates:
        G <-  GI.seminr(latent_mean = latent_mean,
                        infectious_mean = infectmean, 
                        R0 = R0, 
                        nE = nE, 
                        nI = nI,
                        cal.times.fwdbck = t.obs,
                        horizon = hz, 
                        dt = dt)
        
        # Extract the mean bckwd GI :
        b <- G$bck.mean    # ; plot(t.obs,b,typ='b')
        
        # Calculate log likelihood:
        z <- sum(dpois(x = gi.obs, lambda = b, log = TRUE))
        return(-z)
    }
    
    
    # a <- llk(R0=3, infectmean=4, t.obs, gi.obs)
    
    R0.rng <- seq(1.0, 6, by=0.5)
    gimean.rng <- seq(1,10,by=1)
    
    system.time(
        M <- outer(X = R0.rng, 
                   Y = gimean.rng, 
                   FUN = Vectorize(llk, list("R0","infectmean")), 
                   t.obs=t.obs, 
                   gi.obs=gi.obs)
    )
    
    idx <- which(M == min(M, na.rm = TRUE), 
                 arr.ind = TRUE)
    
    # Retrieve neg log likelihood min
    ll.min <- M[idx]
    # Set confidence threshold
    confidence.threshold = 0.90
    # Calculate the level on the likelihoof function
    conf.cutoff <- ll.min+ qchisq(confidence.threshold,2)/2
    
    
    contour(x = R0.rng, 
            y = gimean.rng, 
            z = M, 
            nlevels = 30,
            col='grey',
            xlab = 'R0', ylab = 'Mean GI')
    points(x = R0.true, y=infectious_mean, 
           cex=2, pch=3, col='blue', lwd=3)
    points(x = R0.rng[idx[1]], y = gimean.rng[idx[2]], pch=1, cex=3,lwd=3, col='red')
    # points(x = naive.fit['R0'], y = naive.fit['infectmean'], pch=4, cex=2)
    # Draw the confidence contour
    contour(R0.rng,
            gimean.rng,
            M,
            level = conf.cutoff,
            labels="",
            col="red",
            lwd=3, 
            lty=1,
            add = TRUE)
    
    # Retrieve the extreme values of the confidence contour
    cc = contourLines(R0.rng,gimean.rng,M,level = conf.cutoff)#[[1]]
    
    
    
    # Finally, let’s visually check how good the fit is:
    
    tu <- unique(t.obs)
    Gfit <- GI.resude(cal.times.fwdbck = tu,
                      R0 = fitted.prm['R0'],
                      alpha = 0,
                      kappa = 0,
                      GI_span = 20,
                      GI_mean = fitted.prm['gimean'],
                      GI_var = NA,
                      GI_type = 'pois',
                      horizon = 200)
    
    gbck.fit <- Gfit$bck.mean
    
    plot(x = tu, y = gbck.fit, 
         ylim = range(gbck.fit,gi.obs),
         typ='o', pch=16,lwd=3,
         las = 1, xlab = 'calendar time', ylab='days',
         main = 'RESuDe model fitted to contact tracing data')
    points(t.obs, gi.obs, pch=3, col='red',lwd=6,cex=2)
    grid()
    legend('topleft', legend = c('data','model fit'),col=c('red','black'),pch=c(3,16),lwd=c(NA,3),pt.cex = c(1,1), pt.lwd = c(6,1))
    
}