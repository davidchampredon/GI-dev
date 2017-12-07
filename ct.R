### 
###   Fit to tontact tracing data
###

#' Negative Log-Likelihood given R0 and mean generation interval
#' assuming Poisson observation error from contact tracing data
#' 
#' @param t.obs Numeric vector. Time when the backward generation intervals were observed
#' @param gi.obs Numeric vector. Backward generation intervals observed.
#' @param model.epi String. Epidemic model used. Choice between \code{seminr} or \code{resude}.
#' @param fxd.prm List. Parameters of \code{model.epi} that are fixed.
#' 
nllk <- function(R0, 
                 gimean,
                 t.obs, 
                 gi.obs, 
                 model.epi, 
                 fxd.prm) {
    
    if(model.epi=='seminr'){
        
        im <- 2*(gimean - fxd.prm[['latent_mean']])  # gi.mean ~ latent + infectious/2
        if(im <= 0){
            warning('GI and latent period lengths not consistent!')
            im <- 1
        }
        G <- GI.seminr(latent_mean = fxd.prm[['latent_mean']],
                       infectious_mean = im, 
                       R0 = R0, 
                       nE = fxd.prm[['nE']], 
                       nI = fxd.prm[['nI']],
                       cal.times.fwdbck = t.obs,
                       horizon = fxd.prm[['horizon']], 
                       dt = fxd.prm[['dt']])
    }
    if(model.epi=='resude'){
        G <- GI.resude(cal.times.fwdbck = t.obs,
                       R0 = R0,
                       alpha = fxd.prm[['alpha']], 
                       kappa = fxd.prm[['kappa']], 
                       GI_span = fxd.prm[['GI_span']], 
                       GI_mean = gimean, 
                       GI_var = fxd.prm[['GI_var']], 
                       GI_type = fxd.prm[['GI_type']],
                       horizon = fxd.prm[['horizon']])
    }
    # Extract the mean bckwd GI :
    b <- G$bck.mean    
    
    # Calculate log likelihood:
    z <- sum(dpois(x = gi.obs, lambda = b, log = TRUE))
    return(-z)
}







#' Fit R0 and mean generation interval from contact tracing data
#' 
#' @param t.obs Numeric vector. Time when the backward generation intervals were observed
#' @param gi.obs Numeric vector. Backward generation intervals observed.
#' @param model.epi String. Epidemic model used. Choice between \code{seminr} or \code{resude}.
#' @param fxd.prm List. Parameters of \code{model.epi} that are fixed.
#' @export
gi_ct_fit <- function(t.obs, 
                      gi.obs, 
                      model.epi, 
                      fxd.prm,
                      R0.rng, 
                      gimean.rng,
                      CI = 0.95) {
    
    M <- outer(X = R0.rng, 
               Y = gimean.rng, 
               FUN = Vectorize(nllk, list("R0","gimean")), 
               t.obs = t.obs, 
               gi.obs = gi.obs, 
               model.epi = model.epi, 
               fxd.prm = fxd.prm)
    
    idx <- which(M == min(M, na.rm = TRUE), 
                 arr.ind = TRUE)
    
    # Retrieve neg log likelihood min
    ll.min <- M[idx]
    # Set confidence threshold
    CI = 0.90
    # Calculate the level on the likelihoof function
    conf.cutoff <- ll.min + qchisq(CI,2)/2
    
    # Likelihood surface:
    contour(x = R0.rng, 
            y = gimean.rng, 
            z = M, 
            nlevels = 30,
            col='grey',
            xlab = 'R0', ylab = 'Mean GI')
    
    # Best point estimate:
    points(x = R0.rng[idx[1]], 
           y = gimean.rng[idx[2]],
           pch=1, cex=3,lwd=3, col='red')
    
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
    
}