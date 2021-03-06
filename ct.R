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
#' @export
nllk <- function(R0, 
                 gimean,
                 t.obs, 
                 gi.obs, 
                 model.epi, 
                 fxd.prm) {
    
    message(paste0('Evaluating likelihood R0=',R0, ' gimean=',gimean,' ...'),
            appendLF = FALSE)
    
    if(model.epi=='seminr'){
        nI = fxd.prm[['nI']]
        # gi.mean ~ latent + infectious*(n+1)/(2n)
        im <- 2*nI / (nI+1) *(gimean - fxd.prm[['latent_mean']])  
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
                       dt = fxd.prm[['dt']],
                       calc.fwd = FALSE)
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
    b[b==0] <- 1e-6
    
    # Calculate log likelihood:
    # tmp <- dpois(x = gi.obs, lambda = b, log = TRUE)
    # gi.obs[which(is.na(tmp))]
    z <- -sum(dpois(x = gi.obs, lambda = b, log = TRUE), na.rm = TRUE)
    message(paste0(' nllk = ', z))
    return(z)
}

# Plot neg log likelihood surface
plot_nllk_surf <- function(model.epi, 
                           CI,
                           conf.cutoff,
                           M,
                           R0.rng, gimean.rng,
                           R0.best, gimean.best,
                           R0.ci, gimean.ci) {
    
    title <- paste0('Neg. Log-Likelihood Surface (',model.epi,')\n with ',CI*100,'%CI')
    ylab <- 'Mean intrinsic GI'
    if(model.epi=='seminr') ylab <- 'Mean infectious period'
    
    contour(x = R0.rng, 
            y = gimean.rng, 
            z = M, 
            nlevels = 20,
            col='grey',
            main = title,
            xlab = 'R0', 
            ylab = ylab,
            las = 1)
    
    # Best point estimate:
    col.best <- 'red'
    
    points(x = R0.best, 
           y = gimean.best,
           pch=16, cex=1, col=col.best)
    segments(x0 = R0.best, 
             y0 = 0,
             x1 = R0.best,
             y1 = gimean.best,
             col = col.best, lty=2)
    segments(x0 = 0,
             y0 = gimean.best,
             x1 = R0.best,
             y1 = gimean.best,
             col = col.best, lty=2)
    
    arrows(x0=R0.ci[1], x1=R0.ci[2],
           y0 = gimean.best, y1 = gimean.best, 
           angle = 90, length = 0.1,code = 3, 
           col = col.best,
           lwd=1.3)
    
    arrows(y0=gimean.ci[1], y1=gimean.ci[2],
           x0 = R0.best, x1 = R0.best, 
           angle = 90, length = 0.1, code=3,
           col = col.best,
           lwd=1.3)
    
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


#' Fit R0 and mean generation interval from contact tracing data
#' 
#' @param t.obs Numeric vector. Time when the backward generation intervals were observed
#' @param gi.obs Numeric vector. Backward generation intervals observed.
#' @param model.epi String. Epidemic model used. Choice between \code{seminr} or \code{resude}.
#' @param fxd.prm List. Parameters of \code{model.epi} that are fixed.
#' @param R0.rng Numeric vector. Values of R0 explored to calculate the likelihood surface.
#' @param gimean.rng Numeric vector. Values of mean intrinsic GI explored to calculate the likelihood surface.
#' @param CI Numeric. Confidence interval level. Default = 0.95.
#' @param do.plot Boolean. Fitting diagnostic plots. Default = TRUE.
#' @param R0.true Numeric. Value of the 'true' (simulated) R0
#' @param gimean.true Numeric. Value of the 'true' (simulated) mean GI
#' @export
gi_ct_fit <- function(t.obs, 
                      gi.obs, 
                      model.epi, 
                      fxd.prm,
                      R0.rng, 
                      gimean.rng,
                      CI = 0.95,
                      do.plot = FALSE,
                      R0.true = NULL,
                      gimean.true = NULL) {
    
    t1 <- as.numeric(Sys.time())
    
    # M <- outer(X = R0.rng, 
    #            Y = gimean.rng, 
    #            FUN = Vectorize(nllk, list("R0","gimean")), 
    #            t.obs = t.obs, 
    #            gi.obs = gi.obs, 
    #            model.epi = model.epi, 
    #            fxd.prm = fxd.prm)
    
    ncpus <- parallel::detectCores()
    
    sfInit(parallel = ncpus>1, cpus = ncpus)
    
    nllk_R0 <- function(R0,gimean,t.obs, gi.obs, model.epi, fxd.prm) {
        a <- nllk(R0,gimean,t.obs, gi.obs, model.epi, fxd.prm)
    }
    sfLibrary(GI, lib.loc = './lib')
    sfExportAll()
    x <- list()
    for(i in 1:length(gimean.rng)){
        x[[i]] <- sfSapply(R0.rng,
                 fun = nllk_R0,
                 gimean = gimean.rng[i],
                 t.obs, gi.obs, 
                 model.epi, fxd.prm)
    }
    sfStop()
    
    M <- do.call(cbind, x)
    
    idx <- which(M == min(M, na.rm = TRUE), 
                 arr.ind = TRUE)
    
    # Retrieve neg log likelihood min
    ll.min <- M[idx]
    R0.best     <- R0.rng[idx[1]]
    gimean.best <- gimean.rng[idx[2]]
    
    # Calculate the level on the likelihoof function
    conf.cutoff <- ll.min + qchisq(CI,2)/2
    
    # Estimate the confidence interval from the contour lines
    cc <- contourLines(R0.rng,
                       gimean.rng,
                       M,
                       level = conf.cutoff)
    x.tmp <- list()
    y.tmp <- list()
    
    for(i in seq_along(cc)){
        x.tmp[[i]] <- range(cc[[i]]$x)
        y.tmp[[i]] <- range(cc[[i]]$y)
    }
    R0.ci     <- range(unlist(x.tmp))
    gimean.ci <- range(unlist(y.tmp))
    
    # Likelihood surface:
    if(do.plot){
        
        plot_nllk_surf(model.epi,
                       CI,
                       conf.cutoff,
                       M,
                       R0.rng,
                       gimean.rng,
                       R0.best,
                       gimean.best,
                       R0.ci,
                       gimean.ci)
        
        if(!is.null(R0.true) & !is.null(gimean.true)){
            points(x= R0.true, y = gimean.true, col='blue', cex=3, pch=9)
            text(x= R0.true, y = gimean.true, col='blue', labels = 'true', pos = 3)
        }
        
        # Reorder data, because "GI.resude" expects so:
        idx2 <- order(t.obs)
        tt <- t.obs[idx2]
        
        col.best <- 'red'
        if(model.epi=='seminr'){
            im    <- 2*nI / (nI+1) *(gimean.best - fxd.prm[['latent_mean']])  # gi.mean ~ latent + infectious * (nI+1)/2/nI
            im.lo <- 2*nI / (nI+1) *(gimean.ci[1] - fxd.prm[['latent_mean']])  
            im.hi <- 2*nI / (nI+1) *(gimean.ci[2] - fxd.prm[['latent_mean']])  
            
            if(im <= 0) im <- 1
            if(im.lo <= 0) im.lo <- 1
            if(im.hi <= 0) im.hi <- 1
            
            Gfit <- GI.seminr(latent_mean = fxd.prm[['latent_mean']],
                              infectious_mean = im, 
                              R0 = R0.best, 
                              nE = fxd.prm[['nE']], 
                              nI = fxd.prm[['nI']],
                              cal.times.fwdbck = tt,
                              horizon = fxd.prm[['horizon']], 
                              dt = fxd.prm[['dt']],
                              calc.fwd = FALSE)
            
            Gfit.lo <- GI.seminr(latent_mean = fxd.prm[['latent_mean']],
                                 infectious_mean = im.lo, 
                                 R0 = R0.ci[1], 
                                 nE = fxd.prm[['nE']], 
                                 nI = fxd.prm[['nI']],
                                 cal.times.fwdbck = tt,
                                 horizon = fxd.prm[['horizon']], 
                                 dt = fxd.prm[['dt']],
                                 calc.fwd = FALSE)
            Gfit.hi <- GI.seminr(latent_mean = fxd.prm[['latent_mean']],
                                 infectious_mean = im.hi, 
                                 R0 = R0.ci[2], 
                                 nE = fxd.prm[['nE']], 
                                 nI = fxd.prm[['nI']],
                                 cal.times.fwdbck = tt,
                                 horizon = fxd.prm[['horizon']], 
                                 dt = fxd.prm[['dt']],
                                 calc.fwd = FALSE)
            
        }
        if(model.epi=='resude'){
            Gfit <- GI.resude(cal.times.fwdbck = tt, 
                              R0 = R0.best,
                              alpha = fxd.prm[['alpha']], 
                              kappa = fxd.prm[['kappa']], 
                              GI_span = fxd.prm[['GI_span']], 
                              GI_mean = gimean.best, 
                              GI_var = fxd.prm[['GI_var']], 
                              GI_type = fxd.prm[['GI_type']],
                              horizon = fxd.prm[['horizon']])
            Gfit.lo <- GI.resude(cal.times.fwdbck = tt, 
                                 R0 = R0.ci[1],
                                 alpha = fxd.prm[['alpha']], 
                                 kappa = fxd.prm[['kappa']], 
                                 GI_span = fxd.prm[['GI_span']], 
                                 GI_mean = gimean.ci[1], 
                                 GI_var = fxd.prm[['GI_var']], 
                                 GI_type = fxd.prm[['GI_type']],
                                 horizon = fxd.prm[['horizon']])
            Gfit.hi <- GI.resude(cal.times.fwdbck = tt, 
                                 R0 = R0.ci[2],
                                 alpha = fxd.prm[['alpha']], 
                                 kappa = fxd.prm[['kappa']], 
                                 GI_span = fxd.prm[['GI_span']], 
                                 GI_mean = gimean.ci[2], 
                                 GI_var = fxd.prm[['GI_var']], 
                                 GI_type = fxd.prm[['GI_type']],
                                 horizon = fxd.prm[['horizon']])
        }
        
        gbck.fit <- Gfit$bck.mean
        gbck.fit.lo <- Gfit.lo$bck.mean
        gbck.fit.hi <- Gfit.hi$bck.mean
        
        plot(x = tt, 
             y = gbck.fit, 
             ylim = range(gbck.fit,gi.obs, na.rm = TRUE),
             typ='o', pch=16, lwd=3,
             col = col.best,
             las = 1, 
             xlab = 'calendar time',
             ylab='Mean Backward GI',
             main = paste(model.epi, 'model backward GI\nfitted to contact tracing data'))
        
        
        lines(x=tt, y=gbck.fit.lo, lty=2, col=col.best, lwd=2)
        lines(x=tt, y=gbck.fit.hi, lty=2, col=col.best, lwd=2)
        points(t.obs, gi.obs, pch=1, col=rgb(0,0,0,0.7),lwd=2,cex=1)
        grid()
        legend('topleft', legend = c('data','model fit',paste(CI*100,'%CI')),
               col=c('black',col.best,col.best),pch=c(1,16,NA),lwd=c(NA,3,2),
               pt.cex = c(1,1,NA), pt.lwd = c(2,1,NA), lty=c(1,1,2))
    }
    t2 <- as.numeric(Sys.time())
    dt <- round( (t2-t1)/60, 1)
    msg <- paste('Fit to contact tracing data done in',dt,'minute(s).')
    message(msg)
    return(list(R0.best = R0.best,
                gimean.best = gimean.best,
                R0.ci = R0.ci,
                gimean.ci = gimean.ci))
    
}

#' Fit R0 and mean generation interval from contact tracing data
#' 
#' @param t.obs Numeric vector. Time when the backward generation intervals were observed
#' @param gi.obs Numeric vector. Backward generation intervals observed.
#' @param model.epi String. Epidemic model used. Choice between \code{seminr} or \code{resude}.
#' @param fxd.prm List. Parameters of \code{model.epi} that are fixed.
#' @param R0.rng Numeric vector. Values of R0 explored to calculate the likelihood surface.
#' @param gimean.rng Numeric vector. Values of mean intrinsic GI explored to calculate the likelihood surface.
#' @param CI Numeric. Confidence interval level. Default = 0.95.
#' @param do.plot Boolean. Fitting diagnostic plots. Default = TRUE.
#' @importFrom bbmle mle2 profile confint plot
#' @export
gi_ct_fit_mle2 <- function(t.obs, 
                      gi.obs, 
                      model.epi, 
                      fxd.prm,
                      start.optim,
                      CI = 0.95,
                      do.plot = FALSE ) {
    
    t1 <- as.numeric(Sys.time())
    
    a <- c(list(fxd.prm=fxd.prm), 
           list(model.epi=model.epi),
           list(t.obs = t.obs),
           list(gi.obs = gi.obs))
    
    fit.mle2 <- mle2(minuslogl = nllk, 
                     start = list(R0=start.optim['R0'], 
                                  gimean=start.optim['gimean']), 
                     data = a,
                     control = list(ndeps=c(0.01,0.02)))
    
    pmle2 <- profile(fit.mle2)
    mle.ci <- confint(pmle2, level = CI)
    if(do.plot) plot(pmle2)
    
    t2 <- as.numeric(Sys.time())
    dt <- round( (t2-t1)/60, 1)
    msg <- paste('Fit to contact tracing data done in',dt,'minute(s).')
    message(msg)
    return(list(R0.best = fit.mle2@coef['R0'],
                gimean.best = fit.mle2@coef['gimean'],
                R0.ci = mle.ci['R0',],
                gimean.ci = mle.ci['gimean',]))
}





