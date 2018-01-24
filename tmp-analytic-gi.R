# source('GI.R')
# source('SEmInR_det.R')

library(GI)

x <- GI.seminr(latent_mean = 4,
          infectious_mean = 6,
          R0 = 2.8,
          nE = 5, 
          cal.times.fwdbck = c(1,5,15),
          nI = 5,
          horizon = 50,
          dt = 0.1,
          I.init = 1E-5
)

g <- x$intrinsic
plot(g$tsi, g$density, typ='l', lwd=3)





# library(expm)
# nE <- 4
# nI <- 4
# sig <- 0.1
# gam <- 0.1
# A <- matrix(nrow = nE+nI, ncol=nE+nI)
# 
# for(i in 1:nE) A[i,i] = -sig
# for(i in 2:(nE+1)) A[i,i-1] = sig
# 
# for(i in (nE+1):(nE+nI)) A[i,i] = -gam
# for(i in (nE+2):(nE+nI)) A[i,i-1] = gam
# 
# A[is.na(A)] <- 0
# M <- expm(A) 
# View(M)
