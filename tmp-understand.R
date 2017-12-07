library(GI, lib.loc = './lib')

set.seed(12345)

# Calendar times when the GI were observed:
t.obs <- 1:50 #sort(sample(x = 1:(hz/2), size = 15, replace = TRUE))

G.true <- GI.resude(cal.times.fwdbck = t.obs,
                    R0 = 2.5,
                    alpha = 0, 
                    kappa = 0, 
                    GI_span = 20, 
                    GI_mean = 3, 
                    GI_var = NULL, 
                    GI_type = 'pois',
                    horizon = 100)

gibck.true <-  G.true$bck.mean

par(mfrow=c(1,2))
plot(x=1:length(G.true$incidence) ,
     y=G.true$incidence, 
     typ='l', xlab='day',ylab='incidence') ; 
grid()
plot(x = t.obs, 
     y = gibck.true, 
     typ='o')