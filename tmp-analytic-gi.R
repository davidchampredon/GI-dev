library(GI)
library(expm)

# ----- Build A ----
m <- 6
n <- 4

DOL <- 3
DOI <- 5
sig <- 1 / DOL
gam <- 1 / DOI
A <- matrix(nrow = n+m, ncol=n+m)
D <- matrix(nrow = n+m, ncol=n+m)
N <- matrix(nrow = n+m, ncol=n+m)

for(i in 1:m) D[i,i] = -sig
for(i in 2:(m+1)) N[i,i-1] = sig
for(i in (m+1):(m+n)) D[i,i] = -gam
for(i in (m+2):(m+n)) N[i,i-1] = gam

D[is.na(D)] <- 0
N[is.na(N)] <- 0
A <- D + N
all(D %*% N == N %*% D)


# ----- Test power N ----

Nk_i <- function(sig,gam,m,n,i,k) {
    
    if(k==0) res <- 1
    if(k>0){
        res <- res1 <- res2 <- res3 <- 0
        
        if(k+1 <= i & i <= m+1)   res1 <-  sig^k
        if(2 <= i-m & i-m <=k)    res2 <- gam^(i-m-1) * sig^(k-i+m+1)
        if(m+k+1 <= i & i <= m+n) res3 <- gam^k
        
        # to delete:
        # res2a <- sig^p * ifelse(p+1 <= i & i <= m+1, 1, 0)
        # res2b <- gam^(i-m-1) * sig^(p-i+m+1)
        # res2c <- gam^p * ifelse(m+1+p <= i & i <= m+n, 1, 0)
        # res2 <- res2a+res2b+res2c
        
        res <- res1+res2+res3
    }
    return(res)
}

Nk_i(sig,gam,m,n,i=1,k=0)

Nk_formula <- function(sig,gam,m,n,k){
    Nk_f <- matrix(nrow = n+m, ncol=n+m)
    
    for(i in 1:(m+n)){
        if(i-k>0) Nk_f[i,i-k] <- Nk_i(sig,gam,m,n,i,k)
    }    
    Nk_f[is.na(Nk_f)] <- 0
    return(Nk_f)
}

k <- 5
Nk <- N %^% k
Nk_f <- Nk_formula(sig,gam,m,n,k)
sum( (Nk-Nk_f)^2 )


# ----- Check H:=exp(D) -----

H <- expm(D)

H_formula <- function(sig,gam,m,n) {
    Hf <- matrix(nrow = n+m, ncol=n+m)
    for(i in 1:(m+n)) {
        if(i <= m)   Hf[i,i] <- exp(-sig)
        if(i >= m+1) Hf[i,i] <- exp(-gam)
    }
    Hf[is.na(Hf)] <- 0
    return(Hf)
}

Hf <- H_formula(sig,gam,m,n)
sum( (H-Hf)^2)


# ----- Check G:=exp(N) -----

G_ij <- function(sig,gam,m,n,i,j) {
    
    res <- 0
    if(i>j){
        tmp <- Nk_formula(sig,gam,m,n,k=i-j)
        res <- 1/factorial(i-j) * tmp[i,j]
    }
    if(i==j) res <- 1
    return(res)
}

G_formula <- function(sig,gam,m,n){
    G_f <- matrix(nrow = n+m, ncol=n+m)
    for(i in 1:(m+n)){
        for(j in 1:(m+n)){
            G_f[i,j] <- G_ij(sig,gam,m,n,i,j)
        }
    }
    G_f[is.na(G_f)] <- 0
    return(G_f)
}

G_f <- G_formula(sig,gam,m,n)
G <- expm(N)
sum( (G-G_f)^2 )
#

# M <- expm(A) 
# View(M)


# ----- Eigen values & vectors ----

eg <- eigen(A)
eg$values
round(eg$vectors,4)


# ----- Check M(k,p) -----
# k > p 
k <- 5
p <- 3

M <- (D%^%(k-p)) %*% (N%^%p)


Dk_formula_i <- function(sig,gam,m,n,k,i){
    res1 <- ifelse(i<=m, 1, 0) * (-sig)^k
    res2 <- ifelse(i>=m+1, 1, 0) * (-gam)^k
    return(res1+res2)
}

Dk_formula <- function(sig,gam,m,n,k){
    res <- matrix(ncol=m+n, nrow=m+n, data = 0)
    for(i in 1:(m+n)){
        res[i,i] <- Dk_formula_i(sig,gam,m,n,k,i)
    }
    return(res)
}

Dk <- D%^%k
Dk_f <- Dk_formula(sig,gam,m,n,k)
Dk_f - Dk

M_formula_i <- function(sig,gam,m,n,k,p,i) {
    
    # -- CORRECT!
    # res1 <- Dk_formula_i(sig,gam,m,n,k=k-p,i)
    # res2 <- Nk_i(sig,gam,m,n,i,k=p)
    # return(res1*res2)
    
    
    res1 <- (-sig)^(k-p) * ifelse(i<=m,1,0) + (-gam)^(k-p)*ifelse(i>=m+1,1,0)
    # res1 <- Dk_formula_i(sig,gam,m,n,k=k-p,i)
    # res2 <- Nk_i(sig,gam,m,n,i,k=p)
    res2a <- sig^p * ifelse(p+1 <= i & i <= m+1, 1, 0)
    res2b <- gam^(i-m-1) * sig^(p-i+m+1) * ifelse(m+2 <= i & i <= m+p, 1, 0)
    res2c <- gam^p * ifelse(m+1+p <= i & i <= m+n, 1, 0)
    res2 <- res2a + res2b + res2c
    return(res1*res2)
    
    
    # -- FAUX
    # res1a <- ifelse(p+1 <= i & i<=m, 1,0)
    # res1b <- ifelse(i <= m , 1,0) * (sig/gam)^(m+1-i)
    # res1 <- sig^k *( res1a + res1b )
    # 
    # res2a <- ifelse(i==m+1, 1,0) * (sig/gam)^p
    # res2b <- ifelse(i>=m+1, 1,0) * (sig/gam)^(m+1-i+p)
    # res2c <- ifelse(i>=m+1+p, 1,0)
    # res2 <- gam^k * (res2a + res2b + res2c)
    # 
    # res <- (-1)^(k-p) * (res1 + res2)
    # return(res)
}

M_f <- matrix(ncol = m+n, nrow = m+n)
for(i in (p+1):(m+n)){
    M_f[i,i-p] <- M_formula_i(sig,gam,m,n,k,p,i)
}
M_f[is.na(M_f)] <- 0

M
M - M_f


# ----- Check exp(A) ----

exp_A_i1 <- function(sig,gam,m,n,i){
    res <- 0
    
    if(i <= m) { 
        res <- (sig)^(i-1)/factorial(i-1) * exp(-sig)
    }
    if(i >= m+1){
        #res <- ( (1+(sig/gam)^m) * (gam)^(i-1) / factorial(i-1) * exp(-gam) )
        res <- (gam)^(i-1) / factorial(i-1) * (sig/gam)^m * (1 + ((0))*ifelse(i==m+1,1,0)) * exp(-gam)
    }
    return(res)
}

EXP_A <- expm(A)
EXP_A_col1 <- EXP_A[,1]

exp_a_formula <- vector(length = m+n) 

for(i in 1:(m+n)) 
    exp_a_formula[i] <- exp_A_i1(sig,gam,m,n,i)

EXP_A_col1
exp_a_formula - EXP_A_col1
round(exp_a_formula/EXP_A_col1-1, 3)

plot(EXP_A_col1, log='y', typ='b', cex=2, lwd=2, las=1)
lines(exp_a_formula, pch=16, typ='b', col='blue1', lty=2, lwd=2)
abline(v=m, lty=2)


# ----- Numerical integration ----


hz <- 50
x <- GI.seminr(latent_mean = 1/sig,
               infectious_mean = 1/gam,
               R0 = 2.8,
               nE = m, 
               cal.times.fwdbck = c(1,5,15),
               nI = n,
               horizon = hz,
               dt = 0.1,
               I.init = 1E-5
)

g <- x$intrinsic

g_formula <- function(t,sig,gam,m,n){
    
    s <- 0
    for(i in (m+1):(m+n)) {
        s <- s + (gam*t)^(i-1)/factorial(i-1)
    }
    res <- s * (sig/gam)^m * exp(-gam*t) 
    return(res)}

tt <- seq(1,hz, by=0.5)
gf <- numeric(length = length(tt))
for(i in 1:length(tt)){
    gf[i] <- g_formula(t=tt[i],sig,gam,m,n)
}
gf <- gf / sum(gf)

 plot(g$tsi, g$density, typ='l', lwd=3,
      ylim = range(g$density, gf))
lines(x=tt, y=gf, col='red')
