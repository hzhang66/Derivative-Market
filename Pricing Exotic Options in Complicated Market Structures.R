library(matrixStats)
T <- 1 
N <- 250 #step size
h <- 1/250 #step size
time <- seq(from=0, to=T, by=h)
nSim <- 1000 #number of paths

r <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
r[,1] <- 0.05
S1 <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
S1[1:nSim,1] <- 10
S2 <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
S2[1:nSim,1] <- 10
sigma11 <- 0.1 
sigma12 <- 0.2 
sigma21 <- 0.3 
beta<-0.05
alpha<-0.6
delta<-0.1

#a

#SDE
for(i in 1:nSim){
  for (j in 2:(N+1)){
    B1 <- rnorm(1, mean=0, sd=1)
    B2 <- rnorm(1, mean=0, sd=1)
    S1_t <-S1[i,(j-1)] + S1[i,(j-1)]*r[i,(j-1)]*h+ sigma11*sqrt(S1[i,(j-1)])*sqrt(h)*B1+sigma12*S1[i,(j-1)]*sqrt(h)*B2
    S2_t <-S2[i,(j-1)]+S2[i,(j-1)]*r[i,(j-1)]*h+sigma21*(S1[i,(j-1)]-S2[i,(j-1)])*sqrt(h)*B1
    r_t <- r[i,(j-1)]+alpha*(beta-r[i,(j-1)])*h+delta*sqrt(r[i,(j-1)])*sqrt(h)*B1
    r[i,j]<-r_t
    S1[i,j]<-S1_t 
    S2[i,j]<-S2_t 
    
  }
}

hist(r[1:nSim,N+1], breaks = 20,main="Distribution of r_T at T=1", col="yellow", xlab="(r1 | r0 = 0.05)")

#b

T <- 100 
h <- 1/52
N <- T/h 
time <- seq(from=0, to=T, by=h) 
nSim <- 1 #number of paths
r <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
r[,1] <- 0.05
S1 <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
S1[1:nSim,1] <- 10
S2 <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
S2[1:nSim,1] <- 10

for(i in 1:nSim){
  for (j in 2:(N+1)){
    B1 <- rnorm(1, mean=0, sd=1)
    B2 <- rnorm(1, mean=0, sd=1)
    S1_t <-S1[i,(j-1)] + S1[i,(j-1)]*r[i,(j-1)]*h+ sigma11*sqrt(S1[i,(j-1)])*sqrt(h)*B1+sigma12*S1[i,(j-1)]*sqrt(h)*B2
    S2_t <-S2[i,(j-1)]+S2[i,(j-1)]*r[i,(j-1)]*h+sigma21*(S1[i,(j-1)]-S2[i,(j-1)])*sqrt(h)*B1
    r_t <- r[i,(j-1)]+alpha*(beta-r[i,(j-1)])*h+delta*sqrt(r[i,(j-1)])*sqrt(h)*B1
    r[i,j]<-r_t
    S1[i,j]<-S1_t 
    S2[i,j]<-S2_t 
    
  }
}

plot(time, r, main="Trajectory of r_t from t = 0 to T = 100 with step=1/52", col="blue",ylab="r_t",xlab="time",type="l")


#c

T <- 0.5 
h <- 1/250
N <- T/h 
time <- seq(from=0, to=T, by=h) 
nSim <- 10000
K=10

r <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
r[,1] <- 0.05
S1 <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
S1[1:nSim,1] <- 10
S2 <- matrix(rep(0, (N+1)*nSim), nrow = nSim)
S2[1:nSim,1] <- 10


for(i in 1:nSim){
  for (j in 2:(N+1)){
    B1 <- rnorm(1, mean=0, sd=1)
    B2 <- rnorm(1, mean=0, sd=1)
    S1_t <-S1[i,(j-1)] + S1[i,(j-1)]*r[i,(j-1)]*h+ sigma11*sqrt(S1[i,(j-1)])*sqrt(h)*B1+sigma12*S1[i,(j-1)]*sqrt(h)*B2
    S2_t <-S2[i,(j-1)]+S2[i,(j-1)]*r[i,(j-1)]*h+sigma21*(S1[i,(j-1)]-S2[i,(j-1)])*sqrt(h)*B1
    r_t <- r[i,(j-1)]+alpha*(beta-r[i,(j-1)])*h+delta*sqrt(r[i,(j-1)])*sqrt(h)*B1
    r[i,j]<-r_t
    S1[i,j]<-S1_t 
    S2[i,j]<-S2_t 
    
  }
}
  
Payoff= mean(pmax(S1[,N+1]-K,0))
price=Payoff *exp(-0.05*T)
price
  
#d
maxS <- rep(0, nSim) 
for(i in 1:nSim){
  maxS[i]=max(max(S1[i,]),max(S2[i,]))
}
Payoff= mean(pmax(maxS-K,0))
price=Payoff *exp(-0.05*T)
price
