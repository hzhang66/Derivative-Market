##Question 1: BS and Monte Carlo Simulation for Calls
library(ggplot2)
S = 100
K = 100
sigma = 0.2/sqrt(360*96)
delta = 0
r = 0.05/(360*96)
t = 1/4

stocks = matrix(nrow=5,ncol=96*90)
stocks[,1] = 100
count = ncol(stocks)
for (i in 2:count){
  stocks[,i] = r*stocks[,i-1]+sigma*stocks[,i-1]*rnorm(5)+stocks[,i-1]
}
plot(stocks[1,],type="l",ylim=c(80,120))+lines(stocks[2,],col="blue")+lines(stocks[3,],col="red")+
  lines(stocks[4,],col="yellow")+lines(stocks[5,],col="purple")

Call = function(S, K, r, t, sigma, del) 
{
  d1 = (log(S/K) + (r - del + sigma^2/2)*t) / (sigma*sqrt(t))
  d2 = d1 - sigma*sqrt(t)
  return(list(S*exp(-del*t)*pnorm(d1)- K*exp(-r*t)*pnorm(d2), pnorm(d1), pnorm(d2)))
}
bs_price = Call(S,K,r*360*96,t,sigma*sqrt(360*96),delta)[[1]]

set.seed(42)
MS_Call = function(S,K,r,t,sigma,N){
  prev_stocks = rep(S,N)
  curr_stocks = rep(NA,N)
  count = 96*90
  for (i in 2:count){
    curr_stocks = r*prev_stocks+sigma*prev_stocks*rnorm(N)+prev_stocks
    prev_stocks = curr_stocks
  }
  payoff = pmax(0,curr_stocks-K)
  price =exp(-r*360*96*t)*mean(payoff)
  range = exp(-r*360*96*t)*quantile(payoff,probs=c(.025,.975))
  return(list(price,range))
  
}
MS_Call(S,K,r,t,sigma,100)
MS_Call(S,K,r,t,sigma,1000)
MS_Call(S,K,r,t,sigma,1E6)
MS_Call(S,K,r,t,sigma,1E8)

