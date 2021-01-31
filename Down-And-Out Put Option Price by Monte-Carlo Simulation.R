##2 Down and Out Puts
library(ggplot2)
S = 100
K = 95
Sb = 75
delta = 0
t = 1/4
r = 0.05/(360*96)

out_paths = function(sigma){
  sigma = sigma/sqrt(360*96)
  r = 0.05/(360*96)
  stocks = matrix(nrow=1000,ncol=96*90)
  stocks[,1] = 100
  count = ncol(stocks)
  set.seed(42)
  for (i in 2:count){
    stocks[,i] = r*stocks[,i-1]+sigma*stocks[,i-1]*rnorm(1000)+stocks[,i-1]
    stocks[(stocks[,i] <= 75),i]= 0
  }
  counter = sum(stocks[,count]==0)/1000
  stocks[stocks[,count]==0,count] = 100
  return(list(counter,stocks[,count]))
}
Put = function(S, K, r, t, sigma, del) 
{
  sigma = sigma/sqrt(360*96)
  d1 = (log(S/K) + (r - del + sigma^2/2)*t) / (sigma*sqrt(t))
  d2 = d1 - sigma*sqrt(t)
  return(list(-S*exp(-del*t)*pnorm(-d1)+ K*exp(-r*t)*pnorm(-d2), pnorm(d1), pnorm(d2)))
}

#part a.
price = mean(exp(-r*360*96*t)*pmax(0,K-out_paths(.2)[[2]]))
range = quantile(exp(-r*360*96*t)*pmax(0,K-out_paths(.2)[[2]]),c(.025,.975))
bs_put = Put(S,K,r*360*96,t,.2*sqrt(360*96),delta)[[1]]

#part b.
out_5 = out_paths(.05)[[1]]
out_10 = out_paths(.1)[[1]]
out_15 = out_paths(.15)[[1]]
out_20 = out_paths(.2)[[1]]
out_25 = out_paths(.25)[[1]]
out_30 = out_paths(.3)[[1]]
out_35 = out_paths(.35)[[1]]
out_40 = out_paths(.4)[[1]]

fractions = c(out_5,out_10,out_15,out_20,out_25,out_30,out_35,out_40)*100
plot(seq(5,40,5),fractions,type = "l",xlab="Volatility",ylab = "Proportion Out")

#part c
price_5 = out_paths(.05)[[2]]
price_10 = out_paths(.1)[[2]]
price_15 = out_paths(.15)[[2]]
price_20 = out_paths(.2)[[2]]
price_25 = out_paths(.25)[[2]]
price_30 = out_paths(.3)[[2]]
price_35 = out_paths(.35)[[2]]
price_40 = out_paths(.4)[[2]]


last_cols = cbind(price_5,price_10,price_15,price_20,price_25,price_30,price_35,price_40)
# prices = exp(-r*360*96*t)*pmax(0,K-prices)

prices=vector(length = ncol(last_cols))
conf_int=list()
for(i in (1:ncol(last_cols)))
{
  prices[i] = mean(exp(-r*360*96*t)*pmax(0,K-last_cols[,i]))
  conf_int = append(conf_int,
                    quantile(exp(-r*360*96*t)*pmax(0,K-last_cols[,i]),c(.025,.975)))
}
vols = seq(.05,.4,.05)
bs_puts = Put(S,K,r*360*96,t,vols*sqrt(360*96),delta)[[1]]

plot(vols,bs_puts,type="l",ylim = c(0,15))+lines(vols,prices,col="red")+
  lines(vols,conf_int[seq(1,15,2)],lty=2)+lines(vols,conf_int[seq(2,16,2)],lty=2)



