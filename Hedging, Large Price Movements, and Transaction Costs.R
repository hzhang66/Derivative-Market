S = 100
K = 100
delta = 0
r = 0.05/(360*96)
drift = 0.2/(360*96)
sigma = 0.3/sqrt(360*96)
time_periods = 96*30


#part a
stocks = matrix(nrow=1,ncol=time_periods)
stocks[,1] = S
set.seed(42)
for (i in 2:time_periods)
{
  stocks[,i] = drift*stocks[,i-1]+sigma*stocks[,i-1]*rnorm(1)+stocks[,i-1]
}

#part b
Call = function(S, K, r, t, sigma, del) 
{
  d1 = (log(S/K) + (r - del + sigma^2/2)*t) / (sigma*sqrt(t))
  d2 = d1 - sigma*sqrt(t)
  return(list(S*exp(-del*t)*pnorm(d1)- K*exp(-r*t)*pnorm(d2), pnorm(d1), pnorm(d2)))
}

bs_price = vector(length = time_periods)
bs_deltas = vector(length = time_periods)
for(i in (1:time_periods))
{
  t = 60/365 - (i-1)/(96*365)
  temp = Call(stocks[1,i], K, drift*(360*96), t, sigma*sqrt(360*96), delta)
  bs_price[i] = temp[[1]]
  bs_deltas[i] = temp[[2]]
}


#part c
repl_bs = function(bs_price, bs_deltas, stock_price)
{
  replicating = vector(length = time_periods)
  replicating[1] = bs_price[1]
  bond = vector(length = time_periods)
  bond[1] = replicating[1] - bs_deltas[1]*stock_price[1]
  
  for(i in (2:time_periods))
  {
    replicating[i] = bs_deltas[i-1]*stock_price[i] + bond[i-1]*exp(r)
    bond[i] = replicating[i] - bs_deltas[i]*stock_price[i]
  }
  
  plot((1:length(stock_price)/96),stock_price, type="l", ylab = "Stock Price", xlab = "Day")
  plot((1:length(stock_price)/96),replicating, type="l", col="red", ylab = "Price", xlab = "Day") +
    lines((1:length(stock_price)/96),bs_price, col="green")
  legend("bottomleft", legend=c("Replicating Portfolio", "Black Scholes Price"), col=c("red", "green"), lty=1:1)
  plot((1:length(stock_price)/96),bs_deltas, type="l", ylab = "Hedge Ratio", xlab = "Day")
}

repl_bs(bs_price, bs_deltas, stocks[1,])

#part d
stocks_down = stocks[1,]
stocks_down[time_periods/2] = stocks_down[time_periods/2]*.9
bs_price_down = bs_price
bs_deltas_down = bs_deltas
bs_price_down[time_periods/2] = Call(stocks_down[time_periods/2], K, drift*(360*96), t, sigma*sqrt(360*96), delta)[[1]]
bs_deltas_down[time_periods/2] = Call(stocks_down[time_periods/2], K, drift*(360*96), t, sigma*sqrt(360*96), delta)[[2]]
for(i in ((time_periods/2)+1):time_periods)
{
  stocks_down[i] = drift*stocks_down[i-1]+sigma*stocks_down[i-1]*rnorm(1)+stocks_down[i-1]
  t = 60/365 - (i-1)/(96*365)
  temp = Call(stocks_down[i], K, drift*(360*96), t, sigma*sqrt(360*96), delta)
  bs_price_down[i] = temp[[1]]
  bs_deltas_down[i] = temp[[2]]
}
repl_bs(bs_price_down, bs_deltas_down, stocks_down)

#part e
stocks_up = stocks[1,]
stocks_up[time_periods/2] = stocks_up[time_periods/2]*1.1
bs_price_up = bs_price
bs_deltas_up = bs_deltas
bs_price_up[time_periods/2] = Call(stocks_up[time_periods/2], K, drift*(360*96), t, sigma*sqrt(360*96), delta)[[1]]
bs_deltas_up[time_periods/2] = Call(stocks_up[time_periods/2], K, drift*(360*96), t, sigma*sqrt(360*96), delta)[[2]]
for(i in ((time_periods/2)+1):time_periods)
{
  stocks_up[i] = drift*stocks_up[i-1]+sigma*stocks_up[i-1]*rnorm(1)+stocks_up[i-1]
  t = 60/365 - (i-1)/(96*365)
  temp = Call(stocks_up[i], K, drift*(360*96), t, sigma*sqrt(360*96), delta)
  bs_price_up[i] = temp[[1]]
  bs_deltas_up[i] = temp[[2]]
}
repl_bs(bs_price_up, bs_deltas_up, stocks_up)

#part f
replicating = vector(length = time_periods)
replicating[1] = bs_price[1]
bond = vector(length = time_periods)
bond[1] = replicating[1] - bs_deltas[1]*stocks[1,1]

diff_deltas = abs(diff(bs_deltas))
for(i in (2:time_periods))
{
  replicating[i] = bs_deltas[i-1]*stocks[1,i] + bond[i-1]*exp(r) - diff_deltas[i]*stocks[1,i]*0.002
  bond[i] = replicating[i] - bs_deltas[i]*stocks[1,i]
}


plot((1:length(stocks[1,])/96),stocks[1,], type="l", ylab = "Stock Price", xlab = "Day")
plot((1:length(stocks[1,])/96),replicating, type="l", col="red", ylab = "Price", xlab = "Day") +
  lines((1:length(stocks[1,])/96),bs_price, col="green")
legend("bottomleft", legend=c("Replicating Portfolio", "Black Scholes Price"), col=c("red", "green"), lty=1:1)
plot((1:length(stocks[1,])/96),bs_deltas, type="l", ylab = "Hedge Ratio", xlab = "Day")
