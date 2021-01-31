St = 100
K = 100
r = 0.04
t = .5
sigma = .3
rho = -.5
vt = .01
k = 6
lambda = 0
theta = 0.02
callHestoncf = function(S, K, t, r, q, v0, theta, rho, k, sigma,
                         implVol = FALSE) {
  if (sigma < 0.01)
    sigma = 0.01
  P1 = function(u,S,K,t,r,q,v0,theta,rho,k,sigma) {
    p = Re(exp(-1i * log(K) * u) *
              cfHeston(u - 1i, S, t, r, q, v0, theta, rho, k, sigma) /
              (1i * u * S * exp((r-q) * t)))
    p
  }
  P2 = function(u,S,K,t,r,q,v0,theta,rho,k,sigma) {
    p = Re(exp(-1i * log(K) * u) *
              cfHeston(u  ,S,t,r,q,v0,theta,rho,k,sigma) /
              (1i * u))
    p
  }
  cfHeston = function(u,S,t,r,q,v0,theta,rho,k,sigma) {
    d = sqrt((rho * sigma * 1i * u - k)^2 + sigma^2 *
                (1i * u + u ^ 2))
    g = (k - rho * sigma * 1i * u - d) /
      (k - rho * sigma * 1i * u + d)
    cf1 = 1i * u * (log(S) + (r - q) * t)
    cf2 = theta*k/(sigma^2)*((k - rho * sigma * 1i * u - d) *
                             t - 2 * log((1 - g * exp(-d * t)) / (1 - g)))
    cf3 = v0 / sigma^2 * (k - rho * sigma * 1i * u - d) *
      (1 - exp(-d * t)) / (1 - g * exp(-d * t))
    cf  = exp(cf1 + cf2 + cf3)
    cf
  }
  
  ## pricing
  vP1 = 0.5 + 1/pi * integrate(P1,lower = 0, upper = Inf,
                                S, K, t, r, q, v0, theta, rho, k, sigma)$value
  vP2 = 0.5 + 1/pi * integrate(P2,lower = 0, upper = Inf,
                                S, K, t, r, q, v0, theta, rho, k, sigma)$value
  result = exp(-q * t) * S * vP1 - exp(-r * t) * K * vP2;
  
  ## implied BSM vol
  if (implVol) {
    diffPrice = function(vol,call,S,K,t,r,q){
      d1 = (log(S/K)+(r - q + vol^2/2)*t)/(vol*sqrt(t))
      d2 = d1 - vol*sqrt(t)
      callBSM = S * exp(-q * t) * pnorm(d1) -
        K * exp(-r * t) * pnorm(d2)
      call - callBSM
    }
    impliedVol = uniroot(diffPrice, interval = c(0.0001, 2),
                          call = result, S = S, K = K,
                          t = t, r = r, q = q)[[1L]]
    result = list(value = result, impliedVol = impliedVol)
  }
  result
}

Call = function(S, K, r, t, sigma, del) 
{
  d1 = (log(S/K) + (r - del + sigma^2/2)*t) / (sigma*sqrt(t))
  d2 = d1 - sigma*sqrt(t)
  return(list(S*exp(-del*t)*pnorm(d1)- K*exp(-r*t)*pnorm(d2), pnorm(d1), pnorm(d2)))
}
#Part a
h_price = callHestoncf(St,K,t,r,0,vt,theta,rho,k,sigma)
bs_price = Call(St,K,r,t,sqrt(theta),0)[[1]]

#Part b
stocks = seq(70,130)
diffs = rep(0,length(stocks))
impliedVols =  rep(0,length(stocks))
for (i in 1:length(diffs)){
  diffs[i] = callHestoncf(stocks[i],K,t,r,0,vt,theta,rho,k,sigma)-Call(stocks[i],K,r,t,sqrt(theta),0)[[1]]
  impliedVols[i]= callHestoncf(stocks[i],K,t,r,0,vt,theta,rho,k,sigma,implVol=T)$impliedVol
}
plot(stocks,diffs,type="l",ylab="Difference in Price",xlab="Stock Price")
#part c
plot(stocks,impliedVols,type="l",ylab="Implied Vol",xlab="Stock Price")
#part d
rho = .5
stocks = seq(70,130)
diffs = rep(0,length(stocks))
impliedVols =  rep(0,length(stocks))
for (i in 1:length(diffs)){
  diffs[i] = callHestoncf(stocks[i],K,t,r,0,vt,theta,rho,k,sigma)-Call(stocks[i],K,r,t,sqrt(theta),0)[[1]]
  impliedVols[i]= callHestoncf(stocks[i],K,t,r,0,vt,theta,rho,k,sigma,implVol=T)$impliedVol
}
plot(stocks,impliedVols,type="l",ylab="Implied Vol",xlab="Stock Price")
