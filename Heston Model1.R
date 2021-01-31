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
i = 1i

Call = function(S, K, r, t, sigma, del) 
{
  d1 = (log(S/K) + (r - del + sigma^2/2)*t) / (sigma*sqrt(t))
  d2 = d1 - sigma*sqrt(t)
  return(list(S*exp(-del*t)*pnorm(d1)- K*exp(-r*t)*pnorm(d2), pnorm(d1), pnorm(d2)))
}
bs_price = Call(St,K,r,t,sqrt(theta),0)[[1]]

f1<-function(u){
  b1 = k+lambda-rho*sigma
  b2 = k+lambda
  u1 = 1/2
  u2 = -1/2
  d1 = sqrt((rho*sigma*u*i-b1^2)-sigma^2*(2*u1*u*i-u^2))
  d2 = sqrt((rho*sigma*u*i-b2^2)-sigma^2*(2*u2*u*i-u^2))
  g1 = (b1-rho*sigma*u*i+d1)/(b1-rho*sigma*u*i-d1)
  g2 = (b2-rho*sigma*u*i+d2)/(b2-rho*sigma*u*i-d2)
  A1 = r*u*i*t+k*theta/(sigma^2)*((b1-rho*sigma*u*i+d1)*t-2*log((1-g1*exp(d1*t))/(1-g1)))
  A2 = r*u*i*t+k*theta/(sigma^2)*((b2-rho*sigma*u*i+d2)*t-2*log((1-g2*exp(d2*t))/(1-g2)))
  B1 = (b1-rho*sigma*u*i+d1)/(sigma^2)*(1-exp(d1*t))/(1-g1*exp(d1*t))
  B2 = (b2-rho*sigma*u*i+d2)/(sigma^2)*(1-exp(d2*t))/(1-g2*exp(d2*t))
  phi1 = exp(A1+B1*vt+i*log(St)*u)
  phi2 = exp(A2+B2*vt+i*log(St)*u)
  return(Re(exp(-i*u*log(K))*phi1/(i*u)))
}
f2<-function(u){
  b1 = k+lambda-rho*sigma
  b2 = k+lambda
  u1 = 1/2
  u2 = -1/2
  d1 = sqrt((rho*sigma*u*i-b1^2)-sigma^2*(2*u1*u*i-u^2))
  d2 = sqrt((rho*sigma*u*i-b2^2)-sigma^2*(2*u2*u*i-u^2))
  g1 = (b1-rho*sigma*u*i+d1)/(b1-rho*sigma*u*i-d1)
  g2 = (b2-rho*sigma*u*i+d2)/(b2-rho*sigma*u*i-d2)
  A1 = r*u*i*t+k*theta/(sigma^2)*((b1-rho*sigma*u*i+d1)*t-2*log((1-g1*exp(d1*t))/(1-g1)))
  A2 = r*u*i*t+k*theta/(sigma^2)*((b2-rho*sigma*u*i+d2)*t-2*log((1-g2*exp(d2*t))/(1-g2)))
  B1 = (b1-rho*sigma*u*i+d1)/(sigma^2)*(1-exp(d1*t))/(1-g1*exp(d1*t))
  B2 = (b2-rho*sigma*u*i+d2)/(sigma^2)*(1-exp(d2*t))/(1-g2*exp(d2*t))
  phi1 = exp(A1+B1*vt+i*log(St)*u)
  phi2 = exp(A2+B2*vt+i*log(St)*u)
  return(Re(exp(-i*u*log(K))*phi2/(i*u)))
}

Call_H = function(S,K,r,t,sigma,rho,vt,k,lambda,theta){
  P1 = 1/2+ integrate(f1,lower=0,upper=100)$value/pi
  P2 = 1/2+ integrate(f2,lower=0,upper=100)$value/pi
  return(S*P1-K*exp(-r*t)*P2)
}

Call_H(St,K,r,t,sigma,rho,vt,k,lambda,theta)

