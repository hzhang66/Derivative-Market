library(DataAnalytics) 
library(pracma)
straddle = function(St,K){
  abs(St-K)
}
digital = function(St,K){
  if(St>K) 1
  else 0
}
call = function(St,K){
  max(St-K,0)
}
put = function(St,K){
  max(K-St,0)
}

european = function(S,payoff,r,h,t,K){
  u = exp(r*h+.2*sqrt(h))
  d = exp(r*h-.2*sqrt(h))
  stock_prices = matrix(0,nrow=t+1,ncol=t+1)
  for (i in 1:nrow(stock_prices)){
    for(j in 1:ncol(stock_prices)){
      stock_prices[i,j]=S*u^(j-1)*d^((i-1)-(j-1))
    }
  }
  payoffs = matrix(0,nrow=t+1,ncol=t+1)
   for (i in 1:(t+1)){
    payoffs[t+1,i]=payoff(stock_prices[t+1,i],K)
   }
   deltas = matrix(0,nrow=t+1,ncol=t+1)
   bonds = matrix(0,nrow=t+1,ncol=t+1)
   for (i in 1:t){
     for(j in t:1){
     payoffs[t+1-i,j] = exp(-r*h)*(payoffs[t+2-i,j+1]*((exp(r*h)-d)/(u-d))+payoffs[t+2-i,j]*((u-exp(r*h))/(u-d)))
     deltas[t+1-i,j] = (payoffs[t+2-i,j+1]-payoffs[t+2-i,j])/((stock_prices[t+1-i,j])*(u-d))
     bonds[t+1-i,j] = exp(-r)*((payoffs[t+2-i,j]*u-payoffs[t+2-i,j+1]*d)/(u-d))
     }
   }
   deltas[upper.tri(deltas)] = 0
   payoffs[upper.tri(payoffs)] = 0
   bonds[upper.tri(bonds)] = 0
   stock_prices[upper.tri(stock_prices)] = 0

   print(stock_prices)
   print(payoffs)
   print(bonds)
   print(deltas)
}
straddle_1 = european(100,straddle,.02,1/4,4,90)
straddle_2 = european(100,straddle,.02,1/40,40,90)
digital_price = european(100,digital,.02,1/4,4,90)
put1 = european(100,put,.02,1/4,4,90)


## Quesiton 3
american = function(S,payoff,r,h,t,K){
  u = exp(r*h+.15*sqrt(h))
  d = exp(r*h-.15*sqrt(h))
  p = (exp(r*h)-d)/(u-d)
  stock_prices = matrix(0,nrow=t+1,ncol=t+1)
  for (i in 1:nrow(stock_prices)){
    for(j in 1:ncol(stock_prices)){
      stock_prices[i,j]=S*u^(j-1)*d^((i-1)-(j-1))
    }
  }
  stock_prices[upper.tri(stock_prices)] = 0
  payoffs_no = matrix(0,nrow=t+1,ncol=t+1)
  payoffs_ex = matrix(0,nrow=t+1,ncol=t+1)
  payoffs = matrix(0,nrow=t+1,ncol=t+1)
  deltas = matrix(0,nrow=t+1,ncol=t+1)
  bonds = matrix(0,nrow=t+1,ncol=t+1)
  for (i in 1:(t+1)){
    payoffs_no[t+1,i]=payoff(stock_prices[t+1,i],K)
    payoffs_ex[t+1,i]=payoff(stock_prices[t+1,i],K)
    payoffs[t+1, i] = payoff(stock_prices[t+1,i],K)
  }
  
  ## first row is payoff today, ie option price
  exercise_node = list()
  for (i in seq(t,1,-1)){
    for(j in 1:i){
      payoffs_no[i,j] = exp(-r*h)*(payoffs[i+1,j+1]*((exp(r*h)-d)/(u-d))
                                    +payoffs[i+1,j]*((u-exp(r*h))/(u-d)))
      payoffs_ex[i,j] = payoff(stock_prices[i,j],K)
      payoffs[i, j] = max(payoffs_no[i,j], payoffs_ex[i,j])
      
      if (payoffs_no[i,j] < payoffs_ex[i,j])
        exercise_node = append(exercise_node, paste(i, j))
      
      deltas[i,j] = (payoffs[1+i,j+1]-payoffs[1+i,j])/((stock_prices[i,j])*(u-d))
      bonds[i,j] = exp(-r)*((payoffs[1+i,j]*u-payoffs[1+i,j+1]*d)/(u-d))
    }
  }
  payoffs_no[upper.tri(payoffs_no)] = 0
  payoffs_ex[upper.tri(payoffs_ex)] = 0
  payoffs[upper.tri(payoffs_ex)] = 0
  #print(stock_prices)
  #print(payoffs_no)
  #print(payoffs_ex)
  #print(payoffs)
  return(list(payoffs[1,1], exercise_node, deltas, bonds))
}
american_put = american(10,put,.01,1/365,250,10)
american_call = american(10,call,.01,1/365,250,10)


## Question 4
discrete_dividend = function(S,payoff,nature,t,r,h,K,dividend,dates,sigma){
  u = exp(sigma*sqrt(h))
  d = exp(-sigma*sqrt(h))
  p = (exp(r*h)-d)/(u-d)
  stock_prices = matrix(0,nrow=t+1,ncol=t+1)
  for (i in 1:nrow(stock_prices)){
    for(j in 1:ncol(stock_prices)){
      if (i %in% dates)
        stock_prices[i,j]=(1-dividend)*S*u^(j-1)*d^((i-1)-(j-1))
      else
        stock_prices[i,j]=S*u^(j-1)*d^((i-1)-(j-1))
    }
  }
  payoffs_no = matrix(0,nrow=t+1,ncol=t+1)
  payoffs_ex = matrix(0,nrow=t+1,ncol=t+1)
  payoffs = matrix(0,nrow=t+1,ncol=t+1)
  deltas = matrix(0,nrow=t+1,ncol=t+1)
  bonds = matrix(0,nrow=t+1,ncol=t+1)
  
  if (nature == 'European'){
    payoffs = matrix(0,nrow=t+1,ncol=t+1)
    for (i in 1:(t+1)){
      payoffs[t+1,i]=payoff(stock_prices[t+1,i],K)
    }
  }
  else{
    for (i in 1:(t+1)){
      payoffs_no[t+1,i]=payoff(stock_prices[t+1,i],K)
      payoffs_ex[t+1,i]=payoff(stock_prices[t+1,i],K)
      payoffs[t+1, i] = payoff(stock_prices[t+1,i],K)
    }
    exercise_node = list()
    for (i in seq(t,1,-1)){
      for(j in 1:i){
        payoffs_no[i,j] = exp(-r*h)*(payoffs[i+1,j+1]*((exp(r*h)-d)/(u-d))
                                     +payoffs[i+1,j]*((u-exp(r*h))/(u-d)))
        payoffs_ex[i,j] = payoff(stock_prices[i,j],K)
        payoffs[i, j] = max(payoffs_no[i,j], payoffs_ex[i,j])
        
        if (payoffs_no[i,j] < payoffs_ex[i,j])
          exercise_node = append(exercise_node, paste(i, j))
      }
    }
  
}
  ## first row is payoff today, ie option price
  for (i in seq(t,1,-1)){
    for(j in 1:i){
      deltas[i,j] = (payoffs[1+i,j+1]-payoffs[1+i,j])/((stock_prices[i,j])*(u-d))
      bonds[i,j] = exp(-r)*((payoffs[1+i,j]*u-payoffs[1+i,j+1]*d)/(u-d))
    }
  }
  
  #print(stock_prices)
  #print(payoffs_no)
  #print(payoffs_ex)
  #print(payoffs)
  return(list(payoffs[1,1], exercise_node, deltas, bonds))
}

american_call_dividend = discrete_dividend(10,call,'American', 200, 0.02,1/365,10,0.05,c(50,100,150), 0.2)
american_put_dividend = discrete_dividend(10,put,'American', 200, 0.02,1/365,10,0.05,c(50,100,150), 0.2)
american_straddle_dividend = discrete_dividend(10,straddle,'American', 200, 0.02,1/365,10,0.05,c(50,100,150), 0.2)

##Question 5

asian_opt = function(S,payoff,t,r,h,K,sigma)
{
  t= t/h
  no_paths = 100000
  u = exp(r*h+sigma*sqrt(h))
  d = exp(r*h-sigma*sqrt(h))
  p = (exp(r*h)-d)/(u-d)
  stock_prices = matrix(0,nrow=t+1,ncol=t+1)
  for (i in 1:nrow(stock_prices)){
    for(j in 1:ncol(stock_prices)){
      stock_prices[i,j]=S*u^(j-1)*d^((i-1)-(j-1))
    }
  }
  stock_prices[upper.tri(stock_prices)] = NA
  
  set.seed(9)
  payoffs = vector()
  
  
  for (i in 1:no_paths)
  {
    rands = runif(t)
    rands[rands<=p] = 1
    rands[rands!=1] = 0
    price = vector(length = (t+1))
    price[1] = S
    col = 1
    for(j in (2:(t+1)))
    {
      if(rands[j-1] == 1)
      {
        col=col+1
        price[j] = stock_prices[j, col]
      }
      else
      {
        price[j] = stock_prices[j, col]
      }
    }
    
    payoffs[i]=payoff(mean(price),K)
  }
  opt_price = exp(-r)*mean(payoffs)
  sd = sd(payoffs)
  
  return(list(opt_price,c(exp(-r)*quantile(payoffs,0.025), exp(-r)*quantile(payoffs,0.975))))
}

asian_prop = asian_opt(200,call, 1, 0.02,1/365, 220, 0.2)

##Question 6

lsmput = function(S,t,r,N,K,sigma,n_path){
  dt = t/N
  discount<-exp(-r*dt)
  discountVet<-exp(-r*dt*c(0:N))
  a = array(0, dim = c(3,1))
  nudt = (r-0.5*sigma^2)*dt
  sidt = sigma*sqrt(dt)
  RandMat = matrix(rnorm(round(n_path/2)*N),nrow=round(n_path/2))
  Increments = rbind(nudt + sidt*RandMat, nudt - sidt*RandMat)
  Paths = as.matrix(cbind(log(S)*ones(n_path,1) , Increments))
  LogPaths<-cumsum(Paths[1,])
  for (i in 2:nrow(Paths)){
      a=cumsum(Paths[i,])
      LogPaths<-rbind(LogPaths,a)
  }
  SPaths = exp(LogPaths)
  CashFlows = pmax(K-SPaths,0)[,N+1]
  ExerciseTime = N+1*ones(n_path,1)
  
  for (i in seq(N,1,-1)){
    InMoney = unname(which(SPaths[,i] < K))
    XData = SPaths[InMoney,i]
    RegrMat = unname(as.matrix(cbind(XData, XData^2)))
    RegrMat1 = unname(as.matrix(cbind(ones(length(XData),1), XData, XData^2)))
    YData = unname(CashFlows[InMoney]) * discountVet[ExerciseTime[InMoney] +1  - i]
    out<-lm(YData~RegrMat[,1]+RegrMat[,2])
    a<-out$coefficients
    IntrinsicValue = K - XData 
    ContinuationValue = RegrMat1 %*% a 
    Exercise = which(IntrinsicValue > ContinuationValue) 
    k = InMoney[Exercise]
    CashFlows[k] = IntrinsicValue[Exercise] 
    ExerciseTime[k] = i
  }
  return(mean(CashFlows*discountVet[ExerciseTime]))
}

#1
price=lsmput(200, 1,0.1,250,220,0.3,100000)
price100k = price
price10k = price

#2
num_path=c('10','100','1000','10000','100000')
lsmprice_npath=c(lsmput(200, 1,0.1,250,220,0.3,10),
           lsmput(200, 1,0.1,250,220,0.3,100),
           lsmput(200, 1,0.1,250,220,0.3,1000),
           price10k,
           price100k)
plot(num_path,lsmprice_npath,type="l",xlab="Number of Paths",ylab="Price")

#3
num_period=c('3','10','100','250','1000')
lsmprice_N=c(lsmput(200, 1,0.1,3,220,0.3,100000),
           lsmput(200, 1,0.1,10,220,0.3,100000),
           lsmput(200, 1,0.1,100,220,0.3,100000),
           price100k,
           lsmput(200, 1,0.1,1000,220,0.3,100000))
plot(num_period,lsmprice_N)