library(stats)

# Bachalier
normSpread <- function(Spot, Strike, Vol,  Topt, r, isCall=TRUE){
  d = (Spot - Strike) / (Vol * sqrt(Topt))
  p = (Spot - Strike) * pnorm(d) + Vol * sqrt(Topt) * dnorm(d)
  
  # Call
  if(isCall){
    return(p * exp(-r * Topt))
  }else{
    # Put
    (p - (Spot - Strike)) * exp(-r * Topt)  
  }
}

# Spread option payoff max(S1 - S1 - K, 0) using Monte Carlo
spreadOption_MC <-function(
  Zvec, # Random Number
  S1 = 3800, S2 = 3400,
  rho = 0.9, sig1 = 0.2, sig2 = 0.3, r = 0,
  Tmat = 0.5,  K = 400, isCall = TRUE)
{
  z1 = Zvec[1:(length(Zvec)/2)]
  z2 = rho*z1 + sqrt(1-rho^2)* Zvec[(length(Zvec)/2+1):length(Zvec)]
  S1terminal = S1*exp( (r-0.5*sig1^2)*Tmat + sig1*sqrt(Tmat)*z1)
  S2terminal = S2*exp( (r-0.5*sig2^2)*Tmat + sig2*sqrt(Tmat)*z2)
  payoff = (S1terminal - S2terminal - K)
  if(isCall){
    payoff[payoff < 0] = 0
  }else{
    payoff = -payoff
    payoff[payoff < 0] = 0
  }
  PVvec = exp(-r*Tmat)*payoff
  
  Out = c(PV=mean(PVvec), stderr = sd(PVvec)/sqrt(length(PVvec)))
  Out
}


###### From Murex Asian Spread ########
spreadOptionPV <- function(S1, S2, V1, V2, rho, r, Tmat, K, isCall){
  exp(-r*Tmat)*spreadOptionFV(S1, S2, V1, V2, rho, Tmat, K, isCall)
}

# Future value of Spread Option, ignore Discount Factor
#   @m1,m2  spot price1, 2
#   @V1,V2  
#   @psi,   % Sum of weights corresponding to spent fixings, default = 0
#   @currFixings,  % average Fixings already past
spreadOptionFV <- function(m1, m2, V1, V2, rhoA, T, K, isCall, psi=NULL, currFixings=NULL){
  if( !is.null(psi)) {
    m1 = (1-psi)*m1;
    m2 = (1-psi)*m2;
    K = K - psi*currFixings;
  }
  
  eps = ifelse(isCall == TRUE, 1, -1)
  if( K < 0 ){   # Flip options
    eps = -eps;  # Flip call to put, vice versa
    K = -K;
    # Full exchange
    flippedParam = flipOneToTwo(m1, m2, V1, V2 );
    m1 = flippedParam['m1']
    m2 = flippedParam['m2']
    V1 = flippedParam['V1']
    V2 = flippedParam['V2']
  }
  
  
  inf_limit = Inf;
  #inf_limit = 1000;
  #int1 = quadgk(@(x) integrand1(x, m1, m2, V1, V2, rhoA, T, K, eps), -inf_limit,inf_limit);   
  #int2 = quadgk(@(x) integrand2(x, m1, m2, V1, V2, rhoA, T, K, eps), -inf_limit,inf_limit);   
  #int3 = quadgk(@(x) integrand3(x, m1, m2, V1, V2, rhoA, T, K, eps), -inf_limit,inf_limit);   
  
  #integral(fun, xmin, xmax, method = c("Kronrod","Richardson","Clenshaw","Simpson","Romberg"), vectorized = TRUE, arrayValued = FALSE, reltol = 1e-8, abstol = 0, ...)
  myFun1 = function(x) integrand1(x, m1, m2, V1, V2, rhoA, T, K, eps)
  myFun2 = function(x) integrand2(x, m1, m2, V1, V2, rhoA, T, K, eps)
  myFun3 = function(x) integrand3(x, m1, m2, V1, V2, rhoA, T, K, eps)
  int1 = integrate(myFun1, -inf_limit, inf_limit)$value   
  int2 = integrate(myFun2, -inf_limit, inf_limit)$value   
  int3 = integrate(myFun3, -inf_limit, inf_limit)$value
  
  # My way
  term1 = m1 * int1;
  term2 = -m2 * int2;
  term3 = -K * int3;
  
  PV = eps*(term1 + term2 + term3);
  return(PV)
}

flipOneToTwo<- function(m1, m2, V1, V2 ){
  mTmp = m1;
  m1 = m2;
  m2 = mTmp;
  
  vTmp = V1;
  V1 = V2;
  V2 = vTmp;
  
  Out = c(m1=m1, m2=m2, V1=V1, V2=V2)
  Out
}


integrand1 <- function(y, m1, m2, V1, V2, rho, T, K, eps ){
  X = eps*(-intermD(y, T, m1, m2, K, V1, V2, rho) + V1*sqrt(T)*sqrt(1-rho^2));
  out = dnorm(y-V1*sqrt(T)*rho) * pnorm(X);
  #X = eps*(-intermD(y+V1*sqrt(T)*rho, T, m1, m2, K, V1, V2, rho) + V1*sqrt(T)*sqrt(1-rho^2));
  #out = normpdf(y) .* normcdf(X);
  out
}

integrand2 <- function(y, m1, m2, V1, V2, rho, T, K, eps ){
  X = -eps*intermD(y, T, m1, m2, K, V1, V2, rho);
  out = dnorm(y-V2*sqrt(T)) * pnorm(X);
  #X = -eps*intermD(y+V2*sqrt(T), T, m1, m2, K, V1, V2, rho);
  #out = normpdf(y) .* normcdf(X);        
  out
}

integrand3 <- function(y, m1, m2, V1, V2, rho, T, K, eps ){
  X = -eps*intermD(y, T, m1, m2, K, V1, V2, rho);
  out = dnorm(y) * pnorm(X);
  out
}



intermD <- function(yvec, T, m1, m2, K, V1, V2, rhoA){
  if (K < 0){
    stop('intermD: K -ve !')
  }
  
  num = K + m2 * exp(-0.5*V2^2*T + V2*sqrt(T)*yvec);
  denom = m1 * exp(-0.5*V1^2*T + V1*sqrt(T)*rhoA*yvec);
  dvec = 1/(V1 * sqrt(T) * sqrt(1-rhoA^2))* log(num / denom);
  
  #dvec = 1/(V1 * sqrt(T) * sqrt(1-rhoA^2))*( ...
  #    log(m2) + -0.5*V2^2*T + V2.*sqrt(T).*yvec - log(denom) );
  
  dvec
}
