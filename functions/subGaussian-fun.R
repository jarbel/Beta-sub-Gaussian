
library("gsl")
# Function: double hyperg_1F1 (double a, double b, double x)
# Function: int gsl_sf_hyperg_1F1_e (double a, double b, double x, gsl_sf_result * result)
# These routines compute the confluent hypergeometric function 1F1(a,b,x) = M(a,b,x) for general parameters a, b.

F_conf = function(a,b,x){
  hyperg_1F1(a,b,x)
}
# F_conf = function(a,b,x,log=0){
#   # returns the LOG of the confluent hypergeom fun with log = 1
#   Re(kummerM(x,a,b,lnchf=log))
# }

find_minimizer = function(alp,bet){
  if (alp>bet){
    alp_b = alp
    alp = bet
    bet = alp_b
  }
  if(alp==bet){
    0
  }
  else{
    fun_to_minimize = function(x){
      log(F_conf(alp,alp+bet,x))-
        0.5*alp*x/(alp+bet)*(1+F_conf(alp+1,alp+bet+1,x)/F_conf(alp,alp+bet,x))
    }
    # minimizer
    lower = .01
    # lower = (bet/S-1/2)*10*sqrt(S) #.1 # the pb here can be a postive value of lower bound due to num approx
    upper = 500
    # upper = max((bet/S-1/2)*30*sqrt(S),tan((aa/S-1/2)*pi)*3*sqrt(S))
    uni = uniroot(f = fun_to_minimize, lower = lower, upper = upper)
    # lower = .1 # the pb here can be a postive value of lower bound due to num approx
    # uni = uniroot(f = fun_to_minimize, lower = lower, upper = 100)
    uni$root
  }
}


# variance of a beta
variance_beta = function(alp,bet){
  if(alp==0|bet==0){
    0
  }
  else
  {
    S=alp+bet
    alp*bet/(S^2*(S+1))
  }
}

var_proxy_opt = function(alp,bet) {
  if(alp==bet){
    variance_beta(alp,bet)
  }
  else{
    if(alp==0|bet==0){
      0
    }
    else{
      if (alp>bet){
        alp_b = alp
        alp = bet
        bet = alp_b
      }
      x0 = find_minimizer(alp,bet)
      # optimal proxy variance sigma^2
      alp/(x0*(alp+bet))*(F_conf(alp+1,alp+bet+1,x0)/F_conf(alp,alp+bet,x0)-1)
    }
  }
}

# bound by Bercu et al.
psifun = function(nu){
  (1-nu^2)/abs(log(nu))
}

gfun = function(p){
  if(p==1/2)
    2
  else
    log(p/(1-p))/(2*p-1)
}
gfun = Vectorize(gfun, vectorize.args = "p")


bercu = function(alpha, S){
  nu = variance_beta(alpha, S-alpha)
  psifun(nu)
}
bercu = Vectorize(bercu, vectorize.args = "alpha")

#Hoeffding?
hoef = function(alpha, S){
  mu = alpha/S
  if(mu!=.5)
    (1-2*mu)/2/log(1/mu-1)
  else
    .25
}
hoef = Vectorize(hoef, vectorize.args = "alpha")
