# R script for implementing upper and lower bounds of n (last modified 8/31/2016)
# Originally provided by RSM in C++. Modified by LHU in C++. 
# This version written by AA. 

# tavare's g
g <- function(n,j,t) {
  value <- 0
  if (j > n) return(0)
  else if (n == 0 | j == 0) return(0)
  else if (t > 5000) {
    if (j == 1) return(1)
    else return(0)
  }
  else if (t >= 0.1 & n >= 90 & j >= 50) return(0)
  else if (t >= 1 & n >= 20 & j >= 10) return(0)
  else {
    for (k in j:n) {
      exponent <- exp((-k)*(k-1)*t/2)
      sign <- (-1)^(k-j)
      if (k == 1) {
        val <- exponent * sign
      }
      else {
        num <- (2 * k - 1) * choose(j + k  - 2, j) 
        			* choose (n - 1, k - 1) * choose (k - 1, j - 1)
        den = (n + k - 1) * choose(n + k - 2, n)
        val = exponent * sign * num / den
      }
      value <- value + val
    }
    if (value < 10e-50) value= 0
    return(value) 
  } 
}

# polyphyly function `POLY'
f <- function(x) {
  (x*(x-1)*(x-2)*(3*x-5)-24)/(3*x*(x-1)^2*(x-2))
}

# ---------------------------------
# Doing lower bounds
# ---------------------------------
# first lower bound for n `LOW(q,k,Ttot)'
lowTot <- function(q, k, t) {
  return(max(log(1-q) / log(g(k-1, k-1, t) * f(k)),1))
}

# second lower bound (1) the tavare term
tavare <- function(k, t) {
  x <- 1
  for(i in 2:floor(k / 2)) {
    x = x * g(i, i, t)^(floor(k / i))
  }
  x = x * g(k-2,k-2, t)
  return(x)
}

# second lower bound (1b) the corrected tavare term
fkt <- function(k, t) {
  x <- g(floor(k / 2),floor(k / 2),t)^(floor(k / 2))
  for(i in (floor(k / 2) + 1):k - 2) {
    x = x * g(i,i,t)
  }
  return(x)
}

# second lower bound for n `LOW(q,k,Tmax)'
lowMax <- function(q, k, t) {
  return(max(log(1-q)/log(fkt(k,t) * f(k)),1))
}


# --------------------------------------------------------
# Doing exact calculation for k = 4 (Uricchio et al. 2016)
# --------------------------------------------------------
# k = 4 exact (asymmetric: set T1=Ttot, symmetric: set T1+T2=Ttot) 
n <- function(q, t) {
  log(1 - q) / log(2 / 3 * exp(-t))
}

# ---------------------------------
# Doing upper bound 
# ---------------------------------
# tFac
tFac <- function(A, B) {
  (2 * A + B) * (B - 1) / (A * (A + 1))
}

# Factorial
factorial <- function(A) {
  fact <- 1.0
  while(A > 1) {
    fact <- fact * A
    A <- A - 1
  }
  return(fact)
}

# prMono
prMono <- function(lA, lB, tA, tB) {
  ret <- 0
  for(lAi in 1:lA) {
    for(lBi in 1:lB) {
      ret <- ret + 2 * g(lA, lAi, tA) * g(lB, lBi, tB) * 
        (((factorial(lAi) * factorial(lBi)) / factorial(lAi + lBi))) / (lAi + lBi - 1) 
        					* tFac(lBi,lAi)
    }
  }
  return(ret)
}

# prRecipMono
prRecipMono <- function(lA, lB, tA, tB) {
  ret <- 0
  for(lAi in 1:lA) {
    for(lBi in 1:lB) {
      ret <- ret + g(lA, lAi, tA) * g(lB, lBi, tB) * 
        (((factorial(lAi) * factorial(lBi))/ factorial(lAi + lBi)))/(lAi + lBi - 1) 
    }
  }
  return(ret)
}

# prNotPoly
prNotPoly <- function(lA, lB, tA, tB) {  
  ret <- 0
  for(lAi in 1:lA) {
    for(lBi in 1:lB) {
      ret <- ret + g(lA, lAi, tA) * g(lB, lBi, tB) * 
        (2 * (factorial(lAi) * factorial(lBi) / factorial(lAi+lBi)) * 
           ((lAi+lBi)/(lAi*(lAi+1))+(lAi+lBi)/(lBi*(lBi+1))-1/(lAi+lBi-1)))
    }
  }
  return(ret)
}

# upper bound for n
upp <- function(q, k, tA, tB) {
  min <- 1
  for(i in 2:(k - 2)) { 
    lA = k - i
    lB = i
    
    prM <- prNotPoly(lA, lB, tA, tB)
    if(prM < min) min <- prM
  }
  
  return(log((1-q)/(k - 3))/ log(1 - min))
}

# ---------------------------------------
# Doing old bound (Uricchio et al. 2016)
# ---------------------------------------
# old bound 
oldBound <- function(q, k, t) {
  return(log((1-q)/(k - 3))/ log(1 - g(k-2,1,t)))
} 

# -----------------------------------------
# Doing hypothetical lower bound (guessing)
# -----------------------------------------

# GUESS: improved lower bound for n
test <- function(q, k, t) {
  max <- 0
  for(j in seq(from = 1e-3, to = t, by = 1e-3)) {
    for(i in 2:(k - 2)) {  
      tA = j
      tB = t - j
      lA = k - i
      lB = i
      
      prM <- prNotPoly(lA, lB, tA, tB) 
      if(prM > max) max <- prM
    }
  }
  return(log(1-q)/log(1-max)) 
}
