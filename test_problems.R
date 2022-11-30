
# ********************************************************************************
# One-dimensional test problems (scaled to be used over the probability simplex)
# ********************************************************************************

# from Infinity 77: stationary test function
test.problem02 <- function( weight.vector ) {
  stopifnot( length(weight.vector)==2 && sum(weight.vector)==1 )
  x = weight.vector[1]
  x = x * (7.5-2.7) + 2.7
  y = sin(x) + sin( x * 10/3 )
}

# from Infinity 77: less stationary test function
test.problem13 <- function( weight.vector ) {
  stopifnot( length(weight.vector)==2 && sum(weight.vector)==1 )
  x = weight.vector[1]
  x = x * (0.99-0.001) + 0.001
  y = -x^(2/3) - (1 - x^2)^(1/3)
}

# from Infinity 77: less stationary test function
test.problem15 <- function( weight.vector ) {
  stopifnot( length(weight.vector)==2 && 1-sum(weight.vector)<.Machine$double.eps )
  x = weight.vector[1]
  x = x * (5+5) - 5
  y = (x^2 - 5*x + 6)/(x^2 + 1)
}


# modified Xiong function (from Hebbal et al. 2021): non-stationary test function
test.modifiedXiong <- function( weight.vector ) {
  stopifnot( length(weight.vector)==2 && 1-sum(weight.vector)<.Machine$double.eps )
  x = weight.vector[1]
  term1 = sin(40*(x-0.85)^4)*cos(2.5*(x-0.95))
  term2 = 0.5*(x-0.9)
  test.modifiedXiong = -0.5 * ( term1 + term2 + 1 )
}


# ********************************************************************************
# Two-dimensional test problem (scaled to be used over the probability simplex)
# ********************************************************************************
# Bird function (from Infinity 77)
test.bird <- function( weight.vector ) {
  stopifnot( length(weight.vector)==3 && 1-sum(weight.vector)<.Machine$double.eps )
  x = weight.vector[1:2]
  x = (x * 4* pi) - 2 * pi
  f1 = (x[1]-x[2])^2
  f2 = exp( (1-sin(x[1]))^2 ) * cos(x[2])
  f3 = exp( (1-cos(x[2]))^2 ) * sin(x[1])
  f = f1 + f2 + f3
}