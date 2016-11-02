probpool.default <- function(env.pool, bio.pool, dist.pool ...)
{
  
  
  env.pool <- env.pool
  bio.pool <- bio.pool
  dis.pool <- dis.pool
  
  phi.i    <- 
  
  y <- as.numeric(y)
  est <- linmodEst(x, y)
  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$call <- match.call()
  
  
  class(probpool) <- "probpool"
  est
}



calc.phi.i  <- function (x)
{
  calc(calc(r, fun=function(x){x * 2})
  
}

