dULL <- function(x, mu, sigma, log = FALSE) {
  if (any(mu < 0)) stop("parameter mu must be positive!")
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  if (any(x <= 0 | x >= 1)) stop("x must be in the interval (0, 1)")
  
  log_pdf <- 2 * log(sigma) - log(1 + mu * sigma) + 
    log(mu - log(x)) + (sigma - 1) * log(x)
  
  if (log) {
    return(log_pdf)
  } else {
    return(exp(log_pdf))
  }
}


pULL <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu < 0)) stop("parameter mu must be positive!")
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  if (any(q <= 0 | q >= 1)) stop("parameter q must be positive!")
  
  # Auxiliar function
  auxfun <- function(q) {
    term1 <- (( q ^ sigma ) * (1 + sigma*(mu - log(q))))
    term2 <- 1 + mu * sigma
    res <- term1/term2
    res
  }
  
  cdf <- ifelse(q <= 0, 0,
                ifelse(q >= 1, 1,auxfun(q)))
  
  if (!lower.tail) {
    cdf <- 1 - cdf
  }
  
  if (log.p) {
    cdf <- log(cdf)
  }
  
  return(cdf)
}


qULL <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (any(p <= 0 | p >= 1)) stop("p must be in the interval (0, 1)")
  if (any(mu < 0))        stop("mu must be positive")
  if (any(sigma <= 0))    stop("sigma must be positive")
  
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  
  # Begin auxiliary function
  aux_fun <- function(x, p, mu, sigma){
    res <- p - pULL(x, mu, sigma)
    return(res)
  }
  # End auxiliary function
  
  res <- uniroot(aux_fun, interval=c(.Machine$double.eps, 1 - .Machine$double.eps),
                 p=p, mu=mu, sigma=sigma)$root
  
  return(res)
}
qULL <- Vectorize(qULL)


rULL <- function(n, mu, sigma) {
  if (any(mu < 0)) stop("parameter mu must be positive!")
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  u <- runif(n)
  qULL(u, mu, sigma)
}


