ULL <- function(mu.link = "log", sigma.link = "log") {
  mstats <- checklink("mu.link", "ULL",
                      substitute(mu.link), c("log", "identity"))
  dstats <- checklink("sigma.link", "ULL",
                      substitute(sigma.link), c("log", "identity"))
  
  structure(
    list(
      family = c("ULL", "Unit Log–Lindley"),
      parameters = list(mu = TRUE, sigma = TRUE),
      nopar = 2,
      type = "Continuous",
      mu.link = as.character(substitute(mu.link)),
      sigma.link = as.character(substitute(sigma.link)),
      mu.linkfun = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
      mu.linkinv = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
      mu.dr = mstats$mu.eta,
      sigma.dr = dstats$mu.eta,
      
      # First derivates
      
      dldm = function(y, mu, sigma) {
        dldm <- 1 / (mu - log(y)) - sigma / (1 + mu * sigma)
        return(dldm)
      },
      
      dldd = function(y, mu, sigma) {
        dldd <- 2 / sigma - mu / (1 + mu * sigma) + log(y)
        return(dldd)
      },
      
      # Second derivates 
      
      d2ldm2 = function(y, mu, sigma) {
        dldm <- 1 / (mu - log(y)) - sigma / (1 + mu * sigma)
        d2ldm2 <- -dldm * dldm
        d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
        return(d2ldm2)
      },
      
      d2ldd2 = function(y, mu, sigma) {
        dldd <- 2 / sigma - mu / (1 + mu * sigma) + log(y)
        d2ldd2 <- -dldd * dldd
        d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
        return(d2ldd2)
      },
      
      d2ldmdd = function(y, mu, sigma) {
        dldm <- 1 / (mu - log(y)) - sigma / (1 + mu * sigma)
        dldd <- 2 / sigma - mu / (1 + mu * sigma) + log(y)
        d2ldmdd <- -dldm * dldd
        d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
        return(d2ldmdd)
      },
      
      
      G.dev.incr = function(y, mu, sigma, ...) -2 * dULL(y, mu, sigma, log = TRUE),
      rqres = expression(rqres(pfun = "pULL", type = "Continuous", y = y, mu = mu, sigma = sigma)),
      
      mu.initial = expression(mu <- rep(1, length(y))),
      sigma.initial = expression(sigma <- rep(1, length(y))),
      
      mu.valid = function(mu) all(mu >= 0),
      sigma.valid = function(sigma) all(sigma > 0),
      y.valid = function(y) all(y > 0 & y < 1)
    ),
    class = c("gamlss.family", "family")
  )
}

