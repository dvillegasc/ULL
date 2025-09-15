# Example 1
# Generating some random values with
# known mu
y <- rULL(n=300, mu=0.5, sigma= 0.5)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=ULL)

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ ULL
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu <- exp(-0.5 + 1 * x1)
  sigma <- exp(-0.5 + 1 * x2)
  y <- rULL(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

datos <- gendat(n=300)

mod2 <- gamlss(y~x1, family=ULL, data=datos)
summary(mod2)

# Example 3
# The first dataset measured the concentration of air pollutant CO
# in Alberta, Canada from the Edmonton Central (downtown)
# Monitoring Unit (EDMU) station during 1995.
# Measurements are listed for the period 1976–1995.
# Taken from Bicer et al. (2024) page 12.

data1 <- c(0.19, 0.20, 0.20, 0.27, 0.30,
           0.37, 0.30, 0.25, 0.23, 0.23,
           0.26, 0.23, 0.19, 0.21, 0.20,
           0.22, 0.21, 0.25, 0.25, 0.19)

mod3 <- gamlss(data1 ~ 1, family=ULL)

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod3, what="mu"))
exp(coef(mod3, what="sigma"))
#UMB= 0.8452
#ULL= 1.3758

# Extraction of the log likelihood
logLik(mod3)
#UMB= 19.2862
#ULL= 9.15843

# Example 4
# The second data set measured air quality monitoring of the
# annual average concentration of the pollutant benzo(a)pyrene (BaP).
# The data were obtained from the Edmonton Central (downtown)
# Monitoring Unit (EDMU) location in Alberta, Canada, in 1995.
# Taken from Bicer et al. (2024) page 12.

data2 <- c(0.22, 0.20, 0.25, 0.15, 0.38,
           0.18, 0.52, 0.27, 0.27, 0.27,
           0.13, 0.15, 0.24, 0.37, 0.20)

mod4 <- gamlss(data2 ~ 1, family=ULL)

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod4, what="mu"))
exp(coef(mod4, what="sigma"))
#UMB= 0.8593
#ULL= 1.3864

# Extraction of the log likelihood
logLik(mod4)
#UMB= 12.4300
#ULL= 6.3668

# Replicating figure 5 from Bicer et al. (2024)
# Hist and estimated pdf of Data-I and Data-II
mu1 <- exp(coef(mod3, what="mu"))
mu2 <- exp(coef(mod4, what="mu"))
sigma1 <- exp(coef(mod3, what="sigma"))
sigma2 <- exp(coef(mod4, what="sigma"))

par(mfrow = c(1, 2))

# Data-I
hist(data1, freq = FALSE,
     xlim = c(0, 1.0), ylim = c(0, 10),
     main = "Histogram of Data-I",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")

curve(dULL(x, mu = mu1, sigma = sigma1), add = TRUE,
      from = 0.0000001, to = 0.999,
      col = "blue", lwd = 2)

legend("topright", legend = c("ULL"),
       col = c("blue"), lwd = 2, bty = "n")

# Data-II
hist(data2, freq = FALSE,
     xlim = c(0, 1.0), ylim = c(0, 6),
     main = "Histogram of Data-II",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")

curve(dULL(x,  mu = mu2, sigma = sigma2), add = TRUE,
      from = 0.0000001, to = 0.999,
      col = "blue", lwd = 2)

legend("topright",
       legend = c("ULL"),
       col = c("blue"),
       lwd = 2,
       bty = "n")


par(mfrow = c(1, 1))

# Example 5
# The third dataset measured the concentration of sulphate
# in Calgary from 31 different periods during 1995.
# Taken from Bicer et al. (2024) page 13.

data3 <- c(0.048, 0.013, 0.040, 0.082, 0.073, 0.732, 0.302,
           0.728, 0.305, 0.322,  0.045, 0.261, 0.192,
           0.357, 0.022, 0.143, 0.208, 0.104, 0.330, 0.453,
           0.135, 0.114, 0.049, 0.011, 0.008, 0.037, 0.034,
           0.015, 0.028, 0.069, 0.029)

mod5 <- gamlss(data3 ~ 1, family=ULL)

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod5, what="mu"))
exp(coef(mod5, what="sigma"))
#UMB= 1.5820
#ULL= 0.000001

# Extraction of the log likelihood
logLik(mod5)
#UMB= 23.3506
#ULL= 23.2327

# Example 6
# The fourth dataset measured the concentration of pollutant CO in Alberta, Canada
# from the Calgary northwest (residential) monitoring unit (CRMU) station during 1995.
# Measurements are listed  for the period 1976-95.
# Taken from Bicer et al. (2024) page 13.

data4 <- c(0.16, 0.19, 0.24, 0.25, 0.30, 0.41, 0.40,
           0.33, 0.23, 0.27, 0.30, 0.32, 0.26, 0.25,
           0.22, 0.22, 0.18, 0.18, 0.20, 0.23)

mod6 <- gamlss(data4 ~ 1, family=ULL)

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod6, what="mu"))
exp(coef(mod6, what="sigma"))
#UMB= 0.8161
#ULL= 2.2250

# Extraction of the log likelihood
logLik(mod6)
#UMB= 17.9885
#ULL= 8.5914

# Replicating figure 6 from Bicer et al. (2024)
# Hist and estimated pdf of Data-III and Data-IV
mu3 <- exp(coef(mod5, what="mu"))
mu4 <- exp(coef(mod6, what="mu"))
sigma3 <- exp(coef(mod5, what="sigma"))
sigma4 <- exp(coef(mod6, what="sigma"))

par(mfrow = c(1, 2))

# Data-III
hist(data3, freq = FALSE,
     xlim = c(0, 1.0), ylim = c(0, 10),
     main = "Histogram of Data-III",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")

curve(dULL(x, mu = mu3, sigma= sigma3), add = TRUE,
      from = 0.0000001, to = 0.999,
      col = "blue", lwd = 2)

legend("topright", legend = c("ULL"),
       col = c("blue"), lwd = 2, bty = "n")

# Data-IV
hist(data4, freq = FALSE,
     xlim = c(0, 1.0), ylim = c(0, 6),
     main = "Histogram of Data-IV",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")

curve(dULL(x, mu = mu4, sigma= sigma4), add = TRUE,
      from = 0.0000001, to = 0.999,
      col = "blue", lwd = 2)

legend("topright",
       legend = c("ULL"),
       col = c("blue"),
       lwd = 2,
       bty = "n")


par(mfrow = c(1, 1))


