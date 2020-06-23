##Example code for saving ----
if (F) {
  sims4 <- run_sims(bvals = c(100, 200, 400), rvals = c(0.009, 0.016, 0.02845),
                    avals = c(10**-11, 10**-10, 10**-9), kvals = c(10**9),
                    tauvals = c(22.5, 33.75, 50.625))
  #Save results so they can be re-loaded in future
  write.csv(sims4[[1]], "./Celia/sims4_1.csv", row.names = F)
  if (!is.null(sims4[[2]])) {write.csv(sims4[[2]], "./Celia/sims4_2.csv", row.names = F)}
  if (!is.null(sims4[[3]])) {write.csv(sims4[[3]], "./Celia/sims4_3.csv", row.names = F)}
} else {
  #Load results previously simulated
  temp1 <- read.csv("./Celia/sims4_1.csv", stringsAsFactors = F)
  if ("./Celia/sims4_2.csv" %in% list.files()) {
    temp2 <- read.csv("./Celia/sims4_2.csv", stringsAsFactors = F)
  } else {temp2 <- NULL}
  if ("./Celia/sims4_3.csv" %in% list.files()) {
    temp3 <- read.csv("./Celia/sims4_3.csv", stringsAsFactors = F)
  } else {temp3 <- NULL}
  sims4 <- list(temp1, temp2, temp3)
}

#Okabe and Ito 2008 colorblind-safe qualitative color scale ----
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

## Import libraries ----

library(deSolve)
library(tidyr)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(caret)

## Define derivatives function ----
derivs <- function(t, y, parms) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Issue warning about too small/negative yvals (if warnings is 1)
  if (parms["warnings"]==1 & any(y < parms["thresh_min_dens"])) {
    warning(paste("pop(s)", paste(which(y < parms["thresh_min_dens"]),
                                  collapse = ","), 
                  "below thresh_min_dens, treating as 0"))
  }
  
  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  
  #Create output vector
  dY <- c(S = 0, I = 0, P = 0)
  
  ##Calculate dS
  
  #V3 (logistic dS/dt) (including competition from I pop)
  #dS/dt = rS((K-S-c*I)/K) - aSP
  dY["S"] <- parms["r"] * y["S"] *
    ((parms["K"] - y["S"] - parms["c"] * y["I"])/parms["K"]) -
    parms["a"] * y["S"] * y["P"]
  
  ##Calculate dI
  #dI/dt = aSP - aS(t-tau)P(t-tau)
  if (t < parms["tau"]) {
    dY["I"] <- parms["a"] * y["S"] * y["P"]
  } else {
    dY["I"] <- parms["a"] * y["S"]*y["P"] -
      parms["a"] * lagvalue(t - parms["tau"], 1)*lagvalue(t - parms["tau"], 3)
  }
  
  ##Calculate dP
  #dP/dt = baS(t-tau)P(t-tau) - aSP
  if (t < parms["tau"]) {
    dY["P"] <- -parms["a"] * y["S"] * y["P"]
  } else {
    dY["P"] <- parms["b"] * parms["a"] *
      lagvalue(t-parms["tau"], 1)*lagvalue(t-parms["tau"], 3) -
      parms["a"]*y["S"]*y["P"]
  }
  
  #Issue warning about too large pop (if warnings is TRUE)
  if (parms["warnings"]==1 & any(y > 10**100)) {
    warning(paste("pop(s)",paste(which(y > 10**100), collapse = ","), 
                  "exceed max limit, 10^100, returning dY = 0"))
  }
  dY[y > 10**100] <- 0
  
  #From documentation: The return value of func should be a list, whose first 
  #element is a vector containing the derivatives of y with respect to time
  return(list(dY))
}

#Realistic parameter values: ----
#  
#   r ranges from 0.04 (a 17-minute doubling time) to 0.007 (a 90-minute doubling time).
#   K ranges from 10^6 to 10^10, although typically we focus on 10^8 to 10^9.
#   a ranges from 10^-12 to 10^-8.
#   tau ranges from 10 to 105 (mins).
#   b ranges from 5 to 1000.


#   Also, we'll have to calculate the total bacteria densty by addinc S and I
#   columns into a new one called B. We do thin in "yout", before we pivot_longer it.

## I'll try to plot the results with Density +10 and the wished colors ----
#
#Run simulation with "Density + 10"
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtD <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

## Now, I'l try to change the colors to the colorblind-friendly scale
#
#Run simulation with different colors
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtC <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

##Plot results with different colors
library(tidyr)
youtC$B <- youtC$S+youtC$I
youtC_plot <- pivot_longer(youtC, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtC_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and K = 10**9
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtK <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtK)

## We'll add new rows for each parameter here
youtK$r <- 0.04
youtK$a <- 10**-10
youtK$b <- 50
youtK$tau <- 10
youtK$K <- 10**9
youtK$c <- 1

##Plot results with different colors, Density + 10, and K = 10**9
library(tidyr)
youtK$B <- youtK$S+youtK$I
head(youtK)
youtK_plot <- pivot_longer(youtK, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtK_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and K = 10**8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**8, c = 1,
            warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtK2 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtK2)

## We'll add new rows for each parameter here
youtK2$r <- 0.04
youtK2$a <- 10**-10
youtK2$b <- 50
youtK2$tau <- 10
youtK2$K <- 10**8
youtK2$c <- 1

##Plot results with different colors, Density + 10, and K = 10**8
library(tidyr)
youtK2$B <- youtK2$S+youtK2$I
head(youtK2)
youtK2_plot <- pivot_longer(youtK2, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtK2_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and r = 0.03
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.03, a = 10**-10, b = 50, tau = 10, K = 10**9,
            c = 1,  warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtR <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtR)

## We'll add new rows for each parameter here
youtR$r <- 0.03
youtR$a <- 10**-10
youtR$b <- 50
youtR$tau <- 10
youtR$K <- 10**9
youtR$c <- 1

##Plot results with different colors, Density + 10, and r = 0.03
library(tidyr)
youtR$B <- youtR$S+youtR$I
head(youtR)
youtR_plot <- pivot_longer(youtR, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtR_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) + 
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and r = 0.007
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.007, a = 10**-10, b = 50, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 800, by = 1)
youtR2 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtR2)

## We'll add new rows for each parameter here
youtR2$r <- 0.007
youtR2$a <- 10**-10
youtR2$b <- 50
youtR2$tau <- 10
youtR2$K <- 10**9
youtR2$c <- 1

##Plot results with different colors, Density + 10, and r = 0.007
library(tidyr)
youtR2$B <- youtR2$S+youtR2$I
head(youtR2)
youtR2_plot <- pivot_longer(youtR2, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtR2_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and a = 10**-13
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-13, b = 50, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 11000, by = 1)
youtA <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtA)

## We'll add new rows for each parameter here
youtA$r <- 0.04
youtA$a <- 10**-13
youtA$b <- 50
youtA$tau <- 10
youtA$K <- 10**9
youtA$c <- 1

##Plot results with different colors, Density + 10, and a = 10**-13
library(tidyr)
youtA$B <- youtA$S+youtA$I
head(youtA)
youtA_plot <- pivot_longer(youtA, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtA_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and a = 10**-8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8, b = 50, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtA2 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtA2)

## We'll add new rows for each parameter here
youtA2$r <- 0.04
youtA2$a <- 10**-8
youtA2$b <- 50
youtA2$tau <- 10
youtA2$K <- 10**9
youtA2$c <- 1

##Plot results with different colors, Density + 10, and a = 10**-8
library(tidyr)
youtA2$B <- youtA2$S+youtA2$I
head(youtA2)
youtA2_plot <- pivot_longer(youtA2, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtA2_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and tau = 65
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 65, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 600, by = 1)
youtT <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtT)

## We'll add new rows for each parameter here
youtT$r <- 0.04
youtT$a <- 10**-10
youtT$b <- 50
youtT$tau <- 65
youtT$K <- 10**9
youtT$c <- 1

##Plot results with different colors, Density + 10, and tau = 65
library(tidyr)
youtT$B <- youtT$S+youtT$I
head(youtT)
youtT_plot <- pivot_longer(youtT, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtT_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and tau = 120
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 120, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 800, by = 1)
youtT2 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtT2)

## We'll add new rows for each parameter here
youtT2$r <- 0.04
youtT2$a <- 10**-10
youtT2$b <- 50
youtT2$tau <- 120
youtT2$K <- 10**9
youtT2$c <- 1

##Plot results with different colors, Density + 10, and tau = 120
library(tidyr)
youtT2$B <- youtT2$S+youtT2$I
head(youtT2)
youtT2_plot <- pivot_longer(youtT2, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtT2_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and b = 20
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 20, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtB1 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

## Now, we'll add the new columns for each parameter
youtB1$r <- 0.04
youtB1$a <- 10**-10
youtB1$K <- 10**9
youtB1$c <- 1
head(youtB1)

##Plot results with different colors, Density + 10, and b = 20
library(tidyr)
youtB1$B <- youtB1$S+youtB1$I
youtB1_plot <- pivot_longer(youtB1, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtB1_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and b = 500
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 500, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtB2 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtB2)

## We'll add the new columns for each parameter
youtB2$r <- 0.04
youtB2$a <- 10**-10
youtB2$b <- 500
youtB2$tau <- 10
youtB2$K <- 10**9
youtB2$c <- 1
head(youtB2)

##Plot results with different colors, Density + 10, and b = 500
library(tidyr)
youtB2$B <- youtB2$S+youtB2$I
youtB2_plot <- pivot_longer(youtB2, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtB2_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and b = 850
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 850, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtB3 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtB3)

## We'll add new columns for each parameter
youtB3$r <- 0.04
youtB3$a <- 10**-10
youtB3$b <- 850
youtB3$tau <- 10
youtB3$K <- 10**9
youtB3$c <- 1
head(youtB3)

##Plot results with different colors, Density + 10, and b = 850
library(tidyr)
youtB3$B <- youtB3$S+youtB3$I
youtB3_plot <- pivot_longer(youtB3, c(S, I, P, B), names_to = "Population", values_to = "Density")

tiff("plot.tiff")
ggplot(data = youtB3_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")
dev.off()

## Tiff is the function that converts our graph in an image document, which is
## a ".tiff" one. This is because we saw the legend disrupted here (some problem
## with the code of ggplot). Therefore, I can use this function and close i with
## de.off() to have an imagewhere I CAN clearly see the legend of the plot.


## Finding MAXIMUMS ----
## We have to make sure we also run the command for B!! Since it was introduced
## later than S, I, and P

## To add new columns to a data frame in a quick way, we'd use the following
## command:
youtB1 <- cbind(data.frame(b = 20, tau = 10), youtB1)
head(youtB1)

## Now, we'd like to use all the simulations I've run and combine them into one
## big data frame. Then, summarize them.

## First, I'll try it by adding the simulations for different values of b
bigB <- rbind(youtB1, youtB2, youtB3)
bigB

## Now, I'll do it for ALL the simulations I had run
bigFINAL <- rbind(bigB, youtK, youtK2, youtR, youtR2, youtA, youtA2, youtT, youtT2)
bigFINAL

## Now, I'll try to summarize the data in "bigFINAL"
bigFINAL_plot <- pivot_longer(bigFINAL, c(S, I, P, B), names_to = "Population", values_to = "Density")
ggplot(data = bigFINAL_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_point(lwd = 1) +
  scale_y_continuous(trans = "log10")
## I pivot_longer'd because I remembered that "summarize" works better after doing so
summarise(bigFINAL, maximum_B = max(youtB1$B, na.rm = T), max(youtB2$B, na.rm = T),
          max(youtB3$B, na.rm = T), max(youtK$B, na.rm = T), max(youtK2$B, na.rm = T),
          max(youtR$B, na.rm = T), max(youtR2$B, na.rm = T), max(youtA$B, na.rm = T),
          max(youtA2$B, na.rm = T), max(youtT$B, na.rm = T), max(youtT2$B, na.rm = T))


head(bigFINAL)

## Here, I achived to summarize every maximum of each of the simulations included
## in the big data frame

## Now, I'll try to find the maximum of the maximums

bigFINAL <- dplyr::group_by(bigFINAL, b, tau, a, r, K, c)
bigFINAL
FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), maxtime = time[B == maximum_B])
FINAL

## Let's try find the SLOPE of some simulations ----
## I'll try to find the slope for a = 10**-8
## Run simulation
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8, b = 50, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtA2S <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtA2S)

##Plot results
library(tidyr)
youtA2S$B <- youtA2$S+youtA2$I
head(youtA2S)
youtA2S_plot <- pivot_longer(youtA2S, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtA2S_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous(trans = "log10")

## Find the SLOPE of B in this plot
x <- 1:50
y <- youtA2S$B[1:50]

plot(x, log10(y))
mod <- lm(log10(y)~x) #this can be combined with mod$coefficients[2] within summarize
summary(mod)

cor(x, y)

attributes(mod)
## To pull up some the attributes we can use the $ sign
## For example, we extract the coefficients

plot(x, log10(y))
abline(mod, col = 4) #To include the regression line
mod$coefficients[2]

## Finding the slope with summarize ----

## This wasn't exactly the way to do it. The right way is the following:
## Finding the slope with summarize

## Here, we try to analyze and plot row by row of the data frame bigFINAL

FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), 
                          maxtime = time[B == maximum_B],
                          slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                       time[time < maxtime & B < 0.1*K])$coefficients[2],
                          intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                           time[time < maxtime & B < 0.1*K])$coefficients[1])
FINAL

for (row in 1:nrow(FINAL)) {
  bigfinal_rows <- which(FINAL$b[row] == bigFINAL$b & 
                           FINAL$tau[row] == bigFINAL$tau &
                           FINAL$a[row] == bigFINAL$a &
                           FINAL$r[row] == bigFINAL$r &
                           FINAL$K[row] == bigFINAL$K &
                           FINAL$c[row] == bigFINAL$c)
  print(ggplot(data = bigFINAL[bigfinal_rows, ],
               aes(x = time, y = B)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = FINAL$slope[row], intercept = FINAL$intercept[row],
                      color = "red") +
          geom_point(data = FINAL[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}  

## How to run the simulations with a LOOP ----

mybig <- NA #We do this to erase anything that's in this data frame
myc <- 1
for(myb in c(75, 760)){
  for(mya in c(10**-10.5, 10**-8.25)){
    for(myK in c(10**7.75, 10**8.8)){
      for(mytau in c(33, 95)){
        for(myr in c(0.009, 0.025)){
          timelength <- 250
          timestep <- 1
          runnumber <- 1
          yinit <- c(S = 10**6, I = 0, P = 10**4)
          params <- c(b = myb, a = mya, K = myK, tau = mytau, r = myr, c = myc,
                      warnings = 0, thresh_min_dens = 10**-100)
          keeprunning <- TRUE
          while(keeprunning){
            times <- seq(from = 0, to = timelength, by = timestep)
            myyout <- as.data.frame(dede(y = yinit, times = times, func = derivs, 
                                         parms = params))
            if(tail(myyout$S, 1) < 1000){
              keeprunning <- FALSE
              ## Here, we said that if the density of S is < 1000 we want to stop
              ## the simulation
            }else{
              keeprunning <- TRUE
              timelength <- timelength*2
              timestep <- timestep*2
              ## And, here, that, otherwise, we want it to keep running but doubling
              ## the time length and the timestep
            }
            if(timelength > 11000){
              keeprunning <- FALSE
              ## Finally, we indicate that if the simulation has been running for
              ## 11000, we want it to stop
            }
          }
          
          ## Here, we create the different columns in the data frame and we say 
          ## that they are equal to myb, for example, because we want both values
          ## in myb, not only one.
          myyout$B <- myyout$S + myyout$I
          myyout$burst <- myb
          myyout$infec <- mya
          myyout$K <- myK
          myyout$tau <- mytau
          myyout$r <- myr
          myyout$c <- myc
          myyout$uniq_run <- runnumber
          
          if(is.na(mybig)){
            mybig <- myyout
          }else{
            mybig <- rbind(mybig, myyout)
          }
          ## Here, we specified that, if mybig is empty, we want it to be equal to
          ## myyout (so, we just want to put myyout in it). However, if we have
          ## already run a simulation, and it's inside mybig, what we want is to
          ## add the second simulation to it by using rbind, making mybig bigger.
          
          runnumber <- 1 + runnumber
        }
      }
    }
  }
}

# Let's find the maximums and slopes of the simulations
mybig

# Firts, we summarize the data we have
mybig <- dplyr::group_by(mybig, burst, tau, infec, r, K, c)
mybig

# Let's find the slopes and the inetercepts, too
BIG <- dplyr::summarise(mybig, maximum_B = max(B), 
                              maxtime = time[B == maximum_B],
                              slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                           time[time < maxtime & B < 0.1*K])$coefficients[2],
                              intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                               time[time < maxtime & B < 0.1*K])$coefficients[1])
BIG

## Here, we try to analyze and plot row by row of the data frame mybig
for (row in 1:nrow(BIG)) {
  bigfinal_rows <- which(BIG$burst[row] == mybig$burst & 
                           BIG$tau[row] == mybig$tau &
                           BIG$infec[row] == mybig$infec &
                           BIG$r[row] == mybig$r &
                           BIG$K[row] == mybig$K &
                           BIG$c[row] == mybig$c)
  print(ggplot(data = mybig[bigfinal_rows, ],
               aes(x = time, y = B)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = BIG$slope[row], intercept = BIG$intercept[row],
                      color = "red") +
          geom_point(data = BIG[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}

## To see all the values in BIG:
View(BIG)

## Plotting all the parameters with maxtime in the 32 simulations to see how
## they affect it
ggplot(data = BIG, aes(x = log10(burst), y = maxtime, color = as.factor(K),
                       shape = as.factor(infec))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r)


## Let's run sims1 ----
sims1 <- run_sims(bvals = c(75, 480, 760), avals = c(10**-10.5, 10**-9, 10**-8.25),
                  kvals = c(10**7.75, 10**8.8), tauvals = c(33, 57, 95),
                  rvals = c(0.009, 0.018, 0.025))
length(sims1)
# When we use [] means that the output will be a list, while when we use [[]], means
# that it'll be a dataframe (explanation with trains and their cars).
class(sims1[[1]])
sims1[1]
sims1[[1]]
table(sims1[[1]]$Pop)

# By chencking group 2 and 3 we see if there has been any error: simulations that 
# didn'r reach the equilibrium (2), and simulations that failed (3).
sims1[[2]]
sims1[[3]]

## Notice that the summarize will be slightly different since the data has been
## pivot_longer'd. So, we have to use the density column and make a subset where 
## "Pop" = B

sub_sims1 <- subset(sims1[[1]], Pop == "B")
class(sub_sims1)

## Now, we want to find the maximum_B, the maxtime, and the slope of each simulation
group_sims1 <- dplyr::group_by(sub_sims1, uniq_run, b, tau, a, r, K, c)
group_sims1
class(group_sims1)

sum_sims1 <- dplyr::summarise(group_sims1, maximum_B = max(Density),                
                              maxtime = time[Density == maximum_B],
                              slope = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                           time[time < maxtime & Density < 0.1*K])$coefficients[2],
                              intercept = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                               time[time < maxtime & Density < 0.1*K])$coefficients[1])

for (row in 1:nrow(sum_sims1)) {
  bigfinal_rows <- which(sum_sims1$b[row] == group_sims1$b & 
                           sum_sims1$tau[row] == group_sims1$tau &
                           sum_sims1$a[row] == group_sims1$a &
                           sum_sims1$r[row] == group_sims1$r &
                           sum_sims1$K[row] == group_sims1$K &
                           sum_sims1$c[row] == group_sims1$c)
  tiff(paste("./Celia/Sims1_plots/", sum_sims1[row, "uniq_run"], ".tiff", sep = ""),
       width = 6, height = 6, units = "cm", res = 150)
  
  
  print(ggplot(data = group_sims1[bigfinal_rows, ],
               aes(x = time, y = Density)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims1$slope[row], intercept = sum_sims1$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims1[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
  
  dev.off()
}

group_sims1
sum_sims1
View(sum_sims1)

## Now, we want to make a gglpot that represents all the simulations with the
## summarized data
# We start by analyzing how the paramters affcet maxtime
ggplot(data = sum_sims1, aes(x = log10(b), y = maxtime, color = as.factor(K),
                       shape = as.factor(a))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r) +
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = "lm")
## lm function on our data

# How to see how tau and b affect maxtime as if they were INDEPENDENT from 
# each other?
reg_sims1 <- lm(maxtime ~ tau + b + a + r, data = sum_sims1)
summary(reg_sims1)

confint(reg_sims1)
# Calulate the RSE
sigma(reg_sims1)/mean(sum_sims1$maxtime) # The lower the RSE, the more accurate the model.
# The RSE estimate gives a measure of error of prediction.

# Analyzing the data as if the parameters were CORRELATED
# Split the data into training and test set
set.seed(123)
training.samples <- sum_sims1$maxtime %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- sum_sims1[training.samples, ]
test.data <- sum_sims1[-training.samples, ]

# The standard linear regression model can be computed as follow:
# Buil the model
reg_sims1 <- lm(maxtime ~ tau + b + a + r, data = train.data)
# Summarize the model
summary(reg_sims1)
# Make predictions
predictions <- reg_sims1 %>% predict(test.data)
# Make performance. (a) Prediction error, RMSE
RMSE(predictions, test.data$maxtime)
# (b) R2
R2(predictions, test.data$maxtime)

# INTERACTION EFFECTS
# Build the model
reg2_sims1 <- lm(maxtime ~ b + tau + a + r + b:a + a:r, data = sum_sims1)
summary(reg2_sims1)
# Make predictions
predictions <- reg2_sims1 %>% predict(test.data)
# Model performance (a) Prediction error, RMSE
RMSE(predictions, test.data$maxtime)
# (b) R2
R2(predictions, test.data$maxtime)


# What about maximum_B?
ggplot(data = sum_sims1, aes(x = log10(b), y = maximum_B, color = as.factor(K),
                             shape = as.factor(a))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r) +
  scale_y_continuous(trans = "log10")

# PLotting maxtime and maximum_B with r and K
ggplot(data = sum_sims1, aes(x = maxtime, y = maximum_B,
                             )) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(K ~ r) +
  scale_y_continuous(trans = "log10")



## Let's run sims2 ----
## Where we'll have c, K, and r constant, and we'll be changing between 5
## different values of a, b, and tau
sims2 <- run_sims(bvals = c(50, 100, 200, 400, 800),
                  avals = c(10**-12, 10**-11, 10**-10, 10**-9, 10**-8),
                  kvals = c(10**9), rvals = c(0.04),
                  tauvals = c(15, 22.5, 33.75, 50.625, 75.9375))
length(sims2)                  
sims2[[1]]
sims2[[2]]
sims2[[3]]
table(sims2[[1]]$Pop)

## Now that we're sure that everything went well, we'll strat summarizing the data
sub_sims2 <- subset(sims2[[1]], Pop == "B")
class(sub_sims2)
group_sims2 <- dplyr::group_by(sub_sims2, uniq_run, a, b, c, K, tau, r)
group_sims2

sum_sims2 <- dplyr::summarise(group_sims2, maximum_B = max(Density),                
                              maxtime = time[Density == maximum_B],
                              slope = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                           time[time < maxtime & Density < 0.1*K])$coefficients[2],
                              intercept = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                               time[time < maxtime & Density < 0.1*K])$coefficients[1])

for (row in 1:nrow(sum_sims2)) {
  bigfinal_rows <- which(sum_sims2$b[row] == group_sims2$b & 
                           sum_sims2$tau[row] == group_sims2$tau &
                           sum_sims2$a[row] == group_sims2$a &
                           sum_sims2$r[row] == group_sims2$r &
                           sum_sims2$K[row] == group_sims2$K &
                           sum_sims2$c[row] == group_sims2$c)
  
  print(ggplot(data = group_sims2[bigfinal_rows, ],
               aes(x = time, y = Density)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims2$slope[row], intercept = sum_sims2$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims2[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}

sum_sims2

## Let's make the ggplot for this data
ggplot(data = sum_sims2, aes(x = log10(b), y = maxtime, color = as.factor(a),
                             shape = as.factor(K))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r) +
  geom_smooth(method = "lm")

# How to see how the parameters affect maxtime as if they were INDEPENDENT from 
# each other?
reg_sims2 <- lm(maxtime ~ tau + b + a, data = sum_sims2)
summary(reg_sims2)
plot(reg_sims2)

confint(reg_sims2)
# Calulate the RSE
sigma(reg_sims2)/mean(sum_sims2$maxtime) # The lower the RSE, the more accurate the model.
# The RSE estimate gives a measure of error of prediction.

# Analyzing the data as if the parameters were CORRELATED
# INTERACTION EFFECTS
# Build the model
regi_sims2 <- lm(maxtime ~ tau + b + a + tau:b + tau:a + b:a, data = sum_sims2)
                   
summary(regi_sims2)

## We want to plot the multiple linear regressions
# Read data set
sum_sims2
# Create multiple linear regressions
lm_fit <- lm(maxtime ~ log10(b) + log10(a) + log10(tau), data = sum_sims2)
summary(lm_fit)
# Save predictions of the model in the new data frame together with the variable
# you want to plot against
sum_sims2_predicted <- data.frame(maxtime_pred = predict(lm_fit, sum_sims2),
                                  tau = sum_sims2$tau,
                                  b = sum_sims2$b,
                                  a = sum_sims2$a)
sum_sims2_predicted
# This is the predicted line of multiple linear regressions
ggplot(data = sum_sims2, aes(x = log10(tau), y = maxtime, 
                             color = as.factor(log10(a)))) +
  geom_point() +
  geom_line(data = sum_sims2_predicted, aes(x = log10(tau), y = maxtime_pred, 
                                            color = as.factor(log10(a)))) +
  facet_grid(b ~ .)

# Create multiple linear regressions with interactions
lm_fit2 <- lm(maxtime ~ log10(b)*log10(a)*log10(tau), data = sum_sims2)
summary(lm_fit2)
# Save predictions of the model in the new data frame together with the variable
# you want to plot against
sum_sims2_predicted2 <- data.frame(maxtime_pred = predict(lm_fit2, sum_sims2),
                                  tau = sum_sims2$tau,
                                  b = sum_sims2$b,
                                  a = sum_sims2$a)
sum_sims2_predicted2
# This is the predicted line of multiple linear regressions
ggplot(data = sum_sims2, aes(x = log10(tau), y = maxtime, 
                             color = as.factor(log10(a)))) +
  geom_point() +
  geom_line(data = sum_sims2_predicted2, aes(x = log10(tau), y = maxtime_pred, 
                                            color = as.factor(log10(a)))) +
  facet_grid(b ~ .)



## Let's run sims3 ----
## Where we'll have c, K, a, and r constant, and we'll be changing between 
## different values of b, and tau

# Caluclating the b values with R
tau <- c(30, 45, 62, 87, 102)
intercept <- c(7, 16, 23)
slope <- c(0.932, 0.85, 0.715)
bvals <- as.data.frame(matrix(data = NA, ncol = 4, 
                              nrow = length(tau)*length(intercept)*length(slope)))
bvals
i <- 1

for (tauval in tau){
  for (inter in intercept){
    for (slop in slope){
      b <- slop*(tauval-inter)
      bvals[i,] <- c(tauval, inter, slop, b)
      i <- i+1
    }
  }
}

# This step was made to change the name of the axis and the name of the facets
colnames(bvals) <- c("tau", "tradeintercept", "tradeslope", "b")
bvals
class(bvals)

# Now, we caculate it with the names changed
sims3.1 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                  rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
length(sims3.1)
sims3.1[[1]]
sims3.1[[2]]
sims3.1[[3]]

## Let's group_by these simulations
sub_sims3.1 <- subset(sims3.1[[1]], Pop == "B")
class(sub_sims3.1)

group_sims3.1 <- dplyr::group_by(sub_sims3.1, uniq_run, a, b, c, K, tau, r)
group_sims3.1

sum_sims3.1 <- dplyr::summarise(group_sims3.1, maximum_B = max(Density),
                                maxtime = time[Density == maximum_B])
# Here, we cut the slope because we don't need it in our simulation
View(sum_sims3.1)

# This step is useful to combine to data frames that have interesting columns
# that we want to plot together
joined_sims3.1 <- left_join(sum_sims3.1, bvals)
joined_sims3.1

## Let's make the ggplot for this data
ggplot(data = joined_sims3.1, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

# Caluclating the b values with R
tau <- c(15, 18, 21.59999, 25.92, 31.104)
intercept <- c(7, 16, 23)
slope <- c(0.932, 0.85, 0.715)
bvals <- as.data.frame(matrix(data = NA, ncol = 4, 
                              nrow = length(tau)*length(intercept)*length(slope)))
bvals
i <- 1

for (tauval in tau){
  for (inter in intercept){
    for (slop in slope){
      b <- slop*(tauval-inter)
      bvals[i,] <- c(tauval, inter, slop, b)
      i <- i+1
    }
  }
}

# This step was made to change the name of the axis and the name of the facets
colnames(bvals) <- c("tau", "tradeintercept", "tradeslope", "b")
bvals
class(bvals)

bvals <- bvals[-c(4, 5, 6, 7, 8, 9, 16, 17, 18, 25, 26, 27), ]
bvals
# Now, we caculate it with the names changed
sims3.2 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                  rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
length(sims3.2)
sims3.2[[1]]
sims3.2[[2]]
sims3.2[[3]]

## Let's group_by these simulations
sub_sims3.2 <- subset(sims3.2[[1]], Pop == "B")
class(sub_sims3.2)

group_sims3.2 <- dplyr::group_by(sub_sims3.2, uniq_run, a, b, c, K, tau, r)
group_sims3.2

sum_sims3.2 <- dplyr::summarise(group_sims3.2, maximum_B = max(Density),
                              maxtime = time[Density == maximum_B])
# Here, we cut the slope because we don't need it in our simulation
View(sum_sims3.2)

# This step is useful to combine to data frames that have interesting columns
# that we want to plot together
joined_sims3.2 <- left_join(sum_sims3.2, bvals)
joined_sims3.2

## Let's make the ggplot for this data
ggplot(data = joined_sims3.2, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

# Caluclating the b values with R
tau <- c(22, 25.5, 29.5568, 34.259, 39.709)
intercept <- c(7, 16, 23)
slope <- c(0.932, 0.85, 0.715)
bvals <- as.data.frame(matrix(data = NA, ncol = 4, 
                              nrow = length(tau)*length(intercept)*length(slope)))
bvals
i <- 1

for (tauval in tau){
  for (inter in intercept){
    for (slop in slope){
      b <- slop*(tauval-inter)
      bvals[i,] <- c(tauval, inter, slop, b)
      i <- i+1
    }
  }
}

# This step was made to change the name of the axis and the name of the facets
colnames(bvals) <- c("tau", "tradeintercept", "tradeslope", "b")
bvals
class(bvals)

bvals <- bvals[-c(7, 8, 9), ]
bvals
# Now, we caculate it with the names changed
sims3.3 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                  rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
length(sims3.3)
sims3.3[[1]]
sims3.3[[2]]
sims3.3[[3]]

## Let's group_by these simulations
sub_sims3.3 <- subset(sims3.3[[1]], Pop == "B")
class(sub_sims3.3)

group_sims3.3 <- dplyr::group_by(sub_sims3.3, uniq_run, a, b, c, K, tau, r)
group_sims3.3

sum_sims3.3 <- dplyr::summarise(group_sims3.3, maximum_B = max(Density),
                              maxtime = time[Density == maximum_B])
# Here, we cut the slope because we don't need it in our simulation
View(sum_sims3.3)

# This step is useful to combine to data frames that have interesting columns
# that we want to plot together
joined_sims3.3 <- left_join(sum_sims3.3, bvals)
joined_sims3.3

## Let's make the ggplot for this data
ggplot(data = joined_sims3.3, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

## Take all the summarized data frames in section sims2, rbind them, and make a 
## big grap to compere all of them together.

sims3 <- rbind(joined_sims3.1, joined_sims3.2, joined_sims3.3)
sims3
View(sims3)

## Let's make the ggplot for this data
ggplot(data = sims3, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

# If we wanted to plot a regular graph (which we already have in the Word sheet)
# we should follow the steps followed in the other simulations

# Create multiple linear regressions with interactions
lm_fit3 <- lm(maxtime ~ log10(b)*log10(tau)*log10(a), data = sum_sims3)
summary(lm_fit3)
# Save predictions of the model in the new data frame together with the variable
# you want to plot against
sum_sims3_predicted3 <- data.frame(maxtime_pred = predict(lm_fit3, sum_sims3),
                                   tau = sum_sims3$tau,
                                   b = sum_sims3$b,
                                   a = sum_sims3$a)
sum_sims3_predicted3
# This is the predicted line of multiple linear regressions
ggplot(data = sum_sims3, aes(x = log10(b), y = maxtime, 
                             color = as.factor(log10(a)))) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(data = sum_sims3_predicted3, aes(x = log10(b), y = maxtime_pred, 
                                             color = as.factor(log10(a)))) +
  facet_grid(tau ~ .)




## Let's run sims4 ----
## Where we'll have c and K constant, and we'll be changing between 3 different
## values of b, a, r and tau evenly spaced
sims4 <- run_sims(bvals = c(100, 200, 400), rvals = c(0.009, 0.016, 0.02845),
                  avals = c(10**-11, 10**-10, 10**-9), kvals = c(10**9),
                  tauvals = c(22.5, 33.75, 50.625))
length(sims4)                  
sims4[[1]]
sims4[[2]]
sims4[[3]]
table(sims4[[1]]$Pop)

## Now that we're sure that everything went well, we'll strat summarizing the data
sub_sims4 <- subset(sims4[[1]], Pop == "B")
class(sub_sims4)
group_sims4 <- dplyr::group_by(sub_sims4, uniq_run, a, b, c, K, tau, r)
group_sims4

sum_sims4 <- dplyr::summarise(group_sims4, maximum_B = max(Density),                
                              maxtime = time[Density == maximum_B],
                              slope = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                           time[time < maxtime & Density < 0.1*K])$coefficients[2],
                              intercept = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                               time[time < maxtime & Density < 0.1*K])$coefficients[1])

for (row in 1:nrow(sum_sims4)) {
  bigfinal_rows <- which(sum_sims4$b[row] == group_sims4$b & 
                           sum_sims4$tau[row] == group_sims4$tau &
                           sum_sims4$a[row] == group_sims4$a &
                           sum_sims4$r[row] == group_sims4$r &
                           sum_sims4$K[row] == group_sims4$K &
                           sum_sims4$c[row] == group_sims4$c)
  
  print(ggplot(data = group_sims4[bigfinal_rows, ],
               aes(x = time, y = Density)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims4$slope[row], intercept = sum_sims4$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims4[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}

sum_sims4

## Let's make the ggplot for this data
ggplot(data = sum_sims4, aes(x = log10(b), y = maxtime, color = as.factor(a),
                             shape = as.factor(K))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r) +
  geom_smooth(method = "lm")

# How to see how the parameters affect maxtime as if they were INDEPENDENT from 
# each other?
reg_sims4 <- lm(maxtime ~ tau + log10(b) + log10(a) + log10(r), data = sum_sims4)
summary(reg_sims4)

# Calulate the RSE
sigma(reg_sims4)/mean(sum_sims4$maxtime) # The lower the RSE, the more accurate the model.
# The RSE estimate gives a measure of error of prediction.

# Analyzing the data as if the parameters were CORRELATED
# INTERACTION EFFECTS
# Build the model
regi_sims4 <- lm(maxtime ~ tau + log10(b) + log10(a) + log10(r) + tau:log10(b) +
                 tau:log10(a) + tau:log10(r) + log10(b):log10(a) + log10(b):log10(r) +
                   log10(a):log10(r), data = sum_sims4)

summary(regi_sims4)