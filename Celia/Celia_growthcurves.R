#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

##Celia TODO List
#1. Read code I've put below in derivs. This code should be sufficient to
#   figure out what the model we're using is. Diagram and write out the
#   equations for the model we're simulating. 
#   For the diagram make sure all arrows are labeled.
#   For any symbols you use, make sure to define what the symbol is.
#   (basically I'd like you to do what you did for Problems 2.5 & 3.12
#    but for this model)
#   As some hints:
#   The model includes the following populations:
#     susceptible
#     infected
#     phage
#   The model includes the following symbols:
#     r
#     K
#     c
#     a
#     tau
#     b
#     
#2. Run the simulation with the following conditions:
#       Starting bacterial density = 10**6
#       Starting infected density = 0
#       Starting phage density = 0
#       r = 0.04
#       b = 50
#       a = 10**-10
#       tau = 10
#       c = 1
#       K = 10**9  
#     (hint: I've given you a skeleton code that you can just fill-in)   
#   
#3. Look at how yout is formatted. Use tidyr::pivot_longer to reshape it
#   and ggplot to plot the density of the 3 populations over time
#     hint: the population sizes are quite different, so you'll probably have
#     to use a log10 y-axis to view them. That can be achieved by:
#     + scale_y_continuous(trans = "log10")
#   
#4. Based on the model equations we talked about Tuesday, which parameters
#   do you expect would change this curve?
#
#   ANSWER: Based on what we saw on Tuesday, I’d say that the parameter "a" could change
#           this plot because the three populations depend on it. Plus, "a" is the probability
#           by which the phage successfully infects the bacteria, which means that, if we 
#           increase this parameter, S population will plateau or increase at some point.
#           
#           Then, we have "K", which is the carrying capacity (maximum population size at
#           which the population can sustain itself). Therefore, this parameter will also
#           change the shape of the plot. The curve will be able to grow more if we increase 
#           K, but S population will shrink if we decrease K’s value.
#
#           Also, if we want to change the shape of the susceptible population curve, we could
#           play with "c" (competition coefficient). It measures the strength of competition that
#           does a different specie or, in this case, the infected population, over the susceptible
#           population. Therefore, if "c" is high, the curve will shrink and vice versa.
#    
#           I also think that "r" has something to do with the curve shape. However, I think that
#           "r" depends on the carrying capacity. What I mean is that, if we are in "K" (meaning c
#           = maximum (1?)), "r" will be zero, because competition will be in its maximum, and if 
#           we have no competition (meaning c = 0), "r" will increase. So, the shape of the cure 
#           doesn’t depend on a direct way of "r"’s value. Thus, the value that "r" has is a 
#           consequence of the value that "c" has.
#
#5. Play around with all the parameters and re-run the simulation a few times, 
#   plotting the results each time (you can vary r,b,a,tau,c,K and the starting 
#   bacterial density, but for now keep the starting infected density & 
#   starting phage density at 0). Which parameters affect the density curve?
#   How does each one affect it?
#     note: if you want to simultaneously save the results of multiplt
#     different parameter runs, simply change where you save the results 
#     of dede() to
#     e.g.  yout2 <- as.data.frame(dede(...
#           yout3 <- as.data.frame(dede(...
#   
#   I was wrong!! The value of "r" does change the shape of the curve! If it
#   increases, the curve reaches plateau sooner. 
#
#   I was wrong again, "a" doesn’t affect the shape of the curve!!
#   
#   TO SUM UP, WE ONLY SEE THE PLOT CHANGING WHEN WE VARY THE FOLLOWING
#   PARAMETERS: K, and r.
#
#   Reflection: I think this is because of the conditions we have here.
#   Since P and I = 0, all the parameters that are multiplying them become 0,
#   therefore, they can’t affect the shape of the curve. And, since K and r
#   aren’t multiplying P and I at any point, they don’t become 0, thus are
#   able to change the shape of the cure with the fluctuation of their value.
#   I suppose that if we had different conditions, where P and I weren’t 0,
#   other parameters would make the graph look different.
#
#6. Now run the simulation with the following conditions:
#       Starting bacterial density = 10**6
#       Starting infected density = 0
#       Starting phage density = 10**4
#       r = 0.04
#       b = 50
#       a = 10**-10
#       tau = 10
#       c = 1
#       K = 10**9
#
#       Observations: As we can see, each population has its own starting
#       concentration. We can observe that the infected one increases a lot
#       the first 5 to 8 minutes, and by the 10th it slows down a little bit,
#       even though it's still increasing.
#
#       Phage population doesn’t really grow by time with these conditions,
#       but it maintains its initial concentration, pretty much.
#
#       By these two observations, we can say that phages successfully infect
#       susceptible bacteria and create new phages to maintain their population.
#       However, susceptible population’s rate of growth is higher, because, even
#       though they are infected at some point, they continue to grow and proliferate,
#       more than the phage population.
#
#7. Did both the susceptible and infected populations go extinct?
#   If not, increase the duration of the simulation until they both do go extinct
#   and plot the densities over time.
#   Why do we expect the susceptible and infected populations to always
#   go extinct with this model?
#   (Note that "extinct" can just mean the densities fall below some level,
#   like 10**-5, since the model mathematically allows infinitely small 
#   population sizes)
#
#   ANAWER: None of them went extinct with a duration of 50 (hours?). I tried a
#   longer duration: 250 (hours?) and saw that only the S population went extinct
#   at a time point around 200 (hours?). Then, I tried much longer durations but with
#   none of them I saw I population go extinct!! I was expecting it to go extinct
#   since I thought that with no S population, I oppulation can't be created, thus
#   it's be impossible for them tu survive (if we don't think about resistence to 
#   phages, or sometinh like this).
#
#8. Play around with all the parameters (easiest to leave the starting densities
#   the same for now), plotting the results each time. Are there any
#   parameter combinations that do something unexpected? How do the parameters
#   affect the curve?


## Import libraries ----

library(deSolve)
library(tidyr)
library(ggplot2)
library(dplyr)
library(plyr)

## Define derivatives function ----
derivs <- function(t, y, parms) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Issue warning about too small/negative yvals (if warnings is 1)
  if (parms["warnings"]==1 & any(y < parms["thresh_min_dens"])) {
    warning(paste("pop(s)", paste(which(y < parms["thresh_min_dens"]), collapse = ","), "below thresh_min_dens, treating as 0"))
  }
  
  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  
  #Create output vector
  dY <- c(S = 0, I = 0, P = 0)
  
  ##Calculate dS
  
  #V3 (logistic dS/dt) (including competition from I pop)
  #dS/dt = rS((K-S-c*I)/K) - aSP
  dY["S"] <- parms["r"] * y["S"] * ((parms["K"] - y["S"] - parms["c"] * y["I"])/parms["K"]) - parms["a"] * y["S"] * y["P"]
  
  ##Calculate dI
  #dI/dt = aSP - aS(t-tau)P(t-tau)
  if (t < parms["tau"]) {
    dY["I"] <- parms["a"] * y["S"] * y["P"]
  } else {
    dY["I"] <- parms["a"] * y["S"]*y["P"] - parms["a"] * lagvalue(t - parms["tau"], 1)*lagvalue(t - parms["tau"], 3)
  }
  
  ##Calculate dP
  #dP/dt = baS(t-tau)P(t-tau) - aSP
  if (t < parms["tau"]) {
    dY["P"] <- -parms["a"] * y["S"] * y["P"]
  } else {
    dY["P"] <- parms["b"] * parms["a"] * lagvalue(t-parms["tau"], 1)*lagvalue(t-parms["tau"], 3) - parms["a"]*y["S"]*y["P"]
  }
  
  #Issue warning about too large pop (if warnings is TRUE)
  if (parms["warnings"]==1 & any(y > 10**100)) {
    warning(paste("pop(s)",paste(which(y > 10**100), collapse = ","), "exceed max limit, 10^100, returning dY = 0"))
  }
  dY[y > 10**100] <- 0
  
  #From documentation: The return value of func should be a list, whose first 
  #element is a vector containing the derivatives of y with respect to time
  return(list(dY))
}

#9. Play with the parameters again but knowing the realistic values for each of them:
#  
#   r ranges from 0.04 (a 17-minute doubling time) to 0.007 (a 90-minute doubling time).
#   K ranges from 10^6 to 10^10, although typically we focus on 10^8 to 10^9.
#   a ranges from 10^-12 to 10^-18.
#   tau ranges from 10 to 105 (mins).
#   b ranges from 5 to 1000.
#
#   Also, we'll have to calculate the total bacteria densty by addinc S and I
#   columns into a new one called B. We do thin in "yout", before we pivot_longer it.

## INCORPORATE COLUMN B

##Run simulation with the new column "B"
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtB <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

##Plot results with the new column "B"
library(tidyr)
youtB$B <- youtB$S+youtB$I
youtB_plot <- pivot_longer(youtB, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtB_plot, aes(x = time, y = Density, color = Population)) + geom_line(lwd = 1.5) + scale_y_continuous(trans = "log10")

## I'll try to plot the results with Density +10 and the wished colors
#
#Run simulation with "Density + 10"
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtD <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

##Plot results with "Density + 10"
library(tidyr)
youtD$B <- youtD$S+youtD$I
youtD_plot <- pivot_longer(youtD, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtD_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5) + scale_y_continuous(trans = "log10")

## Now, I'l try to change the colors to the colorblind-friendly scale
#
#Run simulation with different colors
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtC <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

##Plot results with different colors
library(tidyr)
youtC$B <- youtC$S+youtC$I
youtC_plot <- pivot_longer(youtC, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtC_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and K = 10**9
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtK_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and K = 10**8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**8, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtK2_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and r = 0.03
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.03, a = 10**-10, b = 50, tau = 10, K = 10**9, c = 1,  warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtR_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) +  scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and r = 0.007
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.007, a = 10**-10, b = 50, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtR2_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and a = 10**-13
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-13, b = 50, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtA_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and a = 10**-8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8, b = 50, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtA2_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and tau = 65
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 65, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtT_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and tau = 120
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 120, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtT2_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and b = 20
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 20, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtB1_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and b = 500
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 500, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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

ggplot(data = youtB2_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")

#Run simulation with different colors, Density + 10, and b = 850
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 850, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
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
ggplot(data = youtB3_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) + scale_y_continuous(trans = "log10")
dev.off()

## Tiff is the function that converts our graph in an image document, which is
## a ".tiff" one. This is because we saw the legend disrupted here (some problem
## with the code of ggplot). Therefore, I can use this function and close i with
## de.off() to have an imagewhere I CAN clearly see the legend of the plot.


## Finding MAXIMUMS
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

## Now, I'll try to summarizw the data in "bigFINAL"
bigFINAL_plot <- pivot_longer(bigFINAL, c(S, I, P, B), names_to = "Population", values_to = "Density")
ggplot(data = bigFINAL_plot, aes(x = time, y = Density + 10, color = Population)) + geom_point(lwd = 1) + scale_y_continuous(trans = "log10")
## I pivot_longer'd because I remembered that "summarize" works better after doing so
summarise(bigFINAL, maximum_B = max(youtB1$B, na.rm = T), max(youtB2$B, na.rm = T), max(youtB3$B, na.rm = T), max(youtK$B, na.rm = T), max(youtK2$B, na.rm = T), max(youtR$B, na.rm = T), max(youtR2$B, na.rm = T), max(youtA$B, na.rm = T), max(youtA2$B, na.rm = T), max(youtT$B, na.rm = T), max(youtT2$B, na.rm = T))


head(bigFINAL)

## Here, I achived to summarize every maximum of each of the simulations included
## in the big data frame

## Now, I'll try to find the maximum of the maximums
BIG <- summarise(bigFINAL, maximum_B = max(youtB1$B, na.rm = T), max(youtB2$B, na.rm = T), max(youtB3$B, na.rm = T), max(youtK$B, na.rm = T), max(youtK2$B, na.rm = T), max(youtR$B, na.rm = T), max(youtR2$B, na.rm = T), max(youtA$B, na.rm = T), max(youtA2$B, na.rm = T), max(youtT$B, na.rm = T), max(youtT2$B, na.rm = T))
summarize(BIG, maximum_B = max(BIG, na.rm = T))

## There's an easier way to do this!!!

bigFINAL <- dplyr::group_by(bigFINAL, b, tau, a, r, K, c)

FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), maxtime = time[B == maximum_B])
FINAL

head(bigFINAL)
class(bigFINAL)
table(bigFINAL$a)

## Let's try find the SLOPE of some simulations
## I'll try to find the slope for a = 10**-8
## Run simulation
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8, b = 50, tau = 10, K = 10**9, c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtA2S <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

head(youtA2S)

##Plot results
library(tidyr)
youtA2S$B <- youtA2$S+youtA2$I
head(youtA2S)
youtA2S_plot <- pivot_longer(youtA2S, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtA2S_plot, aes(x = time, y = Density + 10, color = Population)) + geom_line(lwd = 1.5, alpha = 1/2) + scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#E69F00")) +scale_y_continuous(trans = "log10")

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

## Finding the slope with summarize.
FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), maxtime = time[B == maximum_B], slope = mod$coefficients[2])
FINAL
