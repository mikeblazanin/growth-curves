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
#   ANSWER: Based on what we saw on Tuesday, I’d say that the parameter "a" could
#           change this plot because the three populations depend on it. Plus, "a" 
#           is the probability by which the phage successfully infects the bacteria,
#           which means that, if we increase this parameter, S population will 
#           plateau or increase at some point.
#           
#           Then, we have "K", which is the carrying capacity (maximum population
#           size at which the population can sustain itself). Therefore, this parameter
#           will also change the shape of the plot. The curve will be able to grow
#           more if we increase K, but S population will shrink if we decrease K’s value.
#
#           Also, if we want to change the shape of the susceptible population curve,
#           we could play with "c" (competition coefficient). It measures the strength
#           of competition that does a different specie or, in this case, the infected
#           population, over the susceptible population. Therefore, if "c" is high,
#           the curve will shrink and vice versa.
#    
#           I also think that "r" has something to do with the curve shape. However,
#           I think that "r" depends on the carrying capacity. What I mean is that,
#           if we are in "K" (meaning c = maximum (1?)), "r" will be zero, because
#           competition will be in its maximum, and if we have no competition 
#           (meaning c = 0), "r" will increase. So, the shape of the cure doesn’t
#           depend on a direct way of "r"’s value. Thus, the value that "r" has 
#           is a consequence of the value that "c" has.
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
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtB <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

##Plot results with the new column "B"
library(tidyr)
youtB$B <- youtB$S+youtB$I
youtB_plot <- pivot_longer(youtB, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtB_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

## I'll try to plot the results with Density +10 and the wished colors
#
#Run simulation with "Density + 10"
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 10, K = 10**9,
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
youtD <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))

##Plot results with "Density + 10"
library(tidyr)
youtD$B <- youtD$S+youtD$I
youtD_plot <- pivot_longer(youtD, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = youtD_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

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
BIG <- summarise(bigFINAL, maximum_B = max(youtB1$B, na.rm = T),
                 max(youtB2$B, na.rm = T), max(youtB3$B, na.rm = T),
                 max(youtK$B, na.rm = T), max(youtK2$B, na.rm = T),
                 max(youtR$B, na.rm = T), max(youtR2$B, na.rm = T),
                 max(youtA$B, na.rm = T), max(youtA2$B, na.rm = T),
                 max(youtT$B, na.rm = T), max(youtT2$B, na.rm = T))
summarize(BIG, maximum_B = max(BIG, na.rm = T))

## There's an easier way to do this!!!

bigFINAL <- dplyr::group_by(bigFINAL, b, tau, a, r, K, c)
bigFINAL
FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), maxtime = time[B == maximum_B])
FINAL

## Let's try find the SLOPE of some simulations
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

## Finding the slope with summarize.
FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), maxtime = time[B == maximum_B],
                          slope = mod$coefficients[2])
FINAL
## This wasn't exactly the way to do it. The right way is the following:
## Finding the slope with summarize
FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), 
                maxtime = time[B == maximum_B],
                slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                time[time < maxtime & B < 0.1*K])$coefficients[2],
                 intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                             time[time < maxtime & B < 0.1*K])$coefficients[1])
FINAL

## Here, we try to analyze and plot row by row of the data frame bigFINAL
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

## Let's run the new simulations
##Run simulation with b = 75 & tau = 33
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 75, tau = 33, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout1 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout1)

##Add the columns for each parameter
yout1$r <- 0.04
yout1$a <- 10**-10
yout1$b <- 75
yout1$tau <- 33
yout1$K <- 10**9
yout1$c <- 1
head(yout1)

##Plot results with b = 75 & tau = 33
library(tidyr)
yout1$B <- yout1$S+yout1$I
yout1_plot <- pivot_longer(yout1, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout1_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 75 & tau = 95
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 75, tau = 95, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout2 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout2)

##Add the columns for each parameter
yout2$r <- 0.04
yout2$a <- 10**-10
yout2$b <- 75
yout2$tau <- 95
yout2$K <- 10**9
yout2$c <- 1
head(yout2)

##Plot results with b = 75 & tau = 95
library(tidyr)
yout2$B <- yout2$S+yout2$I
yout2_plot <- pivot_longer(yout2, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout2_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 75 & r = 0.025
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.025, a = 10**-10, b = 75, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout3 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout3)

##Add the columns for each parameter
yout3$r <- 0.025
yout3$a <- 10**-10
yout3$b <- 75
yout3$tau <- 10
yout3$K <- 10**9
yout3$c <- 1
head(yout3)

##Plot results with b = 75 & r = 0.025
library(tidyr)
yout3$B <- yout3$S+yout3$I
yout3_plot <- pivot_longer(yout3, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout3_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 75 & r = 0.009
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.009, a = 10**-10, b = 75, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout4 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout4)

##Add the columns for each parameter
yout4$r <- 0.009
yout4$a <- 10**-10
yout4$b <- 75
yout4$tau <- 10
yout4$K <- 10**9
yout4$c <- 1
head(yout4)

##Plot results with b = 75 & r = 0.009
library(tidyr)
yout4$B <- yout4$S+yout4$I
yout4_plot <- pivot_longer(yout4, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout4_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 75 & a = 10**-10.5
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10.5, b = 75, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout5 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout5)

##Add the columns for each parameter
yout5$r <- 0.04
yout5$a <- 10**-10.5
yout5$b <- 75
yout5$tau <- 10
yout5$K <- 10**9
yout5$c <- 1
head(yout5)

##Plot results with b = 75 & a = 10**-10.5
library(tidyr)
yout5$B <- yout5$S+yout5$I
yout5_plot <- pivot_longer(yout5, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout5_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 75 & a = 10**-8.25
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8.25, b = 75, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout6 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout6)

##Add the columns for each parameter
yout6$r <- 0.04
yout6$a <- 10**-8.25
yout6$b <- 75
yout6$tau <- 10
yout6$K <- 10**9
yout6$c <- 1
head(yout6)

##Plot results with b = 75 & a = 10**-8.25
library(tidyr)
yout6$B <- yout6$S+yout6$I
yout6_plot <- pivot_longer(yout6, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout6_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 75 & K = 10**7.75
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 75, tau = 10, K = 10**7.75, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout7 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout7)

##Add the columns for each parameter
yout7$r <- 0.04
yout7$a <- 10**-10
yout7$b <- 75
yout7$tau <- 10
yout7$K <- 10**7.75
yout7$c <- 1
head(yout7)

##Plot results with b = 75 & K = 10**7.75
library(tidyr)
yout7$B <- yout7$S+yout7$I
yout7_plot <- pivot_longer(yout7, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout7_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 75 & K = 10**8.8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 75, tau = 10, K = 10**8.8, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout8 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout8)

##Add the columns for each parameter
yout8$r <- 0.04
yout8$a <- 10**-10
yout8$b <- 75
yout8$tau <- 10
yout8$K <- 10**8.8
yout8$c <- 1
head(yout8)

##Plot results with b = 75 & K = 10**8.8
library(tidyr)
yout8$B <- yout8$S+yout8$I
yout8_plot <- pivot_longer(yout8, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout8_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 760 & tau = 33
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 760, tau = 33, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout9 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout9)

##Add the columns for each parameter
yout9$r <- 0.04
yout9$a <- 10**-10
yout9$b <- 760
yout9$tau <- 33
yout9$K <- 10**9
yout9$c <- 1
head(yout9)

##Plot results with b = 760 & tau = 33
library(tidyr)
yout9$B <- yout9$S+yout9$I
yout9_plot <- pivot_longer(yout9, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout9_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 760 & tau = 95
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 760, tau = 95, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout10 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout10)

##Add the columns for each parameter
yout10$r <- 0.04
yout10$a <- 10**-10
yout10$b <- 760
yout10$tau <- 95
yout10$K <- 10**9
yout10$c <- 1
head(yout10)

##Plot results with b = 760 & tau = 95
library(tidyr)
yout10$B <- yout10$S+yout10$I
yout10_plot <- pivot_longer(yout10, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout10_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 760 & r = 0.025
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.025, a = 10**-10, b = 760, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout11 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout11)

##Add the columns for each parameter
yout11$r <- 0.025
yout11$a <- 10**-10
yout11$b <- 760
yout11$tau <- 10
yout11$K <- 10**9
yout11$c <- 1
head(yout11)

##Plot results with b = 760 & r = 0.025
library(tidyr)
yout11$B <- yout11$S+yout11$I
yout11_plot <- pivot_longer(yout11, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout11_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 760 & r = 0.009
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.009, a = 10**-10, b = 760, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout12 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout12)

##Add the columns for each parameter
yout12$r <- 0.009
yout12$a <- 10**-10
yout12$b <- 760
yout12$tau <- 10
yout12$K <- 10**9
yout12$c <- 1
head(yout12)

##Plot results with b = 760 & r = 0.009
library(tidyr)
yout12$B <- yout12$S+yout12$I
yout12_plot <- pivot_longer(yout12, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout12_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 760 & a = 10**-10.5
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10.5, b = 760, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout13 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout13)

##Add the columns for each parameter
yout13$r <- 0.04
yout13$a <- 10**-10.5
yout13$b <- 760
yout13$tau <- 10
yout13$K <- 10**9
yout13$c <- 1
head(yout13)

##Plot results with b = 760 & a = 10**-10.5
library(tidyr)
yout13$B <- yout13$S+yout13$I
yout13_plot <- pivot_longer(yout13, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout13_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 760 & a = 10**-8.25
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8.25, b = 760, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout14 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout14)

##Add the columns for each parameter
yout14$r <- 0.04
yout14$a <- 10**-8.25
yout14$b <- 760
yout14$tau <- 10
yout14$K <- 10**9
yout14$c <- 1
head(yout14)

##Plot results with b = 760 & a = 10**-8.25
library(tidyr)
yout14$B <- yout14$S+yout14$I
yout14_plot <- pivot_longer(yout14, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout14_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 760 & K = 10**7.75
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 760, tau = 10, K = 10**7.75, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout15 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout15)

##Add the columns for each parameter
yout15$r <- 0.04
yout15$a <- 10**-10
yout15$b <- 760
yout15$tau <- 10
yout15$K <- 10**7.75
yout15$c <- 1
head(yout15)

##Plot results with b = 760 & K = 10**7.75
library(tidyr)
yout15$B <- yout15$S+yout15$I
yout15_plot <- pivot_longer(yout15, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout15_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with b = 760 & K = 10**8.8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 760, tau = 10, K = 10**8.8, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout16 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout16)

##Add the columns for each parameter
yout16$r <- 0.04
yout16$a <- 10**-10
yout16$b <- 760
yout16$tau <- 10
yout16$K <- 10**8.8
yout16$c <- 1
head(yout16)

##Plot results with b = 760 & K = 10**9.9
library(tidyr)
yout16$B <- yout16$S+yout16$I
yout16_plot <- pivot_longer(yout16, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout16_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 33 & r = 0.025
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.025, a = 10**-10, b = 50, tau = 33, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout17 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout17)

##Add the columns for each parameter
yout17$r <- 0.025
yout17$a <- 10**-10
yout17$b <- 50
yout17$tau <- 33
yout17$K <- 10**9
yout17$c <- 1
head(yout17)

##Plot results with tau = 33 & r = 0.025
library(tidyr)
yout17$B <- yout17$S+yout17$I
yout17_plot <- pivot_longer(yout17, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout17_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 33 & r = 0.009
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.009, a = 10**-10, b = 50, tau = 33, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout18 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout18)

##Add the columns for each parameter
yout18$r <- 0.009
yout18$a <- 10**-10
yout18$b <- 50
yout18$tau <- 33
yout18$K <- 10**9
yout18$c <- 1
head(yout18)

##Plot results with tau = 33 & r = 0.009
library(tidyr)
yout18$B <- yout18$S+yout18$I
yout18_plot <- pivot_longer(yout18, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout18_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 33 & a = 10**-10.5
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10.5, b = 50, tau = 33, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout19 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout19)

##Add the columns for each parameter
yout19$r <- 0.04
yout19$a <- 10**-10.5
yout19$b <- 50
yout19$tau <- 33
yout19$K <- 10**9
yout19$c <- 1
head(yout19)

##Plot results with tau = 33 & a = 10**-10.5
library(tidyr)
yout19$B <- yout19$S+yout19$I
yout19_plot <- pivot_longer(yout19, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout19_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 33 & a = 10**-8.25
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8.25, b = 50, tau = 33, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout20 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout20)

##Add the columns for each parameter
yout20$r <- 0.04
yout20$a <- 10**-8.25
yout20$b <- 50
yout20$tau <- 33
yout20$K <- 10**9
yout20$c <- 1
head(yout20)

##Plot results with tau = 33 & a = 10**-10.5
library(tidyr)
yout20$B <- yout20$S+yout20$I
yout20_plot <- pivot_longer(yout20, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout20_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 33 & K = 10**7.75
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 33, K = 10**7.75, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout21 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout21)

##Add the columns for each parameter
yout21$r <- 0.04
yout21$a <- 10**-10
yout21$b <- 50
yout21$tau <- 33
yout21$K <- 10**7.75
yout21$c <- 1
head(yout21)

##Plot results with tau = 33 & K = 10**7.75
library(tidyr)
yout21$B <- yout21$S+yout21$I
yout21_plot <- pivot_longer(yout21, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout21_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 33 & K = 10**8.8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 33, K = 10**8.8, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout22 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout22)

##Add the columns for each parameter
yout22$r <- 0.04
yout22$a <- 10**-10
yout22$b <- 50
yout22$tau <- 33
yout22$K <- 10**8.8
yout22$c <- 1
head(yout22)

##Plot results with tau = 33 & K = 10**9.9
library(tidyr)
yout22$B <- yout22$S+yout22$I
yout22_plot <- pivot_longer(yout22, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout22_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 95 & r = 0.025
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.025, a = 10**-10, b = 50, tau = 95, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout23 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout23)

##Add the columns for each parameter
yout23$r <- 0.025
yout23$a <- 10**-10
yout23$b <- 50
yout23$tau <- 95
yout23$K <- 10**9
yout23$c <- 1
head(yout23)

##Plot results with tau = 95 & r = 0.025
library(tidyr)
yout23$B <- yout23$S+yout23$I
yout23_plot <- pivot_longer(yout23, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout23_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 95 & r = 0.009
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.009, a = 10**-10, b = 50, tau = 95, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout24 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout24)

##Add the columns for each parameter
yout24$r <- 0.009
yout24$a <- 10**-10
yout24$b <- 50
yout24$tau <- 95
yout24$K <- 10**9
yout24$c <- 1
head(yout24)

##Plot results with tau = 95 & r = 0.009
library(tidyr)
yout24$B <- yout24$S+yout24$I
yout24_plot <- pivot_longer(yout24, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout24_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 95 & a = 10**-10.5
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10.5, b = 50, tau = 95, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout25 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout25)

##Add the columns for each parameter
yout25$r <- 0.04
yout25$a <- 10**-10.5
yout25$b <- 50
yout25$tau <- 95
yout25$K <- 10**9
yout25$c <- 1
head(yout25)

##Plot results with tau = 95 & a = 10**-10.5
library(tidyr)
yout25$B <- yout25$S+yout25$I
yout25_plot <- pivot_longer(yout25, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout25_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 95 & a = 10**-8.25
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8.25, b = 50, tau = 95, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout26 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout26)

##Add the columns for each parameter
yout26$r <- 0.04
yout26$a <- 10**-8.25
yout26$b <- 50
yout26$tau <- 95
yout26$K <- 10**9
yout26$c <- 1
head(yout26)

##Plot results with tau = 95 & a = 10**-10.5
library(tidyr)
yout26$B <- yout26$S+yout26$I
yout26_plot <- pivot_longer(yout26, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout26_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 95 & K = 10**7.75
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 95, K = 10**7.75, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout27 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout27)

##Add the columns for each parameter
yout27$r <- 0.04
yout27$a <- 10**-10
yout27$b <- 50
yout27$tau <- 95
yout27$K <- 10**7.75
yout27$c <- 1
head(yout27)

##Plot results with tau = 95 & K = 10**7.75
library(tidyr)
yout27$B <- yout27$S+yout27$I
yout27_plot <- pivot_longer(yout27, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout27_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with tau = 95 & K = 10**8.8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10, b = 50, tau = 95, K = 10**8.8, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout28 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout28)

##Add the columns for each parameter
yout28$r <- 0.04
yout28$a <- 10**-10
yout28$b <- 50
yout28$tau <- 95
yout28$K <- 10**8.8
yout28$c <- 1
head(yout28)

##Plot results with tau = 95 & K = 10**9.9
library(tidyr)
yout28$B <- yout28$S+yout28$I
yout28_plot <- pivot_longer(yout28, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout28_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with r = 0.025 & a = 10**-10.5
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.025, a = 10**-10.5, b = 50, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout29 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout29)

##Add the columns for each parameter
yout29$r <- 0.025
yout29$a <- 10**-10.5
yout29$b <- 50
yout29$tau <- 10
yout29$K <- 10**9
yout29$c <- 1
head(yout29)

##Plot results with r = 0.025 & a = 10**-10.5
library(tidyr)
yout29$B <- yout29$S+yout29$I
yout29_plot <- pivot_longer(yout29, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout29_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with r = 0.025 & a = 10**-8.25
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.025, a = 10**-8.25, b = 50, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout30 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout30)

##Add the columns for each parameter
yout30$r <- 0.025
yout30$a <- 10**-8.25
yout30$b <- 50
yout30$tau <- 10
yout30$K <- 10**9
yout30$c <- 1
head(yout30)

##Plot results with r = 0.025 & a = 10**-8.25
library(tidyr)
yout30$B <- yout30$S+yout30$I
yout30_plot <- pivot_longer(yout30, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout30_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with r = 0.025 & K = 10**7.75
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.025, a = 10**-10, b = 50, tau = 10, K = 10**7.75, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout31 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout31)

##Add the columns for each parameter
yout31$r <- 0.025
yout31$a <- 10**-10
yout31$b <- 50
yout31$tau <- 10
yout31$K <- 10**7.75
yout31$c <- 1
head(yout31)

##Plot results with r = 0.025 & K = 10**7.75
library(tidyr)
yout31$B <- yout31$S+yout31$I
yout31_plot <- pivot_longer(yout31, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout31_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with r = 0.025 & K = 10**8.8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.025, a = 10**-10, b = 50, tau = 10, K = 10**8.8, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout32 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout32)

##Add the columns for each parameter
yout32$r <- 0.025
yout32$a <- 10**-10
yout32$b <- 50
yout32$tau <- 10
yout32$K <- 10**8.8
yout32$c <- 1
head(yout32)

##Plot results with r = 0.025 & K = 10**8.8
library(tidyr)
yout32$B <- yout32$S+yout32$I
yout32_plot <- pivot_longer(yout32, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout32_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with r = 0.009 & a = 10**-10.5
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.009, a = 10**-10.5, b = 50, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout33 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout33)

##Add the columns for each parameter
yout33$r <- 0.009
yout33$a <- 10**-10.5
yout33$b <- 50
yout33$tau <- 10
yout33$K <- 10**9
yout33$c <- 1
head(yout33)

##Plot results with r = 0.009 & a = 10**-10.5
library(tidyr)
yout33$B <- yout33$S+yout33$I
yout33_plot <- pivot_longer(yout33, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout33_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with r = 0.009 & a = 10**-8.25
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.009, a = 10**-8.25, b = 50, tau = 10, K = 10**9, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout34 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout34)

##Add the columns for each parameter
yout34$r <- 0.009
yout34$a <- 10**-8.25
yout34$b <- 50
yout34$tau <- 10
yout34$K <- 10**9
yout34$c <- 1
head(yout34)

##Plot results with r = 0.009 & a = 10**-8.25
library(tidyr)
yout34$B <- yout34$S+yout34$I
yout34_plot <- pivot_longer(yout34, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout34_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with r = 0.009 & K = 10**7.75
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.009, a = 10**-10, b = 50, tau = 10, K = 10**7.75, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout35 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout35)

##Add the columns for each parameter
yout35$r <- 0.009
yout35$a <- 10**-10
yout35$b <- 50
yout35$tau <- 10
yout35$K <- 10**7.75
yout35$c <- 1
head(yout35)

##Plot results with r = 0.009 & K = 10**7.75
library(tidyr)
yout35$B <- yout35$S+yout35$I
yout35_plot <- pivot_longer(yout35, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout35_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with r = 0.009 & K = 10**8.8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.009, a = 10**-10, b = 50, tau = 10, K = 10**8.8, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout36 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout36)

##Add the columns for each parameter
yout36$r <- 0.009
yout36$a <- 10**-10
yout36$b <- 50
yout36$tau <- 10
yout36$K <- 10**8.8
yout36$c <- 1
head(yout36)

##Plot results with r = 0.009 & K = 10**8.8
library(tidyr)
yout36$B <- yout36$S+yout36$I
yout36_plot <- pivot_longer(yout36, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout36_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with a = 10**-10.5 & K = 10**7.75
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10.5, b = 50, tau = 10, K = 10**7.75, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout37 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout37)

##Add the columns for each parameter
yout37$r <- 0.04
yout37$a <- 10**-10.5
yout37$b <- 50
yout37$tau <- 10
yout37$K <- 10**7.75
yout37$c <- 1
head(yout37)

##Plot results with a = 10**-10.5 & K = 10**7.75
library(tidyr)
yout37$B <- yout37$S+yout37$I
yout37_plot <- pivot_longer(yout37, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout37_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with a = 10**-10.5 & K = 10**8.8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-10.5, b = 50, tau = 10, K = 10**8.8, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout38 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout38)

##Add the columns for each parameter
yout38$r <- 0.04
yout38$a <- 10**-10.5
yout38$b <- 50
yout38$tau <- 10
yout38$K <- 10**8.8
yout38$c <- 1
head(yout38)

##Plot results with a = 10**-10.5 & K = 10**8.8
library(tidyr)
yout38$B <- yout38$S+yout38$I
yout38_plot <- pivot_longer(yout38, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout38_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with a = 10**-8.25 & K = 10**7.75
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8.25, b = 50, tau = 10, K = 10**7.75, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout39 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout39)

##Add the columns for each parameter
yout39$r <- 0.04
yout39$a <- 10**-8.25
yout39$b <- 50
yout39$tau <- 10
yout39$K <- 10**7.75
yout39$c <- 1
head(yout39)

##Plot results with a = 10**-8.25 & K = 10**7.75
library(tidyr)
yout39$B <- yout39$S+yout39$I
yout39_plot <- pivot_longer(yout39, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout39_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

##Run simulation with a = 10**-8.25 & K = 10**8.8
yinit <- c(S = 10**6, I = 0, P = 10**4)
params <- c(r = 0.04, a = 10**-8.25, b = 50, tau = 10, K = 10**8.8, 
            c = 1, warnings = 0, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 250, by = 1)
yout40 <- as.data.frame(dede(y = yinit, times = times, func = derivs, parms = params))
head(yout40)

##Add the columns for each parameter
yout40$r <- 0.04
yout40$a <- 10**-8.25
yout40$b <- 50
yout40$tau <- 10
yout40$K <- 10**8.8
yout40$c <- 1
head(yout40)

##Plot results with a = 10**-8.25 & K = 10**8.8
library(tidyr)
yout40$B <- yout40$S+yout40$I
yout40_plot <- pivot_longer(yout40, c(S, I, P, B), names_to = "Population", values_to = "Density")

ggplot(data = yout40_plot, aes(x = time, y = Density, color = Population)) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = "log10")

## Let's combine all these simulations
youtF <- rbind(yout1, yout2, yout3, yout4, yout5, yout6, yout7, yout8, yout9, yout10,
               yout11, yout12, yout13, yout14, yout15, yout16, yout17, yout18,
               yout19, yout20, yout21, yout22, yout23, yout24, yout25, yout26, yout27,
               yout28, yout29, yout30, yout31, yout32, yout33, yout34, yout35,
               yout36, yout37, yout38, yout39, yout40)
youtF

## Now, I'll try to summarize the data
youtF <- dplyr::group_by(youtF, b, tau, a, r, K, c)
youtF

# Now, find the maximum and the maxtime for each simulation
youtFINAL <- dplyr::summarise(youtF, maximum_B = max(B), maxtime = time[B == maximum_B])
youtFINAL

# Finding the slope with summarize
youtFINAL <- dplyr::summarise(youtF, maximum_B = max(B), 
                          maxtime = time[B == maximum_B],
                          slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                       time[time < maxtime & B < 0.1*K])$coefficients[2],
                          intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                           time[time < maxtime & B < 0.1*K])$coefficients[1])
youtFINAL

## Here, we try to analyze and plot row by row of the data frame youtF
for (row in 1:nrow(youtFINAL)) {
  bigfinal_rows <- which(youtFINAL$b[row] == youtF$b & 
                           youtFINAL$tau[row] == youtF$tau &
                           youtFINAL$a[row] == youtF$a &
                           youtFINAL$r[row] == youtF$r &
                           youtFINAL$K[row] == youtF$K &
                           youtFINAL$c[row] == youtF$c)
  print(ggplot(data = youtF[bigfinal_rows, ],
               aes(x = time, y = B)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = youtFINAL$slope[row], intercept = youtFINAL$intercept[row],
                      color = "red") +
          geom_point(data = youtFINAL[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}  

## To see all the rows of the data frame youtFINAL:
View(youtFINAL)

## How to run the simulations with a LOOP

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

## Plotting all the parameters with maxtime to see how they affect it
ggplot(data = BIG, aes(x = log10(burst), y = maxtime, color = as.factor(K),
                       shape = as.factor(infec))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r)
