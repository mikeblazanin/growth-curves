## Import libraries ----

library(deSolve)
#library(reshape2)
library(data.table)
library(ggplot2)
library(dplyr)
library(gcplyr)

#Setwd
mywd_split <- strsplit(getwd(), split = "/") 
if (mywd_split[[1]][length(mywd_split[[1]])] == "growth-curves") {
  dir.create("numerical_analysis2", showWarnings = FALSE)
  setwd("./numerical_analysis2/")
} else {
  stop("Not in correct root directory")
}

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

##Testing of plasticity in a ----
x <- 10**(seq(from = -3, to = 0, by = 0.01))
testa_combos <- as.data.frame(expand.grid("f" = c(0, 0.5, 1),
                                         "v1" = c(1/16, 1/4, 1, 4, 16),
                                         "v2" = c(1/16, 1/4, 1, 4, 16)))
testa <- data.frame(x = rep(x, times = nrow(testa_combos)),
                   f = rep(testa_combos$f, each = length(x)),
                   v1 = rep(testa_combos$v1, each = length(x)),
                   v2 = rep(testa_combos$v2, each = length(x)),
                   response = NA)
for (i in 1:nrow(testa_combos)) {
  f <- testa_combos$f[i]
  v1 <- testa_combos$v1[i]
  v2 <- testa_combos$v2[i]
  testa$response[testa$f == f & testa$v1 == v1 & testa$v2 == v2] <- 
    (1 - f + f*(1-x)**v1)**v2
}

ggplot(data = testa, 
       aes(x = x, y = response, color = paste(f))) + 
  geom_line() + 
  scale_x_continuous(trans = "log10") +
  #scale_y_continuous(trans = "log10", limits = c(10**-10, 1)) +
  facet_grid(rows = vars(v1), cols = vars(v2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "1-N/k", y = "a/a_max")


## Define derivatives function ----
derivs <- function(t, y, parms) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Issue warning about too small/negative yvals (if warnings is 1)
  if (parms["warnings"]==1 & any(y < parms["thresh_min_dens"])) {
    warning(paste("pop(s)",
                  paste(which(y < parms["thresh_min_dens"]), collapse = ","),
                  "below thresh_min_dens, treating as 0"))
  }
  
  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  
  #Create output vector
  dY <- c(S = 0, I = 0, P = 0, R = 0, N = 0)
  
  #V5
  # added: plasticity in susceptibility to phage infection
  
  #For all equations, let
  # a_t = a * (1 - f + f*(N/k)^v_a1)^v_a2
  # a_tau = a * (1 - f + f*(N(t-tau)/k)^v_a1)^v_a2
  
  a_t <- 
    parms["a"] *
    (1 - parms["f"] +
       parms["f"]*(y["N"]/parms["k"])**parms["v_a1"])**parms["v_a2"]
  if (t < parms["tau"]) {
    a_tau <- 0
  } else {
    a_tau <- 
      parms["a"] *
      (1 - parms["f"] +
         parms["f"] * (lagvalue(t-parms["tau"], 5)/
            parms["k"])**parms["v_a1"])**parms["v_a2"]
  }
  
  ##Calculate dS
  #dS/dt = u_S*S(N/k) - a_t * SP
  dY["S"] <- parms["u_S"] * y["S"] * (y["N"]/parms["k"]) - 
    a_t * y["S"] * y["P"]
  
  ##Calculate dI
  #dI/dt = a_t*SP - a_tau*S(t-tau)P(t-tau)
  if (t < parms["tau"]) {
    dY["I"] <- a_t * y["S"] * y["P"]
  } else {
    dY["I"] <- a_t * y["S"]*y["P"] - 
      a_tau * lagvalue(t - parms["tau"], 1)*lagvalue(t - parms["tau"], 3)
  }
  
  ##Calculate dP
  #dP/dt = b*a_tau*S(t-tau)P(t-tau) - a_t*SP - a_t*zIP
  if (t < parms["tau"]) {
    dY["P"] <- -a_t * y["S"] * y["P"] -
      a_t * parms["z"] * y["I"] * y["P"]
  } else {
    dY["P"] <- 
      parms["b"]*a_tau*lagvalue(t-parms["tau"], 1)*lagvalue(t-parms["tau"], 3) - 
      a_t*y["S"]*y["P"] -
      a_t * parms["z"] * y["I"] * y["P"]
  }
  
  #Calculate dR
  #dR/dt = u_R*R(N/k) + m*S
  dY["R"] <- parms["u_R"] * y["R"] * (y["N"]/parms["k"]) + 
    parms["m"] * y["S"]
  
  #Calculate dN
  #dN/dt = -u_S*S(N/k) - u_R*R(N/k)
  #        + d*a_tau*S(t-tau)P(t-tau)
  if (t < parms["tau"]) {
    dY["N"] <- -parms["u_S"] * y["S"] * (y["N"]/parms["k"]) -
      parms["u_R"] * y["R"] * (y["N"]/parms["k"])
  } else {
    dY["N"] <- -parms["u_S"] * y["S"] * (y["N"]/parms["k"]) -
      parms["u_R"] * y["R"] * (y["N"]/parms["k"]) +
      parms["d"]*a_tau*lagvalue(t-parms["tau"],1)*lagvalue(t-parms["tau"],3)
  }
  
  #Issue warning about too large pop (if warnings is TRUE)
  if (parms["warnings"]==1 & any(y > 10**100)) {
    warning(paste("pop(s)",
                  paste(which(y > 10**100), collapse = ","),
                  "exceed max limit, 10^100, returning dY = 0"))
  }
  dY[y > 10**100] <- 0
  
  #From documentation: The return value of func should be a list, whose first 
  #element is a vector containing the derivatives of y with respect to time
  return(list(dY))
}

yinit <- c("S" = 10**6, "I" = 0, "P" = 10**4, "R" = 0, "N" = 10**9-10**6)
params <- c(u_S = 0.023,
            u_R = 0.023,
            k = 10**9,
            a = 10**-10,
            tau = 31.6,
            b = 50,
            f = 1,
            d = 0,
            v_a1 = 8,
            v_a2 = 1,
            z = 0,
            m = 0,
            warnings = 1, thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 2000, by = 1)
test <- as.data.frame(dede(y = yinit, times = times, func = derivs, 
                           parms = params))
test2 <- tidyr::pivot_longer(test, cols = -c(time), 
                             names_to = "Pop", values_to = "Density")
ggplot(data = test2, aes(x = time, y = Density, color = Pop)) +
  geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA))

#Plots of a_t vs N/k
