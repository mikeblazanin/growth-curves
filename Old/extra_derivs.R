#This code is just being kept on hand as simplified versions for
# consistency checking of more complex derivs functions

library(deSolve)
library(ggplot2)
library(dplyr)

#This derivs function only has S and R (no phages), and achieves logistic
# growth through one subpopulation being growing and one being non-growing
derivsx <- function(t, y, parms) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  
  #Create output vector
  dY <- c(S = 0, R = 0)
  
  ##Calculate dS
  #dS/dt = u_S*S - 2*u_S*S*(S+R)/k
  dY["S"] <- (
    parms["u_S"] * y["S"] 
    - 2*parms["u_S"] * y["S"] * (y["S"] + y["R"])/parms["k"])
  
  #Calculate dR
  #dR/dt = 2*u_S*S*(S+R)/k
  dY["R"] <- 2*parms["u_S"] * y["S"] * (y["S"] + y["R"])/parms["k"] 
  
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

#This derivs function achieves logistic growth through two subpopulations,
# one growing and one non-growing. The non-growing is also completely
# resistant to phage infection
derivsy <- function(t, y, parms) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  
  #Create output vector
  dY <- c(S = 0, I = 0, P = 0, R = 0, N = 0)
  
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
  #dS/dt = u_S*S - 2*u_S*S*(k-N)/k - a_t * SP
  dY["S"] <- (
    parms["u_S"] * y["S"] 
    - 2*parms["u_S"] * y["S"] * (parms["k"]-y["N"])/parms["k"] 
    - a_t * y["S"] * y["P"])
  
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
  #dR/dt = 2*u_S*S*(k-N)/k
  dY["R"] <- 2*parms["u_S"] * y["S"] * (parms["k"]-y["N"])/parms["k"]
  
  #Calculate dN
  #dN/dt = - u_S*S
  #        + d*a_tau*S(t-tau)P(t-tau)
  if (t < parms["tau"]) {
    dY["N"] <- -parms["u_S"] * y["S"]
  } else {
    dY["N"] <- -parms["u_S"] * y["S"] +
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

times <- seq(from = 0, to = 4000, by = 1)

yinit <- c("S" = 10**6, "R" = 100)
params <- c(u_S = 0.023,
            k = 10**10,
            warnings = 1, thresh_min_dens = 10**-100)

test <- as.data.frame(dede(y = yinit, times = times, func = derivsx, 
                           parms = params, hmax = 0.1))
test$B <- test$S + test$R
test$pred <- params[["k"]]/
  (1+((params[["k"]] - yinit[["S"]])/yinit[["S"]])*
     exp(-params[["u_S"]]*test$time))
test2 <- tidyr::pivot_longer(test, cols = -c(time), 
                             names_to = "Pop", values_to = "Density")
ggplot(data = filter(test2, Pop %in% c("S", "N", "R", "B")), 
       aes(x = time, y = Density, color = Pop)) +
  geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA)) +
  geom_line(data = filter(test2, Pop == "pred"), lty = 2, color = "black")
tail(test)
print(tail(test$B), digits = 10)
print(tail(test$R), digits = 10)

#Weird, there's an error where B is (S_0 - 10^4)^2 higher than it should be
#no idea why
#Also, why is 2* in the derivs function correct?


times <- seq(from = 0, to = 700, by = 1)
yinit <- c("S" = 10**6, "I" = 0, "P" = 1, "R" = 1000, 
           "N" = (10**9 - 10**6 - 1000))
params <- c(u_S = 0.023,
            k = 10**9,
            a = 5*10**-10,
            tau = 31.6,
            b = 50,
            f = 1,
            d = 0,
            v_a1 = 0, v_a2 = 0,
            z = 0,
            warnings = 1, thresh_min_dens = 10**-100)

test <- as.data.frame(dede(y = yinit, times = times, func = derivsy, 
                           parms = params))
test$B <- test$S + test$R
test$pred <- params[["k"]]/
  (1+((params[["k"]] - yinit[["S"]])/yinit[["S"]])*
     exp(-params[["u_S"]]*test$time))
test2 <- tidyr::pivot_longer(test, cols = -c(time), 
                             names_to = "Pop", values_to = "Density")
ggplot(data = filter(test2, Pop %in% c("S", "I", "P", "N", "R", "B")), 
       aes(x = time, y = Density, color = Pop)) +
  geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA)) +
  geom_line(data = filter(test2, Pop == "pred"), lty = 2, color = "black")
tail(test)

ggplot(data = filter(test2, Pop == "B"), 
       aes(x = time, y = Density)) +
  geom_line() +
  geom_line(data = filter(test2, Pop == "pred"), lty = 2, color = "black")