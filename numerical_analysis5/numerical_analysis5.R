## Import libraries ----
library(deSolve)
library(ggplot2)
library(dplyr)
library(gcplyr)

#Setwd
mywd_split <- strsplit(getwd(), split = "/")
if (mywd_split[[1]][length(mywd_split[[1]])] == "growth-curves") {
  dir.create("numerical_analysis5", showWarnings = FALSE)
  setwd("./numerical_analysis5/")
} else if (mywd_split[[1]][length(mywd_split[[1]])] != "numerical_analysis5") {
  stop("Not in correct root directory")
}

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                      "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)
                      
## Define derivatives functions ----
delay_deriv <- function(t, y, parms) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  if (t >= parms["tau"]) {
    lagY <- c("S1" = lagvalue(t-parms["tau"], 1),
              "S2" = lagvalue(t-parms["tau"], 2),
              "I1" = NA,
              "I2" = NA,
              "P" = lagvalue(t-parms["tau"], 5),
              "N" = lagvalue(t-parms["tau"], 6))
    
    lagY[lagY < parms["thresh_min_dens"]] <- 0
  } else {
    lagY <- c("S1" = 0, "S2" = 0, "I1" = NA, "I2" = NA, "P" = 0, "N" = 0)
  }
  
  #Create output vector
  dY <- c(S1 = 0, S2 = 0, I1 = 0, I2 = 0, P = 0, N = 0)
  
  #Function (f) for reduction of parameters depending on N (rd)
  # frac_t = 1 - f + f*(N/k)
  # frac_tau = 1 - f + f(N(t-tau)/k)
  frd <- function(f, N, k) {return(max(0, 1-f+f*(N/k)))}
  
  af_t <- frd(f = parms["f_a"], N = y["N"], k = parms["k"])
  bf_t <- frd(f = parms["f_b"], N = y["N"], k = parms["k"])
  
  if (t >= parms["tau"]) {
    af_tau <- frd(f = parms["f_a"], N = lagY["N"], k = parms["k"])
  }
  
  #Transitions into resistance could increase with growth rate, decrease with
  # growth rate, or be constant
  #To control  these three options, we have two flags: g1 and g2
  # g1 controls if it's constant (g1 = 0) or variable (g1 = 1) with growth
  # g2 controls if it increases (g2 = 1) or decreases (g2 = -1) with growth
  #h*u_S1*S1*(1 - g1 + g1(g2 * (N/k - 1/2) + 1/2))
  trans <- parms["h"] * parms["u_S1"] * y["S1"] *
    (1 - parms["g1"] + parms["g1"]*(parms["g2"]*(y["N"]/parms["k"]-0.5) + 0.5))
  
  ##Calculate dS1
  #dS1/dt = u_S1*S1*N/k - af_t*a_S1*S1*P - trans
  dY["S1"] <- (
    parms["u_S1"] * y["S1"] * y["N"]/parms["k"]
    - af_t * parms["a_S1"] * y["S1"] * y["P"]
    - trans)
  
  #Calculate dS2
  #dS2/dt = u_S2*S2*N/k + trans - afrac_t*a_S2*S2*P
  dY["S2"] <- 
    (parms["u_S2"] * y["S2"] * y["N"]/parms["k"] 
     + trans
     - af_t * parms["a_S2"] * y["S2"] * y["P"])
  
  ##Calculate dI1
  #dI1/dt = af_t * a_S1 * S1*P 
  #        - af_tau * a_S1 * S1(t-tau) * P(t-tau) 
  if (t < parms["tau"]) {
    dY["I1"] <- af_t * parms["a_S1"] * y["S1"] * y["P"] 
  } else {
    dY["I1"] <- 
      (af_t * parms["a_S1"] * y["S1"] * y["P"] 
       - af_tau * parms["a_S1"] * lagY[1] * lagY[5])
  }
  
  ##Calculate dI2
  #dI2/dt = af_t * a_S2 * S2*P
  #        - af_tau * a_S2 * P(t-tau) * S2(t-tau)
  if (t < parms["tau"]) {
    dY["I2"] <- af_t * parms["a_S2"] * y["S2"] * y["P"] 
  } else {
    dY["I2"] <- 
      (af_t * parms["a_S2"] * y["S2"] * y["P"] 
       - af_tau * parms["a_S2"] * lagY[2] * lagY[5])
  }
  
  ##Calculate dP
  #dP/dt = bf_t * b * af_tau * a_S1 * P(t-tau) * S1(t-tau) 
  #        bf_t * b * af_tau * a_S2 * P(t-tau) * S2(t-tau)
  #        - af_t * a_S1 * S1*P
  #        - af_t * a_S2 * S2*P
  #        - z * af_t * a_S1 * I1 * P
  #        - z * af_t * a_S2 * I2 * P
  #        (factored in code for efficiency)
  if (t < parms["tau"]) {
    dY["P"] <-
      (-af_t * y["P"] *
         (parms["a_S1"] * (y["S1"] + parms["z"]*y["I1"])
          + parms["a_S2"] * (y["S2"] + parms["z"]*y["I2"])))
  } else {
    dY["P"] <-
      (bf_t * parms["b"] * af_tau * lagY[5] *
         (parms["a_S1"]*lagY[1] + parms["a_S2"]*lagY[2])
        - af_t * y["P"] *
          (parms["a_S1"] * (y["S1"] + parms["z"]*y["I1"])
           + parms["a_S2"] * (y["S2"] + parms["z"]*y["I2"])))
  }
  
  #Calculate dN
  #dN/dt = - u_S1*S1 * N/k
  #        - u_S2*S2 * N/k
  #        + d*af_tau * a_S1 * S1(t-tau) * P(t-tau)  
  #        + d*af_tau * a_S2 * S2(t-tau) * P(t-tau) 
  #        (factored in code for efficiency)
  if (t < parms["tau"]) {
    dY["N"] <- 
      (-y["N"]/parms["k"] *
         (parms["u_S1"] * y["S1"] + parms["u_S2"] * y["S2"]))
  } else {
    dY["N"] <- 
      (-y["N"]/parms["k"] *
         (parms["u_S1"] * y["S1"] + parms["u_S2"] * y["S2"])
       + parms["d"]*af_tau*lagY[5]* 
         (parms["a_S1"]*lagY[1] + parms["a_S2"]*lagY[2]))
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


}