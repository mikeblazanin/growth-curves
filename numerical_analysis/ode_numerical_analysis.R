library(deSolve)
library(tidyr)
library(ggplot2)
library(gcplyr)
library(dplyr)

## Define derivatives function ----
derivs <- function(t, y, parms) {
  #The derivs function must return a list, whose first element is 
  # the derivative of all the variables at a given time
  
  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  
  #Calculate y["I"] = sum(y[I_i])
  y["I"] <- sum(y[2:(parms["nI"]+1)])
  
  #Create output vector
  I_n <- rep(0, parms["nI"])
  names(I_n) <- paste("I", 1:parms["nI"], sep = "")
  dY <- c(S = 0, I_n, P = 0, R = 0)
  
  ##Calculate dS
  #dS/dt = u_S*S((k_S-S-c_SI*I-c_SR*R)/k_S) - aSP
  dY["S"] <- parms["u_S"] * y["S"] * 
    ((parms["k_S"] - y["S"] - 
        parms["c_SI"]*y["I"] - 
        parms["c_SR"]*y["R"])/parms["k_S"]) - 
    parms["a"] * y["S"] * y["P"]
  
  ##Calculate dI
  #dI_1/dt = aSP - nI / tau * I1
  dY["I1"] <- parms["a"]*y["S"]*y["P"] -
    parms["nI"] / parms["tau"] * y["I1"]
  
  #dI_i/dt = nI / tau * I_(i-1) - nI / tau * I_i
  if(parms["nI"] > 1) {
    #Note S pop is dY[1], so I pops are dY[2: nI+1]
    dY[3:(parms["nI"]+1)] <- 
      parms["nI"] / parms["tau"] * y[2:(parms["nI"])] -
               parms["nI"] / parms["tau"] * y[3:(parms["nI"]+1)]
  }
  
  ##Calculate dP
  #dP/dt = b * nI / tau * I_nI - aSP - azP * sum(I)
  dY["P"] <- parms["b"] * parms["nI"] / parms["tau"] * y[(parms["nI"]+1)] -
    parms["a"]*y["S"]*y["P"] -
    parms["a"] * parms["z"] * y["P"] * y["I"]
  
  #Calculate dR
  #dR/dt = u_R*R((k_R-R-c_RS*S-c_RI*I)/k_R) + m*S
  dY["R"] <- parms["u_R"] * y["R"] * 
    ((parms["k_R"] - y["R"] - 
        parms["c_RS"]*y["S"] - 
        parms["c_RI"]*y["I"])/parms["k_R"]) + 
    parms["m"] * y["S"]
  
  #From documentation: The return value of func should be a list, whose first 
  #element is a vector containing the derivatives of y with respect to time
  return(list(dY))
}

#Units are in mins
params <- c(u_S = 0.023,
            u_R = 0,
            k_S = 10**9,
            k_R = 10**9,
            a = 10**-8,
            tau = 500,
            b = 1,
            c_SI = 1,
            c_SR = 1,
            c_RS = 1,
            c_RI = 1,
            z = 0,
            m = 0,
            nI = 100,
            thresh_min_dens = 10**-10)

I_n <- rep(0, params["nI"])
I_n[1] <- 10**6
names(I_n) <- paste("I", 1:params["nI"], sep = "")
yinit <- c(S = 0,
           I_n,
           P = 0,
           R = 0)

times = seq(from = 0, to = 24*60, by = 5)

out <- as.data.frame(
  ode(y = yinit, times = times, func = derivs, parms = params))
out$nI <- 100
out$I <- rowSums(out[, 3:(params["nI"]+2)])

params["nI"] <- 10
I_n <- rep(0, params["nI"])
I_n[1] <- 10**6
names(I_n) <- paste("I", 1:params["nI"], sep = "")
yinit <- c(S = 0,
           I_n,
           P = 0,
           R = 0)
out2 <- as.data.frame(
  ode(y = yinit, times = times, func = derivs, parms = params))
out2$nI <- 10
out2$I <- rowSums(out2[, 3:(params["nI"]+2)])

params["nI"] <- 1
I_n <- rep(0, params["nI"])
I_n[1] <- 10**6
names(I_n) <- paste("I", 1:params["nI"], sep = "")
yinit <- c(S = 0,
           I_n,
           P = 0,
           R = 0)
out3 <- as.data.frame(
  ode(y = yinit, times = times, func = derivs, parms = params))
out3$nI <- 1
out3$I <- out3[, 3:(params["nI"]+2)]

params["nI"] <- 500
I_n <- rep(0, params["nI"])
I_n[1] <- 10**6
names(I_n) <- paste("I", 1:params["nI"], sep = "")
yinit <- c(S = 0,
           I_n,
           P = 0,
           R = 0)
out4 <- as.data.frame(
  ode(y = yinit, times = times, func = derivs, parms = params))
out4$nI <- 500
out4$I <- rowSums(out4[, 3:(params["nI"]+2)])

out <- select(out, time, S, I, P, R, nI)
out2 <- select(out2, time, S, I, P, R, nI)
out3 <- select(out3, time, S, I, P, R, nI)
out4 <- select(out4, time, S, I, P, R, nI)

out_join <- rbind(out, out2, out3, out4)

out_tdy <- pivot_longer(data = out_join, cols = !time & !nI, 
                        names_to = "Pop", values_to = "Dens")

ggplot(data = dplyr::filter(out_tdy, Pop %in% c("S", "I", "P")),
       aes(x = time, y = Dens, color = Pop)) +
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(nI ~ ., scales = "free_y")

out_tdy$deriv <- calc_deriv(y = out_tdy$Dens, x = out_tdy$time,
                      subset_by = paste(out_tdy$Pop, out_tdy$nI))

ggplot(data = dplyr::filter(out_tdy, Pop %in% c("S", "I", "P")),
       aes(x = time, y = deriv, color = Pop)) +
  geom_line() +
  facet_grid(nI ~ ., scales = "free_y")


