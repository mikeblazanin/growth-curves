## Import libraries ----
library(deSolve)
library(ggplot2)
library(dplyr)
library(gcplyr)
library(cowplot)
library(ggh4x)
library(scales)
library(ggtext)
library(tidyr)

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
deriv_dede <- function(t, y, parms) {
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
  #To control  these three options, we have three flags: h, g1 and g2
  # h controls the rate
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
  #        - P * (v + w * (S1 + S2 + I1 + I2))
  #        (factored in code for efficiency)
  if (t < parms["tau"]) {
    dY["P"] <-
      (-af_t * y["P"] *
         (parms["a_S1"] * (y["S1"] + parms["z"]*y["I1"])
          + parms["a_S2"] * (y["S2"] + parms["z"]*y["I2"]))
       -y["P"] * (parms["v"] + parms["w"] * 
                    (y["S1"] + y["S2"] + y["I1"] + y["I2"])))
  } else {
    dY["P"] <-
      (bf_t * parms["b"] * af_tau * lagY[5] *
         (parms["a_S1"]*lagY[1] + parms["a_S2"]*lagY[2])
        - af_t * y["P"] *
          (parms["a_S1"] * (y["S1"] + parms["z"]*y["I1"])
           + parms["a_S2"] * (y["S2"] + parms["z"]*y["I2"]))
       -y["P"] * (parms["v"] + parms["w"] * 
                    (y["S1"] + y["S2"] + y["I1"] + y["I2"])))
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

deriv_ode <- function(t, y, parms) {
  #The derivs function must return a list, whose first element is 
  # the derivative of all the variables at a given time

  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  
  #Create output vector
  dY <- rep(0, length(y))
  names(dY) <- names(y)
  
  #Calculate increased lysis time depending on nutrients
  # (technically it's the ode rate equivalent of the lysis time)
  #lysrt_t = nI/tau * (1 - f_tau + f_tau * N/k)
  lysrt_t <- 
    (parms["nI"] / parms["tau"] *
     max(0, (1 - parms["f_tau"] + parms["f_tau"]*y["N"]/parms["k"])))
  
  ##Calculate dS
  #dS/dt = u_S*S*N/k - aSP
  dY["S1"] <- 
    (parms["u_S1"] * y["S1"] * y["N"]/parms["k"]
     - parms["a_S1"] * y["S1"] * y["P"])
  
  ##Calculate dI
  #dI_1/dt = aSP - I1 * lysrt_t
  dY["I1"] <- 
    (parms["a_S1"] * y["S1"] * y["P"]
     - y["I1"] * lysrt_t)
  
  #dI_i/dt = lysrt_t * I_(i-1) - lysrt_t * I_i
  if(parms["nI"] > 1) {
    #Note S pop is dY[1], so I pops are dY[2: nI+1]
    dY[3:(parms["nI"]+1)] <- 
      (y[2:(parms["nI"])] * lysrt_t
       - y[3:(parms["nI"]+1)] * lysrt_t)
  }
  
  ##Calculate dP
  #dP/dt = b * I_nI * lysrt_t - aSP - z * aP*sum(I)
  dY["P"] <- 
    (parms["b"] * y[(parms["nI"]+1)] * lysrt_t
     - parms["a_S1"] * y["S1"] * y["P"]
     - parms["z"] * parms["a_S1"] * y["P"] * sum(y[2:(parms["nI"]+1)]))
  
  ##Calculate dN
  #dN/dt = -u_S*S*N/k + d * I_nI * nI / tau
  dY["N"] <-
    (-parms["u_S1"] * y["S1"] * y["N"]/parms["k"] +
      parms["d"] * y[(parms["nI"]+1)] * lysrt_t)
  
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

## Simple test run ----
if(F) {
  #delay differential
  times <- seq(from = 0, to = 48*60, by = 15)
  yinit <- c("S1" = 10**6, "S2" = 0, "I1" = 0, "I2" = 0, 
             "P" = 10**4, "N" = (10**9 - 10**6))
  params <- c(u_S1 = 0.0179, u_S2 = 0,
              k = 10**9,
              a_S1 = 10**-8, a_S2 = 0,
              tau = 100, b = 50,
              z = 1,
              f_a = 0, f_b = 0,
              d = 0,
              h = 0, g1 = 0, g2 = 0,
              v = 0, w = 0,
              warnings = 1, thresh_min_dens = 10**-100)
  
  test <- as.data.frame(dede(y = yinit, times = times, func = deriv_dede, 
                                        parms = params))
  test$S <- test$S1 + test$S2
  test$I <- test$I1 + test$I2
  test$B <- test$S + test$I
  test$pred <- params[["k"]]/
    (1+((params[["k"]] - yinit[["S1"]])/yinit[["S1"]])*
       exp(-params[["u_S1"]]*test$time))
  test2 <- tidyr::pivot_longer(test, cols = -c(time), 
                               names_to = "Pop", values_to = "Density")
  ggplot(data = filter(test2, Pop %in% c("S1", "S2", "I1", "I2", "P")), 
         aes(x = time, y = Density, color = Pop)) +
    #geom_line() + 
    #scale_y_continuous(trans = "log10", limits = c(1, NA)) +
    geom_line(data = filter(test2, Pop == "pred"), lty = 2, color = "black") +
    geom_line(data = filter(test2, Pop == "B"), color = "black", alpha = 0.5) +
    NULL
  
  #ode differential
  times <- seq(from = 0, to = 48*60, by = 15)
  yinit <- c("S1" = 10**6, "I1" = 0, "I2" = 0, "I3" = 0, "I4" = 0, 
             "P" = 10**4, "N" = (10**9 - 10**6))
  params <- c(u_S1 = 0.0179,
              k = 10**9,
              a_S1 = 10**-8,
              tau = 100,
              b = 50,
              z = 1,
              f_tau = 0,
              d = 0,
              nI = 4,
              warnings = 1, thresh_min_dens = 10**-100)
  
  test <- as.data.frame(ode(y = yinit, times = times, func = deriv_ode, 
                             parms = params))
  test$I <- test$I1 + test$I2 + test$I3 + test$I4
  test$B <- test$S1 + test$I
  test$pred <- params[["k"]]/
    (1+((params[["k"]] - yinit[["S1"]])/yinit[["S1"]])*
       exp(-params[["u_S1"]]*test$time))
  test2 <- tidyr::pivot_longer(test, cols = -c(time), 
                               names_to = "Pop", values_to = "Density")
  ggplot(data = filter(test2, Pop %in% c("S1", "I1", "I2", "I3", "I4", "P", "N")), 
         aes(x = time, y = Density, color = Pop)) +
    geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA)) +
    #geom_line(data = filter(test2, Pop == "pred"), lty = 2, color = "black") +
    NULL
  
  #ode w/ nI = 50
  times <- seq(from = 0, to = 48*60, by = 15)
  yinit <- rep(0, 53)
  names(yinit) <- c("S1", paste0("I", 1:50), "P", "N")
  yinit["S1"] <- 10**6
  yinit["P"] <- 10**4
  yinit["N"] <- 10**9 - 10**6
  params <- c(u_S1 = 0.0179,
              k = 10**9,
              a_S1 = 10**-8,
              tau = 100,
              b = 50,
              z = 1,
              f_tau = 0,
              d = 0,
              nI = 50,
              warnings = 1, thresh_min_dens = 10**-100)
  
  test <- as.data.frame(ode(y = yinit, times = times, func = deriv_ode, 
                            parms = params))
  test$I <- rowSums(test[, grep("I", colnames(test))])
  test$B <- test$S1 + test$I
  test$pred <- params[["k"]]/
    (1+((params[["k"]] - yinit[["S1"]])/yinit[["S1"]])*
       exp(-params[["u_S1"]]*test$time))
  test2 <- tidyr::pivot_longer(test, cols = -c(time), 
                               names_to = "Pop", values_to = "Density")
  ggplot(data = filter(test2, Pop %in% c("S1", "I", "P", "N")), 
         aes(x = time, y = Density, color = Pop)) +
    geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA)) +
    #geom_line(data = filter(test2, Pop == "pred"), lty = 2, color = "black") +
    NULL
}

## Define function for running simulations across many parameter values ----

##Save a tryCatch function for later use
##  I didn't write this, see: https://stackoverflow.com/a/24569739/14805829
##  notably, returns list with $value $warning and $error
##  successful will have NULL $warning and $error
##  warning will have NULL $error (and value in $value)
##  error will have NULL $value and $warning
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

#sub-function for checking for equilibrium
check_equil <- function(yout_list, cntrs, fixed_time, equil_cutoff_dens,
                        max_j = 10) {
  #Returns: list(keep_running = TRUE/FALSE,
  #              at_equil = TRUE/FALSE/NA (NA only when fixed_time = TRUE),
  #              cntrs = list([other entries],
  #                           I_only_cntr = new I_only_cntr value,
  #                           j = new j value, 
  #                           k = new k value)
  #              )
  
  #If fixed time, just return to never keep running
  if(fixed_time) {
    return(list(keep_running = FALSE, at_equil = NA, cntrs = cntrs))
  }
  
  #Set placeholders
  keep_running <- TRUE
  at_equil <- FALSE
  
  #Infinite loop prevention check (j = 10 is 24 hrs for init_time 100)
  if (cntrs$j >= max_j || cntrs$k >= 15 || cntrs$j+cntrs$k >= 20) {
    keep_running <- FALSE
  }
  
  #If there was an error, increase k by 1, increase history array size and re-run
  if(!is.null(yout_list$error)) {
    cntrs$k <- cntrs$k+1
    cntrs$cntrl_mxhist <- 10**5
    #If there was a warning, could be several causes, so we
    # generally just halve step size and increase length
  } else if (!is.null(yout_list$warning)) {
    cntrs$j <- cntrs$j+1
    cntrs$k <- cntrs$k+2
    #If it was successful, check for equilibrium
  } else if (is.null(yout_list$warning) & is.null(yout_list$error)) {
    #First drop all rows with nan
    yout_list$value <- 
      yout_list$value[apply(X = yout_list$value, MARGIN = 1,
                            FUN = function(x) {all(!is.nan(x))}), ]
    #keep the last row
    temp <- yout_list$value[nrow(yout_list$value), ]
    #keep the second-to-last row
    temp2 <- yout_list$value[(nrow(yout_list$value)-1), ]
    
    #All I are at equil, and either S1 or N are at equil, we're done
    if (all(temp[grep("I", names(temp))] < equil_cutoff_dens) &&
        (temp$S1 < equil_cutoff_dens || temp$N < equil_cutoff_dens)) {
      keep_running <- FALSE
      at_equil <- TRUE
    #S and N not at equil, need more time
    } else if (temp$S1 >= equil_cutoff_dens && temp$N >= equil_cutoff_dens) {
      cntrs$j <- cntrs$j + 1
    #Any I not at equil (but S or N is because above check failed)
    } else if (any(temp[grep("I", names(temp))] >= equil_cutoff_dens)) {
      #If any I are still changing, lengthen
      if(any(abs(temp[grep("I", names(temp))] - temp2[grep("I", names(temp2))]) > 0)) {
        cntrs$j <- cntrs$j+1
      #If none are changing, shorten step size
      } else {
        cntrs$I_only_cntr <- cntrs$I_only_cntr+1
        cntrs$k <- cntrs$k+1
      }
    } else {stop("check_equil found an unexpected case")}
  } else {stop("tryCatch failed, neither success, warning, nor error detected")}
  
  return(list(keep_running = keep_running, at_equil = at_equil, cntrs = cntrs))
}

run_sims_dede <- function(u_S1, u_S2,
                          k,
                          a_S1, a_S2 = NA,
                          tau, b,
                          z = 0,
                          f_a = 0, f_b = 0,
                          d = 1,
                          h = 0, g1 = 0, g2 = 0,
                          v = 0, w = 0,
                          init_S1 = 10**6,
                          init_S2 = 0,
                          init_moi = 10**-2,
                          init_N = NA,
                          equil_cutoff_dens = 0.1,
                          max_time = 48*60,
                          init_time = 12*60,
                          init_stepsize = 15,
                          combinatorial = TRUE,
                          dynamic_stepsize = TRUE,
                          fixed_time = FALSE,
                          print_info = TRUE) {
  #Inputs: vectors of parameters to be combined factorially to make
  #         all possible combinations & run simulations with
  #       equil_cutoff_dens is threshold density to consider a population at equilibrium
  #       init_time is the length of simulation that will first be tried
  #         (and subsequently lengthened as necessary to find equilibrium)
  #       init_stepsize is the length of stepsize that will first be tried
  #         (step sizes scale as simulation legnth increases so that the 
  #         number of timepoints returned for each simulation is constant)
  #       print_info will print about the simulations (e.g. how many will be run)
  #
  #Output: a list with three entries
  #   1. a dataframe (AKA ybig) with all parameters and densities of each pop 
  #       at each timepoint (melted), along with a boolean for whether it'd 
  #       reached equilibrium
  #   2. a dataframe (AKA y_noequil) with just the parameters for runs where
  #       at least one pop didn't reach equilibrium
  #   3. a dataframe (AKA yfail) with the parameters for runs that failed
  #       to successfully complete

  if(fixed_time == TRUE & dynamic_stepsize == TRUE) {
    warning("fixed_time == TRUE overrides dynamic_stepsize == TRUE")
  }
  
  #placeholder for when derivs changes, currently S1 S2 I1 I2 P N S I B PI
  num_pops <- 10 
  
  if(init_time %% init_stepsize != 0) {
    warning("init_time is not divisible by init_stepsize, this has not been tested")
  }
  
  #Save parameter values provided into dataframe
  # taking all different combinations
  sim_vars <- as.list(environment())
  sim_vars <- sim_vars[names(sim_vars) %in%
                         c("u_S1", "u_S2", "k", "a_S1", "a_S2",
                           "tau", "b", "z", "f_a", "f_b", "d",
                           "h", "g1", "g2", "v", "w",
                           "init_S1", "init_S2", "init_moi", "init_N")]
  
  if (combinatorial) {
    param_combos <- expand.grid(sim_vars, stringsAsFactors = FALSE)
    num_sims <- nrow(param_combos)
  } else { #not combinatorial
    num_sims <- max(sapply(X = sim_vars, FUN = length))
    
    #Check for parameter lengths being non-divisible with num_sims
    if (!all(num_sims %% sapply(X = sim_vars, FUN = length) == 0)) {
      warning("Combinatorial=TRUE but longest param vals length is not a multiple of all other param vals lengths")
    }
    
    #Save parameters into dataframe, replicating as needed
    param_combos <- as.data.frame(sapply(X = sim_vars,
                           FUN = function(x) {rep_len(x, num_sims)}))
  }
  
  #if N is NA, N = k-S-R
  param_combos$init_N[is.na(param_combos$init_N)] <-
    (param_combos$k[is.na(param_combos$init_N)] -
       param_combos$init_S1[is.na(param_combos$init_N)] -
       param_combos$init_S2[is.na(param_combos$init_N)])
  #if a_S2 is NA, a_S2 = a_S1
  param_combos$a_S2[is.na(param_combos$a_S2)] <-
    param_combos$a_S1[is.na(param_combos$a_S2)]
  
  #Print number of simulations that will be run
  if(print_info) {
    print(paste(num_sims, "simulations will be run"))
    
    #Save sequence of 10% cutoff points for later reference
    progress_seq <- round(seq(from = 0, to = num_sims, by = num_sims/10))
  }
  
  #Make placeholders
  yfail <- NULL #for runs that fail
  #for runs that succeed, pre-allocate ybig now to save on memory/speed
  ybig <- 
    data.frame(matrix(as.numeric(NA),
                      nrow = num_pops*(1+init_time/init_stepsize)*num_sims,
                      #(adding cols for uniq_run, equil, time, Pop, Density)
                      ncol = 1+length(sim_vars)+4),
               stringsAsFactors = FALSE)
  colnames(ybig) <- c("uniq_run", names(sim_vars),
                      "equil", "time", "Pop", "Density")
  ybig$Pop <- as.character(ybig$Pop)
  
  #Define counters
  cntrs <- list(
    "start_row" = 1, #row where next sim should start being saved
    #number of rows this sim is (constant when dynamic_stepsize = FALSE)
    "this_run_nrows" = num_pops*(1+init_time/init_stepsize),
    "still_needed_toadd" = 0) #num rows needed to add. Minimized redefining ybig
  
  for (i in 1:nrow(param_combos)) { #i acts as the uniq_run counter
    #Define pops & parameters
    yinit <- c(S1 = param_combos$init_S1[i],
               S2 = param_combos$init_S2[i],
               I1 = 0,
               I2 = 0,
               P = param_combos$init_S1[i]*param_combos$init_moi[i],
               N = param_combos$init_N[i])
    params <- c(unlist(param_combos[i, ]),
                warnings = 0, thresh_min_dens = 10**-100)
    
    #Counters
    cntrs["I_only_cntr"] <- 0 #num times S or N at equil but I is not
    cntrs["j"] <- 0 #length counter (larger is longer times)
    cntrs["k"] <- 0 #step size counter (larger is smaller steps)
    cntrs["cntrl_mxhist"] <- 10**4
    
    #Keep running until meets some quit criteria
    while(TRUE) {
      #Define times
      if (dynamic_stepsize) { #double for ea j so length(times) is constant
        times <- seq(0, init_time*2**cntrs$j, init_stepsize*2**cntrs$j)
      } else { #keep stepsize at init_stepsize
        times <- seq(0, init_time*2**cntrs$j, init_stepsize)
      }
      
      #Calculate hmax (max step size integrator uses)
      if (dynamic_stepsize) {hmax_val <- init_stepsize*2**(cntrs$j-cntrs$k) #halved for ea k
      } else {hmax_val <- min(init_stepsize*2**(cntrs$j-cntrs$k), init_stepsize)}
      
      #Run simulation
      yout_list <- myTryCatch(expr = {
        as.data.frame(dede(y = yinit, times = times, func = deriv_dede, 
                           parms = params, hmax = hmax_val,
                           control = list(mxhist = cntrs["cntrl_mxhist"])))
      })
      
      #Check for equil
      checks <- check_equil(yout_list = yout_list, cntrs = cntrs, 
                            fixed_time = fixed_time, 
                            equil_cutoff_dens = equil_cutoff_dens,
                            max_j = ceiling(log2(max_time/init_time)))
      if(checks$keep_running == FALSE) {break
      } else {cntrs <- checks$cntrs}
    }
    
    #Once end conditions triggered, if run succeeded (or warning-d)
    if(!is.null(yout_list$value)) {
      #Calculate S
      yout_list$value$S <- yout_list$value$S1 + yout_list$value$S2
      #Calculate I
      yout_list$value$I <- yout_list$value$I1 + yout_list$value$I2
      #Calculate all bacteria (B)
      yout_list$value$B <- yout_list$value$S + yout_list$value$I
      #Calculate all phage (PI)
      yout_list$value$PI <- yout_list$value$P + yout_list$value$I
      
      if (!dynamic_stepsize) {
        cntrs$this_run_nrows <- num_pops*nrow(yout_list$value)
        cntrs$still_needed_toadd <- 
          cntrs$still_needed_toadd + 
          (cntrs$this_run_nrows - num_pops*(1+init_time/init_stepsize))
        
        #Add rows if we don't have enough room to save this simulation
        if((cntrs$start_row+cntrs$this_run_nrows-1)>nrow(ybig)) {
          ybig[(nrow(ybig)+1):(nrow(ybig)+cntrs$still_needed_toadd), ] <-
            rep(list(rep(NA, cntrs$still_needed_toadd)), ncol(ybig))
          
          cntrs$still_needed_toadd <- 0
        }
      }
      
      #Reshape, add parameters, and fill into ybig in right rows
      myrows <- cntrs$start_row : (cntrs$start_row + cntrs$this_run_nrows - 1)
      ybig[myrows, ] <-
        cbind(uniq_run = i,
              param_combos[i, ],
              equil = checks$at_equil,
              data.table::melt(data = data.table::as.data.table(yout_list$value), 
                               id.vars = c("time"),
                               value.name = "Density", 
                               variable.name = "Pop",
                               variable.factor = FALSE),
              row.names = myrows)
      #If the run failed
    } else if (!is.null(yout_list$error)) {
      temp <- cbind(uniq_run = i, param_combos[i, ], equil = checks$at_equil)
      if (is.null(yfail)) {yfail <- temp
      } else {yfail <- rbind(yfail, temp)}
    } else {stop("tryCatch failed during saving, neither success nor warning nor error detected")}
    
    #Update cumulative offset (for non-dynamic stepsize runs)
    cntrs$start_row <- cntrs$start_row + cntrs$this_run_nrows
    
    #Print progress update
    if (print_info & i %in% progress_seq) {
      print(paste((which(progress_seq == i)-1)*10,
                  "% completed", sep = ""))
    }
  }
  
  #Pull out all the runs that didn't reach equilibrium (just the params)
  y_noequil <- NULL
  for (run in unique(ybig$uniq_run[which(!ybig$equil)])) {
    if (is.null(y_noequil)) {
      y_noequil <- ybig[min(which(ybig$uniq_run == run)),
                        1:(1+length(sim_vars)+1)] #params & cols for uniq_run, equil
    } else {
      y_noequil <- rbind(y_noequil, ybig[min(which(ybig$uniq_run == run)), 
                                         1:(1+length(sim_vars)+1)])
    }
  }
  
  print(paste(nrow(param_combos), " sims run, ",
              length(unique(ybig$uniq_run)), " succeeded (",
              length(unique(y_noequil$uniq_run)), " did not equil), ",
              length(unique(yfail$uniq_run)), " failed",
              sep = ""))
  
  return(list(ybig, y_noequil, yfail))
}

run_sims_ode <- function(u_S1,
                         k,
                         a_S1,
                         tau, b,
                         z = 0,
                         f_tau = 0,
                         d = 1,
                         nI = 1,
                         init_S1 = 10**6,
                         init_moi = 10**-2,
                         init_N = NA,
                         equil_cutoff_dens = 0.1,
                         max_time = 48*60,
                         init_time = 12*60,
                         init_stepsize = 15,
                         combinatorial = TRUE,
                         dynamic_stepsize = TRUE,
                         fixed_time = FALSE,
                         print_info = TRUE) {
  #Inputs: vectors of parameters to be combined factorially to make
  #         all possible combinations & run simulations with
  #       equil_cutoff_dens is threshold density to consider a population at equilibrium
  #       init_time is the length of simulation that will first be tried
  #         (and subsequently lengthened as necessary to find equilibrium)
  #       init_stepsize is the length of stepsize that will first be tried
  #         (step sizes scale as simulation legnth increases so that the 
  #         number of timepoints returned for each simulation is constant)
  #       print_info will print about the simulations (e.g. how many will be run)
  #
  #Output: a list with three entries
  #   1. a dataframe (AKA ybig) with all parameters and densities of each pop 
  #       at each timepoint (melted), along with a boolean for whether it'd 
  #       reached equilibrium
  #   2. a dataframe (AKA y_noequil) with just the parameters for runs where
  #       at least one pop didn't reach equilibrium
  #   3. a dataframe (AKA yfail) with the parameters for runs that failed
  #       to successfully complete
  
  if(fixed_time == TRUE & dynamic_stepsize == TRUE) {
    warning("fixed_time == TRUE overrides dynamic_stepsize == TRUE")
  }
  
  #placeholder for when derivs changes, currently S1 I1 ... In P N I B PI
  num_pops <- 6 + nI
  
  if(init_time %% init_stepsize != 0) {
    warning("init_time is not divisible by init_stepsize, this has not been tested")
  }
  
  #Save parameter values provided into dataframe
  # taking all different combinations
  sim_vars <- as.list(environment())
  sim_vars <- sim_vars[names(sim_vars) %in%
                         c("u_S1", "k", "a_S1",
                           "tau", "b", "z", "f_tau", "d", "nI",
                           "init_S1", "init_S2", "init_moi", "init_N")]
  
  if (combinatorial) {
    param_combos <- expand.grid(sim_vars, stringsAsFactors = FALSE)
    num_sims <- nrow(param_combos)
  } else { #not combinatorial
    num_sims <- max(sapply(X = sim_vars, FUN = length))
    
    #Check for parameter lengths being non-divisible with num_sims
    if (!all(num_sims %% sapply(X = sim_vars, FUN = length) == 0)) {
      warning("Combinatorial=TRUE but longest param vals length is not a multiple of all other param vals lengths")
    }
    
    #Save parameters into dataframe, replicating as needed
    param_combos <- sapply(X = sim_vars,
                           FUN = function(x) {rep_len(x, num_sims)})
  }
  
  #if N is NA, N = k-S-R
  param_combos$init_N[is.na(param_combos$init_N)] <-
    (param_combos$k[is.na(param_combos$init_N)] 
     - param_combos$init_S1[is.na(param_combos$init_N)])
  
  #Print number of simulations that will be run
  if(print_info) {
    print(paste(num_sims, "simulations will be run"))
    
    #Save sequence of 10% cutoff points for later reference
    progress_seq <- round(seq(from = 0, to = num_sims, by = num_sims/10))
  }
  
  #Make placeholders
  yfail <- NULL #for runs that fail
  #for runs that succeed, pre-allocate ybig now to save on memory/speed
  ybig <- 
    data.frame(matrix(as.numeric(NA),
                      nrow = num_pops*(1+init_time/init_stepsize)*num_sims,
                      #(adding cols for uniq_run, equil, time, Pop, Density)
                      ncol = 1+length(sim_vars)+4),
               stringsAsFactors = FALSE)
  colnames(ybig) <- c("uniq_run", names(sim_vars),
                      "equil", "time", "Pop", "Density")
  ybig$Pop <- as.character(ybig$Pop)
  
  #Define counters
  cntrs <- list(
    "start_row" = 1, #row where next sim should start being saved
    #number of rows this sim is (constant when dynamic_stepsize = FALSE)
    "this_run_nrows" = num_pops*(1+init_time/init_stepsize),
    "still_needed_toadd" = 0) #num rows needed to add. Minimized redefining ybig
  
  for (i in 1:nrow(param_combos)) { #i acts as the uniq_run counter
    #Define pops & parameters
    yinit <- c(S1 = param_combos$init_S1[i],
               rep(0, nI),
               P = param_combos$init_S1[i]*param_combos$init_moi[i],
               N = param_combos$init_N[i])
    names(yinit)[2:(2+nI-1)] <- paste0("I", 1:nI)
    params <- c(unlist(param_combos[i, ]),
                warnings = 0, thresh_min_dens = 10**-100)
    
    #Counters
    cntrs["I_only_cntr"] <- 0 #num times S or N at equil but I is not
    cntrs["j"] <- 0 #length counter (larger is longer times)
    cntrs["k"] <- 0 #step size counter (larger is smaller steps)
    cntrs["cntrl_mxhist"] <- 10**4
    
    #Keep running until meets some quit criteria
    while(TRUE) {
      #Define times
      if (dynamic_stepsize) { #double for ea j so length(times) is constant
        times <- seq(0, init_time*2**cntrs$j, init_stepsize*2**cntrs$j)
      } else { #keep stepsize at init_stepsize
        times <- seq(0, init_time*2**cntrs$j, init_stepsize)
      }
      
      #Calculate hmax (max step size integrator uses)
      if (dynamic_stepsize) {hmax_val <- init_stepsize*2**(cntrs$j-cntrs$k) #halved for ea k
      } else {hmax_val <- min(init_stepsize*2**(cntrs$j-cntrs$k), init_stepsize)}
      
      #Run simulation
      yout_list <- myTryCatch(expr = {
        as.data.frame(ode(y = yinit, times = times, func = deriv_ode, 
                           parms = params, hmax = hmax_val))
      })
      
      #Check for equil
      checks <- check_equil(yout_list = yout_list, cntrs = cntrs, 
                            fixed_time = fixed_time, 
                            equil_cutoff_dens = equil_cutoff_dens,
                            max_j = ceiling(log2(max_time/init_time)))
      if(checks$keep_running == FALSE) {break
      } else {cntrs <- checks$cntrs}
    }
    
    #Once end conditions triggered, if run succeeded (or warning-d)
    if(!is.null(yout_list$value)) {
      #Calculate I
      yout_list$value$I <- 
        rowSums(yout_list$value[, grep("I", colnames(yout_list$value))])
      #Calculate all bacteria (B)
      yout_list$value$B <- yout_list$value$S1 + yout_list$value$I
      #Calculate all phage (PI)
      yout_list$value$PI <- yout_list$value$P + yout_list$value$I
      
      if (!dynamic_stepsize) {
        cntrs$this_run_nrows <- num_pops*nrow(yout_list$value)
        cntrs$still_needed_toadd <- 
          cntrs$still_needed_toadd + 
          (cntrs$this_run_nrows - num_pops*(1+init_time/init_stepsize))
        
        #Add rows if we don't have enough room to save this simulation
        if((cntrs$start_row+cntrs$this_run_nrows-1)>nrow(ybig)) {
          ybig[(nrow(ybig)+1):(nrow(ybig)+cntrs$still_needed_toadd), ] <-
            rep(list(rep(NA, cntrs$still_needed_toadd)), ncol(ybig))
          
          cntrs$still_needed_toadd <- 0
        }
      }
      
      #Reshape, add parameters, and fill into ybig in right rows
      myrows <- cntrs$start_row : (cntrs$start_row + cntrs$this_run_nrows - 1)
      ybig[myrows, ] <-
        cbind(uniq_run = i,
              param_combos[i, ],
              equil = checks$at_equil,
              data.table::melt(data = data.table::as.data.table(yout_list$value), 
                               id.vars = c("time"),
                               value.name = "Density", 
                               variable.name = "Pop",
                               variable.factor = FALSE),
              row.names = myrows)
      #If the run failed
    } else if (!is.null(yout_list$error)) {
      temp <- cbind(uniq_run = i, param_combos[i, ], equil = checks$at_equil)
      if (is.null(yfail)) {yfail <- temp
      } else {yfail <- rbind(yfail, temp)}
    } else {stop("tryCatch failed during saving, neither success nor warning nor error detected")}
    
    #Update cumulative offset (for non-dynamic stepsize runs)
    cntrs$start_row <- cntrs$start_row + cntrs$this_run_nrows
    
    #Print progress update
    if (print_info & i %in% progress_seq) {
      print(paste((which(progress_seq == i)-1)*10,
                  "% completed", sep = ""))
    }
  }

  #Pull out all the runs that didn't reach equilibrium (just the params)
  y_noequil <- NULL
  for (run in unique(ybig$uniq_run[which(!ybig$equil)])) {
    if (is.null(y_noequil)) {
      y_noequil <- ybig[min(which(ybig$uniq_run == run)),
                        1:(1+length(sim_vars)+1)] #params & cols for uniq_run, equil
    } else {
      y_noequil <- rbind(y_noequil, ybig[min(which(ybig$uniq_run == run)), 
                                         1:(1+length(sim_vars)+1)])
    }
  }
  
  print(paste(nrow(param_combos), " sims run, ",
              length(unique(ybig$uniq_run)), " succeeded (",
              length(unique(y_noequil$uniq_run)), " did not equil), ",
              length(unique(yfail$uniq_run)), " failed",
              sep = ""))
  
  return(list(ybig, y_noequil, yfail))
}

run_sims_dede(u_S1 = 0.01, u_S2 = 0,
              k = 10**9,
              a_S1 = 10**-10, a_S2 = NA,
              tau = 100, b = 50)

##Define file save/load wrapper for run_sims ----
run_sims_filewrapper <- function(name, mydir = ".",
                                 read_file = TRUE, write_file = TRUE,
                                 type = "dede", ...) {
  #type = c('dede', 'ode')
  
  #Note: ... can be all the arguments to pass to run sims
  #          or it can be multiple lists, each containing the named arguments
  #          to pass to run sims. In the latter case, each list will be
  #          treated as a separate sub-call to run_sims, and the final results
  #          will be all pasted together into one output
  #          (this facilitates non-factorial combinations of parameters)
  
  if (!read_file & !write_file) {
    warning("simulation will be run and no files will be read or written")
  }
  
  #Check to see if results have previously been simulated & saved
  # if so, load them
  if (read_file == TRUE &
      paste(name, "_1.csv", sep = "") %in% list.files(mydir)) {
    warning("Does not check if current param inputs are same as existing data")
    temp <- list(NULL, NULL, NULL)
    temp[[1]] <- read.csv(paste(mydir, "/", name, "_1.csv", sep = ""),
                          stringsAsFactors = F)
    if (paste(name, "_2.csv", sep = "") %in% list.files(mydir)) {
      temp[[2]] <- read.csv(paste(mydir, "/", name, "_2.csv", sep = ""), 
                            stringsAsFactors = F)
    }
    if (paste(name, "_3.csv", sep = "") %in% list.files(mydir)) {
      temp[[3]] <- read.csv(paste(mydir, "/", name, "_3.csv", sep = ""), 
                            stringsAsFactors = F)
    }
  } else {
    #Run simulations (if files don't exist)
    if(all(lapply(list(...), class) == "list")) {
      #Run multiple sub-simulations
      temp_list <- vector(mode = "list", length(list(...))) 
      for(i in 1:length(temp_list)) {
        if(type == "dede") {
          temp_list[[i]] <- do.call(run_sims_dede, list(...)[[i]])
        } else if (type == "ode") {
          temp_list[[i]] <- do.call(run_sims_ode, list(...)[[i]])
        } else {stop("type must be 'dede' or 'ode'")}
        
        #Modify uniq_run #'s as needed
        if(i > 1) {
          temp_list[[i]][[1]][, "uniq_run"] <-
            temp_list[[i]][[1]][, "uniq_run"] +
            suppressWarnings(max(temp_list[[i-1]][[1]][, "uniq_run"],
                                 temp_list[[i-1]][[2]][, "uniq_run"],
                                 temp_list[[i-1]][[3]][, "uniq_run"],
                                 na.rm = TRUE))
        }
      }
      
      #Paste them together
      temp <- list(do.call(rbind, lapply(temp_list, function(x) x[[1]])),
                   do.call(rbind, lapply(temp_list, function(x) x[[2]])),
                   do.call(rbind, lapply(temp_list, function(x) x[[3]])))
    } else {
      if(type == "dede") {
        temp <- run_sims_dede(...)
      } else if (type == "ode") {
        temp <- run_sims_ode(...)
      } else {stop("type must be 'dede' or 'ode'")}
    }
    
    #Save results so they can be re-loaded in future
    if (write_file) {
      #Save results so they can be re-loaded in future
      write.csv(temp[[1]], row.names = F,
                paste(mydir, "/", name, "_1.csv", sep = ""))
      if (!is.null(temp[[2]])) {
        write.csv(temp[[2]], row.names = F,
                  paste(mydir, "/", name, "_2.csv", sep = ""))
      }
      if (!is.null(temp[[3]])) {
        write.csv(temp[[3]], row.names = F,
                  paste(mydir, "/", name, "_3.csv", sep = ""))
      }
    }
  }
  return(temp)
}

## Global Settings ----
glob_read_files <- TRUE
glob_make_curveplots <- FALSE
glob_make_statplots <- TRUE
dir.create("./statplots", showWarnings = FALSE)

## Useful funcs ----
logis_func <- function(S_0, u_S, k, times) {
  return(k/(1+(((k-S_0)/S_0)*exp(-u_S*times))))
}

logis_integral <- function(S_0, u_S, k, times) {
  return(k/u_S *log(k - S_0 + S_0 * exp(u_S * times)))
}

logis_def_integral <- function(S_0, u_S, k, times) {
  return(logis_integral(S_0 = S_0, u_S = u_S, k = k, times = times)-
           logis_integral(S_0 = S_0, u_S = u_S, k = k, times = min(times)))
}

point_slope <- function(...) {
  #Provide x1 and y1, either m or x2 and y2, and either x or y
  #Provide any 4 of 5: y, y1, x, x1, and either m or x2 and y2
  dots <- list(...)
  if(any(!names(dots) %in% c("y", "y1", "x", "x1", "m", "x2", "y2"))) {
    stop("invalid argument specified")}
  if(!"x1" %in% names(dots) || !"y1" %in% names(dots)) {
    stop("x1 and y1 must be provided")}
  if(!"m" %in% names(dots)) {
    if(!"x2" %in% names(dots) || !"y2" %in% names(dots)) {
      stop("when m is not provided, x2 and y2 must be provided")}
    dots[["m"]] <- (dots[["y2"]] - dots[["y1"]])/(dots[["x2"]] - dots[["x1"]])
  }
  if(!"x" %in% names(dots) && !"y" %in% names(dots)) {
    return(dots[["m"]])
  } else {
    if(!"x" %in% names(dots)) {
      return((dots[["y"]] - dots[["y1"]])/dots[["m"]] + dots[["x1"]])
    } else if (!"y" %in% names(dots)) {
      return(dots[["m"]] * (dots[["x"]] - dots[["x1"]]) + dots[["y1"]])
    }
  }
}

interp_data <- function(df, x, y, subset_by) {
  #This function adds new rows such that all the unique x values in
  #  df are present in each subset_by
  #It interpolates y values as necessary
  #x and y are column names
  #subset_by is a vector
  #Note: this function is slow and could be improved by
  # vectorizing the interpolations

  if(length(subset_by) != nrow(df)) {
    stop("subset_by must be the same length as nrow(df)")}
  if(inherits(df, "tbl_df")) {df <- as.data.frame(df)}
  
  alltimes <- unique(df[, x])
  alltimes <- alltimes[order(alltimes)]
  out <- 
    data.frame(matrix(ncol = ncol(df), 
                      nrow = length(unique(subset_by))*length(alltimes)))
  colnames(out) <- colnames(df)
  for (i in 1:length(unique(subset_by))) {
    if(i %% 10 == 0) {print(paste0(i, "/", length(unique(subset_by))))}
    mygroup <- unique(subset_by)[i]
    mysub <- df[subset_by == mygroup, ]
    mysub <- mysub[order(mysub[, x]), ]
    myrows <- ((i-1)*length(alltimes)+1):(i*length(alltimes))
    
    #Fill in id cols
    for(mycol in colnames(out)) {
      if(!mycol %in% c(x, y)) {out[myrows, mycol] <- mysub[1, mycol]}
    }
    
    #fill in (all) x vals
    out[myrows, x] <- alltimes
    
    #fill in values we already have
    myrows2 <- myrows[out[myrows, x] %in% mysub[, x]]
    out[myrows2, y] <- mysub[match(mysub[, x], out[myrows2, x]), y]
    
    submin <- min(mysub[, x])
    submax <- max(mysub[, x])
    
    for(myrow in myrows[is.na(out[myrows, y])]) {
      if(out[myrow, x] > submin && out[myrow, x] < submax) {
        #then interpolate
        idx1 <- max(which(mysub[, x] < out[myrow, x]))
        idx2 <- min(which(mysub[, x] > out[myrow, x]))
        out[myrow, y] <- 
          point_slope(x1 = mysub[idx1, x], x2 = mysub[idx2, x],
                      y1 = mysub[idx1, y], y2 = mysub[idx2, y],
                      x = out[myrow, x])
      } else if (out[myrow, x] > submax) {
        #project last timepoint out to the end
        out[myrow, y] <- mysub[max(which(mysub[, x] < out[myrow, x])), y]
      }
    }
  }
  return(out)
}

central_diff <- function(x, y, end_behavior = "NA"){
  stopifnot(end_behavior %in% c("NA", "back-forwards"))
  res <- pracma::gradient(F = y, h1 = x)
  if(end_behavior == "NA") {res[1] <- NA; res[length(res)] <- NA}
  return(res)
}

## Run 1: phage traits ----
run1 <- run_sims_filewrapper(
  name = "run1",
  u_S1 = signif(0.04*10**-0.35, 3), u_S2 = 0,
  k = 10**9,
  a_S1 = 10**seq(from = -12, to = -8, length.out = 5),
  a_S2 = 0,
  tau = signif(10**seq(from = 1, to = 2, length.out = 5), 3),
  b = signif(5*10**seq(from = 0, to = 2, length.out = 5), 3),
  z = 1,
  d = 0,
  init_S1 = 10**6,
  init_moi = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig1 <- run1[[1]]

#Set below 0 values to 0
ybig1 <- mutate(group_by(ybig1, uniq_run, Pop),
                Density = ifelse(Density < 0, 0, Density),
                deriv = calc_deriv(y = Density, x = time, x_scale = 60),
                percap_deriv = calc_deriv(y = Density, x = time, percapita = TRUE,
                                          blank = 0, window_width_n = 5))

#Main summarization
ysum1_1 <- summarize(group_by(filter(ybig1, Pop == "B"),
                          uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                          tau, b, z, f_a, f_b, d, h, g1, g2,
                          init_S1, init_S2, init_moi, init_N, equil),
                 peak_dens = max(Density),
                 auc = auc(x = time, y = Density),
                 extin_time_4 = 
                   first_below(y = Density, x = time,
                               threshold = 10**4, return = "x"),
                 run_time = max(time),
                 death_slope = min(deriv, na.rm = TRUE),
                 max_percap = max_gc(percap_deriv),
                 first_above_15106 = first_above(y = Density, x = time, return = "x",
                             threshold = 1.5*10**6))
ysum1_2 <- summarize(group_by(filter(ybig1, Pop == "P"),
                              uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                              tau, b, z, f_a, f_b, d, h, g1, g2,
                              init_S1, init_S2, init_moi, init_N, equil),
                     phage_final = Density[which.max(time)])
ysum1_3 <- summarize(group_by(filter(ybig1, Pop %in% c("B", "P")),
                              uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                              tau, b, z, f_a, f_b, d, h, g1, g2,
                              init_S1, init_S2, init_moi, init_N, equil),
                     peak_time = time[which.max(Density[Pop == "B"])],
                     phage_bactpeak = Density[Pop == "P" & time == peak_time])

ysum1 <- full_join(full_join(ysum1_1, ysum1_2), ysum1_3)

#Add phage growth
ysum1 <- mutate(
  ysum1,
  extin_flag = ifelse(is.na(extin_time_4), "noextin",
                      ifelse(peak_dens >= 0.9*k, "neark", "none")),
  extin_time_4 = ifelse(is.na(extin_time_4), run_time, extin_time_4),
  phage_r = (log(phage_final)-log(init_moi*(init_S1+init_S2)))/
    extin_time_4)

#Run PCA
ybig1_PCA <- filter(ybig1, Pop == "B")
ybig1_PCA <- interp_data(df = ybig1_PCA,
                         x = "time", y = "Density",
                         subset_by = ybig1_PCA$uniq_run)
ybig1_PCA <- mutate(group_by(ybig1_PCA, init_S1, u_S1, k),
                    Dens_norm = Density - logis_func(S_0 = init_S1,
                                                     u_S = u_S1,
                                                     k = k,
                                                     times = time))

ybig1_PCA <- filter(ybig1_PCA, time %% 20 == 0)

ybig1_PCA_wide <- tidyr::pivot_wider(ybig1_PCA,
                                     names_from = time,
                                     names_prefix = "t_",
                                     values_from = c(Density, Dens_norm))

mypca <- prcomp(
  ybig1_PCA_wide[, grep("Density_", colnames(ybig1_PCA_wide))[-1]],
  center = TRUE, scale = TRUE, retx = TRUE)
mypcanorm <- prcomp(
  ybig1_PCA_wide[, grep("Dens_norm", colnames(ybig1_PCA_wide))[-1]],
  center = TRUE, scale = TRUE, retx = TRUE)

#Merge PCA with orig data
colnames(mypcanorm$x) <- paste0("norm_", colnames(mypcanorm$x))
ybig1_PCA_wide <- cbind(ybig1_PCA_wide,
                        as.data.frame(mypca$x),
                        as.data.frame(mypcanorm$x))
ybig1_PCA_wide <- inner_join(ybig1_PCA_wide,
                             ysum1)

ysum1 <- left_join(
  ysum1,
  select(ybig1_PCA_wide, 
         !starts_with("Dens") & !(PC6:PC125) & !(norm_PC6:norm_PC125)))

#Add gradients
ysum1 <- mutate(ungroup(ysum1),
                loga = log10(a_S1),
                logb = log10(b),
                logtau = log10(tau),
                norma = as.vector(scale(loga, center = -10, scale = 4)),
                normb = as.vector(scale(logb, center = 1 + log10(5),
                                        scale = 2)),
                normtau = as.vector(scale(logtau, scale = FALSE)),
                logpeakdens = log10(peak_dens),
                peaktimehr = peak_time/60,
                logauc = log10(auc/60),
                extintimehr = extin_time_4/60)

ysum1 <- 
  mutate(group_by(ysum1, logb, logtau),
         across(.cols = c(logpeakdens, peaktimehr, logauc, extintimehr, PC1),
                list("dloga" = ~central_diff(y = .x, x = loga),
                     "dnorma" = ~central_diff(y = .x, x = norma))))
ysum1 <- 
  mutate(group_by(ysum1, loga, logtau),
         across(.cols = c(logpeakdens, peaktimehr, logauc, extintimehr, PC1),
                list("dlogb" = ~central_diff(y = .x, x = logb),
                     "dnormb" = ~central_diff(y = .x, x = normb))))
ysum1 <- 
  mutate(group_by(ysum1, loga, logb),
         across(.cols = c(logpeakdens, peaktimehr, logauc, extintimehr, PC1),
                list("dlogtau" = ~central_diff(y = .x, x = logtau),
                     "dnormtau" = ~central_diff(y = .x, x = normtau))))

# Run 1: example curves for conceptual figure ----
if(glob_make_statplots) {
  myrun <- 3
  dens_offset <- 1
  
  temp <- filter(ybig1, uniq_run == myrun, Pop %in% c("S", "I", "P", "N"))
  png("./statplots/fig1B_run1.png",
      width = 6, height = 4, units = "in", res = 150)
  print(
    ggplot(data = temp, 
           aes(x = time/60, y = Density+dens_offset, color = Pop)) +
      geom_line(lwd = 1.5, alpha = 1) + 
      scale_y_continuous(trans = "log10",
                         breaks = 10**c(0, 3, 6, 9),
                         labels = scales::trans_format("log10", 
                                                       scales::math_format(10^.x))) +
      scale_x_continuous(breaks = seq(from = 0, to = 10, by = 2.5)) +
      scale_color_manual(limits = c("S", "I", "P", "N"),
                         values = my_cols[c(2, 3, 1, 7)],
                         labels = c("Susceptible\nhosts", 
                                    "Infected\nhosts", 
                                    "Phages", "Nutrients")) +
      geom_hline(yintercept = 1, lty = 2) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 20),
            legend.key.spacing.y = unit(0.1, "in"),
            plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "in")) +
      labs(y = "Density", color = "Population", x = "Time (hr)") +
      guides(color = guide_legend(byrow = TRUE)) +
      NULL
  )
  dev.off()
  
  temp <- filter(ybig1, uniq_run == myrun, Pop == "B")
  p1 <- 
    ggplot(data = temp, 
           aes(x = time/60, y = Density+dens_offset)) +
    geom_line(lwd = 1.5, alpha = 1, color = "black") + 
    scale_y_continuous(trans = "log10",
                       breaks = 10**c(0, 3, 6, 9),
                       labels = scales::trans_format("log10", 
                                                     scales::math_format(10^.x))) +
    scale_x_continuous(breaks = seq(from = 0, to = 10, by = 2.5)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "in")) +
    labs(y = "Density", x = "Time (hr)") +
    NULL
  png("./statplots/run1_example_B.png",
      width = 5, height = 5, units = "in", res = 150)
  print(p1)
  dev.off()
  
  png("./statplots/run1_example_Bpeak.png",
      width = 5, height = 5, units = "in", res = 150)
  p1 <- p1 +
    geom_point(data = filter(ysum1, uniq_run == myrun),
               aes(x = peak_time/60, y = peak_dens+dens_offset), 
               color = "red", size = 3) +
    geom_segment(data = filter(ysum1, uniq_run == myrun),
                 aes(y = peak_dens + dens_offset,
                     yend = peak_dens + dens_offset,
                     x = 0, xend = peak_time/60),
                 lty = 2, size = 1.5, color = "red", alpha = 0.8) +
    geom_segment(data = filter(ysum1, uniq_run == myrun),
                 aes(y = 0 + dens_offset,
                     yend = peak_dens + dens_offset,
                     x = peak_time/60, xend = peak_time/60),
                 lty = 2, size = 1.5, color = "red", alpha = 0.8)
  print(p1)
  dev.off()
  
  png("./statplots/run1_example_Bpeakextin.png",
      width = 5, height = 5, units = "in", res = 150)
  p1 <- p1 +
    geom_point(data = filter(ysum1, uniq_run == myrun),
               aes(x = extin_time_4/60, y = 10**4+dens_offset), 
               color = "red", size = 3) +
    geom_segment(data = filter(ysum1, uniq_run == myrun),
                 aes(y = 0 + dens_offset,
                     yend = 10**4 + dens_offset,
                     x = extin_time_4/60, xend = extin_time_4/60),
                 lty = 2, size = 1.5, color = "red", alpha = 0.8)
  print(p1)
  dev.off()
  
  png("./statplots/fig1C_run1.png",
      width = 4, height = 4, units = "in", res = 150)
  p1 <- p1 +
    geom_area(aes(y = Density+dens_offset),
              fill = "red", alpha = 0.5)
  print(p1)
  dev.off()
}


# Run 1: B curves ----
if(glob_make_statplots) {
  f1a <-
    ggplot(data = filter(ybig1, Pop == "B", b == 50, tau == 31.6),
           aes(x = time/60, y = Density)) +
    geom_line(aes(color = as.factor(a_S1), group = interaction(a_S1, b, tau)),
              lwd = 1.5) +
    labs(x = "Time (hr)", y = "Density (cfu/mL)") +
    scale_x_continuous(limits = c(NA, 24), breaks = c(0, 6, 12, 18, 24)) +
    geom_line(data = data.frame(x = 0:1440,
                                y = logis_func(S_0 = 10**6, u_S = 0.0179,
                                               k = 10**9, times = 0:1440)),
              aes(x = x/60, y = y), lty = 2) +
    scale_color_manual(values = colorRampPalette(c("gray70", "darkblue"))(5),
                       name = "Infection rate\n(/cfu/pfu/mL/min)") +
    theme_bw() +
    theme(axis.title = element_text(size = 17),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    NULL
  
  f1b <-
    ggplot(data = filter(ybig1, Pop == "B", a_S1 == 10**-10, tau == 31.6),
           aes(x = time/60, y = Density)) +
    geom_line(aes(color = as.factor(b), group = interaction(a_S1, b, tau)),
              lwd = 1.5) +
    labs(x = "Time (hr)", y = "Density (cfu/mL)") +
    scale_x_continuous(limits = c(NA, 12), breaks = c(0, 6, 12)) +
    geom_line(data = data.frame(x = 0:1440,
                                y = logis_func(S_0 = 10**6, u_S = 0.0179,
                                               k = 10**9, times = 0:1440)),
              aes(x = x/60, y = y), lty = 2) +
    scale_color_manual(values = colorRampPalette(c("gray70", "darkblue"))(5),
                       name = "Burst size") +
    theme_bw() +
    theme(axis.title = element_text(size = 17),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  
  f1c <-
    ggplot(data = filter(ybig1, Pop == "B", a_S1 == 10**-10, b == 50),
           aes(x = time/60, y = Density)) +
    geom_line(aes(color = as.factor(tau), group = interaction(a_S1, b, tau)),
              lwd = 1.5) +
    theme_bw() +
    labs(x = "Time (hr)", y = "Density (cfu/mL)") +
    scale_x_continuous(limits = c(NA, 13), breaks = c(0, 6, 12)) +
    geom_line(data = data.frame(x = 0:1440,
                                y = logis_func(S_0 = 10**6, u_S = 0.0179,
                                               k = 10**9, times = 0:1440)),
              aes(x = x/60, y = y), lty = 2) +
    scale_color_manual(values = colorRampPalette(c("darkblue", "gray70"))(5),
                       name = "Lysis time (min)") +
    theme_bw() +
    theme(axis.title = element_text(size = 17),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  
  png("./statplots/fig1_run1_Bcurves_a_b_tau.png",
      width = 5.25, height = 9, units = "in", res = 150)
  print(plot_grid(f1a, f1b, f1c, ncol = 1, labels = "AUTO",
                  label_size = 20, align = "hv", axis = "tb"))
  dev.off()
}

# Run 1: stat v trait plots ----
if(glob_make_statplots) {
  #a_S1 x-axis (Fig 2 main text)
  f2a <- ggplot(data = ysum1,
         aes(x = log10(a_S1), y = max_percap)) +
    geom_line(aes(group = paste(b, tau)),
              alpha = 0.5,
              position = position_jitter(width = 0.05, height = 0.0001, seed = 1)) +
    scale_y_continuous(limits = c(0, 0.018), breaks = c(0, 0.009, 0.018)) +
    scale_x_continuous(labels = math_format(10^.x)) +
    labs(x = "Infection rate\n(/cfu/pfu/mL/min)", 
         y = "Maximum cellular\ngrowth rate (/min)") +
    #geom_hline(yintercept = 0.0179, lty = 2) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  f2b <- ggplot(data = ysum1,
         aes(x = log10(a_S1), y = first_above_15106)) +
    geom_line(aes(group = paste(b, tau)),
              alpha = 0.5,
              position = position_jitter(width = 0.05, height = 0.2, seed = 1)) +
    ylim(0, NA) +
    scale_x_continuous(labels = math_format(10^.x)) +
    labs(x = "Infection rate\n(/cfu/pfu/mL/min)", 
         y = expression(atop("Time to reach",
                             paste("1.5×", 10^6, " cfu/mL (min)")))) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
    
  f2c <- ggplot(data = ysum1,
                aes(x = log10(a_S1), y = peak_dens/10**8)) +
    geom_line(aes(group = paste(b, tau))) +
    scale_x_continuous(labels = math_format(10^.x)) +
    scale_y_continuous(labels = c("0", 
                                  "2.5×10<sup>8</sup>", 
                                  "5×10<sup>8</sup>",
                                  "7.5×10<sup>8</sup>", 
                                  "10<sup>9</sup>")) +
    labs(x = "Infection rate\n(/cfu/pfu/mL/min)", 
         y = "Peak Bacterial\nDensity (cfu/mL)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          axis.text.y = element_markdown())
  
  f2d <- ggplot(data = ysum1,
         aes(x = log10(a_S1), y = peak_time/60)) +
    geom_line(aes(group = paste(b, tau))) +
    scale_x_continuous(labels = math_format(10^.x)) +
    labs(x = "Infection rate\n(/cfu/pfu/mL/min)", 
         y = "Time of Peak\nBacterial Density (hr)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  f2e <- ggplot(data = ysum1,
         aes(x = log10(a_S1), y = death_slope)) +
    geom_line(aes(group = paste(b, tau))) +
    scale_x_continuous(labels = math_format(10^.x)) +
    scale_y_continuous(breaks = c(0, -5*10**8, -10**9, -1.5*10**9),
                       labels = c(0, -5, -10, -15)) +
    labs(x = "Infection rate\n(/cfu/pfu/mL/min)", 
         y = expression(atop("Maximum rate of",
                             paste("decline (", 10^8, " cfu/hr)")))) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  f2f <- ggplot(data = ysum1,
         aes(x = log10(a_S1), y = extin_time_4/60)) +
    geom_line(aes(group = paste(b, tau))) +
    scale_x_continuous(labels = math_format(10^.x)) +
    labs(x = "Infection rate\n(/cfu/pfu/mL/min)", 
         y = "Extinction Time (hr)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  f2g <- ggplot(data = ysum1,
         aes(x = log10(a_S1), y = auc/60)) +
    geom_line(aes(group = paste(b, tau))) +
    scale_x_continuous(labels = math_format(10^.x)) +
    scale_y_continuous(breaks = 0:4*10**10,
                       labels = c("0",
                                  "1×10<sup>10</sup>",
                                  "2×10<sup>10</sup>",
                                  "3×10<sup>10</sup>",
                                  "4×10<sup>10</sup>")) +
    labs(x = "Infection rate\n(/cfu/pfu/mL/min)", 
         y = "Area Under the\nCurve (hr cfu/mL)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          axis.text.y = element_markdown())
  
  f2h <- ggplot(data = ybig1_PCA_wide,
         aes(x = log10(a_S1), y = PC1)) +
    geom_line(aes(group = paste(b, tau))) +
    scale_x_continuous(labels = math_format(10^.x)) +
    labs(x = "Infection rate\n(/cfu/pfu/mL/min)", 
         y = "PC1") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  png("./statplots/fig2_run1_metricvaS1.png",
      width = 9, height = 14, units = "in", res = 150)
  print(plot_grid(f2a, f2b, f2c, f2d, f2e, f2f, f2g, f2h,
                  ncol = 2, align = 'hv', axis = 'lr',
                  labels = "AUTO", label_size = 16))
  dev.off()
  
  #b x-axis (Fig S2)
  fs2a <- ggplot(data = ysum1,
                aes(x = b, y = max_percap)) +
    geom_line(aes(group = paste(a_S1, tau)),
              alpha = 0.5,
              position = position_jitter(width = 0.05, height = 0.0001, seed = 1)) +
    scale_y_continuous(limits = c(0, 0.018), breaks = c(0, 0.009, 0.018)) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    labs(x = "Burst size", 
         y = "Maximum cellular\ngrowth rate (/min)") +
    #geom_hline(yintercept = 0.0179, lty = 2) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs2b <- ggplot(data = ysum1,
                aes(x = b, y = first_above_15106)) +
    geom_line(aes(group = paste(a_S1, tau)),
              alpha = 0.5,
              position = position_jitter(width = 0.05, height = 0.2, seed = 1)) +
    ylim(0, NA) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    labs(x = "Burst size", 
         y = expression(atop("Time to reach",
                             paste("1.5×", 10^6, " cfu/mL (min)")))) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs2c <- ggplot(data = ysum1,
                aes(x = b, y = peak_dens)) +
    geom_line(aes(group = paste(a_S1, tau))) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_y_continuous(labels = c("0", 
                                  "2.5×10<sup>8</sup>", 
                                  "5×10<sup>8</sup>",
                                  "7.5×10<sup>8</sup>", 
                                  "10<sup>9</sup>")) +
    labs(x = "Burst size",
         y = "Peak Bacterial\nDensity (cfu/mL)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          axis.text.y = element_markdown())
  
  fs2d <- ggplot(data = ysum1,
                aes(x = b, y = peak_time/60)) +
    geom_line(aes(group = paste(a_S1, tau))) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    labs(x = "Burst size", 
         y = "Time of Peak\nBacterial Density (hr)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs2e <- ggplot(data = ysum1,
                aes(x = b, y = death_slope)) +
    geom_line(aes(group = paste(a_S1, tau))) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_y_continuous(breaks = c(0, -5*10**8, -10**9, -1.5*10**9),
                       labels = c(0, -5, -10, -15)) +
    labs(x = "Burst size", 
         y = expression(atop("Maximum rate of",
                             paste("decline (", 10^8, " cfu/hr)")))) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs2f <- ggplot(data = ysum1,
                aes(x = b, y = extin_time_4/60)) +
    geom_line(aes(group = paste(a_S1, tau))) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    labs(x = "Burst size", 
         y = "Extinction Time (hr)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs2g <- ggplot(data = ysum1,
                aes(x = b, y = auc/60)) +
    geom_line(aes(group = paste(a_S1, tau))) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_y_continuous(breaks = 0:4*10**10,
                       labels = c("0",
                                  "1×10<sup>10</sup>",
                                  "2×10<sup>10</sup>",
                                  "3×10<sup>10</sup>",
                                  "4×10<sup>10</sup>")) +
    labs(x = "Burst size", 
         y = "Area Under the\nCurve (hr cfu/mL)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          axis.text.y = element_markdown())
  
  fs2h <- ggplot(data = ybig1_PCA_wide,
                aes(x = b, y = PC1)) +
    geom_line(aes(group = paste(a_S1, tau))) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    labs(x = "Burst size", 
         y = "PC1") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  png("./statplots/figS2_run1_metricvB.png",
      width = 9, height = 14, units = "in", res = 150)
  print(plot_grid(fs2a, fs2b, fs2c, fs2d, fs2e, fs2f, fs2g, fs2h,
                  ncol = 2, align = 'hv', axis = 'lr',
                  labels = "AUTO", label_size = 16))
  dev.off()
  
  #tau x-axis (Fig S3)
  fs3a <- ggplot(data = ysum1,
                 aes(x = tau, y = max_percap)) +
    geom_line(aes(group = paste(a_S1, b)),
              alpha = 0.5,
              position = position_jitter(width = 0.05, height = 0.0001, seed = 1)) +
    scale_y_continuous(limits = c(0, 0.018), breaks = c(0, 0.009, 0.018)) +
    scale_x_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    labs(x = "Lysis time (min)", 
         y = "Maximum cellular\ngrowth rate (/min)") +
    #geom_hline(yintercept = 0.0179, lty = 2) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs3b <- ggplot(data = ysum1,
                 aes(x = tau, y = first_above_15106)) +
    geom_line(aes(group = paste(a_S1, b)),
              alpha = 0.5,
              position = position_jitter(width = 0.05, height = 0.2, seed = 1)) +
    ylim(0, NA) +
    scale_x_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    labs(x = "Lysis time (min)", 
         y = expression(atop("Time to reach",
                             paste("1.5×", 10^6, " cfu/mL (min)")))) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs3c <- ggplot(data = ysum1,
                 aes(x = tau, y = peak_dens)) +
    geom_line(aes(group = paste(a_S1, b))) +
    scale_x_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    scale_y_continuous(labels = c("0", 
                                  "2.5×10<sup>8</sup>", 
                                  "5×10<sup>8</sup>",
                                  "7.5×10<sup>8</sup>", 
                                  "10<sup>9</sup>")) +
    labs(x = "Lysis time (min)",
         y = "Peak Bacterial\nDensity (cfu/mL)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          axis.text.y = element_markdown())
  
  fs3d <- ggplot(data = ysum1,
                 aes(x = tau, y = peak_time/60)) +
    geom_line(aes(group = paste(a_S1, b))) +
    scale_x_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    labs(x = "Lysis time (min)", 
         y = "Time of Peak\nBacterial Density (hr)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs3e <- ggplot(data = ysum1,
                 aes(x = tau, y = death_slope)) +
    geom_line(aes(group = paste(a_S1, b))) +
    scale_x_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    scale_y_continuous(breaks = c(0, -5*10**8, -10**9, -1.5*10**9),
                       labels = c(0, -5, -10, -15)) +
    labs(x = "Lysis time (min)", 
         y = expression(atop("Maximum rate of",
                             paste("decline (", 10^8, " cfu/hr)")))) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs3f <- ggplot(data = ysum1,
                 aes(x = tau, y = extin_time_4/60)) +
    geom_line(aes(group = paste(a_S1, b))) +
    scale_x_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    labs(x = "Lysis time (min)", 
         y = "Extinction Time (hr)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  fs3g <- ggplot(data = ysum1,
                 aes(x = tau, y = auc/60)) +
    geom_line(aes(group = paste(a_S1, b))) +
    scale_x_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    scale_y_continuous(breaks = 0:4*10**10,
                       labels = c("0",
                                  "1×10<sup>10</sup>",
                                  "2×10<sup>10</sup>",
                                  "3×10<sup>10</sup>",
                                  "4×10<sup>10</sup>")) +
    labs(x = "Lysis time (min)", 
         y = "Area Under the\nCurve (hr cfu/mL)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          axis.text.y = element_markdown())
  
  fs3h <- ggplot(data = ybig1_PCA_wide,
                 aes(x = tau, y = PC1)) +
    geom_line(aes(group = paste(a_S1, b))) +
    scale_x_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    labs(x = "Lysis time (min)", 
         y = "PC1") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  
  png("./statplots/figS3_run1_metricvtau.png",
      width = 9, height = 14, units = "in", res = 150)
  print(plot_grid(fs3a, fs3b, fs3c, fs3d, fs3e, fs3f, fs3g, fs3h,
                  ncol = 2, align = 'hv', axis = 'lr',
                  labels = "AUTO", label_size = 16))
  dev.off()
}
  
  
# Run 1: stat v stat plots ----
if(glob_make_statplots) {  
  png("./statplots/fig3_run1_allmetricsvmetrics_subset.png",
      width = 9, height = 9, units = "in", res = 150)
  GGally::ggpairs(
    data = mutate(ungroup(filter(ybig1_PCA_wide, extin_flag != "noextin")),
                  peak_time_hr = peak_time/60,
                  extin_time_4_hr = extin_time_4/60,
                  auc_hr = auc/60),
    columns = c("peak_dens", "peak_time_hr", "extin_time_4_hr", "auc_hr", "PC1"),
    columnLabels = c("Peak Bacterial\nDensity (cfu/mL)",
                     "Time of Peak\nBacterial\nDensity (hr)",
                     "Extinction\nTime (hr)",
                     "Area Under\nthe Curve\n(hr cfu/mL)",
                     "PC1"),
    upper = list(continuous = "points"), lower = list(continuous = "points"),
    diag = list(continuous = "autopointDiag")) +
    theme_bw() +
    theme(strip.text = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
  dev.off()
  
  png("./statplots/figS4_run1_allmetricsvmetrics_alldata.png",
      width = 14, height = 14, units = "in", res = 150)
  fs4 <- GGally::ggpairs(
    data = mutate(ungroup(ybig1_PCA_wide),
                  peak_time_hr = peak_time/60,
                  extin_time_4_hr = extin_time_4/60,
                  auc_hr = auc/60),
    aes(shape = extin_flag),
    columns = c("peak_dens", "peak_time_hr", 
                "extin_time_4_hr", "auc_hr", "PC1",
                "max_percap", "first_above_15106", "death_slope"),
    columnLabels = c("Peak Bacterial\nDensity (cfu/mL)",
                     "Time of Peak\nBacterial\nDensity (hr)",
                     "Extinction\nTime (hr)",
                     "Area Under\nthe Curve\n(hr cfu/mL)",
                     "PC1",
                     "Maximum\ncellular\ngrowth rate\n(/min)",
                     "Time to reach\n1.5×10^6 cfu/mL\n(min)",
                     "Maximum rate\nof decline\n(10^8 cfu/hr)"),
    upper = list(continuous = "points"), lower = list(continuous = "points"),
    diag = list(continuous = "autopointDiag")) +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    theme_bw() +
    theme(strip.text = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
    guides(shape = "none")
  #Modify x and y limits
  for (i in 1:8) {
    fs4[6, i] <- fs4[6, i] + 
      scale_y_continuous(limits = c(0, 0.018), breaks = c(0, 0.009, 0.018))
    fs4[i, 6] <- fs4[i, 6] +
        scale_x_continuous(limits = c(0, 0.018), breaks = c(0, 0.009, 0.018))
    fs4[7, i] <- fs4[7, i] + scale_y_continuous(limits = c(0, NA))
    fs4[i, 7] <- fs4[i, 7] + scale_x_continuous(limits = c(0, NA))
  }
  print(fs4)
  dev.off()

  # maxtime extintime extra plots
  fexa <- ggplot(data = ysum1,
               aes(x = peak_time/60, y = extin_time_4/60 - peak_time/60,
                   color = as.factor(a_S1), shape = extin_flag)) +
          geom_point() +
          scale_color_viridis_d(
            name = "Infection rate\n(/min)", end = 0.85,
            labels = c(expression(10^-12),
                       expression(10^-11), expression(10^-10),
                       expression(10^-9), expression(10^-8))) +
          scale_shape_manual(breaks = c("none", "neark", "noextin"),
                             values = c(16, 4, 3)) +
          scale_y_log10() +
          labs(x = "Peak density (cfu/mL)",
               y = "Extinction time - peak time\n(hr)") +
          guides(shape = "none") +
          theme_bw() +
          theme(axis.title = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 16))
  
  fexb <- ggplot(data = ysum1,
                 aes(x = peak_time/60, y = extin_time_4/peak_time,
                     color = as.factor(a_S1), shape = extin_flag)) +
    geom_point() +
    scale_color_viridis_d(
      name = "Infection rate\n(/min)", end = 0.85,
      labels = c(expression(10^-12),
                 expression(10^-11), expression(10^-10),
                 expression(10^-9), expression(10^-8))) +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Peak density (cfu/mL)",
         y = "Extinction time:Peak time") +
    guides(shape = "none") +
    theme_bw() +
    theme(axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16))
  
  fexc <- ggplot(data = ysum1,
                 aes(x = phage_bactpeak, y = extin_time_4/peak_time,
                     color = as.factor(a_S1), shape = extin_flag)) +
    geom_point() +
    scale_color_viridis_d(
      name = "Infection rate\n(/min)", end = 0.85,
      labels = c(expression(10^-12),
                 expression(10^-11), expression(10^-10),
                 expression(10^-9), expression(10^-8))) +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    scale_x_log10() +
    labs(x = "Phage density when\nbacteria peak (pfu/mL)",
         y = "Extinction time:Peak time") +
    guides(shape = "none") +
    theme_bw() +
    theme(axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16))
  
  png("./statplots/extrafigure_run1_extintime_peaktime_relationships.png",
      width = 15, height = 5, units = "in", res = 150)
  print(plot_grid(fexa + guides(shape = "none", color = "none"),
                  fexb + guides(shape = "none", color = "none"), 
                  fexc, 
                  nrow = 1, labels = "AUTO", rel_widths = c(1, 1, 1.4),
                  label_size = 20, align = "hv", axis = "tb"))
  dev.off()
}

# Run 1: phage growth plots ----
if (glob_make_statplots) {
  f4a <- 
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = phage_r*60, y = log10(peak_dens))) +
    geom_point() +
    scale_x_log10() + 
    scale_y_continuous(labels = math_format(10^.x)) +
    labs(x = "Average phage\n growth rate (/hr)", 
         y = "Peak Bacterial\nDensity (cfu/mL)") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  f4b <- 
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = phage_r*60, y = peak_time/60)) +
    geom_point() +
    scale_x_log10() + 
    scale_y_log10() +
    labs(x = "Average phage\n growth rate (/hr)", 
         y = "Time of Peak\nBacterial Density (hr)") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  f4c <- 
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = phage_r*60, y = extin_time_4/60)) +
    geom_point() +
    scale_x_log10() + 
    scale_y_log10() +
    labs(x = "Average phage\ngrowth rate (/hr)", 
         y = "Extinction time (hr)") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  png("./statplots/extrafigure_run1_phager_extintime_subset_nocol.png", width = 5, height = 4,
      units = "in", res = 300)
  print(f4c)
  dev.off()
  
  f4d <- 
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = phage_r*60, y = log10(auc/60))) +
    geom_point() +
    scale_x_log10() + 
    scale_y_continuous(labels = math_format(10^.x)) +
    labs(x = "Average phage\n growth rate (/hr)", 
         y = "Area Under the\nCurve (hr cfu/mL)") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  f4e <- 
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = phage_r*60, y = PC1)) +
    geom_point() +
    scale_x_log10() + 
    labs(x = "Average phage\n growth rate (/hr)", 
         y = "PC1") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  png("./statplots/fig4_run1_phager_metrics_subset.png", 
      width = 11, height = 7, units = "in", res = 300)
  print(plot_grid(f4a, f4b, f4c, f4d, f4e,
                  nrow = 2, labels = c("A", "B", "C", "D", "E"),
                  align = "hv", axis = "lr", label_size = 20))
  dev.off()
  
  fs5a <- 
    ggplot(data = ysum1,
           aes(x = phage_r*60, y = log10(peak_dens), shape = extin_flag)) +
    geom_point() +
    scale_x_log10() + 
    scale_y_continuous(labels = math_format(10^.x)) +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Average phage\n growth rate (/hr)", 
         y = "Peak Density (cfu/mL)") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  fs5b <- 
    ggplot(data = ysum1,
           aes(x = phage_r*60, y = peak_time/60, shape = extin_flag)) +
    geom_point() +
    scale_x_log10() + 
    scale_y_log10() +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Average phage\n growth rate (/hr)", 
         y = "Peak Time (hr)") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  fs5c <- 
    ggplot(data = ysum1,
           aes(x = phage_r*60, y = extin_time_4/60, shape = extin_flag)) +
    geom_point() +
    scale_x_log10() + 
    scale_y_log10() +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Average phage\ngrowth rate (/hr)", 
         y = "Extinction time (hr)") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  fs5d <- 
    ggplot(data = ysum1,
           aes(x = phage_r*60, y = log10(auc/60), shape = extin_flag)) +
    geom_point() +
    scale_x_log10() + 
    scale_y_continuous(labels = math_format(10^.x)) +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Average phage\n growth rate (/hr)", 
         y = "Area Under the\nCurve (hr cfu/mL)") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  fs5e <- 
    ggplot(data = ysum1,
           aes(x = phage_r*60, y = PC1, shape = extin_flag)) +
    geom_point() +
    scale_x_log10() + 
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Average phage\n growth rate (/hr)", 
         y = "PC1") +
    theme_bw() +
    guides(shape = "none") +
    theme(axis.title = element_text(size = 16))
  
  png("./statplots/figS5_run1_phager_metrics_alldata.png", 
      width = 11, height = 7, units = "in", res = 300)
  print(plot_grid(fs5a, fs5b, fs5c, fs5d, fs5e,
                  nrow = 2, labels = c("A", "B", "C", "D", "E"),
                  align = "hv", axis = "lr", label_size = 20))
  dev.off()
  
  png("./statplots/extrafigure_run1_finalphagevpeakdens_alldata.png", 
      width = 5.5, height = 3.5, units = "in", res = 300)
  print(
    ggplot(data = ysum1,
           aes(x = peak_dens, y = phage_final, shape = extin_flag,
               color = as.factor(b))) +
      geom_point(size = 2) +
      scale_y_log10() + scale_x_log10() +
      scale_color_viridis_d(end = 0.95, name = "Burst size") +
      scale_shape_manual(breaks = c("none", "neark", "noextin"),
                         values = c(16, 4, 3)) +
      labs(x = "Peak bacterial density (cfu/mL)", 
           y = "Final phage density\n(pfu/mL)") +
      guides(shape = "none") +
      geom_line(aes(y = peak_dens*b)) +
      theme_bw() +
      theme(axis.title = element_text(size = 20),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 14))
  )
  dev.off()
}

#Run 1: contour plots ----
if (glob_make_statplots) {
  f5a <- ggplot(data = filter(ysum1, b == 50),
                aes(x = log10(a_S1), y = tau)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 3) +
    scale_color_viridis_c(name = "Time of\nPeak\nBacterial\nDensity (hr)",
                          limits = c(0, 18), breaks = c(0, 6, 12, 18)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    scale_x_continuous(labels = math_format(10^.x)) +
    xlab("Infection rate\n(/cfu/pfu/mL/min)") +
    ylab("Lysis time (min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 20),
          legend.title = element_text(size = 14, 
                                      margin = margin(0, 0, 0.07, 0, unit = "npc")),
          legend.text = element_text(size = 13)) +
    NULL
  
  f5b <- ggplot(data = filter(ysum1, tau == 31.6),
                aes(x = log10(a_S1), y = b)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 3) +
    scale_color_viridis_c(name = "Time of\nPeak\nBacterial\nDensity (hr)",
                          limits = c(0, 18), breaks = c(0, 6, 12, 18)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(labels = math_format(10^.x)) +
    xlab("Infection rate\n(/cfu/pfu/mL/min)") +
    ylab("Burst size") +
    guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 20),
          legend.title = element_text(size = 14, 
                                      margin = margin(0, 0, 0.07, 0, unit = "npc")),
          legend.text = element_text(size = 13)) +
    NULL
  
  f5c <- ggplot(data = filter(ysum1, a_S1 == 10**-10),
                aes(x = b, y = tau)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 3) +
    scale_color_viridis_c(name = "Time of\nPeak\nBacterial\nDensity (hr)",
                          limits = c(0, 18), breaks = c(0, 6, 12, 18)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_x_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_y_continuous(trans = "log10", breaks = c(10, 32, 100)) +
    xlab("Burst size") +
    ylab("Lysis time (min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 20),
          legend.title = element_text(size = 14, 
                                      margin = margin(0, 0, 0.07, 0, unit = "npc")),
          legend.text = element_text(size = 13)) +
    NULL
  
  f5d <- ggplot(data = pivot_longer(filter(ysum1, extin_flag == "none"),
                             cols = starts_with("peaktimehr_dlog"),
                             names_to = "wrt",
                             names_prefix = "peaktimehr_dlog",
                             values_to = "derivative"),
         aes(x = wrt, y = abs(derivative))) +
    geom_point(position = position_jitter(width = 0.1, seed = 1), alpha = 0.5) +
    labs(y = "Magnitude of derivative\nof peak time (hr/10-fold change)", 
         x = "With respect to") +
    scale_x_discrete(breaks = c("a", "b", "tau"),
                     labels = c("Infection rate", "Burst size", "Lysis time")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  f5e <- ggplot(data = pivot_longer(filter(ysum1, extin_flag == "none"),
                                    cols = starts_with("peaktimehr_dnorm"),
                                    names_to = "wrt",
                                    names_prefix = "peaktimehr_dnorm",
                                    values_to = "derivative"),
                aes(x = wrt, y = abs(derivative))) +
    geom_point(position = position_jitter(width = 0.1, seed = 1), alpha = 0.5) +
    labs(y = "Magnitude of derivative\nof peak time (hr/full range)", 
         x = "With respect to") +
    scale_x_discrete(breaks = c("a", "b", "tau"),
                     labels = c("Infection rate", "Burst size", "Lysis time")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png("./statplots/fig5_run1_contours_derivs.png", width = 10, height = 7.5,
      units = "in", res = 300)
  cowplot::plot_grid(f5a, f5b, f5c,
                     cowplot::plot_grid(f5d, f5e, nrow = 1, labels = c("D", "E")),
                     nrow = 2, labels = c("A", "B", "C", ""),
                     rel_heights = c(1, 0.9))
  dev.off()
  
  ##Supplemental contours
  #Peak time
  p1 <- ggplot(data = ysum1, aes(x = a_S1, y = tau)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 1) +
    facet_grid(~b) +
    scale_color_viridis_c(name = "Peak time (hr)",
                          breaks = c(4, 8, 12, 16)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Lysis time (min)",
         subtitle = "Burst Size") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(data = ysum1, aes(x = a_S1, y = b)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 1) +
    facet_grid(~tau) +
    scale_color_viridis_c(name = "Peak time (hr)",
                          breaks = c(4, 8, 12, 16)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Burst Size",
         subtitle = "Lysis time (min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot(data = ysum1, aes(x = tau, y = b)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 1) +
    facet_grid(~a_S1) +
    scale_color_viridis_c(name = "Peak time (hr)",
                          breaks = c(4, 8, 12, 16)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Lysis time (min)",
         y = "Burst Size",
         subtitle = "Infection rate (/min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png("./statplots/figS6_run1_maxtime_contour_all.png", width = 6, height = 6,
      units = "in", res = 300)
  print(cowplot::plot_grid(
    cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                       p2 + theme(legend.position = "none"), 
                       p3 + theme(legend.position = "none"),
                       ncol = 1, labels = "AUTO"),
    cowplot::get_legend(p1),
    rel_widths = c(1, .2),
    ncol = 2))
  dev.off()
  
  p1 <- ggplot(data = ysum1, aes(x = a_S1, y = tau)) +
    geom_contour_filled(aes(z = log10(peak_dens)), alpha = 0.5) +
    geom_point(aes(color = log10(peak_dens), shape = extin_flag),
               size = 1) +
    facet_grid(~b) +
    scale_color_viridis_c(name = "Peak density\n[log10(cfu/mL)]",
                          breaks = c(6, 7, 8, 9)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Lysis time (min)",
         subtitle = "Burst Size") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(data = ysum1, aes(x = a_S1, y = b)) +
    geom_contour_filled(aes(z = log10(peak_dens)), alpha = 0.5) +
    geom_point(aes(color = log10(peak_dens), shape = extin_flag),
               size = 1) +
    facet_grid(~tau) +
    scale_color_viridis_c(name = "Peak density\n[log10(cfu/mL)]",
                          breaks = c(6, 7, 8, 9)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Burst Size",
         subtitle = "Lysis time (min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot(data = ysum1, aes(x = tau, y = b)) +
    geom_contour_filled(aes(z = log10(peak_dens)), alpha = 0.5) +
    geom_point(aes(color = log10(peak_dens), shape = extin_flag),
               size = 1) +
    facet_grid(~a_S1) +
    scale_color_viridis_c(name = "Peak density\n[log10(cfu/mL)]",
                          breaks = c(6, 7, 8, 9)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Lysis time (min)",
         y = "Burst Size",
         subtitle = "Infection rate (/min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png("./statplots/figS7_run1_maxdens_contour_all.png", width = 6, height = 6,
      units = "in", res = 300)
  print(cowplot::plot_grid(
    cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                       p2 + theme(legend.position = "none"), 
                       p3 + theme(legend.position = "none"),
                       ncol = 1, labels = "AUTO"),
    cowplot::get_legend(p1),
    rel_widths = c(1, .25),
    ncol = 2))
  dev.off()
  
  p1 <- ggplot(data = ysum1, aes(x = a_S1, y = tau)) +
    geom_contour_filled(aes(z = log10(extin_time_4/60)), alpha = 0.5) +
    geom_point(aes(color = log10(extin_time_4/60), shape = extin_flag),
               size = 1) +
    facet_grid(~b) +
    scale_color_viridis_c(name = "Extinction time (hr)",
                          breaks = c(0, 0.5, 1, 1.5),
                          labels = c(1, 3.2, 10, 32)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Lysis time (min)",
         subtitle = "Burst Size") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(data = ysum1, aes(x = a_S1, y = b)) +
    geom_contour_filled(aes(z = log10(extin_time_4/60)), alpha = 0.5) +
    geom_point(aes(color = log10(extin_time_4/60), shape = extin_flag),
               size = 1) +
    facet_grid(~tau) +
    scale_color_viridis_c(name = "Extinction time (hr)",
                          breaks = c(0, 0.5, 1, 1.5),
                          labels = c(1, 3.2, 10, 32)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Burst Size",
         subtitle = "Lysis time (min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot(data = ysum1, aes(x = tau, y = b)) +
    geom_contour_filled(aes(z = log10(extin_time_4/60)), alpha = 0.5) +
    geom_point(aes(color = log10(extin_time_4/60), shape = extin_flag),
               size = 1) +
    facet_grid(~a_S1) +
    scale_color_viridis_c(name = "Extinction time (hr)",
                          breaks = c(0, 0.5, 1, 1.5),
                          labels = c(1, 3.2, 10, 32)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Lysis time (min)",
         y = "Burst Size",
         subtitle = "Infection rate (/min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png("./statplots/figS8_run1_extintime_contour_all.png", width = 6, height = 6,
      units = "in", res = 300)
  print(cowplot::plot_grid(
    cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                       p2 + theme(legend.position = "none"), 
                       p3 + theme(legend.position = "none"),
                       ncol = 1, labels = "AUTO"),
    cowplot::get_legend(p1),
    rel_widths = c(1, .28),
    ncol = 2))
  dev.off()
  
  p1 <- ggplot(data = ysum1, aes(x = a_S1, y = tau)) +
    geom_contour_filled(aes(z = log10(auc/60)), alpha = 0.5) +
    geom_point(aes(color = log10(auc/60), shape = extin_flag),
               size = 1) +
    facet_grid(~b) +
    scale_color_viridis_c(name = "Area under the curve\n(hr cfu/mL)",
                          breaks = c(6, 8, 10),
                          labels = c(expression(10^6),
                                     expression(10^8),
                                     expression(10^10))) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Lysis time (min)",
         subtitle = "Burst Size") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(data = ysum1, aes(x = a_S1, y = b)) +
    geom_contour_filled(aes(z = log10(auc/60)), alpha = 0.5) +
    geom_point(aes(color = log10(auc/60), shape = extin_flag),
               size = 1) +
    facet_grid(~tau) +
    scale_color_viridis_c(name = "Area under the curve\n(hr cfu/mL)",
                          breaks = c(6, 8, 10),
                          labels = c(expression(10^6),
                                     expression(10^8),
                                     expression(10^10))) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Burst Size",
         subtitle = "Lysis time (min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot(data = ysum1, aes(x = tau, y = b)) +
    geom_contour_filled(aes(z = log10(auc/60)), alpha = 0.5) +
    geom_point(aes(color = log10(auc/60), shape = extin_flag),
               size = 1) +
    facet_grid(~a_S1) +
    scale_color_viridis_c(name = "Area under the curve\n(hr cfu/mL)",
                          breaks = c(6, 8, 10),
                          labels = c(expression(10^6),
                                     expression(10^8),
                                     expression(10^10))) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Lysis time (min)",
         y = "Burst Size",
         subtitle = "Infection rate (/min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png("./statplots/figS9_run1_auc_contour_all.png", width = 6.2, height = 6,
      units = "in", res = 300)
  print(cowplot::plot_grid(
    cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                       p2 + theme(legend.position = "none"), 
                       p3 + theme(legend.position = "none"),
                       ncol = 1, labels = "AUTO"),
    cowplot::get_legend(p1),
    rel_widths = c(1, .33),
    ncol = 2))
  dev.off()
  
  p1 <- ggplot(data = ysum1, aes(x = a_S1, y = tau)) +
    geom_contour_filled(aes(z = PC1), alpha = 0.5) +
    geom_point(aes(color = PC1, shape = extin_flag),
               size = 1) +
    facet_grid(~b) +
    scale_color_viridis_c() +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Lysis time (min)",
         subtitle = "Burst Size") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(data = ysum1, aes(x = a_S1, y = b)) +
    geom_contour_filled(aes(z = PC1, alpha = 0.5)) +
    geom_point(aes(color = PC1, shape = extin_flag),
               size = 1) +
    facet_grid(~tau) +
    scale_color_viridis_c() +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Burst Size",
         subtitle = "Lysis time (min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot(data = ysum1, aes(x = tau, y = b)) +
    geom_contour_filled(aes(z = PC1), alpha = 0.5) +
    geom_point(aes(color = PC1, shape = extin_flag),
               size = 1) +
    facet_grid(~a_S1) +
    scale_color_viridis_c() +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10", breaks = c(5, 50, 500)) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Lysis time (min)",
         y = "Burst Size",
         subtitle = "Infection rate (/min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png("./statplots/figS10_run1_PC1_contour_all.png", width = 6.2, height = 6,
      units = "in", res = 300)
  print(cowplot::plot_grid(
    cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                       p2 + theme(legend.position = "none"), 
                       p3 + theme(legend.position = "none"),
                       ncol = 1, labels = "AUTO"),
    cowplot::get_legend(p1),
    rel_widths = c(1, .33),
    ncol = 2))
  dev.off()
  
  #Supplemental diminishing returns plots
  plotlist <- list()
  metrics <- c("logpeakdens", "peaktimehr", "logauc", "extintimehr", "PC1")
  metrics_labels <- c("peak\nbacterial density\n[log10(cfu/mL)\n",
                      "time\nof peak bacterial density\n[hr",
                      "\nlog10(area under the curve)\n[log10(hr cfu/mL)\n",
                      "extinction time\n[hr",
                      "PC1\n[")
  for(metric_i in 1:length(metrics)) {
    metric <- metrics[metric_i]
    against_vals <- c("a", "b", "tau")
    against_labels <- c("infection rate", "burst size", "lysis time")
    against_units <- c("/cfu/pfu/mL/min", "pfu/cfu", "min")
    for(against_i in 1:length(against_vals)) {
      against <- against_vals[against_i]
      not_against <- against_vals[against_vals != against]
      for(form in c("log")) {
        #Build base graph
         p <- ggplot(data = ysum1,
                       aes(x = .data[[paste0(form, against)]],
                           y = .data[[paste0(metric, "_d", form, against)]])) +
                  geom_line(alpha = 0.5,
                    aes(group = paste(
                    .data[[paste0(form, not_against[1])]],
                    .data[[paste0(form, not_against[2])]]))) +
           labs(y = paste0("Derivative of ", metrics_labels[metric_i],
                          "/10-fold change in\n", against_labels[against_i],
                          "]"),
                x = paste0("log10[", against_labels[against_i], "]",
                           "\n(log10[", against_units[against_i], "])")) +
           theme_bw()
        #Modify as appropriate
        if(against == "tau") {
          #Plot x-axis in opposite direction
          p <- p + 
            scale_x_continuous(transform = "reverse") +
            scale_y_continuous(transform = "reverse")
        }
      plotlist[[length(plotlist)+1]] <- p
      }
    }
  }

  png("./statplots/figS11_run1_gradient_slopes.png", width = 10, height = 11,
      units = "in", res = 300)
  print(cowplot::plot_grid(plotlist = plotlist,
                           labels = "AUTO",
                           nrow = 5, align = "hv", axis = "tblr"))
  dev.off()
  
  #Supplemental derivative magnitudes plots
  plotlist <- list()
  metrics <- c("logpeakdens", "peaktimehr", "logauc", "extintimehr", "PC1")
  metrics_labels <- c("peak\nbacterial density\n[log10(cfu/mL)\n",
                      "time\nof peak bacterial density\n[hr",
                      "\nlog10(area under the curve)\n[log10(hr cfu/mL)\n",
                      "extinction time\n[hr",
                      "PC1\n[")
  for(metric_i in 1:length(metrics)) {
    metric <- metrics[metric_i]
    forms <- c("log", "norm")
    form_labels <- c("/10-fold change]", "/full range]")
    for(form_i in 1:length(forms)) {
      form <- forms[form_i]
      #Build base graph
      p <- ggplot(
        data = pivot_longer(
          ysum1,
          cols = starts_with(paste0(metric, "_d", form)),
          names_to = "wrt",
          names_prefix = paste0(metric, "_d", form),
          values_to = "derivative"),
        aes(x = wrt, y = abs(derivative), shape = extin_flag)) +
        geom_point(position = position_jitter(width = 0.1, seed = 1),
                   alpha = 0.5) +
        scale_shape_manual(breaks = c("none", "neark", "noextin"),
                           values = c(16, 4, 3)) +
        guides(shape = "none") +
        labs(y = paste0("Magnitude of\nderivative of ",
                        metrics_labels[metric_i],
                        form_labels[form_i]),
             x = "With respect to") +
        scale_x_discrete(breaks = c("a", "b", "tau"),
                         labels = c("Infection rate", "Burst size", "Lysis time")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      plotlist[[length(plotlist)+1]] <- p
    }
  }
  
  png("./statplots/figS12_run1_gradient_magnitudes.png", width = 8, height = 14,
      units = "in", res = 300)
  print(cowplot::plot_grid(plotlist = plotlist,
                           labels = "AUTO",
                           nrow = 5, align = "hv", axis = "tblr"))
  dev.off()
}

# Run 2: bact traits ----
run2 <- run_sims_filewrapper(
  name = "run2",
  u_S1 = signif(0.04*10**seq(from = 0, to = -0.7, length.out = 5), 3),
  u_S2 = 0,
  k = signif(10**c(8, 8.5, 9, 9.5, 10), 3),
  a_S1 = 10**seq(from = -12, to = -8, length.out = 5),
  a_S2 = 0,
  tau = 31.6,
  b = 50,
  z = 1,
  d = 0,
  init_S1 = 10**6,
  init_moi = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig2 <- run2[[1]]

#Set below 0 values to 0
ybig2 <- mutate(ybig2,
                Density = ifelse(Density < 0, 0, Density))

ysum2 <- summarize(group_by(filter(ybig2, Pop == "B"),
                     uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                     tau, b, z, f_a, f_b, d, h, g1, g2,
                     init_S1, init_S2, init_moi, init_N, equil),
            peak_dens = max(Density),
            peak_time = time[which.max(Density)],
            auc = auc(x = time, y = Density),
            extin_time_4 = 
              first_below(y = Density, x = time,
                          threshold = 10**4, return = "x"),
            run_time = max(time))
ysum2 <- mutate(
  ysum2,
  extin_flag = ifelse(is.na(extin_time_4), "noextin",
                      ifelse(peak_dens >= 0.9*k, "neark", "none")),
  extin_time_4 = ifelse(is.na(extin_time_4), run_time, extin_time_4),
  )

run2_preds <-  expand.grid(time = unique(ybig2$time),
                           u_S1 = unique(ybig2$u_S1),
                           k = unique(ybig2$k),
                           init_S1 = unique(ybig2$init_S1))
run2_preds <- mutate(run2_preds,
                     peak_dens_pred = logis_func(S_0 = init_S1, u_S = u_S1, 
                                                 k = k, times = time),
                     auc_pred = logis_def_integral(S_0 = init_S1, u_S = u_S1, 
                                                   k = k, times = time))
                         

# Run 2: stat v stat plots ----
if(glob_make_statplots) {
  fs3a <- 
    ggplot(data = ysum2,
           aes(x = peak_time/60, y = peak_dens)) +
    geom_point(aes(shape = extin_flag)) +
    geom_line(data = run2_preds,
              aes(x = time/60, y = peak_dens_pred),
              lty = 2) +
    facet_nested("k (cfu/mL)" * signif(k, 2) ~ "u_S1 (/hr)" * signif(u_S1*60, 2), 
                 scales = "free_y") +
    scale_x_continuous(breaks = c(0, 12, 24), limits = c(0, 24)) +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Peak Time (hr)", y = "Peak Density\n(cfu/mL)") +
    guides(shape = "none") + 
    theme_bw() +
    theme(axis.title = element_text(size = 20),
          strip.text = element_text(size = 10))
  
  fs3b <- 
    ggplot(data = ysum2,
           aes(x = peak_time/60, y = auc/60)) +
    geom_point(aes(shape = extin_flag)) +
    geom_line(data = run2_preds,
              aes(x = time/60, y = auc_pred/60),
              lty = 2) +
    facet_nested("k (cfu/mL)" * signif(k, 2) ~ "u_S1 (/hr)" * signif(u_S1*60, 2), 
                 scales = "free_y") +
    scale_x_continuous(breaks = c(0, 12, 24), limits = c(0, 24)) +
    scale_y_log10() +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Peak Time (hr)", y = "Area Under the\nCurve (hr cfu/mL)") +
    guides(shape = "none") +
    theme_bw() +
    theme(axis.title = element_text(size = 20),
          strip.text = element_text(size = 10))
  
  fs3c <- 
    ggplot(data = ysum2,
           aes(x = peak_time/60, y = extin_time_4/60)) +
    geom_point(aes(shape = extin_flag)) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    scale_shape_manual(breaks = c("none", "neark", "noextin"),
                       values = c(16, 4, 3)) +
    labs(x = "Peak Time (hr)", y = "Extinction Time\n(hr)") +
    guides(shape = "none") +
    theme_bw() +
    theme(axis.title = element_text(size = 20))
  
  png("./statplots/figS3_run2_metricvmetric_alldata.png", 
      width = 5, height = 12.5, units = "in", res = 300)
  print(plot_grid(fs3a, fs3b, fs3c, ncol = 1, labels = "AUTO",
                  align = "hv", axis = "tb", label_size = 20))
  dev.off()
}

#Run 2: contour plots ----
if (glob_make_statplots) {
  p1 <- ggplot(data = ysum2, aes(x = a_S1, y = 60*u_S1)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 1) +
    facet_grid(~k) +
    scale_color_viridis_c(name = "Peak time (hr)",
                          breaks = c(6, 12, 18, 24)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Bacterial growth\nrate (/hr)",
         subtitle = "Carrying capacity (cfu/mL)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(data = ysum2, aes(x = a_S1, y = k)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 1) +
    facet_grid(~(u_S1*60)) +
    scale_color_viridis_c(name = "Peak time (hr)",
                          breaks = c(6, 12, 18, 24)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Carrying capacity\n(cfu/mL)",
         subtitle = "Bacterial growth rate (/hr)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot(data = ysum2, aes(x = 60*u_S1, y = k)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 1) +
    facet_grid(~a_S1) +
    scale_color_viridis_c(name = "Peak time (hr)",
                          breaks = c(6, 12, 18, 24)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Bacterial growth rate (/hr)",
         y = "Carrying capacity\ncfu/mL)",
         subtitle = "Infection rate (/min)") +
    guides(fill = "none", shape = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  png("./statplots/figS7_run2_maxtime_contour_all.png", width = 6.2, height = 6,
      units = "in", res = 300)
  print(cowplot::plot_grid(
    cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                       p2 + theme(legend.position = "none"), 
                       p3 + theme(legend.position = "none"),
                       ncol = 1, labels = "AUTO"),
    cowplot::get_legend(p1),
    rel_widths = c(1, .2),
    ncol = 2))
  dev.off()
}

## Run 3: plasticity in a ----
run3 <- run_sims_filewrapper(
  name = "run3",
  u_S1 = signif(0.04*10**-0.35, 3), u_S2 = 0,
  k = 10**9,
  a_S1 = 10**seq(from = -12, to = -8, length.out = 5),
  a_S2 = 0,
  tau = 31.6,
  b = 50,
  z = 1,
  d = 0,
  f_a = round(seq(from = 0, to = 3, length.out = 3), 2),
  init_S1 = 10**6,
  init_moi = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig3 <- run3[[1]]

ybig3 <- mutate(ybig3,
                Density = ifelse(Density <= 0, 0, Density))

ybig3 <- interp_data(ybig3, x = "time", y = "Density",
                     subset_by = paste(ybig3$uniq_run, ybig3$Pop))

ybig3_wide <- tidyr::pivot_wider(filter(ybig3, Pop %in% c("N", "B", "P")),
                           names_from = Pop,
                           values_from = Density)
ybig3_wide <- tidyr::pivot_longer(ybig3_wide,
                                  names_to = "Pop",
                                  values_to = "Density",
                                  cols = c("B", "P"))

#Subsample so point types are visible
ybig3_wide <- filter(ybig3_wide,
                     time %% 20 == 0)

if(glob_make_statplots) {
  png("./statplots/figS14_run3_BvsNk.png", width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(data = ybig3_wide,
               aes(x = (k-N)/k, y = Density, color = time/60)) +
          geom_point(size = 1, aes(pch = Pop)) +
          facet_grid(f_a ~ a_S1, scales = "free_y",
                     labeller = labeller(f_a = function(x){paste("f =", x)})) +
          scale_color_viridis_c(direction = -1, name = "Time (hr)",
                                breaks = c(12, 24, 36, 48)) +
          scale_shape_manual(name = "Population",
                             breaks = c("B", "P"),
                             values = c(19, 4),
                             labels = c("Bacteria", "Phages")) +
          labs(x = "Fraction of N consumed", y = "Density (cfu/mL or pfu/mL)",
               subtitle = "Infection rate (/min)") +
          scale_y_log10() +
          geom_vline(data = filter(ybig3_wide, f_a >= 1), 
                     aes(xintercept = 1/f_a), lty = 2) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          NULL)
  dev.off()
}

## Run 4: plasticity in b ----
run4 <- run_sims_filewrapper(
  name = "run4",
  u_S1 = signif(0.04*10**-0.35, 3), u_S2 = 0,
  k = 10**9,
  a_S1 = 10**-10,
  a_S2 = 0,
  tau = 31.6,
  b = signif(5*10**seq(from = 0, to = 2, length.out = 5), 3),
  z = 1,
  d = 0,
  f_b = round(seq(from = 0, to = 3, length.out = 3), 2),
  init_S1 = 10**6,
  init_moi = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig4 <- run4[[1]]

ybig4 <- mutate(ybig4,
                Density = ifelse(Density <= 0, 0, Density))

ybig4 <- interp_data(ybig4, x = "time", y = "Density",
                     subset_by = paste(ybig4$uniq_run, ybig4$Pop))

ybig4_wide <- tidyr::pivot_wider(filter(ybig4, Pop %in% c("N", "B", "P")),
                                 names_from = Pop,
                                 values_from = Density)
ybig4_wide <- tidyr::pivot_longer(ybig4_wide,
                                  names_to = "Pop",
                                  values_to = "Density",
                                  cols = c("B", "P"))

#Subsample so point types are visible
ybig4_wide <- filter(ybig4_wide,
                     time %% 20 == 0)

if(glob_make_statplots) {
  png("./statplots/figS15_run4_BvsNk.png", width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(data = ybig4_wide,
               aes(x = (k-N)/k, y = Density, color = time/60)) +
          geom_point(size = 1, aes(pch = Pop)) +
          facet_grid(f_b ~ b, scales = "free_y",
                     labeller = labeller(f_b = function(x){paste("f =", x)})) +
          scale_color_viridis_c(direction = -1, name = "Time (hr)",
                                breaks = c(12, 24, 36, 48)) +
          scale_shape_manual(name = "Population",
                             breaks = c("B", "P"),
                             values = c(19, 4),
                             labels = c("Bacteria", "Phages")) +
          labs(x = "Fraction of N consumed", y = "Density (cfu/mL or pfu/mL)",
               subtitle = "Burst size") +
          scale_y_log10() +
          geom_vline(data = filter(ybig4_wide, f_b >= 1), 
                     aes(xintercept = 1/f_b), lty = 2) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          NULL)
  dev.off()
}

## Run 5: plasticity in tau ----
run5 <- run_sims_filewrapper(
  name = "run5",
  type = "ode",
  nI = 50,
  u_S1 = signif(0.04*10**-0.35, 3),
  k = 10**9,
  a_S1 = 10**-10,
  tau = signif(10**seq(from = 1, to = 2, length.out = 5), 3),
  b = 50,
  z = 1,
  d = 0,
  f_tau = round(seq(from = 0, to = 3, length.out = 3), 2),
  init_S1 = 10**6,
  init_moi = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = FALSE)#glob_read_files)

ybig5 <- run5[[1]]

ybig5 <- mutate(ybig5,
                Density = ifelse(Density <= 0, 0, Density))

ybig5_interp <- filter(ybig5, Pop %in% c("N", "B", "P"))
ybig5_interp <- interp_data(ybig5_interp, 
                            x = "time", y = "Density",
                     subset_by = paste(ybig5_interp$uniq_run, ybig5_interp$Pop))

ybig5_wide <- tidyr::pivot_wider(ybig5_interp,
                                 names_from = Pop,
                                 values_from = Density)
ybig5_wide <- tidyr::pivot_longer(ybig5_wide,
                                  names_to = "Pop",
                                  values_to = "Density",
                                  cols = c("B", "P"))

#Subsample so point types are visible
ybig5_wide <- filter(ybig5_wide,
                     time %% 20 == 0)

if(glob_make_statplots) {
  png("./statplots/figS16_run5_BvsNk.png", width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(data = ybig5_wide,
               aes(x = (k-N)/k, y = Density, color = time/60)) +
          geom_point(size = 1, aes(pch = Pop)) +
          facet_grid(f_tau ~ tau, scales = "free_y",
                     labeller = labeller(f_tau = function(x){paste("f =", x)})) +
          scale_color_viridis_c(direction = -1, name = "Time (hr)",
                                breaks = c(12, 24, 36, 48)) +
          scale_shape_manual(name = "Population",
                             breaks = c("B", "P"),
                             values = c(19, 4),
                             labels = c("Bacteria", "Phages")) +
          labs(x = "Fraction of N consumed", y = "Density (cfu/mL or pfu/mL)",
               subtitle = "Lysis time (min)") +
          scale_y_log10() +
          geom_vline(data = filter(ybig5_wide, f_tau >= 1), 
                     aes(xintercept = 1/f_tau), lty = 2) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          NULL)
  dev.off()
}

## Run 6: transitions to resistant subpop ----
run6 <- run_sims_filewrapper(
  name = "run6", read_file = glob_read_files,
  a = list(
    u_S1 = signif(0.04*10**-0.35, 3),
    u_S2 = c(0, signif(0.04*10**-0.35, 3)),
    k = 10**9,
    a_S1 = signif(10**seq(from = -12, to = -8, length.out = 9), 3),
    a_S2 = 0,
    tau = 31.6,
    b = 50,
    z = 1,
    d = 0,
    h = c(0, 0.001, 0.01, 0.1),
    g1 = 0,
    init_S1 = 10**6,
    init_moi = 10**-2,
    equil_cutoff_dens = 0.1,
    init_time = 12*60,
    max_time = 48*60,
    init_stepsize = 5,
    print_info = TRUE),
  b = list(
    u_S1 = signif(0.04*10**-0.35, 3),
    u_S2 = c(0, signif(0.04*10**-0.35, 3)),
    k = 10**9,
    a_S1 = signif(10**seq(from = -12, to = -8, length.out = 9), 3),
    a_S2 = 0,
    tau = 31.6,
    b = 50,
    z = 1,
    d = 0,
    h = c(0, 0.001, 0.01, 0.1),
    g1 = 1,
    g2 = c(1, -1),
    init_S1 = 10**6,
    init_moi = 10**-2,
    equil_cutoff_dens = 0.1,
    init_time = 12*60,
    max_time = 48*60,
    init_stepsize = 5,
    print_info = TRUE)
)

ybig6 <- run6[[1]]

#Add flag for class of transition to resistance
#for plotting extend all to end at same time
#Set Density below 0 to 0
ybig6 <- 
  mutate(group_by(ybig6, uniq_run),
         transition = ifelse(g1 == 0, "Constant",
                             ifelse(g2 == 1, "Decreasing with\nN scarcity", 
                                    "Increasing with\nN scarcity")),
         time = ifelse(time == max(time), max(ybig6$time), time),
         Density = ifelse(Density < 1, 0, Density))

if(glob_make_statplots) {
  png("./statplots/fig6_run6_h_Bcurves_constant_uS2_0.png", width = 6, height = 2.5,
      units = "in", res = 300)
  print(ggplot(data = filter(ybig6, Pop == "B", transition == "Constant",
                             u_S2 == 0),
               aes(x = time/60, y = Density, color = log10(a_S1), group = uniq_run)) +
          geom_line() +
          facet_grid(~ h) +
          scale_y_continuous(trans = "log10", limits = c(1, NA)) +
          coord_cartesian(xlim = c(NA, 30)) +
          theme_bw() +
          scale_color_viridis_c(end = 0.95, name = "log10(infection rate)") +
          labs(x = "Time (hr)", y = "Density (cfu/mL)",
               subtitle = "Resistance Transition Rate") +
          NULL)
  dev.off()
  
  png("./statplots/figS17_run6_h_Bcurves_uS2_0.png", width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(data = filter(ybig6, Pop == "B", u_S2 == 0),
               aes(x = time/60, y = Density, color = log10(a_S1), group = uniq_run)) +
          geom_line() +
          facet_grid(transition ~ h) +
          scale_y_continuous(trans = "log10", limits = c(1, NA)) +
          coord_cartesian(xlim = c(NA, 30)) +
          theme_bw() +
          scale_color_viridis_c(end = 0.95, name = "log10(infection rate)") +
          labs(x = "Time (hr)", y = "Density (cfu/mL)",
               subtitle = "Resistance Transition Rate") +
          NULL)
  dev.off()
  
  png("./statplots/figS18_run6_h_Bcurves_uS2not0.png", width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(data = filter(ybig6, Pop == "B", u_S2 != 0),
               aes(x = time/60, y = Density, color = log10(a_S1), group = uniq_run)) +
          geom_line() +
          facet_grid(transition ~ h) +
          scale_y_continuous(trans = "log10", limits = c(1, NA)) +
          coord_cartesian(xlim = c(NA, 30)) +
          theme_bw() +
          scale_color_viridis_c(end = 0.95, name = "log10(infection rate)") +
          labs(x = "Time (hr)", y = "Density (cfu/mL)",
               subtitle = "Resistance Transition Rate") +
          NULL)
  dev.off()
  
}

## Run 7: init dens, init moi, and a ----
run7 <- run_sims_filewrapper(
  name = "run7",
  u_S1 = signif(0.04*10**-0.35, 3), u_S2 = 0,
  k = 10**9,
  a_S1 = rep(10**seq(from = -12, to = -8, length.out = 5), each = 25),
  a_S2 = 0,
  tau = 31.6,
  b = 50,
  z = 1,
  d = 0,
  init_S1 = rep(10**c(6, 6.5, 7, 7.5, 8), times = 5, each = 5),
  init_moi = rep(10**c(4, 4.5, 5, 5.5, 6), times = 25)/
    rep(10**c(6, 6.5, 7, 7.5, 8), times = 5, each = 5),
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files,
  combinatorial = FALSE)

ybig7 <- run7[[1]]

#Set below 0 values to 0
ybig7 <- mutate(ybig7,
                Density = ifelse(Density < 0, 0, Density),
                init_P = init_moi*(init_S1+init_S2))

ysum7 <- summarize(group_by(filter(ybig7, Pop == "B"),
                            uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                            tau, b, z, f_a, f_b, d, h, g1, g2,
                            init_S1, init_S2, init_moi, init_N, init_P, equil),
                   peak_dens = max(Density),
                   peak_time = time[which.max(Density)],
                   auc = auc(x = time, y = Density),
                   extin_time_4 = 
                     first_below(y = Density, x = time,
                                 threshold = 10**4, return = "x",
                                 return_endpoints = FALSE),
                   run_time = max(time),
                   init_P = (init_S1[1] + init_S2[1])*init_moi[1])
ysum7 <- mutate(
  ysum7,
  extin_flag = ifelse(is.na(extin_time_4), "noextin",
                      ifelse(peak_dens >= 0.9*k, "neark", "none")),
  extin_time_4 = ifelse(is.na(extin_time_4), run_time, extin_time_4))

#Add gradients
ysum7 <- mutate(ungroup(ysum7),
                loga = log10(a_S1),
                loginitS1 = log10(init_S1),
                loginitP = log10(init_P),
                logpeakdens = log10(peak_dens),
                peaktimehr = peak_time/60,
                logauc = log10(auc/60),
                extintimehr = extin_time_4/60)

ysum7 <- 
  mutate(group_by(ysum7, loginitS1, loginitP),
         across(.cols = c(logpeakdens, peaktimehr, logauc, extintimehr),
                list("dloga" = ~central_diff(y = .x, x = loga))))
ysum7 <- 
  mutate(group_by(ysum7, loga, loginitP),
         across(.cols = c(logpeakdens, peaktimehr, logauc, extintimehr),
                list("dlogS1" = ~central_diff(y = .x, x = loginitS1))))
ysum7 <- 
  mutate(group_by(ysum7, loga, loginitS1),
         across(.cols = c(logpeakdens, peaktimehr, logauc, extintimehr),
                list("dlogP" = ~central_diff(y = .x, x = loginitP))))

#Add calculations of stochastic noise
#At each point in the init_S x init_P x a_S1 where we can estimate the
# gradient in all three dimensions (so 27 points, since there are 5 levels
# of each but we can't estimate any edge points so its 3x3x3), we get the
# slope of peak time against log10 changes in a, S1, and P
#We then calculate how Poisson noise at the colony-counting stage or
# the pipetting to inoculate stage would change the initial density
# For the colony counting stage, we set the number of colonies counted and
# the number of replicate plates.
#For both, we set the confidence interval (95%) and use qpois to calculate
# the amount up and down at those 27 values
#(For colony counting only) Multiply by the density to get the density up and down
#Convert via slope of S1 or P to get change in metric
#Then back-convert via slope of a to get change in a

##Here's the original logic I wrote down, which is the same steps but I think
##harder to follow than the new version written above:

#At each point in the grid, we know how much a 10-fold change in initial
# bacterial density or in infection rate will shift the peak time
#We can calculate the amount of error up or down (at 95% confidence interval)
# we could expect based on a given number of colonies counted
#So we take the amount of e.g. upper error, calculate the amount that
# error in initS would shift peak time, then calculate the amount of
# error in a that would shift peak time equivalently.

#Let's say a 10-fold change in initS changes peaktimehr by 1 and
# a 10-fold change in a_S1 changes peaktimehr by 2
#At counting 50 colonies, the errors are -26% and +24%
#That's a log10 change of initS of -0.131 and +0.093
#Which would cause a change in peaktimehr of -0.131 hr and +0.093 hr
#Which is the same effect as a log10 change of a_S1 of -0.0655 and +0.0465
#Which is a raw amount of change in a_S1 as 0.86 and 1.113
#So at 50 colonies using the average effects, the confidence interval of
# errors in a_S1 would be 0.86 and 1.113

#So, for each point with gradients in both directions (27x), we just need to
# run these calculations but for all possible # of colonies

ysum7_gradnoise <- dplyr::filter(ysum7,
                                 !is.na(peaktimehr_dloga) &
                                   !is.na(peaktimehr_dlogP) &
                                   !is.na(peaktimehr_dlogS1))

ysum7_gradnoise <- 
  full_join(ysum7_gradnoise,
            expand.grid(uniq_run = ysum7_gradnoise$uniq_run,
                        num_colony = 1:100,
                        num_replicates = 1:5))
ysum7_gradnoise <- 
  mutate(ungroup(ysum7_gradnoise),
         lower_S =
           (10**((log10(qpois(0.025, num_replicates*num_colony)
                        /num_replicates
                        /num_colony 
                        * init_S1) 
                  - loginitS1)
                 * peaktimehr_dlogS1
                 / peaktimehr_dloga
                 + loga)
            / a_S1),
         upper_S =
           (10**((log10(qpois(0.975, num_replicates*num_colony)
                        /num_replicates
                        /num_colony 
                        * init_S1) 
                  - loginitS1)
                 * peaktimehr_dlogS1
                 / peaktimehr_dloga
                 + loga)
            / a_S1),
         lower_P =
           (10**((log10(qpois(0.025, num_replicates*num_colony)
                        /num_replicates
                        /num_colony 
                        * init_P) 
                  - loginitP)
                 * peaktimehr_dlogP
                 / peaktimehr_dloga
                 + loga)
            / a_S1),
         upper_P =
           (10**((log10(qpois(0.975, num_replicates*num_colony)
                        /num_replicates
                        /num_colony 
                        * init_P) 
                  - loginitP)
                 * peaktimehr_dlogP
                 / peaktimehr_dloga
                 + loga)
            / a_S1)
  )

ysum7_gradnoise2 <- dplyr::filter(ysum7,
                                  !is.na(peaktimehr_dloga) &
                                    !is.na(peaktimehr_dlogP) &
                                    !is.na(peaktimehr_dlogS1))
ysum7_gradnoise2 <- 
  mutate(ungroup(ysum7_gradnoise2),
         lower_S =
           (10**((log10(qpois(0.025, init_S1)) 
                  - loginitS1)
                 * peaktimehr_dlogS1
                 / peaktimehr_dloga
                 + loga)
            / a_S1),
         upper_S =
           (10**((log10(qpois(0.975, init_S1)) 
                  - loginitS1)
                 * peaktimehr_dlogS1
                 / peaktimehr_dloga
                 + loga)
            / a_S1),
         lower_P =
           (10**((log10(qpois(0.025, init_P)) 
                  - loginitP)
                 * peaktimehr_dlogP
                 / peaktimehr_dloga
                 + loga)
            / a_S1),
         upper_P =
           (10**((log10(qpois(0.975, init_P)) 
                  - loginitP)
                 * peaktimehr_dlogP
                 / peaktimehr_dloga
                 + loga)
            / a_S1)
  )

if (glob_make_statplots) {
  p1 <- ggplot(data = filter(ysum7, init_S1 == 10**6), 
               aes(x = log10(a_S1), y = log10(init_P))) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 3) +
    scale_color_viridis_c(name = "Time of\nPeak\nBacterial\nDensity (hr)",
                          breaks = c(3, 6, 9, 12)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(labels = math_format(10^.x)) +
    scale_x_continuous(labels = math_format(10^.x)) +
    labs(x = "Infection rate (/min)",
         y = "Initial phage\ndensity (pfu/mL)") +
    guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 16),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13)) +
    NULL
  
  p2 <- ggplot(data = filter(ysum7, init_P == 10**4), 
               aes(x = log10(a_S1), y = log10(init_S1))) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 3) +
    scale_color_viridis_c(name = "Time of\nPeak\nBacterial\nDensity (hr)",
                          breaks = c(4, 8, 12, 16)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(labels = math_format(10^.x)) +
    scale_x_continuous(labels = math_format(10^.x)) +
    labs(x = "Infection rate (/min)",
         y = "Initial bacterial\ndensity (cfu/mL)") +
    guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 16),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13)) +
    NULL
  
  p3 <- ggplot(data = pivot_longer(filter(ysum7, extin_flag == "none"),
                                          cols = starts_with("peaktimehr_dlog"),
                                          names_to = "wrt",
                                          names_prefix = "peaktimehr_dlog",
                                          values_to = "derivative"),
                      aes(x = wrt, y = derivative)) +
    geom_point(position = position_jitter(width = 0.1, seed = 1), alpha = 0.5) +
    labs(y = "Effect on peak time\n(hr/10-fold increase)", 
         x = "") +
    scale_x_discrete(limits = c("P", "S1", "a"),
                     labels = c("Initial phage density", 
                                "Initial bacterial density", 
                                "Infection rate")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  png("./statplots/fig6_run7_inocdensity_a.png", width = 11, height = 3.5,
      units = "in", res = 300)
  print(cowplot::plot_grid(p1, p2, p3, nrow = 1,
                           rel_widths = c(1, 1, 0.5),
                     labels = "AUTO", align = "h", axis = "tb"))
  dev.off()
  
  png("./statplots/figS13_run7_maxtime_a_moiconst_contour_nobars.png", 
      width = 5, height = 3.2,
      units = "in", res = 300)
  print(ggplot(data = filter(ysum7, init_moi > 0.009, init_moi < 0.011), 
               aes(x = log10(a_S1), y = log10(init_S1))) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          scale_color_viridis_c(name = "Time of\nPeak\nBacterial\nDensity (hr)",
                                breaks = c(3, 6, 9, 12)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(breaks = c(6, 7, 8),
                             labels = c(expression(10^6 * ":" * 10^4), 
                                        expression(10^7 * ":" * 10^5), 
                                        expression(10^8 * ":" * 10^6))) +
          scale_x_continuous(labels = math_format(10^.x)) +
          labs(x = "Infection rate (/cfu/pfu/mL/min)",
               y = "Initial bacterial:initial phage density\n(cfu/mL):(pfu/mL)") +
          guides(fill = "none", shape = "none") +
          theme(axis.title = element_text(size = 12),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11)) +
          NULL)
  dev.off()
  
  p1 <- print(ggplot(data = filter(ysum7, init_P == 10**4), 
                     aes(x = a_S1, y = init_S1)) +
                geom_contour_filled(aes(z = peak_dens), alpha = 0.5) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_S1*qpois(0.025, 10)/10,
                #                  yend = init_S1*qpois(0.975, 10)/10)) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_S1*qpois(0.025, 100)/100,
                #                  yend = init_S1*qpois(0.975, 100)/100),
                #              lwd = 2) +
                geom_point(aes(color = peak_dens, shape = extin_flag),
                           size = 3) +
                scale_color_viridis_c(name = "Peak density\n(cfu/mL)",
                                      breaks = c(0, 5*10**8, 10**9),
                                      labels = c(0,
                                                 expression(5%*%10^8),
                                                 expression(10^9)),
                                      limits = c(0, 10**9)) +
                scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                                   values = c(4, 4, 16)) +
                scale_y_continuous(trans = "log10") +
                scale_x_continuous(trans = "log10") +
                labs(x = "Infection rate (/min)",
                     y = "Initial bacterial\ndensity (cfu/mL)") +
                guides(fill = "none", shape = "none") +
                theme(axis.title = element_text(size = 16),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12)) +
                NULL)
  
  p2 <- print(ggplot(data = filter(ysum7, init_S1 == 10**6), 
                     aes(x = a_S1, y = init_P)) +
                geom_contour_filled(aes(z = peak_dens), alpha = 0.5) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_P*qpois(0.025, 10)/10,
                #                  yend = init_P*qpois(0.975, 10)/10)) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_P*qpois(0.025, 100)/100,
                #                  yend = init_P*qpois(0.975, 100)/100),
                #              lwd = 2) +
                geom_point(aes(color = peak_dens, shape = extin_flag),
                           size = 3) +
                scale_color_viridis_c(name = "Peak density\n(cfu/mL)",
                                      breaks = c(0, 5*10**8, 10**9),
                                      labels = c(0,
                                                 expression(5%*%10^8),
                                                 expression(10^9)),
                                      limits = c(0, 10**9)) +
                scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                                   values = c(4, 4, 16)) +
                scale_y_continuous(trans = "log10") +
                scale_x_continuous(trans = "log10") +
                labs(x = "Infection rate (/min)",
                     y = "Initial phage\ndensity (pfu/mL)") +
                guides(fill = "none", shape = "none") +
                theme(axis.title = element_text(size = 16),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12)) +
                NULL)
  
  p3 <- print(ggplot(data = filter(ysum7, init_P == 10**4), 
                     aes(x = a_S1, y = init_S1)) +
                geom_contour_filled(aes(z = log10(extin_time_4/60)), alpha = 0.5) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_S1*qpois(0.025, 10)/10,
                #                  yend = init_S1*qpois(0.975, 10)/10)) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_S1*qpois(0.025, 100)/100,
                #                  yend = init_S1*qpois(0.975, 100)/100),
                #              lwd = 2) +
                geom_point(aes(color = log10(extin_time_4/60), shape = extin_flag),
                           size = 3) +
                scale_color_viridis_c(name = "Extinction time (hr)",
                                      breaks = c(0, 0.5, 1, 1.5),
                                      labels = c(1, 3.2, 10, 32),
                                      limits = c(NA, 1.5)) +
                scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                                   values = c(4, 4, 16)) +
                scale_y_continuous(trans = "log10") +
                scale_x_continuous(trans = "log10") +
                labs(x = "Infection rate (/min)",
                     y = "Initial bacterial\ndensity (cfu/mL)") +
                guides(fill = "none", shape = "none") +
                theme(axis.title = element_text(size = 16),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12)) +
                NULL)
  
  p4 <- print(ggplot(data = filter(ysum7, init_S1 == 10**6), 
                     aes(x = a_S1, y = init_P)) +
                geom_contour_filled(aes(z = log10(extin_time_4/60)), alpha = 0.5) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_P*qpois(0.025, 10)/10,
                #                  yend = init_P*qpois(0.975, 10)/10)) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_P*qpois(0.025, 100)/100,
                #                  yend = init_P*qpois(0.975, 100)/100),
                #              lwd = 2) +
                geom_point(aes(color = log10(extin_time_4/60), shape = extin_flag),
                           size = 3) +
                scale_color_viridis_c(name = "Extinction time (hr)",
                                      breaks = c(0, 0.5, 1, 1.5),
                                      labels = c(1, 3.2, 10, 32),
                                      limits = c(NA, 1.5)) +
                scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                                   values = c(4, 4, 16)) +
                scale_y_continuous(trans = "log10") +
                scale_x_continuous(trans = "log10") +
                labs(x = "Infection rate (/min)",
                     y = "Initial phage\ndensity (pfu/mL)") +
                guides(fill = "none", shape = "none") +
                theme(axis.title = element_text(size = 16),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12)) +
                NULL)
  
  p5 <- print(ggplot(data = filter(ysum7, init_P == 10**4), 
                     aes(x = a_S1, y = init_S1)) +
                geom_contour_filled(aes(z = log10(auc/60)), alpha = 0.5) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_S1*qpois(0.025, 10)/10,
                #                  yend = init_S1*qpois(0.975, 10)/10)) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_S1*qpois(0.025, 100)/100,
                #                  yend = init_S1*qpois(0.975, 100)/100),
                #              lwd = 2) +
                geom_point(aes(color = log10(auc/60), shape = extin_flag),
                           size = 3) +
                scale_color_viridis_c(name = "Area under the curve\n(hr cfu/mL)",
                                      breaks = 7:10,
                                      labels = c(expression(10^7),
                                                 expression(10^8),
                                                 expression(10^9),
                                                 expression(10^10))) +
                scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                                   values = c(4, 4, 16)) +
                scale_y_continuous(trans = "log10") +
                scale_x_continuous(trans = "log10") +
                labs(x = "Infection rate (/min)",
                     y = "Initial bacterial\ndensity (cfu/mL)") +
                guides(fill = "none", shape = "none") +
                theme(axis.title = element_text(size = 16),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12)) +
                NULL)
  
  p6 <- print(ggplot(data = filter(ysum7, init_S1 == 10**6), 
                     aes(x = a_S1, y = init_P)) +
                geom_contour_filled(aes(z = log10(auc/60)), alpha = 0.5) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_P*qpois(0.025, 10)/10,
                #                  yend = init_P*qpois(0.975, 10)/10)) +
                # geom_segment(aes(x = a_S1, xend = a_S1,
                #                  y = init_P*qpois(0.025, 100)/100,
                #                  yend = init_P*qpois(0.975, 100)/100),
                #              lwd = 2) +
                geom_point(aes(color = log10(auc/60), shape = extin_flag),
                           size = 3) +
                scale_color_viridis_c(name = "Area under the curve\n(hr cfu/mL)",
                                      breaks = 7:10,
                                      labels = c(expression(10^7),
                                                 expression(10^8),
                                                 expression(10^9),
                                                 expression(10^10))) +
                scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                                   values = c(4, 4, 16)) +
                scale_y_continuous(trans = "log10") +
                scale_x_continuous(trans = "log10") +
                labs(x = "Infection rate (/min)",
                     y = "Initial phage\ndensity (pfu/mL)") +
                guides(fill = "none", shape = "none") +
                theme(axis.title = element_text(size = 16),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12)) +
                NULL)
  
  png("./statplots/figS14_run7_othermetrics_a_initP_initS_contour.png", 
      width = 10.5, height = 8,
      units = "in", res = 300)
  print(cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2, 
                           labels = "AUTO", align = "hv", axis = "tblr"))
  dev.off()
  
  #Supplemental derivative magnitudes plots
  plotlist <- list()
  metrics <- c("logpeakdens", "peaktimehr", "logauc", "extintimehr")
  metrics_labels <- c("peak\nbacterial density\n[log10(cfu/mL)\n",
                      "time\nof peak bacterial density\n[hr",
                      "\nlog10(area under the curve)\n[log10(hr cfu/mL)\n",
                      "extinction time\n[hr")
  form <- "log"
  form_label <- "/10-fold increase]"
  for(metric_i in 1:length(metrics)) {
    metric <- metrics[metric_i]
    #Build base graph
    p <- ggplot(
      data = pivot_longer(
        ysum7,
        cols = starts_with(paste0(metric, "_d", form)),
        names_to = "wrt",
        names_prefix = paste0(metric, "_d", form),
        values_to = "derivative"),
      aes(x = wrt, y = derivative, shape = extin_flag)) +
      geom_point(position = position_jitter(width = 0.1, seed = 1),
                 alpha = 0.5) +
      scale_shape_manual(breaks = c("none", "neark", "noextin"),
                         values = c(16, 4, 3)) +
      guides(shape = "none") +
      labs(y = paste0("Effect on ",
                      metrics_labels[metric_i],
                      form_label),
           x = "") +
      scale_x_discrete(limits = c("P", "S1", "a"),
                       labels = c("Initial phage density", 
                                  "Initial bacterial density", 
                                  "Infection rate")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plotlist[[length(plotlist)+1]] <- p
  }
  
  png("./statplots/figS15_run7_gradient_magnitudes.png", width = 6, height = 8,
      units = "in", res = 300)
  print(cowplot::plot_grid(plotlist = plotlist,
                           labels = "AUTO",
                           nrow = 2, align = "hv", axis = "tblr"))
  dev.off()

  
  png("./statplots/figS16_run7_peaktime_sensitivity_avsinitestimate.png", 
      width = 8, height = 5, units = "in", res = 300)
  print(
    ggplot(data = pivot_longer(ysum7_gradnoise,
                               cols = c("lower_S", "upper_S", "lower_P", "upper_P"),
                               names_to = c("bound", "popgradient"),
                               names_sep = "_",
                               values_to = "change_in_a"),
           aes(x = num_colony, y = change_in_a)) +
      geom_point(size = 0.5, alpha = 0.1) +
      facet_grid(popgradient ~ num_replicates,
                 labeller = labeller(
                   num_replicates = function(x){return(x)},
                   popgradient = c("P" = "Initial phage density",
                                   "S" = "Initial bacterial density"))) +
      scale_x_log10() +
      labs(x = "Number of colonies/plaques counted (cfu or pfu)",
           y = "Fold-change in infection\nrate obscured by noise",
           subtitle = "Number of replicate plates/spots") +
      theme_bw() +
      theme(axis.title = element_text(size = 18),
            axis.text = element_text(size = 9),
            plot.subtitle = element_text(size = 18),
            strip.text = element_text(size = 12)) +
      NULL)
  dev.off()
  
  
  p1 <- ggplot(data = dplyr::filter(pivot_longer(ysum7_gradnoise2,
                                                 cols = c("lower_S", "upper_S", "lower_P", "upper_P"),
                                                 names_to = c("bound", "popgradient"),
                                                 names_sep = "_",
                                                 values_to = "change_in_a"),
                                    popgradient == "S"),
               aes(x = log10(init_S1), y = change_in_a)) +
    geom_point(size = 0.9, alpha = 0.5) +
    scale_x_continuous(labels = math_format(10^.x),
                       breaks = c(6.5, 7, 7.5)) +
    labs(x = "Target initial density (cfu/mL)",
         y = "Fold-change in infection\nrate obscured by noise") +
    theme_bw() +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          plot.subtitle = element_text(size = 18)) +
    NULL
  p2 <- ggplot(data = dplyr::filter(pivot_longer(ysum7_gradnoise2,
                                                 cols = c("lower_S", "upper_S", "lower_P", "upper_P"),
                                                 names_to = c("bound", "popgradient"),
                                                 names_sep = "_",
                                                 values_to = "change_in_a"),
                                    popgradient == "P"),
               aes(x = log10(init_P), y = change_in_a)) +
    geom_point(size = 0.9, alpha = 0.5) +
    scale_x_continuous(labels = math_format(10^.x),
                       breaks = c(4.5, 5, 5.5)) +
    labs(x = "Target initial density (pfu/mL)",
         y = "Fold-change in infection\nrate obscured by noise") +
    theme_bw() +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          plot.subtitle = element_text(size = 18)) +
    NULL
  
  png("./statplots/figS17_run7_peaktime_sensitivity_avsinitsample.png", 
      width = 10.5, height = 4, units = "in", res = 300)
  print(cowplot::plot_grid(p1, p2, nrow = 1, labels = "AUTO"))
  dev.off()
  
  #Noise bars from estimating density (old fig S23)
  ggplot(data = filter(ysum7, init_S1 == 10**6), 
               aes(x = a_S1, y = init_P)) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          geom_segment(aes(x = a_S1, xend = a_S1,
                           y = init_P*qpois(0.025, 10)/10,
                           yend = init_P*qpois(0.975, 10)/10)) +
          geom_segment(aes(x = a_S1, xend = a_S1,
                           y = init_P*qpois(0.025, 100)/100,
                           yend = init_P*qpois(0.975, 100)/100),
                       lwd = 2) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          scale_color_viridis_c(name = "Peak time (hr)",
                                breaks = c(3, 6, 9, 12)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Infection rate (/min)",
               y = "Initial phage\ndensity (pfu/mL)") +
          guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
          NULL
  
  ggplot(data = filter(ysum7, init_P == 10**4), 
               aes(x = a_S1, y = init_S1)) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          geom_segment(aes(x = a_S1, xend = a_S1,
                           y = init_S1*qpois(0.025, 10)/10,
                           yend = init_S1*qpois(0.975, 10)/10)) +
          geom_segment(aes(x = a_S1, xend = a_S1,
                           y = init_S1*qpois(0.025, 100)/100,
                           yend = init_S1*qpois(0.975, 100)/100),
                       lwd = 2) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          scale_color_viridis_c(name = "Peak time (hr)",
                                breaks = c(4, 8, 12, 16)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Infection rate (/min)",
               y = "Initial bacterial\ndensity (cfu/mL)") +
          guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
          NULL
  
  #Noise from sampling from flask (old fig S25)
  ggplot(data = filter(ysum7, init_P == 10**4), 
               aes(x = a_S1, y = init_S1)) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          geom_segment(aes(x = a_S1, xend = a_S1,
                           y = qpois(0.025, init_S1),
                           yend = qpois(0.975, init_S1)),
                       lwd = 1) +
          scale_color_viridis_c(name = "Peak time (hr)",
                                breaks = c(6, 12, 18, 24)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Infection rate (/min)",
               y = "Initial bacterial\ndensity (cfu/mL)") +
          guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
          NULL
        
  ggplot(data = filter(ysum7, init_S1 == 10**6), 
               aes(x = a_S1, y = init_P)) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          geom_segment(aes(x = a_S1, xend = a_S1,
                           y = qpois(0.025, init_P),
                           yend = qpois(0.975, init_P)),
                       lwd = 1) +
          scale_color_viridis_c(name = "Peak time (hr)",
                                breaks = c(3, 6, 9, 12)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Infection rate (/min)",
               y = "Initial phage]\ndensity (pfu/mL)") +
          guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
          NULL

  #Noise from multiple replicate plates to estimate density (old Fig S24)
  ggplot(data = filter(ysum7, init_P == 10**4), 
               aes(x = a_S1, y = init_S1)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_segment(aes(x = a_S1, xend = a_S1,
                     y = init_S1*qpois(0.025, 1*10)/1/10,
                     yend = init_S1*qpois(0.975, 1*10)/1/10)) +
    geom_segment(aes(x = a_S1, xend = a_S1,
                     y = init_S1*qpois(0.025, 3*10)/3/10,
                     yend = init_S1*qpois(0.975, 3*10)/3/10),
                 lwd = 2) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 3) +
    scale_color_viridis_c(name = "Peak time (hr)",
                          breaks = c(6, 12, 18, 24)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Initial bacterial\ndensity (cfu/mL)") +
    guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    NULL
  
  ggplot(data = filter(ysum7, init_S1 == 10**6), 
               aes(x = a_S1, y = init_P)) +
    geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
    geom_segment(aes(x = a_S1, xend = a_S1,
                     y = init_P*qpois(0.025, 1*10)/1/10,
                     yend = init_P*qpois(0.975, 1*10)/1/10)) +
    geom_segment(aes(x = a_S1, xend = a_S1,
                     y = init_P*qpois(0.025, 3*10)/3/10,
                     yend = init_P*qpois(0.975, 3*10)/3/10),
                 lwd = 2) +
    geom_point(aes(color = peak_time/60, shape = extin_flag),
               size = 3) +
    scale_color_viridis_c(name = "Peak time (hr)",
                          breaks = c(3, 6, 9, 12)) +
    scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                       values = c(4, 4, 16)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    labs(x = "Infection rate (/min)",
         y = "Initial phage\ndensity (pfu/mL)") +
    guides(fill = "none", shape = "none") +
    theme(axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    NULL
}


## Run 8: init dens, init moi, and b ----
run8 <- run_sims_filewrapper(
  name = "run8",
  u_S1 = signif(0.04*10**-0.35, 3), u_S2 = 0,
  k = 10**9,
  a_S1 = 10**-10,
  a_S2 = 0,
  tau = 31.6,
  b = signif(5*10**seq(from = 0, to = 2, length.out = 5), 3),
  z = 1,
  d = 0,
  init_S1 = 10**(4:8),
  init_moi = 10**(-4:0),
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig8 <- run8[[1]]

#Set below 0 values to 0
ybig8 <- mutate(ybig8,
                Density = ifelse(Density < 0, 0, Density),
                init_P = init_moi*(init_S1+init_S2))

ysum8 <- summarize(group_by(filter(ybig8, Pop == "B"),
                            uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                            tau, b, z, f_a, f_b, d, h, g1, g2,
                            init_S1, init_S2, init_moi, init_N, init_P, equil),
                   peak_dens = max(Density),
                   peak_time = time[which.max(Density)],
                   auc = auc(x = time, y = Density),
                   extin_time_4 = 
                     first_below(y = Density, x = time,
                                 threshold = 10**4, return = "x",
                                 return_endpoints = FALSE),
                   run_time = max(time))
ysum8 <- mutate(
  ysum8,
  extin_flag = ifelse(is.na(extin_time_4), "noextin",
                      ifelse(peak_dens >= 0.9*k, "neark", "none")),
  extin_time_4 = ifelse(is.na(extin_time_4), run_time, extin_time_4))

if (glob_make_statplots) {
  p1 <- ggplot(data = filter(ysum8, init_P == 10**4), 
               aes(x = b, y = init_S1)) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          # geom_segment(aes(x = b, xend = b,
          #                  y = init_S1*qpois(0.025, 10)/10,
          #                  yend = init_S1*qpois(0.975, 10)/10)) +
          # geom_segment(aes(x = b, xend = b,
          #                  y = init_S1*qpois(0.025, 100)/100,
          #                  yend = init_S1*qpois(0.975, 100)/100),
          #              lwd = 2) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          scale_color_viridis_c(name = "Peak time (hr)",
                                breaks = c(6, 12, 18, 24)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Burst size",
               y = "Initial bacterial density\n(cfu/mL)") +
          guides(fill = "none", shape = "none") +
          NULL

  p2 <- ggplot(data = filter(ysum8, init_S1 == 10**6), 
               aes(x = b, y = init_P)) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          # geom_segment(aes(x = b, xend = b,
          #                  y = init_P*qpois(0.025, 10)/10,
          #                  yend = init_P*qpois(0.975, 10)/10)) +
          # geom_segment(aes(x = b, xend = b,
          #                  y = init_P*qpois(0.025, 100)/100,
          #                  yend = init_P*qpois(0.975, 100)/100),
          #              lwd = 2) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          scale_color_viridis_c(name = "Peak time (hr)",
                                breaks = c(3, 6, 9, 12)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Burst size",
               y = "Initial phage density\n(pfu/mL)") +
          guides(fill = "none", shape = "none") +
          NULL
  
  png("./statplots/figS20_run8_peaktime_b_initP_initS_contour.png", 
      width = 8, height = 2.5,
      units = "in", res = 300)
  print(cowplot::plot_grid(p1, p2, ncol = 2, labels = "AUTO"))
  dev.off()
}

## Run 9: init dens, init moi, and tau ----
run9 <- run_sims_filewrapper(
  name = "run9",
  u_S1 = signif(0.04*10**-0.35, 3), u_S2 = 0,
  k = 10**9,
  a_S1 = 10**-10,
  a_S2 = 0,
  tau = signif(10**seq(from = 1, to = 2, length.out = 5), 3),
  b = 50,
  z = 1,
  d = 0,
  init_S1 = 10**(4:8),
  init_moi = 10**(-4:0),
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig9 <- run9[[1]]

#Set below 0 values to 0
ybig9 <- mutate(ybig9,
                Density = ifelse(Density < 0, 0, Density),
                init_P = init_moi*(init_S1+init_S2))

ysum9 <- summarize(group_by(filter(ybig9, Pop == "B"),
                            uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                            tau, b, z, f_a, f_b, d, h, g1, g2,
                            init_S1, init_S2, init_moi, init_N, init_P, equil),
                   peak_dens = max(Density),
                   peak_time = time[which.max(Density)],
                   auc = auc(x = time, y = Density),
                   extin_time_4 = 
                     first_below(y = Density, x = time,
                                 threshold = 10**4, return = "x",
                                 return_endpoints = FALSE),
                   run_time = max(time))
ysum9 <- mutate(
  ysum9,
  extin_flag = ifelse(is.na(extin_time_4), "noextin",
                      ifelse(peak_dens >= 0.9*k, "neark", "none")),
  extin_time_4 = ifelse(is.na(extin_time_4), run_time, extin_time_4))

if (glob_make_statplots) {
  p1 <- ggplot(data = filter(ysum9, init_P == 10**4), 
               aes(x = tau, y = init_S1)) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          # geom_segment(aes(x = tau, xend = tau,
          #                  y = init_S1*qpois(0.025, 10)/10,
          #                  yend = init_S1*qpois(0.975, 10)/10)) +
          # geom_segment(aes(x = tau, xend = tau,
          #                  y = init_S1*qpois(0.025, 100)/100,
          #                  yend = init_S1*qpois(0.975, 100)/100),
          #              lwd = 2) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          scale_color_viridis_c(name = "Peak time (hr)",
                                breaks = c(6, 12, 18, 24)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Lysis time (min)",
               y = "Initial bacterial density\n(cfu/mL)") +
          guides(fill = "none", shape = "none") +
          NULL

  p2 <- ggplot(data = filter(ysum9, init_S1 == 10**6), 
               aes(x = tau, y = init_P)) +
          geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
          # geom_segment(aes(x = tau, xend = tau,
          #                  y = init_P*qpois(0.025, 10)/10,
          #                  yend = init_P*qpois(0.975, 10)/10)) +
          # geom_segment(aes(x = tau, xend = tau,
          #                  y = init_P*qpois(0.025, 100)/100,
          #                  yend = init_P*qpois(0.975, 100)/100),
          #              lwd = 2) +
          geom_point(aes(color = peak_time/60, shape = extin_flag),
                     size = 3) +
          scale_color_viridis_c(name = "Peak time (hr)",
                                breaks = c(3, 6, 9, 12)) +
          scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                             values = c(4, 4, 16)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Lysis time (min)",
               y = "Initial phage density\n(pfu/mL)") +
          guides(fill = "none", shape = "none") +
          NULL
  
  png("./statplots/figS21_run9_peaktime_tau_initP_initS_contour.png", 
      width = 8, height = 2.5,
      units = "in", res = 300)
  print(cowplot::plot_grid(p1, p2, ncol = 2, labels = "AUTO"))
  dev.off()
}


## Run 10: test of metrics across dift bact ----

run10 <- run_sims_filewrapper(
  name = "run10",
  read_file = glob_read_files,
  a = list(
    u_S1 = signif(0.04*10**seq(from = 0, to = -0.7, length.out = 5), 3),
    u_S2 = 0,
    k = signif(10**seq(from = 8, to = 10, length.out = 5)),
    a_S1 = signif(10**seq(from = -12, to = -8, length.out = 5), 3),
    a_S2 = 0,
    tau = 31.6,
    b = 50,
    z = 1,
    d = 0,
    h = c(0, 0.1),
    init_S1 = 10**6,
    init_moi = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
    equil_cutoff_dens = 0.1,
    init_time = 12*60,
    max_time = 48*60,
    init_stepsize = 5,
    print_info = TRUE),
  b = list(
    u_S1 = signif(0.04*10**seq(from = 0, to = -0.7, length.out = 5), 3),
    u_S2 = 0,
    k = signif(10**seq(from = 8, to = 10, length.out = 5)),
    a_S1 = 0,
    a_S2 = 0,
    tau = 31.6,
    b = 50,
    z = 1,
    d = 0,
    h = c(0, 0.1),
    init_S1 = 10**6,
    init_moi = 0,
    equil_cutoff_dens = 0.1,
    init_time = 12*60,
    max_time = 48*60,
    init_stepsize = 5,
    print_info = TRUE))
    
ybig10 <- run10[[1]]

#Set below 0 values to 0
#and set all runs to end at the same time
ybig10 <- mutate(group_by(ybig10, uniq_run),
                Density = ifelse(Density < 0, 0, Density),
                time = ifelse(time == max(time), max(ybig10$time), time))

ysum10 <- summarize(group_by(filter(ybig10, Pop == "B"),
                            uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                            tau, b, z, f_a, f_b, d, h, g1, g2,
                            init_S1, init_S2, init_moi, init_N, equil),
                   peak_dens = max(Density),
                   peak_time = time[which.max(Density)],
                   auc = auc(x = time, y = Density),
                   extin_time_4 = 
                     first_below(y = Density, x = time,
                                 threshold = 10**4, return = "x",
                                 return_endpoints = FALSE),
                   run_time = max(time))
ysum10 <- mutate(
  ysum10,
  extin_flag = ifelse(is.na(extin_time_4), "noextin",
                      ifelse(peak_dens >= 0.9*k, "neark", "none")),
  extin_time_4 = ifelse(is.na(extin_time_4), run_time, extin_time_4),
  bact = paste(u_S1, k, h))

dup_rows <- function(mydf, mut_col, from_vals, to_vals) {
  #this function takes in a data.frame
  #finds all rows where df[, mut_col] == from_vals
  #deletes those rows
  #and adds rows that replicate the original but with to_vals in the
  # mut_col instead
  
  mydf <- as.data.frame(mydf)
  
  #Generate vector of rows to duplicate
  rowsrep <- ifelse(mydf[, mut_col] %in% from_vals,
                    length(to_vals), 1)
  idx <- rep(1:nrow(mydf), rowsrep)
  
  #Duplicate rows
  newdf <- mydf[idx, ]
  
  #replace from_vals with to_vals in duplicated rows
  newdf[which(idx %in% which(rowsrep > 1)), mut_col] <- to_vals
  
  return(newdf)
}

#Change/duplicate rows where a_S1 == 0 to be the 5 a_S1 vals
ysum10 <- dup_rows(ysum10, mut_col = "a_S1", from_vals = 0, 
                   to_vals = 10**c(-8:-12))

#Relative auc
ysum10 <- mutate(group_by(ysum10, u_S1, u_S2, k, z, d, h,
                         init_S1, init_S2, init_N),
                rel_auc = auc/(auc[init_moi == 0][1]),
                ref_auc = auc[init_moi == 0][1])

#Plots
if(glob_make_statplots) {
  p1 <- ggplot(data = filter(ybig10, 
                             a_S1 %in% c(0),
                             u_S1 %in% c(0.0179),
                             k %in% c(10**8, 10**9, 10**10),
                             Pop == "B", h == 0,
                             init_moi %in% c(0)),
               aes(x = time/60, y = Density, 
                   group = paste(u_S1, a_S1, k, init_moi))) +
    geom_area(aes(color = as.factor(k), 
                  fill = as.factor(k)),
              position = "identity", alpha = 0.2) +
    coord_cartesian(xlim = c(0, 24)) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    facet_grid( ~ a_S1, scales = "fixed",
                labeller = labeller("a_S1" = c("0" = "Control"))) +
    scale_color_viridis_d(name = "Bacterial\ncarrying\ncapacity\n(cfu/mL)",
                          end = 0.9) +
    scale_fill_viridis_d(name = "Bacterial\ncarrying\ncapacity\n(cfu/mL)",
                         end = 0.9) +
    theme_bw() +
    labs(x = "Time (hr)", y = "Density (cfu/mL)") +
    guides(color = "none", fill = "none")
  
  p2 <- ggplot(data = filter(ybig10, 
                             a_S1 %in% c(10**-10, 10**-11, 10**-12),
                             u_S1 %in% c(0.0179),
                             k %in% c(10**8, 10**9, 10**10),
                             Pop == "B", h == 0,
                             init_moi %in% c(0.01)),
               aes(x = time/60, y = Density, 
                   group = paste(u_S1, a_S1, k, init_moi))) +
    geom_area(aes(color = as.factor(k), 
                  fill = as.factor(k)),
              position = "identity", alpha = 0.2) +
    coord_cartesian(xlim = c(0, 24)) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    facet_grid( ~ a_S1, scales = "fixed") +
    scale_color_viridis_d(name = "Bacterial\ncarrying\ncapacity\n(cfu/mL)",
                          end = 0.9) +
    scale_fill_viridis_d(name = "Bacterial\ncarrying\ncapacity\n(cfu/mL)",
                         end = 0.9) +
    theme_bw() +
    labs(x = "Time (hr)", y = "Density (cfu/mL)",
         subtitle = "Infection rate (/cfu/pfu/min)") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  p3 <- ggplot(data = filter(ybig10, 
                             a_S1 %in% c(0),
                             k %in% c(10**9),
                             u_S1 %in% c(0.00798, 0.0179, 0.04),
                             Pop == "B", h == 0,
                             init_moi %in% c(0)),
               aes(x = time/60, y = Density, 
                   group = paste(u_S1, a_S1, k, init_moi))) +
    geom_area(aes(color = as.factor(u_S1*60), 
                  fill = as.factor(u_S1*60)),
              position = "identity", alpha = 0.2) +
    coord_cartesian(xlim = c(0, 24)) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    facet_grid( ~ a_S1, scales = "fixed",
                labeller = labeller("a_S1" = c("0" = "Control"))) +
    scale_color_viridis_d(name = "Bacterial\ngrowth\nrate (/hr)", end = 0.9) +
    scale_fill_viridis_d(name = "Bacterial\ngrowth\nrate (/hr)", end = 0.9) +
    theme_bw() +
    labs(x = "Time (hr)", y = "Density (cfu/mL)") +
    guides(color = "none", fill = "none")
  
  p4 <- ggplot(data = filter(ybig10, 
                             a_S1 %in% c(10**-10, 10**-11, 10**-12),
                             k %in% c(10**9),
                             u_S1 %in% c(0.00798, 0.0179, 0.04),
                             Pop == "B", h == 0,
                             init_moi %in% c(0.01)),
               aes(x = time/60, y = Density, 
                   group = paste(u_S1, a_S1, k, init_moi))) +
    geom_area(aes(color = as.factor(u_S1*60), 
                  fill = as.factor(u_S1*60)),
              position = "identity", alpha = 0.2) +
    coord_cartesian(xlim = c(0, 24)) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    facet_grid( ~ a_S1, scales = "fixed") +
    scale_color_viridis_d(name = "Bacterial\ngrowth\nrate (/hr)", end = 0.9) +
    scale_fill_viridis_d(name = "Bacterial\ngrowth\nrate (/hr)", end = 0.9) +
    theme_bw() +
    labs(x = "Time (hr)", y = "Density (cfu/mL)") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  png("./statplots/fig8_run10_auccurves.png", width = 6, height = 4,
      units = "in", res = 300)
  print(cowplot::plot_grid(
    cowplot::plot_grid(p1, p2, rel_widths = c(0.45, 1), 
                       align = "hv", axis = "tb"),
    cowplot::plot_grid(p3, p4, rel_widths = c(0.45, 1), 
                       align = "hv", axis = "tb"),
    nrow = 2, labels = "AUTO"))
  dev.off()
  
  mycolors <- c("black", scales::viridis_pal(end = 0.9)(5))
  
  fs27a <-
    print(ggplot(data = filter(ysum10, h == 0, init_moi == 0.01),
         aes(x = log10(ref_auc), y = log10(auc), 
             fill = as.factor(a_S1), color = as.factor(a_S1))) +
    geom_point(size = 2, aes(shape = as.factor(k))) +
    geom_abline(slope = 1) +
    geom_abline(lty = 3, lwd = 0.5, color = "gray50",
                slope = 0.5, intercept = seq(1, 9, 0.5)) +
    lims(y = log10(c(min(ysum10$auc, ysum10$ref_auc), 
                     max(ysum10$auc, ysum10$ref_auc)))) +
    #geom_smooth(method = "lm", se = FALSE) +
    scale_fill_manual(name = "Infection rate (/min)",
                      breaks = 10**(-12:-8),
                      values = mycolors[1:6],
                      labels = c(expression(10^-12),
                                 expression(10^-11), expression(10^-10),
                                 expression(10^-9), expression(10^-8))) +
    scale_color_manual(name = "Infection rate (/min)",
                       breaks = 10**(-12:-8),
                       values = mycolors[1:6],
                       labels = c(expression(10^-12),
                                  expression(10^-11), expression(10^-10),
                                  expression(10^-9), expression(10^-8))) +
    scale_shape_manual(name = "Carrying capacity\n(cfu/mL)",
                       values = 21:25) +
    guides(shape = guide_legend(override.aes = list(fill = "black"))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
      labs(x = "log10(Control AUC)\n(hr cfu/mL)", 
           y = "log10(AUC)\n(hr cfu/mL)"))
  fs27b <-
    print(ggplot(data = filter(ysum10, h == 0, init_moi == 0.01),
         aes(x = log10(ref_auc), y = log10(rel_auc), 
             color = as.factor(a_S1), fill = as.factor(a_S1))) +
    geom_point(aes(shape = as.factor(k))) +
    geom_abline(slope = 0) +
    geom_abline(lty = 3, lwd = 0.5, color = "gray50",
                slope = -0.5, intercept = seq(0, 8, 0.5)) +
    #geom_smooth(method = "lm", se = FALSE) +
    scale_fill_manual(name = "Infection rate (/min)",
                      breaks = 10**(-12:-8),
                      values = mycolors[1:6],
                      labels = c(expression(10^-12),
                                 expression(10^-11), expression(10^-10),
                                 expression(10^-9), expression(10^-8))) +
    scale_color_manual(name = "Infection rate (/min)",
                       breaks = 10**(-12:-8),
                       values = mycolors[1:6],
                       labels = c(expression(10^-12),
                                  expression(10^-11), expression(10^-10),
                                  expression(10^-9), expression(10^-8))) +
    scale_shape_manual(name = "Carrying capacity\n(cfu/mL)",
                       values = 21:25) +
    guides(shape = guide_legend(override.aes = list(fill = "black"))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
    labs(x = "log10(Control AUC)\n(hr cfu/mL)", 
         y = "log10(Relative AUC)"))
  
  png("./statplots/figS27_run10_relauc_controlauc_all.png", 
      width = 4.5, height = 6,
      units = "in", res = 300)
  print(
    cowplot::plot_grid(
      cowplot::plot_grid(
        fs27a + guides(shape = "none", color = "none", fill = "none"), 
        fs27b + guides(shape = "none", color = "none", fill = "none"), 
        ncol = 1, labels = "AUTO", align = "hv", axis = "tblr"),
      get_legend(fs27a), 
      ncol = 2, rel_widths = c(1, 0.6)))
  dev.off()
  
  png("./statplots/figS28_run10_relauc_moi_VirulenceIndexNull.png", 
      width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(filter(ysum10, h == 0, init_moi != 0),
               aes(x = log10(init_moi), y = 1-rel_auc, color = as.factor(log10(a_S1)))) +
          geom_point(alpha = 0.7) +
          #scale_x_log10() +
          facet_nested("k (cfu/mL)" * signif(k, 2) ~ 
                         "u_S1 (/hr)" * signif(60*u_S1, 2)) +
          labs(x = "log10(initial MOI)", y = "Virulence Index") +
          scale_color_viridis_d(name = "log10(infection rate)",
                                direction = 1,
                                guide = guide_legend(reverse = TRUE)) +
          theme_bw() +
          theme(axis.title = element_text(size = 16),
                axis.text = element_text(size = 7),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                strip.text = element_text(size = 9.5)) +
          NULL)
  dev.off()
}


# Run 10: PCA ----
ybig10_B <- filter(ybig10, Pop == "B")
ybig10_B <- interp_data(df = ybig10_B,
                        x = "time", y = "Density",
                        subset_by = ybig10_B$uniq_run)
ybig10_B <- mutate(group_by(ybig10_B, u_S1, k, h, time),
                   Dens_norm = Density - Density[init_moi == 0])
ybig10_B <- filter(ybig10_B, time %% 20 == 0)

ybig10_B_wide <- tidyr::pivot_wider(ybig10_B,
                                    names_from = time,
                                    names_prefix = "t_",
                                    values_from = c(Density, Dens_norm))

mypca <- prcomp(ybig10_B_wide[, grep("Density_", colnames(ybig10_B_wide))[-1]],
                center = TRUE, scale = TRUE, retx = TRUE)
mypcanorm <- prcomp(ybig10_B_wide[, grep("Dens_norm", colnames(ybig10_B_wide))[-1]],
                center = TRUE, scale = TRUE, retx = TRUE)

#Merge with orig data
colnames(mypcanorm$x) <- paste0("norm_", colnames(mypcanorm$x))
ybig10_B_wide <- cbind(ybig10_B_wide,
                       as.data.frame(mypca$x),
                       as.data.frame(mypcanorm$x))

#Plot
if(glob_make_statplots) {
  fs26a <- 
    ggplot(data = ybig10_B_wide,
           aes(x = PC1, y = PC2)) +
    geom_point(aes(color = as.factor(a_S1)), alpha = 0.7) +
    scale_color_viridis_d(name = "Infection rate\n(/min)", end = 0.85,
                          labels = c("NA", expression(10^-12),
                                     expression(10^-11), expression(10^-10),
                                     expression(10^-9), expression(10^-8))) +
    labs(x = paste("PC1 (",
                   round((100*((mypca$sdev)**2)/
                            sum((mypca$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (",
                   round((100*((mypca$sdev)**2)/
                            sum((mypca$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme_bw() +
    theme(axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    NULL
  
  fs26b <- 
    ggplot(data = ybig10_B_wide,
           aes(x = norm_PC1, y = norm_PC2)) +
    geom_point(aes(color = as.factor(a_S1)), alpha = 0.7) +
    scale_color_viridis_d(name = "Infection rate\n(/min)", end = 0.85,
                          labels = c("NA", expression(10^-12),
                                     expression(10^-11), expression(10^-10),
                                     expression(10^-9), expression(10^-8))) +
    labs(x = paste("PC1 (",
                   round((100*((mypcanorm$sdev)**2)/
                            sum((mypcanorm$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (",
                   round((100*((mypcanorm$sdev)**2)/
                            sum((mypcanorm$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme_bw() +
    theme(axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    NULL
  
  png("./statplots/figS26_run10_pca.png", 
      width = 8, height = 3,
      units = "in", res = 300)
  print(
    cowplot::plot_grid(fs26a + guides(color = "none"), 
                       fs26b + guides(color = "none"),
                       get_legend(fs26a),
                       nrow = 1, align = "hv", axis = "tblr",
                       labels = c("A", "B", ""), rel_widths = c(1, 1, 0.5)))
  dev.off()
}

# Run 11: Transitions to resistant subpop (at mutational rates) ----
run11 <- run_sims_filewrapper(
  name = "run11", read_file = glob_read_files,
  u_S1 = signif(0.04*10**-0.35, 3),
  u_S2 = signif(0.04*10**-0.35, 3),
  k = 10**9,
  a_S1 = signif(10**seq(from = -12, to = -8, length.out = 7), 3),
  a_S2 = 0,
  tau = 31.6,
  b = 50,
  z = 1,
  d = c(0, 1),
  h = 10**(-2:-8),
  g1 = 1,
  g2 = 1,
  init_S1 = 10**6,
  init_moi = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 4*24*60,
  max_time = 4*24*60,
  init_stepsize = 5,
  print_info = TRUE
)

ybig11 <- run11[[1]]

#for plotting extend all to end at same time
#Set Density below 0 to 0
ybig11 <- 
  mutate(group_by(ybig11, uniq_run),
         time = ifelse(time == max(time), max(ybig11$time), time),
         Density = ifelse(Density < 1, 0, Density))

ysum11 <- summarize(group_by(filter(ybig11, Pop == "B"),
                             uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                             tau, b, z, f_a, f_b, d, h, g1, g2,
                             init_S1, init_S2, init_moi, init_N, equil),
                    peak_dens = max(Density),
                    extin_time = 
                      first_minima(y = Density, x = time,
                                   window_width = 2*60, return = "x",
                                  return_endpoints = FALSE),
                    emerg_time_6 =
                      first_above(y = Density[time > extin_time], 
                                  x = time[time > extin_time],
                                  threshold = 10**6, return = "x"))

if(glob_make_curveplots) {
  print(ggplot(data = filter(ybig11, d == 0, Pop == "B"),
               aes(x = time, y = Density)) +
          geom_line() +
          scale_y_log10() +
          facet_grid(a_S1 ~ h) +
          geom_vline(data = filter(ysum11, d == 0), aes(xintercept = extin_time),
                     color = "red") +
          geom_vline(data = filter(ysum11, d == 0), aes(xintercept = emerg_time_6),
                     color = "blue") +
          NULL)
  print(ggplot(data = filter(ybig11, d == 1, Pop == "B"),
               aes(x = time, y = Density)) +
          geom_line() +
          scale_y_log10() +
          facet_grid(a_S1 ~ h) +
          geom_vline(data = filter(ysum11, d == 1), aes(xintercept = extin_time),
                     color = "red") +
          geom_vline(data = filter(ysum11, d == 1), aes(xintercept = emerg_time_6),
                     color = "blue") +
          NULL)
}


if (glob_make_statplots) {
  print(ggplot(data = ysum11, 
               aes(x = a_S1, y = h)) +
    geom_contour_filled(aes(z = extin_time/60), alpha = 0.5) +
    geom_point(aes(color = extin_time/60),
               size = 3) +
    scale_color_viridis_c(name = "Extinction time (hr)",
                          breaks = c(6, 12, 18, 24)) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    guides(fill = "none", shape = "none") +
      facet_grid(~d) +
    NULL)
  
  png("./statplots/figSXX_run11_emergencetime_a_mutrate_contour.png", 
      width = 4.5, height = 5, units = "in", res = 300)
  print(ggplot(data = ysum11, 
               aes(x = a_S1, y = h)) +
          geom_contour_filled(aes(z = emerg_time_6/60), alpha = 0.5) +
          geom_point(aes(color = emerg_time_6/60),
                     size = 3) +
          scale_color_viridis_c(name = "Emergence\ntime (hr)",
                                breaks = c(0, 12, 24, 36, 48)) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          labs(x = "Infection rate (/min)", y = "Resistance Mutation Rate") +
          guides(fill = "none", shape = "none") +
          facet_grid(d ~ .,
                     labeller = labeller(
                       d = c("0" = "No nutrients returned by cell lysis",
                             "1" = "All nutrients returned by cell lysis"))) +
          theme(axis.title = element_text(size = 18),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 14)) +
          NULL)
  dev.off()
}

