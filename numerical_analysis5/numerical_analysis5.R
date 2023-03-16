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
    geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA)) +
    #geom_line(data = filter(test2, Pop == "pred"), lty = 2, color = "black") +
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
                           "h", "g1", "g2",
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
      y_noequil <- ybig[min(which(ybig$uniq_run == run)), 1:21]
    } else {
      y_noequil <- rbind(y_noequil, ybig[min(which(ybig$uniq_run == run)), 1:21])
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
               rep(0, nI),
               P = param_combos$init_S1[i]*param_combos$init_moi[i],
               N = param_combos$init_N[i])
    names(yinit)[3:(3+nI-1)] <- paste0("I", 1:nI)
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
        as.data.frame(ode(y = yinit, times = times, func = deriv_dede, 
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
      y_noequil <- ybig[min(which(ybig$uniq_run == run)), 1:21]
    } else {
      y_noequil <- rbind(y_noequil, ybig[min(which(ybig$uniq_run == run)), 1:21])
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
  init_S1 = 10**6,
  init_moi = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig1 <- run1[[1]]

#Set below 0 values to 0
ybig1 <- mutate(ybig1,
                Density = ifelse(Density < 0, 0, Density))

ysum1 <- full_join(
  summarize(group_by(filter(ybig1, Pop == "B"),
                     uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                     tau, b, z, f_a, f_b, d, h, g1, g2,
                     init_S1, init_S2, init_moi, init_N, equil),
            peak_dens = max(Density),
            peak_time = time[which.max(Density)],
            auc = auc(x = time, y = Density),
            extin_time_4 = 
              first_below(y = Density, x = time,
                          threshold = 10**4, return = "x"),
            run_time = max(time)),
  summarize(group_by(filter(ybig1, Pop == "P"),
                     uniq_run, u_S1, u_S2, k, a_S1, a_S2,
                     tau, b, z, f_a, f_b, d, h, g1, g2,
                     init_S1, init_S2, init_moi, init_N, equil),
            phage_final = Density[which.max(time)]))
ysum1 <- mutate(
  ysum1,
  extin_flag = ifelse(is.na(extin_time_4), "noextin",
                      ifelse(peak_dens >= 0.9*k, "neark", "none")),
  extin_time_4 = ifelse(is.na(extin_time_4), run_time, extin_time_4),
  phage_r = (log(phage_final)-log(init_moi*(init_S1+init_S2)))/
    extin_time_4)


# Run 1: B curves & stat v stat plots ----
dir.create("./statplots", showWarnings = FALSE)
if(glob_make_statplots) {
  png("./statplots/run1_peakdens_peaktime.png",
      width = 5, height = 5, units = "in", res = 150)
  print(
    ggplot(data = ysum1,
           aes(x = peak_time/60, y = peak_dens)) +
      geom_point(aes(shape = extin_flag)) +
      scale_shape_manual(breaks = c("none", "neark", "noextin"),
                         values = c(16, 4, 3)) +
      guides(shape = "none") +    theme_bw() +
      labs(x = "Peak Time (hr)", y = "Peak Density (cfu/mL)") +
      geom_line(data = data.frame(x = 0:1440,
                                  y = logis_func(S_0 = 10**6, u_S = 0.0179,
                                                 k = 10**9, times = 0:1440)),
                aes(x = x/60, y = y), lty = 2)
  )
  dev.off()
  
  png("./statplots/run1_peakdens_peaktime_subset.png",
      width = 5, height = 5, units = "in", res = 150)
  print(
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = peak_time/60, y = peak_dens)) +
      geom_point() +
      theme_bw() +
      labs(x = "Peak Time (hr)", y = "Peak Density (cfu/mL)") +
      geom_line(data = data.frame(x = 0:900,
                                  y = logis_func(S_0 = 10**6, u_S = 0.0179,
                                                 k = 10**9, times = 0:900)),
                aes(x = x/60, y = y), lty = 2)
    + NULL)
  dev.off()
  
  png("./statplots/run1_Bcurves.png",
      width = 5, height = 4, units = "in", res = 150)
  print(
    ggplot(data = filter(ybig1, Pop == "B", b == 50, tau == 31.6),
           aes(x = time/60, y = Density)) +
      geom_line(aes(color = as.factor(a_S1), group = interaction(a_S1, b, tau)),
                lwd = 1.5) +
      theme_bw() +
      labs(x = "Time (hr)", y = "Density (cfu/mL)") +
      scale_x_continuous(limits = c(NA, 24)) +
      geom_line(data = data.frame(x = 0:1440,
                                  y = logis_func(S_0 = 10**6, u_S = 0.0179,
                                                 k = 10**9, times = 0:1440)),
                aes(x = x/60, y = y), lty = 2) +
      scale_color_manual(values = colorRampPalette(c("gray70", "darkblue"))(5),
                         name = "Infection rate\n(/cfu/pfu/min)")
    + NULL)
  dev.off()
  
  png("./statplots/run1_extintime_peaktime.png",
      width = 5, height = 5, units = "in", res = 150)
  print(
    ggplot(data = ysum1,
           aes(x = peak_time/60, y = extin_time_4/60)) +
      geom_point(aes(shape = extin_flag)) +
      scale_shape_manual(breaks = c("none", "neark", "noextin"),
                         values = c(16, 4, 3)) +
      guides(shape = "none") +
      theme_bw() +
      geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
      labs(x = "Peak Time (hr)", y = "Extinction Time (hr)")
    + NULL)
  dev.off()
  
  png("./statplots/run1_extintime_peaktime_subset.png",
      width = 5, height = 5, units = "in", res = 150)
  print(
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = peak_time/60, y = extin_time_4/60)) +
      geom_point() +
      theme_bw() +
      geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
      labs(x = "Peak Time (hr)", y = "Extinction Time (hr)")
    + NULL)
  dev.off()
  
  png("./statplots/run1_auc_peaktime.png",
      width = 5, height = 5, units = "in", res = 150)
  print(
    ggplot(data = ysum1,
           aes(x = peak_time/60, y = auc/60)) +
      geom_point(aes(shape = extin_flag)) +
      scale_shape_manual(breaks = c("none", "neark", "noextin"),
                         values = c(16, 4, 3)) +
      guides(shape = "none") +
      theme_bw() +
      geom_line(data = data.frame(
        x = 0:1080,
        y = logis_def_integral(S_0 = 10**6, u_S = 0.0179,
                               k = 10**9, times = 0:1080)),
        aes(x = x/60, y = y/60), lty = 2) +
      scale_y_log10() +
      labs(x = "Peak Time (hr)", y = "Area Under the Curve (hr cfu/mL)")
    + NULL)
  dev.off()
  
  png("./statplots/run1_auc_peaktime_subset.png",
      width = 5, height = 5, units = "in", res = 150)
  print(
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = peak_time/60, y = auc/60)) +
      geom_point() +
      theme_bw() +
      geom_line(data = data.frame(
        x = 0:540,
        y = logis_def_integral(S_0 = 10**6, u_S = 0.0179,
                               k = 10**9, times = 0:540)),
        aes(x = x/60, y = y/60), lty = 2) +
      scale_y_log10() +
      labs(x = "Peak Time (hr)", y = "Area Under the Curve (hr cfu/mL)")
    + NULL)
  dev.off()
}

#Run 1: contour plots ----
if (glob_make_statplots) {
  png("./statplots/run1_maxtime_a_b_contour.png", width = 5, height = 4,
      units = "in", res = 300)
  print(
    ggplot(data = filter(ysum1, b == 50),
           aes(x = a_S1, y = tau)) +
      geom_contour_filled(aes(z = peak_time/60), alpha = 0.5) +
      geom_point(aes(color = peak_time/60, shape = extin_flag),
                 size = 3) +
      scale_color_viridis_c(name = "Peak time (hr)",
                            breaks = c(4, 8, 12)) +
      scale_shape_manual(breaks = c("neark", "noextin", "none"), 
                         values = c(4, 4, 16)) +
      scale_y_continuous(trans = "log10", breaks = c(16, 40, 100)) +
      scale_x_continuous(trans = "log10") +
      xlab("Infection rate (/min)") +
      ylab("Lysis time (min)") +
      guides(fill = "none", shape = "none") +
      NULL)
  dev.off()
  
  p1 <- ggplot(data = filter(ysum1), aes(x = a_S1, y = tau)) +
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
  
  p2 <- ggplot(data = filter(ysum1), aes(x = a_S1, y = b)) +
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
  
  p3 <- ggplot(data = filter(ysum1), aes(x = tau, y = b)) +
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
  
  png("./statplots/run1_maxtime_contour_all.png", width = 6, height = 6,
      units = "in", res = 300)
  print(cowplot::plot_grid(
    cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                       p2 + theme(legend.position = "none"), 
                       p3 + theme(legend.position = "none"),
                       ncol = 1),
    cowplot::get_legend(p1),
    rel_widths = c(1, .2),
    ncol = 2))
  dev.off()
}

# Run 1: phage growth plots ----
if (glob_make_statplots) {
  png("./statplots/phager_extintime_subset.png", width = 5, height = 4,
      units = "in", res = 300)
  print(
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = extin_time_4/60, y = phage_r*60, color = as.factor(b))) +
      geom_point() +
      scale_color_viridis_d(end = 0.95, name = "Burst Size") +
      scale_x_log10() + 
      scale_y_log10() +
      labs(x = "Extinction time (hr)", 
           y = "Phage Aggregate Growth Rate (e-fold/hour)") +
      theme_bw() +
      NULL)
  dev.off()
  
  png("./statplots/phager_extintime.png", width = 5, height = 4,
      units = "in", res = 300)
  print(
    ggplot(data = ysum1,
           aes(x = extin_time_4/60, y = phage_r*60, color = as.factor(b),
               shape = extin_flag)) +
      geom_point() +
      scale_color_viridis_d(end = 0.95, name = "Burst Size") +
      scale_shape_manual(breaks = c("none", "neark", "noextin"),
                         values = c(16, 4, 3)) +
      scale_x_log10() + 
      scale_y_log10() +
      labs(x = "Extinction time (hr)", 
           y = "Phage Aggregate Growth Rate (e-fold/hour)") +
      theme_bw() +
      guides(shape = "none") +
      NULL)
  dev.off()
  
  png("./statplots/phagefinal_peakdens.png", width = 5, height = 4,
      units = "in", res = 300)
  print(
    ggplot(data = ysum1,
           aes(x = peak_dens, y = phage_final, shape = extin_flag,
               color = as.factor(b))) +
      geom_point(size = 2) +
      scale_y_log10() + scale_x_log10() +
      scale_color_viridis_d(end = 0.95, name = "Burst Size") +
      scale_shape_manual(breaks = c("none", "neark", "noextin"),
                         values = c(16, 4, 3)) +
      labs(x = "Peak Bacterial Density (cfu/mL)", 
           y = "Final Phage Density (pfu/mL)") +
      guides(shape = "none") +
      geom_line(aes(y = peak_dens*b)) +
      theme_bw() +
      NULL)
  dev.off()
  
  png("./statplots/phagefinal_peakdens_subset.png", width = 5, height = 4,
      units = "in", res = 300)
  print(
    ggplot(data = filter(ysum1, extin_flag == "none"),
           aes(x = peak_dens, y = phage_final, color = as.factor(b))) +
      geom_point(size = 2) +
      scale_y_log10() + scale_x_log10() +
      scale_color_viridis_d(end = 0.95, name = "Burst Size") +
      labs(x = "Peak Bacterial Density (cfu/mL)", 
           y = "Final Phage Density (pfu/mL)") +
      geom_line(aes(y = peak_dens*b)) +
      theme_bw() +
      NULL)
  dev.off()
}

