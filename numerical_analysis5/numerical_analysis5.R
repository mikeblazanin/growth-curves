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

ode_deriv <- function(t, y, parms) {
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
  #dS/dt = u_S*S*N/k - aSP
  dY["S"] <- 
    (parms["u_S1"] * y["S1"] * y["N"]/parms["k"]
     - parms["a_S1"] * y["S1"] * y["P"])
  
  ##Calculate dI
  #dI_1/dt = aSP - I1 * nI / tau
  dY["I1"] <- 
    (parms["a_S1"] * y["S1"] * y["P"]
     - y["I1"] * parms["nI"] / parms["tau"])
  
  #dI_i/dt = nI / tau * I_(i-1) - nI / tau * I_i
  if(parms["nI"] > 1) {
    #Note S pop is dY[1], so I pops are dY[2: nI+1]
    dY[3:(parms["nI"]+1)] <- 
      (y[2:(parms["nI"])] * parms["nI"] / parms["tau"]
       - y[3:(parms["nI"]+1)] * parms["nI"] / parms["tau"])
  }
  
  ##Calculate dP
  #dP/dt = b * I_nI * nI / tau - aSP - z * aP*sum(I)
  dY["P"] <- 
    (parms["b"] * y[(parms["nI"]+1)] * parms["nI"] / parms["tau"]
     - parms["a_S1"] * y["S1"] * y["P"]
     - parms["z"] * parms["a_S1"] * y["S1"] * y["P"])
  
  ##Calcunate dN
  #dN/dt = -u_S*S*N/k + d * I_nI * nI / tau
  dY["N"] <-
    (-parms["u_S1"] * y["S1"] * y["N"]/parms["k"] +
      parms["d"] * y[(parms["nI"]+1)] * parms["nI"] / parms["tau"])
  
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
  times <- seq(from = 0, to = 1000, by = 1)
  yinit <- c("S1" = 10**6, "S2" = 0, "I1" = 0, "I2" = 0, 
             "P" = 10**4, "N" = (10**9 - 10**6))
  params <- c(u_S1 = 0.0179,
              k = 10**9,
              a_S1 = 10**-8,
              a_S2 = 0,
              tau = 100,
              b = 50,
              z = 1,
              f = 0,
              d = 0,
              v_a1 = 1, v_a2 = 1,
              g = 0, h = 0,
              warnings = 1, thresh_min_dens = 10**-100)
  
  test <- myTryCatch(as.data.frame(dede(y = yinit, times = times, func = derivs, 
                                        parms = params, hmax = 0.01)))
  test$S <- test$S1 + test$S2
  test$I <- test$I1 + test$I2
  test$B <- test$S + test$I
  test$pred <- params[["k"]]/
    (1+((params[["k"]] - yinit[["S1"]])/yinit[["S1"]])*
       exp(-params[["u_S1"]]*test$time))
  test2 <- tidyr::pivot_longer(test, cols = -c(time), 
                               names_to = "Pop", values_to = "Density")
  ggplot(data = filter(test2, Pop %in% c("S1", "S2", "I1", "I2", "P", "N")), 
         aes(x = time, y = Density, color = Pop)) +
    geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA)) +
    #geom_line(data = filter(test2, Pop == "pred"), lty = 2, color = "black") +
    NULL
}

## Define function for running simulations across many parameter values ----

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
  if (cntrs$j >= max_j | cntrs$k >= 15 | cntrs$j+cntrs$k >= 20) {
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
    
    #I1 and I2 at equil, and either S1 or N are at equil, we're done
    if (yout_list$value$I1[nrow(yout_list$value)] < equil_cutoff_dens &
        yout_list$value$I2[nrow(yout_list$value)] < equil_cutoff_dens &
        (yout_list$value$S1[nrow(yout_list$value)] < equil_cutoff_dens |
         yout_list$value$N[nrow(yout_list$value)] < equil_cutoff_dens)) {
      keep_running <- FALSE
      at_equil <- TRUE
      #S and N not at equil, need more time
    } else if (yout_list$value$S1[nrow(yout_list$value)] >= equil_cutoff_dens &
               yout_list$value$N[nrow(yout_list$value)] >= equil_cutoff_dens) {
      cntrs$j <- cntrs$j + 1
      #I1 or I2 not at equil (but S or N is because above check failed)
    } else if (yout_list$value$I1[nrow(yout_list$value)] >= equil_cutoff_dens |
               yout_list$value$I2[nrow(yout_list$value)] >= equil_cutoff_dens) {
      #If I1 or I2 are still changing, lengthen
      if(abs(yout_list$value$I1[nrow(yout_list$value)] - 
             yout_list$value$I1[nrow(yout_list$value)-1]) > 0 |
         abs(yout_list$value$I2[nrow(yout_list$value)] - 
             yout_list$value$I2[nrow(yout_list$value)-1]) > 0) {
        cntrs$j <- cntrs$j+1
        #If neither are changing, shorten step size
      } else {
        cntrs$I_only_cntr <- cntrs$I_only_cntr+1
        cntrs$k <- cntrs$k+1
      }
    } else {stop("check_equil found an unexpected case")}
  } else {stop("tryCatch failed, niether success, warning, nor error detected")}
  
  return(list(keep_running = keep_running, at_equil = at_equil, cntrs = cntrs))
}

run_sims <- function(u_S1vals,
                     kvals,
                     a_S1vals,
                     a_S2vals = NA,
                     tauvals,
                     bvals,
                     zvals = 0,
                     fvals = 0,
                     dvals = 1,
                     v_a1vals = 1,
                     v_a2vals = 1,
                     gvals = 0,
                     hvals = 0,
                     init_S1_dens_vals = 10**6,
                     init_S2_dens_vals = 0,
                     init_moi_vals = 10**-2,
                     init_N_dens_vals = NA,
                     equil_cutoff_dens = 0.1,
                     max_time = 48*60,
                     init_time = 100,
                     init_stepsize = 1,
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
  sim_vars <- list("u_S1" = u_S1vals, "k" = kvals,
                   "a_S1" = a_S1vals, "a_S2" = a_S2vals,
                   "tau" = tauvals, "b" = bvals,
                   "z" = zvals,
                   "f" = fvals, "d" = dvals,
                   "v_a1" = v_a1vals, "v_a2" = v_a2vals,
                   "g" = gvals, "h" = hvals,
                   "init_S1_dens" = init_S1_dens_vals, 
                   "init_S2_dens" = init_S2_dens_vals,
                   "init_moi" = init_moi_vals,
                   "init_N_dens" = init_N_dens_vals)
  
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
  param_combos$init_N_dens[is.na(param_combos$init_N_dens)] <-
    (param_combos$k[is.na(param_combos$init_N_dens)] -
       param_combos$init_S1_dens[is.na(param_combos$init_N_dens)] -
       param_combos$init_S2_dens[is.na(param_combos$init_N_dens)])
  #if a_S2 is na, a_S2 = a_S1
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
    yinit <- c(S1 = param_combos$init_S1_dens[i],
               S2 = param_combos$init_S2_dens[i],
               I1 = 0,
               I2 = 0,
               P = param_combos$init_S1_dens[i]*param_combos$init_moi[i],
               N = param_combos$init_N_dens[i])
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
        as.data.frame(dede(y = yinit, times = times, func = derivs, 
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

##Define file save/load wrapper for run_sims ----
run_sims_filewrapper <- function(name, dir = ".",
                                 read_file = TRUE, write_file = TRUE,
                                 ...) {
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
      paste(name, "_1.csv", sep = "") %in% list.files(dir)) {
    warning("Does not check if current param inputs are same as existing data")
    temp <- list(NULL, NULL, NULL)
    temp[[1]] <- read.csv(paste(dir, "/", name, "_1.csv", sep = ""),
                          stringsAsFactors = F)
    if (paste(name, "_2.csv", sep = "") %in% list.files(dir)) {
      temp[[2]] <- read.csv(paste(dir, "/", name, "_2.csv", sep = ""), 
                            stringsAsFactors = F)
    }
    if (paste(name, "_3.csv", sep = "") %in% list.files(dir)) {
      temp[[3]] <- read.csv(paste(dir, "/", name, "_3.csv", sep = ""), 
                            stringsAsFactors = F)
    }
  } else {
    #Run simulations (if files don't exist)
    if(all(lapply(list(...), class) == "list")) {
      #Run multiple sub-simulations
      temp_list <- vector(mode = "list", length(list(...))) 
      for(i in 1:length(temp_list)) {
        temp_list[[i]] <- do.call(run_sims, list(...)[[i]])
        
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
      temp <- run_sims(...)
    }
    
    #Save results so they can be re-loaded in future
    if (write_file) {
      #Save results so they can be re-loaded in future
      write.csv(temp[[1]], row.names = F,
                paste(dir, "/", name, "_1.csv", sep = ""))
      if (!is.null(temp[[2]])) {
        write.csv(temp[[2]], row.names = F,
                  paste(dir, "/", name, "_2.csv", sep = ""))
      }
      if (!is.null(temp[[3]])) {
        write.csv(temp[[3]], row.names = F,
                  paste(dir, "/", name, "_3.csv", sep = ""))
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

