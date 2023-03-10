##TODO: parallelize? (with Furrr and future r packages?)

## Import libraries ----
library(deSolve)
library(ggplot2)
library(dplyr)
library(gcplyr)

#Setwd
mywd_split <- strsplit(getwd(), split = "/")
if (mywd_split[[1]][length(mywd_split[[1]])] == "growth-curves") {
  dir.create("numerical_analysis4", showWarnings = FALSE)
  setwd("./numerical_analysis4/")
} else if (mywd_split[[1]][length(mywd_split[[1]])] != "numerical_analysis4") {
  stop("Not in correct root directory")
}

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                      "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)
                      
## Define derivatives function ----
derivs <- function(t, y, parms) {
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
  }

  #Create output vector
  dY <- c(S1 = 0, S2 = 0, I1 = 0, I2 = 0, P = 0, N = 0)
  
  #For all equations, let
  # afrac_t = (1 - f + f*(N/k)^v_a1)^v_a2
  # afrac_tau = (1 - f + f*(N(t-tau)/k)^v_a1)^v_a2
  #Note that when f > 1, afrac will be 0 when N >= k(1 - 1/f)
  afrac_t <- 
    max(0, 
        (1 - parms["f"] +
           parms["f"]*(y["N"]/parms["k"])**parms["v_a1"])**parms["v_a2"])
  if (t < parms["tau"]) {
    afrac_tau <- 0
  } else {
    afrac_tau <- 
      max(0,
          (1 - parms["f"] +
             parms["f"] * 
             (lagY[6]/parms["k"])**parms["v_a1"])**parms["v_a2"])
  }
  
  ##Calculate dS1 (growing subpopulation)
  #dS1/dt = u_S1*S1*N/k - afrac_t*a_S1*S1*P - h*u_S1*S1*(1-(g*N)/k)
  dY["S1"] <- (
    parms["u_S1"] * y["S1"] * y["N"]/parms["k"]
    - afrac_t * parms["a_S1"] * y["S1"] * y["P"]
    - parms["h"] * parms["u_S1"] * y["S1"] * (1 - (parms["g"]*y["N"])/parms["k"]))
  
  #Calculate dS2 (non-growing subpopulation)
  #dS2/dt = h*u_S1*S1*(1-(g*N)/k) - afrac_t*a_S2*S2*P
  dY["S2"] <- 
    (parms["h"] * parms["u_S1"] * y["S1"] * (1 - (parms["g"]*y["N"])/parms["k"])
     - afrac_t * parms["a_S2"] * y["S2"] * y["P"])
  
  ##Calculate dI1
  #dI1/dt = afrac_t * a_S1 * S1*P 
  #        - afrac_tau * a_S1 * S1(t-tau) * P(t-tau) 
  if (t < parms["tau"]) {
    dY["I1"] <- afrac_t * parms["a_S1"] * y["S1"] * y["P"] 
  } else {
    dY["I1"] <- 
      (afrac_t * parms["a_S1"] * y["S1"] * y["P"] 
       - afrac_tau * parms["a_S1"] * lagY[1] * lagY[5]
      )
  }
  
  ##Calculate dI2
  #dI2/dt = afrac_t * a_S2 * S2*P
  #        - afrac_tau * a_S2 * P(t-tau) * S2(t-tau)
  if (t < parms["tau"]) {
    dY["I2"] <- afrac_t * parms["a_S2"] * y["S2"] * y["P"] 
  } else {
    dY["I2"] <- 
      (afrac_t * parms["a_S2"] * y["S2"] * y["P"] 
       - afrac_tau * parms["a_S2"] * lagY[2] * lagY[5])
  }
  
  ##Calculate dP
  #dP/dt = b * afrac_tau * a_S1 * P(t-tau) * S1(t-tau) 
  #        b * afrac_tau * a_S2 * P(t-tau) * S2(t-tau)
  #        - afrac_t * a_S1 * S1*P
  #        - afrac_t * a_S2 * S2*P
  #        - z * afrac_t * a_S1 * I1 * P
  #        - z * afrac_t * a_S2 * I2 * P
  #        (factored in the code for efficiency)
  if (t < parms["tau"]) {
    dY["P"] <- 
      -afrac_t * y["P"] *
      ((parms["a_S1"] * y["S1"] + parms["a_S2"] * y["S2"])
       + parms["z"] * (parms["a_S1"] * y["I1"] + parms["a_S2"] * y["I2"]))
  } else {
    dY["P"] <- 
      (parms["b"]*afrac_tau*lagY[5]*
         (parms["a_S1"]*lagY[1] + parms["a_S2"]*lagY[2])
       -afrac_t * y["P"] *
         ((parms["a_S1"] * y["S1"] + parms["a_S2"] * y["S2"])
          + parms["z"] * (parms["a_S1"] * y["I1"] + parms["a_S2"] * y["I2"]))
      )
  }
  
  #Calculate dN
  #dN/dt = - u_S1*S1 * N/k
  #        + d*afrac_tau * a_S1 * S1(t-tau) * P(t-tau)  
  #        + d*afrac_tau * a_S2 * S2(t-tau) * P(t-tau) 
  #        (factored in code for efficiency)
  if (t < parms["tau"]) {
    dY["N"] <- - parms["u_S1"] * y["S1"] * y["N"]/parms["k"]
  } else {
    dY["N"] <- 
      (- parms["u_S1"] * y["S1"] * y["N"]/parms["k"]
       + parms["d"]*afrac_tau*lagY[5]*
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

## Run 1: phage traits ----
run1 <- run_sims_filewrapper(
  name = "run1",
  u_S1vals = signif(0.04*10**-0.35, 3),
  kvals = 10**9,
  a_S1vals = 10**seq(from = -12, to = -8, length.out = 5),
  tauvals = signif(10**seq(from = 1, to = 2, length.out = 5), 3),
  bvals = signif(5*10**seq(from = 0, to = 2, length.out = 5), 3),
  zvals = 1,
  fvals = 0,
  dvals = 0,
  v_a1vals = 1,
  v_a2vals = 1,
  init_S1_dens_vals = 10**6,
  init_moi_vals = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig1 <- run1[[1]]

ysum1 <- full_join(
  summarize(group_by(filter(ybig1, Pop == "B"),
                     uniq_run, u_S1, k, a_S1, a_S2,
                     tau, b, z, f, d, v_a1, v_a2, g, h,
                     init_S1_dens, init_S2_dens,
                     init_moi, init_N_dens),
            peak_dens = max(Density),
            peak_time = time[which.max(Density)],
            auc = auc(x = time, y = Density),
            extin_time_4 = 
              first_below(y = Density, x = time,
                          threshold = 10**4, return = "x"),
            run_time = max(time)),
  summarize(group_by(filter(ybig1, Pop == "P"),
                     uniq_run, u_S1, k, a_S1, a_S2,
                     tau, b, z, f, d, v_a1, v_a2, g, h,
                     init_S1_dens, init_S2_dens,
                     init_moi, init_N_dens),
            phage_final = Density[which.max(time)]))
ysum1 <- mutate(
  ysum1,
  extin_flag = ifelse(is.na(extin_time_4), "noextin",
                      ifelse(peak_dens >= 0.9*k, "neark", "none")),
  extin_time_4 = ifelse(is.na(extin_time_4), run_time, extin_time_4),
  phage_r = (log(phage_final)-log(init_moi*(init_S1_dens+init_S2_dens)))/
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
      width = 5, height = 5, units = "in", res = 150)
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


## Run 2: bact traits ----
run2 <- run_sims_filewrapper(
  name = "run2",
  u_S1vals = signif(0.04*10**seq(from = 0, to = -0.7, length.out = 5), 3),
  kvals = signif(10**c(8, 8.5, 9, 9.5, 10), 3),
  a_S1vals = 10**seq(from = -12, to = -8, length.out = 5),
  tauvals = signif(10**1.5, 3),
  bvals = 50,
  zvals = 1,
  fvals = 0,
  dvals = 0,
  v_a1vals = 1,
  v_a2vals = 1,
  init_S1_dens_vals = 10**6,
  init_moi_vals = 10**-2,
  equil_cutoff_dens = 0.1,
  init_time = 12*60,
  max_time = 48*60,
  init_stepsize = 5,
  print_info = TRUE, read_file = glob_read_files)

ybig2 <- run2[[1]]

ggplot(data = filter(ybig2, uniq_run %in% run2[[2]]$uniq_run,
                     Pop %in% c("S", "I", "P", "N")),
       aes(x = time/60, y = Density, color = Pop)) +
  geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA)) +
  facet_wrap(~uniq_run)

## Run 3: stationary phase behavior ----
run3 <- run_sims_filewrapper(
  name = "run3", read_file = glob_read_files,
  a = list(
    u_S1vals = signif(0.04*10**-0.35, 3),
    kvals = 10**9,
    a_S1vals = signif(10**seq(from = -12, to = -8, length.out = 9), 3),
    a_S2vals = 0,
    tauvals = signif(10**1.5, 3),
    bvals = 50,
    zvals = 1,
    fvals = c(0, 1, 2),
    dvals = 0,
    v_a1vals = 1,
    v_a2vals = 1,
    g = 0,
    h = c(0, 0.001, 0.01, 0.1),
    init_S1_dens_vals = 10**6,
    init_moi_vals = 10**-2,
    equil_cutoff_dens = 0.1,
    init_time = 12*60,
    max_time = 48*60,
    init_stepsize = 5,
    print_info = TRUE),
  b = list(
    u_S1vals = signif(0.04*10**-0.35, 3),
    kvals = 10**9,
    a_S1vals = signif(10**seq(from = -12, to = -8, length.out = 9), 3),
    a_S2vals = 0,
    tauvals = signif(10**1.5, 3),
    bvals = 50,
    zvals = 1,
    fvals = c(0, 1, 2),
    dvals = 0,
    v_a1vals = 1,
    v_a2vals = 1,
    g = 1,
    h = c(0.01, 0.1, 1),
    init_S1_dens_vals = 10**6,
    init_moi_vals = 10**-2,
    equil_cutoff_dens = 0.1,
    init_time = 12*60,
    max_time = 48*60,
    init_stepsize = 5,
    print_info = TRUE)
)

ybig3 <- run3[[1]]

ysum3 <- full_join(
  summarize(group_by(filter(ybig3, Pop == "B"),
                     uniq_run, u_S1, k, a_S1, a_S2,
                     tau, b, z, f, d, v_a1, v_a2, g, h,
                     init_S1_dens, init_S2_dens,
                     init_moi, init_N_dens),
            peak_dens = max(Density),
            peak_time = time[which.max(Density)],
            auc = auc(x = time, y = Density),
            extin_time_4 = 
              first_below(y = Density, x = time,
                          threshold = 10**4, return = "x"),
            final_dens = Density[which.max(time)]),
  summarize(group_by(filter(ybig3, Pop == "P"),
                     uniq_run, u_S1, k, a_S1, a_S2,
                     tau, b, z, f, d, v_a1, v_a2, g, h,
                     init_S1_dens, init_S2_dens,
                     init_moi, init_N_dens),
            phage_final = Density[which.max(time)]))
ysum3 <- mutate(ysum3,
                final_dens = ifelse(final_dens <= 1, 0, final_dens))

#For the purposes of plotting, extend all runs to same endtime
# and set all densities below 0 to 0
ybig3 <- mutate(group_by(ybig3, uniq_run),
                time = ifelse(time == max(time), max(ybig3$time), time),
                Density = ifelse(Density < 1, 0, Density))

dir.create("run3_dens_curves", showWarnings = FALSE)
if (glob_make_curveplots) {
  for (run in unique(ybig3$uniq_run)) {
    png(paste("./run3_dens_curves/", run, ".png", sep = ""),
         width = 5, height = 5, units = "in", res = 150)
    print(
      ggplot(data = filter(ybig3, uniq_run == run,
                           Pop %in% c("S1", "S2", "I1", "I2", "P", "N")),
             aes(x = time/60, y = Density, color = as.factor(Pop))) + 
        geom_line(lwd = 1.5, alpha = 1) + 
        geom_line(data = filter(ybig3, uniq_run == run, Pop == "B"),
                  color = "black", lwd = 0.6, alpha = 1) +
        scale_y_continuous(trans = "log10", limits = c(1, NA)) +
        scale_color_manual(limits = c("S1", "S2", "I1", "I2", "P", "N"),
                           values = my_cols[c(2, 5, 1, 6, 3, 4)]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              title = element_text(size = 9)) +
        NULL
    )
    dev.off()
  }
}

if(glob_make_statplots) {
  png("./statplots/h_Bcurves.png", width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(data = filter(ybig3, Pop == "B", f == 0, g == 0),
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
  
  png("./statplots/h_Bcurves_g1.png", width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(data = filter(ybig3, Pop == "B", f == 0, h == 0 | g == 1),
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
  
  png("./statplots/f_Bcurves.png", width = 4, height = 7,
      units = "in", res = 300)
  print(ggplot(data = filter(ybig3, Pop == "B", h == 0, g == 0),
               aes(x = time/60, y = Density, group = uniq_run)) +
          geom_line() +
          facet_grid(a_S1 ~ f) +
          scale_y_continuous(trans = "log10", limits = c(1, NA)) +
          coord_cartesian(xlim = c(NA, 24)) +
          theme_bw() +
          labs(x = "Time (hr)", y = "Density (cfu/mL)",
               subtitle = "Reduction in infection rate at N=0") +
          NULL)
  dev.off()
  
  png("./statplots/gh_peakdens_finaldens.png", width = 6, height = 4,
      units = "in", res = 300)
  print(ggplot(data = filter(ysum3, f == 0),
               aes(x = peak_dens, y = final_dens, color = log10(a_S1))) +
          facet_grid(g ~ h,
                     labeller = labeller(g = c("1" = "Inverse logistic\ntransitions",
                                               "0" = "Linear transitions"))) +
          geom_point(alpha = 0.5, size = 2.5) +
          scale_x_log10() + scale_y_log10() +
          scale_color_viridis_c(end = 0.95, name = "log10(infection rate)") +
          labs(x = "Peak Density (cfu/mL)", y = "Final Density (cfu/mL)",
               subtitle = "Resistance Transition Rate") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          NULL)
  dev.off()
  
  ##Make figures to support the point that h has little effect
  ## on other observed patterns, like peak dens-peak time,
  ## peak time at all, peat time - extin time, etc
  
  
}


ggplot(data = ysum3,
       aes(x = peak_time, y = peak_dens)) +
  geom_point(aes(color = as.factor(a_S1))) +
  facet_grid(g*h ~ f) +
  geom_line(data = data.frame(x = seq(from = 0, to = 3000),
                              y = logis_func(S_0 = 10**6,
                                             u_S = signif(0.04*10**-0.35, 3),
                                             k = 10**9,
                                             times = seq(from = 0, to = 3000))),
            aes(x = x, y = y), lty = 2)

ggplot(data = filter(ysum3, f == 0, g == 0),
       aes(y = peak_time, x = h)) +
  facet_wrap(~a_S1) +
  geom_point(alpha = 0.5) +
  scale_x_log10() + scale_y_log10()


#Findings:
  # rate of transition determines how low the pop is after S1 extin
  # a_t modulation must have only a very narrow window where it
  #  produces a partial drop-off. Instead mostly what we see is
  #  either extinction or lack of drop entirely
  # production of resistant cells does produce the partial drop-off
  #  rate and model (logistic vs linear) both appear to control how
  #  low the drop off is. Can't see any other differences between
  #  the two modes
  # drop off from peak is at the same time regardless of g or h,
  #  suggesting that peak time, peak density, and extin time (assuming
  #  the pop drops below that extin threshold) should be unaltered
  #  TODO: why doesn't this hold on the peak dens vs peak time plot
  # logistic transition to resistance leads to phage infectivity
  #  having a larger effect on final bact density bc the resis pop
  #  mostly accumulates right at the end. W/ constant transition rate
  #  phage infectivity has smaller effect on final density
  # f really is just a binary effect - either the population doesn't
  #  drop at all, or the population drops entirely. Which occurs
  #  depends on how high the pop gets before phage would start
  #  killing (and therefore is basically determined by phage
  #  infectivity and bact growth). It can have some small small effect
  #  on the rate of extinction but really tiny effects overall
  # h has only a very weak effect (really basically no effect,
  #  in comparison to even small changes in a) on peak density or
  #  on peak time


##Run 4: test of metrics across dift bact ----
run4 <- run_sims_filewrapper(
  name = "run4", read_file = glob_read_files,
  a = list(
    u_S1vals = signif(0.04*10**seq(from=-0.175, to=-0.525, length.out = 3), 3),
    kvals = signif(10**seq(from = 8.75, to = 9.25, length.out = 3), 3),
    a_S1vals = signif(10**seq(from = -12, to = -8, length.out = 5), 3),
    a_S2vals = 0,
    tauvals = 31.6,
    bvals = 50,
    zvals = 1,
    fvals = 0,
    dvals = 0,
    v_a1vals = 1,
    v_a2vals = 1,
    g = 0,
    h = c(0, 0.1),
    init_S1_dens_vals = 10**6,
    init_moi_vals = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
    equil_cutoff_dens = 0.1,
    init_time = 12*60,
    max_time = 48*60,
    init_stepsize = 5,
    print_info = TRUE),
  b = list(
    u_S1vals = signif(0.04*10**seq(from=-0.175, to=-0.525, length.out = 3), 3),
    kvals = signif(10**seq(from = 8.75, to = 9.25, length.out = 3), 3),
    a_S1vals = 0,
    a_S2vals = 0,
    tauvals = 31.6,
    bvals = 50,
    zvals = 1,
    fvals = 0,
    dvals = 0,
    v_a1vals = 1,
    v_a2vals = 1,
    g = 0,
    h = c(0, 0.1),
    init_S1_dens_vals = 10**6,
    init_moi_vals = 0,
    equil_cutoff_dens = 0.1,
    init_time = 48*60,
    max_time = 48*60,
    init_stepsize = 5,
    print_info = TRUE))

ybig4 <- run4[[1]]

#For the purposes of analyses, set all Dens below 0 to 0
# and set all runs to end at same time
ybig4 <- mutate(group_by(ybig4, uniq_run), 
                time = ifelse(time == max(time), max(ybig4$time), time),
                Density = ifelse(Density < 0, 0, Density))

ysum4 <- 
  summarize(group_by(filter(ybig4, Pop == "B"),
                     uniq_run, u_S1, k, a_S1, a_S2,
                     tau, b, z, f, d, v_a1, v_a2, g, h,
                     init_S1_dens, init_S2_dens,
                     init_moi, init_N_dens),
            peak_dens = max(Density),
            peak_time = time[which.max(Density)],
            auc = auc(x = time, y = Density),
            extin_time_4 = 
              first_below(y = Density, x = time,
                          threshold = 10**4, return = "x"))

#Relative AUC (incl over dift MOIs)
ysum4 <- mutate(group_by(ysum4, u_S1, k, z, f, d, v_a1, v_a2, g, h,
                         init_S1_dens, init_S2_dens, init_N_dens),
                rel_auc = auc/auc[init_moi == 0])

ggplot(data = filter(ysum4, init_moi %in% c(0, 0.01)), 
       aes(x = paste(u_S1, k), y = rel_auc)) +
  geom_point(aes(color = as.factor(init_moi))) +
  #scale_x_log10() +
  scale_y_log10() +
  facet_grid(a_S1 ~ h) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filter(ysum4, init_moi %in% c(0, 0.01)), 
       aes(x = paste(u_S1, k), y = auc)) +
  geom_point(aes(color = as.factor(init_moi))) +
  #scale_x_log10() +
  scale_y_log10() +
  facet_grid(a_S1 ~ h) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


temp <- mutate(filter(ysum4, a_S1 != 0),
               u_S1 = as.factor(u_S1),
               k = as.factor(k),
               h = as.factor(h),
               log_a_S1 = log10(a_S1),
               a_S1 = as.factor(a_S1),
               bact = paste(u_S1, k, h))

library(lme4)
summary(lmer(data = temp,
             log10(auc) ~ a_S1 * (1|bact)))
summary(lm(data = temp,
           log10(auc) ~ a_S1))

summary(lmer(data = temp,
             log10(rel_auc) ~ a_S1 * (1|bact)))
summary(lm(data = temp,
           log10(rel_auc) ~ a_S1))

##when we haven't normalized auc, we want to do a lm of auc ~ bact * phage or not
## thus the "normalization" happens within the lm where it's an assumption
## of additivity or whatever the lm assumes
##when we do normalize auc, we just do an lm of norm auc ~ bact
## in this way, we've already done the normalization ahead of the model
## and it should give us a more biologically meaningful normalization
##this is for the case where we're comparing one phage across multiple bacteria
## if we're imagining our phage has the same a on all the bacteria,
## then the response should vary less with bacteria when we normalize by relative
## auc than when we normalize within the lm
##So what we want to look for is less variance being explained by the
## bacteria? Or rather, since we've already normalized, less residual
## variance
## 

temp <- mutate(
  filter(ysum4, a_S1 %in% c(0, 10**-9, 10**-10), h == 0, 
         init_moi %in% c(0, 0.01)),
  u_S1 = as.factor(u_S1),
  k = as.factor(k),
  h = as.factor(h),
  log_a_S1 = log10(a_S1),
  a_S1 = as.factor(a_S1),
  bact = paste(u_S1, k, h))

ggplot(data = temp,
       aes(x = bact, y = log10(auc), color = as.factor(init_moi),
           shape = as.factor(a_S1))) +
  geom_point() +
  #scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = temp,
       aes(x = bact, y = log10(rel_auc), color = as.factor(init_moi),
           shape = as.factor(a_S1))) +
  geom_point() +
  #scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Need to dup the init_moi == 0 data across all a_S1 vals
temp2 <- temp
mini <- filter(temp2, init_moi == 0)
mini <- rep(mini, length(unique(temp2$a_S1)))
temp2 <- rbind(temp2,
               



lmer(data = temp,
     log10(auc) ~ a_S1 * init_moi * (1|bact))

lmer(data = filter(temp, init_moi != 0),
   log10(rel_auc) ~ a_S1 * (1|bact))



summary(lmer(data = temp,
     log10(auc) ~ log_a_S1 * (log_a_S1|u_S1) * (log_a_S1|k) * (log_a_S1|h)))
summary(lmer(data = temp,
     log10(rel_auc) ~ log_a_S1 * (log_a_S1|u_S1) * (log_a_S1|k) * (log_a_S1|h)))



ggplot(data = filter(temp, a_S1 != 0),
       aes(x = log10(a_S1), y = log10(auc))) +
  geom_point()
ggplot(data = filter(temp, a_S1 != 0),
       aes(x = log10(a_S1), y = log10(rel_auc))) +
  geom_point()



#All OD into PCA
#OD of control minus phage_added into PCA


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

temp <- filter(ybig4, Pop == "B")
ybig4_new <- interp_data(df = temp, 
                         x = "time", y = "Density",
                         subset_by = temp$uniq_run)
ybig4_new <- mutate(group_by(ybig4_new, u_S1, k, h, time),
                    Dens_norm = Density - Density[init_moi == 0])

dir.create("test_ybig4", showWarnings = FALSE)                    
for (run in unique(ybig4$uniq_run)) {
  png(width = 4, height = 6, units = "in", res = 150,
      filename = paste0("./test_ybig4/", run, ".png"))
  print(cowplot::plot_grid(ncol = 1,
    ggplot(data = filter(ybig4_new, uniq_run == run),
           aes(x = time/60, y = Density)) +
      geom_point() +
      scale_y_log10(),
    ggplot(data = filter(ybig4_new, uniq_run == run),
           aes(x = time/60, y = Dens_norm)) +
      geom_point()
  ))
  dev.off()
}

ybig4_wide <- tidyr::pivot_wider(ybig4_new,
                                 id_cols = !Dens_norm,
                                names_from = time,
                                names_prefix = "time_",
                                values_from = Density)
ybig4_wide_norm <- tidyr::pivot_wider(ybig4_new,
                                      id_cols = !Density,
                                      names_from = time,
                                      names_prefix = "time_",
                                      values_from = Dens_norm)

#Define function to make chi-square quantile plots 
# to test for multivariate normality of data or residuals
# (credit to Jonathan Reuning-Scherer)
CSQPlot<-function(vars,label="Chi-Square Quantile Plot"){
  #usually, vars is xxx$residuals or data from one group and label is for plot
  x<-cov(scale(vars),use="pairwise.complete.obs")
  squares<-sort(diag(as.matrix(scale(vars))%*%solve(x)%*%as.matrix(t(scale(vars)))))
  quantiles<-quantile(squares)
  hspr<-quantiles[4]-quantiles[2]
  cumprob<-c(1:length(vars[,1]))/length(vars[,1])-1/(2*length(vars[,1]))
  degf<-dim(x)[1]
  quants<-qchisq(cumprob,df=degf)
  gval<-(quants**(-1+degf/2))/(exp(quants/2)*gamma(degf/2)*(sqrt(2)**degf))
  scale<-hspr / (qchisq(.75,degf)-qchisq(.25,degf))
  se<-(scale/gval)*sqrt(cumprob*(1-cumprob)/length(squares))
  lower<-quants-2*se
  upper<-quants+2*se
  
  plot(quants,squares,col='red',pch=19,cex=1.2,xlab="Chi-Square Quantiles",
       ylab=label,main=paste("Chi-Square Quantiles for",label),ylim=range(upper,lower, squares) , xlim=range(c(0,quants)))
  lines(c(0,100),c(0,100),col=1)
  lines(quants,upper,col="blue",lty=2,lwd=2)
  lines(quants,lower,col="blue",lty=2,lwd=2)
  legend(0,range(upper,lower)[2]*.9,c("Data","95% Conf Limits"),lty=c(0,2),col=c("red","blue"),lwd=c(2,2),
         pch=c(19,NA))
}

#Make CSQ plot for multivariate normality
CSQPlot(ybig4_wide[, grep("time_", colnames(ybig4_wide))[-1]])

mypca <- prcomp(ybig4_wide[, grep("time_", colnames(ybig4_wide))[-1]],
                center = TRUE, scale = TRUE, retx = TRUE)
mypcanorm <- 
  prcomp(ybig4_wide_norm[, grep("time_", colnames(ybig4_wide_norm))[-1]],
         center = TRUE, scale = TRUE, retx = TRUE)

#Merge pca with data
ybig4_wide <- cbind(ybig4_wide, as.data.frame(mypca$x))
ybig4_wide_norm <- cbind(ybig4_wide_norm, as.data.frame(mypcanorm$x))

summary(mypca)

ggplot(ybig4_wide, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = as.factor(a_S1))) +
  # geom_segment(data = as.data.frame(isol_pca_7x$rotation),
  #              aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
  #              arrow = arrow(length = unit(0.02, "npc")),
  #              alpha = .5, lwd = 2, color = "red4") +
  # ggrepel::geom_text_repel(data = as.data.frame(isol_pca_7x$rotation),
  #                          aes(x = arrow_len*PC1, y = arrow_len*PC2,
  #                              label = row.names(isol_pca_7x$rotation)),
  #                          size = 8, alpha = .8, color = "red4", seed = 8,
  #                          min.segment.length = unit(1, "native")) +
  theme_bw() +
  # labs(x = paste("PC1 (", 
  #                round((100*((isol_pca_7x$sdev)**2)/
  #                         sum((isol_pca_7x$sdev)**2))[1], 1),
  #                "%)", sep = ""),
  #      y = paste("PC2 (", 
  #                round((100*((isol_pca_7x$sdev)**2)/
  #                         sum((isol_pca_7x$sdev)**2))[2], 1),
  #                "%)", sep = "")) +
  NULL

ggplot(ybig4_wide_norm, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = as.factor(u_S1), shape = as.factor(k))) +
  # geom_segment(data = as.data.frame(isol_pca_7x$rotation),
  #              aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
  #              arrow = arrow(length = unit(0.02, "npc")),
  #              alpha = .5, lwd = 2, color = "red4") +
  # ggrepel::geom_text_repel(data = as.data.frame(isol_pca_7x$rotation),
  #                          aes(x = arrow_len*PC1, y = arrow_len*PC2,
  #                              label = row.names(isol_pca_7x$rotation)),
  #                          size = 8, alpha = .8, color = "red4", seed = 8,
  #                          min.segment.length = unit(1, "native")) +
  theme_bw() +
  # labs(x = paste("PC1 (", 
  #                round((100*((isol_pca_7x$sdev)**2)/
  #                         sum((isol_pca_7x$sdev)**2))[1], 1),
  #                "%)", sep = ""),
  #      y = paste("PC2 (", 
  #                round((100*((isol_pca_7x$sdev)**2)/
  #                         sum((isol_pca_7x$sdev)**2))[2], 1),
  #                "%)", sep = "")) +
  NULL



temp <- filter(ybig4_wide_norm, k == 10**9, u_S1 == 0.0179,
               init_moi == 0.01)
mypca_tmp <- prcomp(temp[, grep("time_", colnames(temp))[-1]],
  center = TRUE, scale = TRUE, retx = TRUE)

temp <- cbind(temp, as.data.frame(mypca_tmp$x))

ggplot(temp, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = as.factor(a_S1), shape = as.factor(h)))
