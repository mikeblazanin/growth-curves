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

##Some simple runs ----
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
times <- seq(from = 0, to = 1000, by = 1)
test <- as.data.frame(dede(y = yinit, times = times, func = derivs, 
                           parms = params))
test2 <- tidyr::pivot_longer(test, cols = -c(time), 
                             names_to = "Pop", values_to = "Density")
ggplot(data = test2, aes(x = time, y = Density, color = Pop)) +
  geom_line() + scale_y_continuous(trans = "log10", limits = c(1, NA))

## Define function for running simulations across many parameter values ----

#sub-function for checking for equilibrium
check_equil <- function(yout_list, cntrs, fixed_time, equil_cutoff_dens) {
  #Returns: list(keep_running = TRUE/FALSE,
  #              at_equil = TRUE/FALSE/NA (NA only when fixed_time = TRUE),
  #              cntrs = list([other entries],
  #                           i_only_pos_times = new i_only_pos_times value,
  #                           j = new j value, 
  #                           k = new k value)
  #              )
  
  #Infinite loop prevention check (j = 10 is 24 hrs)
  if (cntrs$j >= 10 | k >= 15 | cntrs$j+k >= 20) {
    return(list(keep_running = FALSE, at_equil = FALSE, cntrs))
  }
  
  #If fixed time, don't check for equil
  if(fixed_time) {
    return(list(keep_running = FALSE, at_equil = NA, cntrs))
  }
  
  #If there was an error, increase k by 1 and re-run
  if(!is.null(yout_list$error)) {
    cntrs$k <- cntrs$k+1
    return(list(keep_running = TRUE, at_equil = NA, cntrs))
  #If there was a warning, could be several causes, so we
  # generally just halve step size and increase length
  } else if (!is.null(yout_list$warning)) {
    cntrs$j <- cntrs$j+1
    cntrs$k <- cntrs$k+2
    return(list(keep_running = TRUE, at_equil = NA, cntrs))
  #If it was successful, check for equilibrium
  } else if (is.null(yout_list$warning) & is.null(yout_list$error)) {
    #First drop all rows with nan
    yout_list$value <- 
      yout_list$value[apply(X = yout_list$value, MARGIN = 1,
                            FUN = function(x) {all(!is.nan(x))}), ]

    #S and I both at equil, we're done
    if (yout_list$value$S[nrow(yout_list$value)] < equil_cutoff_dens & 
        yout_list$value$I[nrow(yout_list$value)] < equil_cutoff_dens) {
      return(list(keep_running = FALSE, at_equil = TRUE, cntrs))
    #S not at equil, need more time
    } else if (yout_list$value$S[nrow(yout_list$value)] >= equil_cutoff_dens) {
      cntrs$j <- cntrs$j + 1
      return(list(keep_running = TRUE, at_equil = FALSE, cntrs))
    #I not at equil (but S is because above check failed),
    #   first we'll lengthen the simulation
    #    (to make sure it was long enough to catch the last burst)
    #   then we'll start shrinking our step size
    } else if (yout_list$value$I[nrow(yout_list$value)] >= equil_cutoff_dens) {
      if (cntrs$i_only_pos_times < 1) {
        cntrs$i_only_pos_times <- cntrs$i_only_pos_times+1
        cntrs$j <- cntrs$j+1
        return(list(keep_running = TRUE, at_equil = FALSE, cntrs))
      } else {
        cntrs$k <- cntrs$k+1
        return(list(keep_running = TRUE, at_equil = FALSE, cntrs))
      }
    }
  } else {stop("tryCatch failed, niether success, warning, nor error detected")}
}

run_sims <- function(u_Svals,
                     u_Rvals = 0,
                     kvals,
                     avals,
                     tauvals,
                     bvals,
                     zvals = 0,
                     mvals = 0,
                     fvals = 0,
                     dvals = 1,
                     v_a1vals = 1,
                     v_a2vals = 1,
                     init_S_dens_vals = 10**6,
                     init_R_dens_vals = 0,
                     init_moi_vals = 10**-2,
                     init_N_dens_vals = NA,
                     equil_cutoff_dens = 0.1,
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
  
  num_pops <- 7 #placeholder for when derivs changes, currently SIPRN B PI
  
  if(init_time %% init_stepsize != 0) {
    warning("init_time is not divisible by init_stepsize, this has not been tested")
  }
  
  #Save parameter values provided into dataframe
  # taking all different combinations
  sim_vars <- list("u_S" = u_Svals, "u_R" = u_Rvals, "k" = kvals,
                   "a" = avals, "tau" = tauvals, "b" = bvals,
                   "z" = zvals, "m" = mvals,
                   "f" = fvals, "d" = dvals,
                   "v_a1" = v_a1vals, "v_a2" = v_a2vals,
                   "init_S_dens" = init_S_dens_vals, 
                   "init_R_dens" = init_R_dens_vals,
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
    yinit <- c(S = param_combos$init_S_dens[i],
               I = 0,
               P = param_combos$init_S_dens[i]*param_combos$init_moi[i],
               R = param_combos$init_R_dens[i],
               #if N is NA, N = k-S-R
               N = ifelse(is.na(param_combos$init_N_dens[i]),
                          (param_combos$k[i]
                            - param_combos$init_S_dens[i]
                            - param_combos$init_R_dens[i]),
                          param_combos$init_N_dens[i]))
    params <- c(unlist(param_combos[i, ]),
                warnings = 0, thresh_min_dens = 10**-100)
    
    #Counters
    cntrs["i_only_pos_times"] <- 0 #num times I, but not S, > equil_cutoff_dens
    cntrs["j"] <- 0 #length counter (larger is longer times)
    cntrs["k"] <- 0 #step size counter (larger is smaller steps)
    
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
                           parms = params, hmax = hmax_val))
      })
      
      #Check for equil
      checks <- check_equil(yout_list = yout_list, cntrs = cntrs, 
                            fixed_time = fixed_time, 
                            equil_cutoff_dens = equil_cutoff_dens)
      if(checks$keep_running == FALSE) {break
      } else {cntrs <- checks$cntrs}
    }
    
    #Once end conditions triggered, if run succeeded (or warning-d)
    if(!is.null(yout_list$value)) {
      #Calculate all bacteria (B)
      yout_list$value$B <- yout_list$value$S + yout_list$value$I + yout_list$value$R
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
              equil = at_equil,
              data.table::melt(data = data.table::as.data.table(yout_list$value), 
                               id.vars = c("time"),
                               value.name = "Density", 
                               variable.name = "Pop",
                               variable.factor = FALSE),
              row.names = myrows)
    #If the run failed
    } else if (!is.null(yout_list$error)) {
      temp <- cbind(uniq_run = i, param_combos[i, ], equil = at_equil)
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
  
  return(list(ybig, y_noequil, yfail))
  
  #Code for visualizing while debugging
  # if (F) {
  #   #Code for plotting population sizes over time
  #   ymelt <- reshape2::melt(data = as.data.frame(yout_list[[2]]), 
  #                           id = c("time"),
  #                           value.name = "Density", 
  #                           variable.name = "Pop")
  #   
  #   ggplot(data = ymelt, 
  #          aes(x = time, y = Density+10, color = Pop)) +
  #     geom_line(lwd = 1.5, alpha = 1) + 
  #     scale_y_continuous(trans = "log10") +
  #     scale_x_continuous(breaks = seq(from = 0, to = max(ymelt$time), 
  #                                     by = round(max(ymelt$time)/10))) +
  #     geom_hline(yintercept = 10, lty = 2) +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #     NULL
  # }
  
  #for troubleshooting runs that don't get to equilibrium
  # if (F) {
  #   my_row <- 1
  #   myr <- y_noequil$r[my_row]
  #   mya <- y_noequil$a[my_row]
  #   myb <- y_noequil$b[my_row]
  #   mytau <- y_noequil$tau[my_row]
  #   myk <- y_noequil$K[my_row]
  #   myc <- y_noequil$c[my_row]
  # }
  
  #for troubleshooting runs that fail
  # if (F) {
  #   my_row <- 1
  #   myr <- yfail$r[my_row]
  #   mya <- yfail$a[my_row]
  #   myb <- yfail$b[my_row]
  #   mytau <- yfail$tau[my_row]
  #   myk <- yfail$K[my_row]
  #   myc <- yfail$c[my_row]
  # }
}