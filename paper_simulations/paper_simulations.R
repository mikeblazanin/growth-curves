## Import libraries ----

library(deSolve)
library(ggplot2)
library(dplyr)

#Setwd
mywd_split <- strsplit(getwd(), split = "/") 
if (mywd_split[[1]][length(mywd_split[[1]])] != "paper_simulations") {
  setwd("./paper_simulations/")
}

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

## Define derivatives function for growth curve simulations ----
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
  dY <- c(S = 0, I = 0, P = 0, R = 0)
  
  ##Calculate dS
  #dS/dt = u_S*S((k_S-S-c_SI*I-c_SR*R)/k_S) - aSP
  dY["S"] <- parms["u_S"] * y["S"] * 
    ((parms["k_S"] - y["S"] - 
        parms["c_SI"]*y["I"] - 
        parms["c_SR"]*y["R"])/parms["k_S"]) - 
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
  #dP/dt = baS(t-tau)P(t-tau) - aSP - azIP
  if (t < parms["tau"]) {
    dY["P"] <- -parms["a"] * y["S"] * y["P"] -
      parms["a"] * parms["z"] * y["I"] * y["P"]
  } else {
    dY["P"] <- parms["b"] * parms["a"] * 
      lagvalue(t-parms["tau"], 1)*lagvalue(t-parms["tau"], 3) - 
      parms["a"]*y["S"]*y["P"] -
      parms["a"] * parms["z"] * y["I"] * y["P"]
  }
  
  #Calculate dR
  #dR/dt = u_R*R((k_R-R-c_RS*S-c_RI*I)/k_R) + m*S
  dY["R"] <- parms["u_R"] * y["R"] * 
    ((parms["k_R"] - y["R"] - 
        parms["c_RS"]*y["S"] - 
        parms["c_RI"]*y["I"])/parms["k_R"]) + 
    parms["m"] * y["S"]
  
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

## Define function for running simulations across many parameter values ----
run_sims <- function(u_Svals,
                     u_Rvals = 0,
                     k_Svals,
                     k_Rvals = 1,
                     avals,
                     tauvals,
                     bvals,
                     zvals = 0,
                     mvals = 0,
                     c_SIvals = 1,
                     c_SRvals = 1,
                     c_RSvals = 1,
                     c_RIvals = 1,
                     init_S_dens_vals = 10**6,
                     init_R_dens_vals = 0,
                     init_moi_vals = 10**-2,
                     min_dens = 0.1,
                     init_time = 100,
                     init_stepsize = 1,
                     combinatorial = TRUE,
                     dynamic_stepsize = TRUE,
                     print_info = TRUE) {
  
  #Inputs: vectors of parameters to be combined factorially to make
  #         all possible combinations & run simulations with
  #       min_dens is threshold density to consider a population at equilibrium
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
  
  if(init_time %% init_stepsize != 0) {
    warning("init_time is not divisible by init_stepsize, this has not been tested")
  }
  
  #Save parameter values provided into dataframe
  # taking all different combinations
  if (combinatorial) {
    param_combos <- expand.grid(list("u_Svals" = u_Svals, "u_Rvals" = u_Rvals,
                                     "k_Svals" = k_Svals, "k_Rvals" = k_Rvals,
                                     "avals" = avals, "tauvals" = tauvals, 
                                     "bvals" = bvals, 
                                     "c_SIvals" = c_SIvals, "c_SRvals" = c_SRvals,
                                     "c_RSvals" = c_RSvals, "c_RIvals" = c_RIvals,
                                     "zvals" = zvals, "mvals" = mvals,
                                     "init_S_dens_vals" = init_S_dens_vals, 
                                     "init_R_dens_vals" = init_R_dens_vals,
                                     "init_moi_vals" = init_moi_vals),
                                stringsAsFactors = FALSE)
  } else { #not combinatorial
    num_sims <- max(sapply(X = list(u_Svals, u_Rvals, k_Svals, k_Rvals,
                                    avals, tauvals, bvals, 
                                    c_SIvals, c_SRvals, c_RSvals, c_RIvals,
                                    zvals, mvals,
                                    init_S_dens_vals, init_R_dens_vals,
                                    init_moi_vals), 
                           FUN = length))
    
    #Check for parameter lengths being non-divisible with the
    # number of simulations inferred from the longest parameter length
    if (!all(num_sims %% sapply(X = list(u_Svals, u_Rvals, k_Svals, k_Rvals,
                                         avals, tauvals, bvals, 
                                         c_SIvals, c_SRvals, c_RSvals, c_RIvals,
                                         zvals, mvals,
                                         init_S_dens_vals, init_R_dens_vals,
                                         init_moi_vals), 
                                FUN = length) == 0)) {
      warning("Combinatorial=TRUE but longest param vals length is not a multiple of all other param vals lengths")
    }
    
    #Save parameters into dataframe, replicating shorter parameter
    # vectors as needed to reach # of simulations
    param_combos <- data.frame("u_Svals" = rep_len(u_Svals, num_sims), 
                               "u_Rvals" = rep_len(u_Rvals, num_sims),
                               "k_Svals" = rep_len(k_Svals, num_sims), 
                               "k_Rvals" = rep_len(k_Rvals, num_sims),
                               "avals" = rep_len(avals, num_sims), 
                               "tauvals" = rep_len(tauvals, num_sims), 
                               "bvals" = rep_len(bvals, num_sims), 
                               "c_SIvals" = rep_len(c_SIvals, num_sims), 
                               "c_SRvals" = rep_len(c_SRvals, num_sims),
                               "c_RSvals" = rep_len(c_RSvals, num_sims),
                               "c_RIvals" = rep_len(c_RIvals, num_sims), 
                               "zvals" = rep_len(zvals, num_sims), 
                               "mvals" = rep_len(mvals, num_sims),
                               "init_S_dens_vals" = rep_len(init_S_dens_vals, num_sims), 
                               "init_R_dens_vals" = rep_len(init_R_dens_vals, num_sims),
                               "init_moi_vals" = rep_len(init_moi_vals, num_sims),
                               stringsAsFactors = FALSE)
  }
  
  num_sims <- nrow(param_combos)
  
  #Print number of simulations that will be run
  if(print_info) {
    print(paste(num_sims, "simulations will be run"))
    
    #Save sequence of 10% cutoff points for later reference
    progress_seq <- round(seq(from = 0, to = num_sims, by = num_sims/10))
  }
  
  #Make placeholders
  yfail <- NULL #for runs that fail
  #for runs that succeed, pre-allocate ybig now to save on memory/speed
  ybig <- data.frame("uniq_run" = rep(NA, 6*(1+init_time/init_stepsize)*num_sims),
                     "u_S" = NA, "u_R" = NA, 
                     "k_S" = NA, "k_R" = NA,
                     "a" = NA, "b" = NA, "tau" = NA,
                     "c_SI" = NA, "c_SR" = NA,
                     "c_RS" = NA, "c_RI" = NA, 
                     "z" = NA, "m" = NA,
                     "init_S_dens" = NA, "init_R_dens" = NA,
                     "init_moi" = NA,
                     "equil" = NA, "time" = NA, "Pop" = as.character(NA),
                     "Density" = NA, stringsAsFactors = FALSE)
  
  #Define counters
  rows_tracking <- list(
    #row where data should start being entered for next simulation
    "start_row" = 1,
    #number of rows this simulation is 
    # (default value provided for dynamic_stepsize = FALSE
    #  but when dynamic_stepsize = TRUE this will be overwritten ea time)
    "this_run_nrows" = 6*(1+init_time/init_stepsize),
    #counter for additional rows to add, to minimize number of times
    # ybig has to be re-defined
    "still_needed_toadd" = 0)
  
  for (i in 1:nrow(param_combos)) { #i acts as the uniq_run counter
    #Define pops & parameters
    yinit <- c(S = param_combos$init_S_dens_vals[i],
               I = 0,
               P = param_combos$init_S_dens_vals[i]*param_combos$init_moi_vals[i],
               R = param_combos$init_R_dens_vals[i])
    params <- c(u_S = param_combos$u_Svals[i],
                u_R = param_combos$u_Rvals[i],
                k_S = param_combos$k_Svals[i],
                k_R = param_combos$k_Rvals[i],
                a = param_combos$avals[i],
                tau = param_combos$tauvals[i],
                b = param_combos$bvals[i],
                c_SI = param_combos$c_SIvals[i],
                c_SR = param_combos$c_SRvals[i],
                c_RS = param_combos$c_RSvals[i],
                c_RI = param_combos$c_RIvals[i],
                z = param_combos$zvals[i],
                m = param_combos$mvals[i],
                warnings = 0, thresh_min_dens = 10**-100)
    
    #Run simulation(s) with longer & longer times until equil reached
    #Also, if equil has non-zero I run with shorter steps
    keep_running <- TRUE #placeholder for triggering end of sims
    #Placeholder for the number of times I has been detected above
    # min_dens while S has not been
    i_only_pos_times <- 0
    j <- 0 #length counter (larger is longer times)
    k <- 0 #step size counter (larger is smaller steps)
    while(keep_running) {
      #Define times
      if (dynamic_stepsize) {
        #If dynamic_stepsize true, double lengths & steps for ea j count
        # (so that the number of timepoints returned is constant)
        times <- seq(0, init_time*2**j, init_stepsize*2**j)
      } else {
        #If dynamic_stepsize false, keep stepsize at init_stepsize
        times <- seq(0, init_time*2**j, init_stepsize)
      }
      
      #Calculate hmax (max step size integrator uses)
      if (dynamic_stepsize) {
        #Note that the max step size for the integrator is the 
        # same as our step size except halved for each k count
        hmax_val <- init_stepsize*2**(j-k)
      } else {hmax_val <- min(init_stepsize*2**(j-k), init_stepsize)}
      
      #Run simulation
      yout_list <- myTryCatch(expr = {
        as.data.frame(dede(y = yinit, times = times, func = derivs, 
                           parms = params, hmax = hmax_val))
      })
      
      #Infinite loop prevention check (j = 10 is 24 hrs)
      if (j >= 10 | k >= 15 | j+k >= 20) {
        keep_running <- FALSE
        at_equil <- FALSE
      }
      
      #If there was an error, increase k by 1 and re-run
      if(!is.null(yout_list$error)) {
        k <- k+1
        #If there was a warning, could be several causes, so we
        # generally just halve step size and increase length
      } else if (!is.null(yout_list$warning)) {
        j <- j+1
        k <- k+2
        #If it was successful, check for equilibrium
      } else if (is.null(yout_list$warning) & is.null(yout_list$error)) {
        #First drop all rows with nan
        yout_list$value <- yout_list$value[!(is.nan(yout_list$value$S) |
                                               is.nan(yout_list$value$I) |
                                               is.nan(yout_list$value$P) |
                                               is.nan(yout_list$value$R)), ]
        
        #S and I both at equil, we're done
        if (yout_list$value$S[nrow(yout_list$value)] < min_dens & 
            yout_list$value$I[nrow(yout_list$value)] < min_dens) {
          keep_running <- FALSE
          at_equil <- TRUE
          #S not at equil, need more time
        } else if (yout_list$value$S[nrow(yout_list$value)] >= min_dens) { 
          j <- j+1
          #I not at equil (but S is because above check failed),
          #   first we'll lengthen the simulation
          #    (to make sure it was long enough to catch the last burst)
          #   then we'll start shrinking our step size
        } else if (yout_list$value$I[nrow(yout_list$value)] >= min_dens) {
          if (i_only_pos_times < 1) {
            j <- j+1
            i_only_pos_times <- i_only_pos_times+1
          } else {
            k <- k+1
          }
        }
      } else {stop("tryCatch failed, niether success, warning, nor error detected")}
    }
    
    #Once end conditions triggered, if run succeeded (or warning-d)
    if(!is.null(yout_list$value)) {
      #Calculate all bacteria (B)
      yout_list$value$B <- yout_list$value$S + yout_list$value$I + yout_list$value$R
      #Calculate all phage (PI)
      yout_list$value$PI <- yout_list$value$P + yout_list$value$I
      
      if (!dynamic_stepsize) {
        rows_tracking$this_run_nrows <- 6*nrow(yout_list$value)
        rows_tracking$still_needed_toadd <- 
          rows_tracking$still_needed_toadd + 
          (rows_tracking$this_run_nrows - 6*(1+init_time/init_stepsize))
        
        #If the expected end row of this simulation
        #is larger than the number of rows available, 
        #we need to add more rows
        if((rows_tracking$start_row + rows_tracking$this_run_nrows - 1) >
           nrow(ybig)) {
          ybig <- 
            rbind(ybig,
                  data.frame("uniq_run" = rep(NA, rows_tracking$still_needed_toadd),
                             "u_S" = NA, "u_R" = NA, 
                             "k_S" = NA, "k_R" = NA,
                             "a" = NA, "b" = NA, "tau" = NA,
                             "c_SI" = NA, "c_SR" = NA,
                             "c_RS" = NA, "c_RI" = NA, 
                             "z" = NA, "m" = NA,
                             "init_S_dens" = NA, "init_R_dens" = NA,
                             "init_moi" = NA,
                             "equil" = NA, "time" = NA, "Pop" = as.character(NA),
                             "Density" = NA, stringsAsFactors = FALSE))
          
          rows_tracking$still_needed_toadd <- 0
        }
      }
      
      #Reshape, add parameters, and fill into ybig in right rows
      ybig[rows_tracking$start_row : 
             (rows_tracking$start_row + rows_tracking$this_run_nrows - 1), ] <-
        cbind(data.frame(uniq_run = i, 
                         u_S = param_combos$u_Svals[i], 
                         u_R = param_combos$u_Rvals[i], 
                         k_S = param_combos$k_Svals[i], 
                         k_R = param_combos$k_Rvals[i],
                         a = param_combos$avals[i], 
                         b = param_combos$bvals[i], 
                         tau = param_combos$tauvals[i],
                         c_SI = param_combos$c_SIvals[i],
                         c_SR = param_combos$c_SRvals[i],
                         c_RS = param_combos$c_RSvals[i],
                         c_RI = param_combos$c_RIvals[i],
                         z = param_combos$zvals[i],
                         m = param_combos$mvals[i],
                         init_S_dens = param_combos$init_S_dens_vals[i], 
                         init_R_dens = param_combos$init_R_dens_vals[i], 
                         init_moi = param_combos$init_moi_vals[i],
                         equil = at_equil),
              data.table::melt(data = data.table::as.data.table(yout_list$value), 
                               id.vars = c("time"),
                               value.name = "Density", 
                               variable.name = "Pop",
                               variable.factor = FALSE))
      #If the run failed
    } else if (!is.null(yout_list$error)) {
      temp <- data.frame(uniq_run = i, 
                         u_S = param_combos$u_Svals[i], 
                         u_R = param_combos$u_Rvals[i], 
                         k_S = param_combos$k_Svals[i], 
                         k_R = param_combos$k_Rvals[i],
                         a = param_combos$avals[i], 
                         b = param_combos$bvals[i], 
                         tau = param_combos$tauvals[i],
                         c_SI = param_combos$c_SIvals[i],
                         c_SR = param_combos$c_SRvals[i],
                         c_RS = param_combos$c_RSvals[i],
                         c_RI = param_combos$c_RIvals[i],
                         z = param_combos$zvals[i],
                         m = param_combos$mvals[i],
                         init_S_dens = param_combos$init_S_dens_vals[i], 
                         init_R_dens = param_combos$init_R_dens_vals[i], 
                         init_moi = param_combos$init_moi_vals[i],
                         equil = at_equil)
      if (is.null(yfail)) { #This is the first failed run
        yfail <- temp
      } else { #This is a non-first failed run
        yfail <- rbind(yfail, temp)
      }
    } else {stop("tryCatch failed during saving, neither success nor warning nor error detected")}
    
    #Update cumulative offset (for non-dynamic stepsize runs)
    rows_tracking$start_row <- rows_tracking$start_row + rows_tracking$this_run_nrows
    
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
      y_noequil <- ybig[min(which(ybig$uniq_run == run)), 1:18]
    } else {
      y_noequil <- rbind(y_noequil, ybig[min(which(ybig$uniq_run == run)), 1:18])
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

##Define file save/load wrapper for run_sims ----
run_sims_filewrapper <- function(name, dir = ".",
                                 read_file = TRUE, write_file = TRUE,
                                 ...) {
  
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
    temp <- run_sims(...)
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

##Define peak-finding function ----
find_local_extrema <- function(values, 
                               return_maxima = TRUE,
                               return_minima = TRUE,
                               width_limit = NULL,
                               height_limit = NULL,
                               remove_endpoints = TRUE,
                               na.rm = FALSE) {
  #Takes a vector of values and returns a vector of the indices
  # of all local value extrema (by default, returns both maxima and minima)
  # To only return maxima or minima, change return_maxima/return_minima to FALSE
  
  #width_limit and/or height_limit must be provided
  #Width is how wide the window will be to look for a maxima/minima
  # Narrower width will be more sensitive to narrow local maxima/minima
  # Wider width will be less sensitive to narrow local maxima/minima
  #Height is how high or low a single step is allowed to take
  # e.g. a maxima-finding function will not pass a valley deeper
  # than height_limit
  #Note that this also limits approaches to extrema, so if set too small
  # function may converge on non-peaks
  #If both width_limit and height_limit are provided, steps are limited
  # conservatively (a single step must meet both criteria)
  
  #This function is designed to be compatible with dplyr::group_by and summarize
  
  #Check inputs
  if (!return_maxima & !return_minima) {
    stop("Both return_maxima and return_minima are FALSE, at least one must be TRUE")
  }
  if (is.null(width_limit) & is.null(height_limit)) {
    stop("Either width_limit or height_limit must be provided")
  }
  if (!is.null(width_limit)) {
    if (width_limit%%2 == 0) {
      warning("width_limit must be odd, will use ", width_limit-1, " as width_limit")
      width_limit <- width_limit - 1
    }
  }
  if (is.null(width_limit) & !is.null(height_limit)) {
    warning("height_limit alone tends to be sensitive to height_limit parameter, use with caution")
  }
  if (na.rm == TRUE & sum(is.na(values)) > 0) {
    if (!all(is.na(values[(1+length(values)-sum(is.na(values))):length(values)]))) {
      warning("Removing NAs found within values vector, returned indices will refer to non-NA values")
      print(values)
    }
    values <- values[!is.na(values)]
  } else if(any(is.na(values))) {
    stop("Some provided values are NA and na.rm = FALSE")
  }
  
  #Define sub-function to find limits of the window
  get_window_limits <- function(cnt_pos,
                                width_limit = NULL,
                                height_limit = NULL,
                                looking_for = c("minima", "maxima"),
                                values = NULL) {
    #Check inputs
    if (length(looking_for) > 1) {stop("looking_for must be specified")}
    if (!is.null(height_limit) & is.null(values)) {
      stop("height_limit is specified, but no values are provided")
    }
    if (is.null(width_limit) & is.null(height_limit)) {
      stop("Either width_limit or height_limit must be provided")
    }
    
    #Define window limits
    window_start <- c(NA, NA)
    if (!is.null(width_limit)) { #using width limit
      window_start[1] <- max(c(1, cnt_pos-floor(width_limit/2)))
    }
    if (!is.null(height_limit)) { #using height limig
      #For startpoint height, we want the latest point that is
      #behind of our current point and
      #either:
      # below current height - height limit
      # or above current height + height limit
      #Then we move one place forward 
      # (so it's the last value w/in height limit)
      window_start[2] <- max(c(1,
                               1+which(1:length(values) < cnt_pos &
                                         (values >= (values[cnt_pos] + height_limit) |
                                            values <= (values[cnt_pos] - height_limit)))))
      #Make sure we're going at least 1 point backwards
      if(window_start[2] >= cnt_pos) {window_start[2] <- cnt_pos-1}
    }
    window_end <- c(NA, NA)
    if (!is.null(width_limit)) { #using width limit
      window_end[1] <- min(c(length(values), cnt_pos+floor(width_limit/2)))
    }
    if (!is.null(height_limit)) { #using height limit
      #For endpoint height, we want the earliest point that is
      #forward of our current point and
      #either:
      # below current height - height limit
      # or above current height + height limit
      #Then we move one place back 
      # (so it's the last value w/in height limit)
      window_end[2] <- min(c(length(values),
                             -1+which(1:length(values) > cnt_pos & #not backwards
                                        (values <= (values[cnt_pos] - height_limit) |
                                           values >= (values[cnt_pos] + height_limit)))))
      #Make sure we're going at least one point forwards
      if (window_end[2] <= cnt_pos) {window_end[2] <- cnt_pos+1}
    }
    return(c(max(window_start, na.rm = T), min(window_end, na.rm = T)))
  }
  
  find_next_extrema <- function(cnt_pos, values,
                                width_limit = NULL,
                                height_limit = NULL,
                                looking_for = c("minima", "maxima")) {
    if (cnt_pos == length(values)) {best_pos <- cnt_pos-1
    } else {best_pos <- cnt_pos+1}
    
    #Save the starting position so we never go backwards
    start_pos <- cnt_pos
    
    ##Looking for next maxima
    if(looking_for == "maxima") {
      while (cnt_pos != best_pos) {
        #Move the previous best pointer to current pointer location
        best_pos <- cnt_pos
        #Get next window limits
        window_lims <- get_window_limits(cnt_pos = cnt_pos,
                                         width_limit = width_limit,
                                         height_limit = height_limit,
                                         looking_for = "maxima",
                                         values = values)
        #Make sure we're not going backwards
        window_lims <- c(max(start_pos, window_lims[1]),
                         max(start_pos, window_lims[2]))
        #Then move current pointer to highest point within window
        # (making sure not to check non-integer indices, or indices below 1 or
        #  higher than the length of the vector)
        cnt_pos <- window_lims[1]-1+which.max(values[window_lims[1]:window_lims[2]])
      }
      ##Looking for next minima
    } else if (looking_for == "minima") {
      while (cnt_pos != best_pos) {
        #Move the previous best pointer to current pointer location
        best_pos <- cnt_pos
        #Get next window limits
        window_lims <- get_window_limits(cnt_pos = cnt_pos,
                                         width_limit = width_limit,
                                         height_limit = height_limit,
                                         looking_for = "minima",
                                         values = values)
        #Make sure we're not going backwards
        window_lims <- c(max(start_pos, window_lims[1]),
                         max(start_pos, window_lims[2]))
        #Then move current pointer to lowest point within window
        # (making sure not to check non-integer indices, or indices below 1 or
        #  higher than the length of the vector)
        cnt_pos <- window_lims[1]-1+which.min(values[window_lims[1]:window_lims[2]])
      }
    }
    return(best_pos)
  }
  
  cnt_pos <- 1
  ##Find first maxima
  maxima_list <- c(find_next_extrema(cnt_pos, values,
                                     width_limit = width_limit,
                                     height_limit = height_limit,
                                     looking_for = "maxima"))
  ##Find first minima
  minima_list <- c(find_next_extrema(cnt_pos, values,
                                     width_limit = width_limit,
                                     height_limit = height_limit,
                                     looking_for = "minima"))
  
  ##Check for next extrema until...
  while (TRUE) {
    #we're finding repeats
    if (any(duplicated(c(minima_list, maxima_list)))) {break}
    #or we hit the end of the values
    if (length(values) %in% c(maxima_list, minima_list)) {
      break
    }
    #Since maxima & minima must alternate, always start with furthest one 
    # we've found so far
    cnt_pos <- max(c(minima_list, maxima_list))
    #we're looking for a maxima next
    if (cnt_pos %in% minima_list) {
      maxima_list <- c(maxima_list,
                       find_next_extrema(cnt_pos, values,
                                         width_limit = width_limit,
                                         height_limit = height_limit,
                                         looking_for = "maxima"))
      #we're looking for a minima next
    } else if (cnt_pos %in% maxima_list) {
      minima_list <- c(minima_list,
                       find_next_extrema(cnt_pos, values,
                                         width_limit = width_limit,
                                         height_limit = height_limit,
                                         looking_for = "minima"))
    }
  }
  
  #Combine maxima & minima values & remove duplicates
  output <- c()
  if (return_maxima) {output <- c(output, maxima_list)}
  if (return_minima) {output <- c(output, minima_list)}
  #If remove endpoints is true, remove first or last values from return
  if (remove_endpoints) {
    if (1 %in% output) {output <- output[-which(output == 1)]}
    if (length(values) %in% output) {
      output <- output[-which(output == length(values))]}
  }
  #Remove duplicates
  output <- unique(output)
  #Order
  output <- output[order(output)]
  
  return(output)
}

## Global Settings ----
glob_read_files <- TRUE
glob_make_curveplots <- FALSE
glob_make_statplots <- FALSE

