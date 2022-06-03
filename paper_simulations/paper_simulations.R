## Import libraries ----
library(deSolve)
library(ggplot2)
library(dplyr)
library(plotly)

#Setwd
mywd_split <- strsplit(getwd(), split = "/") 
if (mywd_split[[1]][length(mywd_split[[1]])] == "growth-curves") {
  dir.create("paper_simulations", showWarnings = FALSE)
  setwd("./paper_simulations/")
} else {
  stop("Not in correct root directory")
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
      
      #Infinite loop prevention check 
      # (6400 mins = 4.4 days)
      # (0.001 mins = 0.06 secs)
      if (init_time*2**j > 6400 | hmax_val <= 0.001) {
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

#Parameter values to be used ----
# u_S: 2^(-5, -5.5, -6, -6.5, -7)/min (doub time 22, 31, 44, 63, 88 mins)
# K: 10^(6, 7, 8, 9, 10) cfu/mL
# a: 10^(-8, -9, -10, -11, -12) /min
#       (range of wt adsorption is -8 to -12,
#       and slower infec leads to gc lasting > 3 days or bact near k
#       (discovered from sims including -14 and -16 in a range))
# tau: 10^(1.2, 1.4, 1.6, 1.8, 2)
# b: 10^(1, 1.5, 2, 2.5, 3)
# c: 1
# init_dens: 10^(4, 5, 6, 7, 8)
# init_moi: 10^(-4, -3, -2, -1, 0)
# 
#Units: cfu/mL, pfu/mL, minutes

##Run 1: phage variants simulation & summarization ----
run1 <- run_sims_filewrapper(name = "run1",
                     u_Svals = signif(2**c(-5.5, -6.5), 2),
                     k_Svals = 10**c(7, 9),
                     avals = 10**c(-8, -9, -10, -11, -12),
                     tauvals = signif(10**c(1.2, 1.4, 1.6, 1.8, 2), 2),
                     bvals = signif(10**c(1, 1.5, 2, 2.5, 3), 2),
                     c_SIvals = 1,
                     init_S_dens_vals = 10**6,
                     init_moi_vals = 10**-2,
                     min_dens = 100,
                     init_time = 100,
                     init_stepsize = 1,
                     print_info = TRUE,
                     read_file = glob_read_files)

ybig1 <- run1[[1]]

ybig1 <- group_by_at(ybig1, .vars = 1:17)

#Summarize runs (using ifelse to handle potential non-equil conditions)
y_summarized1 <- 
  dplyr::summarize(
    ybig1,
    max_dens = max(Density[Pop == "B"]),
    max_time = time[Pop == "B" & Density[Pop == "B"] == max_dens],
    #Make references for finding extin point
    extin_dens = 10**4,
    #Technically this is the first index after extinction
    extin_index =
      ifelse(any(Pop == "B" & Density <= extin_dens & time >= max_time),
             min(which(Pop == "B" & Density <= extin_dens & time >= max_time)),
             NA),
    extin_index_back1 = 
      ifelse(is.na(extin_index), NA,
             which(Pop == "B")[match(extin_index, which(Pop == "B")) - 1]),
    #use linear interpolation to find extin time
    extin_time = 
      ifelse(is.na(extin_index), NA,
             time[extin_index] - (Density[extin_index] - extin_dens)*
               (time[extin_index] - time[extin_index_back1])/
               (Density[extin_index] - Density[extin_index_back1])),
    #make references for auc (here indices are within Pop == "B")
    extin_index_winB = 
      ifelse(is.na(extin_time), NA, which(which(Pop == "B") == extin_index)),
    extin_index_back1_winB = 
      ifelse(is.na(extin_time), NA, which(which(Pop == "B") == extin_index)-1),
    extin_index_back2_winB = 
      ifelse(is.na(extin_time), NA, which(which(Pop == "B") == extin_index)-2),
    #using trapezoid rule to find auc
    auc =
      ifelse(is.na(extin_time), NA,
             #trapezoids of all intervals before one ending at extin time
             (time[Pop == "B"][2] - time[Pop == "B"][1])/2 *
               (Density[Pop == "B"][1] +
                  2*sum(Density[Pop == "B"][2:extin_index_back2_winB]) +
                  Density[Pop == "B"][extin_index_back1_winB]) +
               #trapezoid that ends at extin time
               (extin_time-time[extin_index_back1])/2 *
               (Density[extin_index_back1] + extin_dens)),
    phage_final = ifelse(equil[1], max(Density[Pop == "P"]), NA),
    #make references for phage extin dens
    phage_dens_y1 = ifelse(is.na(extin_time), NA,
                           Density[max(which(Pop == "P" & time < extin_time))]),
    phage_dens_y2 = ifelse(is.na(extin_time), NA,
                           Density[which(Pop == "P" & time == time[extin_index])]),
    phage_time_x1 = ifelse(is.na(extin_time), NA,
                           time[max(which(Pop == "P" & time < extin_time))]),
    phage_time_x2 = ifelse(is.na(extin_time), NA,
                           time[which(Pop == "P" & time == time[extin_index])]),
    #using linear interpolation to find phage dens at bact extinction
    phage_extin = ifelse(is.na(extin_time), NA,
                         phage_dens_y1 + 
                           (phage_dens_y2-phage_dens_y1)/(phage_time_x2-phage_time_x1)*
                           (extin_time-phage_time_x1)),
    # phage_r =
    #   (log(phage_final)- log(init_S_dens[1]*init_moi[1]))/extin_time,
    # phage_r_extin =
    #   (log(phage_extin) - log(init_S_dens[1]*init_moi[1]))/extin_time,
    phage_atmaxdens = Density[Pop == "P" & time == max_time],
    maxdens_k_ratio = max_dens/k_S[1],
    run_time = max(time),
    equil = equil[1],
    auc_tomaxtime = gcplyr::auc(x = time[Pop == "B"], 
                                y = Density[Pop == "B"], 
                                xlim = c(NA, max_time)),
    #Integral of the logistic function
    # P(t) = K/(1+((K-P_0)/P_0)e^(-rt))
    # indefinite integral = K((ln(((K - P_0)/P_0)*e^(-rt) + 1))/r + 1) + C
    # here taking integral from 0 to max_time
    pred_auc =
      (k_S[1]*(log((k_S[1]-init_S_dens[1])/init_S_dens[1]*exp(-u_S[1]*max_time) + 1)/u_S[1] + max_time) -
      k_S[1]*(log((k_S[1]-init_S_dens[1])/init_S_dens[1]*exp(-u_S[1]*0) + 1)/u_S[1] + 0))
  )

#Drop unneeded columns
y_summarized1 <- subset(y_summarized1,
                        select = -c(extin_index, extin_index_back1,
                                  extin_index_winB, extin_index_back1_winB,
                                  extin_index_back2_winB, phage_dens_y1,
                                  phage_dens_y2, phage_time_x1,
                                  phage_time_x2))
y_summarized1 <- as.data.frame(y_summarized1)

#Run 1: density dynamics ----
dir.create("./run1_noequil/", showWarnings = F)
if (glob_make_curveplots) {
  for (myrun in unique(ybig1$uniq_run[ybig1$equil == FALSE])) {
    tiff(paste("./run1_noequil/", myrun, ".tiff", sep = ""),
         width = 4, height = 4, units = "in", res = 200)
    print(ggplot(data = ybig1[ybig1$uniq_run == myrun &
                                ybig1$Pop %in% c("S", "I", "P"), ],
                 aes(x = time, y = Density+1, color = Pop)) +
            geom_line(lwd = 1.5, alpha = 0.5) +
            scale_y_continuous(trans = "log10"))
    dev.off()
  }
}

dir.create("./run1_equil/", showWarnings = F)
if (glob_make_curveplots) {
  for (myrun in unique(ybig1$uniq_run[ybig1$equil == TRUE])) {
    tiff(paste("./run1_equil/", myrun, ".tiff", sep = ""),
         width = 4, height = 4, units = "in", res = 200)
    print(ggplot(data = ybig1[ybig1$uniq_run == myrun &
                                ybig1$Pop %in% c("S", "I", "P"), ],
                 aes(x = time, y = Density+1, color = Pop)) +
            geom_line(lwd = 1.5, alpha = 0.5) +
            scale_y_continuous(trans = "log10"))
    dev.off()
  }
}

#Run 1: visualizing which runs didn't equil----
if (glob_make_statplots) {
  for (myu in unique(y_summarized1$u_S)) {
    for (myk in unique(y_summarized1$k_S)) {
      print(ggplot(y_summarized1[y_summarized1$u_S == myu &
                                   y_summarized1$k_S == myk, ],
                   aes(x = a, y = b, 
                       color = as.character(!is.na(extin_time)))) +
              geom_point() +
              facet_wrap(~tau) +
              scale_x_continuous(trans = "log10") +
              scale_y_continuous(trans = "log10") +
              scale_color_manual(name = "Reached\nExtinction",
                                 limits = c("TRUE", "FALSE"),
                                 values = c("black", "gray"),
                                 drop = FALSE) +
              ggtitle(paste("u=", myu, " k=", myk, sep = "")) +
              NULL)
    }
  }
}

###TODO: fix colors here
#Code for Plotly 3D plot of non-equil runs
if (glob_make_statplots) {
  for (myu in unique(y_summarized1$u_S)) {
    for (myk in unique(y_summarized1$k_S)) {
      temp <- y_summarized1[y_summarized1$u_S == myu & y_summarized1$k_S == myk, ]
      plt <- plot_ly(x = log10(temp$a), y = log10(temp$b), 
                     z = log10(temp$tau), type = "scatter3d", 
                     mode = "markers", 
                     color = as.character(!is.na(temp$extin_time)),
                     marker = list(size = 5)) %>%
        layout(title = "Extinction within 4 days?",
               scene = list(xaxis = list(title = "a", tickvals = c(-11, -9),
                                         ticktext = c("10^-11", "10^-9")), 
                            yaxis = list(title = "b", 
                                         tickvals = c(log10(32), log10(320)),
                                         ticktext = c("32", "320")),
                            zaxis = list(title = "tau",
                                         tickvals = c(log10(25), log10(63)),
                                         ticktext = c("25", "63")),
                            camera = list(eye = list(x = -1.77, y = 1.21, z = 1.5))))
      print(plt)
      readline(prompt = paste("Plotting u_S=", myu,
                              ", k_S=", myk, 
                              ". Save, then press any key to continue",
                              sep = ""))
    }
  }
}

#Run 1: max dens-max time plots ----
dir.create("./plots/", showWarnings = F)
if (glob_make_statplots) {
  print(ggplot(data = y_summarized1,
         aes(x = max_time, y = max_dens)) +
    facet_grid(k_S ~ u_S, scales = "free") +
    #scale_y_continuous(trans = "log10") +
    geom_point())
  
  tiff("./plots/run1_Bcurves_tau.tiff", width = 5, height = 4, 
       res = 300, units = "in")
  print(ggplot(data = ybig1[ybig1$a == 10**-10 & ybig1$b == 100 &
                        ybig1$Pop == "B", ],
         aes(x = time/60, y = Density, color = as.factor(tau), group = uniq_run)) +
    geom_line(lwd = 1.5, alpha = 0.7) +
    facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    coord_cartesian(ylim = c(10**5, NA), xlim = c(0, 24)) +
    scale_color_viridis(discrete = TRUE) +
    NULL)
  dev.off()
  
  tiff("./plots/run1_Bcurves_a.tiff", width = 5, height = 4, 
       res = 300, units = "in")
  print(ggplot(data = ybig1[ybig1$b == 100 & ybig1$tau == 40 &
                        ybig1$Pop == "B" & ybig1$Density >= 10**2, ],
         aes(x = time/60, y = Density, color = as.factor(a), group = uniq_run)) +
    geom_line(lwd = 1.5, alpha = 0.7) +
    facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    coord_cartesian(ylim = c(10**5, NA), xlim = c(0, 24)) +
    scale_color_viridis(discrete = TRUE) +
    NULL)
  dev.off()
  
  tiff("./plots/run1_Bcurves_b.tiff", width = 5, height = 4, 
       res = 300, units = "in")
  print(ggplot(data = ybig1[ybig1$a == 10**-10 & ybig1$tau == 40 &
                        ybig1$Pop == "B" & ybig1$Density >= 10**2, ],
         aes(x = time/60, y = Density, color = as.factor(b), group = uniq_run)) +
    geom_line(lwd = 1.5, alpha = 0.7) +
    facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    coord_cartesian(ylim = c(10**5, NA), xlim = c(0, 24)) +
    scale_color_viridis(discrete = TRUE) +
    NULL)
  dev.off()
  
  tiff("./plots/run1_maxdens_maxtime_a.tiff", width = 5, height = 4, 
       res = 300, units = "in")
  print(ggplot(data = y_summarized1,
         aes(x = max_time/60, y = max_dens, color = as.factor(a))) +
    geom_point() +
    facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(breaks = c(0, 12, 24)) +
    #coord_cartesian(xlim = c(0, 36)) +
    scale_color_viridis_d() +
    NULL)
  dev.off()
  
  tiff("./plots/run1_maxdens_maxtime_a_onefacet.tiff", width = 5, height = 4, 
       res = 300, units = "in")
  print(ggplot(data = y_summarized1[y_summarized1$u_S == 0.011 &
                                      y_summarized1$k_S == 10**9, ],
               aes(x = max_time/60, y = max_dens, color = as.factor(a),
                   shape = as.factor(maxdens_k_ratio < 0.95))) +
          geom_point() +
          #facet_grid(k_S ~ u_S, scales = "free") +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(breaks = c(0, 12, 24)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          #coord_cartesian(xlim = c(0, 36)) +
          scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
          labs(x = "Peak Bacterial Density Time (hrs)",
               y = "Peak Bacterial Density (cfu/mL)") +
          theme_bw() +
          guides(shape = FALSE) +
          NULL)
  dev.off()
  
  tiff("./plots/run1_maxdens_maxtime_b.tiff", width = 5, height = 4, 
       res = 300, units = "in")
  print(ggplot(data = y_summarized1,
         aes(x = max_time/60, y = max_dens, color = as.factor(b))) +
    geom_point() +
    facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(breaks = c(0, 12, 24)) +
    #coord_cartesian(xlim = c(0, 36)) +
    scale_color_viridis(discrete = TRUE) +
    NULL)
  dev.off()
  
  tiff("./plots/run1_maxdens_maxtime_tau.tiff", width = 5, height = 4, 
       res = 300, units = "in")
  print(ggplot(data = y_summarized1,
         aes(x = max_time/60, y = max_dens, color = as.factor(tau))) +
    geom_point() +
    facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(breaks = c(0, 12, 24)) +
    #coord_cartesian(xlim = c(0, 36)) +
    scale_color_viridis(discrete = TRUE) +
    NULL)
  dev.off()
}

#Run 1: final phage vs peak bacteria ----
y_summarized1$pred_phage_final <- y_summarized1$b * y_summarized1$max_dens

if (glob_make_statplots) {
  tiff("./plots/run1_maxdens_phagefinal.tiff", width = 5, height = 4, 
       res = 300, units = "in")
  print(ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                y_summarized1$u_S == 0.011, ],
         aes(x = max_dens, y = phage_final, color = as.factor(b))) +
           geom_point() +
    geom_line(aes(y = pred_phage_final)) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    #facet_grid(u_S ~ k_S) +
    scale_color_viridis(discrete = TRUE, end = 0.95) +
    theme_bw() +
    NULL)
  dev.off()
}

##Run 1: landscape of peak time ----
#Code for Plotly 3D plot
if (glob_make_statplots) {
  for (myu in unique(y_summarized1$u_S)) {
    for (myk in unique(y_summarized1$k_S)) {
      temp <- y_summarized1[y_summarized1$u_S == myu & y_summarized1$k_S == myk, ]
      temp <- temp[order(temp$a, temp$b), ]
      temp$maxtime_scld <- c(scale(temp$max_time,
                                center = min(temp$max_time),
                                scale = max(temp$max_time)))
      #Make plot with surfaces by tau
      plt <- plot_ly(x = unique(log10(temp$a)), y = unique(log10(temp$b)))
      for (i in 1:length(unique(y_summarized1$tau))) {
        mytau <- unique(y_summarized1$tau)[i]
        plt <- plt %>% 
          add_surface(z = 0.5*matrix(temp$maxtime_scld[temp$tau == mytau],
                                 nrow = 5, ncol = 5)+(i-1),
                      opacity = 1.1-0.1*i)
      }
      print(plt)
    }
  }
}

#Code for 2D contours
if(glob_make_statplots) {
  tiff("./plots/run1_peaktime_contour_b100_u011_k1e9.tiff", width = 5, height = 4,
       units = "in", res = 300)
  print(ggplot(data = y_summarized1[y_summarized1$b == 100 &
                                      y_summarized1$u_S == 0.011 &
                                      y_summarized1$k_S == 10**9, ],
               aes(x = log10(a), y = tau)) +
          geom_contour_filled(aes(z = max_time/60), alpha = 0.5) +
          geom_point(aes(color = max_time/60, shape = maxdens_k_ratio < 0.95),
                     size = 3) +
          scale_color_viridis(name = "Peak time (hr)",
                              breaks = c(6, 12, 18, 24)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          scale_y_continuous(trans = "log10", breaks = c(16, 40, 100)) +
          xlab("Log10(infection rate) (/min)") +
          ylab("Lysis time (min)") +
          guides(fill = FALSE, shape = FALSE) +
          NULL)
  dev.off()
  
  #b in facets
  tiff("./plots/run1_peaktime_contour_facetb.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(ggplot(data = y_summarized1,
               aes(x = log10(a), y = tau)) +
          geom_contour_filled(aes(z = max_time), alpha = 0.7) +
          geom_point(aes(color = max_time/60, shape = maxdens_k_ratio < 0.95),
                     size = 1.5) +
          scale_color_viridis(name = "Peak time (hr)",
                              breaks = c(6, 12, 18, 24)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          facet_grid(u_S*k_S~b) +
          scale_y_continuous(trans = "log10", breaks = c(16, 40, 100)) +
          xlab("Log10(infection rate) (/min)") +
          ylab("Lysis time (min)") +
          guides(fill = FALSE, shape = FALSE) +
          labs(subtitle = "Burst size") +
          NULL)
  dev.off()
  
  #a in facets
  tiff("./plots/run1_peaktime_contour_faceta.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(ggplot(data = y_summarized1,
               aes(x = b, y = tau)) +
          geom_contour_filled(aes(z = max_time), alpha = 0.7) +
          geom_point(aes(color = max_time/60, shape = maxdens_k_ratio < 0.95),
                     size = 1.5) +
          scale_color_viridis(name = "Peak time (hr)",
                              breaks = c(6, 12, 18, 24)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          facet_grid(u_S*k_S~log10(a)) +
          scale_x_continuous(trans = "log10") +
          scale_y_continuous(trans = "log10", breaks = c(16, 40, 100)) +
          xlab("Burst size") +
          ylab("Lysis time (min)") +
          guides(fill = FALSE, shape = FALSE) +
          labs(subtitle = "Log10(infection rate) (/min)") +
          NULL)
  dev.off()
  
  #tau in facets
  tiff("./plots/run1_peaktime_contour_facettau.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(ggplot(data = y_summarized1,
               aes(x = log10(a), y = b)) +
          geom_contour_filled(aes(z = max_time), alpha = 0.7) +
          geom_point(aes(color = max_time/60, shape = maxdens_k_ratio < 0.95),
                     size = 1.5) +
          scale_color_viridis(name = "Peak time (hr)",
                              breaks = c(6, 12, 18, 24)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          facet_grid(u_S*k_S~tau) +
          scale_y_continuous(trans = "log10") +
          xlab("Log10(infection rate) (/min)") +
          ylab("Burst size") +
          guides(fill = FALSE, shape = FALSE) +
          labs(subtitle = "Lysis time (min)") +
          NULL)
  dev.off()
}

#Run 1: max time - exinction time ----
if (glob_make_statplots) {
  #Extin time vs maxtime (zoomed, one facet)
  tiff("./plots/run1_extintime_maxtime_k1e9_u011.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                y_summarized1$u_S == 0.011, ],
         aes(x = max_time/60, y = extin_time/60,
             color = as.factor(a), shape = maxdens_k_ratio < 0.95)) +
    geom_point(size = 1.5, alpha = 0.6) +
    #facet_grid(u_S ~ k_S) +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis(name = "Infection\nRate\n(/min)", discrete = TRUE) +
    coord_cartesian(ylim = c(0, 22), xlim = c(0, 16)) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Bacterial Extinction Time (hrs)") +
    guides(shape = FALSE)
  dev.off()
  
  #Extin time vs maxtime (all data)
  tiff("./plots/run1_extintime_maxtime.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[!is.na(y_summarized1$max_time) &
                                !is.na(y_summarized1$extin_time), ],
         aes(x = max_time/60, y = extin_time/60,
             color = as.factor(a), shape = maxdens_k_ratio < 0.95)) +
    geom_point() +
    facet_grid(u_S ~ k_S, scales = "free") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis(name = "Infection\nRate\n(/min)", discrete = TRUE) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Bacterial Extinction Time (hrs)") +
    guides(shape = FALSE)
  dev.off()
  
  tiff("./plots/run1_extintime_maxtime_zoomed.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[y_summarized1$maxdens_k_ratio < 0.95, ],
         aes(x = max_time/60, y = extin_time/60,
             color = as.factor(a), shape = maxdens_k_ratio < 0.95)) +
    geom_point() +
    facet_grid(u_S ~ k_S, scales = "free") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis(name = "Infection\nRate\n(/min)", discrete = TRUE) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Bacterial Extinction Time (hrs)") +
    guides(shape = FALSE)
  dev.off()
  
  #Extintime-maxtime ratio vs maxtime (zoomed, one facet)
  tiff("./plots/run1_extinmaxratio_maxtime_k1e9_u011.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                       y_summarized1$u_S == 0.011, ],
         aes(x = max_time/60, y = extin_time/max_time,
             color = as.factor(a), shape = maxdens_k_ratio < 0.95)) +
    geom_point() +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis(name = "Infection\nRate\n(/min)", discrete = TRUE) +
    coord_cartesian(ylim = c(NA, 1.65), xlim = c(0, 16)) +
    guides(shape = FALSE) +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Extinction Time to Peak Density Time Ratio") +
    NULL
  dev.off()
  
  #Extintime-maxtime ratio vs maxtime (all data)
  ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                y_summarized1$u_S == 0.011, ],
         aes(x = max_time/60, y = extin_time/max_time,
             color = as.factor(a), shape = maxdens_k_ratio < 0.95)) +
    geom_point() +
    facet_grid(u_S ~ k_S, scales = "free") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis(name = "Infection\nRate\n(/min)", discrete = TRUE) +
    coord_cartesian(ylim = c(NA, 1.65), xlim = c(0, 16)) +
    guides(shape = FALSE) +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Ratio of Extinction Time to Peak Density Time") +
    NULL
  
  ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                y_summarized1$u_S == 0.011, ],
         aes(x = max_time/60, y = extin_time/max_time,
             color = as.factor(a), shape = maxdens_k_ratio < 0.95)) +
    geom_point() +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    coord_cartesian(ylim = c(NA, 1.8), xlim = c(0, 18)) +
    #facet_grid(u_S ~ k_S, scales = "free") +
    NULL
  
  
  
  
  tiff("./plots/run1_extin-max_maxtime_k1e9_u011.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                y_summarized1$u_S == 0.011, ],
         aes(x = max_time/60, y = (extin_time-max_time)/60,
             color = as.factor(a), shape = maxdens_k_ratio < 0.95)) +
    geom_point() +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis(name = "Infection\nRate\n(/min)", discrete = TRUE) +
    #facet_grid(u_S ~ k_S, scales = "free") +
    coord_cartesian(ylim = c(NA, 6), xlim = c(0, 16)) +
    guides(shape = FALSE) +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Extinction Time - Peak Density Time") +
    NULL
  dev.off()
  
  ggplot(data = y_summarized1[y_summarized1$max_dens < 0.95*y_summarized1$k_S, ],
         aes(x = max_time/60, y = phage_atmaxdens,
             color = as.factor(a), shape = as.factor(b))) +
    geom_point() +
    scale_y_continuous(trans = "log10") +
    facet_grid(u_S ~ k_S, scales = "free") +
    #coord_cartesian(ylim = c(0, 36), xlim = c(0, 24)) +
    NULL
  
  tiff("./plots/run1_extintime_phageatmax_k1e9_u011.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                y_summarized1$u_S == 0.011 &
                                y_summarized1$maxdens_k_ratio < 0.95, ],
         aes(x = phage_atmaxdens, y = extin_time/max_time,
             color = as.factor(a), shape = as.factor(b))) +
    geom_point() +
    scale_x_continuous(trans = "log10") +
    scale_color_viridis(name = "Infection\nRate\n(/min)", discrete = TRUE) +
    scale_shape_discrete(name = "Burst\nSize") +
    #facet_grid(u_S ~ k_S, scales = "free") +
    coord_cartesian(ylim = c(NA, 1.65)) +
    labs(x = "Phage Density at Peak Bacterial Density Time (hrs)",
         y = "Ratio of Extinction Time to Peak Density Time") +
    NULL
  dev.off()
  
  ggplot(data = y_summarized1[y_summarized1$max_dens < 0.95*y_summarized1$k_S, ],
         aes(x = phage_atmaxdens, y = (extin_time - max_time)/60,
             color = as.factor(a), shape = as.factor(b))) +
    geom_point() +
    scale_x_continuous(trans = "log10") +
    facet_grid(u_S ~ k_S, scales = "free") +
    #coord_cartesian(ylim = c(0, 36), xlim = c(0, 24)) +
    NULL
}

#Run 1: auc - max time - extin time ----
if (glob_make_statplots) {
  #Auc vs maxtime (zoomed, one facet)
  # (with prediction line from analytical integral
  # of logistic curve up to peak time)
  tiff("./plots/run1_auc_maxtime_k1e9_u011_zoomed.tiff",
  width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                y_summarized1$u_S == 0.011, ],
         aes(x = max_time/60, y = auc/60)) +
    geom_point(aes(shape = maxdens_k_ratio < 0.95, color = as.factor(a)),
               size = 1.5, alpha = 0.6) +
    geom_line(aes(y = pred_auc/60), color = "black") +
    scale_y_continuous(trans = "log10") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
    coord_cartesian(xlim = c(0, 16), ylim = c(10**6, 10**10)) +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Area Under the Curve (hr cfu/mL)") +
    guides(shape = FALSE) +
    theme_bw() +
    NULL
  dev.off()
  
  #Auc vs maxtime (one facet, not zoomed)
  tiff("./plots/run1_auc_maxtime_k1e9_u011.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[y_summarized1$k_S == 10**9 &
                                y_summarized1$u_S == 0.011, ],
         aes(x = max_time/60, y = auc/60, color = as.factor(a),
             shape = maxdens_k_ratio < 0.95)) +
    geom_point(size = 1.5, alpha = 0.6) +
    #facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Area Under the Curve (hr cfu/mL)") +
    guides(shape = FALSE) +
    theme_bw() +
    NULL
  dev.off()
  
  #Auc vs maxtime (all facets, not zoomed)
  tiff("./plots/run1_auc_maxtime.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1,
         aes(x = max_time/60, y = auc/60)) +
    geom_point(aes(shape = maxdens_k_ratio < 0.95, color = as.factor(a)),
               size = 1.5, alpha = 0.6) +
    geom_line(aes(y = pred_auc/60), color = "black") +
    facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
    labs(x = "Peak Bacterial Density Time (hrs)",
         y = "Area Under the Curve (hr cfu/mL)") +
    guides(shape = FALSE) +
    theme_bw() +
    NULL
  dev.off()
  
  #Auc vs extin time (all facets, not zoomed)
  tiff("./plots/run1_auc_extintime.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1,
         aes(x = extin_time/60, y = auc/60, color = as.factor(a),
             shape = maxdens_k_ratio < 0.95)) +
    geom_point(size = 1.5, alpha = 0.6) +
    facet_grid(k_S ~ u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
    labs(x = "Bacterial Extinction Time (hrs)",
         y = "Area Under the Curve (hr cfu/mL)") +
    guides(shape = FALSE) +
    theme_bw() +
    NULL
  dev.off()
  
  #Auc vs max dens (all facets, not zoomed)
  tiff("./plots/run1_auc_maxdens.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1,
         aes(x = max_dens, y = auc/60, color = as.factor(a),
             shape = maxdens_k_ratio < 0.95)) +
    geom_point(size = 1.5, alpha = 0.6) +
    facet_grid(u_S ~ k_S, scales = "free") +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
    labs(x = "Peak Bacterial Density (cfu/mL)",
         y = "Area Under the Curve (hr cfu/mL)") +
    guides(shape = FALSE) +
    theme_bw() +
    NULL
  dev.off()
  
  #Auc (total) vs empirical auc up to peak
  ggplot(data = y_summarized1[y_summarized1$maxdens_k_ratio < 0.95, ],
         aes(x = auc, y = auc_tomaxtime, color = as.factor(a),
             shape = maxdens_k_ratio < 0.95)) +
    geom_point() +
    facet_grid(u_S ~ k_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    geom_abline(slope = 1, intercept = 0) +
    scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    NULL
  
  #Auc to predicted auc ratio vs auc (all data)
  ggplot(data = y_summarized1,
         aes(x = auc, y = auc/pred_auc,
             shape = maxdens_k_ratio < 0.95)) +
    geom_point() +
    facet_grid(u_S ~ k_S, scales = "free") +
    scale_x_continuous(trans = "log10") +
    scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    NULL
  
  #Auc to predicted auc ratio vs auc (just maxdens < 0.95 * k)
  tiff("./plots/run1_ratio_auc_to_predauc_zoomed.tiff",
  width = 5, height = 4, units = "in", res = 300)
  ggplot(data = y_summarized1[y_summarized1$maxdens_k_ratio < 0.95, ],
         aes(x = auc, y = auc/pred_auc,
             color = as.factor(a))) +
    geom_point() +
    facet_grid(u_S ~ k_S, scales = "free") +
    scale_x_continuous(trans = "log10") +
    scale_color_viridis_d(name = "Infection\nRate\n(/min)") +
    scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
    labs(x = "Area Under the Curve (hr cfu/mL)",
         y = "Ratio of Total AUC: AUC at Peak Bacterial Density") +
    theme_bw() +
    NULL
  dev.off()
  
  summary(
    lm(auc/pred_auc ~ log10(auc) * as.factor(k_S) * as.factor(u_S),
             y_summarized1[y_summarized1$maxdens_k_ratio < 0.95, ]))
}


##Run 2: bact variants simulation & summarization ----
run2 <- run_sims_filewrapper(name = "run2",
                     u_Svals = signif(2**c(-5, -5.5, -6, -6.5, -7), 2),
                     k_Svals = 10**c(7, 8, 9, 10, 11),
                     avals = 10**c(-8, -9, -10, -11, -12),
                     tauvals = signif(10**c(1.4, 1.8), 2),
                     bvals = signif(10**c(1.5, 2.5), 2),
                     c_SIvals = 1,
                     init_S_dens_vals = 10**6,
                     init_moi_vals = 10**-2,
                     min_dens = 100,
                     init_time = 100,
                     init_stepsize = 1,
                     print_info = TRUE,
                     read_file = glob_read_files)

ybig2 <- run2[[1]]

ybig2 <- group_by_at(ybig2, .vars = 1:17)

#Summarize runs (using ifelse to handle potential non-equil conditions)
y_summarized2 <- 
  dplyr::summarize(
    ybig2,
    max_dens = max(Density[Pop == "B"]),
    max_time = time[Pop == "B" & Density[Pop == "B"] == max_dens],
    #Make references for finding extin point
    extin_dens = 10**4,
    #Technically this is the first index after extinction
    extin_index =
      ifelse(any(Pop == "B" & Density <= extin_dens & time >= max_time),
             min(which(Pop == "B" & Density <= extin_dens & time >= max_time)),
             NA),
    extin_index_back1 =
      ifelse(is.na(extin_index), NA,
             which(Pop == "B")[match(extin_index, which(Pop == "B")) - 1]),
    #use linear interpolation to find extin time
    extin_time =
      ifelse(is.na(extin_index), NA,
             time[extin_index] - (Density[extin_index] - extin_dens)*
               (time[extin_index] - time[extin_index_back1])/
               (Density[extin_index] - Density[extin_index_back1])),
    #make references for auc (here indices are within Pop == "B")
    extin_index_winB =
      ifelse(is.na(extin_time), NA, which(which(Pop == "B") == extin_index)),
    extin_index_back1_winB =
      ifelse(is.na(extin_time), NA, which(which(Pop == "B") == extin_index)-1),
    extin_index_back2_winB =
      ifelse(is.na(extin_time), NA, which(which(Pop == "B") == extin_index)-2),
    #using trapezoid rule to find auc
    auc =
      ifelse(is.na(extin_time), NA,
             #trapezoids of all intervals before one ending at extin time
             (time[Pop == "B"][2] - time[Pop == "B"][1])/2 *
               (Density[Pop == "B"][1] +
                  2*sum(Density[Pop == "B"][2:extin_index_back2_winB]) +
                  Density[Pop == "B"][extin_index_back1_winB]) +
               #trapezoid that ends at extin time
               (extin_time-time[extin_index_back1])/2 *
               (Density[extin_index_back1] + extin_dens)),
    phage_final = ifelse(equil[1], max(Density[Pop == "P"]), NA),
    #make references for phage extin dens
    phage_dens_y1 = ifelse(is.na(extin_time), NA,
                           Density[max(which(Pop == "P" & time < extin_time))]),
    phage_dens_y2 = ifelse(is.na(extin_time), NA,
                           Density[which(Pop == "P" & time == time[extin_index])]),
    phage_time_x1 = ifelse(is.na(extin_time), NA,
                           time[max(which(Pop == "P" & time < extin_time))]),
    phage_time_x2 = ifelse(is.na(extin_time), NA,
                           time[which(Pop == "P" & time == time[extin_index])]),
    #using linear interpolation to find phage dens at bact extinction
    phage_extin = ifelse(is.na(extin_time), NA,
                         phage_dens_y1 +
                           (phage_dens_y2-phage_dens_y1)/(phage_time_x2-phage_time_x1)*
                           (extin_time-phage_time_x1)),
    # phage_r =
    #   (log(phage_final)- log(init_S_dens[1]*init_moi[1]))/extin_time,
    # phage_r_extin =
    #   (log(phage_extin) - log(init_S_dens[1]*init_moi[1]))/extin_time,
    phage_atmaxdens = Density[Pop == "P" & time == max_time],
    maxdens_k_ratio = max_dens/k_S[1],
    run_time = max(time),
    equil = equil[1]
  )

#Drop unneeded columns
y_summarized2 <- subset(y_summarized2,
                        select = -c(extin_index, extin_index_back1,
                                    extin_index_winB, extin_index_back1_winB,
                                    extin_index_back2_winB, phage_dens_y1,
                                    phage_dens_y2, phage_time_x1,
                                    phage_time_x2))
y_summarized2 <- as.data.frame(y_summarized2)

#Run 2: density dynamics ----
dir.create("./run2_noequil/", showWarnings = F)
if (glob_make_curveplots) {
  for (myrun in unique(ybig2$uniq_run[ybig2$equil == FALSE])) {
    tiff(paste("./run2_noequil/", myrun, ".tiff", sep = ""),
         width = 4, height = 4, units = "in", res = 200)
    print(ggplot(data = ybig2[ybig2$uniq_run == myrun &
                                ybig2$Pop %in% c("S", "I", "P"), ],
                 aes(x = time, y = Density+1, color = Pop)) +
            geom_line(lwd = 1.5, alpha = 0.5) +
            scale_y_continuous(trans = "log10"))
    dev.off()
  }
}

dir.create("./run2_equil/", showWarnings = F)
if (glob_make_curveplots) {
  for (myrun in unique(ybig2$uniq_run[ybig2$equil == TRUE])) {
    tiff(paste("./run2_equil/", myrun, ".tiff", sep = ""),
         width = 4, height = 4, units = "in", res = 200)
    print(ggplot(data = ybig2[ybig2$uniq_run == myrun &
                                ybig2$Pop %in% c("S", "I", "P"), ],
                 aes(x = time, y = Density+1, color = Pop)) +
            geom_line(lwd = 1.5, alpha = 0.5) +
            scale_y_continuous(trans = "log10"))
    dev.off()
  }
}

#Run 2: landscape of peak time ----
#Code for 2D contours
if(glob_make_statplots) {
  tiff("./plots/run2_peaktime_contour_k1e9_b32_tau63.tiff", width = 5, height = 4,
       units = "in", res = 300)
  print(ggplot(data = y_summarized2[y_summarized2$k_S == 10**9 &
                                      y_summarized2$b == 32 &
                                      y_summarized2$tau == 63, ],
               aes(x = log10(a), y = u_S)) +
          geom_contour_filled(aes(z = max_time/60), alpha = 0.5) +
          geom_point(aes(color = max_time/60, shape = maxdens_k_ratio < 0.95),
                     size = 3) +
          scale_color_viridis(name = "Peak time (hr)", 
                              breaks = c(6, 12, 18, 24, 30)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          scale_y_continuous(trans = "log10") +
          xlab("Log10(infection rate) (/min)") +
          ylab("Bacterial growth rate (/min)") +
          guides(fill = FALSE, shape = FALSE) +
          NULL)
  dev.off()
  
  #k_S in facets
  tiff("./plots/run2_peaktime_contour_facetk_S.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(ggplot(data = y_summarized2,
               aes(x = log10(a), y = u_S)) +
          geom_contour_filled(aes(z = max_time), alpha = 0.7) +
          geom_point(aes(color = max_time/60, shape = maxdens_k_ratio < 0.95),
                     size = 1.5) +
          scale_color_viridis(name = "Peak time (hr)",
                              breaks = c(6, 12, 18, 24, 30, 36)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          facet_grid(b*tau~k_S) +
          scale_y_continuous(trans = "log10") +
          xlab("Log10(infection rate) (/min)") +
          ylab("Bacterial growth rate (/min)") +
          guides(fill = FALSE, shape = FALSE) +
          labs(subtitle = "Carrying capacity (cfu/mL)") +
          NULL)
  dev.off()
  
  #a in facets
  tiff("./plots/run2_peaktime_contour_faceta.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(ggplot(data = y_summarized2,
               aes(x = u_S, y = k_S)) +
          geom_contour_filled(aes(z = max_time), alpha = 0.7) +
          geom_point(aes(color = max_time/60, shape = maxdens_k_ratio < 0.95),
                     size = 1.5) +
          scale_color_viridis(name = "Peak time (hr)",
                              breaks = c(6, 12, 18, 24, 30, 36)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          facet_grid(b*tau~log10(a)) +
          scale_x_continuous(trans = "log10") +
          scale_y_continuous(trans = "log10") +
          xlab("Bacterial growth rate (/min)") +
          ylab("Carrying capacity (cfu/mL)") +
          guides(fill = FALSE, shape = FALSE) +
          labs(subtitle = "Log10(infection rate) (/min)") +
          NULL)
  dev.off()
  
  #u_S in facets
  tiff("./plots/run2_peaktime_contour_facetu_S.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(ggplot(data = y_summarized2,
               aes(x = log10(a), y = k_S)) +
          geom_contour_filled(aes(z = max_time), alpha = 0.7) +
          geom_point(aes(color = max_time/60, shape = maxdens_k_ratio < 0.95),
                     size = 1.5) +
          scale_color_viridis(name = "Peak time (hr)",
                              breaks = c(6, 12, 18, 24, 30, 36)) +
          scale_shape_manual(breaks = c(TRUE, FALSE), values = c(16, 4)) +
          facet_grid(b*tau~u_S) +
          scale_y_continuous(trans = "log10") +
          xlab("Log10(infection rate) (/min)") +
          ylab("Carrying capacity (cfu/mL)") +
          guides(fill = FALSE, shape = FALSE) +
          labs(subtitle = "Bacterial growth rate (/min)") +
          NULL)
  dev.off()
}