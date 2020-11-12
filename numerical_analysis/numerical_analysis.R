#TODO:
# verify that coinfection actually works the way it should!
# should probably come up with a way of not plotting
#   things where bact density reaches ~k
# check whether burst size can be inferred by running curves across 
#   various init conditions then plotting max dens vs final phage
# then think about how that inter-relates to extinction time – phage r 
#   relationship
# Perhaps should measure something like P density at bact peak? 
#   Or P r up to bact peak?
# figure out how to calculate area-under-curve consistently
# use netrd (see Scarpino mtg)
#
# given control curves and curves varying in density & moi,
#   figure out how one could calculate phage parameters
#   or phage fitness
# can dede itself handle a stop-at-equilibrium condition?
#   Yes, they're called roots. However, it's not clear whether
#     it would really make things faster or not
# add in the evolution of resistant bacteria
# fix area-under-curve calculations
#     doing so requires implementing static timestep reporting
# use deriv-percap of B to find max decay rate?
# think about characteristics to be taken when curves are too slow
#   (e.g. time to grow above some threshold density)
# test whether it's better to average replicate wells first, then analyze
#   or analyze each independently then average stats together
# When stats are plotted against each other there are lots of logistic curves
#  fit a logistic eq to them and see what they are
#  (e.g. are they just the bacterial curve when no phage around?)
# Go back and re-learn multivariate calc
# In theory we should be able to predict how much pfu_final
#  is above b*max_dens based on auc before max_dens (or similar)
# Should calculate bacterial decay rate from max_dens to extin_time
#  try as negative exponential growth from max_dens
#  e.g. B(t) = max_dens - e^(r(t-max_time))
# Perhaps the whole B density curve can be reduced to two logistic-like curves?
#  One representing bacterial growth, and one representing bacterial decay
#  We're kind of idealizing bacterial growth as logistic and phage
#   growth as logistic and B(t) is just the integrated difference between them
# Parallelize for running on cluster with much faster walltime
# Measure amount of time I is above S (and so I is more of B) – should be correlated w/ tau
# Run w 2 P pops in competition. Compare outcome to indiv grow curves
# Think about similar approach as Wang et al ’96 for our model
# Think about superinfection parameter being included
#   (simply a rescaling of a for already-infected cells)

## Import libraries ----

library(deSolve)
#library(reshape2)
library(data.table)
library(ggplot2)
library(dplyr)

#Setwd
mywd_split <- strsplit(getwd(), split = "/") 
if (mywd_split[[1]][length(mywd_split[[1]])] != "numerical_analysis") {
  dir.create("numerical_analysis", showWarnings = FALSE)
  setwd("./numerical_analysis/")
}

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

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
  dY <- c(S = 0, I = 0, P = 0, R = 0)
  
  ##Calculate dS
  
  #V1 
  # note: exponential dS/dt
  #dS/dt = rS - aSP
  #dS <- parms["r"] * y["S"] - a * y["S"] * y["P"]
  
  #V2 
  # added: changed dS/dt growth from exp to logistic
  #dS/dt = rS((K-S)/K) - aSP
  # dY["S"] <- parms["r"] * y["S"] * ((parms["K"] - y["S"])/parms["K"]) - 
  #   parms["a"] * y["S"] * y["P"]
  
  #V3
  # added: competition from I pop on S pop
  #dS/dt = rS((K-S-c*I)/K) - aSP
  # dY["S"] <- parms["r"] * y["S"] * 
  #     ((parms["K"] - y["S"] - parms["c"] * y["I"])/parms["K"]) - 
  #   parms["a"] * y["S"] * y["P"]
  #dP/dt = baS(t-tau)P(t-tau) - aSP
  # if (t < parms["tau"]) {
  #   dY["P"] <- -parms["a"] * y["S"] * y["P"]
  # } else {
  #   dY["P"] <- parms["b"] * parms["a"] * 
  #     lagvalue(t-parms["tau"], 1)*lagvalue(t-parms["tau"], 3) - 
  #     parms["a"]*y["S"]*y["P"]
  # }
  
  #V4
  # added:  superinfection
  #         resistant population
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
  
  if(init_time %% init_stepsize != 0) {
    warning("init_time is not divisible by init_stepsize, this has not been tested")
  }
  
  #V3 to V4 argument changes
  # rvals = u_Svals
  # kvals = k_Svals
  # cvals = c_SIvals
  # init_bact_dens_vals = init_S_dens_vals
  # added: zvals (superinfection [0,1]), u_Rvals, k_Rvals, mvals, 
  # c_SRvals, c_RIvals, c_RSvals, init_R_dens_vals
  
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
        stop("dynamic_stepsize = FALSE does not yet work")
        #If dynamic_stepsize false, keep stepsize at init_stepsize
        #times <- seq(0, init_time*2**j, init_stepsize)
      }
      
      #Run simulation (using 1st entry of list as error code:\
      # 0 - success
      # 1 - error
      # 2 - warning
      yout_list <- tryCatch(
        expr = {
          #Note that the max step size for the integrator is the 
          # same as our step size except halved for each k count
          list(0,
               as.data.frame(
                 dede(y = yinit, times = times, func = derivs, 
                      parms = params, hmax = init_stepsize*2**(j-k))))
        },
        error = function(e) {list(1)},
        warning = function(w) {
          list(2,
               as.data.frame(
                 dede(y = yinit, times = times, func = derivs, 
                      parms = params, hmax = 2**(j-k))))
        }
      )
      
      #Infinite loop prevention check (j = 10 is 24 hrs)
      if (j >= 10 | k >= 15 | j+k >= 20) {
        keep_running <- FALSE
        at_equil <- FALSE
      }
      
      #If there was an error, increase k by 1 and re-run
      if(yout_list[[1]] == 1) {
        k <- k+1
        #If there was a warning, could be several causes, so we
        # generally just halve step size and increase length
      } else if (yout_list[[1]] == 2) {
        j <- j+1
        k <- k+2
        #If it was successful, check for equilibrium
      } else if (yout_list[[1]] == 0) {
        #First drop all rows with nan
        yout_list[[2]] <- yout_list[[2]][!(is.nan(yout_list[[2]]$S) |
                                             is.nan(yout_list[[2]]$I) |
                                             is.nan(yout_list[[2]]$P) |
                                             is.nan(yout_list[[2]]$R)), ]
        
        #S and I both at equil, we're done
        if (yout_list[[2]]$S[nrow(yout_list[[2]])] < min_dens & 
            yout_list[[2]]$I[nrow(yout_list[[2]])] < min_dens) {
          keep_running <- FALSE
          at_equil <- TRUE
        #S not at equil, need more time
        } else if (yout_list[[2]]$S[nrow(yout_list[[2]])] >= min_dens) { 
          j <- j+1
        #I not at equil (but S is because above check failed),
        #   first we'll lengthen the simulation
        #    (to make sure it was long enough to catch the last burst)
        #   then we'll start shrinking our step size
        } else if (yout_list[[2]]$I[nrow(yout_list[[2]])] >= min_dens) {
          if (i_only_pos_times < 1) {
            j <- j+1
            i_only_pos_times <- i_only_pos_times+1
          } else {
            k <- k+1
          }
        }
      }
    }
    
    #Once end conditions triggered, if run succeeded
    if(yout_list[[1]] == 0 | yout_list[[1]] == 2) {
      #Calculate all bacteria (B)
      yout_list[[2]]$B <- yout_list[[2]]$S + yout_list[[2]]$I + yout_list[[2]]$R
      #Calculate all phage (PI)
      yout_list[[2]]$PI <- yout_list[[2]]$P + yout_list[[2]]$I
      
      #Reshape, add parameters, and fill into ybig in right rows
      ybig[((i-1)*6*(1+init_time/init_stepsize)+1):
             ((i)*6*(1+init_time/init_stepsize)), ] <- 
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
              data.table::melt(data = data.table::as.data.table(yout_list[[2]]), 
                               id.vars = c("time"),
                               value.name = "Density", 
                               variable.name = "Pop",
                               variable.factor = FALSE))
      #If the run failed
    } else {
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
    }
    
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
  
## Define function that calculates derivatives ----
calc_deriv <- function(density, percapita = FALSE,
                       subset_by = NULL, time = NULL,
                       time_normalize = NULL) {
  #Note! density values must be sorted sequentially into their unique sets already
  
  #Provided a vector of density values, this function returns (by default) the
  # difference between sequential values
  #if percapita = TRUE, the differences of density are divided by density
  #if subset_by is provided, it should be a vector (same length as density),
  # the unique values of which will separate calculations
  #if time_normalize is specified, time should be provided as a simple 
  # numeric (e.g. number of seconds) in some unit
  #Then the difference will be normalized for the time_normalize value
  #(e.g. if time is provided in seconds and the difference per hour is wanted,
  # time_normalize should = 3600)
  
  #Check inputs
  if (!is.numeric(time)) {
    stop("time is not numeric")
  }
  if (!is.null(time_normalize)) {
    if (!is.numeric(time_normalize)) {
      stop("time_normalize is not numeric")
    } else if (is.null(time)) {
      stop("time_normalize is specified, but time is not provided")
    }
  }
  
  #Calc derivative
  ans <- c(density[2:length(density)]-density[1:(length(density)-1)])
  #Percapita (if specified)
  if (percapita) {
    ans <- ans/density[1:(length(density)-1)]
  }
  #Time normalize (if specified)
  if (!is.null(time_normalize)) {
    ans <- ans/
      (c(time[2:length(time)]-time[1:(length(time)-1)])/time_normalize)
  }
  #Subset by (if specified)
  if (!is.null(subset_by)) {
    ans[subset_by[2:length(subset_by)] != subset_by[1:(length(subset_by)-1)]] <- NA
  }
  return(c(ans, NA))
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


## Lit review of parameters ----
#r ranges from .04/min (17 min doubling time)
#         to 0.007/min (90 min doubling time)
#   ln(1/2) = -r * doub_time
#K ranges from say 10^7 to 10^9
#adsorption ranges from 1x10^-12 to 1x10^-8 /min
#Lysis time ranges from 10 to 105 mins
#Burst size ranges from 4.5 to 1000
#Phage resistance mutation rate from 10^-8 to 10^-2
#     with 10^-7 to 10^-4 most common range

## Global Settings ----
glob_read_files <- TRUE
glob_make_curveplots <- TRUE
glob_make_statplots <- TRUE

## Run #1: a, b, tau (phage traits) ----
run1 <- run_sims_filewrapper(name = "run1",
                             u_Svals = c(0.023), #(30 min doubling time)
                             k_Svals = c(10**9),
                             avals = 10**seq(from = -12, to = -8, by = 1),
                             tauvals = signif(10**seq(from = 1, to = 2, by = 0.25), 3),
                             bvals = signif(5*10**seq(from = 0, to = 2, by = 0.5), 3),
                             c_SIvals = 1,
                             init_S_dens_vals = 10**6,
                             init_moi_vals = 10**-2,
                             min_dens = 0.1,
                             init_time = 100,
                             init_stepsize = 1,
                             print_info = TRUE,
                             read_file = glob_read_files)

#Find peaks & extinction via summarize
ybig1 <- group_by_at(run1[[1]], .vars = 1:17)
y_summarized1 <- summarize(ybig1,
                          max_dens = max(Density[Pop == "B"]),
                          max_time = time[Pop == "B" & 
                                            Density[Pop == "B"] == max_dens],
                          extin_index = min(which(Pop == "B" &
                                                    Density <= 10**4)),
                          extin_dens = Density[extin_index],
                          extin_time = time[extin_index],
                          auc = sum(Density[Pop == "B" & time < extin_time])*
                            extin_time,
                          phage_final = max(Density[Pop == "P"]),
                          phage_extin = Density[Pop == "P" & time == extin_time],
                          phage_r = (log(phage_final)-
                                       log(init_S_dens[1]*init_moi[1]))/
                            extin_time,
                          run_time = max(time)
)
                          

#Calculate derivatives
ybig1$deriv <- calc_deriv(density = ybig1$Density, 
                           percapita = FALSE,
                           subset_by = paste(ybig1$uniq_run, ybig1$Pop), 
                           time = ybig1$time,
                           time_normalize = 60)
ybig1$deriv_percap <- calc_deriv(density = ybig1$Density, 
                                percapita = TRUE,
                                subset_by = paste(ybig1$uniq_run, ybig1$Pop), 
                                time = ybig1$time,
                                time_normalize = 60)

#Make plots of density against time ----
dir.create("run1_dens_curves", showWarnings = FALSE)
if (glob_make_curveplots) {
  dens_offset <- 10
  for (run in unique(ybig1$uniq_run)) {
    tiff(paste("./run1_dens_curves/", run, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(
      ggplot(data = ybig1[ybig1$uniq_run == run &
                           ybig1$Pop %in% c("S", "I", "P"),], 
             aes(x = time, y = Density+dens_offset, color = as.factor(Pop))) +
              geom_line(lwd = 1.5, alpha = 1) + 
        geom_line(data = ybig1[ybig1$uniq_run == run &
                                ybig1$Pop == "B",], 
                  aes(x = time, y = Density+dens_offset),
                  color = "black", alpha = 0.5, lwd = 1.1) +
        geom_line(data = ybig1[ybig1$uniq_run == run &
                                ybig1$Pop == "PI",],
                  aes(x = time, y = Density+dens_offset),
                  color = "black", alpha = 0.5, lwd = 1, lty = 3) +
        geom_point(data = y_summarized1[y_summarized1$uniq_run == run, ],
                   aes(x = max_time, y = max_dens+dens_offset), color = "black") +
        geom_point(data = y_summarized1[y_summarized1$uniq_run == run, ],
                   aes(x = extin_time, y = extin_dens+dens_offset), color = "black") +
        scale_y_continuous(trans = "log10") +
        scale_x_continuous(breaks = seq(from = 0, to = max(ybig1$time), 
                                        by = round(max(ybig1[ybig1$uniq_run == run &
                                                            ybig1$Pop != "B", 
                                                            "time"])/10))) +
        scale_color_manual(limits = c("S", "I", "P"),
                           values = my_cols[c(2, 3, 1)]) +
        geom_hline(yintercept = 10, lty = 2) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              title = element_text(size = 9)) +
        ggtitle(paste(ybig1[min(which(ybig1$uniq_run == run)), 6:8],
                      collapse = ", ")) +
        labs(y = paste("Density +", dens_offset)) +
        NULL
    )
    dev.off()
  }
}

#Plot summarized statistics ----
dir.create("run1_statplots", showWarnings = FALSE)
if (glob_make_statplots) {
  for (stat in c("max_dens", "max_time", "extin_time", 
                 "auc", "phage_final", "phage_r")) {
    tiff(paste("./run1_statplots/", stat, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(ggplot(data = y_summarized1,
                 aes(x = a, y = get(stat), color = as.factor(b), 
                     group = as.factor(b))) + 
            geom_point(size = 3, alpha = 0.8) + 
            geom_line(size = 1.1, alpha = 0.6) +
            facet_grid(~tau) +
            labs(y = stat) +
            scale_y_continuous(trans = "log10") +
            scale_x_continuous(trans = "log10") +
            scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggtitle("tau") +
            NULL)
    dev.off()
  }
}

y_sum_melt1 <- reshape2::melt(y_summarized1,
                   id.vars = 1:17,
                   variable.name = "sum_stat",
                   value.name = "stat_val")

if (glob_make_statplots) {
  tiff("./run1_statplots/all_stats.tiff",
       width = 5, height = 6, units = "in", res = 300)
  print(ggplot(data = y_sum_melt1[y_sum_melt1$sum_stat %in%
                              c("max_dens", "max_time", "extin_time", 
                                #"auc", 
                                "phage_final", 
                                "phage_r"
                                ), ],
         aes(x = a, y = stat_val, color = as.factor(b), group = as.factor(b))) +
    geom_point(size = 1.5, alpha = 0.8) + 
    geom_line(size = 1.1, alpha = 0.6) +
    facet_grid(sum_stat~tau, scales = "free_y") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(size = 10)) +
    ggtitle("tau") +
    NULL)
  dev.off()

  tiff("./run1_statplots/all_stats2.tiff",
       width = 5, height = 5, units = "in", res = 300)
  print(ggplot(data = y_sum_melt1[y_sum_melt1$sum_stat != "extin_dens", ],
         aes(x = b, y = stat_val, color = as.factor(tau), group = as.factor(tau))) +
    geom_point(size = 2, alpha = 0.8) + 
    geom_line(size = 1.1, alpha = 0.6) +
    facet_grid(sum_stat~a, scales = "free_y") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("a") +
    NULL)
  dev.off()

  tiff("./run1_statplots/all_stats3.tiff",
       width = 5, height = 5, units = "in", res = 300)
  print(ggplot(data = y_sum_melt1[y_sum_melt1$sum_stat != "extin_dens", ],
         aes(x = tau, y = stat_val, color = as.factor(a), group = as.factor(a))) +
    geom_point(size = 2, alpha = 0.8) + 
    geom_line(size = 1.1, alpha = 0.6) +
    facet_grid(sum_stat~b, scales = "free_y") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("b") +
    NULL)
  dev.off()
}

#Plot stats against ea other ----

#First plot all at once
for (col in c("max_dens", "max_time",
              "extin_time", "auc", "phage_final", "phage_r")) {
  y_summarized1[, paste(col, "_log10", sep = "")] <- log10(y_summarized1[, col])
}

if (glob_make_statplots) {
  tiff("./run1_statplots/stat_cors.tiff", width = 10, height = 10, units = "in", res = 300)
  #Make base figure
  p <- GGally::ggpairs(y_summarized1,
                     aes(color = as.factor(b), shape = as.factor(a)),
                     columns = c("max_dens_log10", "max_time_log10",
                                 "extin_time_log10", 
                                 #"auc_log10", 
                                 "phage_final_log10", 
                                 "phage_r_log10"),
                     lower = list(continuous = "points"),
                     upper = list(continuous = "points")) +
  theme_bw() +
  theme(strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(p)
  dev.off()

  #Then make indiv paired plots
  tiff("./run1_statplots/maxdens_maxtime.tiff",
       width = 5, height = 5, units = "in", res = 300)
  print(ggplot(data = y_summarized1,
         aes(x = max_time, y = max_dens, color = as.factor(b), 
             fill = as.factor(b), shape = as.factor(a))) +
    geom_point(size = 2.5, alpha = 0.5) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    scale_shape_manual(values = 21:25) +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    scale_fill_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    NULL)
  dev.off()
  
  tiff("./run1_statplots/maxdens_maxtime_facet.tiff",
       width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized1,
         aes(x = max_time, y = max_dens, color = as.factor(b), 
             fill = as.factor(b), shape = as.factor(a))) +
    geom_point(size = 2.5, alpha = 0.5) +
    facet_grid(~tau) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    scale_shape_manual(values = 21:25) +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    scale_fill_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("tau") +
    NULL)
  dev.off()
  
  tiff("./run1_statplots/maxdens_extintime.tiff",
       width = 5, height = 5, units = "in", res = 300)
  print(ggplot(data = y_summarized1,
         aes(x = extin_time, y = max_dens, color = as.factor(b), 
             fill = as.factor(tau), shape = as.factor(a))) +
    geom_point(size = 1.5, alpha = 1, stroke = 1) +
    #facet_grid(tau~.) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    scale_shape_manual(values = 21:25) +
    scale_color_manual(values = my_cols) +
    scale_fill_manual(values = my_cols) +
    #scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    #scale_fill_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    NULL)
  dev.off()
  
  tiff("./run1_statplots/maxdens_extintime_facet.tiff",
       width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized1,
         aes(x = extin_time, y = max_dens, color = as.factor(b), 
             fill = as.factor(b), shape = as.factor(a))) +
    geom_point(size = 2.5, alpha = 0.5) +
    facet_grid(~tau) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    scale_shape_manual(values = 21:25) +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    scale_fill_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("tau") +
    NULL)
  dev.off()
  
  tiff("./run1_statplots/maxtime_extintime.tiff",
       width = 5, height = 5, units = "in", res = 300)
  print(ggplot(data = y_summarized1,
         aes(x = max_time, y = extin_time, color = as.factor(b), 
             fill = as.factor(b), shape = as.factor(a))) +
    geom_point(size = 2.5, alpha = 0.5) +
  #  facet_grid(tau~.) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    scale_shape_manual(values = 21:25) +
    theme_bw() +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    scale_fill_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    NULL)
  dev.off()
  
  tiff("./run1_statplots/maxtime_extintime_facet.tiff",
       width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized1,
         aes(x = max_time, y = extin_time, color = as.factor(b), 
             fill = as.factor(b), shape = as.factor(a))) +
    geom_point(size = 2.5, alpha = 0.5) +
    facet_grid(~tau) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    scale_shape_manual(values = 21:25) +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    scale_fill_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("tau") +
    NULL)
  dev.off()
}

#Making contour plots ----
if (glob_make_statplots) {
  tiff("./run1_statplots/maxtime_contour1.tiff", width = 8, height = 4,
       units = "in", res = 300)
  p1 <- ggplot(data = y_summarized1, 
         aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(b)))) +
    geom_contour_filled(aes(z = max_time)) +
    facet_grid(~tau) +
    scale_fill_viridis_d(direction = -1) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Burst Size") +
    xlab("Infection Rate") +
    labs(fill = "Peak Time", subtitle = "Lysis Time") +
    NULL
  print(p1)
  dev.off()
  
  tiff("./run1_statplots/maxtime_contour2.tiff", width = 8, height = 4,
       units = "in", res = 300)
  p2 <- ggplot(data = y_summarized1, 
         aes(x = as.numeric(as.character(tau)), y = as.numeric(as.character(b)))) +
    geom_contour_filled(aes(z = max_time)) +
    facet_grid(~a) +
    scale_fill_viridis_d(direction = -1) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Burst Size") +
    xlab("Lysis Time") +
    labs(fill = "Peak Time", subtitle = "Infection Rate") +
    NULL
  print(p2)
  dev.off()
  
  tiff("./run1_statplots/maxtime_contour3.tiff", width = 8, height = 4,
       units = "in", res = 300)
  p3 <- ggplot(data = y_summarized1, 
         aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(tau)))) +
    geom_contour_filled(aes(z = max_time)) +
    facet_grid(~b) +
    scale_fill_viridis_d(direction = -1) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Lysis Time") +
    xlab("Infection Rate") +
    labs(fill = "Peak Time", subtitle = "Burst Size") +
    NULL
  print(p3)
  dev.off()
  
  tiff("./run1_statplots/maxtime_contour_all.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                     p2 + theme(legend.position = "none"), 
                     p3 + theme(legend.position = "none"),
                     cowplot::get_legend(p1),
                     #rel_widths = c(1, 1, 1, .4),
                     nrow = 2))
  dev.off()
}

#Make plots that include derivs
ybig_melt <- data.table::melt(as.data.table(ybig1),
                              measure.vars = c("Density", "deriv", "deriv_percap"),
                              variable.name = "var_measured",
                              value.name = "value")
dir.create("run1_dens_and_derivs", showWarnings = FALSE)
if (glob_make_curveplots) {
  for (run in unique(ybig_melt$uniq_run)) {
    tiff(paste("./run1_dens_and_derivs/", run, ".tiff", sep = ""),
         width = 4, height = 8, units = "in", res = 300)
    print(ggplot(data = ybig_melt[ybig_melt$uniq_run == run &
                                    ybig_melt$Pop == "B" &
                                    ybig_melt$time > 0, ],
                 aes(x = time, y = value)) +
            geom_point() +
            facet_grid(var_measured~., scales = "free_y") +
            theme_bw())
    dev.off()
  }
}

# Plot multiple B's on same axes ----
if (glob_make_statplots) {
  tiff("./run1_statplots/B_plots.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  print(ggplot(data = ybig1[ybig1$Pop == "B" &
                       ybig1$Density > 0, ],
         aes(x = time, y = Density+10, color = as.factor(a))) +
    geom_line(lwd = 1, alpha = 0.5) +
    geom_hline(yintercept = 10, lty = 2) +
    facet_grid(tau~b, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme_bw() +
    ggtitle("b (top), tau (side)") +
    NULL)
  dev.off()

  tiff("./run1_statplots/B_plots2.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  print(ggplot(data = ybig1[ybig1$Pop == "B" &
                       ybig1$Density > 0, ],
         aes(x = time, y = Density+10, color = as.factor(tau))) +
    geom_line(lwd = 1, alpha = 0.5) +
    geom_hline(yintercept = 10, lty = 2) +
    facet_grid(b~a, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme_bw() +
    ggtitle("a (top), b (side)") +
    NULL)
  dev.off()

  tiff("./run1_statplots/B_plots3.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  print(ggplot(data = ybig1[ybig1$Pop == "B" &
                       ybig1$Density > 0, ],
         aes(x = time, y = Density+10, color = as.factor(b))) +
    geom_line(lwd = 1, alpha = 0.5) +
    geom_hline(yintercept = 10, lty = 2) +
    facet_grid(tau~a, scales = "free") +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    theme_bw() +
    ggtitle("a (top), tau (side)") +
    NULL)
  dev.off()
}

## Run #2: r, a, b, tau ----
run2 <- run_sims_filewrapper(name = "run2",
                             u_Svals = signif(0.04*10**seq(from = 0, to = -0.7, by = -0.175), 3),
                             k_Svals = c(10**9),
                             avals = 10**seq(from = -12, to = -8, by = 1),
                             tauvals = signif(10**seq(from = 1, to = 2, by = 0.25), 3),
                             bvals = signif(5*10**seq(from = 0, to = 2, by = 0.5), 3),
                             c_SIvals = 1,
                             init_S_dens_vals = 10**6,
                             init_moi_vals = 10**-2,
                             min_dens = 0.1,
                             init_time = 100,
                             init_stepsize = 1,
                             print_info = TRUE,
                             read_file = glob_read_files)

#Check fails/no equils
run2[[2]]

run2[[3]]

#Find peaks & extinction via summarize
ybig2 <- group_by_at(run2[[1]], .vars = 1:17)
ybig2 <- ybig2[complete.cases(ybig2), ]
y_summarized2 <- dplyr::summarize(ybig2,
                           max_dens = max(Density[Pop == "B"]),
                           max_time = time[Pop == "B" & 
                                             Density[Pop == "B"] == max_dens],
                           extin_index = min(which(Pop == "B" &
                                                     Density <= 10**4)),
                           extin_dens = Density[extin_index],
                           extin_time = time[extin_index],
                           extin_time_sincemax = extin_time-max_time,
                           auc = sum(Density[Pop == "B" & time < extin_time])*
                             extin_time,
                           phage_final = max(Density[Pop == "P"]),
                           phage_extin = Density[Pop == "P" & time == extin_time],
                           phage_r = (log(phage_final)-
                                        log(init_S_dens[1]*init_moi[1]))/
                             extin_time,
                           phage_atmaxdens = Density[Pop == "P" & time == max_time]
)

## Plot summarized stats ----
y_sum_melt2 <- reshape2::melt(y_summarized2,
                              id.vars = 1:17,
                              variable.name = "sum_stat",
                              value.name = "stat_val")

dir.create("run2_statplots", showWarnings = FALSE)
if (glob_make_statplots) {
  for (myu_S in unique(y_sum_melt2$u_S)) {
    tiff(paste("./run2_statplots/all_stats_r=", 
               formatC(myu_S, digits = 5, format = "f"), 
               ".tiff", sep = ""),
         width = 5, height = 7, units = "in", res = 300)
    print(ggplot(data = y_sum_melt2[y_sum_melt2$u_S == myu_S &
                                y_sum_melt2$sum_stat %in% 
                                c("max_dens", "max_time", 
                                  "extin_time", 
                                  #"extin_time_sincemax",
                                  "phage_final", 
                                  "phage_r",
                                  "phage_atmaxdens"
                                  ), ],
           aes(x = a, y = stat_val, color = as.factor(b), group = as.factor(b))) +
      geom_point(size = 2, alpha = 0.8) + 
      geom_line(size = 1.1, alpha = 0.6) +
      facet_grid(sum_stat~tau, scales = "free_y") +
      scale_y_continuous(trans = "log10") +
      scale_x_continuous(trans = "log10") +
      scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste("r=", myu_S, " tau", sep = "")) +
      NULL
    )
    dev.off()
  }
}

##Just extin_time ----
if (glob_make_statplots) {
  tiff("./run2_statplots/extin_time_1.tiff",
       width = 5, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized2,
         aes(x = a, y = extin_time, color = as.factor(b), group = as.factor(b))) +
    #geom_point() + 
    geom_line(lwd = 1, alpha = 0.65) +
    scale_color_manual(values = my_cols[c(1, 2, 3, 5, 7)]) +
    facet_grid(u_S~tau) +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL)
  dev.off()

  tiff("./run2_statplots/extin_time_2.tiff",
       width = 5, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized2,
         aes(x = b, y = extin_time, color = as.factor(u_S),
             group = as.factor(u_S))) +
    #geom_point() + 
    geom_line(lwd = 1.25, alpha = 0.8) +
    scale_color_manual(values = my_cols[c(1, 2, 3, 5, 7)]) +
    facet_grid(tau~a) +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL)
  dev.off()

  tiff("./run2_statplots/extin_time_3.tiff",
       width = 5, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized2,
         aes(x = u_S, y = extin_time, color = as.factor(tau),
             group = as.factor(tau))) +
    #geom_point() + 
    geom_line(lwd = 1.25, alpha = 0.8) +
    scale_color_manual(values = my_cols[c(1, 2, 3, 5, 7)]) +
    facet_grid(a~b) +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL)
  dev.off()
  
  tiff("./run2_statplots/extin_time_4.tiff",
       width = 5, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized2,
         aes(x = tau, y = extin_time, color = as.factor(a),
             group = as.factor(a))) +
    #geom_point() + 
    geom_line(lwd = 1.25, alpha = 0.8) +
    scale_color_manual(values = my_cols[c(1, 2, 3, 5, 7)]) +
    facet_grid(b~u_S) +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL)
  dev.off()
}

###Plot stats against ea other ----

##First plot all at once

#Calculate log10's
for (col in c("max_dens", "max_time", "extin_time", "extin_time_sincemax",
              "auc", "phage_final", "phage_r")) {
  y_summarized2[, paste(col, "_log10", sep = "")] <- log10(y_summarized2[, col])
}

#Make plots
if (glob_make_statplots) {
  for (myu_S in unique(y_summarized2$u_S)) {
    tiff(paste("./run2_statplots/stat_cors_r=", 
               formatC(myu_S, digits = 5, format = "f"),
               ".tiff", sep = ""),
         width = 15, height = 15, units = "in", res = 300)
    #Make base figure
    p <- GGally::ggpairs(y_summarized2[y_summarized2$u_S == myu_S, ],
                         aes(color = as.factor(b), shape = as.factor(a)),
                         columns = c("max_dens_log10", "max_time_log10",
                                     "extin_time_log10", 
                                     "extin_time_sincemax_log10",
                                     "auc_log10", 
                                     "phage_final_log10", "phage_r_log10"),
                         lower = list(continuous = "points"),
                         upper = list(continuous = "points")) +
      theme_bw() +
      theme(strip.text = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 0))
    print(p)
    dev.off()
  }
}

##Now selected pairs, looking for underlying functions

#Relating max_dens to max_time
max_dens_func <- function(t, K, P_0, u) {K/(1+((K-P_0)/P_0)*exp(-u*t))}

y_summarized2$pred_maxdens <- max_dens_func(t = y_summarized2$max_time,
                                            K = y_summarized2$k_S,
                                            P_0 = y_summarized2$init_S_dens,
                                            u = y_summarized2$u_S)

if (glob_make_statplots) {
  tiff("./run2_statplots/maxdens_maxtime.tiff",
       width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized2,
         aes(x = max_time, y = max_dens, color = as.factor(a), 
             shape = as.factor(b))) +
    geom_point() +
    facet_grid(~u_S) +
    scale_y_continuous(trans = "log10") +
    geom_line(aes(x = max_time, y = pred_maxdens), color = "black", lty = 3) +
    theme_bw() +
    NULL)
  dev.off()

  # for (myu_S in unique(y_summarized2$u_S)) {
  #   tiff(paste("./run2_statplots/maxdens_maxtime_r=", myu_S, ".tiff", sep = ""),
  #        width = 5, height = 5, units = "in", res = 300)
  #   print(ggplot(data = y_summarized2[y_summarized2$r == myu_S, ],
  #                aes(x = max_time, y = max_dens, color = a, shape = b)) +
  #           geom_point() +
  #           scale_y_continuous(trans = "log10") +
  #           stat_function(fun = max_dens_func,
  #                         args = list(K = 10**9, P_0 = 1*10**6, r = myu_S),
  #                         color = "black", lwd = 1, alpha = 0.1) +
  #           ggtitle(paste("u =", myu_S)) +
  #           theme_bw()
  #   )
  #   dev.off()
  # }
}

#Relating phage_final to max_dens and b

#First try fitting a model to the data
temp <- y_summarized2[y_summarized2$u_S == 0.00798 &
                        y_summarized2$max_dens_log10 < 8.95, ]
model1 <- lm(phage_final_log10 ~ max_dens_log10 + as.factor(b),
             temp)
summary(model1)

if (glob_make_statplots) {
  tiff("./run2_statplots/phagefinal_peakdens_fit.tiff",
       width = 5, height = 5, units = "in", res = 300)
  print(ggplot(data = temp,
               aes(x = max_dens_log10, y = phage_final_log10)) +
          geom_point(aes(color = as.factor(b), shape = as.factor(tau))) +
          geom_line(data = fortify(model1),
                    aes(x = max_dens_log10, y = .fitted, color = `as.factor(b)`)) +
          NULL)
  dev.off()
}

#The model seems to suggest the following "true" underlying model:
y_summarized2$pred_phage_final <- y_summarized2$b*y_summarized2$max_dens

if (glob_make_statplots) {
  tiff("./run2_statplots/phagefinal_peakdens_pred.tiff",
       width = 8, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized2, aes(x = max_dens, y = phage_final, 
                                   color = as.factor(b), shape = as.factor(tau))) +
    geom_point(alpha = 0.5) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") +
    scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
    geom_line(aes(x = max_dens, y = pred_phage_final, color = as.factor(b))) +
    facet_grid(~u_S) +
    NULL)
  dev.off()
}

#Since max_dens is related to max_time, and phage_final is related to max_dens
#We can relate phage_final to max_time
phage_final_func2 <- function(max_time, k, P_0, u, b) {
  #phage_final = b*(K/(1+((K-P_0)/P_0)*exp(-u*max_time)))
  b*(k/(1+((k-P_0)/P_0)*exp(-u*max_time)))
}

y_summarized2$pred_phage_final2 <-
  phage_final_func2(max_time = y_summarized2$max_time,
                    k = y_summarized2$k_S,
                    P_0 = y_summarized2$init_S_dens,
                    u = y_summarized2$u_S,
                    b = y_summarized2$b)

if (glob_make_statplots) {
  tiff("./run2_statplots/phagefinal_maxtime_pred.tiff",
       width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized2,
               aes(x = max_time, y = phage_final, 
                   color = as.factor(b), shape = as.factor(a))) +
          geom_point() +
          geom_line(aes(x = max_time, y = pred_phage_final2,
                        color = as.factor(b), group = as.factor(b))) +
          scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
          facet_wrap(~u_S) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          theme_bw() +
          NULL
  )
  dev.off()
}

##Relating max time and extin time
ggplot(data = y_summarized2[y_summarized2$max_dens < 0.8*10**9, ],
       aes(x = max_time, y = extin_time, 
           color = as.factor(tau))) +
           #color = max_dens/k_S)) +
  geom_point(alpha = 0.5) +
  #facet_wrap(~u_S, scales = "free") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  #geom_abline(slope = 1.162, intercept = 26.6) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")

temp <- lm(extin_time ~ max_time,
   data = y_summarized2[y_summarized2$max_dens < 0.8*10**9, ])

summary(temp)

##Run #3: r, a, b, tau, init_dens, init_moi ----
run3 <- run_sims_filewrapper(name = "run3",
                             u_Svals = c(0.04, 0.0179),
                             k_Svals = c(10**9),
                             avals = 10**seq(from = -12, to = -8, by = 2),
                             tauvals = signif(20**seq(from = 1, to = 1.5, by = 0.5), 3),
                             bvals = signif(5*10**seq(from = 1, to = 2, by = 1), 3),
                             c_SIvals = 1,
                             init_S_dens_vals = c(10**4, 10**5, 10**6),
                             init_moi_vals = c(10**-2, 10**-1, 1),
                             min_dens = 0.1,
                             init_time = 100,
                             init_stepsize = 1,
                             print_info = TRUE,
                             read_file = glob_read_files)

#Check fails/no equils
run3[[2]]

run3[[3]]

#Find peaks & extinction via summarize
ybig3 <- group_by_at(run3[[1]], .vars = 1:17)
ybig3 <- ybig3[complete.cases(ybig3), ]
y_summarized3 <- summarize(ybig3,
                           max_dens = max(Density[Pop == "B"]),
                           max_time = time[Pop == "B" & 
                                             Density[Pop == "B"] == max_dens],
                           extin_index = min(which(Pop == "B" &
                                                     Density <= 10**3)),
                           extin_dens = Density[extin_index],
                           extin_time = time[extin_index],
                           extin_time_sincemax = extin_time-max_time,
                           auc = sum(Density[Pop == "B" & time < extin_time])*
                             extin_time,
                           phage_final = max(Density[Pop == "P"]),
                           phage_extin = Density[Pop == "P" & time == extin_time],
                           phage_r = (log(phage_final)-
                                        log(init_S_dens[1]*init_moi[1]))/
                             extin_time
)

#Make plots of density against time ----
dir.create("run3_dens_curves", showWarnings = FALSE)
if (glob_make_curveplots) {
  dens_offset <- 10
  for (run in unique(ybig3$uniq_run)) {
    tiff(paste("./run3_dens_curves/", run, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(
      ggplot(data = ybig3[ybig3$uniq_run == run &
                            ybig3$Pop %in% c("S", "I", "P"),], 
             aes(x = time, y = Density+dens_offset, color = Pop)) +
        geom_line(lwd = 1.5, alpha = 1) + 
        geom_line(data = ybig3[ybig3$uniq_run == run &
                                 ybig3$Pop == "B",], 
                  aes(x = time, y = Density+dens_offset),
                  color = "black", alpha = 0.5, lwd = 1.1) +
        geom_line(data = ybig3[ybig3$uniq_run == run &
                                 ybig3$Pop == "PI",],
                  aes(x = time, y = Density+dens_offset),
                  color = "black", alpha = 0.5, lwd = 1, lty = 3) +
        geom_point(data = y_summarized3[y_summarized3$uniq_run == run, ],
                   aes(x = max_time, y = max_dens+dens_offset), color = "black") +
        geom_point(data = y_summarized3[y_summarized3$uniq_run == run, ],
                   aes(x = extin_time, y = extin_dens+dens_offset), color = "black") +
        scale_y_continuous(trans = "log10") +
        scale_x_continuous(breaks = seq(from = 0, to = max(ybig3$time), 
                                        by = round(max(ybig3[ybig3$uniq_run == run &
                                                               ybig3$Pop != "B", 
                                                             "time"])/10))) +
        scale_color_manual(values = my_cols[c(2, 3, 1)]) +
        geom_hline(yintercept = 10, lty = 2) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              title = element_text(size = 9)) +
        ggtitle(paste(ybig3[min(which(ybig3$uniq_run == run)), 2:9],
                      collapse = ", ")) +
        labs(y = paste("Density +", dens_offset)) +
        NULL
    )
    dev.off()
  }
}

## Plot summarized stats ----
y_sum_melt3 <- reshape2::melt(y_summarized3,
                              id.vars = 1:17,
                              variable.name = "sum_stat",
                              value.name = "stat_val")

dir.create("./run3_statplots", showWarnings = FALSE)
if (glob_make_statplots) {
  for (myu_S in unique(y_sum_melt3$u_S)) {
    for (myb in unique(y_sum_melt3$b)) {
      tiff(paste("./run3_statplots/all_stats_r=", 
                 formatC(myu_S, digits = 5, format = "f"), 
                 ",b=", myb, ".tiff", sep = ""),
           width = 5, height = 7, units = "in", res = 300)
      print(ggplot(data = y_sum_melt3[y_sum_melt3$u_S == myu_S &
                                        y_sum_melt3$b == myb &
                                        y_sum_melt3$sum_stat %in% 
                                        c("max_dens", "max_time", 
                                          "extin_time", "extin_time_sincemax",
                                          "phage_final", "phage_r"), ],
                   aes(x = init_S_dens, y = stat_val, 
                       color = as.factor(init_moi), group = as.factor(init_moi))) +
              geom_point(size = 2, alpha = 0.8) + 
              geom_line(size = 1.1, alpha = 0.6) +
              facet_grid(sum_stat~tau*a, scales = "free_y") +
              scale_y_continuous(trans = "log10") +
              scale_x_continuous(trans = "log10") +
              scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              ggtitle(paste("r=", myu_S, " tau", sep = "")) +
              NULL
      )
      dev.off()
    }
  }

  #Let's focus in on extin_time
  tiff("./run3_statplots/extin_time.tiff",
       width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized3,
         aes(x = init_S_dens, y = extin_time, color = as.factor(init_moi),
             shape = as.factor(a))) +
    geom_point() +
    facet_grid(u_S~b*tau) +
  #  geom_line() +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    NULL)
  dev.off()
}

###Plot stats against ea other ----

##First plot all at once

#Calculate log10's
for (col in c("max_dens", "max_time", "extin_time", "extin_time_sincemax",
              "auc", "phage_final", "phage_r")) {
  y_summarized3[, paste(col, "_log10", sep = "")] <- log10(y_summarized3[, col])
}

#Make plots
if (glob_make_statplots) {
  for (myu_S in unique(y_summarized3$u_S)) {
    tiff(paste("./run3_statplots/stat_cors_r=", 
               formatC(myu_S, digits = 5, format = "f"),
               ".tiff", sep = ""),
         width = 15, height = 15, units = "in", res = 300)
    #Make base figure
    p <- GGally::ggpairs(y_summarized3[y_summarized3$u_S == myu_S, ],
                         aes(color = as.factor(init_moi), 
                             shape = as.factor(init_S_dens)),
                         columns = c("max_dens_log10", "max_time_log10",
                                     "extin_time_log10", 
                                     "extin_time_sincemax_log10",
                                     "auc_log10", 
                                     "phage_final_log10", "phage_r_log10"),
                         lower = list(continuous = "points"),
                         upper = list(continuous = "points")) +
      theme_bw() +
      theme(strip.text = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 0))
    print(p)
    dev.off()
  }
}

###Run #4: r, a, b, tau ----
# run4 <- run_sims_filewrapper(name = "run4",
#                              u_Svals = signif(0.04*10**seq(from = 1, to = -2, by = -0.67), 3),
#                              k_Svals = c(10**9),
#                              avals = 10**seq(from = -14, to = -6, by = 2),
#                              tauvals = signif(10**seq(from = 0, to = 3, by = 0.75), 3),
#                              bvals = signif(5*10**seq(from = -1, to = 3, by = 1), 3),
#                              c_SIvals = 1,
#                              init_S_dens_vals = 10**6,
#                              init_moi_vals = 10**-2,
#                              min_dens = 0.1,
#                              init_time = 100,
#                              init_stepsize = 1,
#                              print_info = TRUE,
#                              read_file = glob_read_files)
#                              
# #Check fails/no equils
# run4[[2]]
# 
# run4[[3]]
# 
# #Find peaks & extinction via summarize
# ybig4 <- group_by_at(run4[[1]], .vars = 1:17)
# ybig4 <- ybig4[complete.cases(ybig4), ]
# y_summarized4 <- summarize(ybig4,
#                            max_dens = max(Density[Pop == "B"]),
#                            max_time = time[Pop == "B" & 
#                                              Density[Pop == "B"] == max_dens],
#                            extin_index = min(which(Pop == "B" &
#                                                      Density <= 10**4)),
#                            extin_dens = Density[extin_index],
#                            extin_time = time[extin_index],
#                            extin_time_sincemax = extin_time-max_time,
#                            auc = sum(Density[Pop == "B" & time < extin_time])*
#                              extin_time,
#                            phage_final = max(Density[Pop == "P"]),
#                            phage_extin = Density[Pop == "P" & time == extin_time],
#                            phage_r = (log(phage_final)-
#                                         log(init_S_dens[1]*init_moi[1]))/
#                              extin_time
# )
# 
# ## Plot summarized stats ----
# y_sum_melt4 <- reshape2::melt(y_summarized4,
#                               id.vars = 1:17,
#                               variable.name = "sum_stat",
#                               value.name = "stat_val")
# 
# if (glob_make_statplots) {
#   for (myu_S in unique(y_sum_melt4$u_S)) {
#     tiff(paste("./run2_statplots/all_stats_r=", 
#                formatC(myu_S, digits = 5, format = "f"), 
#                ".tiff", sep = ""),
#          width = 5, height = 7, units = "in", res = 300)
#     print(ggplot(data = y_sum_melt4[y_sum_melt4$r == myu_S &
#                                       y_sum_melt4$sum_stat %in% 
#                                       c("max_dens", "max_time", 
#                                         "extin_time", "extin_time_sincemax",
#                                         "phage_final", "phage_r"), ],
#                  aes(x = a, y = stat_val, 
#                      color = as.factor(b), group = as.factor(b))) +
#             geom_point(size = 2, alpha = 0.8) + 
#             geom_line(size = 1.1, alpha = 0.6) +
#             facet_grid(sum_stat~tau, scales = "free_y") +
#             scale_y_continuous(trans = "log10") +
#             scale_x_continuous(trans = "log10") +
#             scale_color_manual(values = colorRampPalette(colors = c("gold", "dark red"))(5)) +
#             theme_bw() +
#             theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#             ggtitle(paste("r=", myu_S, " tau", sep = "")) +
#             NULL
#     )
#     dev.off()
#   }
# }

###Run #5: init_dens and init_moi as r,k,a,b,tau indiv ----
run5_params <- data.frame(matrix(NA, ncol = 7, nrow = 25*5*3))
colnames(run5_params) <- c("init_S_dens_vals", "init_moi_vals",
                           "u_Svals", "k_Svals", "avals", "bvals", "tauvals")

run5_params$init_S_dens_vals <- rep(rep(10**seq(from = 4, to = 6, by = 0.5),
                                      times = 5),
                                  times = 5*3)
run5_params$init_moi_vals <- rep(rep(10**seq(from = -3, to = -1, by = 0.5), 
                                each = 5),
                            times = 5*3)
run5_params$u_Svals <- c(rep(0.0179, times = 25*0*3),
                  rep(signif(0.04*10**seq(from = 0, to = -0.7, by = -0.35), 3),
                      each = 25),
                  rep(0.0179, times = 25*4*3))
run5_params$k_Svals <- c(rep(10**9, times = 25*1*3),
                   rep(10**c(8, 9, 10), each = 25),
                   rep(10**9, times = 25*3*3))
run5_params$avals <- c(rep(10**-10, times = 25*2*3),
                   rep(10**seq(from = -12, to = -8, by = 2), each = 25),
                   rep(10**-10, times = 25*2*3))
run5_params$bvals <- c(rep(63.2, times = 25*3*3),
                   rep(signif(2*10**seq(from = 1, to = 2, by = 0.5), 3), each = 25),
                   rep(63.2, times = 25*1*3))
run5_params$tauvals <- c(rep(37.6, times = 25*4*3),
                   rep(signif(10**seq(from = 1.25, to = 2, by = 0.325), 3), each = 25),
                   rep(37.6, times = 25*0*3))
             
run5 <- run_sims_filewrapper(name = "run5",
                             u_Svals = run5_params$u_Svals,
                             k_Svals = run5_params$k_Svals,
                             avals = run5_params$avals,
                             tauvals = run5_params$tauvals,
                             bvals = run5_params$bvals,
                             init_S_dens_vals = run5_params$init_S_dens_vals,
                             init_moi_vals = run5_params$init_moi_vals,
                             c_SIvals = 1,
                             combinatorial = FALSE,
                             min_dens = 0.1,
                             init_time = 100,
                             init_stepsize = 1,
                             print_info = TRUE,
                             read_file = glob_read_files)

run5[[2]]

run5[[3]]

ybig5 <- group_by_at(run5[[1]], .vars = 1:17)
y_summarized5 <- summarize(ybig5,
                           max_dens = max(Density[Pop == "B"]),
                           max_time = time[Pop == "B" & 
                                             Density[Pop == "B"] == max_dens],
                           extin_index = min(which(Pop == "B" &
                                                     Density <= 9.99*10**3)),
                           extin_dens = Density[extin_index],
                           extin_time = time[extin_index],
                           extin_time_sincemax = extin_time-max_time,
                           auc = sum(Density[Pop == "B" & time < extin_time])*
                             extin_time,
                           phage_final = max(Density[Pop == "P"]),
                           phage_extin = Density[Pop == "P" & time == extin_time],
                           phage_r = (log(phage_final)-
                                        log(init_S_dens[1]*init_moi[1]))/
                             extin_time
)

dir.create("./run5_statplots", showWarnings = FALSE)
if (glob_make_statplots) {
  tiff("./run5_statplots/max_time1.tiff", width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized5[y_summarized5$k_S == 10**9 &
                                      y_summarized5$a == 10**-10 &
                                      y_summarized5$b == 63.2 &
                                      y_summarized5$tau == 37.6, ],
               aes(x = log10(init_S_dens), 
                   y = log10(init_moi),
                   z = max_time)) +
          geom_contour_filled() +
          facet_grid(~u_S) +
          ggtitle("r") +
          scale_fill_viridis_d(direction = -1) +
          theme_bw() +
          NULL)
  dev.off()
  
  tiff("./run5_statplots/max_time2.tiff", width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized5[y_summarized5$u_S == 0.0179 &
                                      y_summarized5$a == 10**-10 &
                                      y_summarized5$b == 63.2 &
                                      y_summarized5$tau == 37.6, ],
               aes(x = log10(init_S_dens), 
                   y = log10(init_moi),
                   z = max_time)) +
          geom_contour_filled() +
          facet_grid(~k_S) +
          ggtitle("k") +
          scale_fill_viridis_d(direction = -1) +
          theme_bw() +
          NULL)
  dev.off()
  
  tiff("./run5_statplots/max_time3.tiff", width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized5[y_summarized5$u_S == 0.0179 &
                                      y_summarized5$k_S == 10**9 &
                                      y_summarized5$b == 63.2 &
                                      y_summarized5$tau == 37.6, ],
               aes(x = log10(init_S_dens), 
                   y = log10(init_moi),
                   z = max_time)) +
          geom_contour_filled() +
          facet_grid(~a) +
          ggtitle("a") +
          scale_fill_viridis_d(direction = -1) +
          theme_bw() +
          NULL)
  dev.off()
  
  tiff("./run5_statplots/max_time4.tiff", width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized5[y_summarized5$u_S == 0.0179 &
                                      y_summarized5$k_S == 10**9 &
                                      y_summarized5$a == 10**-10 &
                                      y_summarized5$tau == 37.6, ],
               aes(x = log10(init_S_dens), 
                   y = log10(init_moi),
                   z = max_time)) +
          geom_contour_filled() +
          facet_grid(~b) +
          ggtitle("b") +
          scale_fill_viridis_d(direction = -1) +
          theme_bw() +
          NULL)
  dev.off()
  
  tiff("./run5_statplots/max_time5.tiff", width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized5[y_summarized5$u_S == 0.0179 &
                                      y_summarized5$k_S == 10**9 &
                                      y_summarized5$a == 10**-10 &
                                      y_summarized5$b == 63.2, ],
               aes(x = log10(init_S_dens), 
                   y = log10(init_moi),
                   z = max_time)) +
          geom_contour_filled() +
          facet_grid(~tau) +
          ggtitle("tau") +
          scale_fill_viridis_d(direction = -1) +
          theme_bw() +
          NULL)
  dev.off()
  
  tiff("./run5_statplots/maxdens_phagefinal.tiff", width = 6, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized5[y_summarized5$u_S == 0.0179 &
                                      y_summarized5$k_S == 10**9 &
                                      y_summarized5$a == 10**-10 &
                                      y_summarized5$tau == 37.6, ],
               aes(x = max_dens, y = phage_final, color = as.factor(b),
                   shape = as.factor(init_moi))) +
          geom_point() +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          NULL)
  dev.off()
}


##Run 6: a, u, k (bacterial traits) ----
run6 <- 
  run_sims_filewrapper(name = "run6",
                       u_Svals = signif(0.04*10**seq(from = 0, to = -0.7, by = -0.175), 3),
                       k_Svals = signif(10**c(8, 8.5, 9, 9.5, 10), 3),
                       avals = 10**seq(from = -12, to = -8, by = 1),
                       tauvals = 50,
                       bvals = 50,
                       init_S_dens_vals = 10**6,
                       init_moi_vals = 10**-2,
                       c_SIvals = 1,
                       combinatorial = TRUE,
                       min_dens = 0.1,
                       init_time = 100,
                       init_stepsize = 1,
                       print_info = TRUE,
                       read_file = glob_read_files)

ybig6 <- group_by_at(run6[[1]], .vars = 1:17)
y_summarized6 <- summarize(ybig6,
                           max_dens = max(Density[Pop == "B"]),
                           max_time = time[Pop == "B" & 
                                             Density[Pop == "B"] == max_dens],
                           extin_index = min(which(Pop == "B" &
                                                     Density <= 9.99*10**3)),
                           extin_dens = Density[extin_index],
                           extin_time = time[extin_index],
                           extin_time_sincemax = extin_time-max_time,
                           auc = sum(Density[Pop == "B" & time < extin_time])*
                             extin_time,
                           phage_final = max(Density[Pop == "P"]),
                           phage_extin = Density[Pop == "P" & time == extin_time],
                           phage_r = (log(phage_final)-
                                        log(init_S_dens[1]*init_moi[1]))/
                             extin_time
)

dir.create("run6_statplots", showWarnings = FALSE)
if (glob_make_statplots) {
  #Making contour plots
  tiff("./run6_statplots/maxtime_contour1.tiff", width = 5, height = 5,
       units = "in", res = 300)
  p1 <- ggplot(data = y_summarized6, 
               aes(x = as.numeric(as.character(u_S)), y = as.numeric(as.character(k_S)))) +
    geom_contour_filled(aes(z = max_time)) +
    facet_grid(~a) +
    scale_fill_viridis_d(direction = -1) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Carrying Capacity") +
    xlab("Growth Rate") +
    labs(fill = "Peak Time", subtitle = "Infection Rate") +
    NULL
  print(p1)
  dev.off()

  tiff("./run6_statplots/maxtime_contour2.tiff", width = 5, height = 5,
       units = "in", res = 300)
  p2 <- ggplot(data = y_summarized6, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(k_S)))) +
    geom_contour_filled(aes(z = max_time)) +
    facet_grid(~u_S) +
    scale_fill_viridis_d(direction = -1) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Carrying Capacity") +
    xlab("Infection Rate") +
    labs(fill = "Peak Time", subtitle = "Growth Rate") +
    NULL
  print(p2)
  dev.off()

  tiff("./run6_statplots/maxtime_contour3.tiff", width = 5, height = 5,
       units = "in", res = 300)
  p3 <- ggplot(data = y_summarized6, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(u_S)))) +
    geom_contour_filled(aes(z = max_time)) +
    facet_grid(~k_S) +
    scale_fill_viridis_d(direction = -1) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Growth Rate") +
    xlab("Infection Rate") +
    labs(fill = "Peak Time", subtitle = "Carrying Capacity") +
    NULL
  print(p3)
  dev.off()

  tiff("./run6_statplots/maxtime_contour_all.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                     p2 + theme(legend.position = "none"), 
                     p3 + theme(legend.position = "none"),
                     cowplot::get_legend(p1),
                     #rel_widths = c(1, 1, 1, .4),
                     nrow = 2))
  dev.off()
}

##Run 7: r, k, a, b, tau +/- coinfection ----

###TODO: check all the errors & see what's the cause
###       also figure out why the predict function is making two lines????

run7 <- 
  run_sims_filewrapper(name = "run7",
                       u_Svals = signif(0.04*10**c(0, -0.7), 3),
                       k_Svals = signif(10**c(8, 10), 3),
                       avals = 10**seq(from = -12, to = -8, by = 2),
                       tauvals = signif(10**seq(from = 1, to = 2, by = 0.5), 3),
                       bvals = signif(5*10**seq(from = 0, to = 2, by = 1), 3),
                       init_S_dens_vals = 10**6,
                       init_moi_vals = 10**-2,
                       c_SIvals = 1,
                       zvals = c(0, 1),
                       combinatorial = TRUE,
                       min_dens = 0.1,
                       init_time = 100,
                       init_stepsize = 1,
                       print_info = TRUE,
                       read_file = glob_read_files)

#Find peaks & extinction via summarize
ybig7 <- group_by_at(run7[[1]][!is.na(run7[[1]]$uniq_run), ], .vars = 1:17)
y_summarized7 <- summarize(ybig7,
                           equil = max(equil),
                           max_dens = max(Density[Pop == "B"]),
                           max_time = time[Pop == "B" & 
                                             Density[Pop == "B"] == max_dens],
                           extin_index = min(which(Pop == "B" &
                                                     Density <= 10**4)),
                           extin_dens = Density[extin_index],
                           extin_time = time[extin_index],
                           auc = sum(Density[Pop == "B" & time < extin_time])*
                             extin_time,
                           phage_final = max(Density[Pop == "P"]),
                           phage_extin = max(Density[Pop == "P" & time <= extin_time]),
                           phage_r = (log(phage_final)-
                                        log(init_S_dens[1]*init_moi[1]))/
                             extin_time,
                           run_time = max(time)
)

#Make plots of density against time ----
dir.create("run7_dens_curves", showWarnings = FALSE)
if (glob_make_curveplots) {
  dens_offset <- 10
  for (run in unique(ybig7$uniq_run)) {
    tiff(paste("./run7_dens_curves/", run, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(
      ggplot(data = ybig7[ybig7$uniq_run == run &
                            ybig7$Pop %in% c("S", "I", "P", "R"),], 
             aes(x = time, y = Density+dens_offset, color = Pop)) +
        geom_line(lwd = 1.5, alpha = 1) + 
        geom_line(data = ybig7[ybig7$uniq_run == run &
                                 ybig7$Pop == "B",], 
                  aes(x = time, y = Density+dens_offset),
                  color = "black", alpha = 0.5, lwd = 1.1) +
        geom_line(data = ybig7[ybig7$uniq_run == run &
                                 ybig7$Pop == "PI",],
                  aes(x = time, y = Density+dens_offset),
                  color = "black", alpha = 0.5, lwd = 1, lty = 3) +
        geom_point(data = y_summarized7[y_summarized7$uniq_run == run, ],
                   aes(x = max_time, y = max_dens+dens_offset), color = "black") +
        geom_point(data = y_summarized7[y_summarized7$uniq_run == run, ],
                   aes(x = extin_time, y = extin_dens+dens_offset), color = "black") +
        scale_y_continuous(trans = "log10") +
        # scale_x_continuous(breaks = seq(from = 0, to = max(ybig7$time), 
        #                                 by = round(max(ybig7[ybig7$uniq_run == run &
        #                                                        ybig7$Pop != "B", 
        #                                                      "time"])/10))) +
        scale_color_manual(limits = c("S", "I", "P", "R"),
                           values = my_cols[c(2, 3, 1, 7)]) +
        geom_hline(yintercept = 10, lty = 2) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              title = element_text(size = 9)) +
        ggtitle(paste(ybig7[min(which(ybig7$uniq_run == run)), 6:8],
                      collapse = ", ")) +
        labs(y = paste("Density +", dens_offset)) +
        NULL
    )
    dev.off()
  }
}

dir.create("run7_statplots", showWarnings = FALSE)

#Relating max_dens to max_time
max_dens_func <- function(t, K, P_0, u) {K/(1+((K-P_0)/P_0)*exp(-u*t))}

y_summarized7$pred_maxdens <- max_dens_func(t = y_summarized7$max_time,
                                            K = y_summarized7$k_S,
                                            P_0 = y_summarized7$init_S_dens,
                                            u = y_summarized7$u_S)
y_summarized7 <- y_summarized7[order(y_summarized7$pred_maxdens),]

if (glob_make_statplots) {
  tiff("./run7_statplots/maxdens_maxtime.tiff",
       width = 8, height = 4, units = "in", res = 300)
  print(ggplot(data = y_summarized7[y_summarized7$equil == TRUE &
                                complete.cases(y_summarized7), ],
         aes(x = max_time, y = max_dens, color = as.factor(a), shape = as.factor(b),
             group = k_S*z*u_S)) +
    geom_point() +
    facet_grid(k_S ~ z*u_S, scales = "free") +
    scale_y_continuous(trans = "log10") +
    geom_line(aes(x = max_time, y = pred_maxdens), color = "black", lty = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL)
  dev.off()
}

##Run 8: r, k, a, b, tau +/- (costless) resistance ----
run8_params <- expand.grid(u_vals = signif(0.04*10**c(0, -0.7), 3),
                           k_vals = signif(10**c(8, 10), 3),
                           avals = 10**seq(from = -12, to = -8, by = 4),
                           tauvals = signif(10**seq(from = 1, to = 2, by = 1), 3),
                           bvals = signif(5*10**seq(from = 0, to = 2, by = 2), 3),
                           mvals = signif(10**c(-4, -5.5, -7), 3))
                          
run8 <-
  run_sims_filewrapper(name = "run8",
                       u_Svals = run8_params$u_vals,
                       u_Rvals = run8_params$u_vals,
                       k_Svals = run8_params$k_vals,
                       k_Rvals = run8_params$k_vals,
                       avals = run8_params$avals,
                       tauvals = run8_params$tauvals,
                       bvals = run8_params$bvals,
                       init_S_dens_vals = 10**6,
                       init_moi_vals = 10**-2,
                       c_SIvals = 1,
                       zvals = 0,
                       mvals = run8_params$mvals,
                       combinatorial = FALSE,
                       min_dens = 0.1,
                       init_time = 100,
                       init_stepsize = 1,
                       print_info = TRUE,
                       read_file = glob_read_files)

run8[[2]]
run8[[3]]

#Find peaks & extinction via summarize
ybig8 <- group_by_at(run8[[1]], .vars = 1:17)
y_summarized8 <- summarize(ybig8,
                           max_index = find_local_extrema(values = Density[Pop == "B"],
                                                        return_minima = FALSE,
                                                        width_limit = 3)[1],
                           max_dens = Density[Pop == "B"][max_index],
                           max_time = time[Pop == "B"][max_index],
                           extin_index = find_local_extrema(values = Density[Pop == "B"],
                                                            return_maxima = FALSE,
                                                            width_limit = 3)[1],
                           extin_dens = Density[Pop == "B"][extin_index],
                           extin_time = time[Pop == "B"][extin_index],
                           auc = sum(Density[Pop == "B" & time < extin_time])*
                             extin_time,
                           phage_final = max(Density[Pop == "P"]),
                           phage_extin = Density[Pop == "P" & time == extin_time],
                           phage_r = (log(phage_final)-
                                        log(init_S_dens[1]*init_moi[1]))/
                             extin_time,
                           run_time = max(time)
)

#Make plots of density against time ----
dir.create("run8_dens_curves", showWarnings = FALSE)
if (glob_make_curveplots) {
  dens_offset <- 10
  for (run in unique(ybig8$uniq_run)) {
    tiff(paste("./run8_dens_curves/", run, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(
      ggplot(data = ybig8[ybig8$uniq_run == run &
                            ybig8$Pop %in% c("S", "I", "P", "R"),], 
             aes(x = time, y = Density+dens_offset, color = Pop)) +
        geom_line(lwd = 1.5, alpha = 1) + 
        geom_line(data = ybig8[ybig8$uniq_run == run &
                                 ybig8$Pop == "B",], 
                  aes(x = time, y = Density+dens_offset),
                  color = "black", alpha = 0.5, lwd = 1.1) +
        geom_line(data = ybig8[ybig8$uniq_run == run &
                                 ybig8$Pop == "PI",],
                  aes(x = time, y = Density+dens_offset),
                  color = "black", alpha = 0.5, lwd = 1, lty = 3) +
        geom_point(data = y_summarized8[y_summarized8$uniq_run == run, ],
                   aes(x = max_time, y = max_dens+dens_offset), color = "black") +
        geom_point(data = y_summarized8[y_summarized8$uniq_run == run, ],
                   aes(x = extin_time, y = extin_dens+dens_offset), color = "black") +
        scale_y_continuous(trans = "log10") +
        scale_x_continuous(breaks = seq(from = 0, to = max(ybig8$time), 
                                        by = round(max(ybig8[ybig8$uniq_run == run &
                                                               ybig8$Pop != "B", 
                                                             "time"])/10))) +
        scale_color_manual(limits = c("S", "I", "P", "R"),
                           values = my_cols[c(2, 3, 1, 7)]) +
        geom_hline(yintercept = 10, lty = 2) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              title = element_text(size = 9)) +
        ggtitle(paste(ybig8[min(which(ybig8$uniq_run == run)), 6:8],
                      collapse = ", ")) +
        labs(y = paste("Density +", dens_offset)) +
        NULL
    )
    dev.off()
  }
}

##Run 9: init dens & moi across a few b vals ----



## Run #10: a, b, tau +/- coinfection (phage traits) ----
run10 <- run_sims_filewrapper(name = "run10",
                             u_Svals = c(0.023), #(30 min doubling time)
                             k_Svals = c(10**9),
                             avals = 10**seq(from = -12, to = -8, by = 1),
                             tauvals = signif(10**seq(from = 1, to = 2, by = 0.25), 3),
                             bvals = signif(5*10**seq(from = 0, to = 2, by = 0.5), 3),
                             c_SIvals = 1,
                             init_S_dens_vals = 10**6,
                             init_moi_vals = 10**-2,
                             zvals = c(0, 1),
                             min_dens = 0.1,
                             init_time = 100,
                             init_stepsize = 1,
                             print_info = TRUE,
                             read_file = glob_read_files)

#Find peaks & extinction via summarize
ybig10 <- group_by_at(run10[[1]], .vars = 1:17)
y_summarized10 <- summarize(ybig10,
                           max_dens = max(Density[Pop == "B"]),
                           max_time = time[Pop == "B" & 
                                             Density[Pop == "B"] == max_dens],
                           extin_index = min(which(Pop == "B" &
                                                     Density <= 10**4)),
                           extin_dens = Density[extin_index],
                           extin_time = time[extin_index],
                           auc = sum(Density[Pop == "B" & time < extin_time])*
                             extin_time,
                           phage_final = max(Density[Pop == "P"]),
                           phage_extin = Density[Pop == "P" & time == extin_time],
                           phage_r = (log(phage_final)-
                                        log(init_S_dens[1]*init_moi[1]))/
                             extin_time,
                           run_time = max(time),
                           near_k = if(max_dens >= 0.95*k_S[1]) {1} else{NA}
)

dir.create("./run10_statplots", showWarnings = FALSE)
if (glob_make_statplots) {
  ##Make contour plots of peak time
  tiff("./run10_statplots/maxtime_contour1.tiff", width = 8, height = 4,
       units = "in", res = 300)
  p1 <- ggplot(data = y_summarized10, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(b)))) +
    geom_contour_filled(aes(z = max_time)) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~tau) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Burst Size") +
    xlab("Infection Rate") +
    labs(fill = "Peak Time", subtitle = "Lysis Time") +
    guides(fill = guide_legend(ncol = 2)) +
    NULL
  print(p1)
  dev.off()
  
  tiff("./run10_statplots/maxtime_contour2.tiff", width = 8, height = 4,
       units = "in", res = 300)
  p2 <- ggplot(data = y_summarized10, 
               aes(x = as.numeric(as.character(tau)), y = as.numeric(as.character(b)))) +
    geom_contour_filled(aes(z = max_time)) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~a) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Burst Size") +
    xlab("Lysis Time") +
    labs(fill = "Peak Time", subtitle = "Infection Rate") +
    NULL
  print(p2)
  dev.off()
  
  tiff("./run10_statplots/maxtime_contour3.tiff", width = 8, height = 4,
       units = "in", res = 300)
  p3 <- ggplot(data = y_summarized10, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(tau)))) +
    geom_contour_filled(aes(z = max_time)) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~b) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Lysis Time") +
    xlab("Infection Rate") +
    labs(fill = "Peak Time", subtitle = "Burst Size") +
    NULL
  print(p3)
  dev.off()
  
  tiff("./run10_statplots/maxtime_contour_all.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                           p2 + theme(legend.position = "none"), 
                           p3 + theme(legend.position = "none"),
                           cowplot::get_legend(p1),
                           #rel_widths = c(1, 1, 1, .4),
                           nrow = 2))
  dev.off()
  
  ##Make contour plots of extinction time
  p1 <- ggplot(data = y_summarized10, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(b)))) +
    geom_contour_filled(aes(z = log10(extin_time))) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~tau) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Burst Size") +
    xlab("Infection Rate") +
    labs(fill = "log10(Extinction Time)", subtitle = "Lysis Time") +
    guides(fill = guide_legend(ncol = 2)) +
    NULL
  p2 <- ggplot(data = y_summarized10, 
               aes(x = as.numeric(as.character(tau)), y = as.numeric(as.character(b)))) +
    geom_contour_filled(aes(z = log10(extin_time))) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~a) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Burst Size") +
    xlab("Lysis Time") +
    labs(fill = "Extinction Time", subtitle = "Infection Rate") +
    NULL
  p3 <- ggplot(data = y_summarized10, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(tau)))) +
    geom_contour_filled(aes(z = log10(extin_time))) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~b) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Lysis Time") +
    xlab("Infection Rate") +
    labs(fill = "Extinction Time", subtitle = "Burst Size") +
    NULL
  
  tiff("./run10_statplots/extintime_contour_all.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                           p2 + theme(legend.position = "none"), 
                           p3 + theme(legend.position = "none"),
                           cowplot::get_legend(p1),
                           #rel_widths = c(1, 1, 1, .4),
                           nrow = 2))
  dev.off()
}

##Run 11: a, u, k +/- coinfection (bacterial traits) ----
run11 <- 
  run_sims_filewrapper(name = "run11",
                       u_Svals = signif(0.04*10**seq(from = 0, to = -0.7, by = -0.175), 3),
                       k_Svals = signif(10**c(8, 8.5, 9, 9.5, 10), 3),
                       avals = 10**seq(from = -12, to = -8, by = 1),
                       tauvals = 50,
                       bvals = 50,
                       init_S_dens_vals = 10**6,
                       init_moi_vals = 10**-2,
                       c_SIvals = 1,
                       zvals = c(0, 1),
                       combinatorial = TRUE,
                       min_dens = 0.1,
                       init_time = 100,
                       init_stepsize = 1,
                       print_info = TRUE,
                       read_file = glob_read_files)

ybig11 <- group_by_at(run11[[1]], .vars = 1:17)
y_summarized11 <- summarize(ybig11,
                           max_dens = max(Density[Pop == "B"]),
                           max_time = time[Pop == "B" & 
                                             Density[Pop == "B"] == max_dens],
                           extin_index = min(which(Pop == "B" &
                                                     Density <= 9.99*10**3)),
                           extin_dens = Density[extin_index],
                           extin_time = time[extin_index],
                           extin_time_sincemax = extin_time-max_time,
                           auc = sum(Density[Pop == "B" & time < extin_time])*
                             extin_time,
                           phage_final = max(Density[Pop == "P"]),
                           phage_extin = Density[Pop == "P" & time == extin_time],
                           phage_r = (log(phage_final)-
                                        log(init_S_dens[1]*init_moi[1]))/
                             extin_time,
                           near_k = if(max_dens >= 0.95*k_S[1]) {1} else{NA}
)

dir.create("run11_statplots", showWarnings = FALSE)
if (glob_make_statplots) {
  #Making contour plots for peak time
  tiff("./run11_statplots/maxtime_contour1.tiff", width = 5, height = 5,
       units = "in", res = 300)
  p1 <- ggplot(data = y_summarized11, 
               aes(x = as.numeric(as.character(u_S)), y = as.numeric(as.character(k_S)))) +
    geom_contour_filled(aes(z = max_time)) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~a) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Carrying Capacity") +
    xlab("Growth Rate") +
    labs(fill = "Peak Time", subtitle = "Infection Rate") +
    guides(fill = guide_legend(ncol = 2)) +
    NULL
  print(p1)
  dev.off()
  
  tiff("./run11_statplots/maxtime_contour2.tiff", width = 5, height = 5,
       units = "in", res = 300)
  p2 <- ggplot(data = y_summarized11, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(k_S)))) +
    geom_contour_filled(aes(z = max_time)) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~u_S) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Carrying Capacity") +
    xlab("Infection Rate") +
    labs(fill = "Peak Time", subtitle = "Growth Rate") +
    NULL
  print(p2)
  dev.off()
  
  tiff("./run11_statplots/maxtime_contour3.tiff", width = 5, height = 5,
       units = "in", res = 300)
  p3 <- ggplot(data = y_summarized11, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(u_S)))) +
    geom_contour_filled(aes(z = max_time)) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~k_S) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Growth Rate") +
    xlab("Infection Rate") +
    labs(fill = "Peak Time", subtitle = "Carrying Capacity") +
    NULL
  print(p3)
  dev.off()
  
  tiff("./run11_statplots/maxtime_contour_all.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                           p2 + theme(legend.position = "none"), 
                           p3 + theme(legend.position = "none"),
                           cowplot::get_legend(p1),
                           #rel_widths = c(1, 1, 1, .4),
                           nrow = 2))
  dev.off()
  
  ##Make extinction time contour plots
  p1 <- ggplot(data = y_summarized11, 
               aes(x = as.numeric(as.character(u_S)), y = as.numeric(as.character(k_S)))) +
    geom_contour_filled(aes(z = log10(extin_time))) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~a) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Carrying Capacity") +
    xlab("Growth Rate") +
    labs(fill = "Extinction Time", subtitle = "Infection Rate") +
    guides(fill = guide_legend(ncol = 2)) +
    NULL
  
  p2 <- ggplot(data = y_summarized11, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(k_S)))) +
    geom_contour_filled(aes(z = log10(extin_time))) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~u_S) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Carrying Capacity") +
    xlab("Infection Rate") +
    labs(fill = "Extinction Time", subtitle = "Growth Rate") +
    NULL
  
  p3 <- ggplot(data = y_summarized11, 
               aes(x = as.numeric(as.character(a)), y = as.numeric(as.character(u_S)))) +
    geom_contour_filled(aes(z = log10(extin_time))) +
    scale_fill_viridis_d(direction = -1) +
    geom_point(aes(color = as.factor(near_k))) +
    scale_color_manual(values = c("gray")) +
    facet_grid(z~k_S) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Growth Rate") +
    xlab("Infection Rate") +
    labs(fill = "Extinction Time", subtitle = "Carrying Capacity") +
    NULL
  
  tiff("./run11_statplots/extintime_contour_all.tiff", width = 8, height = 6,
       units = "in", res = 300)
  print(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                           p2 + theme(legend.position = "none"), 
                           p3 + theme(legend.position = "none"),
                           cowplot::get_legend(p1),
                           #rel_widths = c(1, 1, 1, .4),
                           nrow = 2))
  dev.off()
}



###Work in progress below: ----

##Relating max dens to extin_time

#TODO: explore other sigmoid functions
# (e.g. geeralized logistic function)
#Note that, because it is logistic-looking in log-space the underlying shape
# is never going to be logistic.
#Could also look into something where the percap growth rate
# grows then decays with distance to midpoint or something like that

# #Define squared percent errors function
# sq_err_func <- function(params, x_vals, y_vals) {
#   #first input is a vector containing all the parameters
#   # 1 - r
#   # 2 - K
#   # 3 - P_0
#   #second input is the vector of x values
#   r <- params["r"]
#   K <- params["K"]
#   P_0 <- params["P_0"]
#   pred_vals <- K/(1+((K-P_0)/P_0)*exp(-r*x_vals))
#  # return(sum((100*(pred_vals-y_vals))**2))
#   return(sum((100*(log10(pred_vals)-log10(y_vals)))**2))
#   #return(sum((100*(pred_vals-y_vals)/y_vals)**2))
# }
# 
# #this stuff is producing warnings even when not run:
# # Warning messages:
# #   1: Unknown or uninitialised column: 'par'
# 
# # #Find fits
# # fit_output <- data.frame(myr = numeric(), r = numeric(), 
# #                          K = numeric(), P_0 = numeric())
# # # for (myr in unique(y_summarized2$r)) {
# #   myr <- 0.00798
# #   temp <- optim(fn = sq_err_func,
# #                 #par = c(r = myr, K = 10**9, P_0 = 10**6), 
# #                 x_vals = y_summarized2$extin_time[y_summarized2$r == myr],
# #                 y_vals = y_summarized2$max_dens[y_summarized2$r == myr],
# #                 method = "BFGS")
# #   fit_output <- rbind(fit_output,
# #                       data.frame(myr = myr,
# #                                  r = temp$par["r"],
# #                                  K = temp$par["K"],
# #                                  P_0 = temp$par["P_0"]))
# # # }
# 
# #Calc ratio of logistic curve's fit r to the r in the simulations (myr)
# fit_output$r_divby_myr <- fit_output$r/fit_output$myr
# 
# #Relating max_dens to extin_time
# max_dens_func2 <- function(extin_time, K, P_0, r) {
#   log10(K/(1+((K-P_0)/P_0)*exp(-r*extin_time)))
# }
#   
# #for (myr in unique(y_summarized2$r)) {
#   myr <- 0.00798
#   myrow <- which(fit_output$myr == myr)
#   print(ggplot(data = y_summarized2[y_summarized2$r == myr, ],
#                aes(x = extin_time, y = max_dens, color = b, shape = a)) +
#           geom_point() +
#           # stat_function(fun = max_dens_func2,
#           #               args = list(K = fit_output$K[myrow], 
#           #                           P_0 = fit_output$P_0[myrow], 
#           #                           r = fit_output$r[myrow]),
#           #               color = "black", alpha = 0.1, lwd = 1) +
#           # stat_function(fun = max_dens_func2,
#           #               args = list(K = 10**9, 
#           #                           P_0 = 7*10**5, 
#           #                           r = .007),
#           #               color = "black", alpha = 0.1, lwd = 1) +
#           scale_y_continuous(trans = "log10") +
#           scale_x_continuous(trans = "log10") +
#           ggtitle(paste("r =", myr)) +
#           theme_bw() +
#           NULL
#         )
# #}
# 
# ##Relating phage_r to extin_time
# y_summarized2$tau <- as.factor(y_summarized2$tau)
# 
# phage_r_model1 <- lm(phage_r_log10 ~ extin_time_log10 + 
#                        a + a:extin_time_log10 + 
#                        b + b:extin_time_log10 +
#                        tau + tau:extin_time_log10,
#                      y_summarized2)
# anova(phage_r_model1)
# summary(phage_r_model1)
# 
#   #for (myr in unique(y_summarized2$r)) {
#   myr <- 0.00798
#   print(ggplot(data = y_summarized2[y_summarized2$max_dens_log10 < 8.95, ],
#                aes(x = extin_time_log10, y = phage_r_log10, 
#                    color = a, shape = tau)) +
#           geom_point(alpha = 0.5) +
# #          facet_grid(~b) +
#           # stat_function(fun = max_dens_func2,
#           #               args = list(K = 10**9, P_0 = 10**6, r = myr)) +
# #          scale_y_continuous(trans = "log10") +
# #          ggtitle(paste("r =", myr)) +
#           NULL)
#   #}
# 
# for (myr in unique(y_summarized2$r)) {
#   #myr <- 0.00798
#   tiff(paste("./run2_statplots/phager_extintime_r=", myr, ".tiff", sep = ""),
#        width = 6, height = 5, units = "in", res = 300)
#   print(ggplot(data = y_summarized2[y_summarized2$r == myr, ],
#                aes(x = extin_time, y = phage_r, 
#                    color = a, shape = tau)) +
#           geom_point(alpha = 0.8, size = 2) +
#           theme_bw() +
#           scale_y_continuous(trans = "log10") +
#           scale_x_continuous(trans = "log10") +
#           ggtitle(paste("r =", myr)) +
#           NULL
#   )
#   dev.off()
# }
# 
# 
# ##Testing whether B decay after max_time can be fit with a
# # logistic-like curve
# temp1 <- ybig1[ybig1$uniq_run == 4, ]
# temp2 <- y_summarized1[y_summarized1$uniq_run == 4, ]
# 
# func1 <- function(t, max_dens, max_time, r) {
#   #max_dens - exp(r*(t-max_time))
#   max_dens - (max_dens/(1+exp(-r*(t-max_time))))
# }
# 
# func1_log <- function(t, max_dens, max_time, r) {
#   log10(max_dens - (max_dens/(1+exp(-r*(t-max_time)))))
# }
# 
# func1_optim <- function(params, times, density) {
#   max_dens = params[["max_dens"]]
#   max_time = params[["max_time"]]
#   r = params[["r"]]
#   fit_data <- func1(t = times, max_dens = max_dens,
#                     max_time = max_time, r = r)
#   return(sum((density-fit_data)**2))
# }
# 
# fit1 <- optim(par = list(max_dens = temp2$max_dens,
#                         max_time = temp2$max_time,
#                         r = 0.15),
#              fn = func1_optim,
#              times = temp1$time[temp1$time >= temp2$max_time],
#              density = temp1$Density[temp1$time >= temp2$max_time])
#                
# 
# ggplot(data = temp1[temp1$Pop == "B", ],
#        aes(x = time, y = Density+10)) +
#   geom_line(aes(color = Pop, group = Pop),
#             lwd = 1, alpha = 1) +
#   geom_hline(yintercept = 10, lty = 2) +
#   stat_function(mapping = aes(x = time),
#                 fun = func1_log,
#                 args = list(max_dens = 59450913826, #temp2$max_dens,
#                             max_time = 19810832468, #temp2$max_time+60,
#                             r = 14138914293), # 0.15),
#                 color = "black") +
#   scale_y_continuous(trans = "log10", limits = c(10, NA)) +
#   NULL


