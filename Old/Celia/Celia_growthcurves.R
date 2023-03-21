# ONLY-BACTERIA curve
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
  dY <- c(S = 0, I = 0, P = 0)
  
  ##Calculate dS
  
  #V3 (logistic dS/dt) (including competition from I pop)
  #dS/dt = rS((K-S-c*I)/K) - aSP
  dY["S"] <- parms["r"] * y["S"] * 
    ((parms["K"] - y["S"] - parms["c"] * y["I"])/parms["K"]) - 
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
  #dP/dt = baS(t-tau)P(t-tau) - aSP
  if (t < parms["tau"]) {
    dY["P"] <- -parms["a"] * y["S"] * y["P"]
  } else {
    dY["P"] <- parms["b"] * parms["a"] * 
      lagvalue(t-parms["tau"], 1)*lagvalue(t-parms["tau"], 3) - 
      parms["a"]*y["S"]*y["P"]
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

yinit <- c(S = 10**6,
           I = 0,
           P = 0)
params <- c(r = 0.04, 
            a = 10**-10, 
            b = 50, 
            tau = 10,
            K = 10**9,
            c = 1,
            warnings = 0, 
            thresh_min_dens = 10**-100)
times <- seq(from = 0, to = 500, by = 1)
yout <- as.data.frame(
  dede(y = yinit, times = times, func = derivs, parms = params))

##Plot results ----
head(yout)
tidyr::pivot_longer(yout, c(S, I, P), names_to = "Population", values_to = "Density")
# Fisrt, I tried this command above, but R said that Population and Density didn't exist
# when trying to plot. That's why I tried what's below and it worked.

library(tidyr)
yout_plot <- pivot_longer(yout, c(S, I, P), names_to = "Population", values_to = "Density + 10")


##We've reshaped the data. Now, we'll plot the density of the three populations over time
tiff(paste("./Celia/", sep = ""),
     width = 4, height = 4, units = "in", res = 200)
print(ggplot(data = yout_plot, aes(x = time, y = Density + 10, color = Population)) +
  geom_line(lwd = 1.5, alpha = 1/2) +
  scale_color_manual(values = c("#56B4E9", "#009E73", "#E69F00", "#000000")) +
    theme_bw() +
  scale_y_continuous(trans = "log10")) +
xlab("Time")
dev.off()

#Okabe and Ito 2008 colorblind-safe qualitative color scale ----
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

## Import libraries ----
library(deSolve)
library(tidyr)
library(ggplot2)
library(dplyr)
# library(plyr)
# library(tidyverse)
# library(caret)

## Define derivatives function ----
derivs <- function(t, y, parms) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Issue warning about too small/negative yvals (if warnings is 1)
  if (parms["warnings"]==1 & any(y < parms["thresh_min_dens"])) {
    warning(paste("population(s)", paste(which(y < parms["thresh_min_dens"]),
                                  collapse = ","), 
                  "below thresh_min_dens, treating as 0"))
  }
  
  #Set small/negative y values to 0 so they don't affect the dN's
  y[y < parms["thresh_min_dens"]] <- 0
  
  #Create output vector
  dY <- c(Susceptible = 0, Infected = 0, Phage = 0)
  
  ##Calculate dS
  
  #V3 (logistic dS/dt) (including competition from I population)
  #dS/dt = rS((K-S-c*I)/K) - aSP
  dY["Susceptible"] <- parms["r"] * y["Susceptible"] *
    ((parms["K"] - y["Susceptible"] - parms["c"] * y["Infected"])/parms["K"]) -
    parms["a"] * y["Susceptible"] * y["Phage"]
  
  ##Calculate dI
  #dI/dt = aSP - aS(t-tau)P(t-tau)
  if (t < parms["tau"]) {
    dY["Infected"] <- parms["a"] * y["Susceptible"] * y["Phage"]
  } else {
    dY["Infected"] <- parms["a"] * y["Susceptible"]*y["Phage"] -
      parms["a"] * lagvalue(t - parms["tau"], 1)*lagvalue(t - parms["tau"], 3)
  }
  
  ##Calculate dP
  #dP/dt = baS(t-tau)P(t-tau) - aSP
  if (t < parms["tau"]) {
    dY["Phage"] <- -parms["a"] * y["Susceptible"] * y["Phage"]
  } else {
    dY["Phage"] <- parms["b"] * parms["a"] *
      lagvalue(t-parms["tau"], 1)*lagvalue(t-parms["tau"], 3) -
      parms["a"]*y["Susceptible"]*y["Phage"]
  }
  
  #Issue warning about too large population (if warnings is TRUE)
  if (parms["warnings"]==1 & any(y > 10**100)) {
    warning(paste("population(s)",paste(which(y > 10**100), collapse = ","), 
                  "exceed max limit, 10^100, returning dY = 0"))
  }
  dY[y > 10**100] <- 0
  
  #From documentation: The return value of func should be a list, whose first 
  #element is a vector containing the derivatives of y with respect to time
  return(list(dY))
}

#Realistic parameter values ----
#  
#   r ranges from 0.04 (a 17-minute doubling time) to 0.007 (a 90-minute doubling time).
#   K ranges from 10^6 to 10^10, although typically we focus on 10^8 to 10^9.
#   a ranges from 10^-12 to 10^-8.
#   tau ranges from 10 to 105 (mins).
#   b ranges from 5 to 1000.

## Define function for running simulations across many parameter values ----
run_sims <- function(rvals,
                     kvals,
                     avals,
                     tauvals,
                     bvals,
                     cvals = 1,
                     init_bact_dens_vals = 10**6,
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
  
  #Save parameter values provided into dataframe
  # taking all different combinations
  if (combinatorial) {
    param_combos <- expand.grid(list("rvals" = rvals, "kvals" = kvals, 
                                     "avals" = avals, "tauvals" = tauvals, 
                                     "bvals" = bvals, "cvals" = cvals, 
                                     "init_bact_dens_vals" = init_bact_dens_vals, 
                                     "init_moi_vals" = init_moi_vals),
                                stringsAsFactors = FALSE)
    num_sims <- nrow(param_combos)
  } else { #not combinatorial
    num_sims <- max(sapply(X = list(rvals, kvals, avals, tauvals, bvals,
                                    cvals, init_bact_dens_vals, init_moi_vals), 
                           FUN = length))
    
    #Check for parameter lengths being non-divisible with the
    # number of simulations inferred from the longest parameter length
    if (!all(num_sims %% sapply(X = list(rvals, kvals, avals, tauvals, bvals,
                                         cvals, init_bact_dens_vals, init_moi_vals), 
                                FUN = length) == 0)) {
      warning("Combinatorial is true but longest param vals length is not a multiple of all other param vals lengths")
    }
    
    #Save parameters into dataframe, replicating shorter parameter
    # vectors as needed to reach # of simulations
    param_combos <- data.frame("rvals" = rep_len(rvals, num_sims), 
                               "kvals" = rep_len(kvals, num_sims), 
                               "avals" = rep_len(avals, num_sims), 
                               "tauvals" = rep_len(tauvals, num_sims), 
                               "bvals" = rep_len(bvals, num_sims), 
                               "cvals" = rep_len(cvals, num_sims), 
                               "init_bact_dens_vals" = rep_len(init_bact_dens_vals, num_sims), 
                               "init_moi_vals" = rep_len(init_moi_vals, num_sims),
                               stringsAsFactors = FALSE)
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
  ybig <- data.frame("uniq_run" = rep(NA, 5*(1+init_time/init_stepsize)*num_sims),
                     "r" = NA, "a" = NA, "b" = NA, "tau" = NA, "K" = NA,
                     "c" = NA, "init_bact_dens" = NA, "init_moi" = NA,
                     "equil" = NA, "Time" = NA, "Population" = as.character(NA), 
                     "Density" = NA, stringsAsFactors = F)
  
  for (i in 1:nrow(param_combos)) { #i acts as the uniq_run counter
    myr <- param_combos$rvals[i]
    myk <- param_combos$kvals[i]
    mya <- param_combos$avals[i]
    mytau <- param_combos$tauvals[i]
    myb <- param_combos$bvals[i]
    myc <- param_combos$cvals[i]
    my_init_bact <- param_combos$init_bact_dens_vals[i]
    my_moi <- param_combos$init_moi_vals[i]
    
    #Define pops & parameters
    yinit <- c(Susceptible = my_init_bact,
               Infected = 0,
               Phage = my_init_bact*my_moi)
    params <- c(r = myr, a = mya, b = myb, tau = mytau,
                K = myk, c = myc,
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
        yout_list[[2]] <- yout_list[[2]][!(is.nan(yout_list[[2]]$Susceptible) |
                                             is.nan(yout_list[[2]]$Infected) |
                                             is.nan(yout_list[[2]]$Phage)), ]
        
        #S and I both at equil, we're done
        if (yout_list[[2]]$Susceptible[nrow(yout_list[[2]])] < min_dens & 
            yout_list[[2]]$Infected[nrow(yout_list[[2]])] < min_dens) {
          keep_running <- FALSE
          at_equil <- TRUE
          #S not at equil, need more time
        } else if (yout_list[[2]]$Susceptible[nrow(yout_list[[2]])] >= min_dens) { 
          j <- j+1
          #I not at equil (but S is because above check failed),
          #   first we'll lengthen the simulation
          #    (to make sure it was long enough to catch the last burst)
          #   then we'll start shrinking our step size
        } else if (yout_list[[2]]$Infected[nrow(yout_list[[2]])] >= min_dens) {
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
      yout_list[[2]]$All <- yout_list[[2]]$Susceptible + yout_list[[2]]$Infected
      #Calculate all phage (PI)
      yout_list[[2]]$PhageInfected <- yout_list[[2]]$Phage + yout_list[[2]]$Infected
      
      #Reshape, add parameters, and fill into ybig in right rows
      ybig[((i-1)*5*(1+init_time/init_stepsize)+1):
             ((i)*5*(1+init_time/init_stepsize)), ] <- 
        cbind(data.frame(uniq_run = i, r = myr, a = mya, 
                         b = myb, tau = mytau, K = myk, 
                         c = myc, init_bact_dens = my_init_bact, 
                         init_moi = my_moi, equil = at_equil),
              data.table::melt(data = data.table::as.data.table(yout_list[[2]]), 
                               id.vars = c("time"),
                               value.name = "Density", 
                               variable.name = "Population",
                               variable.factor = FALSE))
      
      #If the run failed
    } else {
      if (is.null(yfail)) { #This is the first failed run
        yfail <- data.frame(uniq_run = i, r = myr, a = mya, 
                            b = myb, tau = mytau, K = myk,
                            c = myc, init_bact_dens = my_init_bact, 
                            init_moi = my_moi, equil = at_equil)
      } else { #This is a non-first failed run
        yfail <- rbind(yfail, 
                       data.frame(uniq_run = i, r = myr, a = mya, 
                                  b = myb, tau = mytau, K = myk, 
                                  c = myc, init_bact_dens = my_init_bact, 
                                  init_moi = my_moi, equil = at_equil))
      }
    }
    
    #Print progress update
    if (print_info & i %in% progress_seq) {
      print(paste((which(progress_seq == i)-1)*10,
                  "% completed", sep = ""))
    }
    
  }
  
  #Pull out all the runs that didn't reach equilibrium
  y_noequil <- NULL
  for (run in unique(ybig$uniq_run[which(!ybig$equil)])) {
    if (is.null(y_noequil)) {
      y_noequil <- ybig[min(which(ybig$uniq_run == run)), 1:9]
    } else {
      y_noequil <- rbind(y_noequil, ybig[min(which(ybig$uniq_run == run)), 1:9])
    }
  }
  
  return(list(ybig, y_noequil, yfail))
  
  #Code for visualizing while debugging
  # if (F) {
  #   #Code for plotting population sizes over Time
  #   ymelt <- reshape2::melt(data = as.data.frame(yout_list[[2]]), 
  #                           id = c("Time"),
  #                           value.name = "Density", 
  #                           variable.name = "Population")
  #   
  #   ggplot(data = ymelt, 
  #          aes(x = Time, y = Density+10, color = Population)) +
  #     geom_line(lwd = 1.5, alpha = 1) + 
  #     scale_y_continuous(trans = "log10") +
  #     scale_x_continuous(breaks = seq(from = 0, to = max(ymelt$Time), 
  #                                     by = round(max(ymelt$Time)/10))) +
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

#Initial run that varied one parameter at a time (bigFINAL) ---- 
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sim_bigFINAL_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sim_bigFINAL <- list(NULL, NULL, NULL)
  sim_bigFINAL[[1]] <- read.csv("./Celia/sim_bigFINAL_1.csv", stringsAsFactors = F)
  if ("sim_bigFINAL_2.csv" %in% list.files("./Celia/")) {
    sim_bigFINAL[[2]] <- read.csv("./Celia/sim_bigFINAL_2.csv", stringsAsFactors = F)
  }
  if ("sim_bigFINAL_3.csv" %in% list.files("./Celia/")) {
    sim_bigFINAL[[3]] <- read.csv("./Celia/sim_bigFINAL_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  bigFINAL_params <- as.data.frame(matrix(
    #         r      k     a       tau   b  c  init_S  moi
    data = c(0.04, 10**9, 10**-10, 10, 50, 1, 10**6, 0.01,
             0.04, 10**8, 10**-10, 10, 50, 1, 10**6, 0.01,
             0.03, 10**9, 10**-10, 10, 50, 1, 10**6, 0.01,
             0.007, 10**9, 10**-10, 10, 50, 1, 10**6, 0.01,
             0.04, 10**9, 10**-13, 10, 50, 1, 10**6, 0.01,
             0.04, 10**9, 10**-8, 10, 50, 1, 10**6, 0.01,
             0.04, 10**9, 10**-10, 65, 50, 1, 10**6, 0.01,
             0.04, 10**9, 10**-10, 120, 50, 1, 10**6, 0.01,
             0.04, 10**9, 10**-10, 10, 20, 1, 10**6, 0.01,
             0.04, 10**9, 10**-10, 10, 500, 1, 10**6, 0.01,
             0.04, 10**9, 10**-10, 10, 850, 1, 10**6, 0.01), 
    ncol = 8, byrow = TRUE))
  colnames(bigFINAL_params) <- c("r", "K", "a", "tau", "b", "c", 
                                 "init_bact_dens", "init_moi")
  sim_bigFINAL <- run_sims(rvals = bigFINAL_params$r,
                           kvals = bigFINAL_params$K,
                           avals = bigFINAL_params$a,
                           tauvals = bigFINAL_params$tau,
                           bvals = bigFINAL_params$b,
                           cvals = bigFINAL_params$c,
                           init_bact_dens_vals = bigFINAL_params$init_bact_dens,
                           init_moi_vals = bigFINAL_params$init_moi,
                           init_time = 200, init_stepsize = 0.5,
                           combinatorial = FALSE)
  #Save results so they can be re-loaded in future
  write.csv(sim_bigFINAL[[1]], "./Celia/sim_bigFINAL_1.csv", row.names = F)
  if (!is.null(sim_bigFINAL[[2]])) {write.csv(sim_bigFINAL[[2]], "./Celia/sim_bigFINAL_2.csv", row.names = F)}
  if (!is.null(sim_bigFINAL[[3]])) {write.csv(sim_bigFINAL[[3]], "./Celia/sim_bigFINAL_3.csv", row.names = F)}
}

bigFINAL_plot <- sim_bigFINAL[[1]]
class(bigFINAL_plot)

#This pivot_wider is a stopgap because the original simulations were
# not run using run_sims and were analyzed in the wider format
bigFINAL <- pivot_wider(bigFINAL_plot, names_from = Population, values_from = Density)

#Make plots
for (my_run in unique(bigFINAL_plot$uniq_run)) {
  dir.create("./Celia/bigFINAL_plots/", showWarnings = FALSE)
  tiff(paste("./Celia/bigFINAL_plots/", my_run, ".tiff", sep = ""),
             width = 4, height = 4, units = "in", res = 200)
  print(ggplot(data = bigFINAL_plot[bigFINAL_plot$uniq_run == my_run &
                                      bigFINAL_plot$Population != "PI", ],
               aes(x = time, y = Density + 10, color = Population)) +
          geom_line(lwd = 1.5, alpha = 1/2) +
          scale_color_manual(values = my_cols[c(8, 2, 3, 1)]) +
          scale_y_continuous(trans = "log10") +
          ggtitle(paste("Run #", my_run, sep = "")))
  dev.off()
}

##Summarize bigFINAL (max density & initial slope) ----

## Summarize to find the maximum of each of the simulations included
## in the big data frame
bigFINAL <- dplyr::group_by(bigFINAL, b, tau, a, r, K, c)
FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), 
                          maxtime = time[B == maximum_B])

## Finding the slope with summarize 
FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), 
                          maxtime = time[B == maximum_B],
                          extin_time = time[min(which(time > maxtime & B < 10**4))],
                          slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                       time[time < maxtime & B < 0.1*K])$coefficients[2],
                          intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                           time[time < maxtime & B < 0.1*K])$coefficients[1])
FINAL

## Here, we try to analyze and plot row by row of the data frame bigFINAL
for (row in 1:nrow(FINAL)) {
  bigfinal_rows <- which(FINAL$b[row] == bigFINAL$b & 
                           FINAL$tau[row] == bigFINAL$tau &
                           FINAL$a[row] == bigFINAL$a &
                           FINAL$r[row] == bigFINAL$r &
                           FINAL$K[row] == bigFINAL$K &
                           FINAL$c[row] == bigFINAL$c)
  dir.create("./Celia/bigFINAL_slopeplots/", showWarnings = FALSE)
  tiff(paste("./Celia/bigFINAL_slopeplots/", row, ".tiff", sep = ""),
       width = 4, height = 4, units = "in", res = 200)
  print(ggplot(data = bigFINAL[bigfinal_rows, ],
               aes(x = time, y = B)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = FINAL$slope[row], intercept = FINAL$intercept[row],
                      color = "red") +
          geom_point(data = FINAL[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          geom_point(data = FINAL[row, ], aes(x = extin_time, y = 10**4),
                     col = "green", size = 3) +
          NULL
  )
  dev.off()
}


 
## Let's run sims1 ----
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sims1_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sims1 <- list(NULL, NULL, NULL)
  sims1[[1]] <- read.csv("./Celia/sims1_1.csv", stringsAsFactors = F)
  if ("sims1_2.csv" %in% list.files("./Celia/")) {
    sims1[[2]] <- read.csv("./Celia/sims1_2.csv", stringsAsFactors = F)
  }
  if ("sims1_3.csv" %in% list.files("./Celia/")) {
    sims1[[3]] <- read.csv("./Celia/sims1_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  sims1 <- run_sims(bvals = c(75, 480, 760), avals = c(10**-10.5, 10**-9, 10**-8.25),
                    kvals = c(10**7.75, 10**8.8), tauvals = c(33, 57, 95),
                    rvals = c(0.009, 0.018, 0.025))
  #Save results so they can be re-loaded in future
  write.csv(sims1[[1]], "./Celia/sims1_1.csv", row.names = F)
  if (!is.null(sims1[[2]])) {write.csv(sims1[[2]], "./Celia/sims1_2.csv", row.names = F)}
  if (!is.null(sims1[[3]])) {write.csv(sims1[[3]], "./Celia/sims1_3.csv", row.names = F)}
}

length(sims1)
# When we use [] means that the output will be a list, while when we use [[]], means
# that it'll be a dataframe (explanation with trains and their cars).
class(sims1[[1]])
sims1
table(sims1[[1]]$Population)

# By chencking group 2 and 3 we see if there has been any error: simulations that 
# didn'r reach the equilibrium (2), and simulations that failed (3).
sims1[[2]]
sims1[[3]]

## Notice that the summarize will be slightly different since the data has been
## pivot_longer'd. So, we have to use the density column and make a subset where 
## "Population" = B
sims1_plot <- sims1[[1]]
class(sims1_plot)
#This pivot_wider is a stopgap because the original simulations were
# not run using run_sims and were analyzed in the wider format
sims1 <- pivot_wider(sims1_plot, names_from = Population, values_from = Density)

#Rearrange the legend in the order we want
sims1_plot$Population <- factor(sims1_plot$Population, levels = c("Susceptible",
                                                                  "Infected",
                                                                  "Phage",
                                                                  "All"))
#Make plots
for (my_run in unique(sims1_plot$uniq_run)) {
  dir.create("./Celia/sims1_plots/", showWarnings = FALSE)
  tiff(paste("./Celia/sims1_plots/", my_run, ".tiff", sep = ""),
       width = 4, height = 4, units = "in", res = 200)
  print(ggplot(data = sims1_plot[sims1_plot$uniq_run == my_run &
                                      sims1_plot$Population != "PhageInfected", ],
               aes(x = Time, y = Density + 10, color = Population)) +
          geom_line(lwd = 1.5, alpha = 1/2) +
          scale_color_manual(values = my_cols[c(2, 3, 1, 8)]) +
          scale_y_continuous(trans = "log10") +
          theme_bw())
  dev.off()
}

## Now, we want to find the maximum_B, the maxtime, the extintion time, and the
## slope of each simulation
group_sims1 <- dplyr::group_by(sims1, uniq_run, b, tau, a, r, K, c)
group_sims1
class(group_sims1)

sum_sims1 <- dplyr::summarise(group_sims1, maximum_B = max(All),                
                              maxtime = Time[All == maximum_B],
                              extin_time = Time[min(which(Time > maxtime & All < 10**4))],
                              slope = lm(log10(All[Time < maxtime & All < 0.1*K]) ~ 
                                           Time[Time < maxtime & All < 0.1*K])$coefficients[2],
                              intercept = lm(log10(All[Time < maxtime & All < 0.1*K]) ~ 
                                               Time[Time < maxtime & All < 0.1*K])$coefficients[1])
View(sum_sims1)

for (row in 1:nrow(sum_sims1)) {
  bigfinal_rows <- which(sum_sims1$b[row] == group_sims1$b & 
                           sum_sims1$tau[row] == group_sims1$tau &
                           sum_sims1$a[row] == group_sims1$a &
                           sum_sims1$r[row] == group_sims1$r &
                           sum_sims1$K[row] == group_sims1$K &
                           sum_sims1$c[row] == group_sims1$c)
  dir.create("./Celia/sims1_slopeplots/", showWarnings = FALSE)
  tiff(paste("./Celia/Sims1_slopeplots/", sum_sims1[row, "uniq_run"], ".tiff", sep = ""),
       width = 6, height = 6, units = "cm", res = 150)
  print(ggplot(data = group_sims1[bigfinal_rows, ],
               aes(x = Time, y = All)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims1$slope[row], intercept = sum_sims1$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims1[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          geom_point(data = sum_sims1[row, ], aes(x = extin_time, y = 10**4),
                     cil = "green", size = 3) +
          NULL
  )
  
  dev.off()
}

## Now, we want to make a gglpot that represents all the simulations with the
## summarized data
# We start by analyzing how the paramters affcet maxtime
ggplot(data = sum_sims1, aes(x = log10(b), y = maxtime, color = as.factor(K),
                       shape = as.factor(a))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r) +
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = "lm")
## lm function on our data

# How to see how tau and b affect maxtime as if they were INDEPENDENT from 
# each other?
regi_sims1 <- lm(maxtime ~ tau + b + a + r, data = sum_sims1)
summary(regi_sims1)

# Calulate the RSE
sigma(regi_sims1)/mean(sum_sims1$maxtime) # The lower the RSE, the more accurate the model.
# The RSE estimate gives a measure of error of prediction.

# Analyzing the data as if the parameters were CORRELATED
reg_sims1 <- lm(maxtime ~ tau*b*a*r, data = sum_sims1)
summary(reg_sims1)

# What about maximum_B?
ggplot(data = sum_sims1, aes(x = log10(b), y = maximum_B, color = as.factor(K),
                             shape = as.factor(a))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r) +
  scale_y_continuous(trans = "log10")

# PLotting maxtime and maximum_B with r and K
ggplot(data = sum_sims1, aes(x = maxtime, y = maximum_B,
                             )) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(K ~ r) +
  scale_y_continuous(trans = "log10")



## Let's run sims2 ----
## Where we'll have c, K, and r constant, and we'll be changing between 5
## different values of a, b, and tau
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sims2_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sims2 <- list(NULL, NULL, NULL)
  sims2[[1]] <- read.csv("./Celia/sims2_1.csv", stringsAsFactors = F)
  if ("sims2_2.csv" %in% list.files("./Celia/")) {
    sims2[[2]] <- read.csv("./Celia/sims2_2.csv", stringsAsFactors = F)
  }
  if ("sims2_3.csv" %in% list.files("./Celia/")) {
    sims2[[3]] <- read.csv("./Celia/sims2_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  sims2 <- run_sims(bvals = c(50, 100, 200, 400, 800),
                    avals = c(10**-12, 10**-11, 10**-10, 10**-9, 10**-8),
                    kvals = c(10**9), rvals = c(0.04),
                    tauvals = c(15, 22.5, 33.75, 50.625, 75.9375))
  #Save results so they can be re-loaded in future
  write.csv(sims2[[1]], "./Celia/sims2_1.csv", row.names = F)
  if (!is.null(sims2[[2]])) {write.csv(sims2[[2]], "./Celia/sims2_2.csv", row.names = F)}
  if (!is.null(sims2[[3]])) {write.csv(sims2[[3]], "./Celia/sims2_3.csv", row.names = F)}
}

length(sims2)                  
sims2[[1]]
sims2[[2]]
sims2[[3]]

## Now that we're sure that everything went well, we'll strat summarizing the data
sims2_plot <- sims2[[1]]
class(sims2_plot)
#This pivot_wider is a stopgap because the original simulations were
# not run using run_sims and were analyzed in the wider format
sims2 <- pivot_wider(sims2_plot, names_from = Pop, values_from = Density)

#Make plots
for (my_run in unique(sims2_plot$uniq_run)) {
  dir.create("./Celia/sims2_plots/", showWarnings = FALSE)
  tiff(paste("./Celia/sims2_plots/", my_run, ".tiff", sep = ""),
       width = 4, height = 4, units = "in", res = 200)
  print(ggplot(data = sims2_plot[sims2_plot$uniq_run == my_run &
                                   sims2_plot$Pop != "PI", ],
               aes(x = time, y = Density + 10, color = Pop)) +
          geom_line(lwd = 1.5, alpha = 1/2) +
          scale_color_manual(values = my_cols[c(8, 2, 3, 1)]) +
          scale_y_continuous(trans = "log10") +
          ggtitle(paste("Run #", my_run, sep = "")))
  dev.off()
}


group_sims2 <- dplyr::group_by(sims2, uniq_run, a, b, c, K, tau, r)

sum_sims2 <- dplyr::summarise(group_sims2, maximum_B = max(B),                
                              maxtime = time[B == maximum_B],
                              extin_time = time[min(which(time > maxtime & B < 10**4))],
                              slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                           time[time < maxtime & B < 0.1*K])$coefficients[2],
                              intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                               time[time < maxtime & B < 0.1*K])$coefficients[1])
View(sum_sims2)

for (row in 1:nrow(sum_sims2)) {
  bigfinal_rows <- which(sum_sims2$b[row] == group_sims2$b & 
                           sum_sims2$tau[row] == group_sims2$tau &
                           sum_sims2$a[row] == group_sims2$a &
                           sum_sims2$r[row] == group_sims2$r &
                           sum_sims2$K[row] == group_sims2$K &
                           sum_sims2$c[row] == group_sims2$c)
  dir.create("./Celia/sims2_slopeplots/", showWarnings = FALSE)
  tiff(paste("./Celia/Sims2_slopeplots/", sum_sims2[row, "uniq_run"], ".tiff", sep = ""),
       width = 6, height = 6, units = "cm", res = 150)
  print(ggplot(data = group_sims2[bigfinal_rows, ],
               aes(x = time, y = B)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims2$slope[row], intercept = sum_sims2$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims2[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          geom_point(data = sum_sims2[row, ], aes(x = extin_time, y = 10**4),
                     col = "green", size = 3) +
          NULL
  )
  
  dev.off()
}

## Let's make the ggplot for this data
ggplot(data = sum_sims2, aes(x = tau, y = extin_time, color = as.factor(b))) +
  geom_point(size = 2.5, alpha = 1) +
  facet_grid(. ~ a) +
  theme(legend.title = element_text("Burst Size")) +
  xlab("Lysis Time") +
  ylab("Extinction Time") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(values = my_cols) +
  NULL

# How to see how the parameters affect maxtime as if they were INDEPENDENT from 
# each other?
reg_sims2 <- lm(extin_time ~ tau + b + a, data = sum_sims2)
summary(reg_sims2)

# Calulate the RSE
sigma(reg_sims2)/mean(sum_sims2$maxtime) # The lower the RSE, the more accurate the model.
# The RSE estimate gives a measure of error of prediction.

# Analyzing the data as if the parameters were CORRELATED
# INTERACTION EFFECTS
# Build the model
regi_sims2 <- lm(extin_time ~ tau*b*a, data = sum_sims2)
                   
summary(regi_sims2)

## We want to plot the multiple linear regressions
# Read data set
sum_sims2
# Create multiple linear regressions
lm_fit <- lm(maxtime ~ log10(b)*log10(a)*log10(tau), data = sum_sims2)
summary(lm_fit)
# Save predictions of the model in the new data frame together with the variable
# you want to plot against
sum_sims2_predicted <- data.frame(maxtime_pred = predict(lm_fit, sum_sims2),
                                  tau = sum_sims2$tau,
                                  b = sum_sims2$b,
                                  a = sum_sims2$a)
sum_sims2_predicted

# This is the predicted line of multiple linear regressions
ggplot(data = sum_sims2, aes(x = log10(tau), y = maxtime, 
                             color = as.factor(b))) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_line(data = sum_sims2_predicted, aes(x = log10(tau), y = maxtime_pred, 
                                            color = as.factor(b))) +
  facet_grid(. ~ a) +
  theme(legend.title = element_text("Burst Size")) +
  xlab("Lysis Time") +
  ylab("Maximum Time") +
  theme_bw()

#Plot maxtime using contours

ggplot(data = sum_sims2, aes(y = log10(b), 
                             x = log10(tau), 
                             z = -maxtime)) +
  geom_contour_filled() +
  facet_grid(~a) +
  labs(title = "maxtime", subtitle = "a") +
  #scale_fill_brewer(type = "div", palette = 1)
  NULL

ggplot(data = sum_sims2, aes(y = log10(b), 
                             x = log10(tau), 
                             z = -extin_time)) +
  geom_contour_filled() +
  facet_grid(~a) +
  labs(title = "extin time", subtitle = "a") +
  #scale_fill_brewer(type = "div", palette = 1)
  NULL

#1
ggplot(data = sum_sims2, aes(y = b, 
                             x = tau)) +
  geom_contour_filled(aes(z = maxtime)) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(a~.) +
  theme_bw() +
  labs(title = "Peak Time", subtitle = "Infection Rate") +
  ylab("Burst Size") +
  xlab("Lysis Time") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  NULL
  
#2
ggplot(data = sum_sims2, aes(x = tau, 
                             y = a)) +
  geom_contour_filled(aes(z = maxtime)) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(b~.) +
  theme_bw() +
  labs(title = "Peak Time", subtitle = "Burst Size") +
  ylab("Infectio Rate") +
  xlab("Lysis Time") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  NULL

#3
ggplot(data = sum_sims2, aes(y = b, 
                             x = a)) +
  geom_contour_filled(aes(z = maxtime)) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(tau~.) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Peak Time", subtitle = "Lysis Time") +
  ylab("Busrt Size") +
  xlab("Infection Rate") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  NULL

#Plot extin time using contours
#1
ggplot(data = sum_sims2, aes(y = b, 
                             x = tau)) +
  geom_contour_filled(aes(z = extin_time)) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(a~.) +
  theme_bw() +
  labs(title = "Extinction Time", subtitle = "Infection Rate") +
  ylab("Burst Size") +
  xlab("Lysis Time") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  NULL

#2
ggplot(data = sum_sims2, aes(x = tau, 
                             y = a)) +
  geom_contour_filled(aes(z = extin_time)) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(b~.) +
  theme_bw() +
  labs(title = "Extinction Time", subtitle = "Burst Size") +
  ylab("Infection Rate") +
  xlab("Lysis Time") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  NULL

#3
ggplot(data = sum_sims2, aes(x = a, 
                             y = b)) +
  geom_contour_filled(aes(z = extin_time)) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(tau~.) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Extinction Time", subtitle = "Lysis Time") +
  ylab("Burst Size") +
  xlab("Infection Rate") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  NULL


## Let's run sims3 ----
## Where we'll have c, K, a, and r constant, and we'll be changing between 
## different values of b, and tau

# Caluclating the b values with R
tau <- c(30, 45, 62, 87, 102)
intercept <- c(7, 16, 23)
slope <- c(0.932, 0.85, 0.715)
bvals <- as.data.frame(matrix(data = NA, ncol = 4, 
                              nrow = length(tau)*length(intercept)*length(slope)))
bvals
i <- 1

for (tauval in tau){
  for (inter in intercept){
    for (slop in slope){
      b <- slop*(tauval-inter)
      bvals[i,] <- c(tauval, inter, slop, b)
      i <- i+1
    }
  }
}

# This step was made to change the name of the axis and the name of the facets
colnames(bvals) <- c("tau", "tradeintercept", "tradeslope", "b")
bvals
class(bvals)

# Now, we caculate it with the names changed
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sims3.1_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sims3.1 <- list(NULL, NULL, NULL)
  sims3.1[[1]] <- read.csv("./Celia/sims3.1_1.csv", stringsAsFactors = F)
  if ("sims3.1_2.csv" %in% list.files("./Celia/")) {
    sims3.1[[2]] <- read.csv("./Celia/sims3.1_2.csv", stringsAsFactors = F)
  }
  if ("sims3.1_3.csv" %in% list.files("./Celia/")) {
    sims3.1[[3]] <- read.csv("./Celia/sims3.1_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  sims3.1 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                      rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
  #Save results so they can be re-loaded in future
  write.csv(sims3.1[[1]], "./Celia/sims3.1_1.csv", row.names = F)
  if (!is.null(sims3.1[[2]])) {write.csv(sims3.1[[2]], "./Celia/sims3.1_2.csv", row.names = F)}
  if (!is.null(sims3.1[[3]])) {write.csv(sims3.1[[3]], "./Celia/sims3.1_3.csv", row.names = F)}
}

length(sims3.1)
sims3.1[[1]]
sims3.1[[2]]
sims3.1[[3]]

## Let's group_by these simulations
sub_sims3.1 <- subset(sims3.1[[1]], Population == "B")
class(sub_sims3.1)

group_sims3.1 <- dplyr::group_by(sub_sims3.1, uniq_run, a, b, c, K, tau, r)
group_sims3.1

sum_sims3.1 <- dplyr::summarise(group_sims3.1, maximum_B = max(Density),
                                maxtime = time[Density == maximum_B])

# Here, we cut the slope because we don't need it in our simulation
View(sum_sims3.1)

# This step is useful to combine to data frames that have interesting columns
# that we want to plot together
joined_sims3.1 <- left_join(sum_sims3.1, bvals)
joined_sims3.1

## Let's make the ggplot for this data
ggplot(data = joined_sims3.1, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

# Caluclating the b values with R
tau <- c(15, 18, 21.59999, 25.92, 31.104)
intercept <- c(7, 16, 23)
slope <- c(0.932, 0.85, 0.715)
bvals <- as.data.frame(matrix(data = NA, ncol = 4, 
                              nrow = length(tau)*length(intercept)*length(slope)))
bvals
i <- 1

for (tauval in tau){
  for (inter in intercept){
    for (slop in slope){
      b <- slop*(tauval-inter)
      bvals[i,] <- c(tauval, inter, slop, b)
      i <- i+1
    }
  }
}

# This step was made to change the name of the axis and the name of the facets
colnames(bvals) <- c("tau", "tradeintercept", "tradeslope", "b")
bvals
class(bvals)

bvals <- subset(bvals, b > 0)
bvals
# Now, we caculate it with the names changed
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sims3.2_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sims3.2 <- list(NULL, NULL, NULL)
  sims3.2[[1]] <- read.csv("./Celia/sims3.2_1.csv", stringsAsFactors = F)
  if ("sims3.2_2.csv" %in% list.files("./Celia/")) {
    sims3.2[[2]] <- read.csv("./Celia/sims3.2_2.csv", stringsAsFactors = F)
  }
  if ("sims3.2_3.csv" %in% list.files("./Celia/")) {
    sims3.2[[3]] <- read.csv("./Celia/sims3.2_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  sims3.2 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                      rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
  #Save results so they can be re-loaded in future
  write.csv(sims3.2[[1]], "./Celia/sims3.2_1.csv", row.names = F)
  if (!is.null(sims3.2[[2]])) {write.csv(sims3.2[[2]], "./Celia/sims3.2_2.csv", row.names = F)}
  if (!is.null(sims3.2[[3]])) {write.csv(sims3.2[[3]], "./Celia/sims3.2_3.csv", row.names = F)}
}

length(sims3.2)
sims3.2[[1]]
sims3.2[[2]]
sims3.2[[3]]

## Let's group_by these simulations
sub_sims3.2 <- subset(sims3.2[[1]], Population == "B")
class(sub_sims3.2)

group_sims3.2 <- dplyr::group_by(sub_sims3.2, uniq_run, a, b, c, K, tau, r)
group_sims3.2

sum_sims3.2 <- dplyr::summarise(group_sims3.2, maximum_B = max(Density),
                              maxtime = time[Density == maximum_B])
# Here, we cut the slope because we don't need it in our simulation
View(sum_sims3.2)

# This step is useful to combine to data frames that have interesting columns
# that we want to plot together
joined_sims3.2 <- left_join(sum_sims3.2, bvals)
joined_sims3.2

## Let's make the ggplot for this data
ggplot(data = joined_sims3.2, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

# Caluclating the b values with R
tau <- c(22, 25.5, 29.5568, 34.259, 39.709)
intercept <- c(7, 16, 23)
slope <- c(0.932, 0.85, 0.715)
bvals <- as.data.frame(matrix(data = NA, ncol = 4, 
                              nrow = length(tau)*length(intercept)*length(slope)))
bvals
i <- 1

for (tauval in tau){
  for (inter in intercept){
    for (slop in slope){
      b <- slop*(tauval-inter)
      bvals[i,] <- c(tauval, inter, slop, b)
      i <- i+1
    }
  }
}

# This step was made to change the name of the axis and the name of the facets
colnames(bvals) <- c("tau", "tradeintercept", "tradeslope", "b")
bvals
class(bvals)

bvals <- bvals <- subset(bvals, b > 0)
bvals
# Now, we caculate it with the names changed
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sims3.3_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sims3.3 <- list(NULL, NULL, NULL)
  sims3.3[[1]] <- read.csv("./Celia/sims3.3_1.csv", stringsAsFactors = F)
  if ("sims3.3_2.csv" %in% list.files("./Celia/")) {
    sims3.3[[2]] <- read.csv("./Celia/sims3.3_2.csv", stringsAsFactors = F)
  }
  if ("sims3.3_3.csv" %in% list.files("./Celia/")) {
    sims3.3[[3]] <- read.csv("./Celia/sims3.3_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  sims3.3 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                      rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
  #Save results so they can be re-loaded in future
  write.csv(sims3.3[[1]], "./Celia/sims3.3_1.csv", row.names = F)
  if (!is.null(sims3.3[[2]])) {write.csv(sims3.3[[2]], "./Celia/sims3.3_2.csv", row.names = F)}
  if (!is.null(sims3.3[[3]])) {write.csv(sims3.3[[3]], "./Celia/sims3.3_3.csv", row.names = F)}
}

length(sims3.3)
sims3.3[[1]]
sims3.3[[2]]
sims3.3[[3]]

## Let's group_by these simulations
sub_sims3.3 <- subset(sims3.3[[1]], Population == "B")
class(sub_sims3.3)

group_sims3.3 <- dplyr::group_by(sub_sims3.3, uniq_run, a, b, c, K, tau, r)
group_sims3.3

sum_sims3.3 <- dplyr::summarise(group_sims3.3, maximum_B = max(Density),
                              maxtime = time[Density == maximum_B])
# Here, we cut the slope because we don't need it in our simulation
View(sum_sims3.3)

# This step is useful to combine to data frames that have interesting columns
# that we want to plot together
joined_sims3.3 <- left_join(sum_sims3.3, bvals)
joined_sims3.3

## Let's make the ggplot for this data
ggplot(data = joined_sims3.3, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

## Take all the summarized data frames in section sims3, rbind them, and make a 
## big graph to compare all of them together.
sims3 <- rbind(joined_sims3.1, joined_sims3.2, joined_sims3.3)
sims3
View(sims3)
class(sims3)

## Let's make the ggplot for this data
ggplot(data = sims3, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

# Finding the minimums and slopes in this graph
group_sims3 <- dplyr::group_by(sims3, tradeintercept, tradeslope)

sum_sims3 <- 
  dplyr::summarise(group_sims3, minmaxtime = min(maxtime),
                   average_optimal_tau = mean(tau[maxtime  == minmaxtime]),
                   postoptimal_slope = lm(maxtime[tau > average_optimal_tau] ~
                                            tau[tau > average_optimal_tau])$coefficients[2],
                   postoptimal_intercept = lm(maxtime[tau > average_optimal_tau] ~
                                                tau[tau > average_optimal_tau])$coefficients[1])
sum_sims3

## Plot sum_sims3
ggplot(data = sims3, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  geom_abline(data = sum_sims3, 
              mapping = aes(slope = postoptimal_slope, intercept = postoptimal_intercept)) +
  facet_grid(tradeslope ~ tradeintercept) +
  #scale_y_continuous(trans = "log10") +
  NULL
  
## The following graphs let us see if tradeintercept and/or tradeslope affect the
## value of maxtime
ggplot(data = sum_sims3, aes(y = minmaxtime, x = tradeslope, 
                            color = as.factor(tradeintercept))) +
  geom_point() +
  geom_line()

ggplot(data = sum_sims3, aes(y = average_optimal_tau, x = tradeslope, 
                             color = as.factor(tradeintercept))) +
  geom_point() +
  geom_line()

ggplot(data = sum_sims3, aes(x = tradeslope, y = tradeintercept,
                             color = average_optimal_tau,
                             size = minmaxtime)) +
  geom_point()

# If we wanted to plot a regular graph (which we already have in the Word sheet)
# we should follow the steps followed in the other simulations

# Create multiple linear regressions with interactions
lm_fit3 <- lm(maxtime ~ log10(b)*log10(tau)*log10(a), data = sum_sims3)
summary(lm_fit3)
# Save predictions of the model in the new data frame together with the variable
# you want to plot against
sum_sims3_predicted3 <- data.frame(maxtime_pred = predict(lm_fit3, sum_sims3),
                                   tau = sum_sims3$tau,
                                   b = sum_sims3$b,
                                   a = sum_sims3$a)
sum_sims3_predicted3
# This is the predicted line of multiple linear regressions
ggplot(data = sum_sims3, aes(x = log10(b), y = maxtime, 
                             color = as.factor(log10(a)))) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(data = sum_sims3_predicted3, aes(x = log10(b), y = maxtime_pred, 
                                             color = as.factor(log10(a)))) +
  facet_grid(tau ~ .)



## Let's run a big sims3 with new and smaller values of tau ----
tau <- c(6, 9, 13.5, 15, 18, 20.25, 22, 25.5, 28, 30, 32, 34.259, 39.709, 45,
         56, 62, 87, 100, 115, 130)
intercept <- c(7, 16, 23)
slope <- c(0.932, 0.85, 0.715)
bvals <- as.data.frame(matrix(data = NA, ncol = 4, 
                              nrow = length(tau)*length(intercept)*length(slope)))
bvals
i <- 1

for (tauval in tau){
  for (inter in intercept){
    for (slop in slope){
      b <- slop*(tauval-inter)
      bvals[i,] <- c(tauval, inter, slop, b)
      i <- i+1
    }
  }
}

# This step was made to change the name of the axis and the name of the facets
colnames(bvals) <- c("tau", "tradeintercept", "tradeslope", "b")
bvals
class(bvals)

bvals <- bvals <- subset(bvals, b > 0)
bvals
# Now, we caculate it with the names changed
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sims3BIG2_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sims3BIG2 <- list(NULL, NULL, NULL)
  sims3BIG2[[1]] <- read.csv("./Celia/sims3BIG2_1.csv", stringsAsFactors = F)
  if ("sims3BIG2_2.csv" %in% list.files("./Celia/")) {
    sims3BIG2[[2]] <- read.csv("./Celia/sims3BIG2_2.csv", stringsAsFactors = F)
  }
  if ("sims3BIG2_3.csv" %in% list.files("./Celia/")) {
    sims3BIG2[[3]] <- read.csv("./Celia/sims3BIG2_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  sims3BIG2 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                      rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
  #Save results so they can be re-loaded in future
  write.csv(sims3BIG2[[1]], "./Celia/sims3BIG2_1.csv", row.names = F)
  if (!is.null(sims3BIG2[[2]])) {write.csv(sims3BIG2[[2]], "./Celia/sims3BIG2_2.csv", row.names = F)}
  if (!is.null(sims3BIG2[[3]])) {write.csv(sims3BIG2[[3]], "./Celia/sims3BIG2_3.csv", row.names = F)}
}

length(sims3BIG2)
sims3BIG2[[1]]
sims3BIG2[[2]]
sims3BIG2[[3]]
sims3BIG2

## Let's group_by these simulations
sims3BIG2_plot <- sims3BIG2[[1]]
class(sims3BIG2_plot)
#This pivot_wider is a stopgap because the original simulations were
# not run using run_sims and were analyzed in the wider format
sims3BIG2 <- pivot_wider(sims3BIG2_plot, names_from = Population, values_from = Density)

#Make plots
for (my_run in unique(sims3BIG2_plot$uniq_run)) {
  dir.create("./Celia/sims3BIG2_plots/", showWarnings = FALSE)
  tiff(paste("./Celia/sims3BIG2_plots/", my_run, ".tiff", sep = ""),
       width = 4, height = 4, units = "in", res = 200)
  print(ggplot(data = sims3BIG2_plot[sims3BIG2_plot$uniq_run == my_run &
                                   sims3BIG2_plot$Population != "PI", ],
               aes(x = time, y = Density + 10, color = Population)) +
          geom_line(lwd = 1.5, alpha = 1/2) +
          scale_color_manual(values = my_cols[c(8, 2, 3, 1)]) +
          scale_y_continuous(trans = "log10") +
          ggtitle(paste("Run #", my_run, sep = "")))
  dev.off()
}

group_sims3BIG2 <- dplyr::group_by(sims3BIG2, uniq_run, a, b, c, K, tau, r)
group_sims3BIG2

sum_sims3BIG2 <- dplyr::summarise(group_sims3BIG2, maximum_B = max(B),
                                maxtime = time[B == maximum_B])
## What we do next is to have the same number of digits in the columns tau and b
## from both dataframes bvals and sum_sims3BIG2. Otherwise, there's a rounding
## error
bvals$tau <- round(bvals$tau, digits = 2)
sum_sims3BIG2$tau <- round(sum_sims3BIG2$tau, digits = 2)

bvals$b <- round(bvals$b, digits = 2)
sum_sims3BIG2$b <- round(sum_sims3BIG2$b, digits = 2)

class(sum_sims3BIG2$tau)

# This step is useful to combine to data frames that have interesting columns
# that we want to plot together
joined_sims3BIG2 <- left_join(sum_sims3BIG2, bvals)
View(joined_sims3BIG2)

## Plot joined_sims3BIG2
ggplot(data = joined_sims3BIG2, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept)
  #scale_y_continuous(trans = "log10")

# Let's find the minimum maxtime and the average optimal tau for this plot
group_sims3BIG2 <- dplyr::group_by(joined_sims3BIG2, tradeintercept, tradeslope)
sum_sims3BIG2 <- dplyr::summarise(group_sims3BIG2, minmaxtime = min(maxtime),
                                average_optimal_tau = mean(tau[maxtime  == minmaxtime]),
                                slope = lm(maxtime[tau > 30] ~
                                             tau[tau > 30])$coefficients[2],
                                intercept = lm(maxtime[tau > 30] ~
                                                 tau[tau > 30])$coefficients[1])

sum_sims3BIG2

#Let's plot sum_sims3BIG2 with its slope
ggplot(data = joined_sims3BIG2, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  geom_abline(data = sum_sims3BIG2, 
              mapping = aes(slope = slope, intercept = intercept)) +
  facet_grid(tradeslope ~ tradeintercept)
#scale_y_continuous(trans = "log10") +

#Let's plot the different columns in sum_sims3BIG2
ggplot(data = sum_sims3BIG2, aes(x = tradeslope, y = minmaxtime,
                                 colour = as.factor(tradeintercept))) +
  geom_point(size = 3, alpha = 1/2) +
  geom_line(data = sum_sims3BIG2, aes(x = tradeslope, y = minmaxtime, 
                                             color = as.factor(tradeintercept)))



## Let's run sims4 ----
## Where we'll have c and K constant, and we'll be changing between 3 different
## values of b, a, r and tau evenly spaced
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sims4_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sims4 <- list(NULL, NULL, NULL)
  sims4[[1]] <- read.csv("./Celia/sims4_1.csv", stringsAsFactors = F)
  if ("sims4_2.csv" %in% list.files("./Celia/")) {
    sims4[[2]] <- read.csv("./Celia/sims4_2.csv", stringsAsFactors = F)
  }
  if ("sims4_3.csv" %in% list.files("./Celia/")) {
    sims4[[3]] <- read.csv("./Celia/sims4_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  sims4 <- run_sims(bvals = c(100, 200, 400), rvals = c(0.009, 0.016, 0.02845),
                    avals = c(10**-11, 10**-10, 10**-9), kvals = c(10**9),
                    tauvals = c(22.5, 33.75, 50.625))
  #Save results so they can be re-loaded in future
  write.csv(sims4[[1]], "./Celia/sims4_1.csv", row.names = F)
  if (!is.null(sims4[[2]])) {write.csv(sims4[[2]], "./Celia/sims4_2.csv", row.names = F)}
  if (!is.null(sims4[[3]])) {write.csv(sims4[[3]], "./Celia/sims4_3.csv", row.names = F)}
}
             

length(sims4)                  
sims4[[1]]
sims4[[2]]
sims4[[3]]
table(sims4[[1]]$Population)

## Now that we're sure that everything went well, we'll strat summarizing the data
sims4_plot <- sims4[[1]]
class(sims4_plot)
#This pivot_wider is a stopgap because the original simulations were
# not run using run_sims and were analyzed in the wider format
sims4 <- pivot_wider(sims4_plot, names_from = Population, values_from = Density)

#Make plots
for (my_run in unique(sims4_plot$uniq_run)) {
  dir.create("./Celia/sims4_plots/", showWarnings = FALSE)
  tiff(paste("./Celia/sims4_plots/", my_run, ".tiff", sep = ""),
       width = 4, height = 4, units = "in", res = 200)
  print(ggplot(data = sims4_plot[sims4_plot$uniq_run == my_run &
                                   sims4_plot$Population != "PI", ],
               aes(x = time, y = Density + 10, color = Population)) +
          geom_line(lwd = 1.5, alpha = 1/2) +
          scale_color_manual(values = my_cols[c(8, 2, 3, 1)]) +
          scale_y_continuous(trans = "log10") +
          ggtitle(paste("Run #", my_run, sep = "")))
  dev.off()
}

group_sims4 <- dplyr::group_by(sims4, uniq_run, a, b, c, K, tau, r)
group_sims4

sum_sims4 <- dplyr::summarise(group_sims4, maximum_B = max(B),
                              maxtime = time[B == maximum_B],
                              extin_time = time[min(which(time > maxtime & B < 10**4))],
                              slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                           time[time < maxtime & B < 0.1*K])$coefficients[2],
                              intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                               time[time < maxtime & B < 0.1*K])$coefficients[1])
sum_sims4

for (row in 1:nrow(sum_sims4)) {
  bigfinal_rows <- which(sum_sims4$b[row] == group_sims4$b & 
                           sum_sims4$tau[row] == group_sims4$tau &
                           sum_sims4$a[row] == group_sims4$a &
                           sum_sims4$r[row] == group_sims4$r &
                           sum_sims4$K[row] == group_sims4$K &
                           sum_sims4$c[row] == group_sims4$c)
  dir.create("./Celia/sims4_slopeplots/", showWarnings = FALSE)
  tiff(paste("./Celia/Sims4_slopeplots/", sum_sims4[row, "uniq_run"], ".tiff", sep = ""),
       width = 6, height = 6, units = "cm", res = 150)
  print(ggplot(data = group_sims4[bigfinal_rows, ],
               aes(x = time, y = B)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims4$slope[row], intercept = sum_sims4$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims4[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          geom_point(data = sum_sims4[row, ], aes(x = extin_time, y = 10**4),
                     col = "green", size = 3) +
          NULL
  )
  
  dev.off()
}

## Let's make the ggplot for this data
ggplot(data = sum_sims4, aes(x = log10(b), y = maxtime, color = as.factor(a),
                             shape = as.factor(K))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r) +
  geom_smooth(method = "lm")

# How to see how the parameters affect maxtime as if they were INDEPENDENT from 
# each other?
reg_sims4 <- lm(maxtime ~ tau + log10(b) + log10(a) + log10(r), data = sum_sims4)
summary(reg_sims4)

# Calulate the RSE
sigma(reg_sims4)/mean(sum_sims4$maxtime) # The lower the RSE, the more accurate the model.
# The RSE estimate gives a measure of error of prediction.

# Analyzing the data as if the parameters were CORRELATED
# INTERACTION EFFECTS
# Build the model
regi_sims4 <- lm(maxtime ~ tau + log10(b) + log10(a) + log10(r) + tau:log10(b) +
                 tau:log10(a) + tau:log10(r) + log10(b):log10(a) + log10(b):log10(r) +
                   log10(a):log10(r), data = sum_sims4)

summary(regi_sims4)



## Let's run sims5 ----
## It'll be similar to sims3 but with new values of tradeslope and tradeintercept
# To calculate the values of tau evenly spaced from 11 to 120, we'll use the
# following command:
seq(from = 11, to = 120, length.out = 10)
tau <- c(11.00000, 23.11111,  35.22222,  47.33333,  59.44444,  71.55556,
         83.66667,  95.77778, 107.88889, 120.00000)
intercept <- c(10, 20, 30)
slope <- c(5, 12, 19, 26, 33)
bvals <- as.data.frame(matrix(data = NA, ncol = 4, 
                              nrow = length(tau)*length(intercept)*length(slope)))
bvals
i <- 1

for (tauval in tau){
  for (inter in intercept){
    for (slop in slope){
      b <- slop*(tauval-inter)
      bvals[i,] <- c(tauval, inter, slop, b)
      i <- i+1
    }
  }
}

# This step was made to change the name of the axis and the name of the facets
colnames(bvals) <- c("tau", "tradeintercept", "tradeslope", "b")
bvals
class(bvals)

bvals <- subset(bvals, b > 0)
bvals

# Now, we caculate it with the names changed
#Check if saved simulation results exist. If so, load. If not, run simulation
if ("sims5_1.csv" %in% list.files("./Celia/")) {
  #Load results previously simulated
  sims5 <- list(NULL, NULL, NULL)
  sims5[[1]] <- read.csv("./Celia/sims5_1.csv", stringsAsFactors = F)
  if ("sims5_2.csv" %in% list.files("./Celia/")) {
    sims5[[2]] <- read.csv("./Celia/sims5_2.csv", stringsAsFactors = F)
  }
  if ("sims5_3.csv" %in% list.files("./Celia/")) {
    sims5[[3]] <- read.csv("./Celia/sims5_3.csv", stringsAsFactors = F)
  }
} else {
  #Run simulations (if files don't exist)
  sims5 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                        rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
  #Save results so they can be re-loaded in future
  write.csv(sims5[[1]], "./Celia/sims5_1.csv", row.names = F)
  if (!is.null(sims5[[2]])) {write.csv(sims5[[2]], "./Celia/sims3BIG2_2.csv", row.names = F)}
  if (!is.null(sims5[[3]])) {write.csv(sims5[[3]], "./Celia/sims3BIG2_3.csv", row.names = F)}
}

length(sims5)
sims5[[1]]
sims5[[2]]
sims5[[3]]
sims5

## Let's group_by these simulations
sims5_plot <- sims5[[1]]
class(sims5_plot)
#This pivot_wider is a stopgap because the original simulations were
# not run using run_sims and were analyzed in the wider format
sims5 <- pivot_wider(sims5_plot, names_from = Pop, values_from = Density)

#Make plots
for (my_run in unique(sims5_plot$uniq_run)) {
  dir.create("./Celia/sims5_plots/", showWarnings = FALSE)
  tiff(paste("./Celia/sims5_plots/", my_run, ".tiff", sep = ""),
       width = 4, height = 4, units = "in", res = 200)
  print(ggplot(data = sims5_plot[sims5_plot$uniq_run == my_run &
                                       sims5_plot$Po != "PI", ],
               aes(x = time, y = Density + 10, color = Pop)) +
          geom_line(lwd = 1.5, alpha = 1/2) +
          scale_color_manual(values = my_cols[c(8, 2, 3, 1)]) +
          scale_y_continuous(trans = "log10") +
          ggtitle(paste("Run #", my_run, sep = "")))
  dev.off()
}

group_sims5 <- dplyr::group_by(sims5, uniq_run, a, b, c, K, tau, r)
group_sims5

sum_sims5 <- dplyr::summarise(group_sims5, maximum_B = max(B),
                                  maxtime = time[B == maximum_B],
                              extin_time = time[min(which(time > maxtime & B < 10**4))])
## What we do next is to have the same number of digits in the columns tau and b
## from both dataframes bvals and sum_sims3BIG2. Otherwise, there's a rounding
## error
bvals$tau <- round(bvals$tau, digits = 2)
sum_sims5$tau <- round(sum_sims5$tau, digits = 2)

bvals$b <- round(bvals$b, digits = 3)
sum_sims5$b <- round(sum_sims5$b, digits = 3)

class(sum_sims5$tau)

# This step is useful to combine to data frames that have interesting columns
# that we want to plot together
joined_sims5 <- left_join(sum_sims5, bvals)
View(joined_sims5)

## Plot joined_sims3BIG2
ggplot(data = joined_sims5, aes(x = tau, y = maxtime, colour = b)) +
  scale_color_continuous(trans = "log10") +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  theme_bw() +
  xlab("Lysis Time") +
  ylab("Maximum Time") +
  theme(legend.title = element_text("Burst Size")) +
  NULL

# Let's find the minimum maxtime and the average optimal tau for this plot
group_sims5 <- dplyr::group_by(joined_sims5, tradeintercept, tradeslope)
sum_sims5 <- dplyr::summarise(group_sims5, minmaxtime = min(maxtime),
                                  average_optimal_tau = mean(tau[maxtime == minmaxtime]),
                              slope = lm(extin_time[tau > 40] ~
                                           tau[tau > 40])$coefficients[2],
                              intercept = lm(extin_time[tau > 40] ~
                                               tau[tau > 40])$coefficients[1])

sum_sims5

#Let's plot sum_sims3BIG2 with its slope
ggplot(data = joined_sims5, aes(x = tau, y = extin_time, colour = log10(b), 
                                theme(legend.title = element_text("Burst Size")))) +
  geom_point(size = 3, alpha = 1/2) +
  geom_abline(data = sum_sims5, 
              mapping = aes(slope = slope, intercept = intercept)) +
  facet_grid(tradeslope ~ tradeintercept) +
  theme_bw() +
  xlab("Tau") +
  ylab("Extintion Time") 
#scale_y_continuous(trans = "log10") +

#Let's plot the different columns in sum_sims3BIG2
ggplot(data = sum_sims5, aes(x = tradeslope, y = average_optimal_tau,
                                 colour = as.factor(tradeintercept))) +
  geom_point(size = 3, alpha = 1/2) +
  geom_line(data = sum_sims5, aes(x = tradeslope, y = average_optimal_tau, 
                                      color = as.factor(tradeintercept)))

# Let's add a column to sum_sims5 that gives us the value of average_optimal_b. ----
# We'll calculate it by using the formula b = slope*(tau-intercept).

sum_sims5$average_optimal_b <- 
  sum_sims5$tradeslope*(sum_sims5$average_optimal_tau - sum_sims5$tradeintercept)
sum_sims5

# Plot maxtime using colors + regression lines from sims5
sum_sims6 <- subset(sum_sims2, a == 10**-10)
sum_sims6

ggplot(data = sum_sims6, aes(x = tau, 
                             y = b)) +
  geom_contour_filled(aes(z = maxtime)) +
  scale_fill_viridis_d(direction = -1, alpha = 0.85) +
  ggtitle("maxtime") +
  geom_abline(data = sum_sims5, 
              mapping = aes(slope = tradeslope, intercept = (-tradeintercept)*tradeslope, 
                            color = as.factor(tradeintercept)), size = 1) +
  geom_point(data = sum_sims5, 
            mapping = aes(x = average_optimal_tau, y = average_optimal_b,
                          color = as.factor(tradeintercept)), size = 3) +
  scale_color_manual(values = my_cols[c(8, 6, 4)]) +
  theme_bw() +
  theme(legend.title = element_text("Peak Time")) +
  xlab("Lysis Time") +
  ylab("Burst Size") +
  #ylim(-600, 800) +
  #xlim(0, NA) +
  NULL

# Let's do the same but plotting extintion time using colors + regression lines
# from sims5
ggplot(data = sum_sims6, aes(x = tau, 
                             y = b)) +
  geom_contour_filled(aes(z = extin_time)) +
  scale_fill_viridis_d(direction = -1, alpha = 0.85) +
  ggtitle("extintion time") +
  geom_abline(data = sum_sims5, 
              mapping = aes(slope = tradeslope, intercept = (-tradeintercept)*tradeslope, 
                            color = as.factor(tradeintercept)), size = 1) +
  geom_point(data = sum_sims5, 
             mapping = aes(x = average_optimal_tau, y = average_optimal_b,
                           color = as.factor(tradeintercept)), size = 3) +
  scale_color_manual(values = my_cols[c(8, 6, 4)]) +
  theme_bw() +
  theme(legend.title = element_text("Extintion Time")) +
  xlab("Lysis Time") +
  ylab("Burst Size") +
  #ylim(-600, 800) +
  #xlim(0, NA) +
  NULL
