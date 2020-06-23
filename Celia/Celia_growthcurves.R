##Example code for saving ----
if (F) {
  sims4 <- run_sims(bvals = c(100, 200, 400), rvals = c(0.009, 0.016, 0.02845),
                    avals = c(10**-11, 10**-10, 10**-9), kvals = c(10**9),
                    tauvals = c(22.5, 33.75, 50.625))
  #Save results so they can be re-loaded in future
  write.csv(sims4[[1]], "./Celia/sims4_1.csv", row.names = F)
  if (!is.null(sims4[[2]])) {write.csv(sims4[[2]], "./Celia/sims4_2.csv", row.names = F)}
  if (!is.null(sims4[[3]])) {write.csv(sims4[[3]], "./Celia/sims4_3.csv", row.names = F)}
} else {
  #Load results previously simulated
  temp1 <- read.csv("./Celia/sims4_1.csv", stringsAsFactors = F)
  if ("./Celia/sims4_2.csv" %in% list.files()) {
    temp2 <- read.csv("./Celia/sims4_2.csv", stringsAsFactors = F)
  } else {temp2 <- NULL}
  if ("./Celia/sims4_3.csv" %in% list.files()) {
    temp3 <- read.csv("./Celia/sims4_3.csv", stringsAsFactors = F)
  } else {temp3 <- NULL}
  sims4 <- list(temp1, temp2, temp3)
}

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
    warning(paste("pop(s)", paste(which(y < parms["thresh_min_dens"]),
                                  collapse = ","), 
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
    warning(paste("pop(s)",paste(which(y > 10**100), collapse = ","), 
                  "exceed max limit, 10^100, returning dY = 0"))
  }
  dY[y > 10**100] <- 0
  
  #From documentation: The return value of func should be a list, whose first 
  #element is a vector containing the derivatives of y with respect to time
  return(list(dY))
}

#Realistic parameter values: ----
#  
#   r ranges from 0.04 (a 17-minute doubling time) to 0.007 (a 90-minute doubling time).
#   K ranges from 10^6 to 10^10, although typically we focus on 10^8 to 10^9.
#   a ranges from 10^-12 to 10^-8.
#   tau ranges from 10 to 105 (mins).
#   b ranges from 5 to 1000.

## MIKE'S CODE ----
## Define function for running simulations across many parameter values
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
                     "equil" = NA, "time" = NA, "Pop" = as.character(NA), 
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
    yinit <- c(S = my_init_bact,
               I = 0,
               P = my_init_bact*my_moi)
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
      #Define times, with lengths & steps doubling for ea j count
      # (so that the number of timepoints returned is constant)
      times <- seq(0, init_time*2**j, init_stepsize*2**j)
      
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
                                             is.nan(yout_list[[2]]$P)), ]
        
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
      yout_list[[2]]$B <- yout_list[[2]]$S + yout_list[[2]]$I
      #Calculate all phage (PI)
      yout_list[[2]]$PI <- yout_list[[2]]$P + yout_list[[2]]$I
      
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
                               variable.name = "Pop",
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

#New bigFINAL run using run_sims ----                 
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
                         combinatorial = FALSE)
bigFINAL_plot <- sim_bigFINAL[[1]]
bigFINAL <- pivot_wider(bigFINAL_plot, names_from = Pop, values_from = Density)

#Make plots
for (my_run in unique(bigFINAL_plot$uniq_run)) {
  dir.create("./Celia/bigFINAL_plots/", showWarnings = FALSE)
  tiff(paste("./Celia/bigFINAL_plots/", my_run, ".tiff", sep = ""),
             width = 4, height = 4, units = "in", res = 200)
  print(ggplot(data = bigFINAL_plot[bigFINAL_plot$uniq_run == my_run &
                                      bigFINAL_plot$Pop != "PI", ],
               aes(x = time, y = Density + 10, color = Pop)) +
          geom_line(lwd = 1.5, alpha = 1/2) +
          scale_color_manual(values = my_cols[c(8, 2, 3, 1)]) +
          scale_y_continuous(trans = "log10") +
          ggtitle(paste("Run #", my_run, sep = "")))
  dev.off()
}

## Here, I achived to summarize every maximum of each of the simulations included
## in the big data frame

## Now, I'll try to find the maximum of the maximums

bigFINAL <- dplyr::group_by(bigFINAL, b, tau, a, r, K, c)
bigFINAL
FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), 
                          maxtime = time[B == maximum_B])
FINAL

## Finding the slope with summarize ----

## This wasn't exactly the way to do it. The right way is the following:
## Finding the slope with summarize

## Here, we try to analyze and plot row by row of the data frame bigFINAL

FINAL <- dplyr::summarise(bigFINAL, maximum_B = max(B), 
                          maxtime = time[B == maximum_B],
                          slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                       time[time < maxtime & B < 0.1*K])$coefficients[2],
                          intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                           time[time < maxtime & B < 0.1*K])$coefficients[1])
FINAL

for (row in 1:nrow(FINAL)) {
  bigfinal_rows <- which(FINAL$b[row] == bigFINAL$b & 
                           FINAL$tau[row] == bigFINAL$tau &
                           FINAL$a[row] == bigFINAL$a &
                           FINAL$r[row] == bigFINAL$r &
                           FINAL$K[row] == bigFINAL$K &
                           FINAL$c[row] == bigFINAL$c)
  print(ggplot(data = bigFINAL[bigfinal_rows, ],
               aes(x = time, y = B)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = FINAL$slope[row], intercept = FINAL$intercept[row],
                      color = "red") +
          geom_point(data = FINAL[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}  

## How to run the simulations with a LOOP ----

mybig <- NA #We do this to erase anything that's in this data frame
myc <- 1
for(myb in c(75, 760)){
  for(mya in c(10**-10.5, 10**-8.25)){
    for(myK in c(10**7.75, 10**8.8)){
      for(mytau in c(33, 95)){
        for(myr in c(0.009, 0.025)){
          timelength <- 250
          timestep <- 1
          runnumber <- 1
          yinit <- c(S = 10**6, I = 0, P = 10**4)
          params <- c(b = myb, a = mya, K = myK, tau = mytau, r = myr, c = myc,
                      warnings = 0, thresh_min_dens = 10**-100)
          keeprunning <- TRUE
          while(keeprunning){
            times <- seq(from = 0, to = timelength, by = timestep)
            myyout <- as.data.frame(dede(y = yinit, times = times, func = derivs, 
                                         parms = params))
            if(tail(myyout$S, 1) < 1000){
              keeprunning <- FALSE
              ## Here, we said that if the density of S is < 1000 we want to stop
              ## the simulation
            }else{
              keeprunning <- TRUE
              timelength <- timelength*2
              timestep <- timestep*2
              ## And, here, that, otherwise, we want it to keep running but doubling
              ## the time length and the timestep
            }
            if(timelength > 11000){
              keeprunning <- FALSE
              ## Finally, we indicate that if the simulation has been running for
              ## 11000, we want it to stop
            }
          }
          
          ## Here, we create the different columns in the data frame and we say 
          ## that they are equal to myb, for example, because we want both values
          ## in myb, not only one.
          myyout$B <- myyout$S + myyout$I
          myyout$burst <- myb
          myyout$infec <- mya
          myyout$K <- myK
          myyout$tau <- mytau
          myyout$r <- myr
          myyout$c <- myc
          myyout$uniq_run <- runnumber
          
          if(is.na(mybig)){
            mybig <- myyout
          }else{
            mybig <- rbind(mybig, myyout)
          }
          ## Here, we specified that, if mybig is empty, we want it to be equal to
          ## myyout (so, we just want to put myyout in it). However, if we have
          ## already run a simulation, and it's inside mybig, what we want is to
          ## add the second simulation to it by using rbind, making mybig bigger.
          
          runnumber <- 1 + runnumber
        }
      }
    }
  }
}

# Let's find the maximums and slopes of the simulations
mybig

# Firts, we summarize the data we have
mybig <- dplyr::group_by(mybig, burst, tau, infec, r, K, c)
mybig

# Let's find the slopes and the inetercepts, too
BIG <- dplyr::summarise(mybig, maximum_B = max(B), 
                              maxtime = time[B == maximum_B],
                              slope = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                           time[time < maxtime & B < 0.1*K])$coefficients[2],
                              intercept = lm(log10(B[time < maxtime & B < 0.1*K]) ~ 
                                               time[time < maxtime & B < 0.1*K])$coefficients[1])
BIG

## Here, we try to analyze and plot row by row of the data frame mybig
for (row in 1:nrow(BIG)) {
  bigfinal_rows <- which(BIG$burst[row] == mybig$burst & 
                           BIG$tau[row] == mybig$tau &
                           BIG$infec[row] == mybig$infec &
                           BIG$r[row] == mybig$r &
                           BIG$K[row] == mybig$K &
                           BIG$c[row] == mybig$c)
  print(ggplot(data = mybig[bigfinal_rows, ],
               aes(x = time, y = B)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = BIG$slope[row], intercept = BIG$intercept[row],
                      color = "red") +
          geom_point(data = BIG[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}

## To see all the values in BIG:
View(BIG)

## Plotting all the parameters with maxtime in the 32 simulations to see how
## they affect it
ggplot(data = BIG, aes(x = log10(burst), y = maxtime, color = as.factor(K),
                       shape = as.factor(infec))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r)


## Let's run sims1 ----
sims1 <- run_sims(bvals = c(75, 480, 760), avals = c(10**-10.5, 10**-9, 10**-8.25),
                  kvals = c(10**7.75, 10**8.8), tauvals = c(33, 57, 95),
                  rvals = c(0.009, 0.018, 0.025))
length(sims1)
# When we use [] means that the output will be a list, while when we use [[]], means
# that it'll be a dataframe (explanation with trains and their cars).
class(sims1[[1]])
sims1[1]
sims1[[1]]
table(sims1[[1]]$Pop)

# By chencking group 2 and 3 we see if there has been any error: simulations that 
# didn'r reach the equilibrium (2), and simulations that failed (3).
sims1[[2]]
sims1[[3]]

## Notice that the summarize will be slightly different since the data has been
## pivot_longer'd. So, we have to use the density column and make a subset where 
## "Pop" = B

sub_sims1 <- subset(sims1[[1]], Pop == "B")
class(sub_sims1)

## Now, we want to find the maximum_B, the maxtime, and the slope of each simulation
group_sims1 <- dplyr::group_by(sub_sims1, uniq_run, b, tau, a, r, K, c)
group_sims1
class(group_sims1)

sum_sims1 <- dplyr::summarise(group_sims1, maximum_B = max(Density),                
                              maxtime = time[Density == maximum_B],
                              slope = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                           time[time < maxtime & Density < 0.1*K])$coefficients[2],
                              intercept = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                               time[time < maxtime & Density < 0.1*K])$coefficients[1])

for (row in 1:nrow(sum_sims1)) {
  bigfinal_rows <- which(sum_sims1$b[row] == group_sims1$b & 
                           sum_sims1$tau[row] == group_sims1$tau &
                           sum_sims1$a[row] == group_sims1$a &
                           sum_sims1$r[row] == group_sims1$r &
                           sum_sims1$K[row] == group_sims1$K &
                           sum_sims1$c[row] == group_sims1$c)
  tiff(paste("./Celia/Sims1_plots/", sum_sims1[row, "uniq_run"], ".tiff", sep = ""),
       width = 6, height = 6, units = "cm", res = 150)
  
  
  print(ggplot(data = group_sims1[bigfinal_rows, ],
               aes(x = time, y = Density)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims1$slope[row], intercept = sum_sims1$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims1[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
  
  dev.off()
}

group_sims1
sum_sims1
View(sum_sims1)

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
reg_sims1 <- lm(maxtime ~ tau + b + a + r, data = sum_sims1)
summary(reg_sims1)

confint(reg_sims1)
# Calulate the RSE
sigma(reg_sims1)/mean(sum_sims1$maxtime) # The lower the RSE, the more accurate the model.
# The RSE estimate gives a measure of error of prediction.

# Analyzing the data as if the parameters were CORRELATED
# Split the data into training and test set
set.seed(123)
training.samples <- sum_sims1$maxtime %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- sum_sims1[training.samples, ]
test.data <- sum_sims1[-training.samples, ]

# The standard linear regression model can be computed as follow:
# Buil the model
reg_sims1 <- lm(maxtime ~ tau + b + a + r, data = train.data)
# Summarize the model
summary(reg_sims1)
# Make predictions
predictions <- reg_sims1 %>% predict(test.data)
# Make performance. (a) Prediction error, RMSE
RMSE(predictions, test.data$maxtime)
# (b) R2
R2(predictions, test.data$maxtime)

# INTERACTION EFFECTS
# Build the model
reg2_sims1 <- lm(maxtime ~ b + tau + a + r + b:a + a:r, data = sum_sims1)
summary(reg2_sims1)
# Make predictions
predictions <- reg2_sims1 %>% predict(test.data)
# Model performance (a) Prediction error, RMSE
RMSE(predictions, test.data$maxtime)
# (b) R2
R2(predictions, test.data$maxtime)


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
sims2 <- run_sims(bvals = c(50, 100, 200, 400, 800),
                  avals = c(10**-12, 10**-11, 10**-10, 10**-9, 10**-8),
                  kvals = c(10**9), rvals = c(0.04),
                  tauvals = c(15, 22.5, 33.75, 50.625, 75.9375))
length(sims2)                  
sims2[[1]]
sims2[[2]]
sims2[[3]]
table(sims2[[1]]$Pop)

## Now that we're sure that everything went well, we'll strat summarizing the data
sub_sims2 <- subset(sims2[[1]], Pop == "B")
class(sub_sims2)
group_sims2 <- dplyr::group_by(sub_sims2, uniq_run, a, b, c, K, tau, r)
group_sims2

sum_sims2 <- dplyr::summarise(group_sims2, maximum_B = max(Density),                
                              maxtime = time[Density == maximum_B],
                              slope = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                           time[time < maxtime & Density < 0.1*K])$coefficients[2],
                              intercept = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                               time[time < maxtime & Density < 0.1*K])$coefficients[1])

for (row in 1:nrow(sum_sims2)) {
  bigfinal_rows <- which(sum_sims2$b[row] == group_sims2$b & 
                           sum_sims2$tau[row] == group_sims2$tau &
                           sum_sims2$a[row] == group_sims2$a &
                           sum_sims2$r[row] == group_sims2$r &
                           sum_sims2$K[row] == group_sims2$K &
                           sum_sims2$c[row] == group_sims2$c)
  
  print(ggplot(data = group_sims2[bigfinal_rows, ],
               aes(x = time, y = Density)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims2$slope[row], intercept = sum_sims2$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims2[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}

sum_sims2

## Let's make the ggplot for this data
ggplot(data = sum_sims2, aes(x = log10(b), y = maxtime, color = as.factor(a),
                             shape = as.factor(K))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tau ~ r) +
  geom_smooth(method = "lm")

# How to see how the parameters affect maxtime as if they were INDEPENDENT from 
# each other?
reg_sims2 <- lm(maxtime ~ tau + b + a, data = sum_sims2)
summary(reg_sims2)
plot(reg_sims2)

confint(reg_sims2)
# Calulate the RSE
sigma(reg_sims2)/mean(sum_sims2$maxtime) # The lower the RSE, the more accurate the model.
# The RSE estimate gives a measure of error of prediction.

# Analyzing the data as if the parameters were CORRELATED
# INTERACTION EFFECTS
# Build the model
regi_sims2 <- lm(maxtime ~ tau + b + a + tau:b + tau:a + b:a, data = sum_sims2)
                   
summary(regi_sims2)

## We want to plot the multiple linear regressions
# Read data set
sum_sims2
# Create multiple linear regressions
lm_fit <- lm(maxtime ~ log10(b) + log10(a) + log10(tau), data = sum_sims2)
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
                             color = as.factor(log10(a)))) +
  geom_point() +
  geom_line(data = sum_sims2_predicted, aes(x = log10(tau), y = maxtime_pred, 
                                            color = as.factor(log10(a)))) +
  facet_grid(b ~ .)

# Create multiple linear regressions with interactions
lm_fit2 <- lm(maxtime ~ log10(b)*log10(a)*log10(tau), data = sum_sims2)
summary(lm_fit2)
# Save predictions of the model in the new data frame together with the variable
# you want to plot against
sum_sims2_predicted2 <- data.frame(maxtime_pred = predict(lm_fit2, sum_sims2),
                                  tau = sum_sims2$tau,
                                  b = sum_sims2$b,
                                  a = sum_sims2$a)
sum_sims2_predicted2
# This is the predicted line of multiple linear regressions
ggplot(data = sum_sims2, aes(x = log10(tau), y = maxtime, 
                             color = as.factor(log10(a)))) +
  geom_point() +
  geom_line(data = sum_sims2_predicted2, aes(x = log10(tau), y = maxtime_pred, 
                                            color = as.factor(log10(a)))) +
  facet_grid(b ~ .)



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
sims3.1 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                  rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
length(sims3.1)
sims3.1[[1]]
sims3.1[[2]]
sims3.1[[3]]

## Let's group_by these simulations
sub_sims3.1 <- subset(sims3.1[[1]], Pop == "B")
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

bvals <- bvals[-c(4, 5, 6, 7, 8, 9, 16, 17, 18, 25, 26, 27), ]
bvals
# Now, we caculate it with the names changed
sims3.2 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                  rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
length(sims3.2)
sims3.2[[1]]
sims3.2[[2]]
sims3.2[[3]]

## Let's group_by these simulations
sub_sims3.2 <- subset(sims3.2[[1]], Pop == "B")
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

bvals <- bvals[-c(7, 8, 9), ]
bvals
# Now, we caculate it with the names changed
sims3.3 <- run_sims(bvals = bvals$b, avals = c(10**-10), kvals = c(10**9),
                  rvals = c(0.04), tauvals = bvals$tau, combinatorial = FALSE)
length(sims3.3)
sims3.3[[1]]
sims3.3[[2]]
sims3.3[[3]]

## Let's group_by these simulations
sub_sims3.3 <- subset(sims3.3[[1]], Pop == "B")
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

## Take all the summarized data frames in section sims2, rbind them, and make a 
## big grap to compere all of them together.

sims3 <- rbind(joined_sims3.1, joined_sims3.2, joined_sims3.3)
sims3
View(sims3)

## Let's make the ggplot for this data
ggplot(data = sims3, aes(x = tau, y = maxtime, colour = log10(b))) +
  geom_point(size = 3, alpha = 1/2) +
  facet_grid(tradeslope ~ tradeintercept) +
  scale_y_continuous(trans = "log10")

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




## Let's run sims4 ----
## Where we'll have c and K constant, and we'll be changing between 3 different
## values of b, a, r and tau evenly spaced
sims4 <- run_sims(bvals = c(100, 200, 400), rvals = c(0.009, 0.016, 0.02845),
                  avals = c(10**-11, 10**-10, 10**-9), kvals = c(10**9),
                  tauvals = c(22.5, 33.75, 50.625))
length(sims4)                  
sims4[[1]]
sims4[[2]]
sims4[[3]]
table(sims4[[1]]$Pop)

## Now that we're sure that everything went well, we'll strat summarizing the data
sub_sims4 <- subset(sims4[[1]], Pop == "B")
class(sub_sims4)
group_sims4 <- dplyr::group_by(sub_sims4, uniq_run, a, b, c, K, tau, r)
group_sims4

sum_sims4 <- dplyr::summarise(group_sims4, maximum_B = max(Density),                
                              maxtime = time[Density == maximum_B],
                              slope = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                           time[time < maxtime & Density < 0.1*K])$coefficients[2],
                              intercept = lm(log10(Density[time < maxtime & Density < 0.1*K]) ~ 
                                               time[time < maxtime & Density < 0.1*K])$coefficients[1])

for (row in 1:nrow(sum_sims4)) {
  bigfinal_rows <- which(sum_sims4$b[row] == group_sims4$b & 
                           sum_sims4$tau[row] == group_sims4$tau &
                           sum_sims4$a[row] == group_sims4$a &
                           sum_sims4$r[row] == group_sims4$r &
                           sum_sims4$K[row] == group_sims4$K &
                           sum_sims4$c[row] == group_sims4$c)
  
  print(ggplot(data = group_sims4[bigfinal_rows, ],
               aes(x = time, y = Density)) +
          geom_line() +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = sum_sims4$slope[row], intercept = sum_sims4$intercept[row],
                      color = "red") +
          geom_point(data = sum_sims4[row, ], aes(x = maxtime, y = maximum_B), 
                     col = "blue", size = 3) +
          NULL
  )
}

sum_sims4

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