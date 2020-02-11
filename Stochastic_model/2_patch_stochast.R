#set working directory
setwd("Dropbox/CAMBRIDGE/Part_2/Project/Two_patch_model/Stochastic_model/")

# delete everything
rm(list=ls(all=TRUE))

# import necessary libraries
require(adaptivetau)
require(deSolve)
require(flux)

#PARAMETERS
N <- 10^4*100 #total number of plants

epsilon <- 0.5 #degree of coupling; epsilon = 1 means totally mixed
beta <- 9e-6 # infection rate * when no resistant crops, equilibrium is 50% get infected was 9e-6
rho <- 1 #replanting rate; rho = 1 means planted in seasons *dont bother changing
delta <- 0.3 #fitness cost of rb *
gamma <- 0 #completeness of resistance; gamma = 0 means fully resistant *leave at zero
phi <- 0.5 #cropping ratio; proportion of resistant hosts

#TIMEFRAME
tf <- 40
time_points <- seq(0, tf, length.out = 1001)
num_runs <- 100

#seasons to calculate AUC at
seasons <- seq(1, tf, by = 5)

#GROUPS
group.names <- c("Iwt_S", "Irb_S", "Iwt_R", "Irb_R")

wt_eq_prop <- 1/3
initialIrb_prop <- 3e-6

#INITIAL STATES
initIwt_S <- wt_eq_prop * N - (wt_eq_prop * N * initialIrb_prop)
initIrb_S <- wt_eq_prop * N * initialIrb_prop
  
initIwt_R <- 0
initIrb_R <- 0

initS_S <- (1 - phi)*N - initIwt_S - initIrb_S
initS_R <- phi*N - initIwt_R - initIrb_R

#DEFINE TRANSITION RATES AND SIMULATION FUNCTIONS
SI_transition_rates <- function(x, parms, t) {
  
  #parameters
  N <- parms["N"]
  
  epsilon <- parms["epsilon"]
  p <- 1 / (1 + epsilon) #proportion of inoculum landing on same type of host 
                        #(1-p = amount landing on diff host)
  
  beta <- parms["beta"]
  rho <- parms["rho"]
  delta <- parms["delta"]
  gamma <- parms["gamma"]
  phi <- parms["phi"]

  #states
  Iwt_S <- x["Iwt_S"]
  Irb_S <- x["Irb_S"]
  Iwt_R <- x["Iwt_R"]
  Irb_R <- x["Irb_R"]
  S_S <- x["S_S"]
  S_R <- x["S_R"]
 
  dIwt_S_forward <- beta * S_S * (p * Iwt_S + (1 - p) * Iwt_R)
  dIwt_S_backward <- rho * Iwt_S
  
  dIrb_S_forward <- beta * S_S * (1 - delta) * (p * Irb_S + (1 - p) * Irb_R) 
  dIrb_S_backward <- rho * Irb_S
  
  dIwt_R_forward <- beta * S_R * gamma * (p * Iwt_R + (1 - p) * Iwt_S)
  dIwt_R_backward <- rho * Iwt_R
  
  dIrb_R_forward <- beta * S_R * (1 - delta) * (p * Irb_R + (1 - p) * Irb_S)
  dIrb_R_backward <- rho * Irb_R

  return(c(
    dIwt_S_forward,
    dIwt_S_backward, #wt infection rate of susceptible plants
    dIrb_S_forward,
    dIrb_S_backward,#rb infection rate of susceptible plants
    dIwt_R_forward,
    dIwt_R_backward,#wt infection rate of resistant plants
    dIrb_R_forward,
    dIrb_R_backward #rb infection rate of resistant plants
  ))
}

run_simulation <- function(tf) {
  
  #set.seed?? - done so results are reproducible
  
  SI_transitions <- list(
    c(S_S = -1, Iwt_S = +1),
    c(S_S = +1, Iwt_S = -1), #wt infection rate of susceptible plants
    c(S_S = -1, Irb_S = +1),
    c(S_S = +1, Irb_S = -1), #rb infection rate of susceptible plants
    c(S_R = -1, Iwt_R = +1),
    c(S_R = +1, Iwt_R = -1), #wt infection rate of resistant plants
    c(S_R = -1, Irb_R = +1),
    c(S_R = +1, Irb_R = -1) #rb infection rate of resistant plants
  )
  
  run <- ssa.adaptivetau(init.values = c(Iwt_S = round(initIwt_S), Irb_S = round(initIrb_S), 
                                         Iwt_R = round(initIwt_R), Irb_R = round(initIrb_R),
                                         S_S = round(initS_S), S_R = round(initS_R)),
                         transitions = SI_transitions,
                         rateFunc = SI_transition_rates,
                         params = c(N = N, beta = beta,epsilon = epsilon, beta = beta, 
                                    rho = rho, delta = delta, gamma = gamma, phi = phi),
                         tf = tf#,#Inf = until all transition rates are 0
                         #tl.params = list(epsilon = 0.005/1000)
                          )
  
  run_df <- data.frame(run)
  
  run.times <- vector("list", length(group.names))
  names(run.times) <- group.names
  
  for (group in 1:length(run.times)) {
    run.time <- approx(x = run_df$time,
                    y = run_df[,group+1], #because first column is time
                    xout = time_points,
                    method = "constant")
    run.times[[group]] <- data.frame(run.time[1], run.time[2])
    colnames(run.times[[group]]) <- c("time", names(run.times[group]))
  }
  
  return(run.times)
}

#FUNCTION RUNS MODEL num_runs TIMES AND RETURNS ALL SIMULATIONS (RECORDS AS PROPORTIONS)
simulations <- function(tf) {

  #RUN SIMULATION MANY TIMES AND CREATE MATRICES OF RESULTS FOR I VARIABLES
  simulations <- vector("list", length(group.names))
  names(simulations) <- group.names
  
  extinction_frequencies <- vector("list", length = length(group.names)) #NOTE: doesnt include Iwt_R
  group_end_extinction_count <- c(0,0,0,0)
  extinction_count <- 0
  
  
  for (run in 1:num_runs) {
    trajectory <- run_simulation(tf)
    
    for (group in 1:length(trajectory)) {
      group_trajectory <- trajectory[[group]]
      group_extinction_count <- 0
      
      #if trajectory of group ends at 0 add 1 to end_extinction_count
      if (group_trajectory[group_trajectory["time"] == tf][2] == 0) {
        group_end_extinction_count[group] <- group_end_extinction_count[group] + 1
      }
      
      #iterate through rows and add 1 to extinction count whenever the frequency drops from 1 to 0
      for (row in 1:nrow(group_trajectory)) {
        if (group_trajectory[row, 2] == 1) { #2 refers to column with number of infected of the group
          if (row != nrow(group_trajectory)) { #ensure not at end
            if (group_trajectory[row+1, 2] == 0) { #ie. if extinction event has occured
              group_extinction_count <- group_extinction_count + 1
            }
          }
        }
      }
      #after iterating through all rows, add total extinction count for this run to appropriate vector in extinction list
      extinction_frequencies[[group]] <- c(extinction_frequencies[[group]], group_extinction_count)
      
      #find proportion infected and add to simulations matrix
      group_props_df <- group_trajectory[,-1] / N #divide all columns except col 1 (time) by N
      group_props_df <- cbind(group_trajectory[,1], group_props_df) #add time col back on
      group_matrix <- data.matrix(group_props_df) #data.matrix turns a dataframe to a matrix
      
      if (run == 1) {
        simulations[[group]] <- group_matrix
      }
      else {
        simulations[[group]] <- cbind(simulations[[group]], group_matrix[,2]) #not including time col
      }
    }
    
    #Iwt_S_trajectory <- trajectory[["Iwt_S"]]
    Irb_S_trajectory <- trajectory[["Irb_S"]]
    Irb_R_trajectory <- trajectory[["Irb_R"]]
    
    #if this run led to extinction of rb pathogen add 1 to extinction count
    if (Irb_S_trajectory[Irb_S_trajectory["time"] == tf][2] == 0 #2 gets value of Irb_S
        & Irb_R_trajectory[Irb_R_trajectory["time"] == tf][2] == 0) {
      
      extinction_count <- extinction_count + 1
    }
  }
  
  #after all runs, find average extinction count for each group
  average_extinctions <- c(mean(extinction_frequencies[[1]]), #Iwt_S
                           mean(extinction_frequencies[[2]]), #Irb_S
                           mean(extinction_frequencies[[4]])) #Irb_R
  
  #find extinction proportions
  extinction_prop <- extinction_count / num_runs
  group_end_extinction_prop <- group_end_extinction_count / num_runs
  
  #write csv files
  simulations_output <- c(simulations, extinction_prop)
  write.csv(simulations_output, file = 
            paste0("values_csvs/simulations_output_num=",simulations_called,".csv")) #paste 0 means no spaces put in
  
  write.csv(average_extinctions, file = paste0("values_csvs/average_extinctions_num=",simulations_called,".csv"))
  
  write.csv(group_end_extinction_prop, file = paste0("values_csvs/group_end_extinctions_num=",simulations_called,".csv"))
  
  return(list(simulations, extinction_prop))
}

#COUNT OF HOW MANY TIMES SIMULATIONS FUNCTION HAS BEEN CALLED
simulations_called <- 0
trace(simulations, tracer=function() simulations_called <<- simulations_called + 1)

#FUNCTION THAT CALCULATES MEAN VALUES AND 5th AND 95th PERCENTILES AT EACH TIME POINT 
#OVER ALL THE SIMULATIONS
simulation_means <- function(tf) {
  
  #run simulations
  simulations_output <- simulations(tf)
  simulations <- simulations_output[[1]]
  
  #create calculations vector
  calculations <- vector("list", length(group.names))
  names(calculations) <- group.names
  
  for (group in 1:length(calculations)) {
    group_simulations <- simulations[[group]]
    for (row in 1:nrow(group_simulations)) {
      
      group_percs <- quantile(group_simulations[row, -1], c(0.05, 0.95)) #-1 means not col 1 as it is the times
      
      if (row == 1) {
        group_mean_vals <- mean(group_simulations[row, -1])
        group_fifth_percs <- group_percs["5%"]
        group_ninetyfifth_percs <- group_percs["95%"]
        calculations[[group]] <- data.frame(group_simulations[row,1], #time_points
                                            group_mean_vals, 
                                            group_fifth_percs, 
                                            group_ninetyfifth_percs)
        colnames(calculations[[group]]) <- c("time_points", "mean", "fifth_perc", "ninetyfifth_perc")
        #rownames(calculations[[group]]) <- time_points
      }
      else {
        group_mean_vals <- mean(group_simulations[row, -1])
        group_fifth_percs <- group_percs["5%"]
        group_ninetyfifth_percs <- group_percs["95%"]
        
        calculations[[group]] <- rbind(calculations[[group]], 
                                       c(group_simulations[row,1], #time_points
                                         group_mean_vals, 
                                         group_fifth_percs, 
                                         group_ninetyfifth_percs))
      }
    }
  }
  write.csv(calculations, 
            file = paste0("values_csvs/means+percs_num=",simulation_means_called,".csv"))
  
  return(calculations)
}

#COUNT OF HOW MANY TIMES SIMULATIONS FUNCTION HAS BEEN CALLED
simulation_means_called <- 0
trace(simulation_means, tracer=function() simulation_means_called <<- simulation_means_called + 1)

#FUNCTION THAT CALUCLATES INDIVIDUAL GROUP AUCS AND TOTAL AUC (SUM) - STOCHASTIC
stochastic_aucs <- function(tf) {
  print(1)
  mean_calculations <- simulation_means(tf) #NOTE: infecteds graph produced below will 
                                            #have different means as it calls 
                                            #simulation_means() again
  print(2)
  #CALCULATE AUC FOR EACH CURVE
  #when accessing means do calculations[["Iwt_S"]]$mean etc
  aucs <- vector(length = 4)
  for (group in 1:length(group.names)) {
    group_means <- mean_calculations[[group]]
    group_auc <- auc(group_means$time_points, group_means$mean)
    
    #AUC at end of time frame for group
    aucs[group] <- group_auc
    
    #initialise vector of AUCs at each time point
    seasons_aucs <- vector(length = length(seasons))
    for (season in 1:length(seasons)) {
      
      #calculate AUC up to the season
      seasons_aucs[season] <- auc(group_means[group_means$time_points <= seasons[season],]$time_points, 
                                  group_means[group_means$time_points <= seasons[season],]$mean)
    }
    if (group == 1) {
      all_groups_seasons_aucs <- matrix(seasons_aucs, ncol = 1)
    } else {
      all_groups_seasons_aucs <- cbind(all_groups_seasons_aucs, seasons_aucs)
    }
      
  }
  #total AUC at end of time frame
  total_auc <- sum(aucs)
  
  #calculate total AUC at each season selected
  seasons_total_aucs <- vector(length = nrow(all_groups_seasons_aucs))
  for (row in 1:nrow(all_groups_seasons_aucs)) {
    seasons_total_aucs[row] <- sum(all_groups_seasons_aucs[row,])
  }
  
  print(3)
  #RETURNS VECTOR OF INDIVIDUAL GROUP AUC VALUES AND TOTAL AUC VALUE (SUM) AND A VECTOR OF SEASON AUCS
  return(list(c(aucs, total_auc), seasons_total_aucs))
}

#PLOT INFECTED AGAINST TIME USING MEAN - STOC
calculations <- simulation_means(tf)

par(mar = c(3.1, 4.1, 3.1, 1.1),
    oma = c(3,0,0,0),
    mfrow = c(1,2)) #margins, default is 5.1(bot),4.1(l),4.1(top),2.1(r)
plot(x = time_points, y = calculations[["Iwt_S"]]$mean, type = "l",
     col = "blue",
     ylab = "Proportion infected",
     ylim=c(0, #0.35
            #max(calculations[["Iwt_S"]]$ninetyfifth_perc)
            0.4
            ))
grid()
title("a)", adj = 0, line = 0.5)
# lines(x = time_points, y = calculations[["Iwt_S"]]$fifth_perc,
#       col = "lightblue")
# lines(x = time_points, y = calculations[["Iwt_S"]]$ninetyfifth_perc,
#       col = "lightblue")

lines(x = time_points, y = calculations[["Irb_S"]]$mean,
      col = "red")
# lines(x = time_points, y = calculations[["Irb_S"]]$fifth_perc,
#       col = "pink")
# lines(x = time_points, y = calculations[["Irb_S"]]$ninetyfifth_perc,
#       col = "pink")

###lines(x = time_points, y = calculations[["Iwt_R"]]$mean,
      ###col = "green")
###lines(x = time_points, y = calculations[["Iwt_R"]]$fifth_perc,
      ###col = "lightgreen")
###lines(x = time_points, y = calculations[["Iwt_R"]]$ninetyfifth_perc,
      ###col = "lightgreen")

lines(x = time_points, y = calculations[["Irb_R"]]$mean,
      col = "orange")
# lines(x = time_points, y = calculations[["Irb_R"]]$fifth_perc,
#       col = "yellow")
# lines(x = time_points, y = calculations[["Irb_R"]]$ninetyfifth_perc,
#       col = "yellow")



#RUN DETERMINISTIC MODEL TO COMPARE

SI_ode <- function(times, states, parms) {
  
  #parameters
  epsilon <- parms["epsilon"]
  p <- 1 / (1 + epsilon) #proportion of inoculum landing on same type of host 
  #(1-p = amount landing on diff host)
  N <- parms["N"]
  
  beta <- parms["beta"]
  rho <- parms["rho"]
  delta <- parms["delta"]
  gamma <- parms["gamma"]
  phi <- parms["phi"]
  
  #states
  Iwt_S <- states["Iwt_S"]
  Irb_S <- states["Irb_S"]
  Iwt_R <- states["Iwt_R"]
  Irb_R <- states["Irb_R"]

  
  S_R <- phi * N - Iwt_R - Irb_R
  S_S <- (1-phi) * N - Iwt_S - Irb_S
  
  dIwt_S <- beta * S_S * (p * Iwt_S + (1 - p) * Iwt_R) - rho * Iwt_S
  dIrb_S <- beta * S_S * (1 - delta) * (p * Irb_S + (1 - p) * Irb_R) - rho * Irb_S
  
  dIwt_R <- beta * S_R * gamma * (p * Iwt_R + (1 - p) * Iwt_S) - rho * Iwt_R
  dIrb_R <- beta * S_R * (1 - delta) * (p * Irb_R + (1 - p) * Irb_S) - rho * Irb_R
  
  return(list(c(dIwt_S, dIrb_S, dIwt_R, dIrb_R)))
}

#FIND AUC OF DETERMINISTIC MODEL
deterministic_aucs <- function(tf) {
  print("into determ")
  trajectory <- ode(y = c(Iwt_S = initIwt_S, Irb_S = initIrb_S, Iwt_R = initIwt_R, Irb_R = initIrb_R),
  time = time_points,
  parms = c(N = N, epsilon = epsilon, beta = beta, rho = rho, delta = delta, gamma = gamma, phi = phi),
  func = SI_ode)
  
  trajectory_df <- data.frame(trajectory)
  trajectory_props_df <- trajectory_df[,-1] / N #divide cols except col 1 (time) by N
  trajectory_props_df <- cbind(trajectory_df[,1], trajectory_props_df) #add time col back on
  names(trajectory_props_df) <- c("time", group.names) #rename first column

  #PLOT INFECTED AGAINST TIME - DETERM
  par(mar = c(3.1, 2.1, 3.1, 3.1))
  plot(x = time_points, y = trajectory_props_df$Iwt_S, type = "l",
       col = "blue",
       ylab = "",
       ylim=c(0, #0.35
              #max(trajectory_props_df$Iwt_S)
              0.4
       ))
  grid()
  title("b)", adj = 0, line = 0.5)
  lines(x = time_points, y = trajectory_props_df$Irb_S, col = "red")
  ###lines(x = time_points, y = trajectory_props_df$Iwt_R, col = "green")
  lines(x = time_points, y = trajectory_props_df$Irb_R, col = "orange")

  
  print("determ2")
  #FIND AUC OF INDIVIDUAL GROUPS + TOTAL
  aucs <- vector(length = 4)
  for (group in 1:length(group.names)) { 
    groupname <- group.names[group]
    group_auc <- auc(trajectory_props_df$time, trajectory_props_df[,group+1]) #+1 as first col = time
    aucs[group] <- group_auc
    print("determ3")
    
    #initialise vector of AUCs at each time point
    seasons_aucs <- vector(length = length(seasons))
    for (season in 1:length(seasons)) {
      
      #calculate AUC up to the season
      seasons_aucs[season] <- auc(trajectory_props_df[trajectory_props_df$time <= seasons[season],]$time, 
                                  trajectory_props_df[trajectory_props_df$time <= seasons[season], group + 1])
    }
    
    if (group == 1) {
      all_groups_seasons_aucs <- matrix(seasons_aucs, ncol = 1)
    } else {
      all_groups_seasons_aucs <- cbind(all_groups_seasons_aucs, seasons_aucs)
    }
  }
  #calculate total AUC up to tf
  total_auc <- sum(aucs)
  
  #calculate total AUC at each season selected
  seasons_total_aucs <- vector(length = nrow(all_groups_seasons_aucs))
  for (row in 1:nrow(all_groups_seasons_aucs)) {
    seasons_total_aucs[row] <- sum(all_groups_seasons_aucs[row,])
  }
  
  #RETURNS VECTOR OF INDIVIDUAL GROUP AUC VALUES AND TOTAL AUC VALUE (SUM)
  return(list(c(aucs, total_auc), seasons_total_aucs))
}
deterministic_aucs(tf) #will plot infected vs time graph

#LEGEND - make a new plot and overlay it onto INFECTED VS TIME plot with graphs
mtext("Time (seasons)", side = 1, outer = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",
       horiz = TRUE,
       legend= c("WT infecting S", "RB infecting S", #"R infected by wt",
                 "RB infecting by R"),
       col = c("blue", "red", #"green",
               "orange"),
       xpd = TRUE,   #allows drawing outside plot
       bty = "n",
       lty = c(1,1),
       cex = 0.8) #makes legend font 0.8 size of other font

#PLOT MEAN + 5TH + 95TH PERCENTILES FOR EACH GROUP SEPARATELY
#plot for Iwt_S
par(mar = c(5.1, 4.6, 3.1, 2.1),
    oma = c(0,0,0,0),
    mfrow = c(1,3)) #margins, default is 5.1(bot),4.1(l),4.1(top),2.1(r)
plot(x = time_points, y = calculations[["Iwt_S"]]$mean, type = "l",
     col = "blue",
     ylab = "Proportion infected",
     xlab = "Time (seasons)",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
     ylim=c(0, 0.35
            #max(calculations[["Iwt_S"]]$ninetyfifth_perc)
     ))
polygon(c(calculations[["Iwt_S"]]$time_points,
          rev(calculations[["Iwt_S"]]$time_points)),
        c(calculations[["Iwt_S"]]$fifth_perc,
          rev(calculations[["Iwt_S"]]$ninetyfifth_perc)),
        col = adjustcolor("lightblue",alpha.f=0.5),
        border = NA)
grid()
title("c)", adj = 0, line = 0.5, cex.main = 2)
lines(x = time_points, y = calculations[["Iwt_S"]]$fifth_perc,
      col = "lightblue")
lines(x = time_points, y = calculations[["Iwt_S"]]$ninetyfifth_perc,
      col = "lightblue")

#plot for Irb_S
plot(x = time_points, y = calculations[["Irb_S"]]$mean, type = "l",
     col = "red",
     ylab = "", #<<- not needed as c) has y lab and will be adjacent
     xlab = "Time (seasons)",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
     ylim=c(0, 0.35
            #max(calculations[["Irb_S"]]$ninetyfifth_perc)
     ))
polygon(c(calculations[["Irb_S"]]$time_points,
          rev(calculations[["Irb_S"]]$time_points)),
        c(calculations[["Irb_S"]]$fifth_perc,
          rev(calculations[["Irb_S"]]$ninetyfifth_perc)),
        col = adjustcolor("pink",alpha.f=0.5),
        border = NA)
grid()
title("d)", adj = 0, line = 0.5, cex.main = 2)
lines(x = time_points, y = calculations[["Irb_S"]]$fifth_perc,
      col = "pink")
lines(x = time_points, y = calculations[["Irb_S"]]$ninetyfifth_perc,
      col = "pink")

#plot for Irb_R
plot(x = time_points, y = calculations[["Irb_R"]]$mean, type = "l",
     col = "chocolate3",
     ylab = "", #<<- not needed as c) has y lab and will be adjacent
     xlab = "Time (seasons)",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
     ylim=c(0, 0.35
            #max(calculations[["Irb_R"]]$ninetyfifth_perc)
     ))
polygon(c(calculations[["Irb_R"]]$time_points,
          rev(calculations[["Irb_R"]]$time_points)),
        c(calculations[["Irb_R"]]$fifth_perc,
          rev(calculations[["Irb_R"]]$ninetyfifth_perc)),
        col = adjustcolor("khaki1",alpha.f=0.5),
        border = NA)
grid()
title("e)", adj = 0, line = 0.5, cex.main = 2)
lines(x = time_points, y = calculations[["Irb_R"]]$fifth_perc,
      col = "khaki1")
lines(x = time_points, y = calculations[["Irb_R"]]$ninetyfifth_perc,
      col = "khaki1")

#reset count for simulations called
simulations_called <- 0 
simulation_means_called <- 0