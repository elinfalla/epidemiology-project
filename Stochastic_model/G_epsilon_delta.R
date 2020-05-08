#IMPORTANT: RUN 2_PATCH_STOCHAST FIRST

epsilon_vals <- c(0.95, 0.05)
delta_vals <- seq(0,1, length.out = 5)
epsilon_AUCs <- vector(length = length(epsilon_vals))

stoc_epsilon_aucs <- vector(length = length(epsilon_vals))
determ_epsilon_aucs <- vector(length = length(epsilon_vals))

stoc_epsilon_diff <- vector(length = length(delta_vals))
determ_epsilon_diff <- vector(length = length(delta_vals))

high_stoc_epsilon_aucs <- vector(length = length(delta_vals))
high_determ_epsilon_aucs <- vector(length = length(delta_vals))
low_stoc_epsilon_aucs <- vector(length = length(delta_vals))
low_determ_epsilon_aucs <- vector(length = length(delta_vals))

#ensure simulations() and simulation_means() trace counts are set to 0
simulations_called <- 0
simulation_means_called <- 0

#vector of extinction proportion at each delta value
low_epsilon_extinction_props <- vector(length = length(delta_vals))
high_epsilon_extinction_props <- vector(length = length(delta_vals))
#vector of average extinctions per run for each group
low_epsilon_Iwt_S_extinction <- vector(length = length(delta_vals))
high_epsilon_Iwt_S_extinction <- vector(length = length(delta_vals))
low_epsilon_Irb_S_extinction <- vector(length = length(delta_vals))
high_epsilon_Irb_S_extinction <- vector(length = length(delta_vals))
low_epsilon_Irb_R_extinction <- vector(length = length(delta_vals))
high_epsilon_Irb_R_extinction <- vector(length = length(delta_vals))
#vector of prop extinctions for each group
low_epsilon_Iwt_S_end_extinction <- vector(length = length(delta_vals))
high_epsilon_Iwt_S_end_extinction <- vector(length = length(delta_vals))
low_epsilon_Irb_S_end_extinction <- vector(length = length(delta_vals))
high_epsilon_Irb_S_end_extinction <- vector(length = length(delta_vals))
low_epsilon_Irb_R_end_extinction <- vector(length = length(delta_vals))
high_epsilon_Irb_R_end_extinction <- vector(length = length(delta_vals))

for (d_val in 1:length(delta_vals)) {
  delta <- delta_vals[d_val]
  for (e_val in 1:length(epsilon_vals)) {
    
    epsilon <- epsilon_vals[e_val]

    stoc_aucs <- stochastic_aucs(tf)
    determ_aucs <- deterministic_aucs(tf)
    
    #differentiate between AUC at tf and at each season
    stoc_tf_aucs <- stoc_aucs[[1]]
    determ_tf_aucs <- determ_aucs[[1]]
    
    stoc_epsilon_aucs[e_val] <- stoc_tf_aucs[5] #total AUC
    determ_epsilon_aucs[e_val] <- determ_tf_aucs[5]
  }
  
  #for each delta val, read 2x csv files (one per epsilon val) and add extinction props to vector
  #odd numbers = high epsilon vals, even numbers = low epsilon vals
  low_epsilon_simulations <- read.csv(file = paste0("values_csvs/simulations_output_num=",simulations_called,".csv"))
  high_epsilon_simulations <- read.csv(file = paste0("values_csvs/simulations_output_num=",simulations_called-1,".csv"))
  low_epsilon_extinction_props[d_val] <- low_epsilon_simulations[1,ncol(low_epsilon_simulations)]
  high_epsilon_extinction_props[d_val] <- high_epsilon_simulations[1,ncol(low_epsilon_simulations)]
  
  #same as above but for average num extinctions per run per group
  low_epsilon_group_extinction <- read.csv(file = paste0("values_csvs/average_extinctions_num=",simulations_called,".csv"))
  high_epsilon_group_extinction <- read.csv(file = paste0("values_csvs/average_extinctions_num=",simulations_called-1,".csv"))
  
  low_epsilon_Iwt_S_extinction[d_val] <- low_epsilon_group_extinction[1,2] #row=group, col2=average
  high_epsilon_Iwt_S_extinction[d_val] <- high_epsilon_group_extinction[1,2]
  low_epsilon_Irb_S_extinction[d_val] <- low_epsilon_group_extinction[2,2]
  high_epsilon_Irb_S_extinction[d_val] <- high_epsilon_group_extinction[2,2]
  low_epsilon_Irb_R_extinction[d_val] <- low_epsilon_group_extinction[3,2]
  high_epsilon_Irb_R_extinction[d_val] <- high_epsilon_group_extinction[3,2]
  
  # #same again for mean+perc values
  # low_epsilon_percs <- read.csv(file = paste0("values_csvs/means+percs_num=",simulation_means_called,".csv"))
  # high_epsilon_percs <- read.csv(file = paste0("values_csvs/means+percs_num=",simulation_means_called-1,".csv"))
  # #MAYBE ALL BELOW SHOULD BE DONE IN SIMULATION MEANS FUNCTION
  # #store 5th + 95th in vectors
  # #integrate per column for each group to find their auc
  # #plot aucs against delta
  
  #same again for group extinction props
  low_epsilon_group_end_extinction <- read.csv(file = paste0("values_csvs/group_end_extinctions_num=",simulations_called,".csv"))
  high_epsilon_group_end_extinction <- read.csv(file = paste0("values_csvs/group_end_extinctions_num=",simulations_called-1,".csv"))
  
  low_epsilon_Iwt_S_end_extinction[d_val] <- low_epsilon_group_end_extinction[1,2] #row=group, col2=average
  high_epsilon_Iwt_S_end_extinction[d_val] <- high_epsilon_group_end_extinction[1,2]
  low_epsilon_Irb_S_end_extinction[d_val] <- low_epsilon_group_end_extinction[2,2]
  high_epsilon_Irb_S_end_extinction[d_val] <- high_epsilon_group_end_extinction[2,2] #missed row 3 = Iwt_R
  low_epsilon_Irb_R_end_extinction[d_val] <- low_epsilon_group_end_extinction[4,2]
  high_epsilon_Irb_R_end_extinction[d_val] <- high_epsilon_group_end_extinction[4,2]
  
  #calculate difference in EI between epsilon vals
  stoc_epsilon_diff[d_val] <- stoc_epsilon_aucs[2] - stoc_epsilon_aucs[1] #low - high: effect of increasing spatial hereogenrity (ie. increasing coupling)
  determ_epsilon_diff[d_val] <- determ_epsilon_aucs[2] - determ_epsilon_aucs[1]
  high_stoc_epsilon_aucs[d_val] <- stoc_epsilon_aucs[1]
  low_stoc_epsilon_aucs[d_val] <- stoc_epsilon_aucs[2]
  high_determ_epsilon_aucs[d_val] <- determ_epsilon_aucs[1]
  low_determ_epsilon_aucs[d_val] <- determ_epsilon_aucs[2]
}

#plot difference between epsilon vals against delta vals, 2 lines: stochastic and deterministic
par(mar=c(5.1,4.1,4.1,2.1),
    mfrow=c(1,1),
    oma=c(0,0,0,0))
plot(delta_vals, stoc_epsilon_diff,
     ylim = c(min(stoc_epsilon_diff, determ_epsilon_diff), 
              max(stoc_epsilon_diff, determ_epsilon_diff)),
     xlab = expression(paste("Fitness cost of RB pathogen (",delta,")")),
     ylab = "Benefit of spatial heterogeneity",
     col = "blue",
     type = "l")
grid()
title("a)", adj = 0, line = 0.5)
lines(delta_vals, determ_epsilon_diff, col = "red")
legend("bottomright",
       legend = c("Deterministic", "Stochastic"),
       col = c("red", "blue"),
       bty = "n",
       lty = c(1,1))

###plot epidemic intensity against delta for low and high epsilon, stoc + determ
par(mfrow = c(1,#2),
              1),
    oma = c(0,0,0,0),
    mar = c(5.1,4.1,3.1,1.1))

#plot with y axis starting at 0
plot(delta_vals, high_stoc_epsilon_aucs,
     col = "blue",
     lwd = 1,
     lty = 1,
     type = "l",
     ylim = c(0, max(high_stoc_epsilon_aucs, low_stoc_epsilon_aucs,
                     high_determ_epsilon_aucs, low_determ_epsilon_aucs)),
     xlab = expression(paste("Fitness cost of RB pathogen (",delta,")")),
     ylab = "Epidemic intensity")
grid()
title("b)", adj = 0, line = 0.5)
lines(delta_vals, low_stoc_epsilon_aucs,
      col = "blue",
      lwd = 2,
      lty = 2)
lines(delta_vals, high_determ_epsilon_aucs,
      col = "red",
      lwd = 1,
      lty = 1)
lines(delta_vals, low_determ_epsilon_aucs,
      col = "red",
      lwd = 2,
      lty = 2)
legend("bottomleft",
       legend = c(expression(paste(epsilon," = 0.95 (stochastic)")),
                  expression(paste(epsilon," = 0.05 (stochastic)")),
                  expression(paste(epsilon," = 0.95 (deterministic)")),
                  expression(paste(epsilon," = 0.05 (deterministic)"))),
       col = c("blue", "blue", "red", "red"),
       lty = c(1,2),
       lwd = c(1,2),
       bty = "n")

#plot zoomed in
plot(delta_vals, high_stoc_epsilon_aucs,
     col = "blue",
     lwd = 1,
     lty = 1,
     type = "l",
     ylim = c(min(high_stoc_epsilon_aucs, low_stoc_epsilon_aucs, high_determ_epsilon_aucs, low_determ_epsilon_aucs),
              max(high_stoc_epsilon_aucs, low_stoc_epsilon_aucs, high_determ_epsilon_aucs, low_determ_epsilon_aucs)),
     xlab = expression(paste("Fitness cost of RB pathogen (",delta,")")),
     ylab = "Epidemic intensity")
grid()
title("c)", adj = 0, line = 0.5)
lines(delta_vals, low_stoc_epsilon_aucs,
      col = "blue",
      lwd = 2,
      lty = 2)
lines(delta_vals, high_determ_epsilon_aucs,
      col = "red",
      lwd = 1,
      lty = 1)
lines(delta_vals, low_determ_epsilon_aucs,
      col = "red",
      lwd = 2,
      lty = 2)

#plot proportion of extinctions against delta, high + low epsilon (stoc only) IN RB PATHOGEN ONLY
par(mfrow = c(1,1),
    oma = c(0,0,0,0),
    mar = c(5.1,4.1,4.1,2.1))

plot(delta_vals, low_epsilon_extinction_props,
     col = "blue",
     lty = 2,
     lwd = 2, #normal line width
     ylim = c(0, max(low_epsilon_extinction_props,
                     high_epsilon_extinction_props)),
     ylab = "Proportion of RB extinctions",
     xlab = expression(paste("Fitness cost of RB pathogen (",delta,")")),
     type = "l")
grid()
title("a)", adj = 0, line = 0.5)
lines(delta_vals, high_epsilon_extinction_props,
      col = "blue",
      lwd = 1,
      lty = 1)
legend("bottomright",
       legend = c(expression(paste(epsilon," = 0.95")),
                 expression(paste(epsilon," = 0.05"))),
       col = c("blue", "blue"),
       bty = "n",
       lty = c(1,2),
       lwd = c(1,2))

#plot proportion of extinctions for each group against delta
par(mfrow = c(1,1),
    oma = c(0,0,0,0),
    mar = c(5.1,4.1,3.1,1.1))

plot(delta_vals, low_epsilon_Iwt_S_end_extinction,
     col = "blue",
     lwd = 2,
     lty = 2,
     ylim = c(0, max(low_epsilon_Iwt_S_end_extinction, high_epsilon_Iwt_S_end_extinction,
                     low_epsilon_Irb_S_end_extinction, high_epsilon_Irb_S_end_extinction,
                     low_epsilon_Irb_R_end_extinction, high_epsilon_Irb_R_end_extinction)),
     ylab = "Proportion of extinctions",
     xlab = expression(paste("Fitness cost of RB pathogen (",delta,")")),
       #"",
     type = 'l')
grid()
title("a)", adj = 0, line = 0.5)
# mtext(expression(paste("Fitness cost of RB pathogen (",delta,")")),
#       side = 1, outer = TRUE)
lines(delta_vals, high_epsilon_Iwt_S_end_extinction,
      col = "blue",
      lwd = 1,
      lty = 1)
lines(delta_vals, low_epsilon_Irb_S_end_extinction,
      col = "forestgreen",
      lwd = 2,
      lty = 2)
lines(delta_vals, high_epsilon_Irb_S_end_extinction,
      col = "forestgreen",
      lwd = 1,
      lty = 1)
# lines(delta_vals, low_epsilon_Irb_R_end_extinction,
#       col = "orange",
#       lwd = 2,
#       lty = 2)
# lines(delta_vals, high_epsilon_Irb_R_end_extinction,
#       col = "orange",
#       lwd = 1,
#       lty = 2)

legend("bottomright", inset = c(0.05, 0.05),
       legend = c(expression(paste(epsilon," = 0.95 (WT)")),
                 expression(paste(epsilon," = 0.05 (WT)")),
                 expression(paste(epsilon," = 0.95 (RB)")),
                 expression(paste(epsilon," = 0.05 (RB)"))),
       col = c("blue", "blue", "forestgreen", "forestgreen"),
       bty = "n",
       lty = c(1,2),
       lwd = c(1,2))

#plot mean (average) number of extinctions per run for each group

par(mfrow = c(1,1),
    oma = c(1,0,0,0),
    mar = c(6.5,4.1,3.1,2.1))

plot(delta_vals, low_epsilon_Iwt_S_extinction,
     col = "blue",
     lwd = 2,
     lty = 2,
     ylim = c(0, max(low_epsilon_Iwt_S_extinction, high_epsilon_Iwt_S_extinction,
                     low_epsilon_Irb_S_extinction, high_epsilon_Irb_S_extinction,
                     low_epsilon_Irb_R_extinction, high_epsilon_Irb_R_extinction)),
     ylab = "Mean number of extinctions per run",
     xlab = expression(paste("Fitness cost of RB pathogen (",delta,")")),
       #"",
     type = "l")
grid()
title("b)", adj = 0, line = 0.5)
lines(delta_vals, high_epsilon_Iwt_S_extinction,
      col = "blue",
      lwd = 1,
      lty = 1)
lines(delta_vals, low_epsilon_Irb_S_extinction,
      col = "red",
      lwd = 2,
      lty = 2)
lines(delta_vals, high_epsilon_Irb_S_extinction,
      col = "red",
      lwd = 1,
      lty = 1)
lines(delta_vals, low_epsilon_Irb_R_extinction,
      col = "orange",
      lwd = 2,
      lty = 2)
lines(delta_vals, high_epsilon_Irb_R_extinction,
      col = "orange",
      lwd = 1,
      lty = 1)

#LEGEND - make new plot and overlay it
par(fig = c(0, 1, 0, 1), #coordinates of figure region (x1, x2, y1, y2)
    oma = c(0, 0, 0, 0), 
    mar = c(0, 0, 0, 0), 
    new = TRUE) #needed to ensure it's on existing plot, as fig creates a new plot by default
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottomleft",
       legend= c(expression(paste(epsilon," = 0.95 (WT infecting S)")),
                 expression(paste(epsilon," = 0.05 (WT infecting S)"))),
       col = c("blue", "blue"),
       lty = c(1,2),
       lwd = c(1,2),
       cex = 1,
       bty = "n",
       xpd = TRUE)  #allows drawing outside plot

legend("bottom",
       legend= c(expression(paste(epsilon," = 0.95 (RB infecting S)")),
                 expression(paste(epsilon," = 0.05 (RB infecting S)"))),
       col = c("red", "red"),
       lty = c(1,2),
       lwd = c(1,2),
       cex = 1,
       bty = "n",
       xpd = TRUE)  #allows drawing outside plot

legend("bottomright",
       legend= c(expression(paste(epsilon," = 0.95 (RB infecting R)")),
                 expression(paste(epsilon," = 0.05 (RB infecting R)"))),
       col = c("orange", "orange"),
       lty = c(1,2),
       lwd = c(1,2),
       cex = 1,
       bty = "n",
       xpd = TRUE)  #allows drawing outside plot