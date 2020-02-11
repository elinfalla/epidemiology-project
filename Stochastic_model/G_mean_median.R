#IMPORTANT: RUN 2_PATCH_STOCHAST FIRST

require(matrixStats)

N <- 10^4*100 #total number of plants
init.prop.infected <- wt_eq_prop #assumes only S field are initially infected
init.num.infected <- init.prop.infected * N

initIrb_proportions <- seq(3e-6, 4.5e-5, length.out = 25) 
                        #IF THIS IS CHANGED ALSO CHANGE IF STATEMENT FOR TRAJECOTRY GRAPHS
                       #starts at 1 individual infected by rb: 1/init.num.infected
scale_factor <- 10^5
initIrb_proportions_graph = initIrb_proportions*scale_factor #so numbers on axes aren't too small

initprop_trajectories <- c(initIrb_proportions[1], 
                           initIrb_proportions[length(initIrb_proportions_graph)/3],
                           initIrb_proportions[length(initIrb_proportions_graph)])


#create vectors for total tf aucs for stoc and determ
stoc_total_auc <- vector(length = length(initIrb_proportions))
determ_total_auc <- vector(length = length(initIrb_proportions))

#create vectors for extinctions and simulations
Irb_total_simulations <- vector("list", length = 3) #3 because calculating trajectory at 3 points
rb_extinctions <- vector(length = length(initIrb_proportions))

for (prop in 1:length(initIrb_proportions)) {
  
  #calculate initial numbers in each group
  initIrb_S <- round(init.num.infected * initIrb_proportions[prop])
  initIwt_S <- init.num.infected - initIrb_S
  
  initS_S <- (1 - phi)*N - initIwt_S - initIrb_S
  initS_R <- phi*N - initIwt_R - initIrb_R
  
  #call functions to find AUCs
  stoc_aucs <- stochastic_aucs(tf)
  determ_aucs <- deterministic_aucs(tf)
  
  #differentiate between AUC at tf and at each season
  stoc_tf_aucs <- stoc_aucs[[1]]
  determ_tf_aucs <- determ_aucs[[1]]
  
  #total AUC at tf
  stoc_total_auc[prop] <- stoc_tf_aucs[5] #5 = total auc
  determ_total_auc[prop] <- determ_tf_aucs[5]
  
  initprop_simulations <- read.csv(file = paste0("values_csvs/simulations_output_num=",simulations_called,".csv"))
  rb_extinctions[prop] <- initprop_simulations[1, ncol(initprop_simulations)]
  
  #run trajectory at current init prop
  if (initIrb_proportions[prop] == initprop_trajectories[1] || 
      initIrb_proportions[prop] == initprop_trajectories[2]|| #CHECK RANGE OF VALUES OF INITIRB_PROPS + SET ACCORDINGLY
      initIrb_proportions[prop] == initprop_trajectories[3]) {
    
    if (initIrb_proportions[prop] == initprop_trajectories[1]) {
      count <- 1
    } 
    else if (initIrb_proportions[prop] == initprop_trajectories[2]) {
      count <- 2
    }
    else {
      count <- 3
    }
    
    initprop_simulations[1] <- NULL #all simulations excluding count column at start 
    initprop_simulations[ncol(initprop_simulations)] <- NULL #exclude extinctions calculation
    
    start_Irb_S_simulations <- (num_runs + 1) + 1 #= all Iwt_S cols (inc time) + 1 to get start of Irb_S
    end_Irb_S_simulations <- (num_runs + 1) + (num_runs + 1)
    
    Irb_S_simulations <- initprop_simulations[,start_Irb_S_simulations: 
                                           end_Irb_S_simulations] 
    
    start_Irb_R_simulations <- (num_runs + 1) + (num_runs + 1) + (num_runs + 1) + 1
    end_Irb_R_simulations <-  ncol(initprop_simulations)
    Irb_R_simulations <- initprop_simulations[,start_Irb_R_simulations: # = all Iwt_S, Irb_S, Iwt_R, + 1 = start of Irb_R
                                           end_Irb_R_simulations] 
    Irb_total_simulations[[count]] <- Irb_S_simulations[,-1] + Irb_R_simulations[,-1] #take out column for time
  }
}
 
#PLOT TOTAL AUC (EPIDEMIC INTENSITY) AGAINST INIT PROP WITH VERTICAL LINES FOR WHERE 
#TRAJECTORIES WILL BE TAKEN
par(mar=c(2.1,4.1,3.1,2.1), 
    mfrow=c(1,2),
    oma=c(2.6,0,0,0))
plot(initIrb_proportions_graph, determ_total_auc, type = "l", col = "red",
     xlab = #expression(paste("Initial proportion of RB pathogen (",gamma,")")),
       "",
     ylab = "Epidemic intensity",
     ylim = c(min(determ_total_auc, stoc_total_auc), max(determ_total_auc, stoc_total_auc)),
     xlim = c(0,5))
grid()
title("a)", adj = 0, line = 0.5)
mtext(expression(paste("Initial proportion of RB pathogen (",gamma,") (x",10^{-5},")")), 
      side = 1, outer = TRUE)
lines(initIrb_proportions_graph, stoc_total_auc, col = "blue")
legend("bottomright", inset = c(0.08, 0.01),
       legend = c("Deterministic", "Stochastic"),
       col = c("red", "blue"),
       bty = "n",
       lty = c(1,1))
#vertical lines to mark where trajectory graphs are taken from
initprop_trajectories_graphs <- initprop_trajectories*scale_factor
lines(x = c(initprop_trajectories_graphs[1], initprop_trajectories_graphs[1]), 
      y = c(0,max(determ_total_auc, stoc_total_auc)),
      lty = 2)
lines(x = c(initprop_trajectories_graphs[2], initprop_trajectories_graphs[2]),
      y = c(0,max(determ_total_auc, stoc_total_auc)),
      lty = 2)
lines(x = c(initprop_trajectories_graphs[3], initprop_trajectories_graphs[3]),
      y = c(0,max(determ_total_auc, stoc_total_auc)),
      lty = 2)

plot(initIrb_proportions_graph, rb_extinctions,
     type = "l",
     xlab = #expression(paste("Initial proportion of RB pathogen (",gamma,")")),
       "",
     ylab = "Proportion of extinctions",
     col = "blue")
grid()
title("b)", adj = 0, line = 0.5)
#vertical lines to mark where trajectory graphs are taken from
lines(x = c(initprop_trajectories_graphs[1], initprop_trajectories_graphs[1]), 
      y = c(0,max(determ_total_auc, stoc_total_auc)),
      lty = 2)
lines(x = c(initprop_trajectories_graphs[2], initprop_trajectories_graphs[2]),
      y = c(0,max(determ_total_auc, stoc_total_auc)),
      lty = 2)
lines(x = c(initprop_trajectories_graphs[3], initprop_trajectories_graphs[3]),
      y = c(0,max(determ_total_auc, stoc_total_auc)),
      lty = 2)

#PLOT FOR SIMULATION TRAJECTORY GRAPH
par(mar=c(5.1,4.6,3.1,2.1), 
    mfrow=c(1,3),
    oma=c(2.5,0,0,0))
trajectory_titles <- c("c)", "d)", "e)")
trajectory_extinctions <- c(rb_extinctions[1],
                            rb_extinctions[length(initIrb_proportions)/2],
                            rb_extinctions[length(initIrb_proportions)])
for(count in 1:3) {
  initprop_trajectory <- Irb_total_simulations[[count]]
  r_initprop_trajectories <- round(initprop_trajectories, digits = 7)
  plot(0,0,
       xlim = c(min(time_points), max(time_points)),
       ylim = c(0, max(initprop_trajectory)),
       type = "n",
       xlab = "Time (seasons)",
       ylab = "Proportion infected by RB pathogen",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  grid()
  title(trajectory_titles[count], adj = 0, line = 0.5, cex.main = 2)
  legend("topleft", inset = c(-0.05, 0),
         paste("Initial RB proportion:", as.expression(r_initprop_trajectories[count])),
         bty = "n",
         cex = 1.5)
  if(count == 3) {
    legend("bottomleft", inset = c(0.28, 0.05),
           legend = c("Median", "Mean", "Individual run"),
           col = c("blue", "red", "lightgreen"),
           lty = c(1,1),
           #paste("\n\nExtinctions proportion:", trajectory_extinctions[count]),
           bty = "n",
           cex = 1.5)
  } else {
    legend("bottomleft", inset = c(0.28, 0.05),
           paste("Extinctions proportion:", trajectory_extinctions[count]),
           bty = "n",
           cex = 1.5)
  }
  
  #plot mean and median of all trajectories
  mean_vals <- vector(length = nrow(initprop_trajectory))
  median_vals <-vector(length = nrow(initprop_trajectory))
  
  initprop_trajectory_m <- data.matrix(initprop_trajectory) #convert to matrix to calculate mean + median vals
  
  mean_vals <- rowMeans(initprop_trajectory_m) #rowMeans works for databases
  median_vals <- rowMedians(initprop_trajectory_m)

  #plot each trajectory individually
  for(column in 1:ncol(initprop_trajectory)) {
    lines(time_points, initprop_trajectory[,column],
          col = adjustcolor("lightgreen",alpha.f=0.6))
  }
  
  lines(time_points, mean_vals,
        col = "red")
  lines(time_points, median_vals,
        col = "blue")
}

# #legend under the plot - overlayed
# par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
# plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# legend("bottom",
#        horiz = TRUE,
#        legend= c("Median", "Mean", "Individual run"),
#        col = c("blue", "red", "lightgreen"),
#        xpd = TRUE,   #allows drawing outside plot
#        bty = "n",
#        lty = c(1,1),
#        cex = 1.7) #font size



