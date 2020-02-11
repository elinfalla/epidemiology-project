#all parameters are as specified in 2_patch_stochast.R

#IMPORTANT: RUN 2_PATCH_STOCHAST FIRST

#dev.off() #resets par but clears all previous graphs
par(mar = c(5.1, 4.1, 4.1, 2.1)) #margins, default is 5.1(bot),4.1(l),4.1(top),2.1(r)

N <- 10^4*100 #total number of plants
init.prop.infected <- wt_eq_prop #assumes only S field are initially infected
init.num.infected <- init.prop.infected * N

initIrb_proportions <- seq(3e-6, 4.5e-5, length.out = 25 #starts at 1 individual infected by rb: 1/init.num.infected
                           #1e-3, 1e-1, length.out = 5
                           )
scale_factor <- 10^5
initIrb_proportions_graph <- initIrb_proportions*scale_factor

seasons <- seq(1, tf, length.out = 3)

#create vectors for total tf aucs for stoc and determ
stoc_total_auc <- vector(length = length(initIrb_proportions))
determ_total_auc <- vector(length = length(initIrb_proportions))

stoc_Iwt_S_auc <- vector(length = length(initIrb_proportions))
stoc_Irb_S_auc <- vector(length = length(initIrb_proportions))
stoc_Iwt_R_auc <- vector(length = length(initIrb_proportions))
stoc_Irb_R_auc <- vector(length = length(initIrb_proportions))

determ_Iwt_S_auc <- vector(length = length(initIrb_proportions))
determ_Irb_S_auc <- vector(length = length(initIrb_proportions))
determ_Iwt_R_auc <- vector(length = length(initIrb_proportions))
determ_Irb_R_auc <- vector(length = length(initIrb_proportions))

#create vectors for differences between stoc and determ for each group at each init prop
sd_diff_Iwt_S <- vector(length = length(initIrb_proportions))
sd_diff_Irb_S <- vector(length = length(initIrb_proportions))
sd_diff_Iwt_R <- vector(length = length(initIrb_proportions))
sd_diff_Irb_R <- vector(length = length(initIrb_proportions))

#create vector for difference between total AUC at each init prop
sd_Pdiff_total <- vector(length = length(initIrb_proportions))
sd_Adiff_total <- vector(length = length(initIrb_proportions))

  
col_count <- 1 #to count through colours for seasons vs difference graph
#par(mfrow = c(1,2)) = for when want to plot trajectories graphs: note can't plot num seasons graph at same time

for (prop in 1:length(initIrb_proportions)) {
  
  #calculate initial numbers in each group
  initIrb_S <- init.num.infected * initIrb_proportions[prop]
  initIwt_S <- init.num.infected - initIrb_S
  
  initS_S <- (1 - phi)*N - initIwt_S - initIrb_S
  initS_R <- phi*N - initIwt_R - initIrb_R
  
  #call functions to find AUCs
  stoc_aucs <- stochastic_aucs(tf)
  determ_aucs <- deterministic_aucs(tf)
  
  #differentiate between AUC at tf and at each season
  stoc_tf_aucs <- stoc_aucs[[1]]
  determ_tf_aucs <- determ_aucs[[1]]
  
  stoc_season_aucs <- stoc_aucs[[2]]
  determ_season_aucs <- determ_aucs[[2]]
  
  #AUC AT TF
  stoc_total_auc[prop] <- stoc_tf_aucs[5] #5 = total auc
  determ_total_auc[prop] <- determ_tf_aucs[5]
  
  stoc_Iwt_S_auc[prop] <- stoc_tf_aucs[1]
  stoc_Irb_S_auc[prop] <- stoc_tf_aucs[2]
  stoc_Iwt_R_auc[prop] <- stoc_tf_aucs[3]
  stoc_Irb_R_auc[prop] <- stoc_tf_aucs[4]
  determ_Iwt_S_auc[prop] <- determ_tf_aucs[1]
  determ_Irb_S_auc[prop] <- determ_tf_aucs[2]
  determ_Iwt_R_auc[prop] <- determ_tf_aucs[3]
  determ_Irb_R_auc[prop] <- determ_tf_aucs[4]
  
  #find difference between stoc + determ for total AUC
  sd_Pdiff_total[prop] <- determ_tf_aucs[5] / stoc_tf_aucs[5]
  sd_Adiff_total[prop] <- determ_tf_aucs[5] - stoc_tf_aucs[5]
  
  #find difference between stochastic and determ auc for each group
  sd_diff_Iwt_S[prop] <- determ_tf_aucs[1] - stoc_tf_aucs[1]
  sd_diff_Irb_S[prop] <- determ_tf_aucs[2] - stoc_tf_aucs[2]
  sd_diff_Iwt_R[prop] <- determ_tf_aucs[3] - stoc_tf_aucs[3]
  sd_diff_Irb_R[prop] <- determ_tf_aucs[4] - stoc_tf_aucs[4]
  
  # #run simulations at current init prop
  # if (prop == 1 || prop == 10) {
  #   all_simulations <- simulations(tf)
  #   Irb_S_simulations <- all_simulations[["Irb_S"]]
  #   Irb_R_simulations <- all_simulations[["Irb_R"]]
  #   Irb_total_simulations <- Irb_S_simulations[,-1] + Irb_R_simulations[,-1] #take out column for time
  # 
  #   #PLOT FOR SIMULATION TRAJECTORY GRAPH
  #   #par(mfrow = c(1,1))
  #   plot(0,0,
  #      xlim = c(min(time_points), max(time_points)),
  #      ylim = c(0, max(Irb_total_simulations)),
  #      type = "n",
  #      xlab = "Time",
  #      ylab = "Proportion infected by RB pathogen")
  # 
  #   legend("topleft",
  #          paste("Initial rb proportion:", initIrb_proportions[prop]),
  #          bty = "n")
  # 
  #   for(column in 1:ncol(Irb_total_simulations)) {
  #     lines(time_points, Irb_total_simulations[,column])
  #   }
  # }

  # #AUC FOR SEASONS
  # #if prop is one of desired values:
  # if (prop == 1 || prop == 2 || prop == 3 || prop == 5 || prop == 10) {
  # 
  #   #find difference between stoc + determ at each season
  #   sd_diff_seasons <- vector(length = length(stoc_season_aucs))
  #   for (season in 1:length(stoc_season_aucs)) {
  #     sd_diff_seasons[season] <- determ_season_aucs[season] - stoc_season_aucs[season]
  #   }
  # 
  #   if (prop == 1) {
  #     layout(matrix(1:2, ncol = 2), widths = c(0.8, 0.2))
  #     par(mar=c(4.1,4.1,3.1,0.5),
  #         oma = c(0,0,0,0))
  #     colfunc <- colorRampPalette(c("red", "orange", "yellow", #"green",
  #                                   "blue"))
  # 
  #     #plot empty graph for seasons run against diff between stoc + determ
  #     plot(0,0,
  #          xlim = c(0, max(seasons)),
  #          ylim = c(0, max(sd_diff_seasons)), #maybe 0 wrong??
  #          type = "n",
  #          xlab = "Number of seasons",
  #          ylab = "Difference between stochastic and deterministic")
  #     grid()
  # 
  #   }
  #   cols <- c("blue", "green", "yellow", "orange", "red")
  #   lines(seasons, sd_diff_seasons, col = cols[col_count]) #THIS MAY NEED TO CHANGE LATER
  #   col_count <- col_count + 1
  # }
}

# par(mar=c(4.1,0.1,3.1,0.5))
# legend_image <- as.raster(matrix(colfunc(20), ncol=1))
# plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "")
# text(x=1.5, y = seq(0,1,l=2), labels = seq(min(initIrb_proportions),max(initIrb_proportions),l=2))
# rasterImage(legend_image, 0, 0, 1,1)
# mtext("Initial proportion \n of RB pathogen", cex = 0.8)

#PLOT TOTAL AUC (EPIDEMIC INTENSITY) AGAINST INIT PROP
par(mar=c(2.1,4.1,3.1,1.1), 
    mfrow=c(1,2),
    oma=c(2,0,0,0))

#yaxis from 0
plot(initIrb_proportions_graph, determ_total_auc, type = "l", col = "red",
     xlab = #expression(paste("Initial proportion of RB pathogen (",gamma,")")),
       "",
     ylab = "Epidemic intensity",
     ylim = c(0,
             #min(determ_total_auc, stoc_total_auc), 
             max(determ_total_auc, stoc_total_auc)),
     xlim = c(0,5)
     )
grid()
title("a)", adj = 0, line = 0.5)
mtext(expression(paste("Initial proportion of RB pathogen (",gamma,") (x",10^{-5},")")),
      side = 1, outer = TRUE)
lines(initIrb_proportions_graph, stoc_total_auc, col = "blue")
legend("bottomright", 
       legend = c("Deterministic", "Stochastic"),
       col = c("red", "blue"),
       bty = "n",
       lty = c(1,1))

#zoomed
plot(initIrb_proportions_graph, determ_total_auc, type = "l", col = "red",
     xlab = #expression(paste("Initial proportion of RB pathogen (",gamma,")")),
       "",
     ylab = #"Epidemic intensity",
       "",
     ylim = c(#0,
              min(determ_total_auc, stoc_total_auc), 
              max(determ_total_auc, stoc_total_auc)),
     xlim = c(0,5)
)
grid()
title("b)", adj = 0, line = 0.5)
#mtext(expression(paste("Initial proportion of RB pathogen (",gamma,")")), side = 1, outer = TRUE)
lines(initIrb_proportions_graph, stoc_total_auc, col = "blue")


#PLOT PROP + ABS DIFFERENCE BETWEEN STOC+DETERM FOR ABOVE GRAPH
par(mar=c(2.1,4.1,3.1,1.1),
    mfrow=c(1,2),
    oma=c(2,0,0,0))
plot(initIrb_proportions_graph, sd_Pdiff_total, type = "l",
     xlab = #expression(paste("Initial proportion of RB pathogen (",gamma,")")),
       "",
     ylab = "Proportional Difference",
     #cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2
     xlim = c(0,5)
     )
grid()
title("c)", adj = 0, line = 0.5#, cex.main=2
      )

plot(initIrb_proportions_graph, sd_Adiff_total, type = "l",
     xlab = expression(paste("Initial proportion of RB pathogen (",gamma,")")),
     ylab = "Absolute Difference",
     #cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2
     xlim = c(0,5)
     )
grid()
title("d)", adj = 0, line = 0.5#, cex.main=2
      )
mtext(expression(paste("Initial proportion of RB pathogen (",gamma,") (x",10^{-5},")")), 
      side = 1, outer = TRUE#, cex.main=2
      )

#PLOT EACH GROUP AUC AGAINST INIT PROP - DETERM
par(mar=c(4.1,4.7,3.1,2.1), 
    mfrow=c(1,3),
    oma=c(5.5,0,0,0)) #outer margin, surrounding all plots

plot(initIrb_proportions_graph, determ_Iwt_S_auc, type = "l", col = "blue",
     ylim = c(0,
       #min(determ_Iwt_S_auc, determ_Irb_S_auc, determ_Iwt_R_auc, determ_Irb_R_auc),
              max(determ_Iwt_S_auc, determ_Irb_S_auc, determ_Iwt_R_auc, determ_Irb_R_auc, 11)),
     xlim = c(0,5),
     xlab = "",
     ylab = "Epidemic intensity",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
grid()
title("a)", adj = 0, line = 0.5, cex.main = 2)
mtext(expression(paste("Initial proportion of RB pathogen (",gamma,") (x",10^{-5},")")), 
      side = 1, outer = TRUE, cex=1.5) #outer means in outer margin (oma)
lines(initIrb_proportions_graph, determ_Irb_S_auc, col = "red")
###lines(initIrb_proportions, determ_Iwt_R_auc, col = "green")
lines(initIrb_proportions_graph, determ_Irb_R_auc, col = "orange")

#PLOT EACH GROUP AUC AGAINST INIT PROP - STOC
plot(initIrb_proportions_graph, stoc_Iwt_S_auc, type = "l", col = "blue",
     ylim = c(0,
       #min(stoc_Iwt_S_auc, stoc_Irb_S_auc, stoc_Iwt_R_auc, stoc_Irb_R_auc),
              max(stoc_Iwt_S_auc, stoc_Irb_S_auc, stoc_Iwt_R_auc, stoc_Irb_R_auc, 11)),
     xlim = c(0,5),
     xlab = "",
     ylab = "Epidemic intensity",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
grid()
title("b)", adj = 0, line = 0.5, cex.main = 2)
lines(initIrb_proportions_graph, stoc_Irb_S_auc, col = "red")
###lines(initIrb_proportions, stoc_Iwt_R_auc, col = "green")
lines(initIrb_proportions_graph, stoc_Irb_R_auc, col = "orange")

#PLOT ABS DIFFERENCE BETWEEN STOC + DETERM FOR ABOVE GRAPH(S)
plot(0,0,
     ylim = c(min(sd_diff_Iwt_S, sd_diff_Iwt_R, sd_diff_Irb_S, sd_diff_Irb_R),
              max(sd_diff_Iwt_S, sd_diff_Iwt_R, sd_diff_Irb_S, sd_diff_Irb_R)),
     xlim = #c(initIrb_proportions[1], max(initIrb_proportions)),
       c(0,5),
     type = "n",
     xlab = "",
     ylab = "Difference between stochastic and deterministic",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
grid()
title("c)", adj = 0, line = 0.5, cex.main = 2)
lines(initIrb_proportions_graph, sd_diff_Iwt_S, col = "blue")
lines(initIrb_proportions_graph, sd_diff_Irb_S, col = "red")
###lines(initIrb_proportions, sd_diff_Iwt_R, col = "green")
lines(initIrb_proportions_graph, sd_diff_Irb_R, col = "orange")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",
       horiz = TRUE,
       legend= c("WT infecting S", "RB infecting S", #"R infected by wt", 
                 "RB infecting R"),
       col = c("blue", "red", #"green", 
               "orange"),
       xpd = TRUE,   #allows drawing outside plot
       bty = "n",
       lty = c(1,1),
       cex = 1.7) #font size

