#IMPORTANT: RUN 2_PATCH_STOCHAST FIRST

delta_values <- seq(0, 1, length.out = 25)
seasons <- seq(1, tf, length.out = 5)

#create vectors for total tf aucs for stoc and determ
stoc_total_auc <- vector(length = length(delta_values))
determ_total_auc <- vector(length = length(delta_values))

stoc_Iwt_S_auc <- vector(length = length(delta_values))
stoc_Irb_S_auc <- vector(length = length(delta_values))
stoc_Iwt_R_auc <- vector(length = length(delta_values))
stoc_Irb_R_auc <- vector(length = length(delta_values))

determ_Iwt_S_auc <- vector(length = length(delta_values))
determ_Irb_S_auc <- vector(length = length(delta_values))
determ_Iwt_R_auc <- vector(length = length(delta_values))
determ_Irb_R_auc <- vector(length = length(delta_values))

#create vectors for differences between stoc and determ for each group at each init prop
sd_diff_Iwt_S <- vector(length = length(delta_values))
sd_diff_Irb_S <- vector(length = length(delta_values))
sd_diff_Iwt_R <- vector(length = length(delta_values))
sd_diff_Irb_R <- vector(length = length(delta_values))

#create vector for difference between total AUC at each init prop
sd_Pdiff_total <- vector(length = length(delta_values))
sd_Adiff_total <- vector(length = length(delta_values))


col_count <- 1 #to count through colours for seasons vs difference graph

for (d_val in 1:length(delta_values)) {
  
  delta <- delta_values[d_val]
  
  #call function to calculate AUCs
  stoc_aucs <- stochastic_aucs(tf)
  determ_aucs <- deterministic_aucs(tf)
  
  #differentiate between AUC at tf and at each season
  stoc_tf_aucs <- stoc_aucs[[1]]
  determ_tf_aucs <- determ_aucs[[1]]
  
  stoc_season_aucs <- stoc_aucs[[2]]
  determ_season_aucs <- determ_aucs[[2]]
  
  #AUC AT TF
  stoc_total_auc[d_val] <- stoc_tf_aucs[5] #5 = total auc
  determ_total_auc[d_val] <- determ_tf_aucs[5]
  
  stoc_Iwt_S_auc[d_val] <- stoc_tf_aucs[1]
  stoc_Irb_S_auc[d_val] <- stoc_tf_aucs[2]
  stoc_Iwt_R_auc[d_val] <- stoc_tf_aucs[3]
  stoc_Irb_R_auc[d_val] <- stoc_tf_aucs[4]
  determ_Iwt_S_auc[d_val] <- determ_tf_aucs[1]
  determ_Irb_S_auc[d_val] <- determ_tf_aucs[2]
  determ_Iwt_R_auc[d_val] <- determ_tf_aucs[3]
  determ_Irb_R_auc[d_val] <- determ_tf_aucs[4]
  
  #find difference between stoc + determ for total AUC
  sd_Pdiff_total[d_val] <- determ_tf_aucs[5] / stoc_tf_aucs[5]
  sd_Adiff_total[d_val] <- determ_tf_aucs[5] - stoc_tf_aucs[5]
  
  #find difference between stochastic and determ auc for each group
  sd_diff_Iwt_S[d_val] <- determ_tf_aucs[1] - stoc_tf_aucs[1]
  sd_diff_Irb_S[d_val] <- determ_tf_aucs[2] - stoc_tf_aucs[2]
  sd_diff_Iwt_R[d_val] <- determ_tf_aucs[3] - stoc_tf_aucs[3]
  sd_diff_Irb_R[d_val] <- determ_tf_aucs[4] - stoc_tf_aucs[4]
  
  #AUC FOR SEASONS
  #if delta is one of desired values
  if (d_val == 1 || d_val == 3 || d_val == 5 || d_val == 7 || d_val == 10) {

    #find difference between stoc + determ at each season
    sd_diff_seasons <- vector(length = length(stoc_season_aucs))
    for (season in 1:length(stoc_season_aucs)) {
      sd_diff_seasons[season] <- determ_season_aucs[season] - stoc_season_aucs[season]
    }

    if (d_val == 1) {
      layout(matrix(1:2, ncol = 2), widths = c(0.8, 0.2))
      par(mar=c(4.1,4.1,3.1,0.5),
          oma = c(0,0,0,0))
      colfunc <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"))

      #plot empty graph for seasons run against diff between stoc + determ
      plot(0,0,
           xlim = c(0, max(seasons)),
           ylim = c(0, max(sd_diff_seasons)), #maybe 0 wrong??
           type = "n",
           xlab = "Number of seasons",
           ylab = "Difference between stochastic and deterministic")
      grid()

    }
    cols <- c("blue", "green", "yellow", "orange", "red")
    lines(seasons, sd_diff_seasons, col = cols[col_count]) #THIS MAY NEED TO CHANGE LATER
    col_count <- col_count + 1
  }
}

par(mar=c(4.1,0.1,3.1,0.5))
legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "")
text(x=1.5, y = seq(0,1,l=2), labels = seq(min(delta_values),max(delta_values),l=2))
rasterImage(legend_image, 0, 0, 1,1)
mtext("Fitness cost \n of rb pathogen", cex = 0.8)

#PLOT TOTAL AUC (EPIDEMIC INTENSITY) AGAINST DELTA
par(mar=c(2.1,4.1,3.1,1.1),
    mfrow=c(1,2),
    oma=c(2,0,0,0))

#yaxis from 0
plot(delta_values, determ_total_auc, type = "l", col = "red",
     xlab = #expression(paste("Fitness cost of RB pathogen (",delta,")")),
       "",
     ylab = "Epidemic intensity",
     ylim = c(0,
       #min(determ_total_auc, stoc_total_auc), 
              max(determ_total_auc, stoc_total_auc)),
     xlim = c(0,1))
grid()
title("a)", adj = 0, line = 0.5)
lines(delta_values, stoc_total_auc, col = "blue")
legend("bottomright",
       legend = c("Deterministic", "Stochastic"),
       col = c("red", "blue"),
       bty = "n",
       lty = c(1,1))

#zoomed
plot(delta_values, determ_total_auc, type = "l", col = "red",
     xlab = #expression(paste("Fitness cost of RB pathogen (",delta,")")),
       "",
     ylab = #"Epidemic intensity",
       "",
     ylim = c(#0,
              min(determ_total_auc, stoc_total_auc), 
              max(determ_total_auc, stoc_total_auc)),
     xlim = c(0,1))
grid()
mtext(expression(paste("Fitness cost of RB pathogen (",delta,")")),
      side = 1, outer = TRUE)
title("b)", adj = 0, line = 0.5)
lines(delta_values, stoc_total_auc, col = "blue")

#PLOT PROP + ABS DIFFERENCE BETWEEN STOC+DETERM FOR ABOVE GRAPH
par(mar=c(2.1,4.1,3.1,2.1),
    mfrow=c(1,2),
    oma=c(2,0,0,0))
plot(delta_values, sd_Pdiff_total, type = "l",
     xlab = #expression(paste("Fitness cost of RB pathogen (",delta,")")),
       "",
     ylab = "Proportional Difference")
grid()
title("c)", adj = 0, line = 0.5)

plot(delta_values, sd_Adiff_total, type = "l",
     xlab = #expression(paste("Fitness cost of RB pathogen (",delta,")")),
       "",
     ylab = "Absolute Difference",
     xlim = c(0,1))
grid()
title("d)", adj = 0, line = 0.5)
mtext(expression(paste("Fitness cost of RB pathogen (",delta,")")),
      side = 1, outer = TRUE)


#PLOT EACH GROUP AUC AGAINST INIT PROP - DETERM
par(mar=c(4.1,4.6,3.1,2.1),
    mfrow=c(1,3),
    oma=c(5.5,0,0,0)) #outer margin, surrounding all plots
plot(delta_values, determ_Iwt_S_auc, type = "l", col = "blue",
     ylim = c(0,
       #min(determ_Iwt_S_auc, determ_Irb_S_auc, determ_Iwt_R_auc, determ_Irb_R_auc),
              max(determ_Iwt_S_auc, determ_Irb_S_auc, determ_Iwt_R_auc, determ_Irb_R_auc,
                  stoc_Iwt_S_auc, stoc_Irb_S_auc, stoc_Iwt_R_auc, stoc_Irb_R_auc)),
     xlim = c(0,1),
     xlab = #expression(paste("Fitness cost of RB pathogen (",delta,")")),
       "",
     ylab = "Epidemic intensity",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
grid()
title("a)", adj = 0, line = 0.5, cex.main = 2)
mtext(expression(paste("Fitness cost of RB pathogen (",delta,")")), 
      side = 1, outer = TRUE, cex = 1.5) #outer means in outer margin (oma)
lines(delta_values, determ_Irb_S_auc, col = "red")
###lines(delta_values, determ_Iwt_R_auc, col = "green")
lines(delta_values, determ_Irb_R_auc, col = "orange")


#PLOT EACH GROUP AUC AGAINST INIT PROP - STOC
plot(delta_values, stoc_Iwt_S_auc, type = "l", col = "blue",
     ylim = c(0,
       #min(stoc_Iwt_S_auc, stoc_Irb_S_auc, stoc_Iwt_R_auc, stoc_Irb_R_auc),
              max(stoc_Iwt_S_auc, stoc_Irb_S_auc, stoc_Iwt_R_auc, stoc_Irb_R_auc,
                  determ_Iwt_S_auc, determ_Irb_S_auc, determ_Iwt_R_auc, determ_Irb_R_auc)),
     xlim = c(0,1),
     xlab = #expression(paste("Fitness cost of RB pathogen (",delta,")")),
       "",
     ylab = "Epidemic intensity",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
grid()
title("b)", adj = 0, line = 0.5, cex.main = 2)
lines(delta_values, stoc_Irb_S_auc, col = "red")
###lines(delta_values, stoc_Iwt_R_auc, col = "green")
lines(delta_values, stoc_Irb_R_auc, col = "orange")


#PLOT DIFFERENCE BETWEEN STOC + DETERM FOR ABOVE GRAPH(S)
plot(0,0,
     ylim = c(min(sd_diff_Iwt_S, sd_diff_Irb_S, sd_diff_Irb_R),
              max(sd_diff_Iwt_S, sd_diff_Irb_S, sd_diff_Irb_R)),
     xlim = c(0,1),
     type = "n",
     xlab = #expression(paste("Fitness cost of RB pathogen (",delta,")")),
       "",
     ylab = "Difference between stochastic and deterministic",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
grid()
title("c)", adj = 0, line = 0.5, cex.main = 2)
lines(delta_values, sd_diff_Iwt_S, col = "blue")
lines(delta_values, sd_diff_Irb_S, col = "red")
###lines(delta_values, sd_diff_Iwt_R, col = "green")
lines(delta_values, sd_diff_Irb_R, col = "orange")

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
       cex = 1.7) #makes legend font 1.2x size of other font
