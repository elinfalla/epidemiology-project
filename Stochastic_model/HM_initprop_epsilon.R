#IMPORTANT: RUN 2_PATCH_STOCHAST FIRST
require(ggplot2)

par(mfrow = c(1,1),
    oma = c(0,0,0,0),
    mar = c(5.1,4.1,4.1,2.1))

#range of initial rb proportions to be used
init.prop.infected <- wt_eq_prop #assumes only S field are initially infected
init.num.infected <- init.prop.infected * N
initIrb_proportions <- seq(3e-6, 3e-5, length.out = 7)

#range of epsilon values to be used
epsilon_values <- seq(0,1, length.out = 7)

for(prop in 1:length(initIrb_proportions)) {
  
  #calculate initial numbers in each group
  initIrb_S <- round(init.num.infected * initIrb_proportions[prop])
  initIwt_S <- init.num.infected - initIrb_S
  
  initS_S <- (1 - phi)*N - initIwt_S - initIrb_S
  initS_R <- phi*N - initIwt_R - initIrb_R
  
  for(e_val in 1:length(epsilon_values)) {
    
    epsilon <- epsilon_values[e_val]
    
    stoc_aucs <- stochastic_aucs(tf)
    determ_aucs <- deterministic_aucs(tf)
    
    stoc_tf_aucs <- stoc_aucs[[1]]
    determ_tf_aucs <- determ_aucs[[1]]
    
    stoc_total_auc <- stoc_tf_aucs[5] #5 = total auc
    determ_total_auc <- determ_tf_aucs[5]
    
    #calculate proportional difference
    sd_Pdiff_total <- determ_tf_aucs[5] / stoc_tf_aucs[5]
    
    if(prop == 1 & e_val == 1) {
      heatmap_df <- data.frame(initIrb_proportions[prop], epsilon, sd_Pdiff_total)
      colnames(heatmap_df) <- c("initIrb_prop", "Epsilon", "prop_diff")
    }
    else {
      heatmap_df <- rbind(heatmap_df, c(initIrb_proportions[prop], epsilon, sd_Pdiff_total))
    }
  }
}

plot_heatmap <- ggplot(data = heatmap_df, 
                       mapping = aes(x = initIrb_prop,
                                     y = Epsilon,
                                     fill = prop_diff)) +
  geom_tile() +
  scale_fill_continuous(limits = c(1, 1.4), breaks = c(1, 1.1, 1.2, 1.3, 1.4)) +
  ylab(label = expression(paste("Degree of coupling (",epsilon,")"))) +
  xlab(label = expression(paste("Initial proportion of RB pathogen (",gamma,")"))) +
  labs(fill = "Proportional \n difference")

plot_heatmap
