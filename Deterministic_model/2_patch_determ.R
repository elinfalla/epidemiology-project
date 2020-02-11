# delete everything
rm(list=ls(all=TRUE))

# import library for implementing ode
require(deSolve)

#PARAMETERS
N <- 10^4*100 #total number of plants

epsilon <- 0.5 #degree of coupling; epsilon = 1 means totally mixed
beta <- 9e-6 # infection rate * when no resistant crops, equilibrium is 50% get infected
rho <- 1 #replanting rate; rho = 1 means planted in seasons *dont bother changing
delta <- 0.35 #fitness cost of rb *
gamma <- 0 #completeness of resistance; gamma = 0 means fully resistant *leave at zero
phi <- 0.5 #cropping ratio; proportion of resistant hosts

#TIMEFRAME
tf <- 40
time_points <- seq(0, tf, length.out = 1001)

#INITIAL STATES
initIwt_S <- round(0.33*N - (0.33*N*4e-4))
initIrb_S <- round(0.33*N*4e-4)
initIwt_R <- 0
initIrb_R <- 0

#2 PATCH ODE FUNCTION
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
  #S_S <- states["S_S"]
  #S_R <- states["S_R"]
  
  #dS_S <- rho * N - beta * S_S * (p * ((1 - delta) * Irb_S + Iwt_S) + (1 - p)*((1-delta)*Irb_R +Iwt_R)) - rho*S_S
  S_R <- phi * N - Iwt_R - Irb_R
  S_S <- (1-phi) * N - Iwt_S - Irb_S
  
  dIwt_S <- beta * S_S * (p * Iwt_S + (1 - p) * Iwt_R) - rho * Iwt_S
  dIrb_S <- beta * S_S * (1 - delta) * (p * Irb_S + (1 - p) * Irb_R) - rho * Irb_S
  
  #dS_R <- rho * N - beta * S_R * (p*((1-delta)*Irb_R + gamma*Iwt_R) +(1-p)*((1-delta)*Irb_S+gamma*Iwt_S)) - rho*S_R
  dIwt_R <- beta * S_R * gamma * (p * Iwt_R + (1 - p) * Iwt_S) - rho * Iwt_R
  dIrb_R <- beta * S_R * (1 - delta) * (p * Irb_R + (1 - p) * Irb_S) - rho * Irb_R
  
  return(list(c(dIwt_S, dIrb_S, dIwt_R, dIrb_R#, dS_S, dS_R
                )))
}

trajectory <- ode(y = c(Iwt_S = initIwt_S, Irb_S = initIrb_S, Iwt_R = initIwt_R, Irb_R = initIrb_R#, S_S =800, S_R=800
                        ),
                  time = time_points,
                  parms = c(N = N, epsilon = epsilon, beta = beta, rho = rho, delta = delta, gamma = gamma, phi = phi),
                  func = SI_ode)

trajectory_df <- data.frame(trajectory)

#PLOT NUM INFECTED IN EACH GROUP OVER TIME
par(mar = c(5.1, 4.1, 4.1, 10.5)) #margins, default is 5.1(bot),4.1(l),4.1(top),2.1(r)
plot(x = trajectory_df$time, y = trajectory_df$Irb_R,
     type = "l",
     col = "orange",
     xlab = "Time",
     ylab = "Number infected",
     main = "Numbers of infected in each group - deterministic",
     bty = "l", 
     ylim=c(0,N/2))
lines(x = trajectory_df$time, y = trajectory_df$Irb_S,
      col = "red")
lines(x = trajectory_df$time, y = trajectory_df$Iwt_R,
      col = "green")
lines(x = trajectory_df$time, y = trajectory_df$Iwt_S,
      col = "blue")

#LEGEND
legend("bottomright", inset = c(-0.6,0),
       legend= c("Susceptible by wt", "Susceptible by rb", "Resistant by wt", "Resistant by rb"),
       col = c("blue", "red", "green", "orange"),
       xpd = TRUE,   #allows drawing outside plot
       bty = "n",
       lty = c(1,1))