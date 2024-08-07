 
#RABVITG for a STOCHASTIC Linear Quadratic Model
# Load necessary libraries
library(stats)
library(ggplot2)
library(cowplot)
library(pracma)
library(pspline)
library(gridExtra)
library(nleqslv)
# Define global variables
rho <- 0.1
q1 <- 4
q2 <- 4
r1<-1
r2<-1
a<-(2-0.5*rho)
#a<-2
b1<-1
b2<-1
dt <- 0.3
c<-0

xmin <- 0
xmax <- 1



################################################3

# Define the RBF interpolation function
rbf_interpolate <- function(X, y, s, xnew) {
  # Define the Gaussian RBF function
  gaussian_rbf <- function(r, s) {
    exp(-(r^2) / (2 * s^2))
  }
  
  # Calculate the distance matrix for the input points
  k <- nrow(X)
  N <- ncol(X)
  dist_matrix <- as.matrix(dist(X))
  
  # Compute the RBF matrix for the input points
  K <- gaussian_rbf(dist_matrix, s)
  
  # Solve for the weights
  weights <- solve(K, y)
  
  # Function to calculate the RBF for new points
  interpolate_point <- function(xnew_row) {
    dist_new <- sqrt(rowSums((t(t(X) - xnew_row))^2))
    rbf_new <- gaussian_rbf(dist_new, s)
    sum(rbf_new * weights)
  }
  
  # Apply the interpolation to each row in xnew
  interpolated_values <- apply(xnew, 1, interpolate_point)
  
  return(interpolated_values)
}




# Define the value function for player 1
valfun1 <- function(k, x01, u01, u02, x1, v0) {
  xnew1 <- x01 + dt * (a*x01+b1*k+b2*u02)+ sqrt(dt)*c*x01
  xnew1 <- matrix(pmax(xmin, pmin(xmax*1.20, xnew1)))
  xnew2 <- x01 + dt * (a*x01+b1*k+b2*u02)- sqrt(dt)*c*x01
  xnew2 <- matrix(pmax(xmin, pmin(xmax*1.20, xnew2)))
  
  
  X <- matrix(x1)
  # Define output values
  y <- matrix(v0[,1])
  # Define shape parameter
  s <- 0.5#Define new points for interpolation (m x N matrix)
  
  # Perform interpolation
  val1_interp1  <- rbf_interpolate(X, y, s, xnew1)
  val1_interp2  <- rbf_interpolate(X, y, s, xnew2)
  
  
  val1_interp<- mean(c( val1_interp1,val1_interp2))
  
  val1 <- dt * (q1 * x01**2 + r1*k**2) + (1 - dt * rho) * val1_interp
  return(val1)
}

# Define the value function for player 2
valfun2 <- function(k, x01, u01, u02, x1, v0) {
  xnew1 <- x01 + dt * (a*x01+b1*u01+b2*k)+ sqrt(dt)*c*x01
  xnew1 <- matrix(pmax(xmin, pmin(xmax*1.20, xnew1)))
  xnew2 <- x01 + dt * (a*x01+b1*u01+b2*k)- sqrt(dt)*c*x01
  xnew2 <- matrix(pmax(xmin, pmin(xmax*1.20, xnew2)))
  
  X <- matrix(x1)
  # Define output values
  y <- matrix(v0[,2])
  # Define shape parameter
  s <- 0.5#Define new points for interpolation (m x N matrix)
  
  # Perform interpolation
  val2_interp1  <- rbf_interpolate(X, y, s, xnew1)
  val2_interp2  <- rbf_interpolate(X, y, s, xnew2)
  
  
  
  
  val2_interp<- mean(c( val2_interp1,val2_interp2))
  
  
  val2 <- dt * (q2 * x01**2 + r2*k**2) + (1 - dt * rho) * val2_interp
  return(val2)
}


# Define finite difference function for first derivative
central_difference <- function(x, y) {
  n <- length(x)
  dydx <- numeric(n)
  
  # Forward difference for the first point
  dydx[1] <- (y[2] - y[1]) / (x[2] - x[1])
  
  # Central difference for the interior points
  for (i in 2:(n-1)) {
    dydx[i] <- (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
  }
  
  # Backward difference for the last point
  dydx[n] <- (y[n] - y[n-1]) / (x[n] - x[n-1])
  
  return(dydx)
}
# Define the RABVITG scheme
RABVITG2D <- function(N) {
  tol <- 0.0001
  dif <- tol + 0.1
  its <- 1
  difference<-0
  # Initialize variables
  x1 <- seq(xmin, xmax, length.out = N)
  uold1 <- rep(0, N)
  uold2 <- rep(0, N)
  v0 <- matrix(0, nrow = N, ncol = 2)
  v0[, 1] <- dt * (q1 * x1**2 + r1*uold1**2)
  v0[, 2] <- dt * (q2 * x1**2 + r2*uold2**2)
  
  norma <- 2
  its <- 2
  itmax <- 4000
  tol <- 0.0001
  k11 <- numeric(N)
  k22 <- numeric(N)
  v1 <- matrix(0, nrow = N, ncol = 2, byrow = TRUE)
  C<-0
  M2<-0
  L1<-0
  L2<-0
  
  while (dif > tol) {
    for (i in 1:N) {
      x01 <- x1[i]
      u01 <- uold1[i]
      u02 <- uold2[i]
      
      itsg <- 1
      tolg <- 0.00001
      errorg <- tolg + 1
      while (errorg > tolg & itsg < 25) {
        itsg <- itsg + 1
        k1 <- optimize(valfun1, interval = c(-100, 100), x01 = x01, u01 = u01, u02 = u02, x1 = x1, v0 = v0)$minimum
        k2 <- optimize(valfun2, interval = c(-100, 100), x01 = x01, u01 = u01, u02 = u02, x1 = x1, v0 = v0)$minimum
        policystar <- c(k1, k2)
        policyold <- c(u01, u02)
        errorg[itsg] <- max(abs(policystar - policyold))
        errorg <- tail(errorg, 1)
        u01 <- 0.95 * k1 + 0.05 * u01
        u02 <- 0.95 * k2 + 0.05 * u02
      }
      
      k11[i] <- u01
      k22[i] <- u02
      
      v1[i, 1] <- (valfun1(u01, x01, u01, u02, x1, v0))
      v1[i, 2] <- (valfun2(u02, x01, u01, u02, x1, v0))
      
    }
    
    dif <- max(abs(v1 - v0))
    difference[its]<- dif
    dif <- tail(dif, 1)
    
    ################### Fin the optimal controls per player analitytically##########3
    
    # Define the function to solve
    myfun <- function(x,a,b1,b2,c,q1,q2,r1,r2) {
      
      F <- numeric(4)
      F[1] <- 2 * x[1] * (a + b2 * x[4]) + q1 + c^2 * x[1] + x[3] * b1 * x[1]
      F[2] <- 2 * x[2] * (a + b1 * x[3]) + q2 + c^2 * x[2] + x[4] * b2 * x[2]
      F[3] <- x[3] + (b1 / r1) * x[1]
      F[4] <- x[4] + (b2 / r2) * x[2]
      
      return(F)
    }
    
    # Initial guess
    x0 <- c(1, 1, -1.5, -1.5)
    
    # Solve the system of equations
    solution_NASH <- nleqslv(x0, myfun, control = list(trace = 0),a=a,b1=b1,b2=b2,c=c,q1=q1,q2=q2,r1=r1,r2=r2)
    ##########################################################################
    

    
    # Apply finite difference method to your data
    deriv1st <- central_difference(x1, v1[,1])
    deriv2st <- central_difference(x1, v1[,2])
    
    fderiv1st<-diff(deriv1st)/2   #the inverse of f'_u is y/2
    dvderiv1st<-diff(deriv1st)
    C[its]<-max(fderiv1st/dvderiv1st)
    dx1<-diff(x1)
    M2[its]<-max(dvderiv1st/dx1)
    L1[its]<-abs(max(diff(k11))/max(diff(uold1)))
    L2[its]<-abs(max(diff(k22))/max(diff(uold2)))
  
    #derivada segunda
    deriv21st <- central_difference(x1, deriv1st)
    deriv22st <- central_difference(x1, deriv2st)

    
    
    ##Plots while running##
    plot_data <- data.frame(x1 = x1, k11 = k11, k22 = k22, v1_player1 = v1[,1], v1_player2 = v1[,2],sol1=solution_NASH$x[3],sol2=solution_NASH$x[4])
    plot_difference<-data.frame(iter=seq(1,its),difference,L1,L2)
    # Control plots
      
    
    control_plot1 <- ggplot(plot_data, aes(x = x1)) +
      geom_line(aes(y = k11), col = "blue") +
      geom_point(aes(y = sol1*x1), col = "grey") +
      labs(x = "x", y = expression(u[1]), title = " Player 1") +
      theme_minimal()+
      theme(text = element_text(size = 10))
    
    control_plot2 <- ggplot(plot_data, aes(x = x1)) +
      geom_line(aes(y = k22), col = "red") +
      geom_point(aes(y = sol2*x1), col = "grey") +
      labs(x = "x", y = expression(u[2]), title = "Player 2") +
      theme_minimal()+
      theme(text = element_text(size = 10))
    dmax1<-max(plot_data$dv_player1)*0.9
    dmax2<-max(plot_data$dv_player2)*0.9
    # Derivative Value Function
   
    
    plotL1 <- ggplot(plot_difference, aes(x = iter, y = L1)) +
      geom_line(col = "blue") +
      labs(
        x = "Iteration",
        y = expression(L[1]),
        title = "L Player 1"
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 10, family = "Arial"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")
      )
    
    
    plotL2 <- ggplot(plot_difference, aes(x = iter, y = L2)) +
      geom_line(col = "red") +
      labs(
        x = "Iteration",
        y = expression(L[2]),
        title = "L Player 2"
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 10, family = "Arial"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")
      )    
    # Convergence plot
    convergence_plot <- ggplot(plot_difference, aes(x = iter, y = difference)) +
      geom_line(col = "black") +
      labs(
        x = "Iteration",
        y = expression(paste("||", v[it] - v[it-1], "||")),
        title = "Convergence"
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 10, family = "Arial"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")
      )
    # Value plot
    value_plot1 <- ggplot(plot_data, aes(x = x1)) +
      geom_line(aes(y = v1_player1), col = "blue") +
       labs(x = "x", y = "Value", title = " Player 1") +
      theme_minimal()+
      theme(text = element_text(size = 10))
    value_plot2 <- ggplot(plot_data, aes(x = x1)) +
      geom_line(aes(y = v1_player2), col = "red") +
        labs(x = "x", y = "Value", title = "Player 2") +
      theme_minimal()+
      theme(text = element_text(size = 10))
    
    
    top_row <- plot_grid(control_plot1, control_plot2, plotL1, plotL2, ncol = 4)
    bottom_row <- plot_grid(convergence_plot, NULL, value_plot1, value_plot2, ncol = 4, rel_widths = c(0.5, 0, 0.25, 0.25))
    final_plot <- plot_grid(top_row, bottom_row, nrow = 2)
    
    # Display the final plot
    print(final_plot)
    
    
    
    
    
    v0 <- v1
    uold1<-k11
    uold2<-k22
    
    its <- its + 1
  }
  
  return(list(v0 = v0, k11 = k11, k22 = k22,its=its,dif=difference,x1=x1))
}

# Call the main function
result <- RABVITG2D(8)








