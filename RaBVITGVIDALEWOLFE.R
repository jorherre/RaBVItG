
#RABVITG for a STOCHASTIC Linear Quadratic Model
# Load necessary libraries
library(stats)
library(ggplot2)
library(cowplot)
library(pracma)
library(pspline)
library(gridExtra)
library(nleqslv)
library(akima)
library(interp)
# Define global variables
rho <- 0.1
m1 <- 1
c1 <- 0.1
m2<-1
c2<-0.1
r1<-1
#a<-2
r2<-1
b2<-1
delta<-0.5
dt <- 0.01 #no converge 0.30
sigma<-0.25

xmin <- 0.05
xmax <- 0.975

W<-rho+2*delta
R<-(r1^2)/(4*c1)
beta1<-(sqrt(W^2+12*R*m1)-W)/(6*R)
beta2=beta1
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
valfun1 <- function(k, x01,x02, u01, u02, x1,x2, v0) {

  
  # Compute new state values using the Euler-Maruyama method
  xnew1 <- x01 + dt * (r1*k*sqrt(1-x01)-r2*u02*sqrt(x01)-delta*(2*x01-1)) +sqrt(dt)*sigma*sqrt(x01*(1-x01))
  xnew1 <- matrix(pmax(xmin, pmin(xmax*1.05, xnew1)))
  xnew2 <- x01 + dt * (r1*k*sqrt(1-x01)-r2*u02*sqrt(x01)-delta*(2*x01-1)) -sqrt(dt)*sigma*sqrt(x01*(1-x01))
  xnew2 <- matrix(pmax(xmin, pmin(xmax*1.05, xnew2)))
  
  
  X <- matrix(x1)
  # Define output values
  y <- matrix(v0[,1])
  # Define shape parameter
  s <- 1#Define new points for interpolation (m x N matrix)
  
  # Perform interpolation
  val1_interp1  <- rbf_interpolate(X, y, s, xnew1)
  val1_interp2  <- rbf_interpolate(X, y, s, xnew2)
  
  
  val1_interp<- mean(c( val1_interp1,val1_interp2))
  

  
  # Compute the value function
  val1 <- (dt * (m1 * x01 - c1 * k^2) + (1 - dt * rho) * val1_interp)
  val1=-val1
  return(val1)
}



# Define the value function for player 2
valfun2 <- function(k, x01,x02, u01, u02, x1,x2, v0) {
  # Compute new state values using the Euler-Maruyama method
  xnew1 <- x01 + dt * (r1*u01*sqrt(1-x01)-r2*k*sqrt(x01)-delta*(2*x01-1))+sqrt(dt)*sigma*sqrt(x01*(1-x01))
  xnew1 <- matrix(pmax(xmin, pmin(xmax*1.05, xnew1)))
  xnew2 <- x01 + dt * (r1*u01*sqrt(1-x01)-r2*k*sqrt(x01)-delta*(2*x01-1)) -sqrt(dt)*sigma*sqrt(x01*(1-x01))
  xnew2 <- matrix(pmax(xmin, pmin(xmax*1.05, xnew2)))
  
  X <- matrix(x1)
  # Define output values
  y <- matrix(v0[,2])
  # Define shape parameter
  s <- 1#Define new points for interpolation (m x N matrix)
  
  # Perform interpolation
  val2_interp1  <- rbf_interpolate(X, y, s, xnew1)
  val2_interp2  <- rbf_interpolate(X, y, s, xnew2)
  
  
  val2_interp<- mean(c( val2_interp1,val2_interp2))
  
  
  
 
  
  
 
  
  # Compute the value function
  val2 <- (dt * (m2 * (1-x01) - c2 * k^2) + (1 - dt * rho) * val2_interp)
  val2=-val2
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
  x2 <- 1-x1
  uold1 <- 1+runif( N)
  uold2 <- 1+runif( N)
  v0 <- matrix(0, nrow = N, ncol = 2)
  v0[, 1] <- dt * (m1 * x1 - c1*uold1**2)
  v0[, 2] <- dt * (m2 * (1-x1) - c2*uold2**2)
  
  norma <- 2
  its <- 2
  itmax <- 4000
  tol <- 0.0001
  k11 <- numeric(N)
  k22 <- numeric(N)
  v1 <- matrix(0, nrow = N, ncol = 2, byrow = TRUE)
  
  while (dif > tol) {
    for (i in 1:N) {
      x01 <- x1[i]
      x02 <- x2[i]
      u01 <- uold1[i]
      u02 <- uold2[i]
      
      itsg <- 1
      tolg <- 0.00001
      errorg <- tolg + 1
      while (errorg > tolg & itsg < 25) {
        itsg <- itsg + 1
        k1 <- optimize(valfun1, interval = c(0, 100), x01 = x01, x02=x02,u01 = u01, u02 = u02, x1 = x1, x2=x2,v0 = v0)$minimum
        k2 <- optimize(valfun2, interval = c(0, 100), x01 = x01, x02=x02,u01 = u01, u02 = u02, x1 = x1, x2=x2,v0 = v0)$minimum
        policystar <- c(k1, k2)
        policyold <- c(u01, u02)
        errorg[itsg] <- max(abs(policystar - policyold))
        errorg <- tail(errorg, 1)
        u01 <- 0.95 * k1 + 0.05 * u01
        u02 <- 0.95 * k2 + 0.05 * u02
      }
      
      k11[i] <- u01
      k22[i] <- u02
      
      v1[i, 1] <- -(valfun1(u01, x01,x02, u01, u02, x1,x2, v0))
      v1[i, 2] <- -(valfun2(u02, x01,x02, u01, u02, x1,x2, v0))
      
    }
    
    dif <- max(abs(v1 - v0))
    difference[its]<- dif
    dif <- tail(dif, 1)
    
      ##########################################################################
    
    
    # Define finite difference function for first derivative
    finite_difference <- function(x, y) {
      # Ensure x and y are sorted by x
      ord <- order(x)
      x <- x[ord]
      y <- y[ord]
      
      # Compute differences
      dx <- diff(x)
      dy <- diff(y)
      
      # Compute derivative
      dydx <- dy / dx
      
      # Compute midpoints
      x_mid <- (x[-length(x)] + x[-1]) / 2
      
      return(data.frame(x_mid = x_mid, dydx = dydx))
    }
    
 
    
    # Apply finite difference method to your data
    deriv1st <- central_difference(x1, v1[,1])
    deriv2st <- central_difference(x1, v1[,2])
    
    ##Plots while running##
    plot_data <- data.frame(x1 = x1, k11 = k11, k22 = k22, v1_player1 = v1[,1], v1_player2 = v1[,2],dv_player1=deriv1st,dv_player2=deriv2st)
    plot_difference<-data.frame(iter=seq(1,its),difference)
    # Control plots
    
    
    control_plot1 <- ggplot(plot_data, aes(x = x1)) +
      geom_line(aes(y = k11), col = "blue") +
      geom_point(aes(y = (beta1*r1*sqrt(1-x1))/(2*c1)), col = "grey") +
      labs(x = "x", y = expression(u[1]), title = " Player 1") +
      theme_minimal()+
      theme(text = element_text(size = 10))
    
    control_plot2 <- ggplot(plot_data, aes(x = x1)) +
      geom_line(aes(y = k22), col = "red") +
      geom_point(aes(y = (beta2*r2*sqrt(x1))/(2*c2)), col = "grey") +
      
      labs(x = "x", y = expression(u[2]), title = "Player 2") +
      theme_minimal()+
      theme(text = element_text(size = 10))
    
    # Derivative Value Function
    derivative_1 <- ggplot(plot_data, aes(x = x1)) +
      geom_line(aes(y = dv_player1), col = "blue") +
      labs(x = "x", y = expression(v[1] * "'" * (x)), title = " Player 1") +
      theme_minimal()+
      theme(text = element_text(size = 10))
      
    derivative_2 <- ggplot(plot_data, aes(x = x1)) +
      geom_line(aes(y = dv_player2), col = "red") +
      labs(x = "x", y = expression(v[2] * "'" * (x)), title = "Player 2") +
      theme_minimal()+
      theme(text = element_text(size = 10))
    
    
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
    
    
    top_row <- plot_grid(control_plot1, control_plot2, derivative_1, derivative_2, ncol = 4)
    bottom_row <- plot_grid(convergence_plot, NULL, value_plot1, value_plot2, ncol = 4, rel_widths = c(0.5, 0, 0.25, 0.25))
    final_plot <- plot_grid(top_row, bottom_row, nrow = 2)
    
    # Display the final plot
    print(final_plot)
    
    
    
    
    
    
    v0 <- v1
    
    its <- its + 1
  }
  
  return(list(v0 = v0, k11 = k11, k22 = k22,its=its,dif=difference))
}

# Call the main function
result <- RABVITG2D(6)















# Access individual results from the list
v0_result <- result$v0
k11_result <- result$k11
k22_result <- result$k22
its_result<- result$its
dif_result<-result$dif
deriv1st <- predict(sm.spline(x1, v0_result[,1]), x1, 1) #V'(y)
deriv2st <- predict(sm.spline(x1, v0_result[,2]), x1, 1)





