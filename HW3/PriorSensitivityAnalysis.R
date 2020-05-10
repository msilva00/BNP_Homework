# This is a script for the purposes of prior sensitivity analysis
# THIS IS POORLY WRITTEN 

library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/BNP_Homework/master/HW3/hwk3-data.txt"

### Posterior predictive loss
post_pred_loss = function(y, ypred, k = Inf){
  n = length(y)
  vars = apply(ypred, 2, var)
  means = apply(ypred, 2, mean)
  
  factor = k / (k + 1)
  if (k == Inf)
    factor = 1
  
  return (sum(vars) + factor*sum((y - means)^2))
}

y = read.csv(text = getURL(gitstring), header = F)[,1]

n = length(y)
ord = order(y)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(density(y), ylim = c(0,.2))

### Inverse gamma functions
# parameterization Mean = rate / (shape - 1)
dinvgamma = function(x, shape, rate, log = FALSE){
  out = shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x
  if (log)
    return (out)
  return (exp(out))
}
rinvgamma = function(n, shape, rate){
  1/rgamma(n, shape = shape, rate = rate)
}




post_means = function(prior_phi = c(3,15),
                      prior_tau = c(3,10),
                      prior_mu = c(0,3)){
  ### MCMC
  nburn = 500
  nmcmc = 1000
  
  # Parameter objects
  param_theta = matrix(0, nburn + nmcmc, n)
  param_alpha = double(nburn + nmcmc)
  param_mt = matrix(0, nburn + nmcmc, 2) # param_mt[,1] is mu, param_mt[,2] is tau^2
  param_phi = double(nburn + nmcmc)
  
  # Priors
  prior_phi_a   = prior_phi[1]   # Inverse Gamma (shape)
  prior_phi_b   = prior_phi[2]  # Inverse Gamma (rate)
  prior_alpha_a = 1   # Gamma (shape)           alpha controls n_star (discreteness of G)
  prior_alpha_b = 1   # Gamma (rate)
  prior_mu_a = prior_mu[1]   # Normal (mean)
  prior_mu_b  = prior_mu[2]   # Normal (variance)
  prior_tau2_a  = prior_tau[1]   # Inverse Gamma (shape)
  prior_tau2_b  = prior_tau[2]  # Inverse Gamma (rate)
  
  # Initial values
  #param_theta[1,] = y
  param_phi[1] = rinvgamma(1, prior_phi_a, prior_phi_b)
  param_alpha[1] = rgamma(1, prior_alpha_a, prior_alpha_b)
  param_mt[1,] = c(rnorm(1, prior_mu_a, sqrt(prior_mu_b)),
                   rinvgamma(1, prior_tau2_a, prior_tau2_b))
  
  # Iterations
  for (iter in 2:(nburn + nmcmc)){
    cat("\r", iter, "/", nburn + nmcmc)
    
    # Current parameters (used for convenience)
    curr_theta = param_theta[iter-1,]
    curr_alpha = param_alpha[iter-1]
    curr_mu = param_mt[iter-1,1]
    curr_tau2 = param_mt[iter-1,2]
    curr_phi = param_phi[iter-1]
    
    # Update thetas
    for (i in 1:n){
      n_j_minus = table(curr_theta[-i])
      theta_star_minus = as.numeric(names(n_j_minus))
      n_star_minus = length(theta_star_minus)
      
      # q0 which we derived with the integral of part 1.1
      q0 = dnorm(y[i], curr_mu, sqrt(curr_phi + curr_tau2))
      
      # also derived
      qj = dnorm(y[i], theta_star_minus, sqrt(curr_phi))
      
      # Probabilities of determining which to draw
      # See definitions in pdf
      A = curr_alpha * q0 / (curr_alpha * q0 + sum(n_j_minus * qj))
      
      # See definitions in pdf
      Bj = n_j_minus * qj / (curr_alpha * q0 + sum(n_j_minus * qj))
      
      # Make the update
      draw = sample(n_star_minus + 1, 1, prob = c(A, Bj))
      if (draw == 1){ # Make a draw from h
        curr_theta[i] = rnorm(1, (curr_mu*curr_phi + y[i]+curr_tau2)/(curr_phi + curr_tau2),
                              sqrt(curr_phi*curr_tau2/(curr_phi + curr_tau2)))
      } else { # Make a draw from the existing groups
        curr_theta[i] = theta_star_minus[draw-1]
      }
    }
    n_j = table(curr_theta)
    theta_star = as.numeric(names(n_j))
    n_star = length(theta_star)
    
    # Update alpha
    # Introduce latent variable eta to draw a new alpha (Escobar and West, 1995)
    eta = rbeta(1, curr_alpha + 1, n)
    epsilon = (prior_alpha_a + n_star - 1) /
      (n*(prior_alpha_b - log(eta)) + prior_alpha_a + n_star - 1)
    if (runif(1) < epsilon){
      curr_alpha = rgamma(1, prior_alpha_a + n_star, prior_alpha_b - log(eta))
    } else {
      curr_alpha = rgamma(1, prior_alpha_a + n_star - 1, prior_alpha_b - log(eta))
    }
    
    # Update param_mt: param_mt[,1] = mu, param_mt[,2] = tau^2
    S = sum(theta_star)
    curr_mu = rnorm(1, (prior_mu_a*curr_tau2 + prior_mu_b * S) /
                      (curr_tau2 + prior_mu_b * n_star),
                    sqrt(curr_tau2 * prior_mu_b / (curr_tau2 + prior_mu_b * n_star)))
    curr_tau2 = rinvgamma(1, prior_tau2_a + n_star/2,
                          prior_tau2_b + 1/2 * sum((theta_star - curr_mu)^2))
    
    # Update phi
    curr_phi = rinvgamma(1, prior_phi_a + n/2, prior_phi_b + 1/2 * sum((y - curr_theta)^2))
    
    # Update the clusters (Can be shown to be normally distributed)
    for (j in 1:n_star){
      temp_w = which(curr_theta == theta_star[j])
      curr_theta[temp_w] = rnorm(1,
                                 (curr_phi*curr_mu + curr_tau2*sum(y[temp_w])) / (curr_phi + n_j[j]*curr_tau2),
                                 sqrt(curr_phi*curr_tau2 / (curr_phi + n_j[j]*curr_tau2)))
    }
    
    # Save the current draws
    param_theta[iter,] = curr_theta
    param_alpha[iter] = curr_alpha
    param_mt[iter,1] = curr_mu
    param_mt[iter,2] = curr_tau2
    param_phi[iter] = curr_phi
    
    if (iter == (nburn + nmcmc))
      cat("\n")
  }
  
  ### Remove burn-in
  param_theta = tail(param_theta, nmcmc)
  param_alpha = tail(param_alpha, nmcmc)
  param_phi = tail(param_phi, nmcmc)
  param_mt = tail(param_mt, nmcmc) # param_mt[,1] = mu, param_mt[,2] = tausq
  
  post_means = list(post_mu = mean(param_mt[,1]), 
                    post_tau = mean(param_mt[,2]),
                    post_phi = mean(param_phi))
  post_samples = list(samples_mu = param_mt[,1],
                      samples_tau = param_mt[,2],
                      samples_phi = param_phi)
  # return(post_means)
  return(post_samples)
}
# prior_phi = c(3,15)
# prior_tau = c(3,10)
# prior_mu = c(0,3)
library(parallel)
multiple_tau = list(c(3,10),c(2,0.1), c(3,100))
# NOTE TO SELF DO NOT RUN THE NEXT LINE TWICE
# tests1 = mclapply(multiple_tau, function(x) post_means(prior_tau = x), mc.cores = 3)
par(mfrow = c(3,3), oma = c(2,6,6,2), mar = c(1.5,1.5,1.5,1.5))
for(i in 1:3){
  if(i == 1){
    main1 = expression("Posterior density for " ~ mu)
    main2 = expression("Posterior density for " ~ tau^2)
    main3 = expression("Posterior density for " ~ phi)
  }
  else{
    main1 = ""
    main2 = ""
    main3 = ""
  }
  plot(density(tests1[[i]]$samples_mu), main = main1, xlab = "", ylab = "", cex.main = 1.5)
  
  plot(density(tests1[[i]]$samples_tau), main = main2, xlab = "", ylab = "", cex.main = 1.5)
  plot(density(tests1[[i]]$samples_phi), main = main3, xlab = "", ylab = "", cex.main = 1.5)
}
title(main = expression("Different hyperparameters for"~ tau^2),
      outer = T, line = 4, cex.main = 2)
title(main = expression("Prior specifications: " ~ mu ~"~N(0,3),   " ~ phi ~ "~IG(3,15)"),
      outer = T, line = 2, cex.main = 2)
mtext(text = expression(tau^2 ~ "~ IG(3,100)                 " ~ 
                          tau^2 ~ "~ IG(2,0.1)                " ~
                          tau^2 ~ "~ IG(3,10)"), side=2, line=1.5, outer = T)


multiple_mu = list(c(0,3),c(10,3), c(10,0.1))
# NOTE TO SELF DO NOT RUN THE NEXT LINE TWICE
# tests2 = mclapply(multiple_mu, function(x) post_means(prior_mu = x), mc.cores = 3)
for(i in 1:3){
  if(i == 1){
    main1 = expression("Posterior density for " ~ mu)
    main2 = expression("Posterior density for " ~ tau^2)
    main3 = expression("Posterior density for " ~ phi)
  }
  else{
    main1 = ""
    main2 = ""
    main3 = ""
  }
  plot(density(tests2[[i]]$samples_mu), main = main1, xlab = "", ylab = "", cex.main = 1.5)
  plot(density(tests2[[i]]$samples_tau), main = main2, xlab = "", ylab = "", cex.main = 1.5)
  plot(density(tests2[[i]]$samples_phi), main = main3, xlab = "", ylab = "", cex.main = 1.5)
}
title(main = expression("Different hyperparameters for"~ mu),
      outer = T, line = 4, cex.main = 2)
title(expression("Prior specifications: " ~ tau^2 ~"~IG(3,10),   " ~ phi ~ "~IG(3,15)"), 
      outer = T, line = 2, cex.main = 2)
mtext(text = expression(mu ~ "~ N(10,0.1)                 " ~ 
                          mu ~ "~ N(10,3)                " ~
                          mu ~ "~ N(0,3)"), side=2, line=1.5, outer = T)


multiple_phi = list(c(3,15),c(3,0.1), c(3,100))
# NOTE TO SELF DO NOT RUN THE NEXT LINE TWICE
# tests3 = mclapply(multiple_phi, function(x) post_means(prior_phi = x), mc.cores = 3)
for(i in 1:3){
  if(i == 1){
    main1 = expression("Posterior density for " ~ mu)
    main2 = expression("Posterior density for " ~ tau^2)
    main3 = expression("Posterior density for " ~ phi)
  }
  else{
    main1 = ""
    main2 = ""
    main3 = ""
  }
  plot(density(tests3[[i]]$samples_mu), main = main1, xlab = "", ylab = "", cex.main = 1.5)
  plot(density(tests3[[i]]$samples_tau), main = main2, xlab = "", ylab = "", cex.main = 1.5)
  plot(density(tests3[[i]]$samples_phi), main = main3, xlab = "", ylab = "", cex.main = 1.5)
}
title(main = expression("Different hyperparameters for"~ phi),
      outer = T, line = 4, cex.main = 2) 
title(expression("Prior specifications: " ~ tau^2 ~"~IG(3,10),   " ~ mu ~ "~N(0,3)"), 
      outer = T, line = 2, cex.main = 2)
mtext(text = expression(phi ~ "~ IG(3,100)                 " ~ 
                          phi ~ "~ IG(3,0.1)                " ~
                          phi ~ "~ IG(3,15)"), side=2, line=1.5, outer = T)
