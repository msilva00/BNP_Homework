# Libraries and Functions
library(MASS)
#### Trace, Hist, and Autocorrelation
trace_dens_acf= function(param, prob = 0.95, precision = 1000, main1="main1", main2="main2",
                         xlab2="xlab2", xlab1 = "Iteration", par_settings = T,
                         ylab1 = "ylab1", ylab2 = "Density", breaks = 20, main3 = "main3"){
  # param - the mcmc samples for a parameter
  # prob - high posterior density interval (i.e 95%)
  # precision - controls smoothness of the curve
  if(par_settings == T){
    par(mfrow = c(1,3), mar = c(4.1, 4.1, 3.1, 1.1), oma = c(1,0,1,0))
  }
  range = seq(0, 1-prob, length=precision)
  range = cbind(range, range+prob)
  best = range[which.min(apply(range, 1, function(y)
    diff(quantile(param, y)))),]
  quants = quantile(param, best)
  dens = density(param)
  # trace plot
  plot(param, type = 'l', main = main1,
       xlab = xlab1, ylab = ylab1)
  abline(h = mean(param), col = "forestgreen", lwd = 3)
  # histogram and density
  hist(param, xlab = xlab2, probability = TRUE,ylim = c(0,max(dens$y)), 
       main = main2)
  lines(density(param))
  x1 <- min(which(dens$x >= quants[1]))  
  x2 <- max(which(dens$x <=  quants[2]))
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha("lightgreen", 0.7)))
  abline(v = mean(param), col = "forestgreen", lwd = 3)
  
  acf(param, plot=T, main = main3)
}

autotune = function(accept, target = 0.25, k = 2.5){
  #' @description
  #' Used to adjust candidate sigmas for normal candidate densities
  #' Requires the calculation of acceptance rates within some given
  #' window of mcmc iterations. For example, every 500 draws compute
  #' the acceptance rate (0 <= A <= 1) for each parameter using the
  #' last 500 draws. Multiply each candidate sigma by autotune().
  #' 
  #' @param accept    the acceptance rate within a given window
  #' @param target    the acceptance rate goal
  #' @param k         the maximum the candidate sigma can change by, 1/k is 
  #'     minimum it can change by. For instance, if the current
  #'     sigma is 2, and in some window of mcmc iterations every
  #'     draw was accepted, the new candidate sigma would now
  #'     be 2*k, which should serve to reduce the acceptance rate.
  #'     On the other hand, if no draws are accepted, the new candidate
  #'     sigma would then be 2/k. I recommend k = window/50
  #'     Larger values of k will change sigma by a larger amount,
  #'     and vice versa for smaller values (k -> 1, autotune() -> 1,
  #'     everywhere)
  #' Courtesy of H.Ward, BYU
  (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
                                           ceiling(accept-target))-1))^sign(accept-target)
}

### Posterior predictive loss criterion
criterion = function(y, ypred, k = Inf){
  n = length(y)
  vars = apply(ypred, 2, var)
  means = apply(ypred, 2, mean)
  
  factor = k / (k + 1)
  if (k == Inf)
    factor = 1
  
  return (sum(vars) + factor*sum((y - means)^2))
}
dat = read.table("~/BNP_Homework/HW4/fabric.txt", header = TRUE)
x = dat$length / 100
#x = as.numeric(scale(dat$length))
y = dat$faults
n = length(y)
ord = order(y)
ord_x = order(x)

set.seed(1)
jity = jitter(y)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(x, y, pch = 20)


### Inverse gamma functions
# Mean = rate / (shape - 1)
# What I'm calling rate for the inverse gamma, Wikipedia calls scale
dinvgamma = function(x, shape, rate, log = FALSE){
  out = shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x
  if (log)
    return (out)
  return (exp(out))
}
rinvgamma = function(n, shape, rate){
  1/rgamma(n, shape = shape, rate = rate)
}

### MCMC
nburn = 20000
nmcmc = 100000

# Parameter objects
param_theta = matrix(0, nburn + nmcmc, n)
param_alpha = double(nburn + nmcmc)
param_psi = matrix(0, nburn + nmcmc, 2) # [,1] is zeta, [,2] is mu
param_beta = double(nburn + nmcmc)

# other mcmc
sig_beta = 0.01
acc_beta = double(nburn + nmcmc)
sig_psi = diag(0.01, 2)
acc_psi = double(nburn + nmcmc)

window = 200

# Priors
prior_beta_a  = 0   # Normal (mean)
prior_beta_b  = 1   # Normal (sd)
prior_alpha_a = 1   # Gamma (shape)           alpha controls n_star (discreteness of G)
prior_alpha_b = 1/2 # Gamma (rate)
prior_zeta_a  = 1   # Gamma (shape)
prior_zeta_b  = 1/4 # Gamma (rate)
prior_mu_a    = 1   # Gamma (shape)
prior_mu_b    = 1/4 # Gamma (rate)


# Initial values
#param_theta[1,] = y
param_beta[1] = rnorm(1, prior_beta_a, prior_beta_b)
param_alpha[1] = rgamma(1, prior_alpha_a, prior_alpha_b)
param_psi[1,] = c(rgamma(1, prior_zeta_a, prior_zeta_b),
                  rgamma(1, prior_mu_a, prior_mu_b))
#param_theta[1,] = rep(0, n)
param_theta[1,] = y


# Iterations
for (iter in 2:(nburn + nmcmc)){
  cat("\r", iter, "/", nburn + nmcmc)
  
  # Keep track of current draws
  curr_theta = param_theta[iter-1,]
  curr_alpha = param_alpha[iter-1]
  curr_beta = param_beta[iter-1]
  curr_zeta = param_psi[iter-1,1]
  curr_mu = param_psi[iter-1,2]
  
  # Update thetas
  for (i in 1:n){
    n_j_minus = table(curr_theta[-i])
    theta_star_minus = as.numeric(names(n_j_minus))
    n_star_minus = length(theta_star_minus)
    
    temp_p = (curr_zeta / curr_mu) / (exp(curr_beta*x[i]) + curr_zeta / curr_mu)
    # see report
    q0 = dnbinom(y[i], size = curr_zeta, prob = temp_p)
    
    qj = dpois(y[i], theta_star_minus*exp(curr_beta*x[i]))
    
    # Probabilities of determining which to draw
    # Calculate A (see report)
    A = curr_alpha * q0 / (curr_alpha * q0 + sum(n_j_minus * qj))
    
    # Calculate Bj's (see report)
    Bj = n_j_minus * qj / (curr_alpha * q0 + sum(n_j_minus * qj))
    
    # Make the update
    draw = sample(n_star_minus + 1, 1, prob = c(A, Bj))
    if (draw == 1){ # Make a draw from h
      curr_theta[i] = rgamma(1, y[i] + curr_zeta, exp(curr_beta*x[i]) + curr_zeta/curr_mu)
    } else { # Make a draw from the existing groups
      curr_theta[i] = theta_star_minus[draw-1]
    }
  }
  n_j = table(curr_theta)
  theta_star = as.numeric(names(n_j))
  n_star = length(theta_star)
  
  # Update beta
  cand = rnorm(1, curr_beta, sig_beta)
  temp_curr = dnorm(curr_beta, prior_beta_a, prior_beta_b, log = TRUE) +
    sum(dpois(y, curr_theta*exp(curr_beta*x), log = TRUE))
  temp_cand = dnorm(cand, prior_beta_a, prior_beta_b, log = TRUE) +
    sum(dpois(y, curr_theta*exp(cand*x), log = TRUE))
  if (log(runif(1)) < temp_cand - temp_curr){
    curr_beta = cand
    acc_beta[iter] = 1
  }
  
  # Update psi = (zeta, mu)
  cand = mvrnorm(1, c(curr_zeta, curr_mu), sig_psi)
  if (all(cand > 0)){
    temp_curr = dgamma(curr_zeta, prior_zeta_a, prior_zeta_b, log = TRUE) +
      dgamma(curr_mu, prior_mu_a, prior_mu_b, log = TRUE) +
      sum(dgamma(theta_star, curr_zeta, curr_zeta / curr_mu, log = TRUE))
    temp_cand = dgamma(cand[1], prior_zeta_a, prior_zeta_b, log = TRUE) +
      dgamma(cand[2], prior_mu_a, prior_mu_b, log = TRUE) +
      sum(dgamma(theta_star, cand[1], cand[1] / cand[2], log = TRUE))
    if (log(runif(1)) < temp_cand - temp_curr){
      curr_zeta = cand[1]
      curr_mu = cand[2]
      acc_psi[iter] = 1
    }
  }
  
  # Update alpha
  # Use an auxiliary variable eta to draw a new alpha
  eta = rbeta(1, curr_alpha + 1, n)
  eps = (prior_alpha_a + n_star - 1) /
    (n*(prior_alpha_b - log(eta)) + prior_alpha_a + n_star - 1)
  if (runif(1) < eps){
    curr_alpha = rgamma(1, prior_alpha_a + n_star, prior_alpha_b - log(eta))
  } else {
    curr_alpha = rgamma(1, prior_alpha_a + n_star - 1, prior_alpha_b - log(eta))
  }
  
  # Extra improvement step (only updating the clusters, not individuals)
  for (j in 1:n_star){
    temp_w = which(curr_theta == theta_star[j])
    curr_theta[temp_w] = rgamma(1,
                             curr_zeta + sum(y[temp_w]),
                             curr_zeta / curr_mu + sum(exp(curr_beta*x[temp_w])))
  }
  
  # Update the curr_* objects into the regular ones
  param_theta[iter,] = curr_theta
  param_alpha[iter] = curr_alpha
  param_beta[iter] = curr_beta
  param_psi[iter,1] = curr_zeta
  param_psi[iter,2] = curr_mu
  
  # Tune candidate sigmas
  if ((floor(iter/window) == iter/window) && (iter <= nburn)){
    sig_beta = autotune(mean(acc_beta[(iter-window+1):iter]),
                        target = 0.25, k = window/50) * sig_beta
    sig_psi = autotune(mean(acc_psi[(iter-window+1):iter]),
                       target = 0.25, k = window/50) *
      (sig_psi + window * var(param_psi[(iter-window+1):iter,]) / iter)
  }
  
  if (iter == (nburn + nmcmc))
    cat("\n")
}

### Remove Burn-in
param_theta = tail(param_theta, nmcmc)
param_alpha = tail(param_alpha, nmcmc)
param_beta = tail(param_beta, nmcmc)
param_psi = tail(param_psi, nmcmc)


acc_beta = tail(acc_beta, nmcmc)
acc_psi = tail(acc_psi, nmcmc)

# setwd("~/BNP_Homework/HW4")
# save(param_theta, param_alpha, param_beta, param_psi, acc_beta, acc_psi,
#      file = "prob2.RData")
load("prob2.RData")
mean(acc_beta)
mean(acc_psi)

# Post samples beta, psi[,1:2] = (zeta, mu)
trace_dens_acf(param_psi[,1], main1 = expression("Trace plot for "~ zeta),
               main2 = expression("Density for "~ zeta),
               main3 = expression("ACF for "~ zeta),
               ylab1 = expression(zeta),
               xlab2 = expression(zeta))
trace_dens_acf(param_psi[,2], main1 = expression("Trace plot for "~ mu),
               main2 = expression("Density for "~ mu),
               main3 = expression("ACF for "~ mu),
               ylab1 = expression(mu),
               xlab2 = expression(mu))
trace_dens_acf(param_alpha, main1 = expression("Trace plot for "~ alpha),
               main2 = expression("Density for "~ alpha),
               main3 = expression("ACF for "~ alpha),
               ylab1 = expression(alpha),
               xlab2 = expression(alpha))
trace_dens_acf(param_beta, main1 = expression("Trace plot for "~ beta),
               main2 = expression("Density for "~ beta),
               main3 = expression("ACF for "~ beta),
               ylab1 = expression(beta),
               xlab2 = expression(beta))
title("Acceptance rate = 0.2698083", line = -2)

### Posterior for n*
par(mfrow = c(1,2))
clusters = apply(param_theta, 1, function(x) length(unique(x)))
### Get theta groups
setwd("~/BNP_Homework/HW4/")
# save(param_theta, file = "param_theta.RData")
# load("param_theta.RData")
library(doMC)
library(fields)
# split work across 6 cores in parallel
dim(param_theta)

registerDoMC(6)
theta_tables = function(j){
  out = rep(0, n)
  temp = table(do.call("c", apply(param_theta, 1, function(x) which(x == x[j]))))
  out[as.numeric(names(temp))] = temp
  return (out)
}

group = (foreach(j = 1:n, .combine = rbind) %dopar% theta_tables(j)) / nmcmc
# par(mfrow = c(2,2))
image.plot(group[ord, ord], main = "Pairwise Clustering Probabilities", axes = F, col = heat.colors(250))
ticks_names = seq(0,32, 1)
axis(1, at=seq(0,1,length.out = length(ticks_names)) , labels=ticks_names, las = 2) 
axis(2, at=seq(0,1,length.out = length(ticks_names)) , labels=ticks_names, las = 2) 
box()

plot(table(clusters) / nmcmc, main = expression(n^"*"), cex.main = 2, lwd = 3, ylab="Mass")
a = par("pin")
b = par("plt")

par("pin" = a, "plt" = b)
#

### Box plots of the theta's
qlines = apply(param_theta, 2, quantile, c(0.025, 0.975))
mline = apply(param_theta, 2, mean)

#pdf("./figs/boxplots_1.pdf", width = 8, height = 8)
par(mfrow = c(1,1), mar = c(4.1, 2.1, 2.1, 1.1))
plot(y[ord], col = 'forestgreen', pch = 20, cex = 1,
     xlab = expression("Ordered "~ theta[i] ~ "at index "~ i),
     main = expression("Posterior cluster locations for each" ~ theta[i]))
polygon(c(seq(1,32),rev(seq(1,32))),
        c(qlines[1,ord],rev(qlines[2,ord])),
        col=alpha("lightgreen", 0.7), border = NA)
boxplot(param_theta[,ord], pch = 20, cex = 0.1, add = TRUE)
points(y[ord], col = 'forestgreen', pch = 20, cex = 1)
# points(x, y, pch = 20, ylim = range(c(y, q0)), lwd = 1.5)
lines(mline[ord], col = 'darkgreen', lwd = 1.5)
legend("topleft", legend=c("Data", "Posterior Mean", "95% Posterior CI"), fill=c(NA, NA, alpha("lightgreen", 0.5)),
       bty="n", lwd = c(NA, 3, NA), pch = c(20, NA, NA), border = c(NA, NA, NA), 
       col = c("forestgreen", "forestgreen", "lightgreen"), cex = 1.3) 
#dev.off()

#dev.off()






par(mfrow = c(1, 1), mar = c(3.1, 2.1, 2.1, 1.1))
plot(0, type='n', xlim = range(param_theta), ylim = c(0, 3.0))
for (k in 1:n){
  lines(density(param_theta[,k]), col = k, lty = k)
}

### Posterior predictions
### list of theta_stars
n_j = apply(param_theta, 1, table)
theta_star = lapply(n_j, function(x) as.numeric(names(x)))

# A new cluster (to predict further with this, we would need to specify a particular x)
theta_0 = double(length(param_beta))
for (i in 1:length(param_beta)){
  if (floor(i/window) == i/window)
    cat("\r", i, "/", nmcmc)
  prob = c(param_alpha[i] / (param_alpha[i] + n), n_j[[i]] / (param_alpha[i] + n))
  draw = sample(length(prob), 1, replace = FALSE, prob = prob)
  if (draw == 1){
    theta_0[i] = rgamma(1, param_psi[i,1], param_psi[i,1] / param_psi[i,2])
  } else {
    theta_0[i] = theta_star[[i]][draw - 1]
  }
  if (i == nmcmc)
    cat("\n")
}
y_0 = matrix(0, nmcmc, n)
for (i in 1:n){
  y_0[,i] = rpois(nmcmc, theta_0*exp(param_beta*x[i]))
}


### predictions of observations
q0 = apply(y_0, 2, quantile, c(0.025, 0.975))
m0 = apply(y_0, 2, mean)


# Predictions at each x_i
par(mfrow = c(1,2))
plot(x, y, pch = 20, ylim = range(c(y, q0)), lwd = 1.5, main = "Predictions")
polygon(c(x[ord_x],rev(x[ord_x])),
        c(q0[1,ord_x],rev(q0[2,ord_x])),
        col=alpha("lightgreen", 0.9), border = NA)
points(x, y, pch = 20, ylim = range(c(y, q0)), lwd = 1.5, main = "Predictions")
lines(x[ord_x], m0[ord_x],col = "forestgreen")
legend("topleft", legend=c("Data", "Posterior Mean", "95% Posterior CI"), fill=c(NA, NA, alpha("lightgreen", 0.5)),
       bty="n", lwd = c(NA, 3, NA), pch = c(20, NA, NA), border = c(NA, NA, NA), 
       col = c("black", "forestgreen", "lightgreen"), cex = 1.3) 

### Observed vs. Fitted plot
plot(0, type='n', xlim = range(y), ylim = range(c(y, q0)),
     xlab = "Observed", ylab = "Fitted", main = "Observed vs. Fitted", cex.lab = 1.3)
segments(x0 = jity, y0 = q0[1,], y1 = q0[2,], col = 'lightgreen')
points(jity, m0, col = 'darkgreen', pch = 20)
abline(0, 1, lwd = 3, lty = 2)

criterion(y, y_0, 0)                     # 457.03
criterion(y, y_0, Inf)                   # 590.09
criterion(y, y_0, Inf) - criterion(y, y_0, 0) # 133.06

#### GG ####
(g2=sum((apply(na.omit(y_0),2, mean)-x)^2))

(p2=sum(apply(na.omit(y_0),2, sd)))
(gg_criterion2=g2+p2) # 496.9563
