#### First load some functions
library(scales)
# Function to plot the traces and posterior density
plot_trace_density = function(param, prob = 0.95, precision = 1000, main1="main1", main2="main2",
                              xlab1="xlab1", xlab2 = "Iteration", par_settings = T,
                              ylab1 = "ylab1", ylab2 = "Density", breaks = 20){
  # param - the mcmc samples for a parameter
  # prob - high posterior density interval (i.e 95%)
  # precision - controls smoothness of the curve
  if(par_settings == T){
    par(mfrow = c(1,2), mar = c(4.1, 4.1, 3.1, 1.1), oma = c(0,0,1,0))
  }
  range = seq(0, 1-prob, length=precision)
  range = cbind(range, range+prob)
  best = range[which.min(apply(range, 1, function(y)
    diff(quantile(param, y)))),]
  quants = quantile(param, best)
  dens = density(param)
  # trace plot
  plot(param, type = 'l', main = main1,
       xlab = "Iteration", ylab = ylab1)
  abline(h = mean(param), col = "red", lwd = 3)
  # histogram and density
  hist(param, xlab = xlab1, probability = TRUE,ylim = c(0,max(dens$y)), 
       main = main2)
  lines(density(param))
  x1 <- min(which(dens$x >= quants[1]))  
  x2 <- max(which(dens$x <=  quants[2]))
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha("pink", 0.7)))
  abline(v = mean(param), col = "red", lwd = 3)
}

# Example 
mu_test = rnorm(10000, 0, 1)
plot_trace_density(mu_test, main1 = expression("Traceplot for" ~ mu), 
                   main2 = expression("Posterior density for" ~ mu),
                   ylab1 = expression(mu), xlab1 = expression(mu), 
                   xlab2 = expression(mu), ylab2 = expression(mu))


##################################################################################
# Posterior density plots for multiple component mixtures
post_dens_mix = function(x, dens, prob = 0.95, klength = 5000){
  c.prob = 1
  temp.prob = 0
  k = seq(min(dens$y), max(dens$y), length=klength)
  # i = 2 to prevent certain problems, test.x4 had an issue
  # with the probability 0 region in the middle, (doesn't always
  # occur) perhaps fix by doing f(x) > k, instead of f(x) >= k?
  i = 2
  zeros = function(y, k, return.max.min = FALSE){
    # y is expected to be density(x)$y
    out = NULL
    int = NULL
    for (i in 2:(length(y)-1)){
      # condition to check when the height crosses k
      if ((y[i] > k && y[i-1] < k) || (y[i] < k && y[i-1] > k)){
        # get the x closer to k
        out = c(out, ifelse(which.min(abs(y[c(i,i-1)]-k))==1,
                            i, i-1))
        # -1 if lower interval, +1 if upper
        int = c(int, -sign(y[i] - y[i-1]))
      }
      # check if the height exactly equals k
      if (y[i] == k){
        out = c(out, i)
        # 0 if a maximum or minimum, -1 if lower, +1 if upper
        # y[i] can only be max or min if y[i] = k, so don't
        # check this condition for when height crosses k
        int = c(int, -sign(sign(y[i]-y[i-1]) +
                             sign(y[i+1]-y[i])))
      }
    }
    # ensure that first value is lower end and last is upper end
    if (is.null(int))
      return (NULL)
    if (int[1] == 1){
      int = c(-1, int)
      out = c(1, out)
    }
    if (int[length(int)] == -1){
      int = c(int, 1)
      out = c(out, length(y))
    }
    if (return.max.min)
      return (out)
    return (out[as.logical(int)])
  }
  # repeat until area under curve is <= specified prob
  # (note 14 jun: perhaps do some kind of iterative
  # convergence to reduce the number of iterations;
  # start in the middle i = floor(klength/2), and if
  # temp.prob is too low, set i=floor((i + 0)/2), else
  # set i = floor((i + klength)/2), repeat until value
  # is sufficiently close to prob. need to keep track of
  # previous "lower" and "upper" bounds
  while (c.prob > prob){
    temp.prob = 0
    int = zeros(dens$y, k[i])
    if (!is.null(int)){
      if (length(int) > 0){
        # sum the area in the intervals
        for (j in 1:(length(int)/2))
          temp.prob = temp.prob + mean(x >= dens$x[
            int[2*j-1]] & x <= dens$x[int[2*j]])
        # update (i think this always occurs)
        if (c.prob > temp.prob)
          c.prob = temp.prob
      }
    }
    i = i + 1
  }
  vec = dens$x[int]
  mat_ints = matrix(vec, ncol = 2, byrow = TRUE)
  return (mat_ints)
}



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



### MCMC
nburn = 1000
nmcmc = 8000

# Parameter objects
param_theta = matrix(0, nburn + nmcmc, n)
param_alpha = double(nburn + nmcmc)
param_mt = matrix(0, nburn + nmcmc, 2) # param_mt[,1] is mu, param_mt[,2] is tau^2
param_phi = double(nburn + nmcmc)

# Priors
prior_phi_a   = 3   # Inverse Gamma (shape)
prior_phi_b   = 15  # Inverse Gamma (rate)
prior_alpha_a = 1   # Gamma (shape)           alpha controls n_star (discreteness of G)
prior_alpha_b = 1   # Gamma (rate)
prior_mu_a = 0   # Normal (mean)
prior_mu_b  = 3   # Normal (variance)
prior_tau2_a  = 3   # Inverse Gamma (shape)
prior_tau2_b  = 10  # Inverse Gamma (rate)


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
    
    # q0 here is a normal density (after annoying integration)
    q0 = dnorm(y[i], curr_mu, sqrt(curr_phi + curr_tau2))
    
    # as is qj, by construction
    qj = dnorm(y[i], theta_star_minus, sqrt(curr_phi))
    
    # Probabilities of determining which to draw
    # Calculate A
    A = curr_alpha * q0 / (curr_alpha * q0 + sum(n_j_minus * qj))
    
    # Calculate B's
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

#### Figure 1 ####
plot_trace_density(param_mt[,1], main1 = expression("Traceplot for" ~ mu), 
                   main2 = expression("Posterior density for" ~ mu),
                   ylab1 = expression(mu), xlab1 = expression(mu), 
                   xlab2 = expression(mu), ylab2 = expression(mu))

plot_trace_density(param_mt[,2], main1 = expression("Traceplot for" ~ tau^2),
                   main2 = expression("Posterior density for" ~ tau^2), breaks = 100,
                   ylab1 = expression(tau^2), xlab1 = expression(tau^2), 
                   xlab2 = expression(tau^2), ylab2 = expression(tau^2))

plot_trace_density(param_phi, main1 = expression("Traceplot for" ~ phi),
                   main2 = expression("Posterior density for" ~ phi), breaks = 20,
                   ylab1 = expression(phi), xlab1 = expression(phi), 
                   xlab2 = expression(phi), ylab2 = expression(phi))




#### Figure 2 ####
### Posterior for n*
# hpd.alpha = hpd.uni(param_alpha)
clusters = apply(param_theta, 1, function(x) length(unique(x)))
### Get theta groups
library(doMC)
library(fields)
# split work across 6 cores in parallel
registerDoMC(6)
par.fun = function(j){
  out = double(n)
  temp = table(unlist(apply(param_theta, 1, function(x) which(x == x[j]))))
  out[as.numeric(names(temp))] = temp
  return (out)
}
group = (foreach(j = 1:n, .combine = rbind) %dopar% par.fun(j)) / nmcmc

?foreach
par(mfrow = c(2,2))
plot(table(clusters) / nmcmc, main = expression(n^"*"), cex.main = 2, lwd = 3, ylab="Mass")
a = par("pin")
b = par("plt")
image.plot(group[ord, ord], main = "Cluster groupings heatmap")
par("pin" = a, "plt" = b)
plot_trace_density(param_alpha, main1 = expression("Traceplot for" ~ alpha), par_settings = F,
                   main2 = expression("Posterior density for" ~ alpha),
                   ylab1 = expression(alpha), xlab1 = expression(alpha), 
                   xlab2 = expression(alpha), ylab2 = expression(alpha))
dim(group)
#### Figure 2a ####
par(mfrow = c(1,1))
clusters = apply(param_theta, 1, function(x) length(unique(x)))
plot(table(clusters) / nmcmc, main = expression(n^"*"), cex.main = 2, lwd = 3, ylab="Mass")

#### Figure 3 ####
par(mfrow = c(1,1))
quantile_thetass = apply(param_theta, 2, quantile, c(0.025, 0.975))
median_theta = apply(param_theta, 2, mean)

# Plot each cluster label theta_i, sorted by the order of observations
plot(median_theta[ord], col = 'red', lwd = 3.0, type = 'l', 
     xlab = expression(theta ~ Index), ylab ="", cex.main = 2,
     main = expression("Posterior cluster locations for each" ~ theta[i]))
ci_shade = seq(1,250)
polygon(c(ci_shade,rev(ci_shade)),
        c(quantile_thetass[2,ord],rev(quantile_thetass[1,ord])),
        col=alpha("pink", 0.3), border = NA)
points(y[ord], col = 'black', pch = 20, cex = 0.5)
legend("topleft", legend=c("Data", "Posterior medians", "95% Posterior CI"), fill=c(NA, NA, alpha("pink", 0.5)),
       bty="n", lwd = c(NA, 3, NA), pch = c(20, NA, NA), border = c(NA, NA, NA), 
       col = c("black", "red", "pink"), cex = 1.3) 

#### Figure 4 ####
# Posterior predictions
# list of theta_stars
n_j = apply(param_theta, 1, table)
theta_star = lapply(n_j, function(x) as.numeric(names(x)))

# A new cluster
theta_0 = double(nmcmc)
for (i in 1:nmcmc){
  cat("\r", i, "/", nmcmc)
  prob = c(param_alpha[i] / (param_alpha[i] + n), n_j[[i]] / (param_alpha[i] + n))
  draw = sample(length(prob), 1, replace = FALSE, prob = prob)
  if (draw == 1){
    theta_0[i] = rnorm(1, param_mt[i, 1], sqrt(param_mt[i, 2]))
  } else {
    theta_0[i] = theta_star[[i]][draw - 1]
  }
  if (i == nmcmc)
    cat("\n")
}
y_0 = rnorm(nmcmc, theta_0, sqrt(param_phi))

# post pred theta0
theta0_mat = post_dens_mix(theta_0, density(theta_0))
y0_mat = post_dens_mix(y_0, density(y_0))
plot(density(theta_0), col = "red", main = expression("Posterior predictive" ~ p(theta[0]~ "|data")), xlab = "")
dens = density(theta_0)
theta_lower = apply(theta0_mat, 1, function(mat) min(which(dens$x >= mat))  )
theta_upper = apply(theta0_mat, 1, function(mat) max(which(dens$x <=  mat)))
for(i in 1:dim(theta0_mat)[1]){
  dens = density(theta_0)
  x1 = theta_lower[i]
  x2 = theta_upper[i]
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha("pink", 0.7), border = F))
}
legend("topleft", legend=c("Posterior predictive density", "95% HPD Interval"), fill=c(NA, alpha("pink", 0.7)),
       bty="o", lwd = c(1, NA), border = c(NA, NA), 
       col = c("red", "pink"), cex = 1.3) 



#### Figure 5 ####
# post pred y0
plot(density(y_0),col = "red", main = expression("Posterior predictive density" ~ p(y[0]~ "|data")), xlab = "")
dens = density(y_0)
y0_lower = apply(y0_mat, 1, function(mat) min(which(dens$x >= mat))  )
y0_upper = apply(y0_mat, 1, function(mat) max(which(dens$x <=  mat)))
for(i in 1:dim(y0_mat)[1]){
  dens = density(y_0)
  x1 = y0_lower[i]
  x2 = y0_upper[i]
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha("pink", 0.7), border = NA))
}
lines(density(y))
legend("topleft", legend=c("Data", "Posterior predictive density", "95% HPD Interval"), fill=c(NA, NA, alpha("pink", 0.7)),
       bty="o", lwd = c(1, 1, NA), border = c(NA, NA, NA), 
       col = c("black", "red", "pink"), cex = 1.3) 


### Prior predictive
library(MCMCpack)
# Cluster

prior.theta_0 = double(20000)
prior.y_0 = double(length(prior.theta_0))
tvec = seq(-10, 10, length = 1500)
k = length(tvec)
for (i in 1:length(prior.y_0)){
  cat("\r", i, "/", length(prior.y_0))
  temp.phi = rinvgamma(1, prior_alpha_a, prior_alpha_b)
  temp.alpha = rgamma(1, prior_alpha_a, prior_alpha_b)
  temp.mu = rnorm(1, prior_mu_a, sqrt(prior_mu_b))
  temp.tau = sqrt(rinvgamma(1, prior_tau2_a, prior_tau2_b))
  dG0 = pnorm(tvec, temp.mu, temp.tau) - pnorm(c(-Inf, tvec[-k]), temp.mu, temp.tau)
  flag = TRUE
  while (flag){
    flag = FALSE
    dG = rdirichlet(1, temp.alpha * dG0)
    if (any(is.na(dG)))
      flag = TRUE
  }
  keep = which(dG != 0)
  if (length(keep) == 1){
    prior.theta_0[i] = tvec[keep]
  } else {
    prior.theta_0[i] = sample(tvec[keep], 1, prob = dG[keep])
  }
  prior.y_0[i] = rnorm(1, prior.theta_0[i], sqrt(temp.phi))
  if (i == length(prior.y_0))
    cat("\n")
}

pdf("./figs/prior_1.pdf", height = 6, width = 12)
par(mfrow = c(1,2), mar = c(3.1, 2.1, 2.1, 1.1))
plot(density(prior.theta_0[1:i]), main = expression("Prior predictive for" ~ theta[0]), lwd = 2,
     col = 'blue')
#curve(dnorm(x, mean(prior.theta_0[1:i]), sd(prior.theta_0[1:i])), add = TRUE, col = 'red', lwd = 2)
legend("topleft", legend = "KDE of draws", box.lwd = 0, box.lty = 0, col = 'blue', lwd = 2, cex = 1.3)
#legend("topleft", legend = c("KDE of draws", "Approximate normal"),
#    box.lwd = 0, box.lty = 0, col = c('blue', 'red'), lwd = 2, cex = 1.3)

plot(density(prior.y_0[1:i]), main = expression("Prior predictive for" ~ y[0]), lwd = 2,
     col = 'blue', xlim = c(-10,10), ylim = c(0,0.16))
#curve(dnorm(x, mean(prior.y_0[1:i]), sd(prior.y_0[1:i])), add = TRUE, col = 'red', lwd = 2)
lines(density(y), lwd =3)
legend("topleft", legend = c("KDE of draws", "Data"),
       box.lwd = 0, box.lty = 0, col = c('blue', 'black'), lwd = 2, cex = 1.3)
#legend("topleft", legend = c("KDE of draws", "Approximate normal", "Data"),
#    box.lwd = 0, box.lty = 0, col = c('blue', 'red', 'black'), lwd = 2, cex = 1.3)
dev.off()