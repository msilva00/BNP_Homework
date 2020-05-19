library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/BNP_Homework/master/HW3/hwk3-data.txt"

setwd("~/BNP_Homework/HW3/")
## Installation requires many additional steps, see README file for instructions for mac users
# url = "https://cran.r-project.org/src/contrib/Archive/DPpackage/DPpackage_1.1-7.4.tar.gz"
# pkgFile = "DPpackage_1.1-7.4.tar.gz"
# download.file(url = url, destfile = pkgFile)
# install.packages(c("ada", "ipred", "evd"))
# install.packages(pkgs=pkgFile, type="source", repo=NULL)

library(DPpackage)
### Posterior predictive loss criterion
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
plot(density(y))
curve(0.2*dnorm(x, -5, 1) + 0.5*dnorm(x, 0, 1) + 0.3*dnorm(x, 3.5, 1), add = TRUE, col = 'red')
legend("topleft", legend = c("data", "truth"), col = c("black", "red"), lwd = 1, box.lwd = 0, cex = 1.5)

nburn = 2000
nmcmc = 10000


prior1 = list(m2 = 0, s2 = sqrt(3),    # prior mean and sd for mean of normal part of G0
              tau1 = 2, tau2 = 2,      # twice the shape and rate for kappa
              nu1 = 3,                 # fixed c, shape in inverse gamma part of G0
              nu2 = 5, psiinv2 = 1/10, # rate in inverse gamma part of G0 (an inverse gamma itself)
              a0 = 1, b0 = 1)          # priors on alpha (same as prob 1)
mcmc1 = list(nburn = nburn, nsave = nmcmc, ndisplay = 100)


fit1 = DPdensity(y = y, prior = prior1, mcmc = mcmc1, state = NULL, status = TRUE)

str(fit1$save.state)

plot((fit1$save.state$thetasave[,1]), main="mean of normal part", type='l')
plot((fit1$save.state$thetasave[,2]), main=expression(kappa), type='l')
plot((fit1$save.state$thetasave[,3]), main="rate of ig part", type='l')
plot((fit1$save.state$thetasave[,4]), type='l')
plot((fit1$save.state$thetasave[,5]), main=expression(alpha), type='l')


mu_mean = fit1$save.state$thetasave[,1]
kappa = fit1$save.state$thetasave[,2]
beta_rate = fit1$save.state$thetasave[,3]  
ncluster = fit1$save.state$thetasave[,4]
alpha = fit1$save.state$thetasave[,5]

save_outs = data.frame(mu_mean, kappa, beta_rate, ncluster, alpha)

write.csv(save_outs, "mcmc_outs.csv")

plot(density(fit1$save.state$thetasave[,1]), main="mean of normal part")
plot(density(fit1$save.state$thetasave[,2]), main=expression(kappa))
plot(density(fit1$save.state$thetasave[,3]), main="rate of ig part")
plot(table(fit1$save.state$thetasave[,4])/nmcmc)
plot(density(fit1$save.state$thetasave[,5]), main=expression(alpha))

plot_trace_density(mu_mean, main1 = expression("Traceplot for posterior samples of "~mu), ylab1 = expression(mu),
                   main2 = expression("Histogram for posterior samples of "~mu), xlab2 = expression(mu))

plot_trace_density(kappa, main1 = expression("Traceplot for posterior samples of "~kappa), ylab1 = expression(kappa),
                   main2 = expression("Histogram for posterior samples of "~kappa), xlab2 = expression(kappa))

plot_trace_density(beta_rate, main1 = expression("Traceplot for posterior samples of "~ beta ~" (rate for IG)"), 
                   ylab1 = expression(beta),
                   main2 = expression("Histogram for posterior samples of "~beta), 
                   xlab2 = expression(beta))

plot_trace_density(alpha, main1 = expression("Traceplot for posterior samples of "~alpha), ylab1 = expression(alpha),
                   main2 = expression("Histogram for posterior samples of "~alpha), xlab2 = expression(alpha))


plot(table(fit1$save.state$thetasave[,4])/nmcmc, xlab = "clusters", ylab = "Mass")

?DPdensity
print(fit1)

plot(density(fit1$save.state$randsave[,501]))
plot(density(fit1$save.state$randsave[,502]))
fit1$save.state$randsave[,seq(1, 2*n, by = 2)]
# fit1$save.state$randsave[,seq(2, 2*n, by = 2)], 2, mean)

z1 = apply(fit1$save.state$randsave[,c(2*245-1,2*245+0)], 1,
           function(x) length(unique(x)))
z2 = apply(fit1$save.state$randsave[,c(2*140-1,2*140+0)], 1,
           function(x) length(unique(x)))
z3 = apply(fit1$save.state$randsave[,c(2*245-1,2*245+0,2*140-1,2*140+0)], 1,
           function(x) length(unique(x)))
table(z1); table(z2); table(z3)

### Plot of posterior means for location and scale of eath observation
cluster_loc = apply(fit1$save.state$randsave[,seq(1, 2*n, by = 2)], 2, mean)
cluster_sca = apply(fit1$save.state$randsave[,seq(2, 2*n, by = 2)], 2, mean)
plot(cluster_loc, cluster_sca, pch = 20, cex = 0.5)

plot(0, type = 'n', xlim = range(fit1$save.state$randsave[,seq(1, 2*n, by = 2)]),
     ylim = range(fit1$save.state$randsave[,seq(2, 2*n, by = 2)]))
points(cluster_loc, cluster_sca, pch = 20, cex = 0.6, col = 'red')


### Predictions on mixing parameters
library(MASS)
theta = fit1$save.state$randsave[,c(501, 502)]
d = kde2d(theta[,1], theta[,2], n = 200)
points(theta, pch = 20, cex = 0.5)

points(fit1$save.state$randsave[,c(501, 502)], pch = 20, cex = 0.5,
       xlab = "location", ylab = "scale", main = "clusters")

hpd.theta = apply(theta, 2, function(x) hpd.mult(x, density(x)))
par(mfrow = c(2,1), mar = c(3.1, 2.1, 2.1, 1.1))
hpd.plot(density(theta[,1]), hpd.theta[[1]], main = "new location")
hpd.plot(density(theta[,2]), hpd.theta[[2]], main = "new scale")

### Predictions on new observation
y0 = fit1$save.state$randsave[,503]
hpd.y = hpd.mult(y0, density(y0))

# pdf("./figs/pred_2.pdf", width = 12, height = 6)
par(mfrow=c(1,2),  mar = c(4.1, 4.1, 2.1, 1.1))
contour(d, xlim = c(-7, 7), ylim = c(0, 3), nlevels = 50, cex.lab = 1.3,
        xlab = expression(theta), ylab = expression(phi), main = "Clusters")
#par(mfrow = c(1,1), mar = c(3.1, 2.1, 2.1, 1.1))
hpd.plot(density(y0), hpd.y, lwd = 3, main=expression("Posterior predictions of new" ~ y[0]),
         xlab = "")
lines(density(y), lwd = 3)
legend("topleft", box.lty = 0, legend = "Data", lwd = 3, cex = 1.5)
dev.off()

### Replicate for each observation (don't really need to, but it makes post_pred_loss work)
y1 = matrix(rep(y0, n), ncol = n)

post_pred_loss(y, y1, 0)                      # 2390.80
post_pred_loss(y, y1, Inf)                    # 4745.31
post_pred_loss(y, y1, Inf) - post_pred_loss(y, y1, 0)   # 2354.51
