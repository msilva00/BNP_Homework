#!/usr/bin/env python
# coding: utf-8

# # 4. Simulation of Dirichlet process prior realizations
# ## Consider a $DP(\alpha, G_0)$ prior over the space of distributions (equivalently c.d.f's) $G$ on ${\rm I\!R}$, with $G_0 = N(0,1)$.
# 
# ### (a) Use both Ferguson’s original definition and Sethuraman’s constructive definition to generate(multiple) prior c.d.f.  realizations from the $DP(\alpha,N(0,1))$, for different values ofαrangingfromsmalltolarge.

# In[1]:


import numpy as np
from scipy.stats import norm, dirichlet, beta
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


# In[2]:


def DirPriorSample(samples = 5, N = 50, alpha = 5, mean = 0, sd = 1):
    # Interval bound
    bound = np.maximum(mean-4*sd, mean+4*sd) # not necessary
    
    # Discretization of the space
    x = np.linspace(-bound, bound, N+1) # just (-4,4)
    
    # Probability measure for each interval
    y =np.zeros(N+1)
    y[0] = norm.cdf(x[0], mean, sd)
    for i in range(1,N+1):
        y[i] = norm.cdf(x[i], mean, sd) - norm.cdf(x[i-1], mean, sd)
        
    y = np.append(y, 1 - norm.cdf(x[N], mean, sd))
    
    # Creating the non-negative measures
    param = alpha * y
    
    # Samplind from the dirichlet distribution
    samp_dir = dirichlet.rvs(param, samples)
    
    # Generating the CDF
    draw = np.cumsum(samp_dir.T, axis = 0)
    return np.array([draw, bound])


# In[3]:


test = DirPriorSample(samples = 5, N = 50, alpha = 5, mean = 0, sd = 1)


# Changing $\alpha$ changes the shapes of each realization of the prior distributions. As $\alpha$ grows, the distributions are more concentrated around the theoretical N(0,1) represented by the black line. This is due to the fact that $\alpha$ represents how confident we are that the prior distribution is true. For each simulation using the DP I am generating 10 CDF, each represented with a different color.

# In[4]:


np.random.seed(1)
samples = 10
N = 100


# $\alpha = 0.1$

# In[5]:


draws = DirPriorSample(samples = samples, N = N, alpha = 0.1, mean = 0, sd = 1)
results = draws[0]
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.plot(xx, results)
plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black")
plt.title(label = r'Ferguson`s Definition: $\alpha =0.1$',fontsize=20)
plt.show()


# $\alpha = 1$

# In[6]:


draws = DirPriorSample(samples = samples, N = N, alpha = 1, mean = 0, sd = 1)
results = draws[0]

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.plot(xx, results)
plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black")
plt.title(label = r'Ferguson`s Definition: $\alpha =1$',fontsize=20)
plt.show()


# $\alpha = 10$

# In[7]:


draws = DirPriorSample(samples = samples, N = N, alpha = 10, mean = 0, sd = 1)
results = draws[0]
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.plot(xx, results)
plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black")
plt.title(label = r'Ferguson`s Definition: $\alpha =10$',fontsize=20)
plt.show()


# $\alpha = 100$

# In[8]:


draws = DirPriorSample(samples = samples, N = N, alpha = 100, mean = 0, sd = 1)
results = draws[0]
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.plot(xx, results)
plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black")
plt.title(label = r'Ferguson`s Definition: $\alpha =100$',fontsize=20)
plt.show()


# Now,  we use Sethurman's constructive definition to generate prior CDF realizations. 

# In[9]:


def DirProcessSamplesSethurmans(N = 1000, alpha = 1):
    theta_vec = np.random.normal(0,1, N)
    
    z = np.random.beta(1, alpha, N)
    log_z = np.log(z)
    S_log = np.append(0.0, np.cumsum(np.log((1 - z)))[:N-1])
    log_w = log_z + S_log
    w = np.exp(log_w)
    return np.array([theta_vec,w])


# In[10]:


np.random.seed(1)


# In[11]:


colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
print(colors)


# $\alpha=0.1$

# In[12]:


N = 1000
samples = np.zeros((10, 2, N))
for i in range(10):
    samples[i] = DirProcessSamplesSethurmans(N = N, alpha = 0.1)
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
for i in range(10):
    results = samples[i]
    plt.plot(results[0][np.argsort(results[0])],
            np.cumsum(results[1][np.argsort(results[0])]), color = colors[i])

plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black")
plt.title(label = r'Stick breaking: $\alpha =0.1$',fontsize=20)


# $\alpha = 1$

# In[13]:


samples = np.zeros((10, 2, N))
for i in range(10):
    samples[i] = DirProcessSamplesSethurmans(N = N, alpha = 1)
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
for i in range(10):
    results = samples[i]
    plt.plot(results[0][np.argsort(results[0])],
            np.cumsum(results[1][np.argsort(results[0])]), color = colors[i])
plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black")
plt.title(label = r'Stick breaking: $\alpha =1$',fontsize=20)


# $\alpha = 10$

# In[14]:


samples = np.zeros((10, 2, N))
for i in range(10):
    samples[i] = DirProcessSamplesSethurmans(N = N, alpha = 10)
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
for i in range(10):
    results = samples[i]
    plt.plot(results[0][np.argsort(results[0])],
            np.cumsum(results[1][np.argsort(results[0])]), color = colors[i])
plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black")
plt.title(label = r'Stick breaking: $\alpha =10$',fontsize=20)


# $\alpha = 100$

# In[15]:


samples = np.zeros((10, 2, N))
for i in range(10):
    samples[i] = DirProcessSamplesSethurmans(N = N, alpha = 100)
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
for i in range(10):
    results = samples[i]
    plt.plot(results[0][np.argsort(results[0])],
            np.cumsum(results[1][np.argsort(results[0])]), color = colors[i])
plt.plot(np.linspace(-3,3), norm.cdf(np.linspace(-3,3)), color = "black")
plt.title(label = r'Stick breaking: $\alpha =100$',fontsize=20)


# ### (b) In addition to prior c.d.f. realizations, obtain, for each value of $\alpha$, the corresponding prior distribution for the mean functional $$\mu(G) = \int t dG(t) $$ and for the variance functional  $$\sigma^2(G) = \int t^2 dG(t) - \left\{\int t dG(t) \right\}^2$$

# Using Ferguson's definition of the DP, I can compute the functional mean and the functional variance for each simulation, fixing the value of $\alpha$.

# In[16]:


np.random.seed(1)
samples = 100
N = 100


# In[17]:


xx


# $\alpha = 0.1$

# In[18]:


results[:,1]


# In[19]:


draws = DirPriorSample(samples = samples, N = N, alpha = 0.1, mean = 0, sd = 1)
results = draws[0]

xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(xx, results[:,1], color = "gray", label = "Prior Draws")
plt.plot(xx, results, color = "grey",linewidth=0.5)
mean_fn = np.mean(results, axis = 1)
plt.plot(xx, mean_fn, color = "blue", label = "Estimated prior mean")
var_fn = np.var(results, axis = 1)
plt.plot(xx, var_fn, color = "green", label = "Estimated prior variance")
plt.legend(loc=1, bbox_to_anchor=(1, 0.5), fontsize = 15)
plt.title(label = r'Ferguson`s Definition: $\alpha =0.1$', fontsize = 20)
plt.show()


# $\alpha = 1$

# In[20]:


draws = DirPriorSample(samples = samples, N = N, alpha = 1, mean = 0, sd = 1)
results = draws[0]

xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(xx, results, color = "grey",linewidth=0.5)
plt.plot(xx, results[:,1], color = "gray", label = "Prior Draws")
mean_fn = np.mean(results, axis = 1)
plt.plot(xx, mean_fn, color = "blue", label = "Estimated prior mean")
var_fn = np.var(results, axis = 1)
plt.plot(xx, var_fn, color = "green", label = "Estimated prior variance")
plt.legend(loc=1, bbox_to_anchor=(1, 0.5), fontsize = 15)
plt.title(label = r'Ferguson`s Definition: $\alpha =1$', fontsize = 20)
plt.show()


# $\alpha = 10$

# In[21]:


draws = DirPriorSample(samples = samples, N = N, alpha = 10, mean = 0, sd = 1)
results = draws[0]

xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(xx, results, color = "grey",linewidth=0.5)
plt.plot(xx, results[:,1], color = "gray", label = "Prior Draws")
mean_fn = np.mean(results, axis = 1)
plt.plot(xx, mean_fn, color = "blue", label = "Estimated prior mean")
var_fn = np.var(results, axis = 1)
plt.plot(xx, var_fn, color = "green", label = "Estimated prior variance")
plt.legend(loc=1, bbox_to_anchor=(1, 0.5), fontsize = 15)
plt.title(label = r'Ferguson`s Definition: $\alpha =10$', fontsize = 20)
plt.show()


# The variance decreases as $\alpha$ increases. In each case it grows around the central values and it goes to zero as the distribution reaches values 0 or 1. The functional mean becomes more similar to the c.d.f. as $\alpha$ gets larger.

# ### (c) Consider also a simulation under a mixture of DPs (MDP) prior, which extends the DP above by adding a prior for $\alpha$. Therefore, the MDP prior for G is defined such that, $G|\alpha \sim DP(\alpha, N(0,1)$, with a prior assigned to the precision parameter $\alpha$ from its prior. You can work with a gamma prior for $\alpha$ and 2-3 different choices for the gamma prior parameters.

# Next, we extend the Dirichlet process by adding a gamma prior for M. 
# 
# $\alpha \sim Gamma(3,3)$

# In[22]:


np.random.seed(1)
samples = 100
N = 100
a0 = 3
b0 = 3


# In[23]:


alpha_prior = np.random.gamma(a0,b0, size = 1)
prior_vec = alpha_prior
draws = DirPriorSample(samples = 1, N = N, alpha = alpha_prior, mean = 0, sd = 1)


# In[24]:


results = draws[0]
results_save = np.zeros((results.shape[0],samples))

for i in range(samples):
    alpha_prior = np.random.gamma(a0,b0, size = 1)
    prior_vec = np.append(prior_vec, alpha_prior)
    draws = DirPriorSample(samples = 1, N = N, alpha = alpha_prior, mean = 0, sd = 1)
    new_results = draws[0]
    results_save[:,i] = new_results.flatten()
    
xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.plot(xx, results_save, color = "grey",linewidth=0.5)
mean_fn = np.mean(results_save, axis = 1)
plt.plot(xx, mean_fn, color = "blue", label = "Estimated prior mean")
var_fn = np.var(results_save, axis = 1)
plt.plot(xx, var_fn, color = "green", label = "Estimated prior variance")
plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black", label = "CDF of N(0,1)")
plt.legend(loc=1, bbox_to_anchor=(1, 0.5), fontsize = 15)
plt.title(label = r'$\alpha \sim Gamma(3,3)$', fontsize = 20)


# In[25]:


import seaborn as sns
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
sns.distplot(prior_vec, hist=False, rug=False)
plt.title(label=r'Prior for $\alpha$', fontsize = 20)


# $\alpha \sim Gamma(1,10)$

# In[26]:


np.random.seed(1)
samples = 100
N = 100
a0 = 1
b0 = 10


# In[27]:


alpha_prior = np.random.gamma(a0,b0, size = 1)
prior_vec = alpha_prior
draws = DirPriorSample(samples = 1, N = N, alpha = alpha_prior, mean = 0, sd = 1)

results = draws[0]
results_save = np.zeros((results.shape[0],samples))

for i in range(samples):
    alpha_prior = np.random.gamma(a0,b0, size = 1)
    prior_vec = np.append(prior_vec, alpha_prior)
    draws = DirPriorSample(samples = 1, N = N, alpha = alpha_prior, mean = 0, sd = 1)
    new_results = draws[0]
    results_save[:,i] = new_results.flatten()

xx = np.append(np.linspace(-draws[1], draws[1], N+1),
              np.abs(draws[1]) + (np.abs(draws[1]) + np.abs(draws[1])  )/N )
# Plot the results
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.plot(xx, results_save[:,1], color = "gray", label = "Prior Draws")
plt.plot(xx, results_save, color = "grey",linewidth=0.5)
mean_fn = np.mean(results_save, axis = 1)
plt.plot(xx, mean_fn, color = "blue", label = "Estimated prior mean")
var_fn = np.var(results_save, axis = 1)
plt.plot(xx, var_fn, color = "green", label = "Estimated prior variance")
plt.plot(np.linspace(-4,4), norm.cdf(np.linspace(-4,4)), color = "black", label = "CDF of N(0,1)")
plt.legend(loc=1, bbox_to_anchor=(1, 0.5), fontsize = 15)
plt.title(label = r'$\alpha \sim Gamma(1,10)$', fontsize = 20)


# In[28]:


import seaborn as sns
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
sns.distplot(prior_vec, hist=False, rug=False)
plt.title(label=r'Prior for $\alpha$', fontsize = 20)


# In[ ]:




