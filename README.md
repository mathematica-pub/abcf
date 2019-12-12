# Bayesian Causal Forests

Welcome to the BCF site! This page provides more details on this implementation of Bayesian causal forests (BCF).  

## Why BCF?

BCF is a cutting-edge causal inference methodology that builds on Bayesian Additive Regression Trees (BART, Chipman, George, and McCulloch 2010).  BART and BCF both combine Bayesian regularization with regression trees to provide a highly flexible response surface that, thanks to the Bayesian regularizing priors, is not overfit to the training data.  BCF further extends BART's flexibility by specifying different models for relationships between covariates and the outcome and relationships between covariates and the treatment effect.  For more information, you can find the original BCF paper here: [https://arxiv.org/pdf/1706.09523.pdf](https://arxiv.org/pdf/1706.09523.pdf).

BCF performs remarkably well in simulation and has led the pack at recent rigorous causal inference competitions, such as those held at the Atlantic Causal Inference Conference. This implementation further extends existing BCF functionality by:

- allowing for heteroskedastic error
- automating multi-chain, multi-core implementations
- providing a suite of convergence diagnostic functions via the `{coda}` package
- accelerating some underlying computations, resulting in shorter runtimes

## Getting Started

If you are just getting started with BCF, we recommend starting with the tutorial vignette and the examples throughout the package documentation.

## Installation

This package requires compilation, so make sure you have Rtools properly installed -- details [here](https://cran.r-project.org/bin/windows/Rtools/).

Install the latest release from CRAN:

```{r}
install.packages("bcf")
```

Install the latest development version from GitHub:

```{r}
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jaredsmurray/bcf")
```

## Examples

```{r}
set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 1000
n_burn <- 1000
n_sim <- 1000

x <- matrix(rnorm(n*p), nrow=n)


# create targeted selection, whereby a unit's likelihood of joining the intervention (pi) is related to its expected outcome (mu)
q <- -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) -0.1

# generate treatment variable
pi <- pnorm(q)
z <- rbinom(n,1,pi)

# tau is the true treatment effect. It varies across units as a function of
# X3, the effect moderator
tau <-  1/(1 + exp(-x[,3]))

mu <- q

# generate the response using q, tau and z
y_noiseless <- mu + tau*z

# set the noise level relative to the expected mean function of Y
sigma <- diff(range(mu + tau*pi))/8

# draw the response variable with additive error
y <- y_noiseless + sigma*rnorm(n)

weights <- 1000.0*rep(1, n)

bcf_out <- bcf2::bcf(y            = y,
                 z                = z,
                 x_control        = x,
                 x_moderate       = x,
                 pihat            = pi,
                 nburn            = n_burn,
                 nsim             = n_sim,
                 w                = weights,
                 n_chains         = 4,
                 n_chain_clusters = 2,
                 random_seed      = 1,
                 update_interval  = 1)

```

```{r}
bcf2::summarise_bcf(bcf_out)

mean_square_error <- function (x,y){
  mean((x-y)^2)
}
```

```{r}
assess_closeness <- function(x,y, title){
  cat("Assessing Cloesness of ", title, "\n")
  
  mse = mean_square_error(x,y)
  
  print("MSE")
  print(mse)
  
  print("Error")
  print(sqrt(mse)/abs(mean(x)))
}

mu_correct = bcf_out$yhat - t(t(bcf_out$tau)*z)

print("Y Mean")
print(mean(y))

print("Tau Mean")
print(mean(tau))

print("mu Mean")
print(mean(mu))

assess_closeness(bcf_out$mu,mu_correct,'mu_compare')
```

## Latest package updates: 

### Weights

The original version of {bcf} does not allow for weights, which we often use in practical applications to account for heteroskedasticity. Where the original BCF model was specified as:

$y_i \sim N(\mu(x_i) + \tau(x_i) z_i, \sigma^2)$,

which assumes that all outcomes $y_i$ have the same variance $\sigma^2$, in the extended version we can relax this assumption to allow the variance to reflect the uncertainty in the $y_i$:

`y_i∼N(μ(x_i)+τ(x_i ) z_i,σ^2/w_i )`

We changed several parts of the code to incorporate these weights:

* Code to calculate sufficient statistics
* Code to update leaf node means
* Code to update variance across leaf node means

### Automating multichain processing

It is useful in Bayesian analysis to produce different runs of the same model, with different starting values, as a way of assessing convergence; if the different runs produce drastically different posterior distributions, it is a sign that the model has not converged fully.  In this version of {bcf} we have automated multichain processing and incorporated key MCMC diagnostics from the {coda} package, including effective sample sizes and the Gelman-Rubin statistic ("R hat").

*Question for Peter*: is the multichain processing also multicore?  If so, we should talk it up. 

### Speed-ups

Finally, our implementation parallelizes some steps of the sampling procedure to maximize efficiency in a multicore environment.  Our testing shows that these enhancements have reduced runtimes by [*Peter please fill in here*].
