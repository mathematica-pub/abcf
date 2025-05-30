data {
	int<lower=1> N;										// Number of practices
	real<lower=0> sdy;
	vector[N] y;										// outcome
	vector<lower=0, upper=1>[N] z;						// Treatment indicator
	vector<lower=0>[N] w;								// Number of benes/practice
	int<lower=1> K;										// Number of covariates
	matrix[N,K] X;										// Covariate matrix, standardized
	
	real<lower=0> ate_prior_sd;
	real<lower=0> sigu_hyperprior;
}

parameters {
	real alpha;
	vector[K] beta;
	real<lower=0> sigma_y_raw;
	real<lower=0> sigma_u_raw;
	vector[N] u_raw;
	real tau_bar_raw;
}

transformed parameters {
	vector[N] sigma_lik;
	
	sigma_lik = sigma_y_raw ./ sqrt(w);
}

model {
	alpha ~ std_normal();
	beta ~ std_normal();
	sigma_y_raw ~ std_normal();
	sigma_u_raw ~ normal(0, sigu_hyperprior);
	u_raw ~ std_normal();
	tau_bar_raw ~ normal(0, ate_prior_sd);
	
	y ~ normal(alpha + X * beta + z*tau_bar_raw + u_raw * sigma_u_raw, sigma_lik);
}

generated quantities {
	real sigma_y = sigma_y_raw * sdy;
	real sigma_u = sigma_u_raw * sdy;
	real tau_bar = tau_bar_raw * sdy;
	vector[N] u = u_raw * sigma_u_raw * sdy;
	real rho=0;
	real sigma_v=0;
}