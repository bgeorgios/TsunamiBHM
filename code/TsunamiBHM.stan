/*      ***********************
	© Georgios Boumis, 2023
	***********************
	
	Contact: gboumis@crimson.ua.edu
	Website: www.gboumis.com
*/

functions{
	// Function to return Tapered Pareto log-likelihood for log(b) and log(ac) parameterization
	real tp_lpdf(vector y, real u, real b, real ac) {
		vector[rows(y)] lp;
		real loglike;
		int N;
		
		N = rows(y);
		
		for(i in 1:N){
			lp[i] = log(exp(b) / y[i] + 1 / exp(ac));
		}
		loglike = sum(lp) + exp(b) * N * log(u) - exp(b) * sum(log(y)) + (u * N) / exp(ac) - 1 / exp(ac) * sum(y);
		return(loglike);
	}
	
	// Function to return exponential covariance kernel 
	matrix kernel(int M, vector dist, real alpha, real rho) {
		matrix[M, M] K;
		
		for (i in 1:M){
			for (j in 1:M) {
				K[i, j] = pow(alpha, 2) * exp(-fabs(dist[i] - dist[j]) / rho);
			}
		}
		return K;
	}
}

data {
	int<lower=0> N; // total number of observed tsunami amplitudes
	int<lower=0> M; // number of tide gauges
	int<lower=0> persite[M]; // number of tsunami amplitudes per site
	vector[N] y; // tsunami amplitudes grouped from all tide gauges altogether
	matrix[M, 4] x; // covariates, i.e., intercept (should be 1.0) + latitude + longitude + continental shelf width
	vector[M] dist; // distance of each tide gauge to the baseline "San Diego, CA" tide gauge station
	vector[M] u; // threshold of record completeness for each tide gauge
}

parameters{
	// hyperparameters used for mean functions and covariance matrices of Gaussian processes	
	real w0b;
	real w1b;
	real w2b;
	real w3b;
	
	real<lower=0> ab;
	real<lower=0> rhob;
	
	real w0ac;
	real w1ac;
	real w2ac;
	real w3ac;
	
	real<lower=0> aac;
	real<lower=0> rhoac;
	
	// standard-normal parameters
	vector[M] zb;
	vector[M] zac;
}

transformed parameters{
	// mean function for b
	vector[M] gfb = w0b * x[, 1] + w1b * x[, 2] + w2b * x[, 3] + w3b * x[, 4]; 
	
	// mean function for ac
	vector[M] gfac = w0ac * x[, 1] + w1ac * x[, 2] + w2ac * x[, 3] + w3ac * x[, 4];
	
	// covariance matrix for b and its Cholesky decomposition
	matrix[M, M] Kb = kernel(M, dist, ab, rhob);
	matrix[M, M] Lb = cholesky_decompose(Kb);
	
	// covariance matrix for ac and its Cholesky decomposition
	matrix[M, M] Kac = kernel(M, dist, aac, rhoac);
	matrix[M, M] Lac = cholesky_decompose(Kac);
	
	// non-centered parameterization
	vector[M] b = gfb + Lb * zb;
	vector[M] ac = gfac + Lac * zac;
	
	// Note*: The Cholesky decomposition (L) of the covariance matrix (K) implies that L*L(^T) = K
	// Note**: b ~ multi_normal(gfb, Kb) is equivalent to b = gfb + Lb*zb, where zb ~ normal(0, 1)
	// Note***: Non-centered parameterization is necessary to achieve efficient sampling and convergence due to small amount of data!
}

model {
	int idx;
	
	// priors for wi hyperparameters of the Gaussian process for b
	w0b ~ normal(0, 5);
	w1b ~ normal(0, 5);
	w2b ~ normal(0, 5);
	w3b ~ normal(0, 5);
	
	// priors for rho and alpha hyperparameters of the Gaussian process for b
	ab ~ gamma(1.25, 0.25);
	rhob ~ gamma(112.5, 0.375); // secondary: gamma(36, 0.12)
	
	zb ~ normal(0, 1); // standard-normal parameter used for non-centered parameterization
	
	// // priors for wi hyperparameters of the Gaussian process for ac
	w0ac ~ normal(0, 5);
	w1ac ~ normal(0, 5);
	w2ac ~ normal(0, 5);
	w3ac ~ normal(0, 5);
	
	// priors for rho and alpha hyperparameters of the Gaussian process for ac
	aac ~ gamma(1.25, 0.25);
	rhoac ~ gamma(112.5, 0.375);  // secondary: gamma(36, 0.12)
	
	zac ~ normal(0, 1); // standard-normal parameter used for non-centered parameterization
	
	// increment log-likelihood
	idx = 1;
	for (i in 1:M) {
		segment(y, idx, persite[i]) ~ tp(u[i], b[i], ac[i]);
		idx = idx + persite[i];
	}
}
