functions {
  real kalman_log_lik_1fac_3ind(int T, array[] vector y, vector nu, real phi, real sigma_obs, real sigma_proc, vector lambda) {
    
    real inv_sigma_sq = 1.0 / square(sigma_obs);
    real l2 = lambda[2]; 
    real l3 = lambda[3];
    real nu1 = nu[1]; 
    real nu2 = nu[2]; 
    real nu3 = nu[3];
    
    real c = (1.0 + l2*l2 + l3*l3) * inv_sigma_sq;
    real log_det_Sigma = 6.0 * log(sigma_obs); 
    
    real phi_sq = square(phi);
    real sigma_proc_sq = square(sigma_proc);
    
    array[T] real pre_q;
    array[T] real pre_r;
    for (t in 1:T) {
      real dy1 = y[t, 1] - nu1;
      real dy2 = y[t, 2] - nu2;
      real dy3 = y[t, 3] - nu3;
      pre_q[t] = (dy1 + dy2*l2 + dy3*l3) * inv_sigma_sq;
      pre_r[t] = (dy1*dy1 + dy2*dy2 + dy3*dy3) * inv_sigma_sq;
    }
    
    real a = 0; 
    real P = 1; 
    real sum_ll_dynamic = 0;
    
    for (t in 1:T) {
      real qt = pre_q[t] - a * c;
      real rt = pre_r[t] - 2.0 * a * pre_q[t] + a * a * c;
      
      real F_det_term = 1.0 + P * c;
      real P_over_F = P / F_det_term;
      
      sum_ll_dynamic += log(F_det_term) + rt - P_over_F * qt * qt;
      
      real a_post = a + P_over_F * qt;
      a = phi * a_post; 
      P = phi_sq * P_over_F + sigma_proc_sq;
    }
    
    return -0.5 * (T * 3 * 1.8378770664093453 + T * log_det_Sigma + sum_ll_dynamic);
  }
}

data {
  int<lower=1> N_subj;
  int<lower=1> N_t;
  array[N_subj, N_t] vector[3] y;
}

parameters {
  vector[3] alpha2_nu;
  real<lower=0> tau_nu;
  real alpha2_phi;
  real<lower=0> tau_phi;
  real alpha2_log_sigma;
  real<lower=0> tau_log_sigma;
  real alpha2_log_psi;
  real<lower=0> tau_log_psi;
  vector<lower=0>[2] lambda_free;
  
  matrix[N_subj, 3] z_nu;
  vector[N_subj] z_phi; 
  vector[N_subj] z_log_sigma;
  vector[N_subj] z_log_psi;
}

transformed parameters {
  matrix[N_subj, 3] nu;
  vector[N_subj] phi;
  vector[N_subj] sigma;
  vector[N_subj] psi;
  vector[3] lambda = [1.0, lambda_free[1], lambda_free[2]]';
  
  for (u in 1:3) {
    nu[, u] = alpha2_nu[u] + tau_nu * z_nu[, u];
  }
  phi = tanh(alpha2_phi + tau_phi * z_phi); 
  sigma = exp(alpha2_log_sigma + tau_log_sigma * z_log_sigma);
  psi = exp(alpha2_log_psi + tau_log_psi * z_log_psi);
}

model {
  alpha2_nu[1] ~ normal(0, 5);
  alpha2_nu[2:3] ~ normal(0, 1);
  tau_nu ~ cauchy(0, 2); 
  
  alpha2_phi ~ normal(0, 1);
  tau_phi ~ cauchy(0, 1);
  
  alpha2_log_sigma ~ normal(0, 1);
  tau_log_sigma ~ normal(0, 1);
  
  alpha2_log_psi ~ normal(0, 1);
  tau_log_psi ~ normal(0, 1);
  
  lambda_free ~ normal(1, 2);
  
  to_vector(z_nu) ~ std_normal();
  z_phi ~ std_normal();
  z_log_sigma ~ std_normal();
  z_log_psi ~ std_normal();
  
  for (n in 1:N_subj) {
    target += kalman_log_lik_1fac_3ind(N_t, y[n], nu[n]', phi[n], sigma[n], psi[n], lambda);
  }
}
