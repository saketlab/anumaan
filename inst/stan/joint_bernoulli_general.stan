// joint_bernoulli_general.stan
//
// General D-class joint Bernoulli resistance-profile model with exact
// marginalization of untested classes. This is the same log-linear /
// saturated-interaction likelihood as the production J2-L model
// (anumaan-analysis/Scripts/stewardship/resistance_profiles/, Stan source
// j2_l_higher_order_location.stan), generalized from D = 4 fixed classes to
// any D, and from three separate pairwise/three-way/four-way blocks
// (V, G, W) to one concatenated higher-order design matrix V with a
// matching per-column prior-SD vector -- so adding or removing an
// antibiotic class changes D, S = 2^D, and P (the total interaction-term
// count) automatically, with no hand-edited blocks.
//
// Pair with anumaan::profile_build_design(), which builds every data
// element below from a wide S/R outcome matrix and a covariate design
// matrix.

functions {
  real partial_observed_lpmf(array[] int event_slice, int start, int end,
                             matrix X, matrix B, matrix V,
                             matrix beta, vector eta,
                             array[] int n_compatible,
                             array[,] int compatible_state) {
    real lp = 0;
    int S = rows(B);
    for (i in 1:size(event_slice)) {
      int e = event_slice[i];
      vector[S] score = B * (X[e] * beta)' + V * eta;
      vector[n_compatible[e]] compatible_score;
      for (j in 1:n_compatible[e]) compatible_score[j] = score[compatible_state[e, j]];
      lp += log_sum_exp(compatible_score) - log_sum_exp(score);
    }
    return lp;
  }
}
data {
  int<lower=1> N_events;              // events
  int<lower=2> D;                     // antibiotic classes chosen for the panel
  int<lower=1> K;                     // covariate design columns (caller's choice)
  int<lower=4> S;                     // = 2^D, joint resistance states
  int<lower=0> P;                     // total interaction terms, order 2..max_order
                                       //   = sum_{k=2}^{max_order} choose(D, k)
  matrix[N_events, K] X;              // covariate design (main effects)
  matrix[S, D] B;                     // main-effect state design (+-1 per class)
  matrix[S, P] V;                     // concatenated higher-order state design
  array[N_events] int<lower=1, upper=S> n_compatible;
  array[N_events, S] int<lower=1, upper=S> compatible_state; // only first n_compatible[e] used
  int<lower=1> grainsize;
  vector<lower=0>[K] beta_sd;         // prior SD per covariate column
  vector<lower=0>[P] eta_sd;          // prior SD per interaction term (set by order)
}
transformed data {
  array[N_events] int event_index;
  for (e in 1:N_events) event_index[e] = e;
}
parameters {
  matrix[K, D] beta;   // D main-effect coefficient vectors, one per class
  vector[P] eta;       // every global interaction scalar, orders 2..max_order concatenated
}
model {
  for (k in 1:K) beta[k] ~ normal(0, beta_sd[k]);
  eta ~ normal(0, eta_sd);
  target += reduce_sum(partial_observed_lpmf, event_index, grainsize,
                       X, B, V, beta, eta, n_compatible, compatible_state);
}
generated quantities {
  array[N_events] real log_lik;
  for (e in 1:N_events) {
    vector[S] score = B * (X[e] * beta)' + V * eta;
    vector[n_compatible[e]] compatible_score;
    for (j in 1:n_compatible[e]) compatible_score[j] = score[compatible_state[e, j]];
    log_lik[e] = log_sum_exp(compatible_score) - log_sum_exp(score);
  }
}
