# include <RcppArmadillo.h>
# include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Compares how similar two points are with regards to their clustering across 
// all iterations.
// Works for unsupervised methods (i.e. allows label flipping)
double point_similarity(arma::uword point, 
                        arma::uword comparison_point,
                        arma::umat cluster_record,
                        arma::uword num_iter) {
  double out = 0.0;
  
  for (arma::uword i = 0; i < num_iter; i++){
    if(cluster_record(point, i) == cluster_record(comparison_point, i)){
      out++;
    }
    
  }
  out = out / num_iter;
  return out;
}

// Constructs a similarity matrix comparing all points clustering across the 
// iterations
// [[Rcpp::export]]
arma::mat similarity_mat(arma::umat cluster_record){
  arma::uword sample_size = cluster_record.n_rows;
  arma::uword num_iter = cluster_record.n_cols;
  arma::mat out(sample_size, sample_size);
  
  for (arma::uword point = 0; point < sample_size; point++){ // if not doing diagonal, restrict to sample size - 1
    for (arma::uword comparison_point = point; // + 1; 
         comparison_point < sample_size;
         comparison_point++){
      out(point, comparison_point) = point_similarity(point, 
          comparison_point,
          cluster_record,
          num_iter);
      out(comparison_point, point) = out(point, comparison_point);
    }
  }
  return out;
}

// Calculate the entropy for the current cluster weights
// [[Rcpp::export]]
double entropy(arma::vec class_weights){
  arma::uword n = class_weights.n_elem;
  arma::vec entropy_components(n);
  // std::cout << "\nDeclared\n";
  
  
  for(arma::uword i = 0; i < n; i++){
    entropy_components(i) = - class_weights(i) * log(class_weights(i));
    if (entropy_components.has_nan()){
      entropy_components(i) = 0;
    }
  }
  // std::cout << "Inter";
  double entropy_out = sum(entropy_components);
  return entropy_out;
}

// === Gaussian ================================================================

// Returns a variable involved in updating the scale, it is similar to sample 
// covariance
arma::mat S_n(arma::mat data,
              arma::vec sample_mean,
              arma::uword sample_size,
              arma::uword num_cols){
  arma::mat sample_covariance(num_cols, num_cols);
  sample_covariance.zeros();
  if(sample_size > 0){
    for(arma::uword i = 0; i < sample_size; i++){
      
      arma::vec row_i = trans(data.row(i));
      
      sample_covariance = (sample_covariance 
                             + ((row_i - sample_mean) 
                                  * arma::trans(row_i - sample_mean)
                             )
      );
      
    }
  }
  return sample_covariance;
}

// Returns a vector of the mean of a cluster of size n
arma::vec mean_n(double lambda_0,
                 arma::vec mu_0,
                 arma::uword sample_size,
                 arma::uword num_cols,
                 arma::vec sample_mean){
  arma::vec mu_n(num_cols);
  mu_n = ((lambda_0 * mu_0 + sample_size * sample_mean)
            / (lambda_0 + sample_size));
  return mu_n;
}

// Returns the matrix for the scale hyperparameter for n observations
arma::mat scale_n(arma::mat scale_0,
                  arma::vec mu_0,
                  double lambda_0,
                  arma::mat sample_covariance,
                  arma::uword sample_size,
                  arma::vec sample_mean){
  arma::mat scale_out;
  
  
  if(sample_size > 0){
    scale_out = (scale_0
                   + sample_covariance
                   + ((lambda_0 * sample_size) / (lambda_0 + sample_size))
                   * (sample_mean - mu_0) * arma::trans(sample_mean - mu_0)
    );
    return scale_out;
  }
  
  scale_out = scale_0; 
  return scale_out;
}

// sample a multivariate normal
arma::mat mvrnormArma(arma::uword n,
                      arma::vec mu,
                      arma::mat sigma) {
  arma::uword ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() +
    Y * arma::chol(sigma);
}

// sample from the mean posterior
arma::vec mean_posterior(arma::vec mu_0,
                         arma::mat variance,
                         double lambda_0,
                         arma::mat data){
  arma::uword ncols = data.n_cols;
  arma::vec sample_mean(ncols); sample_mean.zeros();
  arma::vec mu_out(ncols);
  arma::vec mu_n(ncols);
  arma::mat variance_n(ncols, ncols);
  
  arma::uword sample_size = data.n_rows;
  if (sample_size > 0){
    sample_mean = trans(arma::mean(data, 0));
  }
  
  // std::cout << "\nPast initial if\n";
  
  double lambda_n = lambda_0 + sample_size;
  
  mu_n = mean_n(lambda_0, mu_0, sample_size, ncols, sample_mean);
  
  // std::cout << "\nPast initial mu_n\n";
  // std::cout << "\n" << mu_n << "\n";
  
  variance_n = variance / lambda_n;
  
  // std::cout << "\n" << variance_n << "\n";
  
  arma::mat x = mvrnormArma(1, mu_n, variance_n);
  
  mu_out = arma::conv_to<arma::vec>::from(x.row(0));
  
  return mu_out;
  
}


// Sample the variance for after n observations
arma::mat variance_posterior(int df_0,
                             arma::mat scale_0,
                             double lambda_0,
                             arma::vec mu_0,
                             arma::mat data){
  
  arma::uword sample_size = data.n_rows, num_cols = data.n_cols;
  
  // std::cout << "\nCluster data:\n" << data;
  
  arma::vec sample_mean(num_cols); sample_mean.zeros();
  
  int df_n = df_0 + sample_size;
  
  // std::cout << "Reached another dec\n";
  
  arma::mat scale_n_value(num_cols, num_cols);
  arma::mat sample_covariance(num_cols, num_cols);
  arma::mat variance_out(num_cols, num_cols);
  
  // std::cout << sample_size << "\n";
  
  if (sample_size > 0){
    
    sample_mean = arma::trans(arma::mean(data, 0));
    
  } else{
    sample_mean.fill(0.0);
  }
  
  // std::cout << "\nSample covariance reached\n";
  
  sample_covariance = S_n(data, sample_mean, sample_size, num_cols);
  
  // std::cout << "Scale_n reached\n";
  
  arma::mat samp_cov(num_cols, num_cols);
  samp_cov = (sample_size - 1) * arma::cov(data);
  
  // std::cout  << samp_cov << "\n\n";
  
  scale_n_value = scale_n(
    scale_0,
    mu_0,
    lambda_0,
    sample_covariance,
    sample_size,
    sample_mean
  );
  
  
  // std::cout << scale_n_value << "\n";
  
  variance_out = arma::iwishrnd(scale_n_value, df_n);
  
  return variance_out;
  
} 

// Returns a 2-D field of cubes (i.e. a 5D object) for the current iterations 
// means and variances across each cluster (hence a field)
arma::field<arma::cube> mean_variance_sampling(arma::mat data,
                                               arma::uvec cluster_labels,
                                               arma::uword k,
                                               int df_0,
                                               arma::uword num_cols,
                                               arma::mat scale_0,
                                               double lambda_0,
                                               arma::vec mu_0){
  arma::mat cluster_data;
  arma::mat variance(num_cols, num_cols);
  arma::vec mu(num_cols);
  
  arma::field<arma::cube> mean_variance_field(2);
  
  arma::cube var_entry = arma::zeros<arma::cube>(num_cols, num_cols, k);
  arma::cube mu_entry = arma::zeros<arma::cube>(num_cols, 1, k);
  
  mean_variance_field(0) = var_entry;
  mean_variance_field(1) = mu_entry;
  
  for (arma::uword j = 1; j < k + 1; j++) {
    // std::cout << "\nj for loop";
    cluster_data = data.rows(find(cluster_labels == j ));
    
    // std::cout << mean_variance_field(1).slice(j - 1) << "\n";
    
    // std::cout << "\n" << df_0 << "\n";
    
    // if(cluster_data.n_rows < 10){
      // std::cout << cluster_data << "\n\n";
    // }
    
    mean_variance_field(0).slice(j - 1) = variance_posterior(
      df_0,
      scale_0,
      lambda_0,
      mu_0,
      cluster_data
    );
    
    // std::cout << "\nVariance sampled\n";
    
    // std::cout << mean_variance_field(1).slice(j - 1) << "\n";
    
    mean_variance_field(1).slice(j - 1) = mean_posterior(mu_0, 
                        mean_variance_field(0).slice(j - 1), 
                        lambda_0,
                        cluster_data);
    
    // std::cout << "\nAccessed cubes";
  }
  return mean_variance_field;
}

// === Dirichlet ===============================================================

// update the concentration parameter in the Dirichlet distribution
arma::vec concentration_n(arma::vec concentration_0,
                          arma::uvec cluster_labels,
                          arma::uword num_cat){
  
  // arma::uword n = cluster_labels.n_elem;
  arma::uvec class_members;
  arma::uword class_count;
  arma::vec concentration(num_cat);
  
  for (arma::uword i = 1; i < num_cat + 1; i++) {
    class_count = 0;
    class_members = find(cluster_labels == i);
    class_count = class_members.n_elem;
    // for (arma::uword j = 0; j < n; j++ ) {
    //   if (cluster_labels(j) == i) {
    //     class_count++;
    //   }
    // }
    
    concentration(i - 1) = arma::as_scalar(concentration_0(i - 1)) + class_count;
  }
  return concentration;
}

// sample parameters for a dirichlet distribution (normally for the clusters)
// [[Rcpp::export]]
arma::vec dirichlet_posterior(arma::vec concentration_0,
                              arma::uvec cluster_labels,
                              arma::uword num_clusters){
  arma::vec cluster_weight = arma::zeros<arma::vec>(num_clusters);
  
  arma::vec concentration(num_clusters);
  concentration = concentration_n(concentration_0,
                                  cluster_labels,
                                  num_clusters);
  
  
  for (arma::uword i = 0; i < num_clusters; i++) {
    
    // cluster_weight(i) = Rf_rgamma(arma::as_scalar(concentration(i)), 1);
    cluster_weight(i) = arma::randg( arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
  }
  
  double total_cluster_weight = sum(cluster_weight);
  cluster_weight = cluster_weight / total_cluster_weight;
  return cluster_weight;
}

//  TEST FUNCTIONS FOR CLASSES WITH NULL OF 0 //////////////////////////////////

arma::vec concentration_n_class(arma::vec concentration_0,
                          arma::uvec cluster_labels,
                          arma::uword num_cat){
  
  // arma::uword n = cluster_labels.n_elem;
  
  arma::uword class_count = 0;
  arma::vec concentration(num_cat);
  
  arma::uvec class_members;
  
  for (arma::uword i = 0; i < num_cat; i++) {
    // class_count = 0;
    class_members = find(cluster_labels == i);
    class_count = class_members.n_elem;
    // for (arma::uword j = 0; j < n; j++ ) {
    //   if (cluster_labels(j) == i) {
    //     class_count++;
    //   }
    // }
    
    concentration(i) = arma::as_scalar(concentration_0(i)) + class_count;
  }
  return concentration;
}


arma::vec dirichlet_posterior_class(arma::vec concentration_0,
                                    arma::uvec cluster_labels,
                                    arma::uword num_clusters){
  arma::vec cluster_weight = arma::zeros<arma::vec>(num_clusters);
  
  arma::vec concentration(num_clusters);
  concentration = concentration_n_class(concentration_0,
                                        cluster_labels,
                                        num_clusters);
  
  
  for (arma::uword i = 0; i < num_clusters; i++) {
    
    // cluster_weight(i) = Rf_rgamma(arma::as_scalar(concentration(i)), 1);
    // double a = arma::as_scalar(concentration(i));
    cluster_weight(i) = arma::randg(arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
  }
  
  double total_cluster_weight = sum(cluster_weight);
  cluster_weight = cluster_weight / total_cluster_weight;
  return cluster_weight;
}

// Several functions to initialise the 3D array required for the different
// classes for each variable within each cluster

// count unique entries in a vector
arma::uword unique_counter(arma::uvec v){
  std::sort(v.begin(), v.end());
  arma::uword unique_count = std::unique(v.begin(), v.end()) - v.begin();
  return unique_count;
}

// returns a vector of the number of unqiue values in each column
// [[Rcpp::export]]
arma::uvec cat_counter(arma::umat data){
  arma::uword num_cols = data.n_cols;
  arma::uvec num_categories(num_cols);
  for(arma::uword i = 0; i < num_cols; i++){
    num_categories(i) = unique_counter(data.col(i));
  }
  return num_categories;
}

// find the number of categories in each covariate and declare the appropriate
// matrix to record the associated probabilties for each cluster
// [[Rcpp::export]]
arma::field<arma::mat> declare_class_probs_field(arma::uvec cat_count,
                                                 arma::uword num_cols,
                                                 arma::uword num_clusters){
  
  arma::field<arma::mat> class_probabilities(num_cols);
  for(arma::uword i = 0; i < num_cols; i++){
    arma::mat phi_j(num_clusters, cat_count(i));
    phi_j.zeros();
    class_probabilities(i) = phi_j;
  }
  return class_probabilities;
}

// Sample the probabilities for each category across all clusters for each covariate
// [[Rcpp::export]]
arma::field<arma::mat> sample_class_probabilities(arma::umat data,
                                                  arma::field<arma::mat> class_probabilities,
                                                  arma::field<arma::vec> phi_prior,
                                                  arma::uvec cluster_labels,
                                                  arma::uvec cat_count,
                                                  arma::uword num_clusters,
                                                  arma::uword num_cols
){

  arma::umat cluster_data;
  // arma::umat indices;
  for(arma::uword k = 1; k < num_clusters + 1; k++){
    
    // std::cout << "In for loop, k = " << k << "\n";
    
    // indices = find(cluster_labels == k);
    
    // std::cout << "\n" << indices << "\n";
    
    cluster_data = data.rows(find(cluster_labels == k));
    
    // std::cout << "Generic message\n";
    
    for(arma::uword j = 0; j < num_cols; j++){
      
      
      class_probabilities(j).row(k - 1) = arma::trans(dirichlet_posterior_class(phi_prior(j),
                                                                            cluster_data.col(j),
                                                                            cat_count(j)
                                                                            )
      );
      
      // std::cout << "Another message\n";
    }
  }
  return class_probabilities;
}

// Sample the cluster membership of point
// [[Rcpp::export]]
arma::vec categorical_cluster_probabilities(arma::urowvec point,
                                            arma::umat data,
                                            arma::field<arma::mat> class_probabilities,
                                            arma::vec cluster_weights,
                                            arma::uword num_clusters,
                                            arma::uword num_cols){
  
  // std::cout << "In function cluster_probabilities\n";
  arma::vec probabilities = arma::zeros<arma::vec>(num_clusters);

  double curr_weight = 0.0;
  
  // std::cout << "\n\n" << class_probabilities << "\n\n";
  for(arma::uword i = 0; i < num_clusters; i++){
    curr_weight = log(cluster_weights(i));
    for(arma::uword j = 0; j < num_cols; j++){

      probabilities(i) = probabilities(i) + std::log(class_probabilities(j)(i, point(j)));
    }
    probabilities(i) = probabilities(i) + curr_weight;
  }
  probabilities = exp(probabilities - max(probabilities));
  probabilities = probabilities / sum(probabilities);
  
  return probabilities;
}

// Predicts the cluster assignments based on a vector of probabilities using
// the rejection method
// [[Rcpp::export]]
arma::uword cluster_predictor(arma::vec probabilities){
  double u;
  arma::uword pred;
  u = arma::randu<double>( );
  
  // include + 1 if labels centred on 1
  pred = 1 + sum(u > cumsum(probabilities)); // 1 + 
  return pred;
}

// The actual categorical clustering all wrapped up in one function
// [[Rcpp::export]]
Rcpp::List categorical_clustering(arma::umat data,
                                 arma::field<arma::vec> phi_prior,
                                 arma::uvec cluster_labels,
                                 arma::vec fix_vec,
                                 arma::vec cluster_weight_priors,
                                 arma::uword num_clusters,
                                 arma::uword num_iter,
                                 arma::uword burn,
                                 arma::uword thinning){
  
  arma::uword n = data.n_rows;
  arma::uword num_cols = data.n_cols;
  
  arma::uvec cat_count(num_cols);
  cat_count = cat_counter(data);
  
  arma::field<arma::mat> class_probabilities(num_cols);
  
  // std::cout << "Declaring field of matrices for class probs\n";
  
  class_probabilities = declare_class_probs_field(cat_count,
                            num_cols,
                            num_clusters);
  
  // std::cout << "Class probabilities declared\n";
  
  // std::cout << "Num clusters: " << num_clusters << "\n";
  
  arma::vec cluster_weights(num_clusters);
  
  arma::vec curr_cluster_probs(num_clusters);
  
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  // std::cout << "N: " << n << "\nEff count: " << eff_count << "\n";
  arma::umat record(n, eff_count);
  record.zeros();
  
  // std::cout << "Reached for loop\n";
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    cluster_weights = dirichlet_posterior(cluster_weight_priors,
                                          cluster_labels,
                                          num_clusters);
    
    // std::cout << "Cluster weights calculated\n";
    
    class_probabilities = sample_class_probabilities(data,
                                                     class_probabilities,
                                                     phi_prior,
                                                     cluster_labels,
                                                     cat_count,
                                                     num_clusters,
                                                     num_cols
                                 
    );
    
    // std::cout << "Class probs calculated\n";
    // std::cout << class_probabilities << "\n\n";
    
    for(arma::uword j = 0; j < n; j++){
      
      // sample cluster for each point here
      curr_cluster_probs = categorical_cluster_probabilities(data.row(j),
                                               data,
                                               class_probabilities,
                                               cluster_weights,
                                               num_clusters,
                                               num_cols);
      
      // std::cout << "Cluster sampled\n";
      
      if(! fix_vec(j)){
        cluster_labels(j) = cluster_predictor(curr_cluster_probs);
      }
    }
    if (i >= burn && (i - burn) % thinning == 0) {
      // std::cout << "\n" << i << "\n";
      record.col((i - burn) / thinning) = cluster_labels;
    }
  }
  arma::mat sim(n, n);
  sim = similarity_mat(record);
  return List::create(Named("similarity") = sim,
                      Named("class_record") = record);
  // return sim;
}

// === Gaussian clustering =====================================================

arma::vec sample_gaussian_cluster(arma::vec point,
                                  arma::mat data,
                                  arma::uword k,
                                  arma::vec class_weights,
                                  arma::uvec class_labels,
                                  arma::cube mu,
                                  arma::cube variance,
                                  bool outlier = false,
                                  arma::vec global_mean = arma::zeros<arma::vec>(1),
                                  arma::mat global_variance = arma::zeros<arma::mat>(1, 1),
                                  double t_df = 4.0){
  
  double curr_weight;
  double exponent;
  double log_likelihood;
  
  
  double log_det;
  arma::vec prob_vec(k);
  
  arma::uvec count_probs;
  
  arma::uword d = data.n_cols;;
  
  for(arma::uword i = 1; i < k + 1; i++){
    curr_weight = log(class_weights(i - 1));
    
    if(outlier && i == k){
      
      exponent = arma::as_scalar(arma::trans(point - global_mean) 
                                   * arma::inv(global_variance)
                                   * (point - global_mean));
                                   
                                   log_det = arma::log_det(global_variance).real();
                                   
                                   log_likelihood = (lgamma((t_df + d)/2.0) - lgamma(t_df/2.0)
                                                       + d/2.0 * log(t_df * M_PI) + log_det - ((t_df + d)/2.0) *
                                                         log(1 + (1/t_df) * exponent));
                                   
    }
    else {
      exponent = arma::as_scalar(arma::trans(point - mu.slice(i - 1)) 
                                   * arma::inv(variance.slice(i - 1))
                                   * (point - mu.slice(i - 1)));
                                   
                                   log_det = arma::log_det(variance.slice(i - 1)).real();
                                   log_likelihood = -0.5 *(log_det + exponent + d * log(2 * M_PI));
    }
    prob_vec(i - 1) = curr_weight + log_likelihood;
  } 
  prob_vec = exp(prob_vec - max(prob_vec));
  prob_vec = prob_vec / sum(prob_vec);
  
  return prob_vec;
}

// The actual clustering/sampling
// [[Rcpp::export]]
Rcpp::List gaussian_clustering(arma::uword num_iter,
                               arma::vec concentration_0,
                               arma::mat scale_0,
                               arma::uvec class_labels,
                               std::vector<bool> fix_vec,
                               arma::vec mu_0,
                               double lambda_0,
                               arma::mat data,
                               int df_0,
                               arma::uword k,
                               arma::uword burn,
                               arma::uword thinning,
                               bool outlier = false,
                               double t_df = 4.0,
                               bool record_posteriors = false,
                               bool normalise = false){
  
  // std::cout << "In function";
  arma::uword N;
  arma::uword num_cols;
  N = data.n_rows;
  num_cols = data.n_cols;
  
  //  Normalise the continuous data
  if(normalise){
    data = arma::normalise(data);
  }
  
  // for use in the outlier distribution
  arma::mat global_variance(num_cols, num_cols);
  
  // std::cout << arma::cov(data) << "\n";
  
  global_variance = 0.5 * arma::cov(data); // Olly's rec
  
  arma::vec global_mean(num_cols);
  
  global_mean = arma::trans(arma::mean(data, 0));
  
  arma::vec entropy_cw(num_iter);
  
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  arma::uword record_ind;
  
  arma::umat record(N, eff_count);
  record.zeros();
  
  // std::cout << "Record out \n";
  
  arma::mat sim(N, N);
  arma::mat cluster_data;
  
  // Add +1 to k to allow outlier class
  if(outlier){
    k++;
  }
  
  arma::vec class_weights(k);
  
  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration
  ListMatrix variance(eff_count, k);
  ListMatrix mu(eff_count, k);
  
  arma::vec point;
  
  // std::cout << "Faux output sentence\n";
  
  arma::cube class_probs(eff_count, k, N);
  arma::vec curr_class_probs(k);
  
  // std::cout << "Output sentence\n";
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::field<arma::cube> loc_mu_variance(2);
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    class_weights = dirichlet_posterior(concentration_0, class_labels, k);
    
    // std::cout << class_weights << "\n\n";
    // std::cout << "\nENTROPY";
    
    entropy_cw(i) = entropy(class_weights);
    
    // std::cout << "\nBegin sampling parameters\n";
    
    loc_mu_variance = mean_variance_sampling(data,
                                             class_labels,
                                             k,
                                             df_0,
                                             num_cols,
                                             scale_0,
                                             lambda_0,
                                             mu_0);
    
    // std::cout << "\nAccessed cubes\n";
    
    for (arma::uword jj = 0; jj < N; jj++){
      // if the current point is not fixed, sample its class
      point = arma::trans(data.row(jj));
      
      curr_class_probs = sample_gaussian_cluster(point, 
                                                 data,
                                                 k, 
                                                 class_weights, 
                                                 class_labels,
                                                 loc_mu_variance(1),
                                                 loc_mu_variance(0),
                                                 outlier,
                                                 global_mean,
                                                 global_variance,
                                                 t_df
      );
      
      // std::cout << curr_class_probs << "\n\n";
      
      if (i >= burn && (i - burn) % thinning == 0) {
        // std::cout << "record accessed" << "\n";
        record_ind = (i - burn) / thinning;
        class_probs.slice(jj).row(record_ind) = arma::trans(curr_class_probs);
        
      }
      if(! fix_vec[jj]){
        class_labels(jj) = cluster_predictor(curr_class_probs);
        // std::cout << "New label\n" << class_labels(jj) << "\n";
      }
    }
    // std::cout << "Labels\n" << class_labels << "\n";
    // std::cout << "Generic message\n" << "\n";
    
    if (i >= burn && (i - burn) % thinning == 0) {
      // std::cout << i << "\n";
      
      record_ind = (i - burn) / thinning;
      record.col(record_ind) = class_labels;
      // std::cout << "record accessed" << "\n";
      
      if(record_posteriors){
        for(arma::uword j = 0; j < k; j ++){
          
          // std::cout << "Recording params" << j << "\n";
          
          mu(record_ind, j) = loc_mu_variance(1).slice(j);
          variance(record_ind, j) = loc_mu_variance(0).slice(j);
        }
      }
      
    }
  }
  
  // std::cout << "Issue is here";
  sim = similarity_mat(record);
  // std::cout << "Y is here";
  // return sim;
  
  if(record_posteriors){
    
    return List::create(Named("similarity") = sim,
                        Named("class_record") = record,
                        Named("mean_posterior") = mu,
                        Named("variance_posterior") = variance,
                        Named("entropy") = entropy_cw,
                        Named("class_prob") = class_probs);
  }
  
  return List::create(Named("similarity") = sim,
                      Named("class_record") = record,
                      Named("entropy") = entropy_cw,
                      Named("class_prob") = class_probs);
}


// === MDI Code ================================================================

// returns the count of the number of points with the same label in both contexts
arma::uword count_common_cluster(arma::uvec cluster_labels_1,
                                 arma::uvec cluster_labels_2,
                                 arma::uword n){
  arma::uword common_cluster = 0;
  for(arma::uword i = 0; i < n; i++){
    if(cluster_labels_1(i) == cluster_labels_2(i)){
      common_cluster++;
    }
  }
  return common_cluster;
}

// Calculates the rate based on the data
double observed_rate(double v,
                     arma::uword min_num_clusters,
                     arma::vec cluster_weights_gaussian,
                     arma::vec cluster_weights_categorical){
  // declare the rate
  double b = 0.0;
  
  // the rate is proportional to the sum of the product of cluster weights 
  // across contexts (i.e. the product of cluster_weight_i in context 1 times
  // cluster_weight_i in context 2)
  for(arma::uword i = 0; i < min_num_clusters; i++){
    b = b + cluster_weights_gaussian(i) * cluster_weights_categorical(i);
  }
  
  // the rate is equal to this sum times the strategic latent variable, v
  b = v * b;
  
  return b;
}

// Calculate the factorial
int factorial(arma::uword n)
{
  if(n <= 1){
    return 1;
  }
  return n * factorial(n - 1);
  
  // return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


// Returns the log of the factorial of n
double log_factorial(arma::uword n){
  if(n <= 1){
    return 0;
  }
  return log(n) + log_factorial(n - 1);
}

// samples a gamma distn for the current iterations context similarity parameter
// (phi in the original 2012 paper).
double sample_phi(arma::uvec cluster_labels_1,
                  arma::uvec cluster_labels_2,
                  arma::vec cluster_weights_1,
                  arma::vec cluster_weights_2,
                  double v,
                  arma::uword n,
                  arma::uword min_num_clusters,
                  double a0,
                  double b0){
  
  // calculate the shape of the relevant gamma function
  arma::uword count_same_cluster = count_common_cluster(cluster_labels_1,
                                                        cluster_labels_2,
                                                        n);
  
  arma::vec prob_vec(count_same_cluster);
  
  // calculate the rate
  double b = observed_rate(v,
                           min_num_clusters,
                           cluster_weights_1,
                           cluster_weights_2)
    + b0;
  
  double phi = 0.0;
  
  // context similarity is a weighted sum of gammas
  if(count_same_cluster > 0){
    
    for(arma::uword i = 0; i < count_same_cluster; i++){
      prob_vec(i) = log_factorial(count_same_cluster)
      - log_factorial(i)
      - log_factorial(count_same_cluster - i)
      + log_factorial(i + a0 - 1)
      - (i + a0)*log(b);
    }
    
    prob_vec = exp(prob_vec - max(prob_vec));
    prob_vec = prob_vec / sum(prob_vec);
    
    // If move to general null value of 0 can use the prediction function
    double u;
    u = arma::randu<double>( );
    
    // include + 1 if labels centred on 1
    arma::uword pred_ind = 0;
    
    pred_ind = sum(u > cumsum(prob_vec)); // 1 +
    
    // phi_12 = arma::randg( arma::distr_param(count_same_cluster + a0, 1/b) );
    phi = arma::randg( arma::distr_param(pred_ind + a0, 1/b) );
    
  }
  
  return phi;
}


// Calculates the normalising constant for the posterior
double calculate_normalising_constant(arma::vec cluster_weights_1,
                                      arma::vec cluster_weights_2,
                                      double context_similarity,
                                      arma::uword num_clusters_1,
                                      arma::uword num_clusters_2){
  double Z = 0.0;
  
  for(arma::uword i = 0; i < num_clusters_1; i++){
    for(arma::uword j = 0; j < num_clusters_2; j++){
      Z = Z + cluster_weights_1(i) * cluster_weights_2(j) *
        (1 + context_similarity * (i == j));
    }
  }
  return Z;
}

// Sample the cluster membership of a categorical sample for MDI
arma::vec mdi_cat_clust_prob(arma::uword row_index,
                             arma::umat data,
                             arma::field<arma::mat> class_probs,
                             arma::uword num_clusters,
                             arma::uword num_cols_cat,
                             double phi,
                             arma::vec cluster_weights_categorical,
                             arma::uvec clust_labels,
                             arma::uvec clust_labels_comp){
  
  // cluster_labels_comparison is the labels of the data in the other context
  
  arma::vec prob_vec = arma::zeros<arma::vec>(num_clusters);
  arma::urowvec point = data.row(row_index);
  
  // pretty much this is the product of probabilities possibly up-weighted by
  // being in the same cluster in a different context and weighted by the cluster
  // weight in the current context
  double curr_weight = 0.0;
  
  // Upweight for similarity of contexts
  double similarity_upweight = 0.0;
  
  arma::uword common_cluster = 0;
  
  // std::cout << "Reached first for loop in cat function \n";
  
  for(arma::uword i = 0; i < num_clusters; i++){
    
    // calculate the log-weights for the context specific cluster and the across
    // context similarity
    curr_weight = log(cluster_weights_categorical(i));
    
    // Check if in the same cluster in both contexts
    common_cluster = 1 * (clust_labels_comp(row_index) == clust_labels(row_index));
                            
    similarity_upweight = log(1 + phi * common_cluster);
    
    for(arma::uword j = 0; j < num_cols_cat; j++){
      
      prob_vec(i) = prob_vec(i) + std::log(class_probs(j)(i, point(j)));
    }
    
    // As logs can sum rather than multiply the components
    prob_vec(i) = curr_weight + prob_vec(i) + similarity_upweight;
  }
  
  // to handle overflowing
  prob_vec = exp(prob_vec - max(prob_vec));
  
  // normalise
  prob_vec = prob_vec / sum(prob_vec);
  
  return prob_vec;
}


// sample calculate probabilties for cluster allocation in the gaussian data in 
// MDI (hence the presence of the context similarity parameter)
arma::vec mdi_gauss_clust_probs(arma::uword row_index,
                                arma::mat data,
                                arma::uword k,
                                arma::cube mu,
                                arma::cube variance,
                                double context_similarity,
                                arma::vec cluster_weights,
                                arma::uvec cluster_labels,
                                arma::uvec cluster_labels_comp,
                                bool outlier = false,
                                arma::vec global_mean = arma::zeros<arma::vec>(1),
                                arma::mat global_variance = arma::zeros<arma::mat>(1, 1),
                                double t_df = 4.0){
  
  double curr_weight;
  double exponent;
  double log_likelihood;
  
  // Upweight for similarity of contexts
  double similarity_upweight = 0.0;
  
  double log_det;
  arma::vec prob_vec(k);
  
  arma::uvec count_probs;
  
  arma::uword d = data.n_cols;
  
  arma::uword common_cluster = 0;
  
  arma::vec point = arma::trans(data.row(row_index));
  
  for(arma::uword i = 1; i < k + 1; i++){
    curr_weight = log(cluster_weights(i - 1));
    
    common_cluster = 1 * (cluster_labels(row_index) == cluster_labels_comp(row_index));
    
    similarity_upweight = log(1 + context_similarity * common_cluster);
    
    if(outlier && i == k){
      
      exponent = arma::as_scalar(arma::trans(point - global_mean) 
                                   * arma::inv(global_variance)
                                   * (point - global_mean));
                                   
                                   log_det = arma::log_det(global_variance).real();
                                   
                                   log_likelihood = (lgamma((t_df + d)/2.0)
                                                       - lgamma(t_df/2.0)
                                                       + d/2.0 * log(t_df * M_PI) 
                                                       + log_det - ((t_df + d)/2.0) * log(1 + (1/t_df) * exponent)
                                   );
                                   
    }
    else {
      exponent = arma::as_scalar(arma::trans(point - mu.slice(i - 1)) 
                                   * arma::inv(variance.slice(i - 1))
                                   * (point - mu.slice(i - 1)));
                                   
                                   log_det = arma::log_det(variance.slice(i - 1)).real();
                                   log_likelihood = -0.5 *(log_det + exponent + d * log(2 * M_PI));
    }
    prob_vec(i - 1) = curr_weight + log_likelihood + similarity_upweight;
  }
  
  prob_vec = exp(prob_vec - max(prob_vec));
  prob_vec = prob_vec / sum(prob_vec);
  
  return prob_vec;
}


// In a vector changes all values of ``label_1'' to ``label_2'' and vice versa
arma::uvec swap_labels(arma::uvec cluster_labels, 
                       arma::uword label_1, 
                       arma::uword label_2){
  
  arma::uvec label_1_ind = find(cluster_labels == label_1);
  arma::uvec label_2_ind = find(cluster_labels == label_2);
  
  cluster_labels.elem(label_1_ind).fill(label_2);
  cluster_labels.elem(label_2_ind).fill(label_1);
  return cluster_labels;
}

// Swap cluster weights
arma::vec swap_cluster_weights(arma::vec cluster_weights,
                               arma::uword label_1, 
                               arma::uword label_2){
  
  arma::vec decoy_weights = cluster_weights;
  
  cluster_weights(label_1) = decoy_weights(label_2);
  cluster_weights(label_2) = decoy_weights(label_1);
  
  return cluster_weights;
}

// log_likelihood of the model
double log_model_likelihood(double v,
                            double Z,
                            arma::uword n,
                            double context_similarity,
                            arma::uvec cluster_labels_1,
                            arma::uvec cluster_labels_2,
                            arma::vec cluster_weights_1,
                            arma::vec cluster_weights_2){
  double score = 0.0;
  
  for(arma::uword i = 0; i < n; i++){
    score += log(1 + context_similarity * (cluster_labels_1(i) == cluster_labels_2(i)))
    + log(cluster_weights_1(cluster_labels_1(i) - 1)) 
    + log(cluster_weights_2(cluster_labels_2(i) - 1));
  }
  score += (n - 1) * log(v) - v * Z - log_factorial(n - 1);
  
  return score;
}


//  Have to create a function for label swapping
// This will involve comparing likelihoods with and without swap and then 
// a rejection method

// will have to re-order cluster weights vector for dataset 2; easiest to record
// which classes flip and hence the index of the gammas which have to swap
// to generate a vector of the indices use: 
// std::vector<int> v(100) ; // vector with 100 ints.
// std::iota (std::begin(v), std::end(v), 0);
// 
// Will compare improvement in context similarity if cat_cluster_label(j) changed
// to associating with cont_cluster_label(i) and then record swapping of (j) and 
// (i) in the cluster labels in context 2 and the change in the gamma vector

arma::vec cluster_label_update(arma::uvec cluster_labels_1,
                               arma::uvec cluster_labels_2,
                               arma::vec cluster_weights_1,
                               arma::vec cluster_weights_2,
                               arma::uword num_clusters_1,
                               arma::uword num_clusters_2,
                               double phi,
                               arma::uword min_num_clusters,
                               double v,
                               arma::uword n,
                               double a0,
                               double b0,
                               double Z){
  
  arma::uword new_pos = 0;
  
  double new_phi = 0.0;
  
  arma::uvec new_labels(n);
  arma::vec new_weights(num_clusters_2);
  
  double log_accept = 0.0;
  double accept = 0.0;
  
  double old_score = 0.0;
  double new_score = 0.0;
  
  // Should this be bound as the min_num_clusters or min_num_clusters - 1?
  for(arma::uword i = 0; i < num_clusters_2; i++){
    
    // Generate a new position not equal to i
    // multiply a random number from [0,1] by the upper bound of i less 1
    // this less 1 is used so that if we equal i we can add 1 and treat all
    // cases similarly (otherwise we would be biasing the result towards i + 1)
    // Should it be "- 2" rather than "- 1"? Due to C++'s method of accessing
    // elements and the +1 to avoid selecting the same position as i, I think we
    // need to use "- 2".
    new_pos = floor(arma::randu<double>( ) * (num_clusters_2 - (1+1e-12)));
    
    if(new_pos >= i){
      new_pos++;
    }
    
    new_labels = swap_labels(cluster_labels_2, i + 1, new_pos + 1);
    
    new_weights = swap_cluster_weights(cluster_weights_2, i, new_pos);
    
    new_phi = sample_phi(cluster_labels_1,
                         new_labels,
                         cluster_weights_1,
                         new_weights,
                         v,
                         n,
                         min_num_clusters,
                         a0,
                         b0);
    
    old_score = log_model_likelihood(v,
                                     Z,
                                     n,
                                     phi,
                                     cluster_labels_1,
                                     cluster_labels_2,
                                     cluster_weights_1,
                                     cluster_weights_2);
    
    new_score = log_model_likelihood(v,
                                     Z,
                                     n,
                                     new_phi,
                                     cluster_labels_1,
                                     new_labels,
                                     cluster_weights_1,
                                     new_weights);
    
    
    log_accept = new_score - old_score;
    
    accept = 1.0;
    
    if(log_accept < 0){
      accept = exp(log_accept);
    }
    
    if(arma::randu<double>( ) < accept){
      
      // std::cout << "Labels\n";
      cluster_labels_2 = new_labels;
      
      // std::cout << "Weights\n";
      cluster_weights_2 = new_weights;
      
      // std::cout << "Phi\n\n";
      phi = new_phi;
      
      
    }
  }
  
  arma::vec output = arma::zeros<arma::vec>(n + num_clusters_2 + 1);
  output.subvec(0, n - 1) = arma::conv_to<arma::vec>::from(cluster_labels_2);
  
  output.subvec(n, n + num_clusters_2 - 1) = cluster_weights_2;
  output(n + num_clusters_2) = phi;
  return output;
}

// MDI clustering for a gaussian and cateogrical dataset
// [[Rcpp::export]]
Rcpp::List mdi_gauss_cat(arma::mat gaussian_data,
                         arma::umat categorical_data,
                         arma::vec mu_0,
                         double lambda_0,
                         arma::mat scale_0,
                         int df_0,
                         double a0,
                         double b0,
                         arma::vec cluster_weight_priors_gaussian,
                         arma::vec cluster_weight_priors_categorical,
                         arma::field<arma::vec> phi_prior,
                         arma::uvec cluster_labels_gaussian,
                         arma::uvec cluster_labels_categorical,
                         arma::uword num_clusters_gaussian,
                         arma::uword num_clusters_categorical,
                         std::vector<bool> fix_vec,
                         arma::uword num_iter,
                         arma::uword burn,
                         arma::uword thinning,
                         bool outlier = false,
                         double t_df = 4.0,
                         bool record_posteriors = false,
                         bool normalise = false
){
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  
  // std::cout << "In function \n";
  
  arma::uword n = gaussian_data.n_rows;
  arma::uword num_cols_cont = gaussian_data.n_cols;
  
  // arma::uword n_cat = categorical_data.n_rows;
  arma::uword num_cols_cat = categorical_data.n_cols;
  
  // if(n_cont != n_cat){
  //   std::cout << "Issue\n";
  // }
  
  // Add +1 to k to allow outlier class
  if(outlier){
    num_clusters_gaussian++;
  }
  
  // Frequently will compare clusters across contexts so this is a useful bound
  // to iterate to
  arma::uword min_num_clusters = std::min(num_clusters_gaussian,
                                          num_clusters_categorical);
  
  // std::cout << "Nomalising data\n";
  
  //  Normalise the continuous data
  if(normalise){
    gaussian_data = arma::normalise(gaussian_data);
  }
  
  // std::cout << "Data normalised\n";
  
  // Declare global variance and mean - used in outlier t-distribution 
  // (if included)
  arma::mat global_variance(num_cols_cont, num_cols_cont);
  global_variance = 0.5 * arma::cov(gaussian_data); // Olly's rec
  
  arma::vec global_mean(num_cols_cont);
  global_mean = arma::trans(arma::mean(gaussian_data, 0));
  
  
  double v = 0.0; // strategic latent variable
  
  // Cluster weights for each dataset
  arma::vec cluster_weights_gaussian(num_clusters_gaussian);
  arma::vec cluster_weights_categorical(num_clusters_categorical);
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::field<arma::cube> loc_mu_variance(2);
  
  // std::cout << "Declared to mean/variance thing \n";
  
  // Declare the field for the phi variable for the categorical data
  arma::uvec cat_count(num_cols_cat);
  cat_count = cat_counter(categorical_data);
  arma::field<arma::mat> class_probabilities(num_cols_cat);
  
  class_probabilities = declare_class_probs_field(cat_count,
                                                  num_cols_cat,
                                                  num_clusters_categorical);
  
  // double context_similarity = Rf_rgamma(a0, b0);
  
  double context_similarity = arma::randg(arma::distr_param(a0, 1/b0) );
  
  // cluster_weights_gaussian = dirichlet_posterior(cluster_weight_priors_gaussian,
  //                                                cluster_labels_gaussian,
  //                                                num_clusters_gaussian);
  // 
  // cluster_weights_categorical = dirichlet_posterior(cluster_weight_priors_categorical,
  //                                                   cluster_labels_categorical,
  //                                                   num_clusters_categorical);
  
  
  // std::cout << cluster_weights_gaussian << "\n\n";
  // std::cout << cluster_weights_categorical << "\n\n";
  // std::cout << context_similarity << "\n\n";
  
  double Z = 0.0;
  // calculate_normalising_constant(cluster_weights_gaussian,
  //                                           cluster_weights_categorical,
  //                                           context_similarity,
  //                                           num_clusters_gaussian,
  //                                           num_clusters_categorical);
  
  arma::vec curr_gaussian_prob_vec(num_clusters_gaussian);
  arma::vec curr_categorical_prob_vec(num_clusters_categorical);
  
  // std::cout << "Declared prob vectors \n";
  
  // the record for similarity in each clustering
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  arma::umat gaussian_record(n, eff_count);
  gaussian_record.zeros();
  
  arma::umat categorical_record(n, eff_count);
  categorical_record.zeros();
  
  arma::vec context_similarity_record(eff_count);
  
  arma::vec labels_weights_phi(n + num_clusters_categorical + 1);
  
  arma::vec entropy_cw(num_iter);
  
  // std::cout << "All declared \n";
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    
    
    // sample cluster weights for the two datasets
    cluster_weights_gaussian = dirichlet_posterior(cluster_weight_priors_gaussian,
                                                   cluster_labels_gaussian,
                                                   num_clusters_gaussian);
    
    cluster_weights_categorical = dirichlet_posterior(cluster_weight_priors_categorical,
                                                      cluster_labels_categorical,
                                                      num_clusters_categorical);
    
    
    
    // Entropy for graphing convergence
    entropy_cw(i) = entropy(cluster_weights_gaussian);
    
    // std::cout << "Sampled cluster weights \n";
    
    
    // std::cout << "\nSampling variance and mena \n";
    // std::cout << cluster_labels_gaussian << "\n\n";
    
    // Sample the posterior mean and variance for the gaussian data
    loc_mu_variance = mean_variance_sampling(gaussian_data,
                                             cluster_labels_gaussian,
                                             num_clusters_gaussian,
                                             df_0,
                                             num_cols_cont,
                                             scale_0,
                                             lambda_0,
                                             mu_0);
    
    // std::cout << "Variance sampled\n";
    
    // For the categorical data, sample the probabilities for each class
    class_probabilities = sample_class_probabilities(categorical_data,
                                                     class_probabilities,
                                                     phi_prior,
                                                     cluster_labels_categorical,
                                                     cat_count,
                                                     num_clusters_categorical,
                                                     num_cols_cat);
    
    // std::cout << "Sampled class probabilities for categorical data\n";
    
    // sample the context similarity parameter (as only two contexts this is a
    // single double - number not a drink)
    context_similarity = sample_phi(cluster_labels_gaussian,
                                    cluster_labels_categorical,
                                    cluster_weights_gaussian,
                                    cluster_weights_categorical,
                                    v,
                                    n,
                                    min_num_clusters,
                                    a0,
                                    b0);
    
    // std::cout << "Sampled phi\n";
    
    // Calculate the current normalising constant (consider being more clever 
    // about this) 
    Z = calculate_normalising_constant(cluster_weights_gaussian,
                                       cluster_weights_categorical,
                                       context_similarity,
                                       num_clusters_gaussian,
                                       num_clusters_categorical);
    
    // std::cout << "Z calculated \n";
    
    
    
    // sample the strategic latent variable, v
    // v = Rf_rgamma(n, Z);
    v = arma::randg( arma::distr_param(n, 1/Z) );
    
    // std::cout << Z << "\n";
    // std::cout << arma::randg( arma::distr_param(n, 1/Z) ) << "\n";
    // std::cout << arma::randg( arma::distr_param(n, Z) ) << "\n";
    
    // std::cout << "Sampled v \n";
    
    // sample 
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster
      // std::cout << "Gaussian cluser prob vec next\n";
      
      curr_gaussian_prob_vec = mdi_gauss_clust_probs(j,
                                                     gaussian_data,
                                                     num_clusters_gaussian,
                                                     loc_mu_variance(1),
                                                     loc_mu_variance(0),
                                                     context_similarity,
                                                     cluster_weights_gaussian,
                                                     cluster_labels_gaussian,
                                                     cluster_labels_categorical,
                                                     outlier,
                                                     global_mean,
                                                     global_variance,
                                                     t_df);
      
      // std::cout << "Categorical cluser prob vec next\n";
      
      curr_categorical_prob_vec = mdi_cat_clust_prob(j,
                                                     categorical_data,
                                                     class_probabilities,
                                                     num_clusters_categorical,
                                                     num_cols_cat,
                                                     context_similarity,
                                                     cluster_weights_categorical,
                                                     cluster_labels_gaussian,
                                                     cluster_labels_categorical);
      
      // std::cout << "Predict label per point\n";
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(! fix_vec[j]){
        // std::cout << "Gaussian\n";
        cluster_labels_gaussian(j) = cluster_predictor(curr_gaussian_prob_vec);
      }
      
      // std::cout << "Categorical\n";
      cluster_labels_categorical(j) = cluster_predictor(curr_categorical_prob_vec);
      
    }
    
    // std::cout << "All the context update stuff\n";
    
    // std::cout << cluster_labels_categorical.n_elem << "\n"; 
    
    // Update cluster labels in second dataset
    // Return the new labels, weights and similarity in a single vector
    labels_weights_phi = cluster_label_update(cluster_labels_gaussian,
                                              cluster_labels_categorical,
                                              cluster_weights_gaussian,
                                              cluster_weights_categorical,
                                              num_clusters_gaussian,
                                              num_clusters_categorical,
                                              context_similarity,
                                              min_num_clusters,
                                              v,
                                              n,
                                              a0,
                                              b0,
                                              Z);
    
    // labels_weights_phi = cluster_label_update(cluster_labels_gaussian,
    //                                           cluster_labels_categorical,
    //                                           cluster_weights_gaussian,
    //                                           cluster_weights_categorical,
    //                                           num_clusters_categorical,
    //                                           context_similarity,
    //                                           min_num_clusters,
    //                                           v,
    //                                           n,
    //                                           a0,
    //                                           b0);
    
    // Separate the output into the relevant components
    
    // std::cout << "Values calculated now sharing out\n";
    cluster_labels_categorical = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
    
    // std::cout << "cluster labels updated \n";
    
    // std::cout <<"\nCluster weights before:\n" << cluster_weights_categorical << "\n";
    
    cluster_weights_categorical = labels_weights_phi.subvec(n, n + num_clusters_categorical - 1);
    
    // std::cout <<"\nCluster weights after:\n" << cluster_weights_categorical << "\n\n";
    
    // std::cout <<"cluster weights updated \n";
    // std::cout <<"\nContext similarity before checking label swapping:\n" << context_similarity << "\n";
    
    context_similarity = arma::as_scalar(labels_weights_phi(n + num_clusters_categorical));
    // std::cout <<"phi updated \n";
    
    // std::cout <<"\nContext similarity after label swapping:\n" << context_similarity << "\n\n";
    
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thinning == 0) {
      gaussian_record.col((i - burn) / thinning) = cluster_labels_gaussian;
      categorical_record.col((i - burn) / thinning) = cluster_labels_categorical;
      context_similarity_record((i - burn) / thinning) = context_similarity;
    }
  }
  
  // construct similarity matrix
  arma::mat sim(n, n); 
  arma::mat cat_sim(n, n);
  sim = similarity_mat(gaussian_record);
  cat_sim = similarity_mat(categorical_record);
  
  // std::cout << "Context similarity: " << context_similarity << "\n";
  
  return List::create(Named("similarity") = sim,
                      Named("categorical_similarity") = cat_sim,
                      Named("gaussian_class_record") = gaussian_record,
                      Named("categorical_class_record") = categorical_record,
                      Named("context_similarity") = context_similarity_record,
                      Named("entropy") = entropy_cw);
}


///////////////////////////////////////////////////////////////////////////////

// MDI for different types

///////////////////////////////////////////////////////////////////////////////

// MDI clustering for two gaussian datasets
// [[Rcpp::export]]
Rcpp::List mdi_gauss_gauss(arma::mat data_1,
                           arma::mat data_2,
                           arma::vec mu_0_1,
                           double lambda_0_1,
                           arma::mat scale_0_1,
                           int df_0_1,
                           arma::vec mu_0_2,
                           double lambda_0_2,
                           arma::mat scale_0_2,
                           int df_0_2,
                           arma::vec clust_weight_priors_1,
                           arma::vec clust_weight_priors_2,
                           arma::uvec clust_labels_1,
                           arma::uvec clust_labels_2,
                           arma::uword n_clust_1,
                           arma::uword n_clust_2,
                           std::vector<bool> fix_vec_1,
                           std::vector<bool> fix_vec_2,
                           double a0,
                           double b0,
                           arma::uword num_iter,
                           arma::uword burn,
                           arma::uword thinning,
                           bool outlier_1 = false,
                           double t_df_1 = 4.0,
                           bool outlier_2 = false,
                           double t_df_2 = 4.0,
                           bool record_posteriors = false,
                           bool normalise_1 = false,
                           bool normalise_2 = false
){
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  
  // std::cout << "In function \n";
  
  arma::uword n = data_1.n_rows;
  arma::uword n_cols_1 = data_1.n_cols;
  
  // arma::uword n_cat = categorical_data.n_rows;
  arma::uword n_cols_2 = data_2.n_cols;
  
  // Add +1 to k to allow outlier class
  if(outlier_1){
    n_clust_1++;
  }
  
  if(outlier_2){
    n_clust_2++;
  }
  
  // Frequently will compare clusters across contexts so this is a useful bound
  // to iterate to
  arma::uword min_n_clust = std::min(n_clust_1, n_clust_2);
  
  // std::cout << "Nomalising data\n";
  
  //  Normalise the continuous data
  if(normalise_1){
    data_1 = arma::normalise(data_1);
  }
  
  if(normalise_2){
    data_2 = arma::normalise(data_2);
  }
  
  // Declare global variance and mean - used in outlier t-distribution 
  // (if included)
  arma::mat global_variance_1(n_cols_1, n_cols_1);
  global_variance_1 = 0.5 * arma::cov(data_1); // Olly's rec
  
  arma::vec global_mean_1(n_cols_1);
  global_mean_1 = arma::trans(arma::mean(data_1, 0));
  
  arma::mat global_variance_2(n_cols_2, n_cols_2);
  global_variance_2 = 0.5 * arma::cov(data_2); // Olly's rec
  
  arma::vec global_mean_2(n_cols_2);
  global_mean_2 = arma::trans(arma::mean(data_2, 0));
  
  
  double v = 0.0; // strategic latent variable
  
  // Cluster weights for each dataset
  arma::vec clust_weights_1(n_clust_1);
  arma::vec clust_weights_2(n_clust_2);
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::field<arma::cube> loc_mu_variance_1(2);
  arma::field<arma::cube> loc_mu_variance_2(2);
  
  // Context similarity - smaple prior
  double phi = arma::randg(arma::distr_param(a0, 1/b0) );
  
  // Declare normalising constant
  double Z = 0.0;
  
  // Used in each iteration
  arma::vec curr_prob_vec_1(n_clust_1);
  arma::vec curr_prob_vec_2(n_clust_2);
  
  // std::cout << "Declared prob vectors \n";
  
  // the record for similarity in each clustering
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  arma::umat record_1(n, eff_count);
  record_1.zeros();
  
  arma::umat record_2(n, eff_count);
  record_2.zeros();
  
  arma::vec phi_record(eff_count);
  
  arma::vec labels_weights_phi(n + n_clust_2 + 1);
  
  arma::vec entropy_cw(num_iter);
  
  // std::cout << "All declared \n";
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    // sample cluster weights for the two datasets
    clust_weights_1 = dirichlet_posterior(clust_weight_priors_1,
                                          clust_labels_1,
                                          n_clust_1);
    
    clust_weights_2 = dirichlet_posterior(clust_weight_priors_2,
                                          clust_labels_2,
                                          n_clust_2);
    
    // Entropy for graphing convergence
    entropy_cw(i) = entropy(clust_weights_1);
    
    // Sample the posterior mean and variance for the first dataset
    loc_mu_variance_1 = mean_variance_sampling(data_1,
                                               clust_labels_1,
                                               n_clust_1,
                                               df_0_1,
                                               n_cols_1,
                                               scale_0_1,
                                               lambda_0_1,
                                               mu_0_1);
    
    // std::cout << "Variance sampled\n";
    
    // Sample the posterior mean and variance for the second dataset
    loc_mu_variance_2 = mean_variance_sampling(data_2,
                                               clust_labels_2,
                                               n_clust_2,
                                               df_0_2,
                                               n_cols_2,
                                               scale_0_2,
                                               lambda_0_2,
                                               mu_0_2);
    
    // sample the context similarity parameter (as only two contexts this is a
    // single double - number not a drink)
    phi = sample_phi(clust_labels_1,
                     clust_labels_2,
                     clust_weights_1,
                     clust_weights_2,
                     v,
                     n,
                     min_n_clust,
                     a0,
                     b0);
    
    // std::cout << "Sampled phi\n";
    
    // Calculate the current normalising constant (consider being more clever 
    // about this) 
    Z = calculate_normalising_constant(clust_weights_1,
                                       clust_weights_2,
                                       phi,
                                       n_clust_1,
                                       n_clust_2);
    
    // std::cout << "Z calculated \n";
    
    bool data_2_unsupervised = true;
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(n, 1/Z) );
    
    // sample 
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster
      // std::cout << "Gaussian cluser prob vec next\n";
      
      
      curr_prob_vec_1 = mdi_gauss_clust_probs(j,
                                              data_1,
                                              n_clust_1,
                                              loc_mu_variance_1(1),
                                              loc_mu_variance_1(0),
                                              phi,
                                              clust_weights_1,
                                              clust_labels_1,
                                              clust_labels_2,
                                              outlier_1,
                                              global_mean_1,
                                              global_variance_1,
                                              t_df_1);
      
      // std::cout << "Categorical cluser prob vec next\n";
      
      curr_prob_vec_2 = mdi_gauss_clust_probs(j,
                                              data_2,
                                              n_clust_2,
                                              loc_mu_variance_2(1),
                                              loc_mu_variance_2(0),
                                              phi,
                                              clust_weights_2,
                                              clust_labels_2,
                                              clust_labels_1,
                                              outlier_2,
                                              global_mean_2,
                                              global_variance_2,
                                              t_df_2);
      
      // std::cout << "Predict label per point\n";
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(! fix_vec_1[j]){
        // std::cout << "Gaussian\n";
        clust_labels_1(j) = cluster_predictor(curr_prob_vec_1);
      }
      
      if(! fix_vec_2[j]){
        // std::cout << "Categorical\n";
        clust_labels_2(j) = cluster_predictor(curr_prob_vec_2);
      } else {
        data_2_unsupervised = false;
      }
    }
    
    // std::cout << "All the context update stuff\n";
    
    // std::cout << cluster_labels_categorical.n_elem << "\n"; 
    
    // Update cluster labels in second dataset
    // Return the new labels, weights and similarity in a single vector
    
    // Do not allow label flipping if any of context 2 have fixed labels
    
    if(data_2_unsupervised){
      labels_weights_phi = cluster_label_update(clust_labels_1,
                                                clust_labels_2,
                                                clust_weights_1,
                                                clust_weights_2,
                                                n_clust_1,
                                                n_clust_2,
                                                phi,
                                                min_n_clust,
                                                v,
                                                n,
                                                a0,
                                                b0,
                                                Z);
      
      // labels_weights_phi = cluster_label_update(cluster_labels_gaussian,
      //                                           cluster_labels_categorical,
      //                                           cluster_weights_gaussian,
      //                                           cluster_weights_categorical,
      //                                           num_clusters_categorical,
      //                                           context_similarity,
      //                                           min_num_clusters,
      //                                           v,
      //                                           n,
      //                                           a0,
      //                                           b0);
      
      // Separate the output into the relevant components
      
      // std::cout << "Values calculated now sharing out\n";
      clust_labels_2 = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
      
      // std::cout << "cluster labels updated \n";
      
      // std::cout <<"\nCluster weights before:\n" << cluster_weights_categorical << "\n";
      
      clust_weights_2 = labels_weights_phi.subvec(n, n + n_clust_2 - 1);
      
      // std::cout <<"\nCluster weights after:\n" << cluster_weights_categorical << "\n\n";
      
      // std::cout <<"cluster weights updated \n";
      // std::cout <<"\nContext similarity before checking label swapping:\n" << context_similarity << "\n";
      
      phi = arma::as_scalar(labels_weights_phi(n + n_clust_2));
      // std::cout <<"phi updated \n";
      
      // std::cout <<"\nContext similarity after label swapping:\n" << context_similarity << "\n\n";
      
    }
    
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thinning == 0) {
      record_1.col((i - burn) / thinning) = clust_labels_1;
      record_2.col((i - burn) / thinning) = clust_labels_2;
      phi_record((i - burn) / thinning) = phi;
    }
  }
  
  // construct similarity matrix
  arma::mat sim_1(n, n); 
  arma::mat sim_2(n, n);
  sim_1 = similarity_mat(record_1);
  sim_2 = similarity_mat(record_2);
  
  // std::cout << "Context similarity: " << context_similarity << "\n";
  
  return List::create(Named("similarity_1") = sim_1,
                      Named("similarity_2") = sim_2,
                      Named("class_record_1") = record_1,
                      Named("class_record_1") = record_2,
                      Named("context_similarity") = phi_record,
                      Named("entropy") = entropy_cw);
  
}

// MDI clustering for two categorical datasets
// [[Rcpp::export]]
Rcpp::List mdi_cat_cat(arma::umat data_1,
                       arma::umat data_2,
                       arma::field<arma::vec> class_dist_prior_1,
                       arma::field<arma::vec> class_dist_prior_2,
                       arma::vec clust_weight_priors_1,
                       arma::vec clust_weight_priors_2,
                       arma::uvec clust_labels_1,
                       arma::uvec clust_labels_2,
                       arma::uword n_clust_1,
                       arma::uword n_clust_2,
                       std::vector<bool> fix_vec_1,
                       std::vector<bool> fix_vec_2,
                       double a0,
                       double b0,
                       arma::uword num_iter,
                       arma::uword burn,
                       arma::uword thinning
){
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  
  // std::cout << "In function \n";
  
  arma::uword n = data_1.n_rows;
  arma::uword n_cols_1 = data_1.n_cols;
  
  // arma::uword n_cat = categorical_data.n_rows;
  arma::uword n_cols_2 = data_2.n_cols;
  
  // Frequently will compare clusters across contexts so this is a useful bound
  // to iterate to
  arma::uword min_n_clust = std::min(n_clust_1, n_clust_2);
  
  double v = 0.0; // strategic latent variable
  
  // Cluster weights for each dataset
  arma::vec clust_weights_1(n_clust_1);
  arma::vec clust_weights_2(n_clust_2);
  
  // Declare the field for the phi variable for the categorical data
  arma::uvec cat_count_1(n_cols_1);
  cat_count_1 = cat_counter(data_1);
  arma::field<arma::mat> class_prob_1(n_cols_1);
  
  class_prob_1 = declare_class_probs_field(cat_count_1,
                                           n_cols_1,
                                           n_clust_1);
  
  
  arma::uvec cat_count_2(n_cols_2);
  cat_count_2 = cat_counter(data_2);
  arma::field<arma::mat> class_prob_2(n_cols_2);
  
  class_prob_2 = declare_class_probs_field(cat_count_2,
                                           n_cols_2,
                                           n_clust_2);
  
  // Context similarity - smaple prior
  double phi = arma::randg(arma::distr_param(a0, 1/b0) );
  
  // Declare normalising constant
  double Z = 0.0;
  
  // Used in each iteration
  arma::vec curr_prob_vec_1(n_clust_1);
  arma::vec curr_prob_vec_2(n_clust_2);
  
  // std::cout << "Declared prob vectors \n";
  
  // the record for similarity in each clustering
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  arma::umat record_1(n, eff_count);
  record_1.zeros();
  
  arma::umat record_2(n, eff_count);
  record_2.zeros();
  
  arma::vec phi_record(eff_count);
  
  arma::vec labels_weights_phi(n + n_clust_2 + 1);
  
  arma::vec entropy_cw(num_iter);
  
  bool data_2_unsupervised = true;
  
  // std::cout << "All declared \n";
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    // sample cluster weights for the two datasets
    clust_weights_1 = dirichlet_posterior(clust_weight_priors_1,
                                          clust_labels_1,
                                          n_clust_1);
    
    clust_weights_2 = dirichlet_posterior(clust_weight_priors_2,
                                          clust_labels_2,
                                          n_clust_2);
    
    // Entropy for graphing convergence
    entropy_cw(i) = entropy(clust_weights_1);
    
    // For the categorical data, sample the probabilities for each class
    class_prob_1 = sample_class_probabilities(data_1,
                                              class_prob_1,
                                              class_dist_prior_1,
                                              clust_labels_1,
                                              cat_count_1,
                                              n_clust_1,
                                              n_cols_1);
    
    class_prob_2 = sample_class_probabilities(data_2,
                                              class_prob_2,
                                              class_dist_prior_2,
                                              clust_labels_2,
                                              cat_count_2,
                                              n_clust_2,
                                              n_cols_2);
    
    // sample the context similarity parameter (as only two contexts this is a
    // single double - number not a drink)
    phi = sample_phi(clust_labels_1,
                     clust_labels_2,
                     clust_weights_1,
                     clust_weights_2,
                     v,
                     n,
                     min_n_clust,
                     a0,
                     b0);
    
    // std::cout << "Sampled phi\n";
    
    // Calculate the current normalising constant (consider being more clever 
    // about this) 
    Z = calculate_normalising_constant(clust_weights_1,
                                       clust_weights_2,
                                       phi,
                                       n_clust_1,
                                       n_clust_2);
    
    // std::cout << "Z calculated \n";
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(n, 1/Z) );
    
    // sample 
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster
      // std::cout << "Gaussian cluser prob vec next\n";
      
      
      curr_prob_vec_1 = mdi_cat_clust_prob(j, 
                                           data_1,
                                           class_prob_1,
                                           n_clust_1,
                                           n_cols_1,
                                           phi,
                                           clust_weights_1,
                                           clust_labels_1,
                                           clust_labels_2);
      
      curr_prob_vec_2 = mdi_cat_clust_prob(j, 
                                           data_2,
                                           class_prob_2,
                                           n_clust_2,
                                           n_cols_2,
                                           phi,
                                           clust_weights_2,
                                           clust_labels_2,
                                           clust_labels_1);
      
      // std::cout << "Predict label per point\n";
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(! fix_vec_1[j]){
        // std::cout << "Gaussian\n";
        clust_labels_1(j) = cluster_predictor(curr_prob_vec_1);
      }
      
      if(! fix_vec_2[j]){
        // std::cout << "Categorical\n";
        clust_labels_2(j) = cluster_predictor(curr_prob_vec_2);
      } else {
        data_2_unsupervised = false;
      }
    }
    
    // std::cout << "All the context update stuff\n";
    
    // std::cout << cluster_labels_categorical.n_elem << "\n"; 
    
    if(data_2_unsupervised){
      
      // Update cluster labels in second dataset
      // Return the new labels, weights and similarity in a single vector
      labels_weights_phi = cluster_label_update(clust_labels_1,
                                                clust_labels_2,
                                                clust_weights_1,
                                                clust_weights_2,
                                                n_clust_1,
                                                n_clust_2,
                                                phi,
                                                min_n_clust,
                                                v,
                                                n,
                                                a0,
                                                b0,
                                                Z);
      
      // labels_weights_phi = cluster_label_update(cluster_labels_gaussian,
      //                                           cluster_labels_categorical,
      //                                           cluster_weights_gaussian,
      //                                           cluster_weights_categorical,
      //                                           num_clusters_categorical,
      //                                           context_similarity,
      //                                           min_num_clusters,
      //                                           v,
      //                                           n,
      //                                           a0,
      //                                           b0);
      
      // Separate the output into the relevant components
      
      // std::cout << "Values calculated now sharing out\n";
      clust_labels_2 = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
      
      // std::cout << "cluster labels updated \n";
      
      // std::cout <<"\nCluster weights before:\n" << cluster_weights_categorical << "\n";
      
      clust_weights_2 = labels_weights_phi.subvec(n, n + n_clust_2 - 1);
      
      // std::cout <<"\nCluster weights after:\n" << cluster_weights_categorical << "\n\n";
      
      // std::cout <<"cluster weights updated \n";
      // std::cout <<"\nContext similarity before checking label swapping:\n" << context_similarity << "\n";
      
      phi = arma::as_scalar(labels_weights_phi(n + n_clust_2));
      // std::cout <<"phi updated \n";
      
      // std::cout <<"\nContext similarity after label swapping:\n" << context_similarity << "\n\n";
    }
      
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thinning == 0) {
      record_1.col((i - burn) / thinning) = clust_labels_1;
      record_2.col((i - burn) / thinning) = clust_labels_2;
      phi_record((i - burn) / thinning) = phi;
    }
  }
  
  // construct similarity matrix
  arma::mat sim_1(n, n); 
  arma::mat sim_2(n, n);
  sim_1 = similarity_mat(record_1);
  sim_2 = similarity_mat(record_2);
  
  // std::cout << "Context similarity: " << context_similarity << "\n";
  
  return List::create(Named("similarity_1") = sim_1,
                      Named("similarity_2") = sim_2,
                      Named("class_record_1") = record_1,
                      Named("class_record_1") = record_2,
                      Named("context_similarity") = phi_record,
                      Named("entropy") = entropy_cw);
  
}
