// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// similarity_mat
arma::mat similarity_mat(arma::Mat<int> cluster_record);
RcppExport SEXP _BayesicGibbs_similarity_mat(SEXP cluster_recordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type cluster_record(cluster_recordSEXP);
    rcpp_result_gen = Rcpp::wrap(similarity_mat(cluster_record));
    return rcpp_result_gen;
END_RCPP
}
// cluster_weight_posterior
arma::vec cluster_weight_posterior(arma::vec concentration_0, arma::Col<int> cluster_labels, int num_clusters);
RcppExport SEXP _BayesicGibbs_cluster_weight_posterior(SEXP concentration_0SEXP, SEXP cluster_labelsSEXP, SEXP num_clustersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type concentration_0(concentration_0SEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type cluster_labels(cluster_labelsSEXP);
    Rcpp::traits::input_parameter< int >::type num_clusters(num_clustersSEXP);
    rcpp_result_gen = Rcpp::wrap(cluster_weight_posterior(concentration_0, cluster_labels, num_clusters));
    return rcpp_result_gen;
END_RCPP
}
// cat_counter
arma::uvec cat_counter(arma::mat data);
RcppExport SEXP _BayesicGibbs_cat_counter(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(cat_counter(data));
    return rcpp_result_gen;
END_RCPP
}
// sampling
arma::mat sampling(arma::mat data, arma::field<arma::vec> phi_prior, arma::Col<int> cluster_labels, arma::vec fix_vec, arma::vec cluster_weight_priors, arma::uword num_clusters, arma::uword num_iter, arma::uword burn, arma::uword thinning);
RcppExport SEXP _BayesicGibbs_sampling(SEXP dataSEXP, SEXP phi_priorSEXP, SEXP cluster_labelsSEXP, SEXP fix_vecSEXP, SEXP cluster_weight_priorsSEXP, SEXP num_clustersSEXP, SEXP num_iterSEXP, SEXP burnSEXP, SEXP thinningSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type phi_prior(phi_priorSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type cluster_labels(cluster_labelsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fix_vec(fix_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cluster_weight_priors(cluster_weight_priorsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type num_clusters(num_clustersSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thinning(thinningSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling(data, phi_prior, cluster_labels, fix_vec, cluster_weight_priors, num_clusters, num_iter, burn, thinning));
    return rcpp_result_gen;
END_RCPP
}
// entropy
double entropy(arma::vec class_weights);
RcppExport SEXP _BayesicGibbs_entropy(SEXP class_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type class_weights(class_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(entropy(class_weights));
    return rcpp_result_gen;
END_RCPP
}
// point_comparison
Rcpp::List point_comparison(arma::uword num_iter, arma::vec concentration_0, arma::mat scale_0, arma::uvec class_labels, std::vector<bool> fix_vec, arma::vec mu_0, double lambda_0, arma::mat data, int df_0, arma::uword k, arma::uword burn, arma::uword thinning, bool outlier, double t_df, bool record_posteriors);
RcppExport SEXP _BayesicGibbs_point_comparison(SEXP num_iterSEXP, SEXP concentration_0SEXP, SEXP scale_0SEXP, SEXP class_labelsSEXP, SEXP fix_vecSEXP, SEXP mu_0SEXP, SEXP lambda_0SEXP, SEXP dataSEXP, SEXP df_0SEXP, SEXP kSEXP, SEXP burnSEXP, SEXP thinningSEXP, SEXP outlierSEXP, SEXP t_dfSEXP, SEXP record_posteriorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type concentration_0(concentration_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale_0(scale_0SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type class_labels(class_labelsSEXP);
    Rcpp::traits::input_parameter< std::vector<bool> >::type fix_vec(fix_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_0(mu_0SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_0(lambda_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type df_0(df_0SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thinning(thinningSEXP);
    Rcpp::traits::input_parameter< bool >::type outlier(outlierSEXP);
    Rcpp::traits::input_parameter< double >::type t_df(t_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type record_posteriors(record_posteriorsSEXP);
    rcpp_result_gen = Rcpp::wrap(point_comparison(num_iter, concentration_0, scale_0, class_labels, fix_vec, mu_0, lambda_0, data, df_0, k, burn, thinning, outlier, t_df, record_posteriors));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _BayesicGibbs_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _BayesicGibbs_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _BayesicGibbs_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _BayesicGibbs_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesicGibbs_similarity_mat", (DL_FUNC) &_BayesicGibbs_similarity_mat, 1},
    {"_BayesicGibbs_cluster_weight_posterior", (DL_FUNC) &_BayesicGibbs_cluster_weight_posterior, 3},
    {"_BayesicGibbs_cat_counter", (DL_FUNC) &_BayesicGibbs_cat_counter, 1},
    {"_BayesicGibbs_sampling", (DL_FUNC) &_BayesicGibbs_sampling, 9},
    {"_BayesicGibbs_entropy", (DL_FUNC) &_BayesicGibbs_entropy, 1},
    {"_BayesicGibbs_point_comparison", (DL_FUNC) &_BayesicGibbs_point_comparison, 15},
    {"_BayesicGibbs_rcpparma_hello_world", (DL_FUNC) &_BayesicGibbs_rcpparma_hello_world, 0},
    {"_BayesicGibbs_rcpparma_outerproduct", (DL_FUNC) &_BayesicGibbs_rcpparma_outerproduct, 1},
    {"_BayesicGibbs_rcpparma_innerproduct", (DL_FUNC) &_BayesicGibbs_rcpparma_innerproduct, 1},
    {"_BayesicGibbs_rcpparma_bothproducts", (DL_FUNC) &_BayesicGibbs_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesicGibbs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
