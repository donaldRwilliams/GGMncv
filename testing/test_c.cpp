#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat remove_row(arma::mat x, int which){
  x.shed_row(which);
  return(x);
}

arma::mat remove_col(arma::mat x, int index){
  x.shed_col(index);
  return(x);
}

arma::mat Sigma_i_not_i(arma::mat x, int index) {
  arma::mat sub_x = x.row(index);
  sub_x.shed_col(index);
  return(sub_x);
}

// [[Rcpp::export]]
Rcpp::List hft(arma::mat Sigma, arma::mat adj, double tol) {

  arma::mat S = Sigma;
  arma::mat W = S;

  arma::uvec upper_indices = trimatu_ind( size(S) );
  arma::mat W_previous = S;
  double p = S.n_cols;
  arma::mat iter(1,1, arma::fill::zeros);
  double max_diff = 100;
  arma::mat w_12(1,p-1);

  while(max_diff > tol){

    for(int i = 0; i < p; ++i){

      arma::mat beta(1,p-1, arma::fill::zeros);
      arma::uvec pad_index =  find(Sigma_i_not_i(adj, i) == 1);

      if(pad_index.n_elem == 0 ){
        w_12 = beta;
      } else {

        arma::mat W_11 = remove_row(remove_col(W , i), i);
        arma::mat s_12 = Sigma_i_not_i(S,i);

        arma::mat W_11_star = W_11.submat(pad_index, pad_index);
        arma::mat s_12_star = s_12(pad_index);

        beta(pad_index) = inv(W_11_star) * s_12_star;
        arma::mat w_12 = W_11 * beta.t();
        arma::mat temp = W.col(i).row(i);
        w_12.insert_rows(i, temp);

        for(int k = 0; k < p; ++k){
          W.row(i).col(k) = w_12(k);
          W.row(k).col(i) = w_12(k);
        }

        max_diff = max(W.elem(upper_indices) -  W_previous.elem(upper_indices));

        W_previous = W;
      }


    }

    iter(0,0) = iter(0,0) + 1;
  }

  arma::mat Theta = inv(W) % adj;
  arma::mat W_return = inv(Theta);

 List ret;
 ret["Theta"] = Theta;
 ret["Sigma"] = W_return;
 ret["iter"]  =  iter;
 return ret;
}

/*** R
# diag(adj) <- 1
# # system.time({
# # replicate(100,
# #  test <- hft(cor(Y), adj,  0.000001)
# # qpgraph::qpHTF(cor(Y), adj, tol = 0.000001)
# # # test2 <- GGMncv:::htf(cor(Y), adj = adj, tol = 0.000001)
# # )})
# test <- hft(cor(Y), adj,  0.001)
# test2 <- GGMncv:::htf(cor(Y), adj = adj, tol = 0.001)
#
# round(test$Sigma, 3) -round(test2$Sigma,3)
# t <- qpgraph::qpHTF(cor(Y), adj, tol = 0.0001)
# round(solve(t),3) - round(solve(test2$Sigma),3)
#
# #
# #
# # solve(solve(test$W) * adj)
# #
# # test$W[,7]
# # round(solve(test$W), 3)[1,]
# #
# # round(test$W[2,], 3)
# # test2 <- GGMncv:::htf(cor(Y), adj = adj)
# # test2$Theta
# # round(test2$Sigma[2,], 3)
# # qpgraph::qpHTF(cor(Y), adj)
# # solve(qpgraph::qpHTF(cor(Y), adj))[1,]
*/
