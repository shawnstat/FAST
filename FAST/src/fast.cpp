// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

//References
// John St. John. (2014). GrNMF:A Network Constrained Version of Non-negative Matrix Factorization.R package version 0.2
// Cai, D., He, X., Wu, X., & Han, J. (2008). Non-negative Matrix Factorization on Manifold. 2008 Eighth IEEE International Conference on Data Mining (ICDM), 63-72. doi:10.1109/ICDM.2008.57
// Xu, W., Liu, X., & Gong, Y. (2003). Document clustering based on non-negative matrix factorization. the 26th annual international ACM SIGIR conference (pp. 267-273). New York, New York, USA: ACM. doi:10.1145/860435.860485
// Hofree, M., Shen, J. P., Carter, H., Gross, A., & Ideker, T. (2013). Network-based stratification of tumor mutations. Nature Methods, 10(11), 1108-1115. doi:10.1038/nmeth.2651

//' @export
// [[Rcpp::export]]
List fast(NumericMatrix Xr, NumericMatrix Gr, int r=5, double lambda_1r=1.0, double lambda_2r=1.0, int n_iter=1000, double converge=1e-6) {

    int p = Xr.nrow(), n = Xr.ncol();

    if(Gr.nrow() != n || Gr.ncol() != n){
      throw std::invalid_argument("The edge-weight matrix G should be n by n (where n is the number of spots)!");
    }
    // r must be no more than the number of features
    if(r > p)
      throw std::invalid_argument("r must be less than or equal to p, the number of genes");

    double lastFit=std::numeric_limits<double>::max();

    int convergenceCheckFrequency=10;
 //   int iterations_to_modify_lambda=100;

    double lambda_1 = lambda_1r;
    double lambda_2 = lambda_2r;

    if(lambda_1<0 | lambda_2<0){
      throw std::invalid_argument("lambdas must be >= 0!");
    }

    arma::mat X(Xr.begin(), p, n, false);       // reuses memory and avoids extra copy
    arma::mat G(Gr.begin(), n, n, false);
    arma::mat W(p, r, arma::fill::randu); // initialize W,H to random values in 0-1
    arma::mat H(n, r, arma::fill::randu);
    arma::mat Eo(r,n);
    arma::mat Fo(n,n);
    Eo.ones();
    Fo.ones();


    double minX = X.min(), maxX = X.max();

    if(minX < 0){
      throw std::invalid_argument("The input matrices should be non-negative!");
    }


    // set W/H to be random values in the range
    // of X, the NMF package in R recommends this
//    W*=(maxX-minX);
//    W+=minX+1e-4;
//    H*=(maxX-minX);
//    H+=minX+1e-4;

    // the D matrix is the diagonal matrix of row (or column since G is symmetric)
    // sums of G.
    arma::mat D = arma::diagmat(arma::sum(G));
    // L is the graph laplacian matrix
    arma::mat L = arma::mat(D - G);

    //loop over iterations
    int it=0;
    for(it = 0; it < n_iter; it++){

      if(it % convergenceCheckFrequency == 0){
        double fit =  arma::as_scalar(arma::norm(X-W*H.t(),"fro")) + arma::as_scalar(lambda_1 * arma::trace(H.t() * L * H)) + arma::as_scalar(lambda_2*arma::accu(W));
        //Rcout << "The value of fit : " << fit << "\n";
        //allow the user to turn off by setting converge to a negative
        if(std::abs(lastFit-fit) <= converge && converge >= 0)
          break;
        lastFit=fit;
      }
      // The following does element wize multiplicative updates
      // on the element-wise division

      W %= (X * H) / (W * H.t() * H);
      H %= (X.t()* W + lambda_1 * G * H + lambda_2*Fo*Eo.t()) / (H * W.t() * W + lambda_1 * D * H + lambda_2*H*Eo*Eo.t());
    }

    return Rcpp::List::create(
        Rcpp::Named("W") = W,
        Rcpp::Named("H") = H,
        Rcpp::Named("Max.iter") = it,
        Rcpp::Named("ObjectiveMain") =  arma::as_scalar(arma::norm(X-W*H.t(),"fro")),
        Rcpp::Named("ObjectiveGraph") =   arma::as_scalar(lambda_1 * arma::trace(H.t() * L * H)),
//        Rcpp::Named("ObjectiveSparse") =  arma::as_scalar(lambda_2 * arma::accu(W)),
        Rcpp::Named("ObjectiveFitNMF") =  arma::as_scalar(arma::norm(X-W*H.t(),"fro")) +  arma::as_scalar(lambda_1 * arma::trace(H.t() * L * H)) + arma::as_scalar(lambda_2 * arma::norm(H*Eo-Fo,"fro"))
    ) ;
}
