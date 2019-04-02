// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <set>
#include <utility>
#include <algorithm>
#include <vector>

using namespace Rcpp;
using namespace std;



Eigen::VectorXd ColDispersion_sp(Eigen::SparseMatrix<double> mat){
  int nrows = mat.rows();

  Eigen::VectorXd dispersion(mat.cols());

  for (int i=0; i<mat.outerSize(); ++i){

    double mean = 0;
    double var = 0;
    int counter = 0;

    for (Eigen::SparseMatrix<double>::InnerIterator iter(mat,i); iter; ++iter){
      mean += iter.value();
    }

    mean = mean / nrows;
    for (Eigen::SparseMatrix<double>::InnerIterator iter(mat,i); iter; ++iter){
      var += pow(iter.value() - mean, 2);
      counter += 1;
    }
    var = (var + (nrows - counter) * pow(mean, 2)) / (nrows - 1);
    dispersion[i] = var/mean;

  }
  return(dispersion);
}


Eigen::VectorXd ColDispersion_dense(Eigen::MatrixXd mat){
  int nrows = mat.rows();

  Eigen::VectorXd dispersion(mat.cols());

  for (int i=0; i<mat.outerSize(); ++i){

    double mean = 0;
    double var = 0;
    int counter = 0;

    for (Eigen::MatrixXd::InnerIterator iter(mat,i); iter; ++iter){
      mean += iter.value();
    }


    mean = mean / nrows;
    for (Eigen::MatrixXd::InnerIterator iter(mat,i); iter; ++iter){
      var += pow(iter.value() - mean, 2);
      counter += 1;
    }
    var = (var + (nrows - counter) * pow(mean, 2)) / (nrows - 1);
    dispersion[i] = var/mean;

  }
  return(dispersion);
}


//[[Rcpp::export]]
Eigen::VectorXd ColDispersion(SEXP mat) {
  if (Rf_isS4(mat)) {
    if(Rf_inherits(mat, "dgCMatrix")) {
      return ColDispersion_sp(as<Eigen::SparseMatrix<double>>(mat)) ;
    } ;
    stop("unsupported matrix class.") ;
  } else {
    return ColDispersion_dense(as<Eigen::MatrixXd>(mat)) ;
  }
}

inline int check_close(double v1, double v2, double close_thresh) {
  if(v1 > v2){
    swap(v1, v2);
  }
  if(abs(v1) < 0.001){
    return 0;
  }

  if((v2 - v1)/v1 < close_thresh){
    return 1;
  } else {
    return 0;
  }
}

int too_close_total(NumericVector & gval1, NumericVector & gval2, double close_thresh) {
  int n = gval1.length();

  int result = 0;
  for ( int i = 0; i < n; ++i ) {
    result += check_close(gval1(i), gval2(i), close_thresh);
  }
  return result;
}

// [[Rcpp::export]]
RcppExport SEXP reduce_genes_cpp(List & l) {
  NumericMatrix mat = as<NumericMatrix>(l["matrix"]);
  double close_thresh = as<double>(l["close_threshold"]);
  double cells_thresh = as<double>(l["cells_threshold"]);
  double nreduce = as<int>(l["nreduce"]);
  int ncells = mat.nrow();
  int ngenes = mat.ncol();

  set<pair<int, int> > check_pairs;

  for(int i = 0; i < ncells; i ++) {
    vector<pair<double, int> > temp_vec;
    for(int j = 0; j < ngenes; j ++) {
      temp_vec.push_back(make_pair(mat(i, j), j));
    }
    sort(temp_vec.begin(), temp_vec.end());
    for(int j = 1; j < ngenes; j ++) {
      check_pairs.insert(make_pair(min(temp_vec[j - 1].second, temp_vec[j].second), max(temp_vec[j - 1].second, temp_vec[j].second)));
    }
  }

  int total_removed = 0;
  vector<int> removed_genes(ngenes, 0);
  for(pair<int, int> pa : check_pairs) {
    int gene1 = pa.first, gene2 = pa.second;
    if(!removed_genes[gene1] and !removed_genes[gene2]) {
      NumericVector col1 = mat(_, gene1), col2 = mat(_, gene2);
      if(too_close_total(col1, col2, close_thresh) > cells_thresh * ncells) {
        removed_genes[gene1] = 1;
        total_removed ++;
      }
    }
    if(total_removed >= nreduce) break;
  }
  NumericVector ret(removed_genes.begin(), removed_genes.end());
  return ret;
}
