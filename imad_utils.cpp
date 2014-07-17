#include "imad.h"
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions.hpp> //for cdf()

//typedef Map<Matrix<float,Dynamic,Dynamic,RowMajor> > MapRMMatrixXf;
//defined in imad.h

namespace imad_utils{

  void reorder_eigens(VectorXf& lambda, MatrixXf& A, MatrixXf& B){
    /* This function takes a vector of eigenvalues and a matrix
     * of eigenvectors (in columns). It sorts the eigenvalues in
     * descending order and then reorders the eigenvectors appropriately */
    //NOTE: This function is probably super inefficient --K.Song, 2014-7-14

    //TODO: Rewrite function in accord with structures using eigen structs

    int N = lambda.rows();

    /* We create a vector of (eigenvalue, eigenvector) pairs, then sort it
     * based on the eigenvalue (using templating of "<" which is defined below).
     * We then sort this vector based on the eigenvalues */
    std::vector<Eigentup> pairs = std::vector<Eigentup>();  //"Eigentuple"
    for(int i = 0; i < N; i++){     //eigenvalue  //eigenvec   //eigenvec
      Eigentup thistup = Eigentup(  lambda(i),     A.col(i),    B.col(i));
      pairs.push_back(thistup);
    }
    //Sort in order and reverse. Should probably fix this later.
    std::sort(pairs.begin(), pairs.end());
    std::reverse(pairs.begin(), pairs.end());

    /* Now populate the original matrix and eigenvalue vector. Eigenvectors
     * are stored in columns, so this may be backwards compared to what
     * you are used to seeing in C codes. */
    for(int i = 0; i < N; i++){
      lambda(i) = pairs[i].eigenval;
      for(int j = 0; j < N; j++){
        A(j,i) = pairs[i].eigenvec_a(j);
        B(j,i) = pairs[i].eigenvec_b(j);
      }
    }
    return;
  }

  /*************************************************************************/

  void colwise_subtract(MapRMMatrixXf& A, VectorXf& toSubtract){
    /* Subtracts a column vector from each column of the matrix. One hell of a
     * lot more efficient than the Python solution (repeat the vector as many
     * times as necessary to get a matrix, then do matrix subtraction) */
    if(A.rows() != toSubtract.rows()){
      throw std::invalid_argument("Matrix and vector must have same number of rows!");
    }

    for(int i = 0; i < A.cols(); i++){
      A.col(i) = A.col(i) - toSubtract;
    }
  }

  void rowwise_divide(MatrixXf& A, VectorXf& toDivide){
    /* Takes a matrix and a vector and does an element-wise division
     * NOTE: toDivide is still a column vector, we just treat it as a row*/

     if(A.cols() != toDivide.rows()){
       throw std::invalid_argument("Matrix ncol must be same as vector length!");
     }
     for(int i = 0; i < A.rows(); i++){
       A.row(i) = (A.row(i).array() / toDivide.transpose().array()).matrix();
     }
  }

  /*************************************************************************/

  RowVectorXf& getWeights(RowVectorXf& input, RowVectorXf& output, std::vector<int>& bands){
    /* Gets CDF of Chisquared and places into weights array. CDF function is
     * backwards compared to Python (in python, param to chi2.cdf is second arg)
     * so cdf(rng(a), b) = stats.chi2.cdf(b,a) */
    if(input.cols() != (int)bands.size()){
      throw std::invalid_argument("Must have same number of inputs as bands!");
    }
    for(size_t i = 0; i < bands.size(); i++){
      boost::math::chi_squared_distribution<double> rng(bands[i]);
      double tmpout = cdf(rng, input(i));
      output(i) = tmpout;
    }
    return output;
  }

  /*************************************************************************/

  bool Eigentup::operator< (Eigentup const &other) const{
    return (eigenval < other.eigenval);
  }

  Eigentup::Eigentup(double val, VectorXf vec_a, VectorXf vec_b){
    eigenval = val;
    eigenvec_a = vec_a;
    eigenvec_b = vec_b;
  }
}
