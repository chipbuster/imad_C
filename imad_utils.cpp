#include "imad.h"
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <vector>

namespace imad_util{
  void reorder_eigens(VectorXf * lambda, MatrixXf* A, MatrixXf* B){
    /* This function takes a vector of eigenvalues and a matrix
     * of eigenvectors (in columns). It sorts the eigenvalues in
     * descending order and then reorders the eigenvectors appropriately */
    //NOTE: This function is probably super inefficient --K.Song, 2014-7-14

    int N = lambda->rows();

    /* We create a vector of (eigenvalue, eigenvector) pairs, then sort it
     * based on the eigenvalue (using templating of "<" which is defined below).
     * We then sort this vector based on the eigenvalues */
    std::vector<Eigentup> pairs = std::vector<Eigentup>();  //"Eigentuple"
    for(int i = 0; i < N; i++){       //eigenvalue    //eigenvec   //eigenvec
      Eigentup thistup = Eigentup( (*lambda)(i),   A->col(i),    B->col(i));
      pairs.push_back(thistup);
    }
    //Sort in order and reverse. Should probably fix this later.
    std::sort(pairs.begin(), pairs.end());
    std::reverse(pairs.begin(), pairs.end());

    /* Now populate the original matrix and eigenvalue vector. Eigenvectors
     * are stored in columns, so the idiom is a little backwards. */
    for(int i = 0; i < N; i++){
      (*lambda)(i) = pairs[i].eigenval;
      for(int j = 0; j < N; j++){
        (*A)(j,i) = pairs[i].eigenvec_a(j);
        (*B)(j,i) = pairs[i].eigenvec_b(j);
      }
    }
    return;
  }

  bool Eigentup::operator < (Eigentup const &other) const{
    return (eigenval < other.eigenval);
  }

  Eigentup::Eigentup(double val, VectorXf vec_a, VectorXf vec_b){
    eigenval = val;
    eigenvec_a = vec_a;
    eigenvec_b = vec_b;
  }
}
