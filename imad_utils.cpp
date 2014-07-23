#include "imad.h"
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <limits>
#include <assert.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions.hpp> //for cdf()

//typedef Map<Matrix<double,Dynamic,Dynamic,RowMajor> > MapRMMatrixXd;
//defined in imad.h

using std::cout;
using std::endl;

namespace imad_utils{

  /* This function takes a vector of eigenvalues and a matrix
   * of eigenvectors (in columns). It sorts the eigenvalues in
   * descending order and then reorders the eigenvectors appropriately */
  //NOTE: This function is probably super inefficient --K.Song, 2014-7-14

  void reorder_eigens(VectorXd& lambda, MatrixXd& A, MatrixXd& B){

    int N = lambda.rows();

    /* We create a vector of (eigenvalue, eigenvector) pairs, then sort it
     * based on the eigenvalue (using templating of "<" which is defined below).
     * We then sort this vector based on the eigenvalues */
    std::vector<Eigentup> pairs = std::vector<Eigentup>();  //"Eigentuple"
    for(int i = 0; i < N; i++){     //eigenvalue  //eigenvec   //eigenvec
      Eigentup thistup = Eigentup(  lambda(i),     A.col(i),    B.col(i));
      pairs.push_back(thistup);
    }
    //Sort in descending order by eigenvalue
    std::sort(pairs.begin(), pairs.end());

    /* Now populate the original matrix and eigenvalue vector. Eigenvectors
     * are stored in columns, so this may be backwards compared to what
     * you are used to seeing in C codes. */
    for(int i = 0; i < N; i++){
      bool neg_eig = (pairs[i].eigenval < 0);
      //Make sure eigenvalues are positive
      lambda(i) = neg_eig ? (-1) * pairs[i].eigenval : pairs[i].eigenval;
      for(int j = 0; j < N; j++){
        //If eigenvalue was made positive, flip the eigenvector
        A(j,i) = neg_eig ? (-1) * pairs[i].eigenvec_a(j) : pairs[i].eigenvec_a(j);
        B(j,i) = neg_eig ? (-1) * pairs[i].eigenvec_b(j) : pairs[i].eigenvec_b(j);
      }
    }
    return;
  }

  bool Eigentup::operator< (Eigentup const &other) const{
    return (eigenval > other.eigenval);
  }

  Eigentup::Eigentup(double val, VectorXd vec_a, VectorXd vec_b){
    eigenval = val;
    eigenvec_a = vec_a;
    eigenvec_b = vec_b;
  }



  /*************************************************************************/

  /* Subtracts a column vector from each column of the matrix. One hell of a
   * lot more efficient than the Python solution (repeat the vector as many
   * times as necessary to get a matrix, then do matrix subtraction). Does
   * two subtractions at once to make the code a little cleaner */

  void rowwise_subtract(MatrixRXd& A, VectorXd& toSubtract1,
                        MatrixRXd& B, VectorXd& toSubtract2){
    assert(A.cols() == toSubtract1.rows());
    assert(B.cols() == toSubtract2.rows());
    for(int i = 0; i < A.rows(); i++){
      A.row(i) = A.row(i) - toSubtract1.transpose();
      B.row(i) = B.row(i) - toSubtract2.transpose();
    }
  }

  /* Takes a matrix and a vector and does an element-wise division.
   * NOTE: toDivide is still a column vector, we just treat it as a row*/

  void rowwise_divide(MatrixXd& A, VectorXd& toDivide){
     assert(A.cols() == toDivide.rows());
     for(int i = 0; i < A.rows(); i++){
       A.row(i) = (A.row(i).array() / toDivide.transpose().array()).matrix();
     }
  }

  /* Takes a matrix and a vector and does an element-wise multiplication. */

  void colwise_multiply(MatrixXd& A, VectorXd& toMult){
    assert(A.rows() == toMult.rows());
    for(int i = 0; i < A.cols(); i++){
      A.col(i) = (A.col(i).array() * toMult.array()).matrix();
    }
  }

  /*************************************************************************/

   /* Modifies the buffersize in-place, return value is number of chunks n
   needed in order to competely cover the input row. Can theoretically
   handle an image the size of the solar system...hopefully that won't be
   needed anytime soon. */

  int find_chunksize(int& buffersize, int ncol, int nBands){
    const int max_memory_in_bytes = 1000000000; //1GB in bytes
    const int max_bufsize = max_memory_in_bytes / ncol / nBands / sizeof(double);

    if(ncol < max_bufsize){ //Easy!
      buffersize = ncol;
      return 1;
    }
    else{
      for(int i = 1; i < std::numeric_limits<int>::max() - 2; i++){
        if(ncol / i < max_bufsize){ //Search through possible ways to split
          buffersize = ncol / i + 1; //rounding up
          return i;
        }
      }
    }
    throw std::invalid_argument("find_chunksize: Image is too large to chunk...what the hell do you have?");
  }

/*************************************************************************/

/* What follows are some of the most poorly-commented lines from the original
 * Python implementation, and my understanding of them is still fuzzy.
 * I have included modified versions of the original comments, and the
 * lines are each numbered, with comment-footnotes at the end of this file.
 * They have been moved out of the main function because they don't actually
 * contribute to the math and are incredibly cryptic. */

 //You are not expected to understand these lines. I certainly don't.

  void math_cleanup(MatrixXd& A, MatrixXd& B, MatrixXd& s11, MatrixXd& s12){
/* 1 */ MatrixXd D = s11.diagonal().array()      //Bookkeeping
                        .sqrt().cwiseInverse()   //Inverse sqrt
                        .matrix().asDiagonal();  //More bookkeeping
/* 2 */ VectorXd s = (D*s11*A).rowwise().sum();
        VectorXd tmp1 = (s.array() / (s.array().abs())).matrix();
/* 3 */ imad_utils::colwise_multiply(A, tmp1);
/* 4 */ VectorXd cov = (A.transpose() * s12 * B).diagonal();
/* 5 */ MatrixXd tmp2 = (cov.array() / cov.array().abs()).matrix().asDiagonal();
        B = B * tmp2;
  }
}


  /*************************************************************************/

/**** Notes on tile matrix ****/
/* tile is a strange name for a matrix (borrowed from the python code).
 * It holds a single row of both images. The top half of the matrix holds
 * the row from file1, with all of its bands in separate rows (e.g. for any
 * given row, row 1 of tile will hold image1, band1, row 2 will hold image1,
 * band2, etc.). The same is true of the second image in the bottom
 * half of tile. Note that this tile is transposed compared to the one in
 * Python: Pyversion has rows in cols and files split left/right. */

/**** Notes on numbered lines in imad() ****/
/* While going through these, it is helpful to remember that left multiplying by
 * a diagonal matrix (D*A) multiplies each row by the diagonal element, i.e.
 * row1(A) = row1(A) * D(1,1), row2(A) = row2(A) * D(2,2), etc. and right-multiplying
 * by a diagonal matrix multiplies the columns instead of the rows. Off we go!*/

/* 1 */ /* Completely unintelligible. First, we take the diagonal of s11. We then
 * cast it to an array (to deal with some issues with eigen), then take an
 * element-wise inverse square root. We then revert it back to a diagonal matrix.
 * The end effect is identical to masking all but the diagonal of s11, then
 * taking an element-wise inverse square root. */

/* 2 */ /* s is actually a row-vector. Not sure what this actually does */

/* 3 */ /* If S(i,i) was negative, we set it to -1. Otherwise, we set it to +1.
 * We then multiply A by this matrix to ensure the correlation sum is positive
 * (Yes, this part feels like magic to me too). */

/* 4,5 */ /* A similar process for B */
