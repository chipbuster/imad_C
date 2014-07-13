#include <Eigen/Dense>

//Bad practice, yada yada. This does keep the build system rather clean.
//We'll move executables to another file if time allows.

using Eigen::VectorXf;
using Eigen::MatrixXf;

class ImageStats{
  int n2Bands;
  /* Because the number of bands passed to ImageStats is 2x the number of
   * bands defined in the main iMad algorithm, it has a different name
   * to avoid confusion. */
  double sum_weights;
  VectorXf means;
  MatrixXf covar;

public:
  ImageStats(int bands);
  VectorXf get_means();
  MatrixXf get_covar();

  void zero();
  void update(float* inputs, float* weights, size_t nrow, size_t ncol);
};

ImageStats::ImageStats(int input){
    n2Bands = input;
    means(n2Bands);
    covar(n2Bands,n2Bands);
    sum_weights = 0.0;
}

VectorXf ImageStats::get_means(){
  return means;
}

MatrixXf ImageStats::get_covar(){
  covar = covar / (sum_weights - 1.0);
  MatrixXf diagonal = covar.diagonal().asDiagonal();
  return covar + covar.transpose() - diagonal;
}

void ImageStats::zero(){
  means.setZero(n2Bands);
  covar.setZero(n2Bands,n2Bands);
  sum_weights = 0.0;
}

//TODO: Separate updates for mean and covar?

void ImageStats::update(float* inputs, float* weights, size_t nrow, size_t ncol){
  double weight, ratio;
  double* diff = new double[ncol]; //Difference between element and mean

  for(size_t thisrow = 0; thisrow < nrow; thisrow++){
    /*Eww ternary statement. If we have no weights vector, weights are all 1.
      Otherwise, weights are the values that are given to us */
    weight = weights==NULL ? 1 : weights[thisrow];
    this->sum_weights += weight;
    ratio = weight / this->sum_weights;

    //Calculate mean
    for(size_t index = 0; index < ncol; index++){
      diff[index] = inputs[thisrow*ncol + index] - means(index);
      means(index) += diff[index] * ratio;
    }

    //Fill in upper triangular matrix of covariances
    for(size_t j = 0; j < ncol; j++){
      for(size_t k = j; k < ncol; k++){
        covar(thisrow,k) += diff[j]*diff[k]*(1-ratio)*weight;
      }
    }
  }

  delete[] diff;
}
