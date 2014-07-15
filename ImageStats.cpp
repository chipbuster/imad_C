#include "imad.h"

//typedef Map<Matrix<float,Dynamic,Dynamic,RowMajor> > MapRMMatrixXf;
//defined in imad.h

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
  covar = (covar.array() / (sum_weights - 1.0)).matrix();
  MatrixXf diagonal = covar.diagonal().asDiagonal();
  return covar + covar.transpose() - diagonal;
}

void ImageStats::zero(){
  means.setZero(n2Bands);
  covar.setZero(n2Bands,n2Bands);
  sum_weights = 0.0;
}

//TODO: Separate updates for mean and covar?

void ImageStats::update(MapRMMatrixXf& input,
                        MatrixXf& weights, size_t nrow, size_t ncol){
  /* inputs has rows in rows, each row is a diff. band. This is a transpose
   * of the Python algorithm, where we have rows in columns and each column
   * is a difference band. */
  double weight, ratio;
  double* diff = new double[nrow]; //Difference between element and mean
  bool no_weights = (weights.cols() == 0 && weights.rows() == 0);

  for(int pix = 0; pix < input.cols(); pix++ ){
    /* We traverse down columns to get the averages. Each column in input is the
     * same pixel, but in different bands/images. */

    weight = no_weights ? 1 : weights(pix); //If no weights, weight default to 1
    sum_weights += weight;
    ratio = weight / sum_weights;

    //Calculate mean of band via provisional means algorithm
    for(int band = 0; band < input.cols(); band++){
      diff[band] = input(band,pix) - means(band);
      means(band) += diff[band] * ratio;
    }
    /* Calculate covariance similarly, fill in upper tri covar matrix
     * B2 is the entry for band 1, b2 is the entry for band 2. */
    for(int b1 = 0; b1 < input.cols(); b1++){
      for(int b2 = b1; b2 < input.rows(); b2++){
        covar(b1,b2) = diff[b1] * diff[b2] * (1-ratio) * weight;
      }
    }
  }

  delete[] diff;
}

/* Python stores an image row in a column, with different bands in
 * different columns. The C++ code stores images in rows, with different
 * bands in distinct rows (that is, in C++, you read right to find more
 * pixels from the same band, and read down to find the pixel value from
 * a different band in the same location, whereas in Python, you read
 * down to find more pixels in the same band, and right for different
 * bands). */
