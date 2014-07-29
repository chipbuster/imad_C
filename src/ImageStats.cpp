#include "imad.h"

//typedef Map<Matrix<float,Dynamic,Dynamic,RowMajor> > MapRMMatrixXd;
//defined in imad.h

using std::cout; //debugging purposes
using std::endl;

ImageStats::ImageStats(size_t input){
    n2Bands = input;
    means = VectorXd(n2Bands);
    covar = MatrixXd(n2Bands,n2Bands);
    means.setZero(n2Bands);
    covar.setZero(n2Bands,n2Bands);
    sum_weights = 0.00001; //Not zero due to issues with div by zero
}

ImageStats::~ImageStats(){
  means.resize(0);
  covar.resize(0,0);
}

VectorXd ImageStats::get_means(){
  return means;
}

MatrixXd ImageStats::get_covar(){
  if(covar == MatrixXd::Zero(n2Bands,n2Bands)){
    cout << "Error: iMad was unable to find pixel values that had data in the "
         << "region specified. If you are using automatic overlap detection, "
         << "this probably means that the detected overlap area has no data. "
         << "Please try selecting the overlapping region manually. Also, note "
         << "that a single no-data value in any layer causes that pixel to be "
         << "treated as no data in ALL the layers." << endl;
    throw std::runtime_error("No data in the image region specified");
  }
  covar = (covar / (sum_weights - 1.0));
  MatrixXd diagonal = covar.diagonal().asDiagonal();
  return covar + covar.transpose() - diagonal;
}

void ImageStats::reset(){
  means.setZero(n2Bands);
  covar.setZero(n2Bands,n2Bands);
  sum_weights = 0.00001;
}

//TODO: Separate updates for mean and covar?

void ImageStats::update(double* input,
                        double* weights, int nrow, int ncol){

  double weight, ratio;
  double* diff = new double[nrow]; //Difference between element and mean

  for(int row = 0; row < nrow; row++ ){
    /* Check for NODATA values. If any of the values are nodata, we increment
     * the row counter to skip that pixel. Note that this means that the outer
     * loop is NOT guaranteed to run nrow times---rows may be skipped here.*/
    bool has_nodata = false;
    for(int i = 0 ; i < ncol; i++){
      if(input[i + ncol * row] == -9999){
        has_nodata = true;
        break;
      }
    }
    if(has_nodata) continue;

    //If no weights, weight defaults to 1
    weight = weights == NULL ? 1 : weights[row];
    sum_weights += weight;
    ratio = weight / sum_weights;

    //Calculate mean of band via provisional means algorithm
    for(int band = 0; band < ncol; band++){
      diff[band] = input[band + ncol * row] - means(band);
      means(band) += diff[band] * ratio;
    }
    /* Calculate covariance similarly, fill in upper tri covar matrix
     * B2 is the entry for band 1, b2 is the entry for band 2. */
    for(int b1 = 0; b1 < ncol; b1++){
      for(int b2 = b1; b2 < ncol; b2++){
        covar(b2,b1) += diff[b1] * diff[b2] * (1-ratio) * weight;
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
