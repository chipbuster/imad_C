#ifndef IMAGE_STATS_H
#define IMAGE_STATS_H

#include <Eigen/Dense>

using Eigen::VectorXf;
using Eigen::MatrixXf;

class ImageStats{
  int n2Bands;
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

#endif
