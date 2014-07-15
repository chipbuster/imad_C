#ifndef HELPERS_H
#define HELPERS_H

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>


using namespace Eigen;

/* This file contains the headers for the iMad project. */
/* All GDAL functions derive from gdal_priv.h           */

typedef Map<Matrix<float,Dynamic,Dynamic,RowMajor> > MapRMMatrixXf;

namespace GdalFileIO{

/*
  struct ImageDims{
    int xoffset;
    int yoffset;
    int nrow;
    int ncol;

    bool operator==(const ImageDims& other) const{
      if(xoffset == other.xoffset &&
         yoffset == other.yoffset &&
         nrow    == other.nrow    &&
         ncol    == other.ncol ){
           return true;
         }
      else return false;
    }
  };
  */

  //Found in openGdalFiles.cpp
  GDALDataset* openFile(std::string filename);
  bool dimensionsmatch(GDALDataset* dataset1, GDALDataset* dataset2);
  std::vector<int>* selectBands();

}

namespace imad_utils{

  struct Eigentup{
    double eigenval;
    VectorXf eigenvec_a;
    VectorXf eigenvec_b;
    bool operator < (Eigentup const &other) const;
    Eigentup(double val, VectorXf vec_a, VectorXf vec_b);
  };

  void reorder_eigens(VectorXf& lambda, MatrixXf& A, MatrixXf& B);
  void colwise_subtract(MapRMMatrixXf& A, VectorXf& toSubtract);
}



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
  void update(MapRMMatrixXf& input, MatrixXf& weights, size_t nrow, size_t ncol);
};

#endif
