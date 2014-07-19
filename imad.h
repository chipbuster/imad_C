#ifndef HELPERS_H
#define HELPERS_H

#include "gdal_priv.h"
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>


using namespace Eigen;

/* This file contains the headers for the iMad project. */
/* All GDAL functions derive from gdal_priv.h           */
/* Namespaces are named after their .cpp files, e.g. all functions in namespace
 * GdalFileIO are implemented in GdalFileIO.cpp. */

typedef Map<Matrix<double,Dynamic,Dynamic,RowMajor> > MapRMMatrixXd;
typedef Matrix<double,Dynamic,Dynamic,RowMajor> MatrixRXd;

namespace GdalFileIO{

  GDALDataset* openFile(std::string filename);
  bool dimensionsmatch(GDALDataset* dataset1, GDALDataset* dataset2);
  std::vector<int>* selectBands();

}

namespace imad_utils{

  struct Eigentup{
    bool operator < (Eigentup const &other) const;
    Eigentup(double val, VectorXd vec_a, VectorXd vec_b);

    double eigenval;
    VectorXd eigenvec_a;
    VectorXd eigenvec_b;
  };

  void reorder_eigens(VectorXd& lambda, MatrixXd& A, MatrixXd& B);
  void rowwise_subtract(MatrixRXd& A, VectorXd& toSubtract1,
                        MatrixRXd& B, VectorXd& toSubtract2);
  void rowwise_divide(MatrixXd& A, VectorXd& toDivide);
  void colwise_multiply(MatrixXd& A, VectorXd& toMult);
  VectorXd& getWeights(VectorXd& inputs, VectorXd& outputs, size_t num_bands);
}

namespace geo_utils{

  struct ImageInfo{
    int ncol;
    int nrow;
    int nBands;
    char* projection;
    double* geotransform;
    ImageInfo(GDALDataset* target);
    ~ImageInfo();

    bool operator == (const ImageInfo& other) const;
  };

  class CoordTransform{
     Matrix3d Img2Geo;
     Matrix3d Geo2Img;
     Vector3d  input;
     Vector3d  output;

   public:
     //Needed for vectorizable ops
     EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CoordTransform(double* transform_coeff);
    CoordTransform(GDALDataset* input_file);
    ~CoordTransform();

    //return one coordinate
    double ImgtoGeo_X(double imgX, double imgY);
    double ImgtoGeo_Y(double imgX, double imgY);
    double GeotoImg_X(double geoX, double geoY);
    double GeotoImg_Y(double geoX, double geoY);

    //Modify coordinates in-place
    void GeotoImg(double& imgX, double& imgY);
    void ImgtoGeo(double& geoX, double& geoY);
  };
}

class ImageStats{
  int n2Bands;
  double sum_weights;
  VectorXd means;
  MatrixXd covar;

 public:
  ImageStats(size_t bands);
  ~ImageStats();
  VectorXd get_means();
  MatrixXd get_covar();

  void reset();
  void update(double* input, double* weights, int nrow, int ncol);
};

#endif
