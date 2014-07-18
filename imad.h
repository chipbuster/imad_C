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
/* Namespaces are named after their .cpp files, e.g. all functions in namespace
 * GdalFileIO are implemented in GdalFileIO.cpp. */

typedef Map<Matrix<float,Dynamic,Dynamic,RowMajor> > MapRMMatrixXf;
typedef Matrix<float,Dynamic,Dynamic,RowMajor> MatrixRXf;

namespace GdalFileIO{

  struct CoordTransform{
    /* Takes the Geotransform of an image and uses them to transform coordinate
     * Only one CoordTransform object per image projection, please. */
    //Look up affine transforms for more details on implementation.
    /* in the transform constants, X,Y refer to geospatial coordinates, while
     * P,L refer to the raster (x,y) coordinates. Raster coordinates start
     * in the top left and increase going right and down, so a pixel at
     * (P,L) = (2,5) is in the second column and fifth row of that image.
     * A derivative (e.g. dA/dB) is denoted as dA_dB in these names.
     * See docs for GDAL's GetGeoTransform() for more information. */
     Matrix3d Img2Geo;
     Matrix3d Geo2Img;
     Vector3d  input;
     Vector3d  output;

     //Needed for vectorizable ops
     EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CoordTransform(double* transform_coeff);
    double ImgtoGeo_X(double imgX, double imgY);
    double ImgtoGeo_Y(double imgX, double imgY);
    double GeotoImg_X(double geoX, double geoY);
    double GeotoImg_Y(double geoX, double geoY);
    void GeotoImg(double& imgX, double& imgY);
    void ImgtoGeo(double& geoX, double& geoY);
  };

  //Found in openGdalFiles.cpp
  GDALDataset* openFile(std::string filename);
  bool dimensionsmatch(GDALDataset* dataset1, GDALDataset* dataset2);
  std::vector<int>* selectBands();

}

namespace imad_utils{

  struct Eigentup{
    bool operator < (Eigentup const &other) const;
    Eigentup(double val, VectorXf vec_a, VectorXf vec_b);

      double eigenval;
      VectorXf eigenvec_a;
      VectorXf eigenvec_b;
  };

  void reorder_eigens(VectorXf& lambda, MatrixXf& A, MatrixXf& B);
  void colwise_subtract(MapRMMatrixXf& A, VectorXf& toSubtract);
  void rowwise_divide(MatrixXf& A, VectorXf& toDivide);
  void colwise_multiply(MatrixXf& A, VectorXf& toMult);
  RowVectorXf& getWeights(RowVectorXf& inputs, RowVectorXf& outputs, std::vector<int>& bands);
}

class ImageStats{
public:
  int n2Bands;
  double sum_weights;
  VectorXf means;
  MatrixXf covar;


  ImageStats(size_t bands);
  VectorXf get_means();
  MatrixXf get_covar();

  void zero();
  void update(float* input, float* weights, int nrow, int ncol);
};

#endif
