#ifndef HELPERS_H
#define HELPERS_H

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "ImageStats.h" //Includes <Eigen/Dense>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>

/* This file contains the headers for the iMad project. */
/* All GDAL functions derive from gdal_priv.h           */

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

namespace imad_util{

  struct Eigentup{
    double eigenval;
    VectorXf eigenvec_a;
    VectorXf eigenvec_b;
    bool operator < (Eigentup const &other) const;
    Eigentup(double val, VectorXf vec_a, VectorXf vec_b);
  };

  void reorder_eigens(VectorXf * lambda, MatrixXf* A, MatrixXf* B);

}

#endif
