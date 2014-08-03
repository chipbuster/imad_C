#include "imad.h"
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstring>
#include <assert.h>

namespace geo_utils{
  /** begin class CoordTransform, defined in header **/

  CoordTransform::CoordTransform(double* transform_coeff){
    //Construct affine transform matrix
    Img2Geo(0,2) = transform_coeff[0];
    Img2Geo(1,2) = transform_coeff[3];
    Img2Geo(0,0) = transform_coeff[1];
    Img2Geo(0,1) = transform_coeff[2];
    Img2Geo(1,0) = transform_coeff[4];
    Img2Geo(1,1) = transform_coeff[5];
    Img2Geo(2,0) = 0;
    Img2Geo(2,1) = 0;
    Img2Geo(2,2) = 1;
    Geo2Img = Img2Geo.inverse();

    //Last element of affine transform vectors is always 1
    input(2) = 1;
  }
  //Alternate constructor: get geotransform for user from a dataset handle
  CoordTransform::CoordTransform(GDALDataset* input_file){
    double* transform_coeff = new double[6];
    input_file->GetGeoTransform(transform_coeff);
    Img2Geo(0,2) = transform_coeff[0];
    Img2Geo(1,2) = transform_coeff[3];
    Img2Geo(0,0) = transform_coeff[1];
    Img2Geo(0,1) = transform_coeff[2];
    Img2Geo(1,0) = transform_coeff[4];
    Img2Geo(1,1) = transform_coeff[5];
    Img2Geo(2,0) = 0;
    Img2Geo(2,1) = 0;
    Img2Geo(2,2) = 1;
    Geo2Img = Img2Geo.inverse();

    //Last element of affine transform vectors is always 1
    input(2) = 1;
    delete[] transform_coeff;
  }
  CoordTransform::~CoordTransform(){
    Img2Geo.resize(0,0);
    Geo2Img.resize(0,0);
  }

  double CoordTransform::ImgtoGeo_X(double imgP, double imgL){
    input(0) = imgP;
    input(1) = imgL;
    output = Img2Geo * input;
    return output(0);
  }
  double CoordTransform::ImgtoGeo_Y(double imgP, double imgL){
    input(0) = imgP;
    input(1) = imgL;
    output = Img2Geo * input;
    return output(1);
  }
  double CoordTransform::GeotoImg_X(double geoX, double geoY){
    input(0) = geoX;
    input(1) = geoY;
    output = Geo2Img * input;
    return round(output(0));
  }
  double CoordTransform::GeotoImg_Y(double geoX, double geoY){
    input(0) = geoX;
    input(1) = geoY;
    output = Geo2Img * input;
    return round(output(1));
  }
  void CoordTransform::GeotoImg(double& X, double& Y){
    double tmpX = GeotoImg_X(X,Y);
    double tmpY = GeotoImg_Y(X,Y);
    X = tmpX; Y = tmpY;
  }
  void CoordTransform::ImgtoGeo(double& X, double& Y){
    double tmpX = ImgtoGeo_X(X,Y);
    double tmpY = ImgtoGeo_Y(X,Y);
    X = tmpX; Y = tmpY;
  }

  void CoordTransform::GeotoImg(Coord& inp){
    double tmpX = ImgtoGeo_X(inp.x,inp.y);
    double tmpY = ImgtoGeo_Y(inp.x,inp.y);
    inp.x = tmpX; inp.y = tmpY;
  }

  void CoordTransform::ImgtoGeo(Coord& inp){
    double tmpX = ImgtoGeo_X(inp.x,inp.y);
    double tmpY = ImgtoGeo_Y(inp.x,inp.y);
    inp.x = tmpX; inp.y = tmpY;
  }

  /** begin class ImageInfo **/

  ImageInfo::ImageInfo(GDALDataset* input){
     ncol = input->GetRasterXSize();
     nrow = input->GetRasterYSize();
   nBands = input->GetRasterCount();

     //GetProjectionRef's pointer can be gc'ed. Make our own copy for records.
     const char* tmp = input->GetProjectionRef();
     size_t len = strlen(tmp);
     projection = new char[len + 1]; //Don't forget the null terminator!
     strcpy(projection,tmp);

     //Projection copied. Let's get the geotransform
     geotransform = new double[6];
     input->GetGeoTransform(geotransform);
  }
  ImageInfo::~ImageInfo(){
    delete[] projection;
    delete[] geotransform;
  }
  bool ImageInfo::operator == (const ImageInfo& other) const{
    if( other.ncol != ncol ||
        other.nrow != nrow ||
        other.nBands != nBands
      ) return false;
    for(size_t i = 0; i < strlen(projection); i++){
      if(other.projection[i] != projection[i]) return false;
    }
    for(size_t i = 0; i < 6; i++){
      if(other.geotransform[i] != geotransform[i]) return false;
    }
    return true;
  }

  bool ImageInfo::compatible(const ImageInfo& other) const{
    /* Check that dimensions of the two files match in X, Y, and num. bands
     * Compatibility does not require that the origins match, only that
     * the projections and pixel sizes do. */
    if(other.ncol != ncol){
      std::cout<< "X-dimension mismatch in images. Exiting..." << std::endl;
      return false;
    }
    if(other.nrow != nrow){
      std::cout<< "Y-dimension mismatch in images. Exiting..." << std::endl;
      return false;
    }
    if(other.nBands != nBands){
      std::cout<< "Raster count mismatch in images. Exiting..." << std::endl;
      return false;
    }
    //If projections are not equal
    if(strcmp(other.projection,projection)){
      std::cout<< "Projections of files do not match. Exiting..." << std::endl;
      return false;
    }
    //Check pixel sizes
    if(abs(geotransform[1] - other.geotransform[4]) > 0.00001){
      std::cout<< "Pixel sizes do not match. Exiting..." << std::endl;
      return false;
    }
    if(abs(geotransform[2] - other.geotransform[5]) > 0.00001){
      std::cout<< "Pixel sizes do not match. Exiting..." << std::endl;
      return false;
    }
    return true;
  }


}
