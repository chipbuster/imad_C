#include "gdal_priv.h"
#include "imad.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <math.h>

using namespace std;

namespace GdalFileIO{

    GDALDataset* openFile(string filename){
      /* Returns a GDALDataset handle to the file given by file1.
       * If the file cannot be opened, an exception is thrown */

      if(filename.empty()){
        cout << "Please enter a file to be opened: ";
        getline(cin, filename);
      }

      GDALDataset* dataset1 = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);

      if(dataset1 == NULL){
        throw std::invalid_argument("Error while opening file" + filename);
      }
      return dataset1;
    }

    /*************************************************************************/

    bool dimensionsmatch(GDALDataset* dataset1, GDALDataset* dataset2){
      /* Check that dimensions of the two files match in X, Y, and num. bands */
      if(dataset1->GetRasterXSize() != dataset2->GetRasterXSize()){
        cout << "X-dimension mismatch in images. Exiting..." << endl;
        return false;
      }
      if(dataset1->GetRasterYSize() != dataset2->GetRasterYSize()){
        cout << "Y-dimension mismatch in images. Exiting..." << endl;
        return false;
      }
      if(dataset1->GetRasterCount() != dataset2->GetRasterCount()){
        cout << "Raster count mismatch in images. Exiting..." << endl;
        return false;
      }
      //All dimensions match!
      return true;
    }

    /*************************************************************************/

    vector<int>* selectBands(){
      vector<int>* bands = new vector<int>();
      cout << "Please enter the band numbers you would like to use, followed by a return"
           << endl << "When you are finished, please enter 0. Yes this interface sucks." << endl;

      string input;
      cin >> input;
      while(input != "0"){
        bands->push_back(atof(input.c_str()));
        cin >> input;
      }
      if (bands->size() < 1){
        throw std::out_of_range("Need to specify at least one band.");
      }
      return bands;
    }

    /*************************************************************************/

    CoordTransform::CoordTransform(double* transform_coeff){
      //Notes can be found in the header file. Construct affine transform matrix
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

    double CoordTransform::ImgtoGeo_X(double imgP, double imgL){
      input(0) = imgP;
      input(1) = imgL;
      input(2) = 1;
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

    /*************************************************************************/
}
