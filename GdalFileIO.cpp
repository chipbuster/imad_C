#include "gdal_priv.h"
#include "imad.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <math.h>
#include <string.h> //Gdal libraries use const char*, need strcmp()

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
      return geo_utils::ImageInfo(dataset1).compatible(geo_utils::ImageInfo(dataset2));
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
    void getOutputFileInfo(string& output_file, string& format){
      if(format.empty()){
        std::cout<< "Valid formats: geotransformiff, PCIDSK, HFA, ENVI" << std::endl;
        std::cout<< "Please enter output format: ";
        getline(std::cin, format);
      }
      if(output_file.empty()){
        std::cout<< "Please enter output file name";
        getline(std::cin, output_file);
      }
    }

    /*************************************************************************/
    void writeOutputToFile(GDALDataset* outfile, double* tile,
                           MatrixXd& A, MatrixXd& B, //Eigenvector matrices
                           int xoffset, int yoffset, int ncol, int nrow,
                           GDALRasterBand** bands_1,
                           GDALRasterBand** bands_2,
                           GDALDataset* reference_file,
                           VectorXd& sigMADs){
      int nBands = A.cols();
      double* geotransform   = new double[6];
      reference_file->GetGeoTransform(geotransform);
        //Move the origins of the picture to the overlap zone. Not implemented yet.
          //  geotransform[0] = geotransform[0] + xof*geotransform[1]
          //  geotransform[3] = geotransform[3] + y10*geotransform[5]
      const char* projection = reference_file->GetProjectionRef();
      outfile->SetGeoTransform(geotransform);
      outfile->SetProjection(projection);

      vector<GDALRasterBand*> outbands  = vector<GDALRasterBand*>(nBands);

      for(int i = 0; i < nBands; i++){
        outbands[i] = outfile->GetRasterBand(i+1);
      }

      for(int row = 0; row < nrow; row++){
        for(int k = 0; k < nBands; k++){
          bands_1[k]->RasterIO(GF_Read, xoffset, yoffset + row, ncol, 1,
                               &tile[k], ncol, 1,
                               GDT_Float64, sizeof(double)*(nBands*2), 0);
          bands_2[k]->RasterIO(GF_Read, xoffset, yoffset + row, ncol, 1,
                               &tile[k + nBands], ncol, 1,
                               GDT_Float64, sizeof(double)*(nBands*2), 0);
        }
        MapRMMatrixXd tileMat(tile, ncol, 2*nBands);
        MatrixRXd top = tileMat.block(0,0,ncol,nBands);
        MatrixRXd bot = tileMat.block(0,nBands,ncol,nBands);
        MatrixXd mads = (top * A) - (bot * B);
        for(int k = 1; k < nBands; k++){
          outbands[k]->RasterIO(GF_Write, xoffset, yoffset + row, ncol, 1,
                               &(mads.data()[k*ncol]), ncol, 1,
                               GDT_Float64, 0,0 );
          //implement chisqr to last band
          outbands[k]->FlushCache();
        }
        imad_utils::rowwise_divide(mads,sigMADs);
        //Take columnwise sum of squares, result is (1 x nBands) row vector
        VectorXd chisqr = mads.array().square().rowwise().sum().matrix();
      }
      delete[] geotransform;
    }
}
