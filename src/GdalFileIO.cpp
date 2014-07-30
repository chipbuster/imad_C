#include "gdal_priv.h"
#include "imad.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <math.h>
#include <string.h> //Gdal libraries use const char*, need strcmp() to compare

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

    // Writes output to file. By default, in the C code, reference_file is file1

    void writeOutputToFile(GDALDataset* outfile, double* tile,
                           MatrixXd& A, MatrixXd& B, //Eigenvector matrices
                           int x10, int y10, int x20, int y20,
                           int ncol, int nrow,int nb_pr, int bufsize,
                           GDALRasterBand** bands_1,
                           GDALRasterBand** bands_2,
                           GDALDataset* reference_file,
                           VectorXd& sigMADs){
      int nBands = A.cols();
      double* geotransform   = new double[6];
      reference_file->GetGeoTransform(geotransform);
      //Move the origins of the picture to the overlap zone.
      geotransform[0] = geotransform[0] + x10*geotransform[1]
      geotransform[3] = geotransform[3] + y10*geotransform[5]
      const char* projection = reference_file->GetProjectionRef();
      outfile->SetGeoTransform(geotransform);
      outfile->SetProjection(projection);

      vector<GDALRasterBand*> outbands  = vector<GDALRasterBand*>(nBands + 2);

      for(int i = 0; i <= nBands + 1; i++){
        outbands[i] = outfile->GetRasterBand(i+1);
      }

      //Read array into memory
      int this_bufsize = 0;
      int xstart = 0;
      for(int row = 0; row < nrow; row++){
        for(int bufnum = 0; bufnum < nb_pr; bufnum++){
          for(int k = 0; k < nBands; k++){
            xstart = bufsize * bufnum;
            this_bufsize = min(bufsize, ncol - xstart);
            imad_bigfun::readToBuf((tile + k), bands_1[k],
                                   x10 + xstart, y10 + row,
                                   this_bufsize, nBands);
            imad_bigfun::readToBuf((tile + nBands + k), bands_2[k],
                                   x20 + xstart, y20 + row,
                                   this_bufsize, nBands);
          }
          //Compute MAD variates
          MapRMMatrixXd tileMat(tile, bufsize, 2*nBands);
          MatrixRXd top = tileMat.block(0,0,bufsize,nBands);
          MatrixRXd bot = tileMat.block(0,nBands,bufsize,nBands);
          MatrixXd mads = (top * A) - (bot * B);

          //Write MAD variates to file
          for(int k = 0; k < nBands; k++){
            outbands[k]->RasterIO(GF_Write, xstart, row, this_bufsize, 1,
                                 &(mads.data()[k*this_bufsize]), this_bufsize, 1,
                                 GDT_Float64, 0,0 );
            outbands[k]->FlushCache();
          }
          imad_utils::rowwise_divide(mads,sigMADs);
          //Take columnwise sum of squares, result is (1 x nBands) row vector
          VectorXd chisqr = mads.array().square().rowwise().sum().matrix();

          //Write ChiSqr values to extra band #1
          outbands[nBands]->RasterIO(GF_Write, xstart, row, this_bufsize, 1,
                                     chisqr.data(), this_bufsize, 1,
                                     GDT_Float64, 0, 0);
        }
      }
      delete[] geotransform;
    }

    /*************************************************************************/
/*   === The remainder of code in this file is error checking which is     ===
     === pathetically incomplete and full of holes. If you plan on using   ===
     === this code in another application (e.g. R, python), you should try ===
     === to write error checking code in that language (osgeo bindings in  ===
     === Python, rgdal in R). This code is not well maintained and should  ===
     === NOT be relied on to catch all errors in input.                    ===
*/


    /* If no output file, get an output file and format */
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

    // Called from fix_missing_band_data()
    int* selectBands(int nBands){
      cout << "You have specified " << nBands << " bands will be used." << endl;
      cout << "Please enter the band numbers you would like to use, separated "
           << "by a space. Press 'enter' when done." << endl;

      int* retval = new int[nBands];
      string input;
      cin >> input;
      for(int i = 0; i < nBands; i++){
        retval[i] = atof(input.c_str());
        cin >> input;
      }
      return retval;
    }

    /*************************************************************************/

    /* If explicit values are not given for bands (i.e. bands to use and number
     * of bands), this function will prompt the user to fill them */
    void fix_missing_band_data(int** bands1, int** bands2, int& nBands){
      if(nBands == -1){
        cout << "Please enter the number of bands you would like to test." << endl;
        cin >> nBands;
      }
      if (nBands < 1){
        throw std::out_of_range("Invalid number of bands.");
      }
      if(*bands1 == NULL){
        cout << "Selecting bands for image 1." << endl;
        *bands1 = selectBands(nBands);
      }
      if(*bands2 == NULL){
        cout << "Selecting bands for image 2." << endl;
        *bands2 = selectBands(nBands);
      }
    }
    /*************************************************************************/

    /* If explicit values are not given for dims (i.e. size of window and
     * image offsets), this function will prompt the user to fill them */
    void fix_missing_dims_data(int** win_size, int** offset_1,int** offset_2){
      //First, get window size
      int nrow,ncol;
      if(*win_size == NULL){
        cout << "Please enter the number of columns in the image (X-size): ";
        cin >> ncol;
        cout << "Please enter the number of rows in the image (Y-size): ";
        cin >> nrow;
        *win_size = new int[2];
        (*win_size)[0] = ncol; //use XSize x YSize instead of nrow x ncol
        (*win_size)[1] = nrow;
      }
      if(*offset_1 == NULL){
        (*offset_1) = new int[2];
        cout << "Please enter the X-offset for image 1: ";
        cin >> (*offset_1)[0];
        cout << "Please enter the Y-offset for image 1: ";
        cin >> (*offset_1)[1];
      }
      if(*offset_2 == NULL){
        (*offset_2) = new int[2];
        cout << "Please enter the X-offset for image 2: ";
        cin >> (*offset_2)[0];
        cout << "Please enter the Y-offset for image 2: ";
        cin >> (*offset_2)[1];
      }
    }

    /*************************************************************************/

    /* Checks input paramenters for all potential errors. This function should
     * never trigger in R or Python versions */

    bool has_errors(GDALDataset* file1, GDALDataset* file2,
                   int* bands1_arg, int* bands2_arg, int nBands,
                   int* win_size, int* offsets_1, int* offsets_2, int inp_pen ){

      // Check whether input bands are within limits.

      //offsets cannot be negative
      if(offsets_1[0] < 0 || offsets_2[0] < 0 ||
         offsets_1[1] < 0 || offsets_2[1] < 0 ){
        cout << "Error: image offsets cannot be negative!" << endl;
        return true;
      }

      //Image dims can't be negative either
      if(win_size[0] < 0 || offsets_2[1] < 0){
        cout << "Error: image dimensions cannot be negative!" << endl;
        return true;
      }

      //Xoffset + Xwidth cannot exceed image Xsize
      if(offsets_1[0] + win_size[0] > file1->GetRasterXSize()){
        cout << "Error in Image 1: Dimensions too large." << endl;
        cout << "X-offset + window width exceeds image x-dimension!" << endl;
        return true;
      }

      //Perform the same check for image 2
      if(offsets_2[0] + win_size[0] > file2->GetRasterXSize()){
        cout << "Error in Image 2: Dimensions too large." << endl;
        cout << "X-offset + window width exceeds image x-dimension!" << endl;
        return true;
      }

      //Yoffset + Ywidth cannot exceed image Xsize
      if(offsets_1[1] + win_size[1] > file1->GetRasterYSize()){
        cout << "Error in Image 1: Dimensions too large." << endl;
        cout << "Y-offset + window height exceeds image y-dimension!" << endl;
        return true;
      }

      //Perform the same check for image 2
      if(offsets_2[1] + win_size[1] > file2->GetRasterYSize()){
        cout << "Error in Image 2: Dimensions too large." << endl;
        cout << "Y-offset + window height exceeds image y-dimension!" << endl;
        return true;
      }

      //Image penalty must be between 0 and 1
      if(inp_pen > 1 || inp_pen < 0){
        cout << "Error: penalization value must be between 0 and 1!" << endl;
        return true;
      }

      if(nBands < 1){
        cout << "Error: must have more than zero bands!" << endl;
        return true;
      }

      return false;
    }


}
