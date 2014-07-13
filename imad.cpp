#include "gdal_priv.h"
#include "imad.h"
#include <string>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXf;
using Eigen::VectorXf;

void imad(string filename1="", string filename2=""){
  GDALAllRegister(); //Must be called before GDAL functions

  //openFile() will ask us for filenames if the names are blank
  GDALDataset* file1 = GdalFileIO::openFile(filename1);
  GDALDataset* file2 = GdalFileIO::openFile(filename2);

  if(!GdalFileIO::dimensionsmatch(file1,file2)) return; //function reports exact error

  int nrow = file1->GetRasterXSize();
  int ncol = file1->GetRasterYSize();

  //Temporary hardcoding. Eventually, will be user-specified
  int xoffset = 0;
  int yoffset = 0;
  int maxiter = 10;
  double tolerance = 0.001;
  double pen = 0.0; //penalty

  /* selectBands() gets a series of bands from the user. For example, for
  a LandSat7 image, the vector might return [1,2,3,4,5,7] */
  vector<int>& bandnums = *GdalFileIO::selectBands();
  int nBands = bandnums.size();

  vector<GDALRasterBand*> bands_1 = vector<GDALRasterBand*>(nBands);
  vector<GDALRasterBand*> bands_2 = vector<GDALRasterBand*>(nBands);
  vector<double> noDataValues      = vector<double>(nBands);

  //Get easy handles to raster bands (pointers)
  for(size_t i = 0; i < bandnums.size(); i++){
    bands_1[i]     = (file1->GetRasterBand(bandnums[i]));
    bands_2[i]     = (file1->GetRasterBand(bandnums[i]));
    noDataValues[i] = bands_1[i]->GetNoDataValue();
  }

  //Input state is done. Set up for calculations.

  /* tile is a strange name for a matrix (borrowed from the python code).
  It holds a single row of both images. The top half of the matrix holds
  the row from file1, with all of its bands in separate rows (e.g. for any
  given row, row 1 of tile will hold image1, band1, row 2 will hold image1,
  band2, etc.). The same is true of the second image in the bottom half of the
  matrix. */

  float* tile = new float[2 * nBands * ncol];
  ImageStats cpm = ImageStats(nBands * 2);
  double delta = 1.0;
  VectorXf rho, oldrho;
  MatrixXf sigMADs, means1, means2, A, B, cov;

  //Start calculations (note double conditional in for-loop)
  for(int iter = 0; iter < maxiter && delta > tolerance; iter++){
    for(int row = 0; row < nrow; row++){
      for(int k = 0; k < nBands; k++){
        //Read into the appropriate location in tile.
        //The Python code for this might be easier to understand.
        bands_1[k]->RasterIO(GF_Read, xoffset, yoffset, ncol, 1,
                             &tile[k * ncol], ncol, 1,
                             GDT_Float32, 0,0 );
        bands_1[k]->RasterIO(GF_Read, xoffset, yoffset, ncol, 1,
                             &tile[(k + nBands) * ncol], ncol, 1,
                             GDT_Float32, 0,0 );
      }
      if(iter > 0){} //Not sure how to implement this yet
      else cpm.update(tile, NULL, ncol, 2*nBands);
    }
    MatrixXf S = cpm.get_covar();
    VectorXf means = cpm.get_means();
    MatrixXf s11 = S.block(0,0,nBands,nBands);
    MatrixXf s12 = S.block(0,nBands,nBands,nBands);
    MatrixXf s21 = S.block(nBands,0,nBands,nBands);
    MatrixXf s22 = S.block(nBands,nBands,nBands,nBands);
    MatrixXf b1  = s11;
    MatrixXf c1  = s12 * s22.inverse() * s21;
    MatrixXf b2  = s22;
    MatrixXf c2  = s21 * s11.inverse() * s22;

    


  }

  delete &bandnums;
  delete[] tile;
}

int main(){ //dummy main
  return 0;
}
