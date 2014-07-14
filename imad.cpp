#include "gdal_priv.h"
#include "imad.h"
#include <string>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::ArrayXXf;
using Eigen::ArrayXf;

void imad(string filename1="", string filename2=""){
  GDALAllRegister(); //Must be called at start, see GDAL API for details

  //openFile() will ask us for filenames if the names are blank
  GDALDataset* file1 = GdalFileIO::openFile(filename1);
  GDALDataset* file2 = GdalFileIO::openFile(filename2);

  if(!GdalFileIO::dimensionsmatch(file1,file2)) return; //function reports error

  int ncol = file1->GetRasterXSize();
  int nrow = file1->GetRasterYSize();

  //Temporary hardcoding of values. Eventually, will be user-specified
  int xoffset = 0;
  int yoffset = 0;
  int maxiter = 10;
  double tolerance = 0.001;
  double pen = 0.0; //penalty

  /* selectBands() gets a series of bands from the user. For example, for
  a LandSat7 image, the vector might return [1,2,3,4,5,7] */
  vector<int>& bandnums = *GdalFileIO::selectBands();
  int nBands = bandnums.size();

  vector<GDALRasterBand*> bands_1  = vector<GDALRasterBand*>(nBands);
  vector<GDALRasterBand*> bands_2  = vector<GDALRasterBand*>(nBands);
  vector<double> noDataValues      = vector<double>(nBands);

  //Get easy handles to raster bands (pointers)
  for(size_t i = 0; i < bandnums.size(); i++){
    bands_1[i]     = (file1->GetRasterBand(bandnums[i]));
    bands_2[i]     = (file1->GetRasterBand(bandnums[i]));
    noDataValues[i] = bands_1[i]->GetNoDataValue();
  }

  //Input state is done. Set up for calculations.

  /* tile is a strange name for a matrix (borrowed from the python code).
   * It holds a single row of both images. The top half of the matrix holds
   * the row from file1, with all of its bands in separate rows (e.g. for any
   * given row, row 1 of tile will hold image1, band1, row 2 will hold image1,
   * band2, etc.). The same is true of the second image in the bottom
   * half of the matrix. */

  float* tile = new float[2 * nBands * ncol];
  ImageStats cpm = ImageStats(nBands * 2);
  double delta = 1.0;
  VectorXf rho, oldrho = VectorXf::Zero(nBands);
  MatrixXf sigMADs, means1, means2, A, B, cov;
  //Declare two special vectors
  const VectorXf two = VectorXf::Constant(nBands, 2); //col vector of all "2"
  const VectorXf one = VectorXf::Constant(nBands, 1); //col vector of all "1"

  /* Start calculations (note double conditional in for-loop: we will continue
   * until the delta is smaller than the tolerance or to maxiter iterations) */
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
    MatrixXf s11 = S.block(0,0,nBands,nBands);           //Upper left
    MatrixXf s12 = S.block(0,nBands,nBands,nBands);      //Upper right
    MatrixXf s21 = S.block(nBands,0,nBands,nBands);      //Lower left
    MatrixXf s22 = S.block(nBands,nBands,nBands,nBands); //Lower right
    MatrixXf b1  = s11;
    MatrixXf c1  = s12 * s22.inverse() * s21;
    MatrixXf b2  = s22;
    MatrixXf c2  = s21 * s11.inverse() * s22;

    //Solve the eigenproblem
    VectorXf mu2, mu;
    if(nBands == 1){
      mu2(0) = c1(0,0)/b1(0,0);
      A(0,0) = 1 / sqrt(b1(0,0));
      B(0,0) = 1 / sqrt(b2(0,0));
    }
    else{
      Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXf> solver1(c1,b1);
      Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXf> solver2(c2,b2);
      if (solver1.info() != Eigen::Success || solver2.info() != Eigen::Success){
        throw std::runtime_error("Eigensolver did not converge!");
      }
      mu2 = solver1.eigenvalues(); //Eigenvalues should be identical
      A = solver1.eigenvectors();
      B = solver2.eigenvectors();
    }

      /* We now need to sort the eigenvalues by their eigenvectors and return.
       * This is handled by a utility function. At the end of the function,
       * mu2 should be sorted in descending order, and the eigenvecs should be
       * in the appropriate order */

      imad_util::reorder_eigens(&mu2, &A, &B);

      mu = mu2.array().sqrt().matrix(); //All aboard the crazy train!

      //Calculate the variances of the eigenvectors (Var(e1), Var(e2), etc.)
      VectorXf eigvarA = (A * A.transpose()).diagonal();
      VectorXf eigvarB = (B * B.transpose()).diagonal();

      //Calculate the penalized MAD significances
      //If you do not understand the math, try mentally setting pen to zero
      VectorXf sigma = ((two - pen*(eigvarA + eigvarB)) / (1 - pen) - 2 * mu);

      //A lot of nasty dancing between array and matrix types
      VectorXf rho = (
        mu.array() * (1 - pen)                                       //numerator
                   /                                                 //---------
        ((one - pen*eigvarA) * (one - pen*eigvarB)).array().sqrt()  //denominator
                     ).matrix();

      delta = (rho - oldrho).maxCoeff();
      oldrho = rho;

      /*We now "tile" the sigmas and means */

      /* These last five lines were basically copied over verbatim from Python.
       * I stopped keeping track of what everything was for. Sorry :\
       * Note: Do expect overly-verbose commenting here. */

      /* We make a vector that is just the diagonal elements of s11. Then we
       * take the inverse square root (1/sqrt(x)) of those diagonal elements
       * and recast to a matrix. */

       //TODO: Rewrite and understand this section
        MatrixXf D = s11.diagonal().array().sqrt().cwiseInverse().matrix().asDiagonal();
        MatrixXf s = (D*s11*A).colwise().sum();
        A = A * (s.array()/s.array().abs()).matrix();
        cov = (A.transpose() * s12 * B).diagonal();
        B = B * (cov.array() / cov.array().abs()).matrix().asDiagonal();

    }

  //End iterations. Gear up to write final result to file.

  delete &bandnums;
  delete[] tile;
}

int main(){ //dummy main
  vector<int>* asd = GdalFileIO::selectBands();
  return 0;
}
