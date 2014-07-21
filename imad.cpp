#include "gdal_priv.h"
#include "imad.h"
#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <Eigen/Dense>
#include <math.h>

inline int max(int a, int b){  return a > b ? a : b; }

using namespace Eigen;
using std::cout; //debugging purposes
using std::endl;

//typedef Map<Matrix<double,Dynamic,Dynamic,RowMajor> > MapRMMatrixXd;
//defined in imad.h

void imad(std::string filename1="", std::string filename2="",
          std::string output_file="", std::string format="",
          int* bands1_arg=NULL, int* bands2_arg=NULL){
  GDALAllRegister(); //Must be called at start, see GDAL API for details

  //openFile() will ask us for filenames if the names are blank
  GDALDataset* file1 = GdalFileIO::openFile(filename1);
  GDALDataset* file2 = GdalFileIO::openFile(filename2);
  if(!GdalFileIO::dimensionsmatch(file1,file2)) return; //function reports error

  // int* bands1, bands2;
  // int nBands = GdalFileIO::getBandInfo(bands1,bands2,bands1_arg, bands2_arg);

  GdalFileIO::getOutputFileInfo(output_file, format);

  int ncol = file1->GetRasterXSize();
  int nrow = file1->GetRasterYSize();

  //Temporary hardcoding of values. Eventually, will be user-specified
  int xoffset = 0;
  int yoffset = 0;
  int maxiter = 100; //Set to 1 to test output
  double tolerance = 0.0001;
  double pen = 0.0; //penalty

  /* selectBands() gets a series of bands from the user. For example, for
  a LandSat7 image, the vector might return [1,2,3,4,5,7] */
  std::vector<int>& bandnums = *GdalFileIO::selectBands();
  int nBands = bandnums.size();

  GDALRasterBand** bands_1  = new GDALRasterBand*[nBands];
  GDALRasterBand** bands_2  = new GDALRasterBand*[nBands];
  double* noDataValues     = new double[nBands];

  //Get easy handles to raster bands (pointers)
  for(int i = 0; i < nBands; i++){
    bands_1[i]      = (file1->GetRasterBand(bandnums[i]));
    bands_2[i]      = (file2->GetRasterBand(bandnums[i]));
    noDataValues[i] = bands_1[i]->GetNoDataValue();
  }

  //Input state is done. Set up for calculations.

  /* See notes at bottom for explanation of tile matrix. We set up the chunking
   * procedure for large images. In most cases, bufsize will be ncol and nb_pr
   * (number of buffers per row) will be 1. For large, large images, you might
   * need more than that. find_chunksize() modifies bufsize in-place (by ref)*/

/* Part of chunking---not yet implemented
  float* tile;
  int bufrsize = -1; //size of one row in the tile buffer
  int nb_pr = find_chunksize(bufsize, ncol, nBands);
*/
  double* tile = new double[2 * nBands * ncol];

  ImageStats cpm = ImageStats(nBands * 2);
  double delta = 1.0;
  VectorXd rho, oldrho = VectorXd::Zero(nBands);
  MatrixXd A, B;
  VectorXd  sigMADs, means1, means2;
  //Declare a special vector for calculations later on
  const VectorXd one = VectorXd::Constant(nBands, 1); //col vector of all "1"

  /* Start calculations (note double conditional in for-loop: we will continue
   * until the delta is smaller than the tolerance or to maxiter iterations) */
  for(int iter = 0; iter < maxiter && delta > tolerance; iter++){
    for(int row = 0; row < nrow; row++){
      for(int k = 0; k < nBands; k++){
        //Read into the appropriate location in tile.
        //The Python code for this might be easier to understand.
        bands_1[k]->RasterIO(GF_Read, xoffset, yoffset + row, ncol, 1,
                             &tile[k], ncol, 1,
                             GDT_Float64, sizeof(double)*(nBands*2), 0);
        bands_2[k]->RasterIO(GF_Read, xoffset, yoffset + row, ncol, 1,
                             &tile[k + nBands], ncol, 1,
                             GDT_Float64, sizeof(double)*(nBands*2), 0);
      }
      //The image data for a single row of all bands is now in tile.
      //Update the mean and covariance matrix by a provisional algorithm

      if(iter > 0){ //Repeated iterations use the previous iterations as weights
        VectorXd weights(ncol);
        weights = imad_utils::calc_weights(tile, weights, A, B, means1, means2,
                                                        sigMADs, ncol, nBands);
        cpm.update(tile, weights.data(), ncol, 2*nBands);
      }
      else{
        cpm.update(tile, NULL, ncol, 2*nBands);
      }
    }
    MatrixXd S = cpm.get_covar();
    VectorXd means = cpm.get_means();
    cpm.reset();
    MatrixXd s11 = S.block(0,0,nBands,nBands);           //Upper left
    MatrixXd s12 = S.block(0,nBands,nBands,nBands);      //Upper right
    MatrixXd s21 = S.block(nBands,0,nBands,nBands);      //Lower left
    MatrixXd s22 = S.block(nBands,nBands,nBands,nBands); //Lower right
    s11 = (1 - pen) * s11  + pen * MatrixXd::Identity(nBands,nBands);
    s22 = (1 - pen) * s22  + pen * MatrixXd::Identity(nBands,nBands);
    MatrixXd b1  = s11;
    MatrixXd c1  = s12 * s22.inverse() * s21;
    MatrixXd b2  = s22;
    MatrixXd c2  = s21 * s11.inverse() * s12;

    //Solve the eigenproblem
    VectorXd mu2, mu;
    if(nBands == 1){
      mu2(0) = c1(0,0)/b1(0,0);
      A(0,0) = 1 / sqrt(b1(0,0));
      B(0,0) = 1 / sqrt(b2(0,0));
    }
    else{
      Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> solver1(c1,b1);
      Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> solver2(c2,b2);
      if (solver1.info() != Eigen::Success || solver2.info() != Eigen::Success){
        throw std::runtime_error("Eigensolver did not converge!");
      }
      mu2 = solver2.eigenvalues(); //Eigenvalues should be identical, use either
      A = solver1.eigenvectors();
      B = solver2.eigenvectors();
    }

    /* We now need to sort the eigenvals and their eigenvecs. mu2 should be
     * sorted in desc. order and eigenvecs (A and B) reordered appropriately*/

    imad_utils::reorder_eigens(mu2, A, B);

    mu = mu2.array().sqrt().matrix(); //Element-wise square root (mu2 = rho**2)

    //Calculate the variances of the eigenvectors (Var(e1), Var(e2), etc.)
    //Note that in Python, these are row vectors, while they are column in C
    VectorXd eigvarA = (A.transpose() * A).diagonal();
    VectorXd eigvarB = (B.transpose() * B).diagonal();

    /* Calculate the penalized MAD significances. Check the README for more
     * information on the source of the variable names */

    //if pen == 0, sigma <-- 2*(1 - mu) and rho <-- mu

    VectorXd sigma = ((2*one - pen*(eigvarA + eigvarB))          //numerator
                                    /                            //---------
                    (1 - pen) - 2 * mu).array().sqrt().matrix(); //denominator

    VectorXd rho = (
      mu.array() * (1 - pen)                                       //numerator
                 /                                                 //---------
((one - pen*eigvarA).array() * (one - pen*eigvarB).array()).sqrt() //denominator
                   ).matrix();

    delta = (rho - oldrho).maxCoeff(); //The max of diffs between correlations
    oldrho = rho;

    cout << "Iteration: " << iter << "  --   Delta: " << delta << endl;

    sigMADs = sigma;
    means1 = means.block(0,0,nBands,1);
    means2 = means.block(nBands,0,nBands,1);

    imad_utils::math_cleanup(A,B,s11,s12);

    }

  //End iterations. Gear up to write final result to file.

  GDALDriver* outdriver=GetGDALDriverManager()->GetDriverByName(format.c_str());
  GDALDataset *outfile = outdriver->Create(output_file.c_str(),ncol,nrow,
                                           nBands + 1, GDT_Float64, NULL);

  GdalFileIO::writeOutputToFile(outfile, tile, A, B, xoffset, yoffset, ncol,
                                nrow, bands_1, bands_2, file1, sigMADs);

  //Cleanup in aisle 7!
  GDALClose( (GDALDatasetH) file1 );
  GDALClose( (GDALDatasetH) file2 );
  GDALClose( (GDALDatasetH) outfile);

  delete &bandnums;
  delete[] tile;
}

int main(){ //dummy main
  std::string file1 = std::string("/home/chipbuster/POP_2014/lndsr.LE71960531999293EDC00.tif");
  std::string file2 = std::string("/home/chipbuster/POP_2014/lndsr.LE71960532000120EDC00.tif");
  std::string out = std::string("/home/chipbuster/POP_2014/C++_asd.tif");
  std::string fmt = std::string("GTiff");
  imad(file1,file2,out,fmt);
  return 0;
}
