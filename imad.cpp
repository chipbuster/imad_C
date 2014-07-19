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
          std::string output_file="", std::string format=""){
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
  int maxiter = 2; //Set to 1 to test output
  double tolerance = 0.00001;
  double pen = 0.0; //penalty

  /* selectBands() gets a series of bands from the user. For example, for
  a LandSat7 image, the vector might return [1,2,3,4,5,7] */
  std::vector<int>& bandnums = *GdalFileIO::selectBands();
  int nBands = bandnums.size();

  std::vector<GDALRasterBand*> bands_1  = std::vector<GDALRasterBand*>(nBands);
  std::vector<GDALRasterBand*> bands_2  = std::vector<GDALRasterBand*>(nBands);
  std::vector<double> noDataValues      = std::vector<double>(nBands);

  //Get easy handles to raster bands (pointers)
  for(size_t i = 0; i < bandnums.size(); i++){
    bands_1[i]     = (file1->GetRasterBand(bandnums[i]));
    bands_2[i]     = (file2->GetRasterBand(bandnums[i]));
    noDataValues[i] = bands_1[i]->GetNoDataValue();
  }

  if(output_file.empty()){
    std::cout<< "Valid formats: GTiff, PCIDSK, HFA, ENVI" << std::endl;
    std::cout<< "Please enter output format: ";
    getline(std::cin, format);
    std::cout<< "Please enter output file name";
    getline(std::cin, output_file);
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
  std::cout<< "Image is " << nrow << " by " << ncol << std::endl;

  ImageStats cpm = ImageStats(nBands * 2);
  double delta = 1.0;
  VectorXd rho, oldrho = VectorXd::Zero(nBands);
  MatrixXd A, B;
  VectorXd  sigMADs, means1, means2;
  //Declare two special vectors for calculations later on
  const VectorXd two = VectorXd::Constant(nBands, 2); //col vector of all "2"
  const VectorXd one = VectorXd::Constant(nBands, 1); //col vector of all "1"

  /* Start calculations (note double conditional in for-loop: we will continue
   * until the delta is smaller than the tolerance or to maxiter iterations) */
  for(int iter = 0; iter < maxiter && delta > tolerance; iter++){
    for(int row = 0; row < nrow; row++){
      std::cout<< "Row " << row << " of " << nrow << std::endl;
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

      if(iter > 0){
        //Mapped row-major matrix, typedef'ed in imad.h, see Eigen API on "Map"
        MapRMMatrixXd tileMat(tile, ncol, 2*nBands);
        MatrixRXd top = tileMat.block(0,0,ncol,nBands);
        MatrixRXd bot = tileMat.block(0,nBands,ncol,nBands);

        imad_utils::rowwise_subtract(top, means1, bot, means2);
        //Subtract off the means from the previous iteration (utility function)
        /* top and bot are now zero-mean matrices. Mads must be row-major for
         * output purposes later on. */

        MatrixXd mads = (top * A) - (bot * B);
        imad_utils::rowwise_divide(mads,sigMADs);
        //Take columnwise sum of squares, result is (1 x nBands) row vector
        VectorXd chisqr = mads.array().square().rowwise().sum().matrix();
        VectorXd tmp(ncol,1); //Pass this into getWeights()
        VectorXd weights = imad_utils::getWeights(chisqr, tmp, nBands);
        cpm.update(tile, weights.data(), ncol, 2*nBands);
      }
      else{
        cpm.update(tile, NULL, ncol, 2*nBands);
      }
    }
    std::cout<< "Done updating!" << std::endl;
    MatrixXd S = cpm.get_covar();
    VectorXd means = cpm.get_means();
    cpm.reset();
    cout << means.transpose() << endl << S << endl;
    MatrixXd s11 = S.block(0,0,nBands,nBands);           //Upper left
    MatrixXd s12 = S.block(0,nBands,nBands,nBands);      //Upper right
    MatrixXd s21 = S.block(nBands,0,nBands,nBands);      //Lower left
    MatrixXd s22 = S.block(nBands,nBands,nBands,nBands); //Lower right
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
      mu2 = solver1.eigenvalues(); //Eigenvalues should be identical
      A = solver1.eigenvectors();
      B = solver2.eigenvectors();
    }

      /* We now need to sort the eigenvalues by their eigenvectors and return.
       * This is handled by a utility function. At the end of the function,
       * mu2 should be sorted in descending order, and the eigenvecs should be
       * in the appropriate order */

      imad_utils::reorder_eigens(mu2, A, B);

      mu = mu2.array().sqrt().matrix(); //All aboard the crazy train!

      //Calculate the variances of the eigenvectors (Var(e1), Var(e2), etc.)
      //Note that in Python, these are row vectors, while they are column in C
      VectorXd eigvarA = (A.transpose() * A).diagonal();
      VectorXd eigvarB = (B.transpose() * B).diagonal();

      /* Calculate the penalized MAD significances. If you do not understand the
       * math, try mentally setting pen to zero and then matching this with the
       * equations from DOI:10.1016/j.rse.2003.10.024. The Greek variable names
       * have been chosen to line up with the variable names in that paper. */

      VectorXd sigma = ((two - pen*(eigvarA + eigvarB))
                                      /
                            (1 - pen) - 2 * mu).array().sqrt().matrix();
      //A lot of nasty dancing between array and matrix types
      VectorXd rho = (
        mu.array() * (1 - pen)                                       //numerator
                   /                                                 //---------
((one - pen*eigvarA).array() * (one - pen*eigvarB).array()).sqrt()  //denominator
                     ).matrix();



      delta = (rho - oldrho).maxCoeff(); //The max of diffs between correltaions
      oldrho = rho;

      /* In the python code, you tile (repeat) these vectors to get an array. To
       * save memory for large images, we will use the funcs in imad_utils to do
       * column/rowwise ops, saving the cost of huge repeating array in memory*/

      sigMADs = sigma;
      means1 = means.block(0,0,nBands,1);
      means2 = means.block(nBands,0,nBands,1);

  /* What follows are some of the most poorly-commented lines from the original
   * Python implementation, and my understanding of them is still fuzzy.
   * I have included modified versions of the original comments, and the
   * lines are each numbered, with comment-footnotes at the end of the file */

   //You are not expected to understand these lines. I certainly don't.

/* 1 */ MatrixXd D = s11.diagonal().array()      //Bookkeeping
                        .sqrt().cwiseInverse()   //Inverse sqrt
                        .matrix().asDiagonal();  //More bookkeeping
/* 2 */ VectorXd s = (D*s11*A).rowwise().sum();
        VectorXd tmp1 = (s.array() / (s.array().abs())).matrix();
/* 3 */ imad_utils::colwise_multiply(A, tmp1);
/* 4 */ VectorXd cov = (A.transpose() * s12 * B).diagonal();
/* 5 */ MatrixXd tmp2 = (cov.array() / cov.array().abs()).matrix().asDiagonal();
        B = B * tmp2;

    }

  //FINISHED WITH THE COMPUTATION CODE! LET'S GO GET A DRINK!

  //End iterations. Gear up to write final result to file.

  GDALDriver* outdriver=GetGDALDriverManager()->GetDriverByName(format.c_str());
  GDALDataset *outfile = outdriver->Create(output_file.c_str(),ncol,nrow,
                                           nBands + 1, GDT_Float64, NULL);
  const char* projection = file1->GetProjectionRef();
  double* geotransform   = new double[6];
  file1->GetGeoTransform(geotransform);
    //Move the origins of the picture to the overlap zone. Not implemented yet.
        // gt[0] = gt[0] + x10*gt[1]
        // gt[3] = gt[3] + y10*gt[5]
  outfile->SetGeoTransform(geotransform);
  outfile->SetProjection(projection);

  std::vector<GDALRasterBand*> outbands  = std::vector<GDALRasterBand*>(nBands);
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

  //Cleanup in aisle 7!
  GDALClose( (GDALDatasetH) file1 );
  GDALClose( (GDALDatasetH) file2 );
  GDALClose( (GDALDatasetH) outfile);
  delete &bandnums;
  delete[] geotransform;
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

/**** Notes on tile matrix ****/
/* tile is a strange name for a matrix (borrowed from the python code).
 * It holds a single row of both images. The top half of the matrix holds
 * the row from file1, with all of its bands in separate rows (e.g. for any
 * given row, row 1 of tile will hold image1, band1, row 2 will hold image1,
 * band2, etc.). The same is true of the second image in the bottom
 * half of tile. Note that this tile is transposed compared to the one in
 * Python: Pyversion has rows in cols and files split left/right. */

/**** Notes on numbered lines in imad() ****/
/* While going through these, it is helpful to remember that left multiplying by
 * a diagonal matrix (D*A) multiplies each row by the diagonal element, i.e.
 * row1(A) = row1(A) * D(1,1), row2(A) = row2(A) * D(2,2), etc. and right-multiplying
 * by a diagonal matrix multiplies the columns instead of the rows. Off we go!*/

/* 1 */ /* Completely unintelligible. First, we take the diagonal of s11. We then
 * cast it to an array (to deal with some issues with eigen), then take an
 * element-wise inverse square root. We then revert it back to a diagonal matrix.
 * The end effect is identical to masking all but the diagonal of s11, then
 * taking an element-wise inverse square root. */

/* 2 */ /* s is actually a row-vector. Not sure what this actually does */

/* 3 */ /* If S(i,i) was negative, we set it to -1. Otherwise, we set it to +1.
 * We then multiply A by this matrix to ensure the correlation sum is positive
 * (Yes, this part feels like magic to me too). */

/* 4,5 */ /* A similar process for B */
