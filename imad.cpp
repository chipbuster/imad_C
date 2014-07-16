#include "gdal_priv.h"
#include "imad.h"
#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;

//typedef Map<Matrix<float,Dynamic,Dynamic,RowMajor> > MapRMMatrixXf;
//defined in imad.h

void imad(std::string filename1="", std::string filename2=""){
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
  std::vector<int>& bandnums = *GdalFileIO::selectBands();
  int nBands = bandnums.size();

  std::vector<GDALRasterBand*> bands_1  = std::vector<GDALRasterBand*>(nBands);
  std::vector<GDALRasterBand*> bands_2  = std::vector<GDALRasterBand*>(nBands);
  std::vector<double> noDataValues      = std::vector<double>(nBands);

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
   * half of tile. Note that this tile is transposed compared to the one in
   * Python: Pyversion has rows in cols and files split left/right. */

  float* tile = new float[2 * nBands * ncol];
  ImageStats cpm = ImageStats(nBands * 2);
  double delta = 1.0;
  VectorXf rho, oldrho = VectorXf::Zero(nBands);
  MatrixXf A, B, cov;
  VectorXf  sigMADs, means1, means2;
  //Declare two special vectors for calculations later on
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

      //Mapped row-major matrix, typedef'ed in imad.h, see Eigen API on "Map"
      MapRMMatrixXf tileMat(tile,2*nBands,ncol);

      if(iter > 0){
        MapRMMatrixXf top(tile, nBands, ncol);
        MapRMMatrixXf bot(tile + nBands*ncol, nBands, ncol);

        //Subtract off the means from the previous iteration (utility function)
        imad_utils::colwise_subtract(top, means1);
        imad_utils::colwise_subtract(bot, means2);
        //top and bot are now zero-mean matrices. Do some math to them.
        MatrixXf mads = top.transpose() * A - bot.transpose() * B;
        imad_utils::rowwise_divide(mads,sigMADs);
        //Take columnwise sum of squares, result is (1 x nBands) row vector
        RowVectorXf chisqr = mads.array().square().colwise().sum().matrix();
        RowVectorXf tmp(1,bandnums.size()); //Pass this into getWeights()
        RowVectorXf weights = one - imad_utils::getWeights(chisqr, tmp, bandnums);
        cpm.update(tileMat, weights, ncol, 2*nBands);
      }
      else{
        RowVectorXf tmp(0,0); //No-dims indicates we do not want to use weights
        cpm.update(tileMat, tmp, ncol, 2*nBands);
      }
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

      imad_utils::reorder_eigens(mu2, A, B);

      mu = mu2.array().sqrt().matrix(); //All aboard the crazy train!

      //Calculate the variances of the eigenvectors (Var(e1), Var(e2), etc.)
      //Note that in Python, these are row vectors, while they are column in C
      VectorXf eigvarA = (A * A.transpose()).diagonal();
      VectorXf eigvarB = (B * B.transpose()).diagonal();

      /* Calculate the penalized MAD significances. If you do not understand the
       * math, try mentally setting pen to zero and then matching this with the
       * equations from DOI:10.1016/j.rse.2003.10.024. The Greek variable names
       * have been chosen to line up with the variable names in that paper. */

      VectorXf sigma = ((two - pen*(eigvarA + eigvarB)) / (1 - pen) - 2 * mu);

      //A lot of nasty dancing between array and matrix types
      VectorXf rho = (
        mu.array() * (1 - pen)                                       //numerator
                   /                                                 //---------
        ((one - pen*eigvarA) * (one - pen*eigvarB)).array().sqrt()  //denominator
                     ).matrix();

      delta = (rho - oldrho).maxCoeff();
      oldrho = rho;

      /* In the python code, you tile (repeat) these vectors to get an array. To
       * save memory for large images, we will use the funcs in imad_utils to do
       * column/rowwise ops, saving the cost of huge repeating array in memory*/

      sigMADs = sigma;
      means1 = means.block(0,0,nBands,1);
      means2 = means.block(nBands,0,nBands,1);

      /* What follows are some of the most poorly-commented lines from the
       * original Python implementation, and my understanding of them is still
       * fuzzy. I have included modified versions of the original comments, and
       * the lines are each numbered, with comment-footnotes at the end of the file */

/* 1 */ MatrixXf D = s11.diagonal().array().sqrt().cwiseInverse().matrix().asDiagonal();
/* 2 */ MatrixXf s = (D*s11*A).colwise().sum();
/* 3 */ A = A * (s.array()/s.array().abs()).matrix();
/* 4 */ cov = (A.transpose() * s12 * B).diagonal();
/* 5 */ B = B * (cov.array() / cov.array().abs()).matrix().asDiagonal();

    }

    //FINISHED WITH THE COMPUTATION CODE! LET'S GO GET A DRINK!

  //End iterations. Gear up to write final result to file.

  delete &bandnums;
  delete[] tile;
}

int main(){ //dummy main
  std::vector<int>* asd = GdalFileIO::selectBands();
  delete asd;
  return 0;
}
/**** Notes on lines 1-5 in imad() ****/
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
