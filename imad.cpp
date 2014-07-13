#include "gdal_priv.h"
#include "imad.h"
#include <string>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

void imad(std::string filename1="", std::string filename2=""){
  GDALAllRegister(); //Must be called before GDAL functions

  //openFile() will ask us for filenames if the names are blank
  GDALDataset* file1 = GdalFileIO::openFile(filename1);
  GDALDataset* file2 = GdalFileIO::openFile(filename2);

  if(!GdalFileIO::dimensionsmatch(file1,file2)) return; //function reports exact error

  int nrow = file1->GetRasterXSize();
  int ncol = file1->GetRasterYSize();

  //Temporary hardcoding. Eventually, offsets will be user-specified
  int xoffset = 0;
  int yoffset = 0;

  std::vector<int>& bandnums = *GdalFileIO::selectBands();
  int nBands = bandnums.size();

//   GDALRasterBand& bands_1 = new GDALRasterBand[bandnums.size()];
//   GDALRasterBand& bands_2 = new GDALRasterBand[bandnums.size()];
//
// //Get Raster band objects and place in array
//   for(size_t i = 0; i < bandnums.size(); i++){
//     bands_1[i] = file1.GetRasterBand(bandnums[i]);
//     bands_2[i] = file2.GetRasterBand(bandnums[i]);
//   }

}

int main(){ //dummy main
  return 0;
}
