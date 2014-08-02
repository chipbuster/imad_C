#include "gdal_priv.h"
#include "imad.h"
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <boost/math/distributions.hpp> //for cdf()

using namespace std;

//BoundingBox is only ever used locally, so keep def'n in-file
struct BoundingBox{
  Coord UL;
  Coord LR;
  BoundingBox(GDALDataset* input){
    CoordTransform myTrans = CoordTransform(input);
    int xsize = input->GetRasterXSize();
    int ysize = input->GetRasterYSize();

    //List image corners in Pixel,Line coordinates
    UL = Coord(0,0);
    LR = Coord(xsize+1,ysize+1);

    //Transform corners to geographic coordinates
    myTrans.ImgtoGeo(UL);
    myTrans.ImgtoGeo(LR);
    //Done!
  }
}

namespace imad_ImageOverlap{

  void ImageOverlap(GDALDataset* input_file1,GDALDataset* input_file2,
                    int* winsize, int* offsets_1, int* offsets_2){

    CoordTransform img_1 = CoordTransform(input_file1);
    CoordTransform img_2 = CoordTransform(input_file2);




  }


}
