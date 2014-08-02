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
  Coord UL; // coord has x and y
  Coord UR;
  Coord LL;
  Coord LR;
  BoundingBox(GDALDataset* input){ // get coordinates of the corners easily
    CoordTransform myTrans = CoordTransform(input); // will change coordinate in place
    int xsize = input->GetRasterXSize();
    int ysize = input->GetRasterYSize();

    //List image corners in Pixel,Line coordinates
    UL = Coord(0,0);
    UR = Coord(xsize+1,0);
    LR = Coord(xsize+1,ysize+1);
    LL = Coord(0,ysize+1);

    //Transform corners to geographic coordinates
    myTrans.ImgtoGeo(UL);
    myTrans.ImgtoGeo(UR);
    myTrans.ImgtoGeo(LL);
    myTrans.ImgtoGeo(LR);
    //Done
  }

  bool inside(const BoundingBox& other){

  }
}

namespace imad_ImageOverlap{

  void ImageOverlap(GDALDataset* input_file1,GDALDataset* input_file2,
                    int* winsize, int* offsets_1, int* offsets_2){

    CoordTransform img_1 = CoordTransform(input_file1);
    CoordTransform img_2 = CoordTransform(input_file2);
    ImageInfo     info_1 = ImageInfo(input_file1);
    ImageInfo     info_2 = ImageInfo(input_file2);
    //Compatible() will give details if it fails
    if(! info_1.compatible(info_2)){
      throw std::invalid_argument("Images are not compatible!")
    }
    // Case 1, the two images do not overlap, should throw an exception
    if(((box1.UL.x <= box2.UR.x) ||(box1.UR.x >= box2.UL.x))&&((box1.LL.y>= box2.LR.y)
      ||(box1.LR.y <= box2.LL.y))){
        throw std::invalid_argument("Error while opening file" + filename);
      }


  }


}
