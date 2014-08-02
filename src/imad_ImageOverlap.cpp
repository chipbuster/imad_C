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
    UL = Coord(myTrans.)
  }
}

namespace imad_ImageOverlap{

  void ImageOverlap(GDALDataset* input_file1,GDALDataset* input_file2,
                    int* winsize, int* offsets_1, int* offsets_2){

    CoordTransform img_1 = CoordTransform(input_file1);
    CoordTransform img_2 = CoordTransform(input_file2);

    BoundingBox box1 = BoundingBox(input_file1);
    BoundingBox box2 = BoundingBox(input_file2);

    // Case 1, the two images do not overlap, should throw an exception
    if(((box1.UL.x <= box2.UR.x) ||(box1.UR.x >= box2.UL.x))&&((box1.LL.y>= box2.LR.y)
      ||(box1.LR.y <= box2.LL.y))){
        throw std::invalid_argument("Error while opening file" + filename);
      }




  }


}
