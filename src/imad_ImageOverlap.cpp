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
  BoundingBox(GDALDataset* input){ //Get coordinate in pixel/line
    int xsize = input->GetRasterXSize();
    int ysize = input->GetRasterYSize();
    //List image corners in Pixel,Line coordinates
    UL = Coord(0,0);
    UR = Coord(xsize,0);
    LR = Coord(xsize,ysize);
    LL = Coord(0,ysize);
  }
  //This constructor places the coordinates relative to the first (i.e. assumes
  //that the origin of the reference file is 0,0). These are still P,L coords
  //and some of the coords may be negative.
  BoundingBox(GDALDataset* input, GDALDataset* reference){
    //Files need to have identical geotransforms for this to work
    CoordTransform myT_inp = CoordTransform(input);
    CoordTransform myT_ref = CoordTransform(reference);

    //Transform input (P,L) coords to geo, then to ref (P,L) coords
    double ipxg = myT_inp.ImgtoGeo_X(0,0); //Input's x coord in geo
    double ipyg = myT_inp.ImgtoGeo_Y(0,0); //Input's y coord in geo

    int x_00 = static_cast<int>(myT_inp.GeotoImg_X(ipxg,ipyg));
    int y_00 = static_cast<int>(myT_inp.GeotoImg_Y(ipxg,ipyg));

    int xsize = input->GetRasterXSize();
    int ysize = input->GetRasterYSize();

    UL = Coord(x_00, y_00);
    UR = Coord(xsize,0) + UL;
    LR = Coord(xsize,ysize) + UL;
    LL = Coord(0,ysize) + UL;

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
