#include "helperfunctions.h"

using namespace std;

GDALDataset* openFile(string filename){
  /* Returns a GDALDataset handle to the file given by file1.
   * If the file cannot be opened, an exception is thrown */

  if(filename.empty()){
    cout << "Please enter a file to be opened: ";
    getline(cin, filename);
  }

  GDALDataset* dataset1 = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);

  if(dataset1 == NULL){
    throw std::invalid_argument("Error while opening file" + filename);
  }
  return dataset1;
}

bool dimensionsmatch(GDALDataset* dataset1, GDALDataset* dataset2){
  /* Check that dimensions of the two files match in X, Y, and num. bands */
  if(dataset1->GetRasterXSize() != dataset2->GetRasterXSize()){
    cout << "X-dimension mismatch in images. Exiting..." << endl;
    return false;
  }
  if(dataset1->GetRasterYSize() != dataset2->GetRasterYSize()){
    cout << "Y-dimension mismatch in images. Exiting..." << endl;
    return false;
  }
  if(dataset1->GetRasterCount() != dataset2->GetRasterCount()){
    cout << "Raster count mismatch in images. Exiting..." << endl;
    return false;
  }
  //All dimensions match!
  return true;
}
