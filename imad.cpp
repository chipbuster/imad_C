#include "helperfunctions.h"
#include "gdal_priv.h"

void imad(std::string filename1="", std::string filename2=""){
  //openFile() will ask us for filenames if the names are blank
  GDALDataset* file1 = openFile(filename1);
  GDALDataset* file2 = openFile(filename2);

  

}
