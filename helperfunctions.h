#ifndef HELPERS_H
#define HELPERS_H

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "ImageStats.h" //Includes <Eigen/Dense>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>

/* This file contains the headers for the iMad project. */
/* All GDAL functions derive from gdal_priv.h           */

//Found in openGdalFiles.cpp
GDALDataset* openFile();
GDALDataset* openFile(std::string filename);
bool dimensionsmatch(GDALDataset* dataset1, GDALDataset* dataset2);




#endif
