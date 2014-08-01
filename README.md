Iteratively-reweighted Multivariate Alteration Detection (iMAD)
===============================================================

This package implements Nielsen's iMAD algorithm. It is a C++ program
ported from Python, with the eventual goal of making an R package.

## Required Packages

See Building/Install section for details on how to acquire these packages.

To build iMad, you will need the GDAL libraries, which can be acquired from
the [OSGEO website](http://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries).

You will also need the Eigen libraries, which can be found
[here](http://eigen.tuxfamily.org/index.php?title=Main_Page). This
is the primary linear algebra library in this package.

The [Boost Libraries](http://www.boost.org/) provide some statistical functions
needed for iMad.

If you cannot get the Boost and Eigen libraries, this package comes with a
version preinstalled in the src/ directory. However, these may be out-of-date,
and you will still need to make sure to have the GDAL libraries installed.

## Building/Install

###On all systems

Make sure you have OSGEO/GDAL installed. Instructions for doing this on
specific systems can be found below.

If you want to install and use your own versions of the Boost and Eigen
libraries, edit the INCLUDEDIRS flag in the Makefile (experts only).

### On Linux

To install GDAL on an Ubuntu system, run
`sudo apt-get install libgdal1-1.7.0 libgdal1-dev` from the command line.

For RedHat systems (Fedora, RHEL, and CentOS), use
`yum install gdal-libs gdal-devel`.

On Arch Linux and derivatives, `sudo pacman -S gdal`.

I don't have any other systems to test on right now.

cd to the src/ directory and run `make`. There is no `make install`
at the moment.

### On other systems

Work in progress, hang tight!

## Use

I dunno. You guys are the community that came up with the algorithm, not me.
