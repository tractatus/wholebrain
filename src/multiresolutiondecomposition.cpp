// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

// FFTW3
#include "fftw3.h"

//for wavelet
#include <ctime>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <algorithm>
#include "wavelet2d.h"
#include <math.h>
#include <iostream>

//for sleep
#include <unistd.h>

#include <Rcpp.h>

#include "wholebrain.h"

using namespace cv;
using namespace Rcpp;




/* function to compute the maximum from an image array */
void* maxv(vector<vector<double> > &arr, double &max){
    max = 0;


    for (unsigned int i =0; i < arr.size(); i++) {
        for (unsigned int j =0; j < arr[0].size(); j++) {
            if (max <= arr[i][j]){
                max = arr[i][j];
            }
        }
    }
    return 0;
}

//find nearest value in a vector
//std::lower_bound: returns the first value that does not compare less
//std::upper_bound: returns the first value that compares strictly greater
int nearest(std::vector<int> const& vec, int value) {
 //auto const it = std::upper_bound(vec.begin(), vec.end(), value);
 
 auto const upper = std::lower_bound(vec.begin(), vec.end(), value);
 auto const lower = upper-1;
    
 if (upper == vec.end()) { return -1; }

 return *lower;
}


RcppExport SEXP getMaxMin(SEXP input) {
  BEGIN_RCPP
  Rcpp::RNGScope __rngScope; 

  Rcpp::CharacterVector fname(input);
  
  std::string ffname(fname[0]);
  Rcpp::Rcout << "Loading image:" << ffname << std::endl;
  Mat src = imread(ffname, -1); // -1 tag means "load as is"
  double minVal;
  double maxVal;
  minMaxLoc(src, &minVal, &maxVal);

  return List::create(
    _["max"] = maxVal,
    _["min"] = minVal
  );
  
  END_RCPP  
}

/*
template<> SEXP wrap(const cv::Mat &obj) {

}
*/
/* show a image */
RcppExport SEXP multiresolutiondecomposition(SEXP input, SEXP scales, SEXP family, SEXP outputfile, SEXP cellBodies, SEXP computeenergy, SEXP computecoherency, SEXP computeorientation, SEXP maxV, SEXP minV, SEXP maskorig, SEXP roiOnly) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //define wavelet scales
  int J = Rcpp::as<int>(scales);
  int compEnergy = Rcpp::as<int>(computeenergy); //used as Sigma, 10 is recommended
  int cellB = Rcpp::as<int>(cellBodies);

  int maxNorm = Rcpp::as<int>(maxV);
  int minNorm = Rcpp::as<int>(minV);
  int maskOrig = Rcpp::as<int>(maskorig);
  bool saveroiOnly = Rcpp::as<bool>(roiOnly);

  //fluorescence threshold to minimize false positives
  //int fluorThresh = Rcpp::as<int>(fluorescenceThreshold);
  //std::vector<int> fluorThresh  = Rcpp::as< std::vector<int> >(fluorescenceThreshold);
  double maxVal, minVal;

  Mat tobefiltered; // no local
Mat cellImg;
Mat processImg;

//to handle execution time logging
clock_t t0, t1;
 

  Rcpp::CharacterVector f(input);
  std::string ff(f[0]);
  Rcpp::Rcout << "Loading image:" << ff << std::endl;
  t0 = clock();
  tobefiltered = imread(ff, -1); //CV_LOAD_IMAGE_GRAYSCALE
  //for debugging  imgSWT('/Users/danielfurth/Documents/mouse_atlas/lighsheet_figure/analysis/Resliceofanalysis.tif', cell.bodies=4, fluorescence.threshold=0, areaSize=0)->mupp2
  t1 = clock();
  double time_elapsed = (double)((t1 - t0) / CLOCKS_PER_SEC);
  Rcpp::Rcout << "LOADED." << " loading took " <<  time_elapsed << " seconds." << std::endl;

  //get width and height of original image.
  int rowsWzero = tobefiltered.rows;
  int colsWzero = tobefiltered.cols;

  //check that the image is symmetric. If not then subsample the image by bilinear interpolation.
  bool sumbsample;
  if(rowsWzero!=colsWzero){
    if(rowsWzero>colsWzero){
      Rcpp::Rcout << "Image size "<< colsWzero << "x" << rowsWzero << " is not symmetric, subsampling will be performed." << std::endl;
      Size downsampledsize(static_cast<int>(colsWzero), static_cast<int>(colsWzero));
      resize(tobefiltered, tobefiltered, downsampledsize, INTER_LINEAR);
      sumbsample = true;
    }else{
      Rcpp::Rcout << "Image size "<< colsWzero << "x" << rowsWzero << " is not symmetric, subsampling will be performed." << std::endl;
      Size downsampledsize(static_cast<int>(rowsWzero), static_cast<int>(rowsWzero));
      resize(tobefiltered, tobefiltered, downsampledsize, INTER_LINEAR);
      sumbsample = true;
    }
  }else{
    sumbsample = false;
  }

  //get width and height of image.
  int rowsW = tobefiltered.rows;
  int colsW = tobefiltered.cols;

  //set the different wavelet resolutions for cases where the image is not sampled by a factor of 2
  int myints[] = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
  std::vector<int> vec (myints, myints + sizeof(myints) / sizeof(int) );  

  //resize image to a J size if it is not contained with 2^J sampling.
  if(std::find(vec.begin(), vec.end(), rowsW) != vec.end()) {
    /* resolutions contains rowsW */
  } else {
    /* resolutions does not contain rowsW */
    Rcpp::Rcout << "Image size "<< rowsW << "x" << rowsW << " needs to be adjused to 2^J sampling." << std::endl;
    //save the old row for future reference.
    int oldRowW = rowsW;
    //function nearest defined on line 60.
    int newwidth = nearest(vec, rowsW);
    Size oldsize = tobefiltered.size();
    Size newsize(static_cast<int>(newwidth), static_cast<int>(newwidth));
    resize(tobefiltered, tobefiltered, newsize, INTER_LINEAR);
    //get width and height of image.
    rowsW = tobefiltered.rows;
    colsW = tobefiltered.cols;
     Rcpp::Rcout << "New size: "<< rowsW << "x" << rowsW << "." << std::endl;
      sumbsample = true;
  }

  

  //make a 1D vector of the image.
  vector<vector<double> > vec1(rowsW, vector<double>(colsW));

  int k =1;
  for (int i=0; i < rowsW; i++) {
    for (int j =0; j < colsW; j++){
      unsigned char temp;
      temp = ((uchar*) tobefiltered.data + i * tobefiltered.step)[j  * tobefiltered.elemSize() + k ];
      vec1[i][j] = (double) temp;
    }
  }

  //create output vector
  vector<double> output;
 
  //select wavelet
  Rcpp::CharacterVector fam(family);
  std::string famwav(fam[0]);
  string nm = famwav; // default "db2"

  //define wavelet scales
  //int J = Rcpp::as<int>(scales);// Rcpp::as<int >(scales); //default is J=6

  //this is the computationally heavy part. Improve swt_2d with OpenCL
  Rcpp::Rcout << "Beginning wavelet transform" << std::endl;
  t0 = clock();
  swt_2d(vec1,J,nm,output);
  t1 = clock();
  time_elapsed = (double)((t1 - t0) / CLOCKS_PER_SEC);
  Rcpp::Rcout << "DONE." << " wavelet took " <<  time_elapsed << " seconds." << std::endl;


  // get the output dimensions
  int rowW,colW;
  dwt_output_dim(vec1, rowW, colW );

  //blur is a matrix containing the appoximations coefficients of the last filter stage.
  Rcpp::Rcout << "Extracting approximation coefficients" << std::endl;
  vector<vector<double> > blur(rowW, vector<double>(colW));

  for (int i=0;i < rowW; i++){
    for (int j=0; j < colW;j++){
      double temp = output[i*colW + j];
      blur[i][j]= temp;
    }
  }

  double max;
  maxv(blur,max);

  //create OpenCv matrix to store the approximation coefficients in.
  Size imgSize = tobefiltered.size(); // size of output image

  Mat cvImg(imgSize, tobefiltered.depth()); 


  //extract approximations coefficients from blur martix
  for (int i = 0; i < imgSize.height; i++ ) {
    for (int j = 0; j < imgSize.width; j++ ){
        if ( blur[i][j] <= 0.0){
          blur[i][j] = 0.0;
        }

        ((uchar*)(cvImg.data + i * cvImg.step))[j] = (char) ( (blur[i][j] / max) * 255.0);
    }
  }

  //save approximation coefficients
  //get output filename suffix
  Rcpp::CharacterVector of(outputfile);
  std::string off(of[0]);

  string input;
  input = "approximations_coefficents/a_" + off + ".tif";
  try {
    if(!saveroiOnly){
      imwrite(input, cvImg);
      Rcpp::Rcout << input << " SAVED" << endl;
    }

  }
  catch (runtime_error& ex) {
    Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
    return(R_NilValue);
  }

  //define detail coefficients
  vector<vector<double> >  detail(1*rowW,vector<double>(J * colW));
  //extract detail coefficients from output martix
   Rcpp::Rcout << "Extracting detail coefficients" << std::endl;
  for (int k=0; k < J; k++) {
    for (int i=0; i < rowW; i++) {
      for(int j=0+ k*colW; j < (k+1)*colW; j++) {
        double temp = (output[(3*k+1)*rowW*colW+ i * colW +j - k*colW]+output[(3*k+2)*rowW*colW+ i * colW +j - k*colW]+output[(3*k+3)*rowW*colW+ i * colW +j - k*colW]);
        detail[i][j]= temp;
      }
    }
  }

  IplImage *dvImg; // image used for output
  CvSize imgSz; // size of output image

  imgSz.width = J*colW;
  imgSz.height = 1*rowW;

  dvImg = cvCreateImage( imgSz, 16, 1 );

  /*for (int i = 0; i < imgSz.height; i++ ) {
    for (int j = 0; j < imgSz.width; j++ ){
      if ( detail[i][j] <= 0.0){
        detail[i][j] = 0.0;
      }

      ((ushort*)(dvImg->imageData + dvImg->widthStep*i))[j] = (short) detail[i][j];

    }
  }*/

  maxv(detail, max);
  /*double imageMin;
  double imageMax;
  minMaxLoc(tobefiltered, &imageMin, &imageMax);*/
  for (int i = 0; i < imgSz.height; i++ ) {
    for (int j = 0; j < imgSz.width; j++ ){
      if ( detail[i][j] <= 0.0){
        detail[i][j] = 0.0;
      }

      ((ushort*)(dvImg->imageData + dvImg->widthStep*i))[j] =
       (short) (detail[i][j]); // use "max) * 65536.0" instead of max) * maxvue if normalie on single tile.  (short) ( (detail[i][j]/ max) * imageMax);
    }
  }


  /* Size detailSz; // size of output image

  detailSz.width = J*colW;
  detailSz.height = 1*rowW;

  Mat dvImg = Mat(detailSz, CV_16SC1);  //tobefiltered.depth()
  //this has to be done better
  for (int i = 0; i < detailSz.height; i++ ) {
    for (int j = 0; j < detailSz.width; j++ ){
      if ( detail[i][j] <= 0.0){
        detail[i][j] = 0.0;
      }
      ((ushort*)(dvImg.data + i * dvImg.step))[j] = (short) detail[i][j];
    }
  }

  maxv(detail, max);
  for (int i = 0; i < detailSz.height; i++ ) {
    for (int j = 0; j < detailSz.width; j++ ){
      if ( detail[i][j] <= 0.0){
        detail[i][j] = 0.0;
      }
      ((ushort*)(dvImg.data + i * dvImg.step))[j] = (short) ( (detail[i][j]/ max) * 65536.0);
    }
  }
  */


  //imgSWT('/Users/danielfurth/Documents/mouse_atlas/lighsheet_figure/8bit_cropped_x0.25_unspmask3-0.6_s_0607.tif')

  //make adjustments to the original image.
  if(sumbsample==true){
      Size upsampledsize(static_cast<int>(colsWzero), static_cast<int>(rowsWzero));
      resize(tobefiltered, tobefiltered, upsampledsize, INTER_LINEAR);
  }


  Mat muppo = cvarrToMat(dvImg);

  //take out each separate detail coefficient.
  std::vector< int > x1;
  int x2 = 0;
  for (int j=0; j < J; j++) {
    x1.push_back(colW*(j+1));
    Mat subImg = muppo(cv::Range(0, rowW), cv::Range(x2,  x1[j]));

    Mat dest(subImg);
    Mat dest2(subImg);

    //check if subsampling was done due to unequal dimesnions, if so correct for it.
    if(sumbsample==true){
      Size upsampledsize(static_cast<int>(colsWzero), static_cast<int>(rowsWzero));
      resize(dest2, dest2, upsampledsize, INTER_LINEAR);
    }

    if(j==cellB){
      cellImg = dest2.clone();
    }
    //adjust brightness contrast

    /*double minS, maxS;
    cv::minMaxLoc(subImg, &minS, &maxS);

    double alpha = (double)65536/(double)maxS;
    int beta = minS;
    // Do the operation new_image(i,j) = alpha*image(i,j) + beta
    for( int y = 0; y < dest.rows; y++ ){ 
      for( int x = 0; x < dest.cols; x++ ){ 
        for( int c = 0; c < 3; c++ ){
          dest2.at<Vec3b>(y,x)[c] = saturate_cast<uchar>( alpha*( dest.at<Vec3b>(y,x)[c] ) - beta );
        }
      }
    } */

    Mat displayImage(dest2);
    cv::resize(dest2, displayImage, Size(), 0.25, 0.25);

    string String = static_cast<ostringstream*>( &(ostringstream() << j) )->str();

    //imshow( String, displayImage ); // image visualisation
    x2 = x1[j];

    input = "d" + String + "/" + "d" + String + "_" + off + ".tif";

    if(maskOrig){
      double minSrcImg;
      double maxSrcImg;
      minMaxLoc(tobefiltered, &minSrcImg, &maxSrcImg);
      minMaxLoc(dest2, &minVal, &maxVal);
      dest2.convertTo(dest2, CV_32F); //       dest2.convertTo(dest2, CV_32F, 1.00/(maxVal - minVal), -minVal * 1.00/(maxVal - minVal));
      Mat temptobefiltered;
      tobefiltered.convertTo(temptobefiltered, CV_32F);
      dest2 = dest2.mul(temptobefiltered);
      minMaxLoc(dest2, &minVal, &maxVal);
      dest2.convertTo(dest2, CV_16U, maxSrcImg/(maxVal - minVal), -minVal * maxSrcImg/(maxVal - minVal));
    
    }

    try {
      if(!saveroiOnly){
        imwrite(input, dest2);
        Rcpp::Rcout << input << " SAVED" << endl;
      }else{
        if(j==cellB){
          imwrite(input, dest2);
          Rcpp::Rcout << input << " SAVED" << endl;
        }
      }

    }
    catch (runtime_error& ex) {
      Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
      return(R_NilValue);
    }

  }

  if(compEnergy>0){
    //compute optional tensor energy
    Mat gx, gy;
    Scharr(cellImg, gx, CV_32F, 1, 0); //CV_32F
    Scharr(cellImg, gy, CV_32F, 0, 1);


    /// Generate grad_x and grad_y
    Mat grad_x, grad_y;
    Mat abs_grad_x, abs_grad_y;

    /// Gradient X
    //Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
    Scharr( cellImg, grad_x, CV_32F, 1, 0, 1, 0, BORDER_DEFAULT );
    convertScaleAbs( grad_x, abs_grad_x );

    /// Gradient Y
    //Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
    Scharr( cellImg, grad_y, CV_32F, 0, 1, 1, 0, BORDER_DEFAULT );
    convertScaleAbs( grad_y, abs_grad_y );

    /// Total Gradient (approximate)
    Mat grad;
    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
    //resize(grad, grad, Size(), 0.25, 0.25);

    //imshow( "Sobel test", grad );
  

    // Compute the structure tensor, and from it, anisotropy
    //int Sigma = 10;
    Mat gx2, gxy, gy2;
    GaussianBlur(gx.mul(gx), gx2, Size(0, 0), compEnergy); 
    GaussianBlur(gx.mul(gy), gxy, Size(0, 0), compEnergy);
    GaussianBlur(gy.mul(gy), gy2, Size(0, 0), compEnergy);
    MatExpr trace = gx2 + gy2;
    MatExpr det = gx2.mul(gy2) - gxy.mul(gxy);
    Mat second_term;
    sqrt(trace.mul(trace) / 4 - det, second_term);
    MatExpr eig1 = trace / 2 + second_term;
    MatExpr eig2 = trace / 2 - second_term;
    //compute the tensors
    //Mat anisotropy = eig1 / eig2;

    Mat energy = trace;
  minMaxLoc(energy, &minVal, &maxVal);
  Rcpp::Rcout << minVal << " MinVALUE" << endl;
    Rcpp::Rcout << maxVal << " MAXVALUE" << endl;

    if(maxNorm!=0){
      maxVal = (double)maxNorm;
      minVal = (double)minNorm;
    }

    //output mask with original image
    if(maskOrig){
      /*double minSrcImg;
      double maxSrcImg;
      minMaxLoc(tobefiltered, &minSrcImg, &maxSrcImg);
      energy.convertTo(energy, CV_32F, 1.00/(maxVal - minVal), -minVal * 1.00/(maxVal - minVal));
      Mat temptobefiltered;
      tobefiltered.convertTo(temptobefiltered, CV_32F);
      energy = energy.mul(temptobefiltered);
      minMaxLoc(energy, &minVal, &maxVal);
      energy.convertTo(energy, CV_16U, maxSrcImg/(maxVal - minVal), -minVal * maxSrcImg/(maxVal - minVal)); */
      energy.convertTo(energy, CV_8U, 65535.0/(maxVal - minVal), -minVal * 65535.0/(maxVal - minVal));

    }else{
        energy.convertTo(energy, CV_8U, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));

    }
    
    string cellBodyfilename = static_cast<ostringstream*>( &(ostringstream() << cellB) )->str();
    input = "trace/d" + cellBodyfilename + "_" + "Energy_" + off + ".tif";
    try {
      imwrite(input, energy);
      Rcpp::Rcout << input << " SAVED" << endl;
    }
    catch (runtime_error& ex) {
      Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
      return(R_NilValue);
    }
  }
     /*
  END
  */
END_RCPP  
}



