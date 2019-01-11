// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/imgproc/types_c.h>

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

/*
template<> SEXP wrap(const cv::Mat &obj) {

}
*/

Mat tobefiltered; // no local
Mat cellImg;
Mat processImg;

//to handle execution time logging
clock_t t0, t1;




/* function to compute the maximum from an image array */
void* maxval(vector<vector<double> > &arr, double &max){
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

//find closest value in a vector
//std::lower_bound: returns the first value that does not compare less
//std::upper_bound: returns the first value that compares strictly greater
int closest(std::vector<int> const& vec, int value) {
 //auto const it = std::upper_bound(vec.begin(), vec.end(), value);
 
 auto const upper = std::lower_bound(vec.begin(), vec.end(), value);
 auto const lower = upper-1;
    
 if (upper == vec.end()) { return -1; }

 return *lower;
}

// take number image type number (from cv::Mat.type()), get OpenCV's enum string.
string getImgType(int imgTypeInt)
{
    int numImgTypes = 35; // 7 base types, with five channel options each (none or C1, ..., C4)

    int enum_ints[] =       {CV_8U,  CV_8UC1,  CV_8UC2,  CV_8UC3,  CV_8UC4,
                             CV_8S,  CV_8SC1,  CV_8SC2,  CV_8SC3,  CV_8SC4,
                             CV_16U, CV_16UC1, CV_16UC2, CV_16UC3, CV_16UC4,
                             CV_16S, CV_16SC1, CV_16SC2, CV_16SC3, CV_16SC4,
                             CV_32S, CV_32SC1, CV_32SC2, CV_32SC3, CV_32SC4,
                             CV_32F, CV_32FC1, CV_32FC2, CV_32FC3, CV_32FC4,
                             CV_64F, CV_64FC1, CV_64FC2, CV_64FC3, CV_64FC4};

    string enum_strings[] = {"CV_8U",  "CV_8UC1",  "CV_8UC2",  "CV_8UC3",  "CV_8UC4",
                             "CV_8S",  "CV_8SC1",  "CV_8SC2",  "CV_8SC3",  "CV_8SC4",
                             "CV_16U", "CV_16UC1", "CV_16UC2", "CV_16UC3", "CV_16UC4",
                             "CV_16S", "CV_16SC1", "CV_16SC2", "CV_16SC3", "CV_16SC4",
                             "CV_32S", "CV_32SC1", "CV_32SC2", "CV_32SC3", "CV_32SC4",
                             "CV_32F", "CV_32FC1", "CV_32FC2", "CV_32FC3", "CV_32FC4",
                             "CV_64F", "CV_64FC1", "CV_64FC2", "CV_64FC3", "CV_64FC4"};

    for(int i=0; i<numImgTypes; i++)
    {
        if(imgTypeInt == enum_ints[i]) return enum_strings[i];
    }
    return "unknown image type";
}

void FindBlobs(const cv::Mat &binary, std::vector < std::vector<cv::Point2i> > &blobs, int lengthlimitation)
{
    blobs.clear();

    // Fill the label_image with the blobs
    // 0  - background
    // 1  - unlabelled foreground
    // 2+ - labelled foreground

    cv::Mat label_image;
    binary.convertTo(label_image, CV_32SC1);

    int label_count = 2; // starts at 2 because 0,1 are used already

    for(int y=0; y < label_image.rows; y++) {
        int *row = (int*)label_image.ptr(y);
        for(int x=0; x < label_image.cols; x++) {
            if(row[x] != 1) {
                continue;
            }

            cv::Rect rect;
            cv::floodFill(label_image, cv::Point(x,y), label_count, &rect, 0, 0, 8); //flag of four only the four pixels (non-diagonal) flag of 4 only 

            std::vector <cv::Point2i> blob;

            for(int i=rect.y; i < (rect.y+rect.height); i++) {
                int *row2 = (int*)label_image.ptr(i);
                for(int j=rect.x; j < (rect.x+rect.width); j++) {
                    if(row2[j] != label_count) {
                        continue;
                    }

                    blob.push_back(cv::Point2i(j,i));
                }
            }
            if(blob.size()>lengthlimitation){
              blobs.push_back(blob);

              label_count++;
            }
        }
    }
}

/**
 * Perform one thinning iteration.
 * Normally you wouldn't call this function directly from your code. image using Zhang-Suen algorithm. http://answers.opencv.org/question/3207/what-is-a-good-thinning-algorithm-for-getting-the-skeleton-of-characters-for-ocr/
 *
 * Parameters:
 *    im    Binary image with range = [0,1]
 *    iter  0=even, 1=odd
 */
void thinningIteration(cv::Mat& img, int iter)
{
    CV_Assert(img.channels() == 1);
    CV_Assert(img.depth() != sizeof(uchar));
    CV_Assert(img.rows > 3 && img.cols > 3);

    cv::Mat marker = cv::Mat::zeros(img.size(), CV_8UC1);

    int nRows = img.rows;
    int nCols = img.cols;

    if (img.isContinuous()) {
        nCols *= nRows;
        nRows = 1;
    }

    int x, y;
    uchar *pAbove;
    uchar *pCurr;
    uchar *pBelow;
    uchar *nw, *no, *ne;    // north (pAbove)
    uchar *we, *me, *ea;
    uchar *sw, *so, *se;    // south (pBelow)

    uchar *pDst;

    // initialize row pointers
    pAbove = NULL;
    pCurr  = img.ptr<uchar>(0);
    pBelow = img.ptr<uchar>(1);

    for (y = 1; y < img.rows-1; ++y) {
        // shift the rows up by one
        pAbove = pCurr;
        pCurr  = pBelow;
        pBelow = img.ptr<uchar>(y+1);

        pDst = marker.ptr<uchar>(y);

        // initialize col pointers
        no = &(pAbove[0]);
        ne = &(pAbove[1]);
        me = &(pCurr[0]);
        ea = &(pCurr[1]);
        so = &(pBelow[0]);
        se = &(pBelow[1]);

        for (x = 1; x < img.cols-1; ++x) {
            // shift col pointers left by one (scan left to right)
            nw = no;
            no = ne;
            ne = &(pAbove[x+1]);
            we = me;
            me = ea;
            ea = &(pCurr[x+1]);
            sw = so;
            so = se;
            se = &(pBelow[x+1]);

            int A  = (*no == 0 && *ne == 1) + (*ne == 0 && *ea == 1) + 
                     (*ea == 0 && *se == 1) + (*se == 0 && *so == 1) + 
                     (*so == 0 && *sw == 1) + (*sw == 0 && *we == 1) +
                     (*we == 0 && *nw == 1) + (*nw == 0 && *no == 1);
            int B  = *no + *ne + *ea + *se + *so + *sw + *we + *nw;
            int m1 = iter == 0 ? (*no * *ea * *so) : (*no * *ea * *we);
            int m2 = iter == 0 ? (*ea * *so * *we) : (*no * *so * *we);

            if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
                pDst[x] = 1;
        }
    }

    img &= ~marker;
}

/**
 * Function for thinning the given binary image
 *
 * Parameters:
 *    src  The source image, binary with range = [0,255]
 *    dst  The destination image
 */
void thinning(const cv::Mat& src, cv::Mat& dst)
{
    dst = src.clone();
    dst /= 255;         // convert to binary image

    cv::Mat prev = cv::Mat::zeros(dst.size(), CV_8UC1);
    cv::Mat diff;

    do {
        thinningIteration(dst, 0);
        thinningIteration(dst, 1);
        cv::absdiff(dst, prev, diff);
        dst.copyTo(prev);
    } 
    while (cv::countNonZero(diff) > 0);

    dst *= 255;
}

void stepwisethreshold(const cv::Mat& src, vector<int> thresholds, bool displayImage, double minArea, double maxArea, vector<vector<Point> >& contours,  vector<Vec4i>& hierarchy,  cv::Mat& dst)
{
  Mat tmp = Mat::zeros(src.size(), src.type());
  //create a lineary spaced array with spacing of eight.
  vector<int> linearspaced(8); 
  int p = thresholds[0];
  int increment = (thresholds[1] - thresholds[0])/8;
  for (unsigned int j=0; j<linearspaced.size(); j++){
    linearspaced[j]=p;
    p+= increment;
    Rcpp::Rcout << j << " = " << p << std::endl;
  }
  
  for (unsigned int i=0; i<linearspaced.size(); i++) {
    int thresh = linearspaced[i];
    threshold(src, tmp, thresh, 255, CV_THRESH_BINARY );
    findContours(tmp, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
    Mat mask = Mat::zeros(src.size(), src.type());

    for (size_t k = 0; k < contours.size(); ++k)
    {
      const vector<Point>& contour = contours[k];
      double area0 = contourArea(contour);
      //if area falls within the area limits given by R via "alim" parameter then draw the contour.
      if( (area0 > minArea) && (area0 < maxArea) ){
        Scalar color(255, 255, 255);
        drawContours( mask, contours, k, color, cv::FILLED, 8, hierarchy );
      }
    }
    //add the contours to the destination Mat.
    bitwise_or(dst, mask, dst);
    
  } 

  findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
  Mat mask = Mat::zeros(src.size(), src.type());
  for (size_t k = 0; k < contours.size(); ++k)
  {
    const vector<Point>& contour = contours[k];
    double area0 = contourArea(contour);
    //if area falls within the area limits given by R via "alim" parameter then draw the contour.
    if( (area0 > minArea) && (area0 < maxArea) ){
        Scalar color(255, 255, 255);
        drawContours( mask, contours, k, color, cv::FILLED, 8, hierarchy );
    }
  }
  dst = mask;



  if(displayImage){
    Mat displayImage;
    src.copyTo(displayImage);
    
    Scalar red(0, 0, 255);
    cvtColor(displayImage, displayImage, CV_GRAY2RGB);
    findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
    for( int i = 0; i< contours.size(); i++ ){
      drawContours( displayImage, contours, i, red, 1.5, 8, hierarchy, 0, Point() );
    }

    double ScaleFactor;
    if(src.rows>600){
      ScaleFactor = 600/(double)src.rows;
      resize(displayImage, displayImage, Size(), ScaleFactor, ScaleFactor, INTER_LINEAR);
    }

    Rcpp::Rcout << "ScaleFactor: " << ScaleFactor << std::endl;
    namedWindow("Image", WINDOW_NORMAL);
    imshow("Image", displayImage);
  }
}

class StepThreshold
{
public:
    Mat src; 
    Mat dst;
    Mat dsp;
    Mat out;
    bool displayImage;
    int alpha = 10;
    double alphaADJ;
    double Max;
    double Min;
    int beta = 0;
    int imgdepth;
    int numThresholds;
    int minArea;
    int maxArea; 
    int minThresh;
    int maxThresh;
    int eccentricityThresh;
    int hideFilter = 0;
    vector<vector<Point> > contours;  
    vector<Vec4i> hierarchy;
    vector<float> eccentricity;
    Point2f centroid;
    vector<float> centroidX;
    vector<float> centroidY;

    unsigned int p;
    unsigned int increment;
    unsigned int maxValue;

    void runthreshold(void)
    {
      centroidX.clear();
      centroidY.clear();
      Mat tmp = Mat::zeros(src.size(), CV_8UC1);
      dst = Mat::zeros(src.size(), CV_8UC1);

      minMaxIdx(src, &Min, &Max);

      if(maxArea<minArea){
        maxArea = minArea+1;
      }

      int linearspaced[numThresholds];
      int p = minThresh;
      double range = ((double)maxThresh - (double)minThresh);
      int increment = range/numThresholds;

      for (unsigned int j=0; j<numThresholds; j++){
        linearspaced[j]=p;
        p+= increment;
      }
  
      for (unsigned int i=0; i<numThresholds; i++) {
        int thresh;
        if(imgdepth==1000){
          thresh = 65535*((double)linearspaced[i]/1000);
        }else{
          thresh = linearspaced[i];
        }
        threshold(src, tmp, thresh, imgdepth, CV_THRESH_BINARY );
        tmp.convertTo(tmp, CV_8UC1);
        findContours(tmp, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
        Mat mask = Mat::zeros(src.size(),CV_8UC1);
        for (size_t k = 0; k < contours.size(); ++k)
        {
          const vector<Point>& contour = contours[k];
          double area0 = contourArea(contour);
          double eccentricity;
          if(contour.size()>5){
          RotatedRect ellipse = fitEllipse( contour ); 
            eccentricity = ellipse.size.height / ellipse.size.width;
          }else{
            eccentricity = 1;
          }
          //if area falls within the area limits given by R via "alim" parameter then draw the contour.
          if( (area0 > minArea) && (area0 < maxArea) && (eccentricity< (3*((double)eccentricityThresh/1000)+1) ) ){
            Scalar color(255, 255, 255);
            drawContours( mask, contours, k, color, cv::FILLED, 8, hierarchy );
      
          }
        }
        //add the contours to the destination Mat.
        bitwise_or(dst, mask, dst);
      } 

      findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
      Mat mask = Mat::zeros(src.size(), CV_8UC1);
      vector<Moments> mu(contours.size());
      for (size_t k = 0; k < contours.size(); ++k)
      {
        const vector<Point>& contour = contours[k];
        double area0 = contourArea(contour);
        
        //if area falls within the area limits given by R via "alim" parameter then draw the contour.
        if( (area0 > minArea) && (area0 < maxArea) ){
          Scalar color(255, 255, 255);
          drawContours( mask, contours, k, color, cv::FILLED, 8, hierarchy, 0  );
           mu[k] = moments( contours[k], false );
           centroidX.push_back( mu[k].m10/mu[k].m00 ); //
           centroidY.push_back( mu[k].m01/mu[k].m00 ); //
        }
      }
      dst = mask.clone();



      if(displayImage){
        dsp = src.clone();
        dsp -= Min;
        dsp.convertTo(dsp,CV_8UC1,255.0/(Max-Min)); //-255.0/Min
        alphaADJ = (double)alpha/10;
        dsp.convertTo(dsp, -1, alphaADJ, beta);
        Scalar red(0, 0, 255);
        cvtColor(dsp, dsp, CV_GRAY2RGB);
        findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
        if(hideFilter){}else{
        for( int i = 0; i< contours.size(); i++ ){
          drawContours( dsp, contours, i, red, cv::FILLED, 8, hierarchy, 0 );
        }
        }
        //REMOVE
        out = dst.clone();
        out.convertTo(out, CV_8UC1);
        Scalar color(255, 255, 255);
        findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
        for( int i = 0; i< contours.size(); i++ ){
          drawContours( out, contours, i, color, cv::FILLED, 8, hierarchy, 0 );
        }
        //REMOVE

        double ScaleFactor;
        if(src.rows>600){
          ScaleFactor = 600/(double)src.rows;
          resize(dsp, dsp, Size(), ScaleFactor, ScaleFactor, INTER_LINEAR);
        }
      }else{
        out = dst.clone();
        out.convertTo(out, CV_8UC1);
        Scalar color(255, 255, 255);
        findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
        for( int i = 0; i< contours.size(); i++ ){
          drawContours( out, contours, i, color, cv::FILLED, 8, hierarchy, 0 );
        }
      }
    }
};


void onTrackbarLT(int v, void *vp)
{
    StepThreshold *pd = static_cast<StepThreshold *>(vp);
    (*pd).runthreshold();
    imshow("display", pd->dsp);
    ////imshow("DST", pd->dst);
    ////imshow("OUT", pd->dst);
    if(pd->minThresh > pd->maxThresh){
      int rangcorr =  pd->maxThresh;
      pd->minThresh = rangcorr;
      setTrackbarPos("upper threshold", "Controls", pd->minThresh);
    }
    //Rprintf("Lower threshold: \v%d \r", pd->thresholds[0]);
    Rcpp::Rcout << " Lower threshold: " << (pd->minThresh) << "\r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
}

void onTrackbarUT(int v, void *vp)
{
    StepThreshold *pd = static_cast<StepThreshold *>(vp);
    (*pd).runthreshold();
    //imshow("DST", pd->dst);
    //imshow("OUT", pd->dst);
    imshow("display", pd->dsp);
    if(pd->maxThresh < pd->minThresh){
      int rangcorr = pd->minThresh;
      pd->maxThresh = rangcorr;
      setTrackbarPos("upper threshold", "Controls", pd->maxThresh);
    }
    //Rprintf("Lower threshold: \v%d \r", pd->thresholds[0]);
    Rcpp::Rcout << "Upper threshold: " << (pd->maxThresh) << "\r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
}

void onTrackbarMinA(int v, void *vp)
{
    StepThreshold *pd = static_cast<StepThreshold *>(vp);
    (*pd).runthreshold();
    imshow("display", pd->dsp);
    //imshow("DST", pd->dst);
    //imshow("OUT", pd->dst);
    //Rprintf("Lower threshold: \v%d \r", pd->thresholds[0]);
    Rcpp::Rcout << "Min. area: " << (pd->minArea) << "\r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
}

void onTrackbarMaxA(int v, void *vp)
{
    StepThreshold *pd = static_cast<StepThreshold *>(vp);
    (*pd).runthreshold();
    imshow("display", pd->dsp);
    //imshow("DST", pd->dst);
    //imshow("OUT", pd->dst);
    //Rprintf("Lower threshold: \v%d \r", pd->thresholds[0]);
    Rcpp::Rcout << "Max. area: " << (pd->maxArea) << "\r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
}

void onTrackbarEcc(int v, void *vp)
{
    StepThreshold *pd = static_cast<StepThreshold *>(vp);
    (*pd).runthreshold();
    imshow("display", pd->dsp);
    //imshow("DST", pd->dst);
    //imshow("OUT", pd->dst);
    //Rprintf("Lower threshold: \v%d \r", pd->thresholds[0]);
    Rcpp::Rcout << "Eccentricity: " << (pd->eccentricityThresh) << "\r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
}


RcppExport SEXP interactiveThreshold(SEXP input, SEXP numthresh){
BEGIN_RCPP
  Rcpp::RNGScope __rngScope;  
  StepThreshold pd;


  Rcpp::CharacterVector f(input);
  std::string ff(f[0]);
  Rcpp::Rcout << "Loading image:" << ff << std::endl;
  t0 = clock();
  pd.src = imread(ff, -1); // -1 tag means "load as is"
  t1 = clock();
  double time_elapsed = (double)((t1 - t0) / CLOCKS_PER_SEC);
  Rcpp::Rcout << "LOADED." << " loading took " <<  time_elapsed << " seconds." << std::endl;
  Rcpp::Rcout << "Image type: " <<  getImgType(pd.src.type()) << "__" << pd.src.type()  << std::endl;
  int depth;
  if(pd.src.type()==0){
    pd.imgdepth = 255;  
    depth = 255;
  }else if(pd.src.type()==2){
    pd.imgdepth = 1000; //
    depth = 65535;
    pd.src.convertTo(pd.src, CV_16S);
    Rcpp::Rcout << "Changed image type to: " <<  getImgType(pd.src.type()) << "__" << pd.src.type() << std::endl;
  }
 
  pd.minThresh = 0;
  pd.maxThresh = pd.imgdepth;
  pd.displayImage = true;
  pd.minArea = 0;
  pd.maxArea = 1000;
  pd.numThresholds = Rcpp::as<int>(numthresh);
  pd.eccentricityThresh = 1000;



  if(pd.src.cols>4096){
    resize(pd.src, pd.src, Size(), 0.25, 0.25, INTER_LINEAR);
  }
    
    // Initialize parameters
    int histSize = 256;    // bin size
    float histRange[] = { 0, 255 };
    const float *ranges[] = { histRange };
    Mat srcFloat;
    pd.src.convertTo(srcFloat,CV_32F,255/(double)depth);
    
    // Calculate histogram
    MatND hist;
    calcHist( &srcFloat, 1, 0, Mat(), hist, 1, &histSize, ranges, true, false );
    
    // Plot the histogram
    int hist_w = 200; int hist_h = 125;
    int bin_w = cvRound( (double) hist_w/histSize );
    
    Mat histImage( hist_h, hist_w, CV_8UC1, Scalar( 0,0,0) );
    normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
    
    for( int i = 1; i < histSize; i++ )
    {
        line( histImage, Point( bin_w*(i-1), hist_h - cvRound(hist.at<float>(i-1)) ) ,
             Point( bin_w*(i), hist_h - cvRound(hist.at<float>(i)) ),
             Scalar( 255, 0, 0), 2, 8, 0  );
    }
    
    namedWindow( "Histogram", 1 );    imshow( "Histogram", histImage );


  pd.runthreshold();
  imshow("display", pd.dsp);
  moveWindow("display", 100, 300);

  namedWindow("Controls", WINDOW_AUTOSIZE);  
 

  namedWindow("View", WINDOW_AUTOSIZE);
  imshow("View", Mat(1, 200, CV_8UC1));
  moveWindow("View", 800, 300);

  imshow("Controls", Mat(1, 800, CV_8UC1));
     

  Rcpp::Rcout << "The finish the command and return to R command prompt select 'Display Image' and press any key." << std::endl;


  waitKey(0);
  vector<float> arealimits;
  vector<int> range;
  vector<float> eccentricity;

  arealimits.push_back(pd.minArea);
  arealimits.push_back(pd.maxArea);
  range.push_back(pd.minThresh);
  range.push_back(pd.maxThresh);

  destroyAllWindows();
   return List::create(
    _["alim"] = arealimits,
    _["threshold.range"] = range,
    _["eccentricity"] = pd.eccentricityThresh,
    _["alpha"] = pd.alphaADJ,
    _["beta"] = pd.beta,
    _["Max"] = pd.Max,
    _["Min"] = pd.Min,
    _["x"] = pd.centroidX,
    _["y"] = pd.centroidY
  );

  //return(R_NilValue);
END_RCPP  

}

void copySourceTile(const cv::Mat& src, cv::Mat& srcTile, cv::Rect &tile)
{
    //top left corner
    auto tl = tile.tl();
    //bottom left corner
    auto br = tile.br();
    
    cv::Point tloffset, broffset;
    
    //Take care of border cases
    if (tile.x < 0)
    {
      
        tloffset.x = -tile.x;
        tile.x = 0;
    }
    
    if (tile.y < 0)
    {
        tloffset.y = -tile.y;
        tile.y = 0;
    }
    
    if (br.x >= src.cols)
    {
      
        broffset.x = br.x - src.cols + 1;
        tile.width -= broffset.x; 
    }
    
    if (br.y >= src.rows)
    {

        broffset.y = br.y - src.rows + 1;
        tile.height -= broffset.y;
    }
    
    // If any of the tile sides exceed source image boundary we must use copyMakeBorder to make proper paddings for this side
    if (tloffset.x > 0 || tloffset.y > 0 || broffset.x > 0 || broffset.y > 0)
    {
        cv::Rect paddedTile(tile.tl(), tile.br());
        assert(paddedTile.x >= 0);
        assert(paddedTile.y >= 0);
        assert(paddedTile.br().x < src.cols);
        assert(paddedTile.br().y < src.rows);
        
        cv::copyMakeBorder(src(paddedTile), srcTile, tloffset.y, broffset.y, tloffset.x, broffset.x, BORDER_CONSTANT);
    }
    else
    {
        // Entire tile (with paddings lies inside image and it's safe to just take a region:
        src(tile).copyTo(srcTile);
    }
}

RcppExport SEXP createTiles(SEXP input, SEXP tilesize, SEXP overlap, SEXP position, SEXP outputfile) {
    BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    //convert to Cpp int
    int tileSize = Rcpp::as<int>(tilesize);
    int overlapPerc = Rcpp::as<int>(overlap);
    
    Rcpp::CharacterVector of(outputfile);
    std::string off(of[0]);

    int pos = Rcpp::as<int>(position);

    
    Rcpp::CharacterVector f(input);
    std::string ff(f[0]);
    Rcpp::Rcout << "Loading image:" << ff << std::endl;
    Mat sourceImage = imread(ff, -1); //CV_LOAD_IMAGE_GRAYSCALE
    Mat finalImage = sourceImage.clone();
    Rcpp::Rcout << "LOADED." << std::endl;
       R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
    int rows;
    int cols;
    int overlapPixels;
    int newHeight;
    int newWidth; 
    if(overlapPerc==-999){
      rows = (int)( ((double)sourceImage.rows / (double)tileSize)+0.5);
      cols = (int)( ((double)sourceImage.cols / (double)tileSize)+0.5);
  
      Rcpp::Rcout << "Rows: " << rows << std::endl;
      Rcpp::Rcout << "Cols: " << cols << std::endl; 

      int overlapRows = abs( (sourceImage.rows - rows*tileSize) )/rows;
      int overlapCols = abs( (sourceImage.cols - cols*tileSize) )/cols;


      if(overlapCols==overlapRows){
        overlapPixels = overlapRows/2;
      }else if(overlapCols>overlapRows){
        //correct the number of rows
        tileSize = tileSize - overlapCols;
        newHeight = (rows)*(tileSize+overlapCols)-overlapCols*(rows-1);
        newWidth = (cols)*(tileSize+overlapCols)-overlapCols*(cols-1);
        overlapPixels = overlapCols/2;

        

      }else{
        //correct the number of cols
        newHeight = (rows)*(tileSize)-overlapRows*(rows-2);
        newWidth = (cols)*(tileSize)-overlapRows*(cols-2);
        overlapPixels = overlapRows/2;

        
      }

    }else{
      overlapPixels = (tileSize*((double)overlapPerc/100))/2;

      rows = (sourceImage.rows / tileSize) + (sourceImage.rows % tileSize ? 1 : 0);
      cols = (sourceImage.cols / tileSize) + (sourceImage.cols % tileSize ? 1 : 0);
    }

    finalImage = Mat::zeros(newHeight, newWidth, sourceImage.type());

    int x0;
    int y0;

    switch (pos) {
        case 1:
            //topleft
             x0 =   0;
             y0 = 0;
            break;
        case 2:
            //topright
             x0 =   0;
             y0 = 0;
            break;
        case 3:
            //bottomright
             x0 =   (sourceImage.cols - finalImage.cols);
             y0 = (sourceImage.rows - finalImage.rows);
            break;
        case 4:
            //bottomleft
             x0 = 0;
             y0 = (sourceImage.rows - finalImage.rows);
            break;
        }


    Mat dst_roi = sourceImage(Rect(x0,y0,finalImage.cols,finalImage.rows));
    dst_roi.copyTo(finalImage); 
        

    tileSize = tileSize - 2 * overlapPixels;


    Rcpp::Rcout << "Padding:" << overlapPixels << std::endl;
    Rcpp::Rcout << "Center tileSize:" << tileSize << std::endl;
    
    double empiricalOverlap = (double)(overlapPixels*2)/(double)(overlapPixels*2+tileSize);
 
    cv::Mat tileInput, tileOutput;

    int barWidth = 70;
    float progress = 0.0;
    
    int k = 0;
    string filepath;
    for (int rowTile = 0; rowTile < rows; rowTile++)
    {
        for (int colTile = 0; colTile < cols; colTile++)
        {
            cv::Rect srcTile(colTile * tileSize + overlapPixels,
                             rowTile * tileSize + overlapPixels,
                             tileSize + 2 * overlapPixels,
                             tileSize + 2 * overlapPixels);
            
            
            copySourceTile(finalImage, tileInput, srcTile);
            
            string String = std::to_string(k);//static_cast<ostringstream*>( &(ostringstream() << k) )->str();
            
            //imshow( String, displayImage ); // image visualisation
            string filepath;
            filepath = "Tiled_" + off + "/" + off + "_Tile" + String + ".tif";
            try {
                imwrite(filepath, tileInput);
                
            }
            catch (runtime_error& ex) {
                Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
                return(R_NilValue);
            }
            k++;
            
            
        }

        Rcpp::Rcout << "  [";
        int pos = barWidth * progress;
        for (int j = 0; j < barWidth; ++j) {
          if (j < pos) Rcpp::Rcout << "=";
            else if (j == pos) Rcpp::Rcout << ">";
          else Rcpp::Rcout << " ";
        }
        Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::flush;
          R_FlushConsole();
          R_ProcessEvents();
          R_CheckUserInterrupt();
        progress += (float)1/(rows-1);
    }
    
    return List::create(
                        _["rows"] = rows,
                        _["cols"] = cols,
                        _["tilesize"] = tileSize,
                        _["padding"] = overlapPixels,
                        _["overlap"] = empiricalOverlap
                        );
    
    END_RCPP
}

RcppExport SEXP getContours(SEXP input, SEXP thresh, SEXP numobj) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    
    //define integers
    
    int threshcpp = Rcpp::as<int>(thresh);
    int numobjcpp = Rcpp::as<int>(numobj);
    
    Rcpp::CharacterVector f(input);
    std::string ff(f[0]);
    Rcpp::Rcout << "Loading image:" << ff << std::endl;
    Mat image = imread(ff, -1); //CV_LOAD_IMAGE_GRAYSCALE
    Rcpp::Rcout << "LOADED." << std::endl;
    
    
    //canny edge detection to also locate the ventricles
    Mat canny_output;
    vector<vector<Point> > contours;
    
    vector<Vec4i> hierarchy;
    
    vector<int> contourX;
    vector<int> contourY;
    vector<int> id;
    
    /// Detect edges using canny
    Canny( image, canny_output, threshcpp, threshcpp*2, 3 );
    /// Find contours
    findContours( canny_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );
    
    vector<float>contourAreaVector;
    for( int i = 0; i< contours.size(); i++ )
    {
        contourArea(contours[i]);
        double area = contourArea(contours[i]);
        contourAreaVector.push_back(area);
    }
    
    vector<int>indexLarge;
    
    auto largestArea = std::max_element(std::begin(contourAreaVector), std::end(contourAreaVector));
    auto indexLargest = std::distance(std::begin(contourAreaVector), largestArea);
    indexLarge.push_back(indexLargest);
    
    
    //get the largest objects base don area. numobj is the numbe rof objects the user asked for. Usually defined by how many
    for(int i=1; i<numobjcpp;i++){
        contourAreaVector[indexLargest] = 0;
    
        largestArea = std::max_element(std::begin(contourAreaVector), std::end(contourAreaVector));
        indexLargest = std::distance(std::begin(contourAreaVector), largestArea);
        indexLarge.push_back(indexLargest);
        
        
    }
    
    for(size_t k1 = 0; k1 < indexLarge.size(); ++k1){
        const vector<Point>& contour = contours[indexLarge[k1]];

    //after checking for largest objects then extract contours
    for (size_t k2 = 0; k2 < contour.size(); ++k2)
    {
        const Point& pt = contour[k2];
        contourX.push_back(pt.x);
        contourY.push_back(pt.y);
        id.push_back((k1+1));
        
    }
    
    }
    
    
    return List::create(
                        _["x"] = contourX,
                        _["y"] = contourY,
                        _["id"] = id
                        );
    
    //return(R_NilValue);
    END_RCPP
}


/* show a image */
RcppExport SEXP filterImage(SEXP input, SEXP noFilter, SEXP filter_minArea, SEXP filter_maxArea, SEXP filter_minThresh, SEXP filter_maxThresh, SEXP filter_eccentricity, SEXP filter_alpha, SEXP filter_beta, SEXP filter_Max, SEXP filter_Min, SEXP scales, SEXP cellBodies, SEXP processes, SEXP family, SEXP sigmaR, SEXP areaSize, SEXP processLength, SEXP outputfile) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //define wavelet scales
  int J = Rcpp::as<int>(scales);
  //define scale with cell bodies
  int cellB = Rcpp::as<int>(cellBodies);

  //define scale with processes
  int process = Rcpp::as<int>(processes);

  std::vector<double> minArea = Rcpp::as< std::vector<double> >(areaSize); 

  int fiberlength = Rcpp::as<int>(processLength);

  //fluorescence threshold to minimize false positives
  //int fluorThresh = Rcpp::as<int>(fluorescenceThreshold);
  //std::vector<int> fluorThresh  = Rcpp::as< std::vector<int> >(fluorescenceThreshold);

  // BEGIN INITILIZATION OF FILTER LIST
  
  int nofilter = Rcpp::as<int>(noFilter); //change to bool
  int filterminArea = Rcpp::as<int>(filter_minArea);
  int filtermaxArea = Rcpp::as<int>(filter_maxArea);
  int filterminThresh = Rcpp::as<int>(filter_minThresh);
  int filtermaxThresh = Rcpp::as<int>(filter_maxThresh);
  int filtereccentricity = Rcpp::as<int>(filter_eccentricity);
  double alpha = Rcpp::as<double>(filter_alpha);
  int beta = Rcpp::as<int>(filter_beta);
  double Max = Rcpp::as<double>(filter_Max);
  double Min = Rcpp::as<double>(filter_Min);
  // END INITILIZATION OF FILTER LIST 

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
    //function closest defined on line 60.
    int newwidth = closest(vec, rowsW);
    Size oldsize = tobefiltered.size();
    Size newsize(static_cast<int>(newwidth), static_cast<int>(newwidth));
    resize(tobefiltered, tobefiltered, newsize, INTER_LINEAR);
    //get width and height of image.
    rowsW = tobefiltered.rows;
    colsW = tobefiltered.cols;
     Rcpp::Rcout << "New size: "<< rowsW << "x" << rowsW << "." << std::endl;
      sumbsample = true;
  }

  //set up the threshold for fluorescence baseline.
  Mat fluorescenceMask = Mat::zeros(tobefiltered.size(), tobefiltered.type());
  //threshold(tobefiltered, fluorescenceMask, fluorThresh, 255, CV_THRESH_BINARY);
  //vector<vector<Point> > contoursStepwise;
  //vector<Vec4i> hierarchyStepwise;

  /*
  stepwisethreshold(tobefiltered, fluorThresh, true, minArea[0], minArea[1], contoursStepwise, hierarchyStepwise, fluorescenceMask);
  */
  // INITIALIZE A THRESHOLD OBJECT
 /* Rcpp::Rcout << "Initialize threshold:" << std::endl;
  int depth;
  if(!nofilter){
    StepThreshold pd;
    pd.src = tobefiltered.clone();
    if(pd.src.type()==0){
      pd.imgdepth = 255;  
      depth = 255;
    }else if(pd.src.type()==2){
      pd.imgdepth = 1000; //
      depth = 65535;
      pd.src.convertTo(pd.src, CV_16S);
    }
 
    pd.minThresh = filterminThresh;
    pd.maxThresh = filtermaxThresh;
    pd.minArea = filterminArea;
    pd.maxArea = filtermaxArea;
    pd.eccentricityThresh = filtereccentricity;
    pd.displayImage = false; //need to turn displayImage off otherwise interactive mode will start.
    pd.numThresholds = 8;
    pd.runthreshold();

    fluorescenceMask = pd.out.clone();
  }
  Rcpp::Rcout << "Threshold done:" << std::endl;
  // END OF A THRESHOLD OBJECT */

  //bitwise_not( fluorescenceMask, fluorescenceMask ); //invert for background subtraction

  // Set up the detector with default parameters.
  //SimpleBlobDetector detector;

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
  maxval(blur,max);

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
    imwrite(input, cvImg);
    Rcpp::Rcout << input << " SAVED" << endl;

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

  maxval(detail, max);
  /*double imageMin;
  double imageMax;
  minMaxLoc(tobefiltered, &imageMin, &imageMax);*/
  for (int i = 0; i < imgSz.height; i++ ) {
    for (int j = 0; j < imgSz.width; j++ ){
      if ( detail[i][j] <= 0.0){
        detail[i][j] = 0.0;
      }

      ((ushort*)(dvImg->imageData + dvImg->widthStep*i))[j] =
       (short) (detail[i][j]); // use "max) * 65536.0" instead of max) * maxvalue if normalie on single tile.  (short) ( (detail[i][j]/ max) * imageMax);
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

  maxval(detail, max);
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
      resize(fluorescenceMask, fluorescenceMask, upsampledsize, INTER_LINEAR);
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

    //check if wavelet scale is the one chosen for segmenting cell bodies
    if(j==cellB){
      cellImg = dest2.clone();
    }
    //check if wavelet scale is the one chosen for segmenting processes
    if(j==process){
      processImg = dest2.clone();
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

    string String = std::to_string(j);//static_cast<ostringstream*>( &(ostringstream() << j) )->str();

    //imshow( String, displayImage ); // image visualisation
    x2 = x1[j];

    input = "d" + String + "/" + "d" + String + "_" + off + ".tif";
    try {
      imwrite(input, dest2);
      Rcpp::Rcout << input << " SAVED" << endl;

    }
    catch (runtime_error& ex) {
      Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
      return(R_NilValue);
    }

  }

  //START SEGMENTATION OF CELL BODIES
  // Compute images gradients
  resize(cellImg, cellImg, Size(), 1, 1);

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
  int Sigma = Rcpp::as<int>(sigmaR);
  Mat gx2, gxy, gy2;
  GaussianBlur(gx.mul(gx), gx2, Size(0, 0), Sigma); 
  GaussianBlur(gx.mul(gy), gxy, Size(0, 0), Sigma);
  GaussianBlur(gy.mul(gy), gy2, Size(0, 0), Sigma);
  MatExpr trace = gx2 + gy2;
  MatExpr det = gx2.mul(gy2) - gxy.mul(gxy);
  Mat second_term;
  sqrt(trace.mul(trace) / 4 - det, second_term);
  MatExpr eig1 = trace / 2 + second_term;
  MatExpr eig2 = trace / 2 - second_term;
  //compute the tensors
  //Mat anisotropy = eig1 / eig2;
  Mat energy = trace;
  Mat coherency = (eig1 - eig2)/(eig2 + eig1);
  Mat orientation = Mat(gx2.rows, gx2.cols, CV_32F);
  for(int i = 0; i < gx2.rows; i++){
    for(int j = 0; j < gx2.cols; j++){
      // Retrieve a single value
      float valueX = gx2.at<float>(i,j);
      float valueY = gy2.at<float>(i,j);
      float valueXY = gxy.at<float>(i,j);
      // Calculate the corresponding single direction, done by applying the arctangens function
      float result = 0.5*atan2(2*valueXY, (valueY-valueX) );
      // Store in orientation matrix element
      orientation.at<float>(i,j) = result;
    }
  }

  //convert Energy, Coherence, and Orientation to 8-bit.

  double maxVal, minVal;

  minMaxLoc(energy, &minVal, &maxVal);
  energy.convertTo(energy, CV_8U, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));

  //minMaxLoc(coherency, &minVal, &maxVal); //find minimum and maximum intensities
  coherency.convertTo(coherency, CV_8U, 255.0/(1 - 0), -0 * 255.0/(1 - 0));
  
  minMaxLoc(orientation, &minVal, &maxVal);
  orientation.convertTo(orientation, CV_8U, 179.0/(maxVal*1.2 - minVal), -minVal * 179.0/(maxVal - minVal));

    
    
    bool autorange = false;
    
    Mat displaypurposeAdjust;
    if(autorange)
    {
        minMaxLoc(tobefiltered, &minVal, &maxVal);
        tobefiltered.convertTo(tobefiltered, -1, alpha, beta); //if 16 bit make adjustment for display purposes.
        minMaxLoc(tobefiltered, &minVal, &maxVal);
        tobefiltered.convertTo(tobefiltered, CV_8U, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));
        displaypurposeAdjust = tobefiltered.clone();
        minMaxLoc(displaypurposeAdjust, &minVal, &maxVal);
        displaypurposeAdjust -= minVal;
        displaypurposeAdjust.convertTo(displaypurposeAdjust,CV_8UC1,255.0/(maxVal-minVal)); //-255.0/Min       255.0/(maxVal-minVal)
        displaypurposeAdjust.convertTo(displaypurposeAdjust, -1, alpha, beta);
    }else{
        tobefiltered -= Min;
        tobefiltered.convertTo(tobefiltered,CV_8UC1,255.0/(Max-Min)); //-255.0/Min
        tobefiltered.convertTo(tobefiltered, -1, alpha, beta);
        
        
        displaypurposeAdjust = tobefiltered.clone();
        
    }
    
    //minMaxLoc(cellImg, &minVal, &maxVal); //find minimum and maximum intensities
    


    cellImg.convertTo(cellImg,CV_8U,255.0/800); //-255.0/Min

  
    //minMaxLoc(cellImg, &minVal, &maxVal); //find minimum and maximum intensities

    //cellImg.convertTo(cellImg, CV_8U, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));
    
    
    //COMPUTE COHERENCY FOR PROCESS
    Mat gxProcess, gyProcess;
    Scharr(processImg, gxProcess, CV_32F, 1, 0); //CV_32F
    Scharr(processImg, gyProcess, CV_32F, 0, 1);
    /// Generate grad_x and grad_y
    Mat grad_xProcess, grad_yProcess;
    Mat abs_grad_xProcess, abs_grad_yProcess;
    
    /// Gradient X
    Scharr( processImg, grad_xProcess, CV_32F, 1, 0, 1, 0, BORDER_DEFAULT );
    convertScaleAbs( grad_xProcess, abs_grad_xProcess );
    /// Gradient Y
    Scharr( processImg, grad_yProcess, CV_32F, 0, 1, 1, 0, BORDER_DEFAULT );
    convertScaleAbs( grad_yProcess, abs_grad_yProcess );
    
    /// Total Gradient (approximate)
    Mat gradProcess;
    addWeighted( abs_grad_xProcess, 0.5, abs_grad_yProcess, 0.5, 0, grad );
    
    // Compute the structure tensor, and from it, anisotropy
    Mat gx2Process, gxyProcess, gy2Process;
    GaussianBlur(gxProcess.mul(gxProcess), gx2Process, Size(0, 0), Sigma);
    GaussianBlur(gxProcess.mul(gyProcess), gxyProcess, Size(0, 0), Sigma);
    GaussianBlur(gyProcess.mul(gyProcess), gy2Process, Size(0, 0), Sigma);
    MatExpr traceProcess = gx2Process + gy2Process;
    MatExpr detProcess = gx2Process.mul(gy2Process) - gxyProcess.mul(gxyProcess);
    Mat second_termProcess;
    sqrt(traceProcess.mul(traceProcess) / 4 - detProcess, second_termProcess);
    MatExpr eig1Process = traceProcess / 2 + second_termProcess;
    MatExpr eig2Process = traceProcess / 2 - second_termProcess;
    //compute the tensors
    Mat coherencyProcess = (eig1Process - eig2Process)/(eig2Process + eig1Process);
    coherencyProcess.convertTo(coherencyProcess, CV_8U, 255.0/(1 - 0), -0 * 255.0/(1 - 0));
    /*
  resize(energy, energy, Size(), 0.25, 0.25);
  resize(coherency, coherency, Size(), 0.25, 0.25);
  resize(orientation, orientation, Size(), 0.25, 0.25);
  resize(cellImg, cellImg, Size(), 0.25, 0.25);
  */
    
    string cellBodyfilename = std::to_string(cellB);//static_cast<ostringstream*>( &(ostringstream() << cellB) )->str();
    input = "scharr/d" + cellBodyfilename + "Scharr_"   + off + ".tif";
    try {
        imwrite(input, grad);
        Rcpp::Rcout << input << " SAVED" << endl;
        
    }
    catch (runtime_error& ex) {
        Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
        return(R_NilValue);
    }

  input = "trace/d" + cellBodyfilename + "_" + "Energy_" + off + ".tif";
  try {
    imwrite(input, energy);
    Rcpp::Rcout << input << " SAVED" << endl;

  }
  catch (runtime_error& ex) {
    Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
    return(R_NilValue);
  }

  input = "coherency/d" + cellBodyfilename + "_" + "Coherency_"  + off + ".tif";
  try {
    imwrite(input, coherency);
    Rcpp::Rcout << input << " SAVED" << endl;

  }
  catch (runtime_error& ex) {
    Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
    return(R_NilValue);
  }

    input = "orientation/d" + cellBodyfilename + "_" + "Orientation_" + off + ".tif";
  try {
    imwrite(input, orientation);
    Rcpp::Rcout << input << " SAVED" << endl;

  }
  catch (runtime_error& ex) {
    Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
    return(R_NilValue);
  }



 
  

  vector<Mat> array_to_merge;
  vector<Mat> hue;
  vector<Mat> outputimage;
  coherency.convertTo(coherency, CV_8UC1);
  coherencyProcess.convertTo(coherencyProcess, CV_8UC1);
  tobefiltered.convertTo(tobefiltered, CV_8UC1);
  orientation.convertTo(orientation, CV_8UC1);
  cellImg.convertTo(cellImg, CV_8UC1);

  Mat g = Mat::zeros(Size(cellImg.cols, cellImg.rows), CV_8UC1);
  Mat white;
  bitwise_not( g, white);

  array_to_merge.push_back(orientation); //H
  array_to_merge.push_back(cellImg); //L
  array_to_merge.push_back(coherency); //S

  hue.push_back(orientation); //H
  hue.push_back(white); //L
  hue.push_back(white); //S

  outputimage.push_back(orientation); //H
  outputimage.push_back(displaypurposeAdjust); //L
  outputimage.push_back(coherencyProcess); //S

  //imgSWT('/Users/danielfurth/Documents/mouse_atlas/lighsheet_figure/analysis/Resliceofanalysis.tif', cell.bodies=4, fluorescence.threshold=0, areaSize=0)->mupp2

  
  Mat hslImg;

  merge(array_to_merge, hslImg);

  Mat hueMerge;

  merge(hue, hueMerge);

  Mat hslOutput;

  merge(outputimage, hslOutput);

  imwrite("cellImg.tif", cellImg);
    
    input = "rangecorrected/d" + cellBodyfilename + "_" + "rangecorrected_" + off + ".tif";
    try {
        imwrite(input, tobefiltered);
        Rcpp::Rcout << input << " SAVED" << endl;
        
    }
    catch (runtime_error& ex) {
        Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
        return(R_NilValue);
    }

  /* for(int i=0; i<hslImg.rows; i++){
   for(int j=0; j<hslImg.cols; j++){
           hslImg.data[hslImg.step[0]*i + hslImg.step[1]* j + 0] = coherency.data[coherency.step[0]*i + orientation.step[1]* j];
           hslImg.data[hslImg.step[0]*i + hslImg.step[1]* j + 1] = cellImg.data[cellImg.step[0]*i + tobefiltered.step[1]* j];
           hslImg.data[hslImg.step[0]*i + hslImg.step[1]* j + 2] = orientation.data[orientation.step[0]*i + coherency.step[1]* j];
    }
  } */
  Mat final;
  Mat finaloutput;
  cvtColor(hslImg, final, CV_HLS2BGR);
  cvtColor(hslOutput, finaloutput, CV_HLS2BGR);
  cvtColor(tobefiltered, tobefiltered, CV_GRAY2BGR);
  double alphaBlended = 0.7; 
  double betaBlended;
  betaBlended = ( 1.0 - alphaBlended );
  Mat blended;
  addWeighted( final, alphaBlended, tobefiltered, betaBlended, 0.0, blended);

  //only after blending we actually do the segmentation of cells in the Saturation Channel
  Mat cellsegmentation;
  cvtColor(blended, cellsegmentation, CV_BGR2HSV); //has to be HSV and not HLS
  vector<Mat> spl;
  split(cellsegmentation,spl);

  double thres_val = threshold(spl[2], spl[2], 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
  threshold(spl[1], spl[1], thres_val, 255, CV_THRESH_BINARY );
  Rcpp::Rcout << thres_val << " THRESHOLD" << endl;
  //imshow("spl1 H",spl[0]);//h
  //imshow("spl2 L",spl[1]);//l
  //imshow("spl3 S",spl[2]);//s
  Mat segmented;
  //segmented = spl[2].mul(spl[1]);
    bool blendedSegmentation = true;
    if(blendedSegmentation){
  segmented = spl[2] - spl[1];
  
  //bitwise_or(segmented, fluorescenceMask, segmented);
    }
    else{
        segmented=fluorescenceMask;
    }
  //imshow("Fluoroescence Mask", fluorescenceMask);
  //imshow("Segmented",segmented);
  //imshow("DNAIEL", segmented);

  //BLUR A SECOND TIME TO AVOID SMALL fluorescent DEBRIS.
  GaussianBlur(segmented, segmented, Size(0, 0), 2); //set parameter for Sigma2
  threshold(segmented, segmented, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);

  Mat binarycellbodies;
  segmented.copyTo(binarycellbodies); //to be subtracted from Mat final
  //cvtColor(binarycellbodies, binarycellbodies,CV_GRAY2BGR);
  //cvtColor(binarycellbodies, binarycellbodies,CV_BGR2HSV);
  Mat watershedMat;
  segmented.convertTo(watershedMat, CV_8U); 

  distanceTransform(segmented, segmented, CV_DIST_L2, 3);
  normalize(segmented, segmented, 0, 1., cv::NORM_MINMAX);
  threshold(segmented, segmented, .35, 1., CV_THRESH_BINARY); //had to have 0.3

  // Create the CV_8U version of the distance image
  // It is needed for cv::findContours()
  Mat dist_8u;
  segmented.convertTo(dist_8u, CV_8U);

  // Find total markers
  vector<vector<Point> > contours;
  findContours(dist_8u, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

  // Total numbe rof cell bodies
  int ncomp = contours.size();

  cv::Mat markers = Mat::zeros(segmented.size(), CV_32SC1);
  for (int i = 0; i < ncomp; i++)
    drawContours(markers, contours, i, Scalar::all(i+1), -1);

  circle(markers, Point(5,5), 3, CV_RGB(255,255,255), -1); //background marker

  Mat finalseg;
  cvtColor(watershedMat, finalseg,CV_GRAY2BGR,3);
  watershed(finalseg, markers);

  markers.convertTo(markers, CV_8U);
  threshold(markers, markers, 254, 255, 4);

  double min, ncompWatershed;
  minMaxIdx(markers, &min, &ncompWatershed);
  Rcpp::Rcout << "Cell bodies found: " << ncompWatershed << endl;

  Mat matEdge = markers==1;
  Mat matContour(matEdge.size(), CV_8UC1); 

  vector<int> contourPointsX;
  vector<int> contourPointsY;
  vector<int> contourID;
  vector<float> contourSomaArea;
  vector<float> eccentricity;

    //get properties of the cell bodies
  vector<float> somaX;
  vector<float> somaY;

  for (int m = 1; m < ((int)ncompWatershed+1); m++){
  matEdge = markers==m;
  findContours(matEdge, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
      /// Get the moments
  vector<Moments> mu(contours.size() );
    ///  Get the mass centers:
  vector<Point2f> mc( contours.size() );

  // ---- Code below is only used for visualizing the result. ----
  for (size_t k = 0; k < contours.size(); ++k)
  {
      const vector<Point>& contour = contours[k];
      double area0 = contourArea(contour);
      contourSomaArea.push_back(area0);
      if(contour.size()>=5){
        RotatedRect ellipse = fitEllipse( contour ); 
        eccentricity.push_back( ellipse.size.height / ellipse.size.width );
      }else{
        eccentricity.push_back( 0 );
      }
      mu[k] = moments( contours[k], false );
      mc[k] = Point2f( mu[k].m10/mu[k].m00 , mu[k].m01/mu[k].m00 ); 
      somaX.push_back(mc[k].x);
      somaY.push_back(mc[k].y); 
      for (size_t k2 = 0; k2 < contour.size(); ++k2)
      {
          const Point& pt = contour[k2];
          //cout << pt.x << ", " << pt.y << ", "<< endl;
          contourPointsX.push_back(pt.x);
          contourPointsY.push_back(pt.y);
          contourID.push_back(m);

          matContour.at<uint8_t>(pt) = 255;
      }
  }
  }



  //compute intensity
  Mat labels = Mat::zeros(tobefiltered.size(), CV_8UC1);     
  vector<float> contourAvgIntensity(contours.size(), 0.f); // This contains the averages of each contour

  for (size_t i = 0; i < contours.size(); ++i)
  {
    drawContours(labels, contours, i, Scalar(i), cv::FILLED);
    Rect roi = boundingRect(contours[i]);
    Scalar Avg = cv::mean( tobefiltered(roi), labels(roi) == i);
    contourAvgIntensity[i] = Avg[0];
  }

  //imshow("contour", matContour);

  //imshow("final", final);
  //imshow("blended", blended);
  input = "processes_anisotropic/d" + cellBodyfilename + "_anisotropic_fluorescence"  + "_" + off + ".tif";
  try {
    imwrite(input, final);
    input = "output/d" + cellBodyfilename + "_output"  + "_" + off + ".png";
    imwrite(input, finaloutput);
    input = "blended/d" + cellBodyfilename + "blended"  + "_" + off + ".tif";
    imwrite(input, blended);
    Rcpp::Rcout << input << " SAVED" << endl;

  }
  catch (runtime_error& ex) {
    Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
    return(R_NilValue);
  }

  // Gradient magnitude
  Mat magnitude, abs_gx, abs_gy;
  //sqrt(gx.mul(gx) + gy.mul(gy), magnitude);
  convertScaleAbs( gx, abs_gx );
  convertScaleAbs( gy, abs_gy );
  addWeighted( abs_gx, 0.5, abs_gy, 0.5, 0, magnitude );  
  //imshow( "Magnitude", magnitude );

  Mat processes;
  final.copyTo(processes);
  cvtColor(processes, processes, CV_BGR2HSV);
  vector<Mat> processesSplit;
  split(processes, processesSplit);
  cvtColor(processes, processes, CV_HSV2BGR);
  //imshow("Processes",processes);

  thres_val = threshold(processesSplit[2], processesSplit[2], 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
  threshold(processesSplit[1], processesSplit[1], thres_val, 255, CV_THRESH_BINARY );
  Rcpp::Rcout << thres_val << " THRESHOLD" << endl;
  Mat segmentedProcesses;
  segmentedProcesses = processesSplit[1] -binarycellbodies ; //- processesSplit[1]

  Mat processesFiltered;
  medianBlur(segmentedProcesses, processesFiltered, 5);  //change this to 5
  //threshold(segmentedProcesses, segmentedProcesses, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
  //imshow("Processes filtered", processesFiltered);
  Mat skel;
  threshold(processesFiltered, skel, 1, 255, cv::THRESH_BINARY);
  thinning(skel, skel);


  /* //perform skeletonization
  Mat skel(processesFiltered.size(), CV_8UC1, cv::Scalar(0));
  Mat temp;
  Mat eroded;
 
  Mat element = cv::getStructuringElement(MORPH_CROSS, cv::Size(3, 3));
 
  bool done;    
  do
  {
    cv::erode(processesFiltered, eroded, element);
    cv::dilate(eroded, temp, element); // temp = open(img)
    cv::subtract(processesFiltered, temp, temp);
    cv::bitwise_or(skel, temp, skel);
    eroded.copyTo(processesFiltered);
 
    done = (cv::countNonZero(processesFiltered) == 0);
  }while (!done);
  //imshow("Skeleton", skel);

  morphologyEx( skel, skel, 3, element ); */

  Mat labelledSkel = cv::Mat::zeros(skel.size(), CV_8UC3);
  Mat binary;
  std::vector < std::vector<cv::Point2i > > blobs;

  skel.copyTo(binary);
  threshold(binary, binary, 0.0, 1.0, cv::THRESH_BINARY);
  FindBlobs(binary, blobs, fiberlength);
  vector<int> processesPointsX;
  vector<int> processesPointsY;
  vector<int> processesID;
  vector<int> processesRED;
  vector<int> processesGREEN;
  vector<int> processesBLUE;
  vector<int> processesTheta;
  vector<int> processesCoherence;
  vector<int> processType;
  uchar pix;  // To store a pixel intensity for endpoint determination

  Rcpp::Rcout << blobs.size() << " BLOBS" << endl;
  // Randomy color the blobs
  cvtColor(processes, processes, CV_BGR2HSV);
  cvtColor(hueMerge, hueMerge, CV_HSV2BGR);
  for(size_t i=0; i < blobs.size(); i++) {
      unsigned char r = 255 * (rand()/(1.0 + RAND_MAX));
      unsigned char g = 255 * (rand()/(1.0 + RAND_MAX));
      unsigned char b = 255 * (rand()/(1.0 + RAND_MAX));

      for(size_t j=0; j < blobs[i].size(); j++) {
          int x = blobs[i][j].x;
          int y = blobs[i][j].y;

          //check what type of pixel in the skeleton
          // Reset counter that then determnes type
          int count = 0;     

          // For each pixel in the neighbourhood
          // centered at this skeleton location...
          for (int r = -1; r <= 1; r++) {
            for (int c = -1; c <= 1; c++) {

                // Get the pixel in the neighbourhood, but first make sure your within the image boundaries.
                if( ((x+c)>=0)&&((y+r)>=0)&&((c+x)<binary.cols)&&((r+y)<binary.rows) ){
                  pix = binary.at<uchar>(r+y,c+x);
                }else{
                  pix=0;
                }
                // Count if non-zero
                if (pix != 0)
                  count++;
            }
          }
          //Check for vertices
          /*
            
          REPLACE ALL THIS "else if"-SHIT WITH A DATA STRUCTURE LOOK UP

          */
          //if count is equal to 3 it means the pixel is part of an edge and hence no further processing
          if(count==3){
            processType.push_back(1); //its an edge
          }else if(count==2){
            processType.push_back(2); //its a end point vertex
          }else{
            int north, northeast, west, east, northwest, south, southeast, southwest;

            //if its not an endpoint and not an edge then it has to be a joint. Find out which type.
            if( ((x+1) > binary.cols) ){
              east = 0;
              northeast = 0;
              southeast = 0;
            }else{
              east = binary.at<uchar>(y,x+1);
            }

            if( ((x-1) < 0) ){
              west = 0;
              northwest = 0;
              southwest = 0;
            }else{
              west = binary.at<uchar>(y,x-1);
            }

            if( ((y-1) < 0)  ){
              north = 0;
              northeast = 0;
              northwest = 0;
            }else{
              north = binary.at<uchar>(y-1,x);
              if( ((x+1) <= binary.cols) ){
                northeast = binary.at<uchar>(y-1,x+1);
              }
              if( ((x-1) >= 0) ){
                northwest = binary.at<uchar>(y-1,x-1);
              }
            }

            if( ((y+1) > binary.rows) ){
              south = 0;
              southeast = 0;
              southwest = 0;
            }else{
              south = binary.at<uchar>(y+1,x);
              if( ((x+1) <= binary.cols) ){
                southeast = binary.at<uchar>(y+1,x+1);
              }
              if( ((x-1) >= 0) ){
                southwest = binary.at<uchar>(y+1,x-1);
              }
            }

            //define which type it is.
            if( (count==5)&&( (northwest&&northeast&&southwest&&southeast)||(west&&north&&east&&south) ) ){
              processType.push_back(4); //its a vertex with 4 joints
            }else if( (count==4) && ( 
                                      (northwest&&northeast&&southeast)
                                      ||
                                      (northwest&&southeast&&southwest)
                                      ||
                                      (west&&north&&east)
                                      ||
                                      (west&&south&&east)
                                      ||
                                      (west&&north&&south)
                                      ||
                                      (north&&east&&south)
                                      ||
                                      (north&&east&&southwest)
                                      ||
                                      (northwest&&east&&south)
                                      ||
                                      (west&&northeast&&south)
                                      ||
                                      (west&&north&&southeast)
                                    )
                    ){
              processType.push_back(3); //its a vertex with 3 joints
            }else{
              processType.push_back(5); //undefined type
            }
          }


          /* switch(count){
            case 2 :
              processType.push_back(2); //its a end point
              break;
            case 3:
              processType.push_back(1); //its a slab
              break;
            case 4:
              processType.push_back(3); //its a joint with 3 joints
              break;
            case 5:
              processType.push_back(4); //its a joint with 4 joints
              break;    
          } */


          processesPointsX.push_back(x);
          processesPointsY.push_back(y);
          processesID.push_back(i);
          processesBLUE.push_back(hueMerge.at<Vec3b>(y,x)[0]);
          processesGREEN.push_back(hueMerge.at<Vec3b>(y,x)[1]);
          processesRED.push_back(hueMerge.at<Vec3b>(y,x)[2]);

          processesTheta.push_back(processes.at<Vec3b>(y,x)[0]);
          
          processesCoherence.push_back(coherencyProcess.at<char>(y,x));

          labelledSkel.at<Vec3b>(y,x)[0] = b;
          labelledSkel.at<Vec3b>(y,x)[1] = g;
          labelledSkel.at<Vec3b>(y,x)[2] = r;
      }
  }

  //imshow("labelled", labelledSkel);
   

  //save segmented cell bodies
  input = "cellbodies/d" + cellBodyfilename + "segmented"   + off + ".tif";
  try {
    imwrite(input, markers);
    input = "processes_labeled/d" + cellBodyfilename + "labeled"   + off + ".tif";
    imwrite(input, labelledSkel);
    Rcpp::Rcout << input << " SAVED" << endl;

  }
  catch (runtime_error& ex) {
    Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
    return(R_NilValue);
  }


  return List::create(
    _["area"] = contourSomaArea,
    _["centroid.x"] = somaX,
    _["centroid.y"] = somaY,
    _["intensity"] = contourAvgIntensity,
    _["eccentricity"] = eccentricity,
    _["x"] = contourPointsX, 
    _["y"] = contourPointsY,
    _["id"] = contourID,
    _["process.x"] = processesPointsX,
    _["process.y"] = processesPointsY,
    _["process.id"] = processesID,
    _["process.r"] = processesRED,
    _["process.g"] = processesGREEN,
    _["process.b"] = processesBLUE,
    _["process.theta"] = processesTheta,
    _["process.coherence"] = processesCoherence,
    _["process.type"] = processType
  );

  //return(R_NilValue);
END_RCPP  
}
