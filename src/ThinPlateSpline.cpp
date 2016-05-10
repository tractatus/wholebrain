// ThinPlateSpline.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

//for sleep
#include <unistd.h>

#include <Rcpp.h>

#include "wholebrain.h"

using namespace cv;
using namespace Rcpp;

//#include "stdafx.h"
//#include <cv.h>
//#include <cxcore.h>
//#include <highgui.h>
#include "CThinPlateSpline.h"

RcppExport SEXP ThinPlateRegistration(SEXP input, SEXP thresh, SEXP verbose){
    BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::CharacterVector fname(input);
    std::string ffname(fname[0]);
    Rcpp::Rcout << "Loading image:" << ffname << std::endl;
	// load a nice picture
    Mat img = imread(ffname, -1);
    imshow("original",img);
	// generate some generic points
	// usually you would use a interest point detector such as SURF or SIFT
	std::vector<cv::Point> iP, iiP;

	// push some points into the vector for the source image
	iP.push_back(cv::Point(50,50));
	iP.push_back(cv::Point(400,50));
	iP.push_back(cv::Point(50,400));
	iP.push_back(cv::Point(400,400));
	iP.push_back(cv::Point(256,256));
	iP.push_back(cv::Point(150,256));

	// push some point into the vector for the dst image
	iiP.push_back(cv::Point(70,70));
	iiP.push_back(cv::Point(430,60));
	iiP.push_back(cv::Point(60,410));
	iiP.push_back(cv::Point(430,420));
	iiP.push_back(cv::Point(220,280));
	iiP.push_back(cv::Point(180,240));

	// create thin plate spline object and put the vectors into the constructor
	CThinPlateSpline tps(iP,iiP);
	
	// warp the image to dst
	Mat dst;
	tps.warpImage(img,dst,0.01,INTER_CUBIC,BACK_WARP);

	// show images
	
	cv::imshow("distorted",dst);
    
	cv::waitKey(0);

    END_RCPP
    
}

