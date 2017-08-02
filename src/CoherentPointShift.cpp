
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

//for sleep
#include <unistd.h>
//Rcpp
#include <Rcpp.h>

using namespace cv;
using namespace Rcpp;

#include <cpd/utils.hpp>
#include <cpd/matrix.hpp>
//#include <Eigen/Dense>

RcppExport SEXP CoherentPointDriftRegistration(){
    BEGIN_RCPP
    	// SEXP input, SEXP srcX, SEXP srcY, SEXP dstX, SEXP dstY, SEXP resizeP, SEXP MaxDisp, SEXP MinDisp, SEXP outputfile

   Rcpp::RNGScope __rngScope;
   /*
    Rcpp::CharacterVector fname(input);
    std::string ffname(fname[0]);
    Rcpp::Rcout << "Loading image:" << ffname << std::endl;
    Mat img = imread(ffname, -1);
    
    //original image size
    int width = img.cols;
    int height = img.rows;

    double resizeParam = Rcpp::as<double>(resizeP);

    double Max = Rcpp::as<int>(MaxDisp);
    double Min = Rcpp::as<int>(MinDisp);

  	resize(img, img, Size(), resizeParam, resizeParam, INTER_LINEAR);


  	Rcpp::CharacterVector of(outputfile);
    std::string off(of[0]);
    */

    cpd::Matrix fixed(2,2);
    cpd::Matrix moving(2,2);
    
    fixed(0,0) = 3;
  	fixed(1,0) = 2.5;
  	fixed(0,1) = -1;
  	fixed(1,1) = 2;
  	moving = fixed;
  	moving(1,1) = fixed(1,0) + fixed(0,1);

  	Rcpp::Rcout << "Moving:\n" << moving << std::endl;
  	Rcpp::Rcout << "Fixed:\n" << fixed << std::endl;

  	

  	cpd::Matrix fixedFish = cpd::matrix_from_path("/Users/danielfurth/Downloads/cpd-master/tests/fixtures/fish_distorted.csv");
  	cpd::Matrix movingFish = cpd::matrix_from_path("/Users/danielfurth/Downloads/cpd-master/tests/fixtures/fish.csv");


  	Rcpp::Rcout << "FISH:\n" << movingFish << std::endl;

	//cpd::RigidResult result = cpd::rigid(fixed, moving);
	//cpd::RigidResult result = cpd::rigid(fixed, moving);

    /*cpd::Matrix fixed = load_points_from_somewhere();
    cpd::Matrix moving = load_points_from_somewhere();

    
    registration.correspondence(true).outliers(0.2);
    cpd::RegiResult result = registration.run(fixed, moving);
    */

  	return R_NilValue;

    END_RCPP
    
}