// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"

//for sleep
#include <unistd.h>

#include <Rcpp.h>

#include "wholebrain.h"

using namespace cv;

/*
template<> SEXP wrap(const cv::Mat &obj) {

}
*/

Mat img; // no local

/* show a image */
RcppExport SEXP loadImage(SEXP filename) {
  Rcpp::CharacterVector f(filename);
  std::string ff(f[0]);
  img = imread(ff);

    namedWindow("image", cv::WINDOW_AUTOSIZE);
  imshow("image", img);
  waitKey(10000);
  destroyWindow("image");

  Rcpp::XPtr<cv::Mat> p(&img, false);
  return(p);
}


RcppExport SEXP showImage(SEXP filename) {
  Rcpp::XPtr<cv::Mat> p(filename);

  namedWindow("image", cv::WINDOW_AUTOSIZE);
  imshow("image", *p);
  waitKey(10000);
  destroyWindow("image");
  return(R_NilValue);
}
