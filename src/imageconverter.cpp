#include <Rcpp.h>
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;
using namespace Rcpp;


/*
string ImgTypes(int imgTypeInt)
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
*/


/* apply operation to stack */
RcppExport SEXP rgbTogray(SEXP input, SEXP writetoconsole, SEXP saveoutput, SEXP invertimg) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);


  Rcpp::CharacterVector outputfilename(saveoutput);
  std::string ffoutputfilename(outputfilename[0]);

  bool verbose = Rcpp::as<bool>(writetoconsole);
  bool invert = Rcpp::as<bool>(invertimg);

  if(verbose){Rcpp::Rcout << "Loading image:" << ffname << std::endl;}
  Mat img = imread(ffname, -1);
  if(verbose){Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;}

   
  
  if(img.type()==16){
    cvtColor(img , img , CV_BGR2GRAY);
    if(invert){
      bitwise_not ( img, img );
    }
    img.convertTo(img, CV_16S);

    //Rcpp::Rcout << "Image type: " <<  ImgTypes(img.type()) << "_" << img.type()  << std::endl;
    imwrite(ffoutputfilename, img);

  }

  /*
  return List::create(
    _["x"] = xPoint,
    _["y"] = yPoint,
    _["contour.ID"] = contourID
  );
  */
return R_NilValue;

   /*
  END
  */
END_RCPP  
}



