#include <Rcpp.h>
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;
using namespace Rcpp;


// take number image type number (from cv::Mat.type()), get OpenCV's enum string.
string ImgTyper(int imgTypeInt)
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


/* apply operation to stack */
RcppExport SEXP getAtlasNissl(SEXP maskfile, SEXP ageafile) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector fname(maskfile);
  std::string ffname(fname[0]);
    
  Rcpp::CharacterVector ageaname(ageafile);
  std::string ageafname(ageaname[0]);

  Mat mask = imread(ffname, -1);
  Mat agea = imread(ageafname, -1);

  //resize(img, img, Size(), resizeParam, resizeParam, INTER_LINEAR);
  cvtColor(mask , mask , CV_BGRA2GRAY);
  cvtColor(agea , agea , CV_BGR2GRAY);  
  bitwise_not( agea, agea );
  //GaussianBlur(agea, agea, 15, 2, 2)
  GaussianBlur(agea, agea, Size(15,15), 2);
  double resizeParam = (double)(7900/395);
  Mat ageaResize;
  resize(agea, ageaResize, Size(11300, 7900));
  GaussianBlur(ageaResize, ageaResize, Size(301,301), 40);

  Mat outputMat;
  ageaResize.copyTo(outputMat, mask);
  normalize(outputMat, outputMat, 0, 255, NORM_MINMAX, CV_8UC1);
  imwrite(ffname, outputMat);
 
  return R_NilValue;
   /*
  END
  */
END_RCPP  
}



