#include <Rcpp.h>
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/imgproc/types_c.h>

using namespace cv;
using namespace std;
using namespace Rcpp;


string ImgTypesConv(int imgTypeInt)
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

void rotation(cv::Mat &img, int rotflag){
  //1=CW, 2=CCW, 3=180
  if (rotflag == 90){
    Mat tmp = img.t();
    flip(tmp, img, 1); //transpose+flip(1)=CW
  } else if (rotflag == -90) {
    Mat tmp = img.t();
    flip(tmp, img, 0); //transpose+flip(0)=CCW     
  } else if (rotflag == 180){
    flip(img, img,-1);    //flip(-1)=180          
  } else if (rotflag != 0){ //if not 0,1,2,3:
    Point2f src_center(img.cols/2.0F, img.rows/2.0F);
    Mat rot_mat = getRotationMatrix2D(src_center, (double)rotflag, 1.0);
    warpAffine(img, img, rot_mat, img.size());
  }
}

/* apply operation to stack */
RcppExport SEXP rgbTogray(SEXP input, SEXP writetoconsole, SEXP saveoutput, SEXP invertimg, SEXP rotate) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);


  Rcpp::CharacterVector outputfilename(saveoutput);
  std::string ffoutputfilename(outputfilename[0]);

  bool verbose = Rcpp::as<bool>(writetoconsole);
  bool invert = Rcpp::as<bool>(invertimg);
  int angle = Rcpp::as<int>(rotate);

  if(verbose){Rcpp::Rcout << "Loading image:" << ffname << std::endl;}
  Mat img = imread(ffname, -1);
  if(verbose){Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;}

  if(rotate!=0){
    if(verbose){Rcpp::Rcout << "Rotating..." << angle << " degrees" << std::endl;}
    rotation(img, angle);
    if(verbose){Rcpp::Rcout << "====== ROTATION DONE ======" << std::endl;}
  }
   
  
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


RcppExport SEXP invertImg(SEXP input, SEXP writetoconsole, SEXP saveoutput) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);


  Rcpp::CharacterVector outputfilename(saveoutput);
  std::string ffoutputfilename(outputfilename[0]);

  bool verbose = Rcpp::as<bool>(writetoconsole);

  if(verbose){Rcpp::Rcout << "Loading image:" << ffname << std::endl;}
  Mat img = imread(ffname, -1);
  if(verbose){Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;}

   
  
  if(img.type()==16){
    cvtColor(img , img , CV_BGR2GRAY);
      bitwise_not ( img, img );

    img.convertTo(img, CV_16S);

    //Rcpp::Rcout << "Image type: " <<  ImgTypes(img.type()) << "_" << img.type()  << std::endl;
    imwrite(ffoutputfilename, img);

  }else{
    bitwise_not ( img, img );
    //img.convertTo(img, CV_16S);
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

RcppExport SEXP morphologyEx(SEXP input, SEXP morphElem, SEXP morphSize, SEXP morphOperator, SEXP saveuchar, SEXP writetoconsole, SEXP saveoutput) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);

  Rcpp::CharacterVector outputfilename(saveoutput);
  std::string ffoutputfilename(outputfilename[0]);

  bool eightbit = Rcpp::as<bool>(saveuchar);
  bool verbose = Rcpp::as<bool>(writetoconsole);



  Mat src, dst;

  int morph_elem = Rcpp::as<int>(morphElem);// 0;
  int morph_size = Rcpp::as<int>(morphSize);// 0;
  int morph_operator = Rcpp::as<int>(morphOperator);// 0;  "Operator:\n 0: Opening - 1: Closing \n 2: Gradient - 3: Top Hat \n 4: Black Hat",
  int const max_operator = 4;
  int const max_elem = 2;
  int const max_kernel_size = 21; //"Kernel size:\n 2n +1"



  if(verbose){Rcpp::Rcout << "Loading image:" << ffname << std::endl;}
   src = imread(ffname, -1);
  if(verbose){Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;}

   if(verbose){
    Rcpp::Rcout << "Image type: " <<  ImgTypesConv(src.type()) << "_" << src.type()  << std::endl;
  }
  //resize(img, img, Size(), resizeParam, resizeParam, INTER_LINEAR);
  //Mat originalImag;
  //img.copyTo(originalImag);

    /// Apply the specified morphology operation
    if(verbose){Rcpp::Rcout << "====== SETTING OPERATOR ======" << std::endl;}

    int operation = morph_operator + 2;

  Mat element = getStructuringElement( MORPH_RECT, Size( 2*morph_size + 1, 2*morph_size+1 ), Point( morph_size, morph_size ) );
    if(verbose){Rcpp::Rcout << "====== APPLYING FILTER ======" << std::endl;}

  morphologyEx( src, dst, operation, element );
    if(verbose){Rcpp::Rcout << "====== FILTER DONE ======" << std::endl;}

    if(eightbit){
      double min, max;
      cv::minMaxLoc(dst, &min, &max);
      dst -= min;
      dst.convertTo(dst,CV_8UC1,255.0/(max-min)); //-255.0/Min
    }


    imwrite(ffoutputfilename, dst);

  return R_NilValue;

END_RCPP  
}


