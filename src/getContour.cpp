#include <Rcpp.h>
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;
using namespace Rcpp;
RNG rng(12345);


// take number image type number (from cv::Mat.type()), get OpenCV's enum string.
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


/* apply operation to stack */
RcppExport SEXP getContour(SEXP input, SEXP thresho,  SEXP getLarge, SEXP getNested, SEXP display, SEXP resizeP, SEXP blurImg, SEXP writetoconsole, SEXP saveoutput, SEXP writeout) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);

  int saveout = Rcpp::as<int>(writeout);

  Rcpp::CharacterVector outputfilename(saveoutput);
  std::string ffoutputfilename(outputfilename[0]);

  int thresh = Rcpp::as<int>(thresho);
  int DISPLAY = Rcpp::as<int>(display);
  int numNested = Rcpp::as<int>(getNested);

  bool verbose = Rcpp::as<bool>(writetoconsole);

  double resizeParam = Rcpp::as<int>(resizeP);
  resizeParam = resizeParam/100;
  int blurImage = Rcpp::as<int>(blurImg);



  if(verbose){Rcpp::Rcout << "Loading image:" << ffname << std::endl;}
  Mat img = imread(ffname, -1);
  if(verbose){Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;}

   if(verbose){
    Rcpp::Rcout << "Resizing to: " <<  resizeParam*100 << "% of original size." << std::endl;
    Rcpp::Rcout << "Image type: " <<  ImgTypes(img.type()) << "_" << img.type()  << std::endl;
  }
  resize(img, img, Size(), resizeParam, resizeParam, INTER_LINEAR);

  if(DISPLAY){
    imshow("Original", img);
  }
  if(img.type()==24){
    cvtColor(img , img , CV_BGRA2GRAY);
    Rcpp::Rcout << "Image type: " <<  ImgTypes(img.type()) << "_" << img.type()  << std::endl;
  }
  int depth;
  if(img.type()==0){
    depth = 255;
  }else if(img.type()==2){
    depth = 65535;
    img.convertTo(img, CV_16S);
  }

   if(saveout){
     imwrite(ffoutputfilename, img);
   }

  Mat img_bw;

  if(blurImage!=0){
    blur( img, img, Size(blurImage,blurImage) );
  }

  if(thresh==0){
    threshold(img, img_bw, 0, depth, CV_THRESH_BINARY | CV_THRESH_OTSU);
  }else{
    threshold(img, img_bw, thresh, depth, CV_THRESH_BINARY);
    img.convertTo(img, CV_16S);
  }

  if(DISPLAY){
    imshow("Binary segmentation", img_bw);
  }


  vector<vector<Point> > contours;
  vector<Vec4i> hierarchy;


  img_bw.convertTo(img_bw, CV_8U);

  //cv::cvtColor(colorMat, greyMat, cv::COLOR_BGR2GRAY);

  /// Find contours
  findContours( img_bw, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );

  /// Get the moments
  vector<Moments> mu(contours.size() );
  for( int i = 0; i < contours.size(); i++ )
     { mu[i] = moments( contours[i], false ); }

  ///  Get the mass centers:
  vector<Point2f> mc( contours.size() );
  for( int i = 0; i < contours.size(); i++ )
     { mc[i] = Point2f( mu[i].m10/mu[i].m00 , mu[i].m01/mu[i].m00 ); }

  /// Draw contours
  Mat drawing = Mat::zeros( img.size(), CV_8UC3 );
  for( int i = 0; i< contours.size(); i++ )
     {
       Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
       drawContours( drawing, contours, i, color, 2, 8, hierarchy, 0, Point() );
       circle( drawing, mc[i], 4, color, -1, 8, 0 );
     }

  /// Show in a window
  if(DISPLAY){
    namedWindow( "Contours", CV_WINDOW_AUTOSIZE );
    imshow( "Contours", drawing );
  }
  int largest_area=0;
  int largest_contour_index=0;
 Rect bounding_rect;
 Point largest_centroid;

  /// Calculate the area with the moments 00 and compare with the result of the OpenCV function
  printf("\t Info: Area and Contour Length \n");
  for( int i = 0; i< contours.size(); i++ )
     {

      double a=contourArea( contours[i],false);  //  Find the area of contour
       if(a>largest_area){
       largest_area=a;
       largest_contour_index=i;                //Store the index of largest contour
       bounding_rect=boundingRect(contours[i]); // Find the bounding rectangle for biggest contour
       }


       printf(" * Contour[%d] - Area: %.2f - Length: %.2f  - Parent: %d   - Nested: %d  \n", i, a, arcLength( contours[i], true ), hierarchy[i][3], hierarchy[i][2]);
       Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
       drawContours( drawing, contours, i, color, 2, 8, hierarchy, 0, Point() );
       circle( drawing, mc[i], 4, color, -1, 8, 0 );
     }
     largest_centroid = mc[largest_contour_index];
      printf("Largest contour: %d, # contours detected: %d \n", largest_contour_index, contours.size());

   std::vector<int> ventricles_indexes; 
   std::vector<double> ventricles_size;

      std::vector<int> xPoint;
   std::vector<int> yPoint;
   std::vector<int> contourID;
   vector<int>::iterator it;
   int counter = 0;  // counter.


   if(contours.size()>1){    
   for( int i = 0; i< contours.size(); i++ ){
      if(hierarchy[i][3]==largest_contour_index){
        ventricles_indexes.push_back(i);
        ventricles_size.push_back(contourArea( contours[i],false));
      }
   }


   std::vector<int> ventricles_selected; 
   ventricles_selected.push_back(largest_contour_index);

   std::partial_sort(ventricles_size.begin(), ventricles_size.begin()+2,ventricles_size.end(), std::greater<int>());

   for(int j=0; j < numNested; j++){
    for (unsigned i=0; i < ventricles_indexes.size(); i++) {
      if( ventricles_size[j] == contourArea( contours[ventricles_indexes[i]],false) ){
            ventricles_selected.push_back(ventricles_indexes[i]);
            //printf("%d, ", ventricles_indexes[i]);
        }
      }
    }
   //printf("\n"); 
   std::sort( ventricles_selected.begin(), ventricles_selected.end() );
   ventricles_selected.erase( unique( ventricles_selected.begin(), ventricles_selected.end() ), ventricles_selected.end() );


   

       printf("Ventricles: "); 
  for (it=ventricles_selected.begin() ; it < ventricles_selected.end(); it++,counter++ ){
    if(counter>0){
      printf("%d, ",*it);
    }
   for (int j = 0; j < contours[*it].size(); ++j){
      // Do whatever you need to do with the points in the ith contour
      xPoint.push_back(contours[*it][j].x);
      yPoint.push_back(contours[*it][j].y);
      contourID.push_back(counter);
   }
      xPoint.push_back(contours[*it][0].x); //make sure polygon is closed
      yPoint.push_back(contours[*it][0].y);
      contourID.push_back(counter);
  }
  printf("\n"); 
  }else{
    
      for (int j = 0; j < contours[0].size(); ++j){
      // Do whatever you need to do with the points in the ith contour
        xPoint.push_back(contours[0][j].x);
        yPoint.push_back(contours[0][j].y);
        contourID.push_back(counter);
      }
      xPoint.push_back(contours[0][0].x); //make sure polygon is closed
      yPoint.push_back(contours[0][0].y);
      contourID.push_back(counter);
  }
  //int rows = img.rows;
  //int cols = img.cols;
  if(DISPLAY){
     int k;
  while(1){

  k=waitKey(0);
  //cout << '\n' << "KEY PRESSED: " << k << endl;
  if( (k == 27)|(k == -1)|(k == 115) ){
    cout << '\n' << "Assembling output list" << endl;
  
    destroyWindow("Original");
    destroyWindow("Binary segmentation");
    destroyWindow("Contours");
    break;
  }
  }
  }
  
  return List::create(
    _["x"] = xPoint,
    _["y"] = yPoint,
    _["contour.ID"] = contourID
  );

   /*
  END
  */
END_RCPP  
}



