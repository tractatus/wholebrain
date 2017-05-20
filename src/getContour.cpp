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
RcppExport SEXP getContour(SEXP input, SEXP thresho, SEXP invertimg, SEXP getLarge, SEXP getNested, SEXP display, SEXP resizeP, SEXP blurImg, SEXP writetoconsole, SEXP saveoutput, SEXP writeout) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);

  int saveout = Rcpp::as<int>(writeout);

  Rcpp::CharacterVector outputfilename(saveoutput);
  std::string ffoutputfilename(outputfilename[0]);

  int thresh = Rcpp::as<int>(thresho);
  int invert = Rcpp::as<int>(invertimg);
  int DISPLAY = Rcpp::as<int>(display);
  int numNested = Rcpp::as<int>(getNested);

  bool verbose = Rcpp::as<bool>(writetoconsole);

  double resizeParam = Rcpp::as<double>(resizeP);
  int blurImage = Rcpp::as<int>(blurImg);



  if(verbose){Rcpp::Rcout << "Loading image:" << ffname << std::endl;}
  Mat img = imread(ffname, -1);
  if(verbose){Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;}

   if(verbose){
    Rcpp::Rcout << "Resizing to: " <<  resizeParam << "% of original size." << std::endl;
    Rcpp::Rcout << "Image type: " <<  ImgTypes(img.type()) << "_" << img.type()  << std::endl;
  }
  resize(img, img, Size(), resizeParam, resizeParam, INTER_LINEAR);

  if(DISPLAY){
    imshow("Original", img);
  }
  if(img.type()==16){
    cvtColor(img , img , CV_BGR2GRAY);
    if(invert){
      bitwise_not ( img, img );
    }
    Rcpp::Rcout << "Image type: " <<  ImgTypes(img.type()) << "_" << img.type()  << std::endl;
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

    printf("LOOP DONE: %d \n", contours.size() ); 

   std::vector<int> ventricles_selected; 
   ventricles_selected.push_back(largest_contour_index);


   if(numNested>0){

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
  }

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
      printf("Contours less than 1: "); 
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

Mat createSegmentationDisplay(Mat& segments, int numOfSegments, Mat const& image)
{
  // Create a new image
  Mat wshed(segments.size(), CV_8UC3);
 
  //Create color tab for coloring the segments
  vector<Vec3b> colorTab;
  for(int i = 0; i < numOfSegments; i++) {
    int b = theRNG().uniform(0, 255);
    int g = theRNG().uniform(0, 255);
    int r = theRNG().uniform(0, 255);
    colorTab.push_back(Vec3b((uchar)b, (uchar)g, (uchar)r));
  }
 
  // Assign different color to different segments
  for(int i = 0; i < segments.rows; i++) {
    for(int j = 0; j < segments.cols; j++) {
      int index = segments.at<int>(i,j);
      if(index == -1) {
        wshed.at<Vec3b>(i,j) = Vec3b(255,255,255);
      }
      else if(index <= 0 || index > numOfSegments) {
        wshed.at<Vec3b>(i,j) = Vec3b(0,0,0);
      }
      else {
        wshed.at<Vec3b>(i,j) = colorTab[index - 1];
      }
    }
  }
 
  //If the original image available then merge with the colors of segments
  if(image.dims > 0) {
    wshed = (wshed * 0.5) + (image * 0.5);
  }
 
  return wshed;
}


class WatershedSegmenter{
private:
    Mat markers;
public:
    void setMarkers(Mat& markerImage)
    {
        markerImage.convertTo(markers, CV_32S);
    }

    Mat process(Mat &image)
    {
        watershed(image, markers);
        markers.convertTo(markers, CV_8U);
        return markers;
    }
};
// NUCLER SEGMENT
RcppExport SEXP nuclearSegment(SEXP input, SEXP distanceThresh, SEXP kernel, SEXP iter, SEXP thresho, SEXP invertimg, SEXP getLarge, SEXP getNested, SEXP display, SEXP resizeP, SEXP blurImg, SEXP writetoconsole, SEXP saveoutput, SEXP writeout) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);

  int saveout = Rcpp::as<int>(writeout);

  Rcpp::CharacterVector outputfilename(saveoutput);
  std::string ffoutputfilename(outputfilename[0]);

  int thresh = Rcpp::as<int>(thresho);
  int invert = Rcpp::as<int>(invertimg);
  int DISPLAY = Rcpp::as<int>(display);
  int numNested = Rcpp::as<int>(getNested);

  bool verbose = Rcpp::as<bool>(writetoconsole);

  double resizeParam = Rcpp::as<double>(resizeP);

  double distansThresh = Rcpp::as<double>(distanceThresh);
  int blurImage = Rcpp::as<int>(blurImg);

  int morph_size = Rcpp::as<int>(kernel);
  int ITERATIONS = Rcpp::as<int>(iter);;


  if(verbose){Rcpp::Rcout << "Loading image:" << ffname << std::endl;}
  Mat img = imread(ffname, -1);
  if(verbose){Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;}

   if(verbose){
    Rcpp::Rcout << "Resizing to: " <<  resizeParam << "% of original size." << std::endl;
    Rcpp::Rcout << "Image type: " <<  ImgTypes(img.type()) << "_" << img.type()  << std::endl;
  }
  resize(img, img, Size(), resizeParam, resizeParam, INTER_LINEAR);
  Mat originalImag;
  img.copyTo(originalImag);
  if(DISPLAY){
    imshow("Original", img);
  }
  if(img.type()==16){
    cvtColor(img , img , CV_BGR2GRAY);
    if(invert){
      bitwise_not ( img, img );
    }
    Rcpp::Rcout << "Image type: " <<  ImgTypes(img.type()) << "_" << img.type()  << std::endl;
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


  Mat gx, gy;
  Scharr(img, gx, CV_32F, 1, 0); //CV_32F
  Scharr(img, gy, CV_32F, 0, 1);


  /// Generate grad_x and grad_y
  Mat grad_x, grad_y;
  Mat abs_grad_x, abs_grad_y;

  /// Gradient X
  //Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
  Scharr( img, grad_x, CV_32F, 1, 0, 1, 0, BORDER_DEFAULT );
  convertScaleAbs( grad_x, abs_grad_x );

  /// Gradient Y
  //Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
  Scharr( img, grad_y, CV_32F, 0, 1, 1, 0, BORDER_DEFAULT );
  convertScaleAbs( grad_y, abs_grad_y );

  /// Total Gradient (approximate)
  Mat grad;
  addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
  //resize(grad, grad, Size(), 0.25, 0.25);

  //imshow( "Sobel test", grad );
  

  // Compute the structure tensor, and from it, anisotropy
  int Sigma =10;
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
    double maxVal, minVal;

  minMaxLoc(energy, &minVal, &maxVal);
  energy.convertTo(img, CV_16S, 65536.0/(maxVal - minVal), -minVal * 65536.0/(maxVal - minVal));
if(DISPLAY){
    imshow("energy", img);
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


       //imwrite(ffoutputfilename, img);

  // Perform the distance transform algorithm


  Mat element = getStructuringElement( MORPH_ELLIPSE, Size( 2*morph_size + 1, 2*morph_size+1 ), Point( morph_size, morph_size ) );
  Mat sureBG;
  morphologyEx(img_bw, sureBG, MORPH_CLOSE, element );
   for( int i = 0; i< ITERATIONS; i++ ){
      dilate(sureBG,sureBG,element);

   }
   if(DISPLAY){
      imshow("Sure BG", sureBG);
    }


  Mat distanceTest;
  Mat distantrans;
  img_bw.convertTo(distanceTest, CV_8UC1);

  distanceTransform(distanceTest, distantrans, CV_DIST_L2, 5);
  // Normalize the distance image on the unit interval
    normalize(distantrans, distantrans, 0, 1., NORM_MINMAX);
  if(DISPLAY){
    imshow("Distance Transform Image", distantrans);
  }
threshold(distantrans, distantrans, distansThresh, 1., CV_THRESH_BINARY);
  distantrans.convertTo(distantrans, CV_8U);

  // Find total markers
  vector<vector<Point> > contoursDist;
  findContours(distantrans, contoursDist, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
 
  // Create the marker image for the watershed algorithm
  Mat markers = Mat::zeros(distantrans.size(), CV_32SC1);
 
  // Draw the foreground markers
  for (size_t i = 0; i < contoursDist.size(); i++)
    drawContours(markers, contoursDist, static_cast<int>(i), Scalar::all(static_cast<int>(i)+1), -1);
 
  // Draw the background marker
  circle(markers, Point(5,5), 3, CV_RGB(255,255,255), -1);
  //imshow("Markers", markers*10000);
  Mat dist_CV_8UC3;
  img.convertTo(dist_CV_8UC3,CV_8UC1,255.0/(depth)); //-255.0/Min
    cvtColor(dist_CV_8UC3,dist_CV_8UC3,CV_GRAY2RGB);

  // Perform the watershed algorithm
  watershed(dist_CV_8UC3, markers);
  Mat mark = Mat::zeros(markers.size(), CV_8UC1);
  markers.convertTo(mark, CV_8UC1);
  bitwise_not(mark, mark);
  if(DISPLAY){
  imshow("Markers_v2", mark); // uncomment this if you want to see how the mark image looks like at that point
}

  int numOfSegments = contoursDist.size();
  Mat wshed = createSegmentationDisplay(markers, numOfSegments, dist_CV_8UC3);
  if(DISPLAY){
  imshow("wshed", wshed); // uncomment this if you want to see how the mark image looks like at that point
}
  // Get the borders of the objects (beans)
Mat distanceTest2;
  Mat distantransT;
  img_bw.convertTo(distanceTest2, CV_8UC1);

  distanceTransform(distanceTest2, distantransT, CV_DIST_L2, 5);
  // Normalize the distance image on the unit interval
    normalize(distantransT, distantransT, 0, 1., NORM_MINMAX);
    if(DISPLAY){
    imshow("Distance Transform Image2", distantransT);
  }
threshold(distantransT, distantransT, distansThresh, 1., CV_THRESH_BINARY);
  //distantransT.convertTo(distantransT, CV_8U);


  minMaxLoc(distantransT, &minVal, &maxVal);
  distantransT.convertTo(distantransT, CV_8U, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));
  if(DISPLAY){
    imshow("Distance Transform Image3", distantransT);
  }

  Mat breans, mupp;
  img_bw.convertTo(breans, CV_8UC1);
    threshold(img, mupp, thresh, depth, 3);

    cv::Mat objects,peaks;
    threshold(mupp, objects, 0, depth, CV_THRESH_BINARY);
  objects.convertTo(objects, CV_8UC1);

  cv::dilate(mupp,peaks,cv::Mat(),cv::Point(-1,-1),3);
    cv::dilate(objects,objects,cv::Mat(),cv::Point(-1,-1),3);

    /* Now all the peaks should be exactely 0*/
    peaks=peaks-mupp;

    /* And the non-peaks 255*/
    cv::threshold(peaks,peaks,0,255,cv::THRESH_BINARY);
    peaks.convertTo(peaks,CV_8U);

    /* Only the zero values of "peaks" that are non-zero
     * in "objects" are the real peaks*/
    cv::bitwise_xor(peaks,objects,peaks);

    /* The peaks that are distant from less than
     * 2 pixels are merged by dilatation */
    cv::dilate(peaks,peaks,cv::Mat(),cv::Point(-1,-1),1);
        std::vector<std::vector<cv::Point> > localMaxima;

    cv::findContours(peaks, localMaxima, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE);
    /* just draw them and save the image */
    cv::Mat dt(mupp.size(), CV_8UC3);

    cv::drawContours(dt,localMaxima,-1,cv::Scalar(255,255,255),-1);
    if(DISPLAY){
        cv::imshow(std::string("dist_CV_8UC3"), dist_CV_8UC3);
      }

    cv::Mat dilation, erosion;
    cv::dilate(sureBG, dilation, cv::Mat(), cv::Point(-1,-1), 3);
    cv::erode(dilation, erosion, cv::Mat());
    cv::Mat border(sureBG.size(), CV_8U, cv::Scalar(0));
    border = dilation - erosion;
    if(DISPLAY){
        cv::imshow(std::string("border"), border);
      }

/*
// Get the distance transform and normalize the result to [0,255]
    cv::Mat dt;
    //bitwise_not(breans,breans);
    cv::distanceTransform(breans, dt, cv::DIST_L2, 3);
    normalize(dt, dt, 0, 1., NORM_MINMAX);

        cv::imshow(std::string("distance"), dt);


    // Threshold it to isolate the peaks
    cv::threshold(dt, dt, distansThresh, 1., cv::THRESH_BINARY);
        cv::imshow(std::string("dt"), dt);

    // Use connectedComponents() to isolate the objects
    dt.convertTo(dt, CV_8U);
            cv::imshow(std::string("dt2"), dt);
*/
              cv::cvtColor(dt, dt, CV_RGB2GRAY);

    cv::Mat lbl(dt.size(), CV_32S);
    int ncc = cv::connectedComponents(dt, lbl, 8, CV_32S);
    std::cout << "Number of Connected Components: " << ncc << std::endl;


     // Create the marker image for watershed
    //cv::Mat markersB = cv::Mat::zeros(dt.size(), CV_8U);
    //markersB = dt + border;
    border.convertTo(border, CV_32SC1);
    cv::normalize(border, border, 0, ncc+1, cv::NORM_MINMAX);

    cv::Mat markersB = cv::Mat::zeros(dt.size(), CV_32SC1);
    markersB = lbl + border; 

    //markersB.convertTo(markersB, CV_32SC1);
      minMaxLoc(originalImag, &minVal, &maxVal);
  originalImag.convertTo(originalImag, CV_8U, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));

    equalizeHist( distantransT, distantransT );

        cvtColor(originalImag,originalImag,CV_GRAY2RGB);

    cv::cvtColor(distantransT, distantransT, CV_GRAY2RGB);
    cv::watershed(distantransT, markersB);

    markersB.convertTo(markersB, CV_8U);
    markersB = 255 - markersB; // cv::bitwise_not()

    // Create a green mask with the sillhuete of the detected objects
    cv::Mat green_mask = cv::Mat::zeros(markersB.size(), CV_8UC3);
    cv::Mat white_mask = cv::Mat::zeros(markersB.size(), CV_8UC3);

    for (int i = 0; i < markersB.rows; i++)
    {
        for (int j = 0; j < markersB.cols; j++)
        {
            // Draw the marker's silhuete (white) as green in the mask
            if (markersB.at<unsigned char>(i,j) == 255){
                green_mask.at<cv::Vec3b>(i,j) = cv::Vec3b(0, 255, 0);
                white_mask.at<cv::Vec3b>(i,j) = cv::Vec3b(255, 255, 255);
              }
        }
    }
    if(DISPLAY){
  imshow("whitemask", white_mask); // uncomment this if you want to see how the mark image looks like at that point
}
    cv::dilate(white_mask, white_mask, cv::Mat());

Mat nuclei = distantransT - white_mask;
if(DISPLAY){
  imshow("nuclei", nuclei); // uncomment this if you want to see how the mark image looks like at that point
}


  vector<vector<Point> > contours;
  vector<Vec4i> hierarchy;
  cv::cvtColor(nuclei, nuclei, CV_RGB2GRAY);

  /// Find contours
  findContours( nuclei, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );
   std::vector<int> xPoint;
   std::vector<int> yPoint;
   std::vector<int> contourID;
   int whichIsPerimeter;
   int counter = 0;

 for (int i = 0; i< contours.size(); ++i){
   for (int j = 0; j < contours[i].size(); ++j){
      // Do whatever you need to do with the points in the ith contour
      xPoint.push_back(contours[i][j].x);
      yPoint.push_back(contours[i][j].y);
      contourID.push_back(counter);
   }
      xPoint.push_back(contours[i][0].x); //make sure polygon is closed
      yPoint.push_back(contours[i][0].y);
      contourID.push_back(counter);
      counter++;
  }

  whichIsPerimeter = counter;
minMaxLoc(dilation, &minVal, &maxVal);

  dilation.convertTo(dilation, CV_8U, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));


  findContours( dilation, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );
  for (int i = 0; i< contours.size(); ++i){
   for (int j = 0; j < contours[i].size(); ++j){
      // Do whatever you need to do with the points in the ith contour
      xPoint.push_back(contours[i][j].x);
      yPoint.push_back(contours[i][j].y);
      contourID.push_back(counter);
   }
      xPoint.push_back(contours[i][0].x); //make sure polygon is closed
      yPoint.push_back(contours[i][0].y);
      contourID.push_back(counter);
      counter++;
  }


int numOfSegments2 = contoursDist.size();
  Mat wshed2 = createSegmentationDisplay(markers, numOfSegments, distantransT);//dist_CV_8UC3);
  if(DISPLAY){
  imshow("output", wshed2); // uncomment this if you want to see how the mark image looks like at that point
  }
    cv::dilate(green_mask, green_mask, cv::Mat());
    if(DISPLAY){
    cv::imshow(std::string("green_mask"), green_mask);
  }


minMaxLoc(border, &minVal, &maxVal);
  border.convertTo(border, CV_8U, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));


    cv::cvtColor(border, border, CV_GRAY2RGB);

    cv::Mat result = originalImag + green_mask + border;

    cv::imshow(std::string("result"), result);
  
/*
Mat labeledImage;
connectedComponents(distantrans, labeledImage, 8, CV_32S);   
    imshow("labeledImage", labeledImage*10000);
     imwrite(ffoutputfilename, labeledImage);

     //watershed only takes 8UC3 as input
     Mat dist_CV_8UC3;
    img_bw.convertTo(dist_CV_8UC3, CV_8UC3);
     labeledImage.convertTo(labeledImage, CV_32SC1);

    watershed(dist_CV_8UC3, labeledImage);
    Mat mark = Mat::zeros(labeledImage.size(), CV_8UC1);
    labeledImage.convertTo(mark, CV_8UC1);
    bitwise_not(mark, mark);
    imshow("Markers_v2", mark);
*/

  
  //int rows = img.rows;
  //int cols = img.cols;
  if(DISPLAY){
     int k;
  while(1){

  k=waitKey(0);
  //cout << '\n' << "KEY PRESSED: " << k << endl;
  if( (k == 27)|(k == -1)|(k == 115) ){
    cout << '\n' << "Assembling output list" << endl;
  

    destroyAllWindows();
    break;
  }
  }
  }else{
     cout << '\n' << "Assembling output list" << endl;
    destroyAllWindows();

  }

  //return R_NilValue;

  return List::create(
    _["x"] = xPoint,
    _["y"] = yPoint,
    _["contour.ID"] = contourID,
    _["perim"] = whichIsPerimeter
  );


END_RCPP  
}








