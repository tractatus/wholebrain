#include <Rcpp.h>
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;

class LaplacianBlending {
private:
    Mat_<float> left;
    Mat_<float> right;
    Mat_<float> blendMask;
     
    vector<Mat_<float> > leftLapPyr,rightLapPyr,resultLapPyr;
    Mat leftSmallestLevel, rightSmallestLevel, resultSmallestLevel;
    vector<Mat_<float> > maskGaussianPyramid; //masks are 3-channels for easier multiplication with RGB
     
    int levels;
     
     
    void buildPyramids() {
        buildLaplacianPyramid(left,leftLapPyr,leftSmallestLevel);
        buildLaplacianPyramid(right,rightLapPyr,rightSmallestLevel);
        buildGaussianPyramid();
    }
     
    void buildGaussianPyramid() {
        assert(leftLapPyr.size()>0);
         
        maskGaussianPyramid.clear();
        Mat currentImg;
        //cvtColor(blendMask, currentImg, CV_GRAY2BGR);
        currentImg = blendMask.clone();
        maskGaussianPyramid.push_back(currentImg); //highest level
         
        currentImg = blendMask;
        for (int l=1; l<levels+1; l++) {
            Mat _down; 
            if (leftLapPyr.size() > l) {
                pyrDown(currentImg, _down, leftLapPyr[l].size());
            } else {
                pyrDown(currentImg, _down, leftSmallestLevel.size()); //smallest level
            }
 
            Mat down; 
            //cvtColor(_down, down, CV_GRAY2BGR);
            down = _down.clone();
            maskGaussianPyramid.push_back(down);
            currentImg = _down;
        }
    }
     
    void buildLaplacianPyramid(const Mat& img, vector<Mat_<float> >& lapPyr, Mat& smallestLevel) {
        lapPyr.clear();
        Mat currentImg = img;
        for (int l=0; l<levels; l++) {
            Mat down,up;
            pyrDown(currentImg, down);
            pyrUp(down, up, currentImg.size());
            Mat lap = currentImg - up;
            lapPyr.push_back(lap);
            currentImg = down;
        }
        currentImg.copyTo(smallestLevel);
    }
     
    Mat_<float> reconstructImgFromLapPyramid() {
        Mat currentImg = resultSmallestLevel;
        for (int l=levels-1; l>=0; l--) {
            Mat up;
                         
            pyrUp(currentImg, up, resultLapPyr[l].size());
            currentImg = up + resultLapPyr[l];
        }
        return currentImg;
    }
     
    void blendLapPyrs() {
        resultSmallestLevel = leftSmallestLevel.mul(maskGaussianPyramid.back())  + 
        rightSmallestLevel.mul(1.0 - maskGaussianPyramid.back());  //rightSmallestLevel.mul(Scalar(1.0,1.0,1.0) - maskGaussianPyramid.back()); */
        for (int l=0; l<levels; l++) {
            Mat A = leftLapPyr[l].mul(maskGaussianPyramid[l]);
            Mat antiMask = 1.0 - maskGaussianPyramid[l]; //Scalar(1.0,1.0,1.0) - maskGaussianPyramid[l];
            Mat B = rightLapPyr[l].mul(antiMask);
            Mat blendedLevel = A + B;
            //Mat_<Vec3f> blendedLevel = A + B;
             
            resultLapPyr.push_back(blendedLevel);

        }
    }
 
public:
    LaplacianBlending(const Mat_<float>& _left, const Mat_<float>& _right, const Mat_<float>& _blendMask, int _levels):
    left(_left),right(_right),blendMask(_blendMask),levels(_levels) 
    { 
        assert(_left.size() == _right.size());
        assert(_left.size() == _blendMask.size());

        buildPyramids();

        blendLapPyrs();
    };
     
    Mat_<float> blend() {

        return reconstructImgFromLapPyramid();
    }   
};
 
Mat_<float> LaplacianBlend(const Mat_<float>& l, const Mat_<float>& r, const Mat_<float>& m) {
    LaplacianBlending lb(l,r,m,4);
    return lb.blend();
}


Mat_<ushort> multiBandBlending(Mat src1, Mat src2, bool vertical) {
  if ( src1.empty() || src2.empty() )
    {
     cerr << "cannot get image overlap" << endl;
    }
   Mat_<float> dst1; src1.convertTo(dst1, CV_32F, 1.0/65535.0);
   Mat_<float> dst2; src2.convertTo(dst2, CV_32F, 1.0/65535.0);
   Mat_<float> m(dst1.rows,dst2.cols,0.0);
   

   if(vertical){
     m(Range(0,m.rows/2),Range::all()) = 1.0;
   }else{
     m(Range::all(),Range(0,m.cols/2)) = 1.0;
   }

   Mat_<float> blend = LaplacianBlend(dst1, dst2, m);
   Mat output;
   blend.convertTo(output,CV_16U, 65535.0);   
   return output;
}



RcppExport SEXP testLaplace() {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>

   Mat l8u = imread("/Users/danielfurth/Documents/bottom_FFC_test0033.tif", -1);
   Mat r8u = imread("/Users/danielfurth/Documents/top_FFC_test0038.tif", -1);

   if ( l8u.empty() || l8u.empty() )
    {
     cerr << "cannot find file, check image path of: " << "/Users/danielfurth/Documents/bottom_FFC_test0033.tif" << endl;
    }

    if ( r8u.empty() || l8u.empty() )
    {
     cerr << "cannot find file, check image path of: " << "/Users/danielfurth/Documents/top_FFC_test0038.tif" << endl;
    }

   Mat_<float> l; l8u.convertTo(l,CV_32F,1.0/65535.0);
   Mat_<float> r; r8u.convertTo(r,CV_32F,1.0/65535.0);
   Rcpp::Rcout << "Converted " << std::endl;
 


     Mat_<float> m(l.rows,l.cols,0.0);
     m(Range(0,m.rows/2),Range::all()) = 1.0;


   Mat_<float> blend = LaplacianBlend(l, r, m);
   imshow("blended",blend);
   waitKey(0);
   Mat output;
   blend.convertTo(output,CV_16U, 65535.0);
   imwrite("/Users/danielfurth/Documents/blended_FFC.tif",output);
END_RCPP  
}

void placeCenterTiles ( int l, int m, int n, int overlap, Rcpp::IntegerVector tilePosTop, Rcpp::IntegerVector tilePosBottom, Rcpp::IntegerVector tilePosLeft, Rcpp::IntegerVector tilePosRight, Rcpp::IntegerVector x0, Rcpp::IntegerVector x1, Rcpp::IntegerVector y0, Rcpp::IntegerVector y1, vector<Mat>& src, Mat& dst, bool verbose)

//****************************************************************************80
//
//  Purpose:
//
//    getGain computes the z-stack wise quantile of a image position and average across these samples.
//
//  Discussion:
//
//    On input, the A array contains values of 0 or 1.
//
//    The 0 pixels are to be ignored.  The 1 pixels are to be grouped

//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 November 2015
//
//  Author:
//
//    Daniel Furth
//
//  Parameters:
//
//    Input, int L, M, N, the order of the array.
//
//
{



  int STACKS = l; //number of tiles
  

//  "Read" the array one pixel at a time.
//
int barWidth = 70;
float progress = 0.0;


          for(size_t s = 0; s != STACKS; ++s){
            int x;
            int y;
            int w;
            int h;

            
            if(tilePosTop[s]==1){
              //top tile
              y = 0;
              h = src.at(s).rows - overlap;
            }else{
              y = overlap;
              h = src.at(s).rows - 2*overlap;
            }
            if(tilePosBottom[s]==1){
              //bottom tile
              h = src.at(s).rows - overlap;
            }
            if(tilePosLeft[s]==1){
              //left tile
              x = 0;
              w = src.at(s).cols - overlap;
            }else{
              x = overlap;
              w = src.at(s).cols - 2*overlap;
            }
            if(tilePosRight[s]==1){
              //right tile
              w = src.at(s).cols - overlap;
            }
            
            Rect region_without_borders = Rect(x, y, w, h);
            Mat image_roi = src.at(s)(region_without_borders);

            image_roi.copyTo(dst(cv::Rect(x0[s] + x, y0[s] + y, w, h)));
                //print progress bar in console
            if(verbose){
    Rcpp::Rcout << "  [";
    int pos = barWidth * progress;
    for (int p = 0; p < barWidth; ++p) {
        if (p < pos) Rcpp::Rcout << "=";
        else if (p == pos) Rcpp::Rcout << ">";
        else Rcpp::Rcout << " ";
    }
        Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
    progress += (float)1/(STACKS-1);

          }
        }
          

}

void placeHorizontalOverlap ( int l, int m, int n, int overlap, Rcpp::IntegerVector tilePosTop, Rcpp::IntegerVector horizLeftIndex, Rcpp::IntegerVector x0, Rcpp::IntegerVector y0, Rcpp::IntegerVector y1, vector<Mat>& src, Mat& dst, bool verbose)

//****************************************************************************80
//
//  Purpose:
//
//    getGain computes the z-stack wise quantile of a image position and average across these samples.
//
//  Discussion:
//
//    On input, the A array contains values of 0 or 1.
//
//    The 0 pixels are to be ignored.  The 1 pixels are to be grouped

//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 November 2015
//
//  Author:
//
//    Daniel Furth
//
//  Parameters:
//
//    Input, int L, M, N, the order of the array.
//
//
{



  int STACKS = horizLeftIndex.size(); //number of tiles
  

//  "Read" the array one pixel at a time.
//
int barWidth = 70;
float progress = 0.0;


          for(size_t s = 0; s < STACKS; ++s){
            
            int y;
            int h = y1[s] - y0[s];

            
            if(tilePosTop[horizLeftIndex[s]-1]==1){
              //top tile
              y = 0;
            }else{
              y = overlap;
            }
            

            Rect left_ROI = Rect(m-overlap, y, overlap, h);
            Rect right_ROI = Rect(0, y, overlap, h);
            Mat src1 = src.at(horizLeftIndex[s]-1)(left_ROI);
            Mat src2 = src.at(horizLeftIndex[s])(right_ROI);

            Mat blended = multiBandBlending(src1, src2, false);
            //Rcpp::Rcout << "blended.cols: " <<  blended.cols << ", blended.rows: " << blended.rows << std::endl;
            //Rcpp::Rcout << "S: " <<  s << ", x0: " << x0[s] << ", y0: " << y0[s] << ", overlap: " << overlap << ", h: " << h  << std::endl;

            blended.copyTo(dst(cv::Rect(x0[s], y0[s], overlap, h)));
             // Rcpp::Rcout << "placed" << std::endl;

                //print progress bar in console
            if(verbose){
    Rcpp::Rcout << "  [";
    int pos = barWidth * progress;
    for (int p = 0; p < barWidth; ++p) {
        if (p < pos) Rcpp::Rcout << "=";
        else if (p == pos) Rcpp::Rcout << ">";
        else Rcpp::Rcout << " ";
    }
        Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
    progress += (float)1/(STACKS-1);

          }
        }
          

}

void placeVerticalOverlap ( int l, int m, int n, int overlap, Rcpp::IntegerVector tilePosLeft, Rcpp::IntegerVector verticTopIndex, Rcpp::IntegerVector verticBottomIndex, Rcpp::IntegerVector x0, Rcpp::IntegerVector x1, Rcpp::IntegerVector y0, vector<Mat>& src, Mat& dst, bool verbose)

//****************************************************************************80 tilePosLeft, verticTopIndex, verticBottomIndex, verticX0, verticX1, verticY0
//
//  Purpose:
//
//    getGain computes the z-stack wise quantile of a image position and average across these samples.
//
//  Discussion:
//
//    On input, the A array contains values of 0 or 1.
//
//    The 0 pixels are to be ignored.  The 1 pixels are to be grouped

//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 November 2015
//
//  Author:
//
//    Daniel Furth
//
//  Parameters:
//
//    Input, int L, M, N, the order of the array.
//
//
{



  int STACKS = verticTopIndex.size(); //number of tiles
  

//  "Read" the array one pixel at a time.
//
int barWidth = 70;
float progress = 0.0;


          for(size_t s = 0; s < STACKS; ++s){
            int x;
            int w = x1[s] - x0[s];

            if(tilePosLeft[verticTopIndex[s]-1]==1){
              //top tile
              x = 0;
            }else{
              x = overlap;
            }



            Rect top_ROI = Rect(x, n-overlap, w, overlap);
            Rect bottom_ROI = Rect(x, 0, w, overlap);
            Mat src1 = src.at(verticTopIndex[s]-1)(top_ROI);
            Mat src2 = src.at(verticBottomIndex[s]-1)(bottom_ROI);

            Mat blended = multiBandBlending(src1, src2, true);


            blended.copyTo(dst(cv::Rect(x0[s], y0[s], w, overlap)));
                //print progress bar in console
    if(verbose){
    Rcpp::Rcout << "  [";
    int pos = barWidth * progress;
    for (int p = 0; p < barWidth; ++p) {
        if (p < pos) Rcpp::Rcout << "=";
        else if (p == pos) Rcpp::Rcout << ">";
        else Rcpp::Rcout << " ";
    }
        Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
    progress += (float)1/(STACKS-1);

          }
        }
          

}

void placeSmallOverlap ( int l, int m, int n, int overlap, Rcpp::IntegerVector smallTopLeft, Rcpp::IntegerVector smallTopRight, Rcpp::IntegerVector smallBottomLeft, Rcpp::IntegerVector smallBottomRight, Rcpp::IntegerVector x0, Rcpp::IntegerVector y0, vector<Mat>& src, Mat& dst, bool verbose)

//****************************************************************************80 tilePosLeft, verticTopIndex, verticBottomIndex, verticX0, verticX1, verticY0
//
//  Purpose:
//
//    getGain computes the z-stack wise quantile of a image position and average across these samples.
//
//  Discussion:
//
//    On input, the A array contains values of 0 or 1.
//
//    The 0 pixels are to be ignored.  The 1 pixels are to be grouped

//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 November 2015
//
//  Author:
//
//    Daniel Furth
//
//  Parameters:
//
//    Input, int L, M, N, the order of the array.
//
//
{

  int STACKS = smallTopLeft.size(); //number of tiles
  

//  "Read" the array one pixel at a time.
//
int barWidth = 70;
float progress = 0.0;


          for(size_t s = 0; s < STACKS; ++s){
          


            Rect top_left_ROI = Rect(m-overlap, n-overlap, overlap, overlap);
            Rect top_right_ROI = Rect(0, n-overlap, overlap, overlap);
            Mat src1 = src.at(smallTopLeft[s]-1)(top_left_ROI);
            Mat src2 = src.at(smallTopRight[s]-1)(top_right_ROI);

            Mat top = multiBandBlending(src1, src2, false);

            Rect bottom_left_ROI = Rect(m-overlap, 0, overlap, overlap);
            Rect bottom_right_ROI = Rect(0, 0, overlap, overlap);
            Mat src3 = src.at(smallBottomLeft[s]-1)(bottom_left_ROI);
            Mat src4 = src.at(smallBottomRight[s]-1)(bottom_right_ROI);

            Mat bottom = multiBandBlending(src3, src4, false);

            Mat blended = multiBandBlending(top, bottom, true);


            blended.copyTo(dst(cv::Rect(x0[s], y0[s], overlap, overlap)));
                //print progress bar in console
            if(verbose){
    Rcpp::Rcout << "  [";
    int pos = barWidth * progress;
    for (int p = 0; p < barWidth; ++p) {
        if (p < pos) Rcpp::Rcout << "=";
        else if (p == pos) Rcpp::Rcout << ">";
        else Rcpp::Rcout << " ";
    }
        Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::flush;
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
    progress += (float)1/(STACKS-1);

          }
        }
          

}

/* apply operation to stack */
RcppExport SEXP LaplacianBlendPipe(SEXP input, SEXP outname, SEXP nrows, SEXP ncols, SEXP overlappixels, SEXP widthPixels, SEXP heightPixels, SEXP topGridLayout, SEXP bottomGridLayout, SEXP leftGridLayout, SEXP rightGridLayout, SEXP x0GridLayout,   SEXP x1GridLayout,   SEXP y0GridLayout, SEXP y1GridLayout, SEXP horizleftImage, SEXP x0horiz, SEXP y0horiz, SEXP y1horiz, SEXP vertictopImage, SEXP verticbottomImage, SEXP x0vertic,  SEXP  x1vertic , SEXP  y0vertic, SEXP topleftImage, SEXP toprightImage, SEXP bottomleftImage, SEXP bottomrightImage, SEXP x0small, SEXP y0small, SEXP imagedisplay, SEXP writetoconsole, SEXP contrast, SEXP brightness, SEXP matching, SEXP rotate,SEXP outputfile) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector f(input);
  int num_files = f.size();

  Rcpp::CharacterVector outf(outname);
  std::string outfolder(outf[0]);

  int verticaltiles = Rcpp::as<int>(nrows);
  int horizontaltiles = Rcpp::as<int>(ncols);
  int pixeloverlap = Rcpp::as<int>(overlappixels);
  int stitchedimageWidth = Rcpp::as<int>(widthPixels);
  int stitchedimageHeight = Rcpp::as<int>(heightPixels);

  int shouldImageBeShown = Rcpp::as<int>(imagedisplay);
  bool verbose = Rcpp::as<bool>(writetoconsole);

  double alpha = Rcpp::as<double>(contrast);
  int beta = Rcpp::as<int>(brightness);
  double angle = Rcpp::as<double>(rotate);

  Rcpp::IntegerVector tilePosTop(topGridLayout);
  Rcpp::IntegerVector tilePosBottom(bottomGridLayout);
  Rcpp::IntegerVector tilePosLeft(leftGridLayout);
  Rcpp::IntegerVector tilePosRight(rightGridLayout);

  Rcpp::IntegerVector tilePosx0(x0GridLayout);
  Rcpp::IntegerVector tilePosx1(x1GridLayout);
  Rcpp::IntegerVector tilePosy0(y0GridLayout);
  Rcpp::IntegerVector tilePosy1(y1GridLayout);

  Rcpp::IntegerVector horizLeftIndex(horizleftImage);
  Rcpp::IntegerVector horizX0(x0horiz);
  Rcpp::IntegerVector horizY0(y0horiz);
  Rcpp::IntegerVector horizY1(y1horiz);

  Rcpp::IntegerVector verticTopIndex(vertictopImage);
  Rcpp::IntegerVector verticBottomIndex(verticbottomImage);  
  Rcpp::IntegerVector verticX0(x0vertic);
  Rcpp::IntegerVector verticX1(x1vertic);
  Rcpp::IntegerVector verticY0(y0vertic);

  Rcpp::IntegerVector smallTopLeft(topleftImage);
  Rcpp::IntegerVector smallTopRight(toprightImage); 
  Rcpp::IntegerVector smallBottomLeft(bottomleftImage);  
  Rcpp::IntegerVector smallBottomRight(bottomrightImage);   
  Rcpp::IntegerVector smallX0(x0small);
  Rcpp::IntegerVector smallY0(y0small);

  int featurematching = Rcpp::as<int>(matching);
  

  Rcpp::CharacterVector of(outputfile);
  std::string off(of[0]);

  int barWidth = 70;
  float progress = 0.0;

  std::string ff(f[0]);
  Mat remove = imread(ff, -1);
  int rows = remove.rows;
  int cols = remove.cols;
  vector<Mat> zpositions(num_files,Mat(cols,rows,remove.type()));
  vector<Mat> imgs(num_files,Mat(cols,rows,CV_8UC1));
  vector<Mat> destination(num_files,Mat(cols,rows,remove.type()));
  Mat src;
  
  Mat dst = Mat::zeros(Size(stitchedimageWidth, stitchedimageHeight), remove.type());
  dst = dst + 65535;
  if(verbose){Rcpp::Rcout << "Loading " <<  num_files << " tiles into RAM." << std::endl;}
  for (int i = 0; i < f.size(); i++){
    std::string filename(f[i]);

    src = imread(filename, -1);
    if(featurematching){
      double Min, Max;
      src.convertTo(src, -1, alpha, beta);  
      minMaxLoc(src, &Min, &Max);
      src.convertTo(src,CV_8UC1,255.0/(Max-Min)); 
      imgs.push_back(src);
    }else{
    zpositions.at(i) = src;
    }

    


    //print progress bar in console
    if(verbose){
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
    progress += (float)1/(num_files-1);
    }

  }


  if(featurematching){
    //Mat pano;
    //bool try_use_gpu = false;
    if(verbose){Rcpp::Rcout << "\n" << std::endl;
  Rcpp::Rcout << "====== RUNNING FEATURE MATCHING STITCHER ======" << std::endl;}
    /*Remove THIS BECAUSE NOT SUPPORTED MODULE Stitcher stitcher = Stitcher::createDefault(try_use_gpu);
    Stitcher::Status status = stitcher.stitch(imgs, pano);
 
    if (status != Stitcher::OK)
    {
           Rcpp::Rcout <<  "Can't stitch images, error code = " << status << std::endl;
           return R_NilValue;
    }
      if(verbose){Rcpp::Rcout << "====== FEATURE MATCHING DONE ======" << std::endl;
  */
  Rcpp::Rcout << "feature matching not supported anymore..." << std::endl;
    //}
/*
  string outputname;
  outputname =  outfolder + "/" + off;
  if(angle!=0){

    // get rotation matrix for rotating the image around its center
    cv::Point2f center(pano.cols/2.0, pano.rows/2.0);
    cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1.0);
    // determine bounding rectangle
    cv::Rect bbox = cv::RotatedRect(center,pano.size(), angle).boundingRect();
    // adjust transformation matrix
    rot.at<double>(0,2) += bbox.width/2.0 - center.x;
    rot.at<double>(1,2) += bbox.height/2.0 - center.y;

    cv::Mat dst;
    cv::warpAffine(src, dst, rot, bbox.size());
  }else{
    imwrite(outputname, pano);
  }

  if(verbose){Rcpp::Rcout << "====== SAVING DONE ======" << std::endl;}
  
  if(shouldImageBeShown){
  int k;  
  string displaywindow = off;
  Mat normalized;
  double Min, Max;
  pano.convertTo(pano, -1, alpha, beta);  
  minMaxLoc(pano, &Min, &Max);
  pano.convertTo(normalized,CV_8UC1,255.0/(Max-Min)); 
  if(normalized.rows>600){
          double ScaleFactor = 600/(double)normalized.rows;
          resize(normalized, normalized, Size(), ScaleFactor, ScaleFactor, INTER_LINEAR);
        }
  namedWindow( displaywindow, cv::WINDOW_AUTOSIZE);
  imshow(displaywindow, normalized);
  Rcpp::Rcout << '\n' << "Press ESC to close Display window of "<< displaywindow << std::endl;
  while(k > 0){
  k=waitKey(0);
  Rcpp::Rcout << '\n' << "KEY PRESSED: " << k << std::endl;
  if( (k == 27)|(k == -1) ){
  destroyWindow(displaywindow);
  break;
  }
  }
  }
*/

  }else{


  if(verbose){Rcpp::Rcout << "\n" << std::endl;
  Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;

  Rcpp::Rcout << "Running command. " << "place center tiles" << std::endl;}
  placeCenterTiles(num_files, cols, rows, pixeloverlap, tilePosTop, tilePosBottom, tilePosLeft, tilePosRight, tilePosx0, tilePosx1, tilePosy0, tilePosy1, zpositions, dst, verbose);
    // Gaussian smoothing

  if(verbose){Rcpp::Rcout << "\n" << std::endl;
  Rcpp::Rcout << "====== CENTER TILES DONE ======" << std::endl;

  Rcpp::Rcout << "Running command. " << "blend horizontal overlap" << std::endl;}

  placeHorizontalOverlap(num_files, cols, rows, pixeloverlap, tilePosTop, horizLeftIndex, horizX0, horizY0, horizY1, zpositions, dst, verbose);

  if(verbose){Rcpp::Rcout << "====== HORIZONTAL OVERLAP DONE ======" << std::endl;

  Rcpp::Rcout << "Running command. " << "blend vertical overlap" << std::endl;}

  placeVerticalOverlap(num_files, cols, rows, pixeloverlap, tilePosLeft, verticTopIndex, verticBottomIndex, verticX0, verticX1, verticY0, zpositions, dst, verbose);

  if(verbose){Rcpp::Rcout << "====== VERTICAL OVERLAP DONE ======" << std::endl;

  Rcpp::Rcout << "Running command. " << "blend small overlap" << std::endl;}

  placeSmallOverlap(num_files, cols, rows, pixeloverlap, smallTopLeft, smallTopRight, smallBottomLeft, smallBottomRight, smallX0, smallY0, zpositions, dst, verbose);

  if(verbose){Rcpp::Rcout << "====== SMALL OVERLAP DONE ======" << std::endl;

  Rcpp::Rcout << "saving stitched image..." << std::endl;
   R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
  }

  cv::Mat rotated;
  if(angle!=0){
    // get rotation matrix for rotating the image around its center
    cv::Point2f center(dst.cols/2.0, dst.rows/2.0);
    cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1.0);
    // determine bounding rectangle
    cv::Rect bbox = cv::RotatedRect(center,dst.size(), angle).boundingRect();
    // adjust transformation matrix
    rot.at<double>(0,2) += bbox.width/2.0 - center.x;
    rot.at<double>(1,2) += bbox.height/2.0 - center.y;

    
    cv::warpAffine(dst, rotated, rot, bbox.size());
    string outputname;
    outputname =  outfolder + "/" + off;
    imwrite(outputname, rotated);
  }else{

    dst.copyTo(rotated);
    string outputname;
    outputname =  outfolder + "/" + off;
    imwrite(outputname, rotated);
  }

  if(verbose){Rcpp::Rcout << "====== SAVING DONE ======" << std::endl;}
  
  if(shouldImageBeShown){
  int k;  
  string displaywindow = off;
  Mat normalized;
  double Min, Max;
  rotated.convertTo(rotated, -1, alpha, beta);  
  minMaxLoc(rotated, &Min, &Max);
  rotated.convertTo(normalized,CV_8UC1,255.0/(Max-Min)); 
  if(normalized.rows>600){
          double ScaleFactor = 600/(double)normalized.rows;
          resize(normalized, normalized, Size(), ScaleFactor, ScaleFactor, INTER_LINEAR);
        }
  namedWindow( displaywindow, cv::WINDOW_AUTOSIZE);
  imshow(displaywindow, normalized);
  Rcpp::Rcout << '\n' << "Press ESC to close Display window of "<< displaywindow << std::endl;
  while(k > 0){
  k=waitKey(0);
  Rcpp::Rcout << '\n' << "KEY PRESSED: " << k << std::endl;
  if( (k == 27)|(k == -1) ){
  destroyWindow(displaywindow);
  break;
  }
  }
  }
   
   //else feature matching ended
  }

  return R_NilValue;
END_RCPP  
}



