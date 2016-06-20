// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>

//for sleep
#include <unistd.h>

#include <Rcpp.h>

//#include "ropencv.h"
#include "DoubleSlider.h"
#include "Widget.h"

using namespace cv;
using namespace std;
using namespace Rcpp;
/*
template<> SEXP wrap(const cv::Mat &obj) {

}
*/
#include <ctime>


Mat ui_backgroundImage, ui_display, img2;
void releaseImg(int x,int y);
void showImage();
void findImg(int x,int y, vector <Widget> widgets);

vector <Widget> imshowwidget(1); //this is for command imshow
vector <Widget> widgets(5);
DoubleSlider intensitySliderImage;
DoubleSlider imshowSliderImage;


int sliderMin = 110;
int sliderMax = 692;
int globalCoordinateX[]={sliderMin,sliderMax};	//point coordinates
int globalCoordinateY[]={36,36};
int selectedWidget=-1;			//currently selected point
int numWidgets=2;

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

// take number image type number (from cv::Mat.type()), get OpenCV's enum string.
string getImgTypes(int imgTypeInt)
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

class ImageShow
{
public:
    Mat src; 
    Mat dst;
    Mat dsp;
    Mat out;
    Mat dspHighRes;
    bool displayImage;
    double Max;
    double Min;
    double ScaleFactor;
    int imgdepth;
  
    void updateImage(void)
    {
      
      dst = Mat::zeros(src.size(), CV_8UC1);



      if(displayImage){
        dsp = src.clone();
        dsp -= Min;
        dsp.convertTo(dsp,CV_8UC1,255.0/(Max-Min)); //-255.0/Min

        //REMOVE
        out = dst.clone();
        out.convertTo(out, CV_8UC1);

        dspHighRes = dsp.clone();
        if(src.rows>600){
          ScaleFactor = 600/(double)src.rows;
          resize(dsp, dsp, Size(), ScaleFactor, ScaleFactor, INTER_LINEAR);
        }
      }else{
        out = dst.clone();
        out.convertTo(out, CV_8UC1);
      }
    }
};



class Segmentation
{
public:
    Mat src; 
    Mat dst;
    Mat dsp;
    Mat out;
    Mat dspHighRes;
    bool displayImage;
    int alpha = 10;
    double alphaADJ;
    double Max;
    double Min;
    double ScaleFactor;
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
    vector<float> contourSomaArea;
    vector<float> intensitySoma;

    bool endSegment;

    unsigned int p;
    unsigned int increment;
    unsigned int maxValue;

    void runthreshold(void)
    {
      centroidX.clear();
      centroidY.clear();
      contourSomaArea.clear();
      intensitySoma.clear();
      Mat tmp = Mat::zeros(src.size(), CV_8UC1);
      dst = Mat::zeros(src.size(), CV_8UC1);

      //minMaxIdx(src, &Min, &Max);

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
          thresh = linearspaced[i];//65535*((double)linearspaced[i]/1000);
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
            drawContours( mask, contours, k, color, CV_FILLED, 8, hierarchy );
      
          }
        }
        //add the contours to the destination Mat.
        bitwise_or(dst, mask, dst);
      } 

      findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
      Mat mask = Mat::zeros(src.size(), CV_8UC1);
      Mat labels;
      vector<Moments> mu(contours.size());
      for (size_t k = 0; k < contours.size(); ++k)
      {
        const vector<Point>& contour = contours[k];
        double area0 = contourArea(contour);
        
        //if area falls within the area limits given by R via "alim" parameter then draw the contour.
        if( (area0 > minArea) && (area0 < maxArea) ){
          Scalar color(255, 255, 255);
          drawContours( mask, contours, k, color, CV_FILLED, 8, hierarchy, 0  );

          if(endSegment){
           mu[k] = moments( contours[k], false );
           centroidX.push_back( mu[k].m10/mu[k].m00 ); //
           centroidY.push_back( mu[k].m01/mu[k].m00 ); //
           contourSomaArea.push_back(area0);
            mask.copyTo(labels);
          drawContours(labels, contours, k, Scalar(k), CV_FILLED);
          Rect roi = boundingRect(contours[k]);
          Scalar Avg = cv::mean( src(roi), labels(roi) == k);
          intensitySoma.push_back(Avg[0]);
          }


        }
      }
      dst = mask.clone();



      if(displayImage){
        dsp = src.clone();
        dsp -= Min;
        dsp.convertTo(dsp,CV_8UC1,255.0/(Max-Min)); //-255.0/Min
        
        Scalar red(0, 0, 255);
        cvtColor(dsp, dsp, CV_GRAY2RGB);
        findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
        if(hideFilter){
        }else{
        for( int i = 0; i< contours.size(); i++ ){
          drawContours( dsp, contours, i, red, CV_FILLED, 8, hierarchy, 0 );
        }
        }
        //REMOVE
        out = dst.clone();
        out.convertTo(out, CV_8UC1);
        Scalar color(255, 255, 255);
        findContours(dst, contours, hierarchy, RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE );
        for( int i = 0; i< contours.size(); i++ ){
          drawContours( out, contours, i, color, CV_FILLED, 8, hierarchy, 0 );
        }
        //REMOVE

        dspHighRes = dsp.clone();
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
          drawContours( out, contours, i, color, CV_FILLED, 8, hierarchy, 0 );
        }
      }
    }
};


//template to convert number to string for display in UI
template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

//generate a linear spacing between points
vector<double> linspace(double min, double max, int n)
{
 vector<double> result;
 // vector iterator
 int iterator = 0;

for (int i = 0; i <= n-2; i++)
 {
 double temp = min + i*(max-min)/(floor((double)n) - 1);
 result.insert(result.begin() + iterator, temp);
 iterator += 1;
 }

//iterator += 1;

result.insert(result.begin() + iterator, max);
 return result;
}

void mergewithBackground(Mat* src, Mat* dst, Mat* overlay, const Point& location)
{
    for (int y = max(location.y, 0); y < src->rows; ++y)
    {
        int fY = y - location.y;

        if (fY >= overlay->rows)
            break;

        for (int x = max(location.x, 0); x < src->cols; ++x)
        {
            int fX = x - location.x;

            if (fX >= overlay->cols)
                break;

            double opacity = ((double)overlay->data[fY * overlay->step + fX * overlay->channels() + 3]) / 255;

            for (int c = 0; opacity > 0 && c < src->channels(); ++c)
            {
                unsigned char overlayPx = overlay->data[fY * overlay->step + fX * overlay->channels() + c];
                unsigned char srcPx = src->data[y * src->step + x * src->channels() + c];
                dst->data[y * src->step + src->channels() * x + c] = srcPx * (1. - opacity) + overlayPx * opacity;
            }
        }
    }
}

int insetX0;
int insetY0;
int insetWidth;
int insetHeight;
bool zoomOn = false;
void CallBackZoom(int event, int x, int y, int flags, void* userdata)
{
  switch(event) {
    case EVENT_LBUTTONDOWN:
      cout << '\r' << "LEFT MOUSE DOUBLECLICK (" << x << ", " << y << ")" << flush;
      Segmentation *pd = static_cast<Segmentation *>(userdata);

      Mat zoomImage;
      pd->dsp.copyTo(zoomImage);

      int zoomX0;
      int zoomX1;
      int zoomY0;
      int zoomY1;
      int zoomWidth = zoomImage.cols/10;
      int zoomHeight = zoomImage.rows/10;
      zoomX0 = ((x-zoomWidth)<0) ? 0 : x-zoomWidth;
      zoomY0 = ((y-zoomHeight)<0) ? 0 : y-zoomHeight;
      zoomX1 = ((x+zoomWidth)>zoomImage.cols) ? zoomImage.cols : x+zoomWidth;
      zoomY1 = ((y+zoomHeight)>zoomImage.rows) ? zoomImage.rows : y+zoomHeight;
      rectangle(zoomImage, Point(zoomX0, zoomY0), Point(zoomX1, zoomY1), Scalar( 100, 163, 117), 2);
      imshow("display", zoomImage);

      insetX0 = pd->dspHighRes.cols*((double)zoomX0/zoomImage.cols);
      insetY0 = pd->dspHighRes.rows*((double)zoomY0/zoomImage.rows);
      insetWidth = pd->dspHighRes.cols*2*((double)zoomWidth/zoomImage.cols);
      insetHeight = pd->dspHighRes.rows*2*((double)zoomHeight/zoomImage.rows);


      Mat inset = pd->dspHighRes( Rect(insetX0, insetY0, insetWidth, insetHeight) );
      imshow("zoom", inset);
      zoomOn = true;
      break;
  }

}

//Mat img; // no local
void CallBackFunc(int event, int x, int y, int flags, void* userdata)
{

   switch(event) {
    case EVENT_LBUTTONDBLCLK:
          /*cout << '\r' << "LEFT MOUSE DOUBLECLICK (" << x << ", " << y << ")" << flush;
          for(unsigned int i=0;i<widgets.size();i++){
            //check which widget the user pressed
            if((x>=(widgets[i].x)) && (x<=(widgets[i].x+widgets[i].width ))&& (y>=(widgets[i].y ))&& (y<=(widgets[i].y+widgets[i].height ))){
              selectedWidget=i;
              intensitySliderImage.globalCoordinateX[0] = widgets[selectedWidget].guiPixelValue.at(0);
              intensitySliderImage.globalCoordinateX[1] = widgets[selectedWidget].guiPixelValue.at(1);
              //check if double click is within the slider bar
              int xTmp = x-widgets[selectedWidget].x;
              int yTmp = y-widgets[selectedWidget].y;
              if( (xTmp>(intensitySliderImage.globalCoordinateX[0]+16))&&(xTmp<(intensitySliderImage.globalCoordinateX[1]-4)) &&( (yTmp>=41) && (yTmp<=45) ) ){
                intensitySliderImage.feature = 4;
                cout << '\r' << "INSIDE DOUBLECLICK (" << x << ", " << y << ")" << flush;

                intensitySliderImage.conversionFactorLast = intensitySliderImage.conversionFactor;
                intensitySliderImage.minValue = intensitySliderImage.conversionFactorLast*(intensitySliderImage.globalCoordinateX[0]-intensitySliderImage.sliderMin)+intensitySliderImage.intercept;
                intensitySliderImage.maxValue = intensitySliderImage.conversionFactorLast*(intensitySliderImage.globalCoordinateX[1]-intensitySliderImage.sliderMin)+intensitySliderImage.intercept;

                intensitySliderImage.conversionFactor = (double)((intensitySliderImage.maxValue-intensitySliderImage.minValue)/(intensitySliderImage.sliderMax-intensitySliderImage.sliderMin));
                intensitySliderImage.globalCoordinateX[0] = intensitySliderImage.sliderMin;
                intensitySliderImage.globalCoordinateX[1] = intensitySliderImage.sliderMax;

              }

              //intensitySliderImage.findFeature(x-widgets[selectedWidget].x,y-widgets[selectedWidget].y);
              break;
            }
          } */
        /*if( (x>(globalCoordinateX[0]+16))&&(x<(globalCoordinateX[1]-4)) &&( (y>=41) && (y<=45) ) ){
                double conversionFactorLast = conversionFactor;
                double intercept2 = conversionFactorLast*(globalCoordinateX[0]-sliderMin)+intercept;
                double upplimit = conversionFactorLast*(globalCoordinateX[1]-sliderMin)+intercept;

                conversionFactor = (double)((upplimit-intercept2)/(sliderMax-sliderMin));
                intercept = intercept2;
                globalCoordinateX[0] = sliderMin;
                globalCoordinateX[1] = sliderMax;

                axisLabelsString = generateAxes(intercept2, upplimit, 7);

             //cout << '\r' << "FSAFSAFSAFSAF (" << x << ", " << y << ")" << flush;
             //showImage();
        } */
        //cout << '\r' << "LEFT MOUSE (" << x << ", " << y << ")" << flush;
        break;

    case EVENT_RBUTTONDBLCLK:
        //cout << '\r' << "Right button of the mouse is doublclicked - position (" << x << ", " << y << ")" << flush;
        /* if( (x>(globalCoordinateX[0]+16))&&(x<(globalCoordinateX[1]-4)) &&( (y>=41) && (y<=45) ) ){

                double intercept2 = 0;
                conversionFactor = (double)((65536.0-intercept2)/(sliderMax-sliderMin));
                intercept = intercept2;
                globalCoordinateX[0] = sliderMin;
                globalCoordinateX[1] = sliderMax;

                axisLabelsString = generateAxes(intercept, 65536, 7);

             showImage();
        } */
        break;
    case EVENT_LBUTTONDOWN:
        /*findImg( x, y, img0,selectedImg);
        cout << '\r' << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << flush;
        */
        //cout << '\r' << "LEFT MOUSE (" << x << ", " << y << ")" << flush;
        findImg( x, y, widgets);
        break;
    case EVENT_LBUTTONUP:
    if(selectedWidget!=-1){
			releaseImg(x,y);
      widgets[selectedWidget].value[0] = intensitySliderImage.lowerLimit;
      widgets[selectedWidget].value[1] = intensitySliderImage.upperLimit;

      Segmentation *pd = static_cast<Segmentation *>(userdata);
      if(selectedWidget==1){
        pd->minThresh = widgets[selectedWidget].value[0];
        pd->maxThresh = widgets[selectedWidget].value[1];
      }
      if(selectedWidget==2){
        pd->minArea = widgets[selectedWidget].value[0];
        pd->maxArea = widgets[selectedWidget].value[1];
      }
      if(selectedWidget==3){
        pd->eccentricityThresh = widgets[selectedWidget].value[1];
      }
      if(selectedWidget==4){
        pd->Min = widgets[selectedWidget].value[0];
        pd->Max = widgets[selectedWidget].value[1];
      }
      (*pd).runthreshold();
      imshow("display", pd->dsp);
      if(zoomOn){
        Mat inset = pd->dspHighRes( Rect(insetX0, insetY0, insetWidth, insetHeight) );
        imshow("zoom", inset);
      }
      R_FlushConsole();
      R_ProcessEvents();


      
      /*pd.Min = widgets[selectedWidget].value[0];
      //pd.maxThresh = intensitySliderImage.upperLimit;
      pd.Max = widgets[selectedWidget].value[1];
      pd.runthreshold(); */

      /*pd.runthreshold();
      cout << '\n' << "Min" << pd.Min << " Max:" << pd.Max << endl;
      imshow("display", pd->dsp); */

			selectedWidget = -1;
			intensitySliderImage.feature = -1;
		}
		//cout << '\n' << "Left mouse button is released - position" << x << ", " << y << ")" << endl;

		break;

    case EVENT_RBUTTONDOWN:
        //cout << '\r' << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << flush;
        break;
    case EVENT_MBUTTONDOWN:
        //cout << '\r' << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" << flush;
        break;
    case EVENT_MOUSEMOVE:
        //cout << '\r' << "Mouse move over the window - position (" << x << ", " << y << ")" << flush;

		if(selectedWidget!=-1){
                releaseImg(x,y);
			break;
		}

        break;
    }
}

void findImg(int x,int y, vector <Widget> widgets){

	for(unsigned int i=0;i<widgets.size();i++){
    //check which widget the user pressed
		if((x>=(widgets[i].x)) && (x<=(widgets[i].x+widgets[i].width ))&& (y>=(widgets[i].y ))&& (y<=(widgets[i].y+widgets[i].height ))){
			selectedWidget=i;
			intensitySliderImage.globalCoordinateX[0] = widgets[selectedWidget].guiPixelValue.at(0);
      intensitySliderImage.globalCoordinateX[1] = widgets[selectedWidget].guiPixelValue.at(1);
      intensitySliderImage.findFeature(x-widgets[selectedWidget].x,y-widgets[selectedWidget].y);
			break;
		}
	}

}

//set new cordiinates and redeaw the scene
void releaseImg(int x,int y){


    intensitySliderImage.updateCoordinates(x-widgets[selectedWidget].x,y-widgets[selectedWidget].y);
    intensitySliderImage.textlabel = widgets[selectedWidget].name;
    intensitySliderImage.maxValue = widgets[selectedWidget].balue[1];
    intensitySliderImage.updateAxes();
    intensitySliderImage.updateImage();
    widgets[selectedWidget].updateImage(intensitySliderImage.displayImage);

    widgets[selectedWidget].guiPixelValue.at(0) = intensitySliderImage.globalCoordinateX[0];
    widgets[selectedWidget].guiPixelValue.at(1) = intensitySliderImage.globalCoordinateX[1];

    /*
    if(point<3){
    //if point is not equal to three it means a specific point is being dragged.

    globalCoordinateX[point]=x;
    if( (globalCoordinateX[0]>=globalCoordinateX[1])|(globalCoordinateX[1]<=globalCoordinateX[0])  ){
        if(point==0){
        globalCoordinateX[0] = globalCoordinateX[1]-20;
        }else{
        globalCoordinateX[1] = globalCoordinateX[0]+20;
        }
    }else if(x>sliderMax){
        globalCoordinateX[point]=sliderMax;
    }
    else if(x<sliderMin){
        globalCoordinateX[point]=sliderMin;
    }

    }else if(point==3){
    //if point is equal too 3 then drag the the whole bar thing
    if( ((x-bardistanceLow)<sliderMin) ){
        //if slider bar is out or the lower bound then don't let it slide any further
        globalCoordinateX[0]=sliderMin;
        globalCoordinateX[1]=sliderMin+(bardistanceLow+bardistanceHigh);

    }else if( ((x+bardistanceHigh)>sliderMax) ){
        //if slider bar is out or the upper bound then don't let it slide any further
        globalCoordinateX[1]=sliderMax;
        globalCoordinateX[0]=sliderMax-(bardistanceLow+bardistanceHigh);
    }else{
        //normal sliude behavior puts the sliders at +/- from the point the user clicked.
        globalCoordinateX[0]=x-bardistanceLow;
        globalCoordinateX[1]=x+bardistanceHigh;
    }
    }else if(point==4){

    }
    */
	showImage();
}

//draw the scene components
void showImage(){

    /*

	Mat img3=img2.clone();

	//draw the points
	for(int j=0;j<nop;j++){
		overlayImage( &img3, &slider, Point( globalCoordinateX[j], 36 ) );
	}
    rectangle(
    img3,
    cv::Point(globalCoordinateX[0]+20, 41),
    cv::Point(globalCoordinateX[1]-4, 45),
    cv::Scalar(255, 255, 255),
    CV_FILLED
     );

    string lowlimit = NumberToString( (int)(conversionFactor*(globalCoordinateX[0]-sliderMin)+intercept)  );
    string upplimit = NumberToString( (int)(conversionFactor*(globalCoordinateX[1]-sliderMin)+intercept) );
    Point posLL(99-charLength*(lowlimit.size()), 48);
    Point posUL(720, 48);
    putText(img3, lowlimit, posLL, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);
    putText(img3, upplimit, posUL, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);

      putText(img3, axisLabelsString.at(0), pos1, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);
  putText(img3, axisLabelsString.at(1), pos2, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);
  putText(img3, axisLabelsString.at(2), pos3, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);
  putText(img3, axisLabelsString.at(3), pos4, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);
  putText(img3, axisLabelsString.at(4), pos5, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);
  putText(img3, axisLabelsString.at(5), pos6, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);
  putText(img3, axisLabelsString.at(6), pos7, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);

	imshow("controls", img3);
	img3.release();

	*/



	//for(int i=0; i < widgets.size(); i++){
        mergewithBackground( &ui_backgroundImage, &ui_display, &widgets[selectedWidget].storedImage, Point(widgets[selectedWidget].x, widgets[selectedWidget].y ) );
    //}

    imshow("controls", ui_display);

}


RcppExport SEXP imageshow(SEXP input, SEXP resizeP, SEXP filename, SEXP sliderFilename, SEXP backgroundFilename){
  double resizeParam = Rcpp::as<int>(resizeP);
  resizeParam = resizeParam/100;

  Rcpp::CharacterVector f(backgroundFilename);
  std::string ff(f[0]);
  ui_backgroundImage = imread(ff, IMREAD_UNCHANGED);

  ui_display = ui_backgroundImage.clone();

   Rcpp::CharacterVector sf(filename);
  std::string sff(sf[0]);
  //slider = imread(sff, IMREAD_UNCHANGED);

  Rcpp::CharacterVector sliderfile(sliderFilename);
  std::string sliderfilename(sliderfile[0]);

  imshowSliderImage.backgroundImage = imread(sff, IMREAD_UNCHANGED);
  imshowSliderImage.slider = imread(sliderfilename, IMREAD_UNCHANGED);
  imshowSliderImage.updateImage();

  std::string labels[] = {"8-bit render"};
  int rangeValues[] = {65536}; 
  for(unsigned int i=0; i < imshowwidget.size(); i++){
        //widgets[i].name.push_back(labels[i]);
        imshowwidget[i].name = labels[i];
        imshowwidget[i].balue.push_back(0);
        imshowwidget[i].balue.push_back(rangeValues[i]);

        imshowwidget[i].assignImage(imshowSliderImage.displayImage);
        imshowwidget[i].assignPos(0, ui_backgroundImage.rows -  (1*i + 1)*imshowwidget[i].height);
        //assign values
        if( imshowwidget[i].guiPixelValue.size()==0 ){
        imshowwidget[i].value.push_back(0);
        imshowwidget[i].value.push_back(rangeValues[i]);
        imshowwidget[i].guiPixelValue.push_back(110);
        imshowwidget[i].guiPixelValue.push_back(692);
        }
  }

    //DRAW THE GUI
  for(unsigned int i=0; i < imshowwidget.size(); i++){
    imshowSliderImage.globalCoordinateX[0] =  imshowwidget[i].guiPixelValue[0];
    imshowSliderImage.globalCoordinateX[1] =  imshowwidget[i].guiPixelValue[1];
    imshowSliderImage.textlabel = widgets[i].name;
    imshowSliderImage.maxValue = rangeValues[i];
    imshowSliderImage.updateAxes();
    imshowSliderImage.updateImage();
    imshowwidget[i].updateImage(imshowSliderImage.displayImage);

    mergewithBackground( &ui_backgroundImage, &ui_display, &imshowwidget[i].storedImage, Point(imshowwidget[i].x, imshowwidget[i].y ) );
  }

  namedWindow("controls", CV_WINDOW_AUTOSIZE);
  imshow("controls", ui_display);
  moveWindow("controls", 600, 200);

  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);
  Rcpp::Rcout << "Loading image:" << ffname << std::endl;
  ImageShow dispayimage;
  dispayimage.src = imread(ffname, -1); // -1 tag means "load as is"
  Rcpp::Rcout << "LOADED." << std::endl;
  Rcpp::Rcout << "Image type: " <<  getImgTypes(dispayimage.src.type()) << "_" << dispayimage.src.type()  << std::endl;

  int depth;
  if(dispayimage.src.type()==0){
    dispayimage.imgdepth = 255;  
    depth = 255;
  }else if(dispayimage.src.type()==2){
    dispayimage.imgdepth = 1000; //
    depth = 65535;
    dispayimage.src.convertTo(dispayimage.src, CV_16S);
    Rcpp::Rcout << "Changed image type to: " <<  getImgTypes(dispayimage.src.type()) << "_" << dispayimage.src.type() << std::endl;
  }
 

  dispayimage.Min = 0;
  dispayimage.Max = 65535;

    Rcpp::Rcout << "Resizing to: " <<  resizeParam*100 << "% of original size." << std::endl;
    resize(dispayimage.src, dispayimage.src, Size(), resizeParam, resizeParam, INTER_LINEAR);

  dispayimage.updateImage();
  imshow("display", dispayimage.dsp);
  moveWindow("display", 100, 300);

   //set the callback function for any mouse event
  setMouseCallback("controls", CallBackFunc, &dispayimage); // NULL
  setMouseCallback("display", CallBackZoom, &dispayimage); // NULL
  //INITIALIZE LOOP

  int k;
  while(1){


  k=waitKey(0);
  //cout << '\n' << "KEY PRESSED: " << k << endl;
  if(k == 104){
    dispayimage.updateImage();
    imshow("display", dispayimage.dsp);
    if(zoomOn){
      Mat inset = dispayimage.dspHighRes( Rect(insetX0, insetY0, insetWidth, insetHeight) );
      imshow("zoom", inset);
    }
  }
  if(k == 122){
    zoomOn = false;
    destroyWindow("zoom");
  }
  if( (k == 27)|(k == -1)|(k == 115) ){
  destroyWindow("controls");
  destroyWindow("display");
  destroyWindow("zoom");
  break;
  }
  }


}


RcppExport SEXP GUI(SEXP input, SEXP numthresh, SEXP resizeP, SEXP filename, SEXP sliderFilename, SEXP backgroundFilename) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; 

  //to handle execution time logging
clock_t tStart, tStop;
  double resizeParam = Rcpp::as<int>(resizeP);
  resizeParam = resizeParam/100;

  Rcpp::CharacterVector f(backgroundFilename);
  std::string ff(f[0]);
  ui_backgroundImage = imread(ff, IMREAD_UNCHANGED);

  ui_display = ui_backgroundImage.clone();

   Rcpp::CharacterVector sf(filename);
  std::string sff(sf[0]);
  //slider = imread(sff, IMREAD_UNCHANGED);

  Rcpp::CharacterVector sliderfile(sliderFilename);
  std::string sliderfilename(sliderfile[0]);

  intensitySliderImage.backgroundImage = imread(sff, IMREAD_UNCHANGED);
  intensitySliderImage.slider = imread(sliderfilename, IMREAD_UNCHANGED);
  intensitySliderImage.updateImage();

  std::string labels[] = {"brain outline", "intensity", "soma area", "eccentricity", "8-bit render"};
  int rangeValues[] = {65536, 65536, 1000, 1000, 65536}; 
  for(unsigned int i=0; i < widgets.size(); i++){
        //widgets[i].name.push_back(labels[i]);
        widgets[i].name = labels[i];
        widgets[i].balue.push_back(0);
        widgets[i].balue.push_back(rangeValues[i]);

        widgets[i].assignImage(intensitySliderImage.displayImage);
        widgets[i].assignPos(0, ui_backgroundImage.rows -  (1*i + 1)*widgets[i].height);
        //assign values
        if( widgets[i].guiPixelValue.size()==0 ){
        widgets[i].value.push_back(0);
        widgets[i].value.push_back(rangeValues[i]);
        widgets[i].guiPixelValue.push_back(110);
        widgets[i].guiPixelValue.push_back(692);
        }
  }

  // USE THIS LATE RON IN DOUBLE SLIDER
  /* string text = "intensity";
  Point textOrg(603, 22);
  putText(img2, text, textOrg, fontFace, fontScale, Scalar::all(255), thickness, CV_AA);
  */


  //widgets[3].updateImage()

  //DRAW THE GUI
  for(unsigned int i=0; i < widgets.size(); i++){
    intensitySliderImage.globalCoordinateX[0] =  widgets[i].guiPixelValue[0];
    intensitySliderImage.globalCoordinateX[1] =  widgets[i].guiPixelValue[1];
    intensitySliderImage.textlabel = widgets[i].name;
    intensitySliderImage.maxValue = rangeValues[i];
    intensitySliderImage.updateAxes();
    intensitySliderImage.updateImage();
    widgets[i].updateImage(intensitySliderImage.displayImage);

    mergewithBackground( &ui_backgroundImage, &ui_display, &widgets[i].storedImage, Point(widgets[i].x, widgets[i].y ) );
  }

  namedWindow("controls", CV_WINDOW_AUTOSIZE);
  imshow("controls", ui_display);
  moveWindow("controls", 600, 200);
  //showImage();

  


  //load the image 
  /*
  Rcpp::CharacterVector mf(input);
  std::string mff(mf[0]);
  Rcpp::Rcout << "Loading image:" << mff << std::endl;
  t0 = clock();
  pd.src = imread(mff, -1); // -1 tag means "load as is"
  t1 = clock();
  double time_elapsed = (double)((t1 - t0) / CLOCKS_PER_SEC);
  Rcpp::Rcout << "LOADED." << " loading took " <<  time_elapsed << " seconds." << std::endl;
  Rcpp::Rcout << "Image type: " <<  getImgType(pd.src.type()) << "__" << pd.src.type()  << std::endl;
  */


  Rcpp::CharacterVector fname(input);
  std::string ffname(fname[0]);
  Rcpp::Rcout << "Loading image:" << ffname << std::endl;
    Segmentation pd;
  tStart = clock();
  pd.src = imread(ffname, -1); // -1 tag means "load as is"
  tStop = clock();
  double time_elapsed = (double)((tStop - tStart) / CLOCKS_PER_SEC);
  Rcpp::Rcout << "LOADED." << " loading took " <<  time_elapsed << " seconds." << std::endl;
  Rcpp::Rcout << "Image type: " <<  getImgTypes(pd.src.type()) << "_" << pd.src.type()  << std::endl;

  int depth;
  if(pd.src.type()==0){
    pd.imgdepth = 255;  
    depth = 255;
  }else if(pd.src.type()==2){
    pd.imgdepth = 1000; //
    depth = 65535;
    pd.src.convertTo(pd.src, CV_16S);
    Rcpp::Rcout << "Changed image type to: " <<  getImgTypes(pd.src.type()) << "_" << pd.src.type() << std::endl;
  }
 
  pd.minThresh = 0;
  pd.maxThresh = pd.imgdepth;
  pd.displayImage = true;
  pd.minArea = 0;
  pd.maxArea = 1000;
  pd.numThresholds = Rcpp::as<int>(numthresh);
    pd.eccentricityThresh = 1000;
  pd.Min = 0;
  pd.Max = 65535;
  pd.minThresh = 0;
  pd.maxThresh = 65535;
  pd.endSegment = false;


    Rcpp::Rcout << "Resizing to: " <<  resizeParam*100 << "% of original size." << std::endl;
    resize(pd.src, pd.src, Size(), resizeParam, resizeParam, INTER_LINEAR);

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
    
    Mat histImage( hist_h, hist_w, CV_8UC3, Scalar( 0,0,0) );
    normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
    
    for( int i = 1; i < histSize; i++ )
    {
        line( histImage, Point( bin_w*(i-1), hist_h - cvRound(hist.at<float>(i-1)) ) ,
             Point( bin_w*(i), hist_h - cvRound(hist.at<float>(i)) ),
             Scalar( 100, 163, 117), 2, 8, 0  );
    }
    
    namedWindow( "histogram", 1 );    imshow( "histogram", histImage );

  pd.runthreshold();
  imshow("display", pd.dsp);
  moveWindow("display", 100, 300);

   //set the callback function for any mouse event
  setMouseCallback("controls", CallBackFunc, &pd); // NULL
  setMouseCallback("display", CallBackZoom, &pd); // NULL
  //INITIALIZE LOOP

  int k;
  while(1){


  k=waitKey(0);
  //cout << '\n' << "KEY PRESSED: " << k << endl;
  if(k == 104){
    if(pd.hideFilter){
      pd.hideFilter = 0;
    }else{
      pd.hideFilter = 1;
    }
    pd.runthreshold();
    imshow("display", pd.dsp);
    if(zoomOn){
      Mat inset = pd.dspHighRes( Rect(insetX0, insetY0, insetWidth, insetHeight) );
      imshow("zoom", inset);
    }
  }
  if(k == 122){
    zoomOn = false;
    destroyWindow("zoom");
  }
  if( (k == 27)|(k == -1)|(k == 115) ){
    cout << '\n' << "Assembling output list" << endl;
  pd.endSegment = true;
  pd.runthreshold();
  cout << '\n' << "OUTPUT SEGMENTED CELLS: " << pd.centroidX.size() << endl;
  destroyWindow("controls");
  destroyWindow("display");
  destroyWindow("histogram");
  destroyWindow("zoom");
  break;
  }
  }
    //return R_NilValue;



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
    _["Max"] = pd.Max,
    _["Min"] = pd.Min,
    _["x"] = pd.centroidX,
    _["y"] = pd.centroidY,
    _["intensity"] = pd.intensitySoma,
    _["soma.area"] = pd.contourSomaArea
  );

END_RCPP  
}
