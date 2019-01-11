#include "DoubleSlider.h"

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>

//for sleep
#include <iostream>
#include <unistd.h>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

using namespace cv;
using namespace std;

DoubleSlider::DoubleSlider()
{
    sliderMin = 110; // lower bound in pixels of the slider
    sliderMax = 692; // upper bound in pixels of the slider
    bardistanceLow = 0; //distance from the lower slider to the feature where the suer decided to drag the slider bar
    bardistanceHigh = 0; //distance from the upper slider to the feature where the suer decided to drag the slider bar
    const static int newglobalCoordinateX[] = {sliderMin,sliderMax};
    copy(newglobalCoordinateX, newglobalCoordinateX+sizeof(newglobalCoordinateX)/sizeof(newglobalCoordinateX[0]), globalCoordinateX);//x coordinates of the sliders, default if min and max
    const static int newglobalCoordinateY[] = {36,36};
    copy(newglobalCoordinateY, newglobalCoordinateY+sizeof(newglobalCoordinateY)/sizeof(newglobalCoordinateY[0]), globalCoordinateY);
    feature=-1;			//currently selected feature
    nop=2;
    intercept = 0;
    maxValue = 65536;
    minValue = 0;
    conversionFactor = (double)((maxValue-intercept)/(sliderMax-sliderMin));

    Mat slider;

    axisLabelsString = generateAxes(minValue, maxValue, 7);

    fontFace = FONT_HERSHEY_SIMPLEX;
    fontScale = 0.5;
    thickness = 1;
    linePos = 77;
    charLength = 10;
    pos1 = Point(118-charLength*(axisLabelsString.at(0).size()/2), linePos);
    pos2 = Point(215-charLength*(axisLabelsString.at(1).size()/2), linePos);
    pos3 = Point(312-charLength*(axisLabelsString.at(2).size()/2), linePos);
    pos4 = Point(409-charLength*(axisLabelsString.at(3).size()/2), linePos);
    pos5 = Point(506-charLength*(axisLabelsString.at(4).size()/2), linePos);
    pos6 = Point(603-charLength*(axisLabelsString.at(5).size()/2), linePos);
    pos7 = Point(700-charLength*(axisLabelsString.at(6).size()/2), linePos);

    textlabel = "brain";
}

DoubleSlider::~DoubleSlider()
{
    //dtor
}

template <typename T>
  string DoubleSlider :: NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

vector<double> DoubleSlider :: linspace(double min, double max, int n)
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

vector<string>  DoubleSlider::generateAxes(double min, double max, int n){
      vector<double> axisLabelsDouble =  linspace(min, max, n);
      vector<int> axisLabelsInteger(axisLabelsDouble.begin(), axisLabelsDouble.end());

      vector<string> axisLabelsString;

      for(vector<int>::size_type i = 0; i != axisLabelsInteger.size(); i++) {
        axisLabelsString.push_back(NumberToString(axisLabelsInteger[i]));
      }
      return(axisLabelsString);
}

void DoubleSlider :: overlayImage(Mat* src, Mat* overlay, const Point& location)
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
                src->data[y * src->step + src->channels() * x + c] = srcPx * (1. - opacity) + overlayPx * opacity;
            }
        }
    }
}


//a function in order to find exactly which feature in the scene the user interacted with
void DoubleSlider :: findFeature(int x,int y)
{

	for(int i=0;i<nop;i++){
		if((x>=(globalCoordinateX[i])) && (x<=(globalCoordinateX[i]+16 ))&& (y>=(globalCoordinateY[i] ))&& (y<=(globalCoordinateY[i]+15 ))){
			feature=i;
			break;
		}
	}

	  if( (x>(globalCoordinateX[0]+16))&&(x<(globalCoordinateX[1]-4)) &&( (y>=41) && (y<=55) ) ){
        feature = 3;
        bardistanceLow = x-globalCoordinateX[0];
        bardistanceHigh = globalCoordinateX[1]-x;
      }

      if(feature!=-1){this->updateCoordinates(x,y);}

}

//a function to update the coordinates of the feature the user is interacting with
void DoubleSlider :: updateCoordinates(int x,int y)
{ // Mat &backgroundImage, Mat &displayImage
        if(feature<3){
    //if feature is not equal to three it means a specific feature is being dragged.

    globalCoordinateX[feature]=x;
    if( (globalCoordinateX[0]>=globalCoordinateX[1])|(globalCoordinateX[1]<=globalCoordinateX[0])  ){
        if(feature==0){
        globalCoordinateX[0] = globalCoordinateX[1]-20;
        }else{
        globalCoordinateX[1] = globalCoordinateX[0]+20;
        }
    }else if(x>sliderMax){
        globalCoordinateX[feature]=sliderMax;
    }
    else if(x<sliderMin){
        globalCoordinateX[feature]=sliderMin;
    }

    }else if(feature==3){
    //if feature is equal too 3 then drag the the whole bar thing
    if( ((x-bardistanceLow)<sliderMin) ){
        //if slider bar is out or the lower bound then don't let it slide any further
        globalCoordinateX[0]=sliderMin;
        globalCoordinateX[1]=sliderMin+(bardistanceLow+bardistanceHigh);

    }else if( ((x+bardistanceHigh)>sliderMax) ){
        //if slider bar is out or the upper bound then don't let it slide any further
        globalCoordinateX[1]=sliderMax;
        globalCoordinateX[0]=sliderMax-(bardistanceLow+bardistanceHigh);
    }else{
        //normal sliude behavior puts the sliders at +/- from the feature the user clicked.
        globalCoordinateX[0]=x-bardistanceLow;
        globalCoordinateX[1]=x+bardistanceHigh;
    }
    }else if(feature==4){

    }

	this->updateImage(); //this->updateImage(backgroundImage, displayImage);

}

void DoubleSlider :: updateAxes(){
  conversionFactor = (double)((maxValue-intercept)/(sliderMax-sliderMin));
  axisLabelsString = generateAxes(intercept, maxValue, 7);
  pos1 = Point(118-charLength*(axisLabelsString.at(0).size()/2), linePos);
    pos2 = Point(215-charLength*(axisLabelsString.at(1).size()/2), linePos);
    pos3 = Point(312-charLength*(axisLabelsString.at(2).size()/2), linePos);
    pos4 = Point(409-charLength*(axisLabelsString.at(3).size()/2), linePos);
    pos5 = Point(506-charLength*(axisLabelsString.at(4).size()/2), linePos);
    pos6 = Point(603-charLength*(axisLabelsString.at(5).size()/2), linePos);
    pos7 = Point(700-charLength*(axisLabelsString.at(6).size()/2), linePos);
}

void DoubleSlider :: updateImage(){ //intensitySliderImage.backgroundImage
  conversionFactor = (double)((maxValue-intercept)/(sliderMax-sliderMin));
	Mat img3=backgroundImage.clone();

	//draw the features
	for(int j=0;j<nop;j++){
		overlayImage( &img3, &slider, Point( globalCoordinateX[j], 36 ) );
	}

	//draw bar
    rectangle(
    img3,
    cv::Point(globalCoordinateX[0]+20, 41),
    cv::Point(globalCoordinateX[1]-4, 45),
    cv::Scalar(255, 255, 255),
    cv::FILLED
     );
    lowerLimit = (int)(conversionFactor*(globalCoordinateX[0]-sliderMin)+intercept);
    upperLimit = (int)(conversionFactor*(globalCoordinateX[1]-sliderMin)+intercept);
    string lowlimit = NumberToString(lowerLimit  );
    string upplimit = NumberToString( upperLimit );
    Point posLL(99-charLength*(lowlimit.size()), 48);
    Point posUL(720, 48);
    putText(img3, lowlimit, posLL, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);
    putText(img3, upplimit, posUL, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);

      putText(img3, axisLabelsString.at(0), pos1, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);
      putText(img3, axisLabelsString.at(1), pos2, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);
      putText(img3, axisLabelsString.at(2), pos3, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);
      putText(img3, axisLabelsString.at(3), pos4, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);
      putText(img3, axisLabelsString.at(4), pos5, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);
      putText(img3, axisLabelsString.at(5), pos6, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);
      putText(img3, axisLabelsString.at(6), pos7, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);


     Point textOrg(598, 22);
     putText(img3, textlabel, textOrg, fontFace, fontScale, Scalar::all(255), thickness, cv::LINE_AA);

	displayImage = img3.clone(); //imshow("controls", img3);
	img3.release();
}
