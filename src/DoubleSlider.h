#ifndef DOUBLESLIDER_H
#define DOUBLESLIDER_H

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>

//for sleep
#include <unistd.h>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

class DoubleSlider
{
    public:
        DoubleSlider();
        virtual ~DoubleSlider();

        template <typename T>
        std::string NumberToString ( T Number );
        //functions for finding, updating and displaying features
        void findFeature(int x,int y);
        void updateCoordinates(int x,int y); //, cv::Mat &backgroundImage, cv::Mat &displayImage
        void updateAxes();
        void updateImage();

        std::vector<double> linspace(double min, double max, int n);
        void overlayImage(cv::Mat* src, cv::Mat* overlay, const cv::Point& location);

        cv::Mat backgroundImage;
        cv::Mat slider;
        cv::Mat displayImage;

        std::string textlabel;

        //for exeses generated
        std::vector<std::string> generateAxes(double min, double max, int n); //function for generating axes labels

        std::vector<std::string> axisLabelsString;

        int fontFace;
        double fontScale;
        int thickness;
        int linePos;
        int charLength;
        cv::Point pos1;
        cv::Point pos2;
        cv::Point pos3;
        cv::Point pos4;
        cv::Point pos5;
        cv::Point pos6;
        cv::Point pos7;

        int sliderMin; // lower bound in pixels of the slider
        int sliderMax; // upper bound in pixels of the slider
        int lowerLimit;
        int upperLimit;
        int bardistanceLow; //distance from the lower slider to the point where the suer decided to drag the slider bar
        int bardistanceHigh; //distance from the upper slider to the point where the suer decided to drag the slider bar
        int globalCoordinateX[2];	//x coordinates of the sliders, default if min and max
        int globalCoordinateY[2]; //y coordinates of the same
        int feature;			//currently selected feature
        int nop; //number of slider points (usually 2)
        double intercept;
        double conversionFactor;
        double conversionFactorLast;
        int maxValue; // max value for the slider 65536 for intensity based siders
        int minValue; // max value for the slider 65536 for intensity based siders

    protected:
    private:



};

#endif // DOUBLESLIDER_H
