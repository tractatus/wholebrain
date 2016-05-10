#include "Widget.h"

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

Widget::Widget()
{
     name = "";

     x = 0;
     y  = 0;
     width = 0;
     height  = 0;
}

Widget::~Widget()
{
    //dtor
}

void Widget::assignImage(cv::Mat &src)
{
    storedImage = src;
    width = storedImage.cols;
    height  = storedImage.rows;
}

void Widget::assignPos(int xPos, int yPos)
{
    x = xPos;
    y  = yPos;
}

void Widget::updateImage(cv::Mat &src)
{
    storedImage = src.clone();
}
