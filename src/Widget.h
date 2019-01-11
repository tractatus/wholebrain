#ifndef WIDGET_H
#define WIDGET_H

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>

//for sleep
#include <unistd.h>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>


class Widget
{
    public:
        Widget();
        ~Widget();

        cv::Mat storedImage;
        int x;
        int y;
        int width;
        int height;


        std::string name;
        std::vector<int> guiPixelValue;
        std::vector<int> value;
        std::vector<int> balue;
        int range[];

        void assignImage(cv::Mat &src);
        void assignPos(int xPos, int yPos);
        void updateImage(cv::Mat &src);
    protected:
    private:
};

#endif // WIDGET_H
