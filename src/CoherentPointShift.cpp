
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

//for sleep
#include <unistd.h>
//Rcpp
#include <Rcpp.h>

using namespace cv;
using namespace Rcpp;

#include <cpd/nonrigid.hpp>
//#include <Eigen/Dense>

RcppExport SEXP CoherentPointDriftRegistration(SEXP input, SEXP srcX, SEXP srcY, SEXP dstX, SEXP dstY, SEXP resizeP, SEXP MaxDisp, SEXP MinDisp, SEXP outputfile){
    BEGIN_RCPP

   Rcpp::RNGScope __rngScope;

    Rcpp::CharacterVector fname(input);
    std::string ffname(fname[0]);
    Rcpp::Rcout << "Loading image:" << ffname << std::endl;
    Mat img = imread(ffname, -1);
    
    //original image size
    int origwidth = img.cols; 
    int origheight = img.rows;

    double resizeParam = Rcpp::as<double>(resizeP);

    double Max = Rcpp::as<int>(MaxDisp);
    double Min = Rcpp::as<int>(MinDisp);

  	resize(img, img, Size(), resizeParam, resizeParam, INTER_LINEAR);
  	//resized used for registration
    int width = img.cols; // N
    int height = img.rows; // M

  	Rcpp::CharacterVector of(outputfile);
    std::string off(of[0]);


	Mat displayoutput;
	Mat img8bit;
	Mat displayoutputDST;
	Mat img8bitDST;
	double min, max;
	if(Max>0){
		min = Min;
		max = Max;
	}else{
		minMaxLoc(img, &min, &max);
	}
	img.convertTo(img8bit,CV_8UC1,255.0/(max-min) );
	normalize(img8bit, displayoutput, 0, 255, NORM_MINMAX, CV_8UC1);

	std::string filepath;
    filepath = off + "_undistorted.png";
	imwrite(filepath, displayoutput);
	filepath = off + "_undistorted.tif";
    imwrite(filepath, img);

    //define meshgrid in y-coordinates could be done more elegant with Eigen but quick and dirty for now
    cpd::Matrix y(height*width, 1);
    cpd::Matrix x(y.rows(),1);
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            y(j+height*i, 0) = j;
            x(j+height*i, 0) = i;
        }
    }

    //column bind to a 2D matrix
    cpd::Matrix grid(y.rows(), x.cols()+y.cols());
	grid << x, y;


	//reading in the correspondence points
	std::vector<int> srX = as<std::vector<int> >(srcX);
	std::vector<int> srY = as<std::vector<int> >(srcY); 
	std::vector<int> dtX = as<std::vector<int> >(dstX); 
 	std::vector<int> dtY = as<std::vector<int> >(dstY); 


	cpd::Matrix iP(srX.size(), 2),  iiP(srX.size(), 2);
	//iP << srX, srY;
	//iiP << dtX, dtY;
	for (unsigned i=0; i < srX.size(); i++) {
	   iP(i, 0) = srX[i];
	   iP(i, 1) = srY[i];
	   iiP(i, 0) = dtX[i];
	   iiP(i, 1) = dtY[i];   
	}



	//initialize the runner
  	cpd::Nonrigid runner;
  	//set the parameters
    runner.outliers(0.7).sigma2(0.0).tolerance(1e-5).max_iterations(150);
    cpd::NonrigidResult result = runner.run(iP, iiP);

  	//Rcpp::Rcout << "Results:\n" << result.points << std::endl;
  	//Rcpp::Rcout << "Reference:\n" << m_fish << std::endl;
  	Rcpp::Rcout << "Iterations: " << result.iterations << std::endl;

  	/* cpd::GaussTransformDirect direct;
    cpd::Probabilities probabilities = direct.compute(iP, result.points, result.sigma2, 0.1);

    Rcpp::Rcout << "Correspondance:\n" << probabilities.correspondence << std::endl;
   	*/

   	Rcpp::Rcout << "\nComputing Coherent Point Shift (CPD) transformation\n" << std::endl;
   	cpd::Matrix transform = result.transformation_grid(iiP, grid);

   	//from the eigen vector to the std vector
    std::vector<double> trans(transform.data(), transform.data() + transform.rows() * transform.cols());
    //std::vector<double> correspondance(probabilities.correspondence.data(), probabilities.correspondence.data() + probabilities.correspondence.rows() * probabilities.correspondence.cols());

  	return List::create(
    	_["trans"] = trans,
    	_["width"] = origwidth,
    	_["height"] = origheight,
       	_["dwp.width"] = width,
    	_["dwp.height"] = height	
  	);

    END_RCPP
    
}