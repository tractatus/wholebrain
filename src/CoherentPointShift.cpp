
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

RcppExport SEXP CoherentPointDriftRegistration(){
    BEGIN_RCPP
    	// SEXP input, SEXP srcX, SEXP srcY, SEXP dstX, SEXP dstY, SEXP resizeP, SEXP MaxDisp, SEXP MinDisp, SEXP outputfile

   Rcpp::RNGScope __rngScope;
   /*
    Rcpp::CharacterVector fname(input);
    std::string ffname(fname[0]);
    Rcpp::Rcout << "Loading image:" << ffname << std::endl;
    Mat img = imread(ffname, -1);
    
    //original image size
    int width = img.cols; // N
    int height = img.rows; // M

    double resizeParam = Rcpp::as<double>(resizeP);

    double Max = Rcpp::as<int>(MaxDisp);
    double Min = Rcpp::as<int>(MinDisp);

  	resize(img, img, Size(), resizeParam, resizeParam, INTER_LINEAR);


  	Rcpp::CharacterVector of(outputfile);
    std::string off(of[0]);
    */

    int width = 50;
    int height = 20;

    //define meshgrid in y-coordinates
    cpd::Matrix y(height, 1);
    y.col(0) = cpd::Array::LinSpaced(height, 0, height-1);
   	y = y.replicate(width,1);

   	Rcpp::Rcout << "Y:\n" << y.rows() << std::endl;


   	//define meshgrid in x-coordinates
   	cpd::Matrix x(y.rows(),1);
   	for (cpd::Matrix::Index i = 0; i < width; ++i) {
        for (cpd::Matrix::Index j = 0; j < height; ++j) {
            x(i * , 1) = i;
        }
    }
    /*

  	cpd::Matrix m_fish_distorted = cpd::matrix_from_path("/Users/danielfurth/Downloads/cpd-master/tests/fixtures/fish_distorted.csv");
  	cpd::Matrix m_fish = cpd::matrix_from_path("/Users/danielfurth/Downloads/cpd-master/tests/fixtures/fish.csv");

  	Rcpp::Rcout << "File:\n" << "/Users/danielfurth/Downloads/cpd-master/tests/fixtures/fish_distorted2.csv" << std::endl;



  	cpd::Nonrigid runner;
    runner.outliers(0.7).sigma2(0.0).tolerance(1e-5).max_iterations(150);
    cpd::NonrigidResult result = runner.run(m_fish, m_fish_distorted);

  	Rcpp::Rcout << "Results:\n" << result.points << std::endl;
  	Rcpp::Rcout << "Reference:\n" << m_fish << std::endl;
  	Rcpp::Rcout << "Iter:\n" << result.iterations << std::endl;

  	cpd::GaussTransformDirect direct;
    cpd::Probabilities probabilities = direct.compute(m_fish, result.points, result.sigma2, 0.1);

    Rcpp::Rcout << "Correspondance:\n" << probabilities.correspondence << std::endl;

   	Rcpp::Rcout << "Scale:\n" << result.scale << std::endl;
   	Rcpp::Rcout << "Translation:\n" << result.translation << std::endl;

   	Rcpp::Rcout << "m_W size:\n" << result.m_W.rows() << " | " << result.m_W.cols() << std::endl;

   	result.grid = m_fish.replicate(2,1);

   	Rcpp::Rcout << "Grid size:\n" << result.grid.rows() << " | " << result.grid.cols() << std::endl;


   	cpd::Matrix transform = result.matrix(m_fish_distorted);

   	Rcpp::Rcout << "T size:\n" << transform.rows() << " | " << transform.cols() << std::endl;


   	Rcpp::Rcout << "X:\n" << x.rows() << std::endl;
   	Rcpp::Rcout << "Y:\n" << y.rows() << std::endl;


    //cpd::Matrix transform = result.matrix();

    //cpd::Matrix fish = cpd::apply_transformation_matrix(m_fish_distorted, transform);

    //NumericMatrix xx(wrap(fish));
    //Rcpp::Rcout << "Results:\n" << fish << std::endl;

	//cpd::RigidResult result = cpd::rigid(fixed, moving);
	//cpd::RigidResult result = cpd::rigid(fixed, moving);

    /*cpd::Matrix fixed = load_points_from_somewhere();
    cpd::Matrix moving = load_points_from_somewhere();

    
    registration.correspondence(true).outliers(0.2);
    cpd::RegiResult result = registration.run(fixed, moving);
    */

	/*return List::create(
    	_["transform"] = xx
  	);*/
  	return R_NilValue;

    END_RCPP
    
}