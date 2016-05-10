#include <Rcpp.h>
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;


/* apply operation to stack */
RcppExport SEXP getContour(SEXP input, SEXP outputfolder, SEXP outname, SEXP kernel, SEXP display,  SEXP outputfile, SEXP writetoconsole) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector f(input);
  int num_files = f.size();

  Rcpp::CharacterVector outputfilenames(outname);

  Rcpp::CharacterVector outf(outputfolder);
  std::string outfolder(outf[0]);

  int KERNEL = Rcpp::as<int>(kernel);
  int DISPLAY = Rcpp::as<int>(display);
  Rcpp::CharacterVector of(outputfile);
  std::string off(of[0]);

  bool verbose = Rcpp::as<bool>(writetoconsole);

  int barWidth = 70;
  float progress = 0.0;

  std::string ff(f[0]);
  Mat remove = imread(ff, -1);
  int rows = remove.rows;
  int cols = remove.cols;
  vector<Mat> zpositions(num_files,Mat(cols,rows,remove.type()));
  vector<Mat> destination(num_files,Mat(cols,rows,remove.type()));
  Mat src;
  
  Mat dst = Mat::zeros(remove.size(), remove.type());
  if(verbose){Rcpp::Rcout << "Loading " <<  num_files << " stacks into RAM." << std::endl;}
  for (int i = 0; i < f.size(); i++){
    std::string filename(f[i]);

    src = imread(filename, -1);
    zpositions.at(i) = src;

    


    //print progress bar in console
    if(verbose){Rcpp::Rcout << "  [";
    int pos = barWidth * progress;
    for (int j = 0; j < barWidth; ++j) {
        if (j < pos) Rcpp::Rcout << "=";
        else if (j == pos) Rcpp::Rcout << ">";
        else Rcpp::Rcout << " ";
    }
    Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::cout.flush();
    R_FlushConsole();
    R_ProcessEvents();
    progress += (float)1/(num_files-1);
  }

  }
  if(verbose){Rcpp::Rcout << "\n" << std::endl;
  Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;

  Rcpp::Rcout << "Running command. " << "generate gain image" << " on stack" << std::endl;}
  getGain(num_files, cols, rows, zpositions, dst);
    // Gaussian smoothing
  GaussianBlur( dst, dst, Size( KERNEL, KERNEL ), 0, 0 );
      imwrite(off, dst);
  if(verbose){  Rcpp::Rcout << "\n" << std::endl;
    Rcpp::Rcout << "Gain image saved as: " << off << std::endl;
  Rcpp::Rcout << "====== GAIN IMAGE GENERATED ======" << std::endl;

  Rcpp::Rcout << "Running command. " << "illumination correction" << " on stack" << std::endl;}

  Mat output;

  barWidth = 70;
  progress = 0.0;
 for(size_t s = 0; s != num_files; ++s){
  Scalar meanGain = mean(dst);
  divide(zpositions.at(s), dst, output, meanGain.val[0]);
  
  //multiply(output, meanGain, output);
  string filename(outputfilenames[s]);
  string outputname;
  outputname =  outfolder + "/" + "FFC_" + filename;
  imwrite(outputname, output);

  //print progress bar in console
  if(verbose){  Rcpp::Rcout << "  [";
    int pos = barWidth * progress;
    for (int j = 0; j < barWidth; ++j) {
        if (j < pos) Rcpp::Rcout << "=";
        else if (j == pos) Rcpp::Rcout << ">";
        else Rcpp::Rcout << " ";
    }
    Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::cout.flush();
    R_FlushConsole();
    R_ProcessEvents();
    progress += (float)1/(num_files-1); }


 }
 if(verbose){Rcpp::Rcout << "\n" << std::endl;
 Rcpp::Rcout << "====== OUTPUT SAVED ======" << std::endl;}
   /*
  END
  */
END_RCPP  
}



