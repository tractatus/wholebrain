
#include <Rcpp.h>
// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"


#include <math.h>       /* log2 */

using namespace cv;
using namespace std;
using namespace Rcpp;

void copySrcTile(const cv::Mat& src, cv::Mat& srcTile, cv::Rect &tile)
{
    //top left corner
    auto tl = tile.tl();
    //bottom left corner
    auto br = tile.br();
    
    cv::Point tloffset, broffset;
    
    //Take care of border cases
    if (tile.x < 0)
    {
      
        tloffset.x = -tile.x;
        tile.x = 0;
    }
    
    if (tile.y < 0)
    {
        tloffset.y = -tile.y;
        tile.y = 0;
    }
    
    if (br.x >= src.cols)
    {
      
        broffset.x = br.x - src.cols + 1;
        tile.width -= broffset.x; 
    }
    
    if (br.y >= src.rows)
    {

        broffset.y = br.y - src.rows + 1;
        tile.height -= broffset.y;
    }
    
    // If any of the tile sides exceed source image boundary we must use copyMakeBorder to make proper paddings for this side
    if (tloffset.x > 0 || tloffset.y > 0 || broffset.x > 0 || broffset.y > 0)
    {
        cv::Rect paddedTile(tile.tl(), tile.br());
        assert(paddedTile.x >= 0);
        assert(paddedTile.y >= 0);
        assert(paddedTile.br().x < src.cols);
        assert(paddedTile.br().y < src.rows);
        
        cv::copyMakeBorder(src(paddedTile), srcTile, tloffset.y, broffset.y, tloffset.x, broffset.x, BORDER_CONSTANT);
    }
    else
    {
        // Entire tile (with paddings lies inside image and it's safe to just take a region:
        src(tile).copyTo(srcTile);
    }
}



/**
* @brief generic algorithm for other channel types except of uchar
* @param input   the input image
* @param output  the output image
* @param smin    total number of minimum pixels
* @param smax    total number maximum pixels
* @param channel the channel used to compute the histogram
*
* This algorithm only support uchar channel and float channel by now
*/
void get_quantile(cv::Mat &input, size_t smin, size_t smax,  int &alpha, int &beta ) 
{

    std::vector<float> temp_input((float*)input.data, (float*)input.data + input.rows * input.cols);

    std::sort(std::begin(temp_input), std::end(temp_input));

    beta = temp_input[smin];
    alpha =  temp_input[smax];
}



RcppExport SEXP createWeb(SEXP input, SEXP alpha, SEXP beta, SEXP verbose, SEXP outputfile) {
    BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    //convert to Cpp int
    int alphaInt = Rcpp::as<int>(alpha);
    int betaADJ = Rcpp::as<int>(beta);

    int printOutput = Rcpp::as<int>(verbose);
    
    Rcpp::CharacterVector of(outputfile);
    std::string off(of[0]);

    
    Rcpp::CharacterVector f(input);
    std::string ff(f[0]);
    if(printOutput){
    Rcpp::Rcout << "Loading image:" << ff << std::endl;
    }
    Mat sourceImage = imread(ff, -1); //CV_LOAD_IMAGE_GRAYSCALE
    Mat finalImage = sourceImage.clone();
    if(printOutput){
    Rcpp::Rcout << "LOADED." << std::endl;
    }

    int tiers;
    int height = sourceImage.rows;
    int width = sourceImage.cols;
    int maxdim = (width > height ? width : height);
    tiers = (int)log2(maxdim/256);

    //Range correct the image
    /*double MaxRange;
    double MinRange;
    minMaxIdx(sourceImage, &MinRange, &MaxRange);

    sourceImage -= MinRange;
    sourceImage.convertTo(sourceImage,CV_8UC1,255.0/(MaxRange-MinRange));
    double alphaADJ = (double)alphaInt/10;
    sourceImage.convertTo(sourceImage, -1, alphaADJ, betaADJ);
    */
    std::pair<float, float> quantiles;
    Mat sourceImageFloat;
    if(alphaInt==0){
        
        sourceImage.convertTo(sourceImageFloat,CV_32F, 1.0); //-255.0/Min
        
        get_quantile(sourceImageFloat, sourceImageFloat.total() * 0.001, sourceImageFloat.total() * 0.999, alphaInt, betaADJ); 
            if(printOutput){
            Rcpp::Rcout << "Auto range clipping with 0.1%: "<< std::endl;
            Rcpp::Rcout << "Max: " << alphaInt << std::endl;
            Rcpp::Rcout << "Min: " << betaADJ << std::endl;
        }

    }

    sourceImage -= betaADJ;
        sourceImage.convertTo(sourceImage,CV_8UC1,255.0/(alphaInt-betaADJ)); //-255.0/Min
        
        //cvtColor(sourceImage, sourceImage, CV_GRAY2RGB);
        //

    Mat dstImage;
    sourceImage.copyTo(dstImage);

    if((maxdim/pow(2, tiers))>256){
      tiers++;
    }

    int rows;
    int cols;
    Mat tileInput;

    int barWidth = 70;
    float progress = 0.0;


    if(printOutput){
    Rcpp::Rcout << "Zoom levels: " << tiers << std::endl;
    }

    int numberOfTiles = 0;
    int tileGroupNumber = -1;

    tiers++;
    for (int tierInd = 0; tierInd < tiers; tierInd++)
    {

      double downsamplingfactor = 1/pow(2, tiers-1-tierInd);
      resize(sourceImage, dstImage, Size(), downsamplingfactor, downsamplingfactor);


      rows = (dstImage.rows / 256) + (dstImage.rows % 256 ? 1 : 0);
      cols = (dstImage.cols / 256) + (dstImage.cols % 256 ? 1 : 0);


      for (int rowTile = 0; rowTile < rows; rowTile++) 
      {
        for(int colTile = 0; colTile < cols; colTile++)
        {

          if(numberOfTiles % 256 == 0){
            tileGroupNumber++;
            string tilegrouppath;
            tilegrouppath = "mkdir \"Web_" + off + "/" + "Tiles_" + off + "/TileGroup" + to_string(tileGroupNumber) + "\"";
            system(tilegrouppath.c_str());
          }

            int width = 256;
            int height = 256;

            if(((colTile * 256)+width)>dstImage.cols ){
              width = dstImage.cols - (colTile * 256);
            }
            if(((rowTile * 256)+height)>dstImage.rows ){
              height = dstImage.rows - (rowTile * 256);
            }

            cv::Rect srcTile(colTile * 256,
                             rowTile * 256,
                             width,
                             height );
            
            tileInput = dstImage(srcTile);
            
            string tier_number = static_cast<ostringstream*>( &(ostringstream() << tierInd) )->str();
            string col_number = static_cast<ostringstream*>( &(ostringstream() << colTile) )->str();
            string row_number = static_cast<ostringstream*>( &(ostringstream() << rowTile) )->str();
            
            //imshow( String, displayImage ); // image visualisation
            string filepath;
            filepath = "Web_" + off + "/" + "Tiles_" + off + "/TileGroup" + to_string(tileGroupNumber) + "/" + tier_number + "-" + col_number + "-" + row_number + ".jpg";
            try {
                imwrite(filepath, tileInput);
                
            }
            catch (runtime_error& ex) {
                Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
                return(R_NilValue);
            }
            
            numberOfTiles++;
            
        }
    }

    if(printOutput){
    Rcpp::Rcout << "  [";
    int pos = barWidth * progress;
    for (int p = 0; p < barWidth; ++p) {
        if (p < pos) Rcpp::Rcout << "=";
        else if (p == pos) Rcpp::Rcout << ">";
        else Rcpp::Rcout << " ";
    }
        Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::cout.flush();
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
    progress += (float)1/(tiers-1);
    }


      

    }  


    /*
    int rows;
    int cols;
    int overlapPixels;
    int newHeight;
    int newWidth; 
    if(overlapPerc==-999){
      rows = (int)( ((double)sourceImage.rows / (double)tileSize)+0.5);
      cols = (int)( ((double)sourceImage.cols / (double)tileSize)+0.5);
  
      Rcpp::Rcout << "Rows: " << rows << std::endl;
      Rcpp::Rcout << "Cols: " << cols << std::endl; 

      int overlapRows = abs( (sourceImage.rows - rows*tileSize) )/rows;
      int overlapCols = abs( (sourceImage.cols - cols*tileSize) )/cols;


      if(overlapCols==overlapRows){
        overlapPixels = overlapRows/2;
      }else if(overlapCols>overlapRows){
        //correct the number of rows
        tileSize = tileSize - overlapCols;
        newHeight = (rows)*(tileSize+overlapCols)-overlapCols*(rows-1);
        newWidth = (cols)*(tileSize+overlapCols)-overlapCols*(cols-1);
        overlapPixels = overlapCols/2;

        

      }else{
        //correct the number of cols
        newHeight = (rows)*(tileSize)-overlapRows*(rows-2);
        newWidth = (cols)*(tileSize)-overlapRows*(cols-2);
        overlapPixels = overlapRows/2;

        
      }

    }else{
      overlapPixels = (tileSize*((double)overlapPerc/100))/2;

      rows = (sourceImage.rows / tileSize) + (sourceImage.rows % tileSize ? 1 : 0);
      cols = (sourceImage.cols / tileSize) + (sourceImage.cols % tileSize ? 1 : 0);
    }

    finalImage = Mat::zeros(newHeight, newWidth, sourceImage.type());

    int x0;
    int y0;

    switch (pos) {
        case 1:
            //topleft
             x0 =   0;
             y0 = 0;
            break;
        case 2:
            //topright
             x0 =   0;
             y0 = 0;
            break;
        case 3:
            //bottomright
             x0 =   (sourceImage.cols - finalImage.cols);
             y0 = (sourceImage.rows - finalImage.rows);
            break;
        case 4:
            //bottomleft
             x0 = 0;
             y0 = (sourceImage.rows - finalImage.rows);
            break;
        }


    Mat dst_roi = sourceImage(Rect(x0,y0,finalImage.cols,finalImage.rows));
    dst_roi.copyTo(finalImage); 
        

    tileSize = tileSize - 2 * overlapPixels;


    Rcpp::Rcout << "Overlap:" << overlapPixels << std::endl;
    Rcpp::Rcout << "tileSize:" << tileSize << std::endl;
    
    double empiricalOverlap = (double)(overlapPixels*2)/(double)(overlapPixels*2+tileSize);
 
    cv::Mat tileInput, tileOutput;
    
    int k = 0;
    string filepath;
    for (int rowTile = 0; rowTile < rows; rowTile++)
    {
        for (int colTile = 0; colTile < cols; colTile++)
        {
            cv::Rect srcTile(colTile * tileSize + overlapPixels,
                             rowTile * tileSize + overlapPixels,
                             tileSize + 2 * overlapPixels,
                             tileSize + 2 * overlapPixels);
            
            
            copySrcTile(finalImage, tileInput, srcTile);
            
            string String = static_cast<ostringstream*>( &(ostringstream() << k) )->str();
            
            //imshow( String, displayImage ); // image visualisation
            string filepath;
            filepath = "Tiled_" + off + "/" + off + "_Tile" + String + ".tif";
            try {
                imwrite(filepath, tileInput);
                
            }
            catch (runtime_error& ex) {
                Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
                return(R_NilValue);
            }
            k++;
            
            
        }
    }
    */
    return List::create(
                        _["width"] = sourceImage.cols,
                        _["height"] = sourceImage.rows
                        );

    return(R_NilValue);
    END_RCPP
}