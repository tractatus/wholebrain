#include <Rcpp.h>

// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;

//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}

int connectedcomponents3D ( int l, int m, int n, vector<Mat>& src, vector<Mat>& dst, vector<double>&  centerofmassZ,  vector<double>&  centerofmassX,  vector<double>&  centerofmassY) // int a[], int c[]

//****************************************************************************80
//
//  Purpose:
//
//    connectedcomponents3D assigns contiguous nonzero pixels to a common component.
//
//  Discussion:
//
//    On input, the A array contains values of 0 or 1.
//
//    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
//    into connected components.
//
//    The pixel A(I,J,K) is "connected" to the pixels:
//
//      A(I-1,J,  K  ),  A(I+1,J,  K  ),
//      A(I,  J-1,K  ),  A(I,  J+1,K  ),
//      A(I,  J,  K-1),  A(I,  J,  K+1),
//
//    so most pixels have 6 neighbors.
//
//    On output, COMPONENT_NUM reports the number of components of nonzero
//    data, and the array C contains the component assignment for
//    each nonzero pixel, and is 0 for zero pixels.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L, M, N, the order of the array.
//
//    Input, int A[L*M*N], the pixel array.
//
//    Output, int C[L*M*N], the component array.
//
//    Output, int I4BLOCK_COMPONENTS, the number of components
//    of nonzero data.
//
{
  int b;
  int c1;
  int component;
  int component_num;
  int i;
  int j;
  int k;
  int north;
  int *p;
  int *q;
  int up;
  int west;
//
//  Initialization.  src.at<uchar>(i,j)
//
  /*for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        c[i+j*l+k*l*m] = 0;
      }
    }
  } */
  component_num = 0;
//
//  P is simply used to store the component labels.  The dimension used
//  here is, of course, usually an absurd overestimate.
//
  p = new int[l*m*n+1];
  for ( i = 0; i <= l * m * n; i++ )
  {
    p[i] = i;
  }
//
//  "Read" the array one pixel at a time.  If a (nonzero) pixel has a north or
//  west neighbor with a label, the current pixel inherits it.
//  In case the labels disagree, we need to adjust the P array so we can
//  later deal with the fact that the two labels need to be merged.
//
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        //get the array content in each direction N-S, W-E, U-D
        //NORTH-SOUTH
        if ( i == 0 )
        {
          north = 0;
        }
        else
        {
          north = dst.at(i-1).at<ushort>(j,k); //north = c[i-1+j*l+k*l*m];
        }
        //WEST-EAST
        if ( j == 0 )
        {
          west = 0;
        }
        else
        {
          west = dst.at(i).at<ushort>(j-1,k); //west = c[i+(j-1)*l+k*l*m];
        }
        //UP-DOWN
        if ( k == 0 )
        {
          up = 0;
        }
        else
        {
          up = dst.at(i).at<ushort>(j,k-1);//up = c[i+j*l+(k-1)*l*m];
        }
        //check if target pixel is empty  a[i+j*l+k*l*m] != 0
        if ( src.at(i).at<ushort>(j,k) != 0 )
        {
//
//  New component?
//
          if ( north == 0 && west == 0 && up == 0 )
          {
            component_num = component_num + 1;
            dst.at(i).at<ushort>(j,k) = component_num;//  c[i+j*l+k*l*m] = component_num;
          }
//
//  One predecessor is labeled.
//
          else if ( north != 0 && west == 0 && up == 0 )
          {
            dst.at(i).at<ushort>(j,k) = north; // c[i+j*l+k*l*m] = north;
          }
          else if ( north == 0 && west != 0 && up == 0 )
          {
            dst.at(i).at<ushort>(j,k) = west; //c[i+j*l+k*l*m] = west;
          }
          else if ( north == 0 && west == 0 && up != 0 )
          {
            dst.at(i).at<ushort>(j,k) = up; //c[i+j*l+k*l*m] = up;
          }
//
//  Two predecessors are labeled.
//
          else if ( north == 0 && west != 0 && up != 0 )
          {
            dst.at(i).at<ushort>(j,k) = i4_min ( west, up );//c[i+j*l+k*l*m] = i4_min ( west, up );
            c1 = i4_min ( p[west], p[up] );
            p[west] = c1;
            p[up] = c1;
          }
          else if ( north != 0 && west == 0 && up != 0 )
          {
            dst.at(i).at<ushort>(j,k) = i4_min ( north, up );//c[i+j*l+k*l*m] = i4_min ( north, up );
            c1 = i4_min ( p[north], p[up] );
            p[north] = c1;
            p[up] = c1;
          }
          else if ( north != 0 && west != 0 && up == 0 )
          {
            dst.at(i).at<ushort>(j,k) = i4_min ( north, west );//c[i+j*l+k*l*m] = i4_min ( north, west );
            c1 = i4_min ( p[north], p[west] );
            p[north] = c1;
            p[west] = c1;
          }
//
//  Three predecessors are labeled.
//
          else if ( north != 0 && west != 0 && up != 0 )
          {
            dst.at(i).at<ushort>(j,k) = i4_min ( north, i4_min ( west, up ) ); //c[i+j*l+k*l*m] = i4_min ( north, i4_min ( west, up ) );
            c1 = i4_min ( p[north], i4_min ( p[west], p[up] ) );
            p[north] = c1;
            p[west] = c1;
            p[up] = c1;
          }
        }
      }
    }
  }
//
//  When a component has multiple labels, have the higher labels
//  point to the lowest one.
//
  for ( component = component_num; 1 <= component; component-- )
  {
    b = component;
    while ( p[b] != b )
    {
      b = p[b];
    }
    p[component] = b;
  }
//
//  Locate the minimum label for each component.
//  Assign these mininum labels new consecutive indices.
//
  q = new int[component_num+1];

  for ( j = 0; j <= component_num; j++ )
  {
    q[j] = 0;
  }

  i = 0;
  for ( component = 1; component <= component_num; component++ )
  {
    if ( p[component] == component )
    {
      i = i + 1;
      q[component] = i;
    }
  }

  component_num = i;

//for statistics
  double X[component_num];
  double Y[component_num];
  double Z[component_num];
  uint number_voxels[component_num];
  for(i = 0; i < component_num; i++){
    X[i] = 0;
    Y[i] = 0;
    Z[i] = 0;
    number_voxels[i] = 0;
  }
//
//  Replace the labels by consecutive labels.
//
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        int index = dst.at(i).at<ushort>(j,k);
        int value = q [ p [  index ] ];
        dst.at(i).at<ushort>(j,k) = value; //q [ p [  index ] ]; // c[i+j*l+k*l*m] = q [ p [ c[i+j*l+k*l*m] ] ];
        Z[value] += (double)i/l;
        X[value] += (double)j/m;
        Y[value] += (double)k/n;
        number_voxels[value]++;
      }
    }
  }

  for(i = 0; i < component_num; i++){
    centerofmassZ.push_back( ((double)Z[i]/number_voxels[i])  );
    centerofmassX.push_back( ((double)X[i]/number_voxels[i])  );
    centerofmassY.push_back( ((double)Y[i]/number_voxels[i])  );
  }

  delete [] p;
  delete [] q; 

  return component_num;
}



/* apply operation to stack */
RcppExport SEXP runonstack(SEXP input, SEXP operation, SEXP outputfile) {
BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
  //Rcpp::CharacterVector std::vector< std::string >
  Rcpp::CharacterVector f(input);
  int num_files = f.size();
 // int size[3] = { src.cols, src.rows, f.size() };//a strack is a three dimensional array.
 // cv::Mat M(3, size, CV_8UC1, cv::Scalar(0));
  int barWidth = 70;
  float progress = 0.0;

  std::string ff(f[0]);
  Mat remove = imread(ff, CV_LOAD_IMAGE_GRAYSCALE);
  //vector<Mat> zpositions(num_files);
  //int size[3] = { 561, 561, num_files };
  //Mat stack(3, size, CV_8UC1, cv::Scalar(0));
  int rows = remove.rows;
  int cols = remove.cols;
  vector<Mat> zpositions(num_files,Mat(cols,rows,CV_16U));
  vector<Mat> destination(num_files,Mat(cols,rows,CV_16U));
  Rcpp::Rcout << "Loading " <<  num_files << " stacks into RAM." << std::endl;
  for (int i = 0; i < f.size(); i++){
    std::string filename(f[i]);

    Mat src = imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
    Mat tmp;
    threshold(src, tmp, 0, 1, CV_THRESH_BINARY | CV_THRESH_OTSU );
    tmp.convertTo(tmp, CV_16U);
    zpositions.at(i) = tmp;

    Mat dst = Mat::zeros(src.size(), CV_16U);

    destination.at(i) = dst;

    //print progress bar in console
    Rcpp::Rcout << "  [";
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
  Rcpp::Rcout << "\n" << std::endl;
  Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;

  /* merge(zpositions, stack);
 
  Rcpp::Rcout << "---- Dimensions: " << stack.dims << std::endl;
  Rcpp::Rcout << "---- Array size: ";
  for (int i=0; i<stack.dims; ++i)
      Rcpp::Rcout << stack.size[i] << ' x ';
  Rcpp::Rcout << std::endl; */

  Rcpp::Rcout << "Running command. " << "[INSERT COMMAND HERE]" << " on stack" << std::endl;
  vector<double> zpos;
  vector<double> xpos;
  vector<double> ypos;
  int mupp = connectedcomponents3D ( num_files, cols, rows, zpositions, destination, zpos, xpos, ypos );
  Rcpp::Rcout << "====== PROCESSING DONE ======" << std::endl;
  Rcpp::Rcout << "MUPP = "<< mupp << std::endl;

   /*
  imshow("input423", zpositions.at(423));
  imshow("input500", zpositions.at(500));
  imshow("input599", zpositions.at(599));

  imshow("output423", destination.at(423));
  imshow("output500", destination.at(500));
  imshow("output599", destination.at(599));
  */

  //Mat reference;
  //destination.at(499).convertTo(reference, CV_16U);

  return Rcpp::List::create(
    Rcpp::_["components"] = mupp,
    Rcpp::_["x"] = xpos,
    Rcpp::_["y"] = ypos,
    Rcpp::_["z"] = zpos
  );

END_RCPP  
}


// imshow('/Users/danielfurth/Documents/D2_1222_03_02_TileScan_005_Merging001_ch00.tif')

/* 
files<-c('/Users/danielfurth/Documents/D2_1222_03_02_TileScan_005_Merging001_ch00.tif', '/Users/danielfurth/Documents/mouse_atlas/lighsheet_figure/cropped_x0.25_unspmask3-0.6_s_0607.tif') 
stackapply(files) 
#~/Documents/mouse_atlas/lighsheet_figure/analysis')


library(wholebrain)
setwd('~/Documents/mouse_atlas/lighsheet_figure/soma_for_rendering')  
files.to.be.processed<-dir()


istif<-substr(files.to.be.processed, (nchar(files.to.be.processed)+1)-4, nchar(files.to.be.processed) ) 

files.to.be.processed <- files.to.be.processed[which(istif==".tif")] 
stackapply(files.to.be.processed) 


plot3d(df$x, df$y, df$z, size=2, type='s')

for (int i = 0; i < 100; i++)
  for (int j = 0; j < 100; j++)
    for (int k = 0; k < 3; k++) 
      std::cout << M.at<cv::Vec3f>(i,j)[k] << ", ";

*/

