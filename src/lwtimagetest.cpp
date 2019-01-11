//============================================================================
// Name        : imagedemo_sym.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   :
// Description : DWT of arbitrary size image using symmetric extension
//============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <algorithm>
#include "lwave.h"
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <ctime>
#include <Rcpp.h>
#include "wholebrain.h"

using namespace std;
using namespace cv;
using namespace Rcpp;

void* maxvalLWT(vector<vector<double> > &arr, double &max){
	max = 0;
	for (unsigned int i =0; i < arr.size(); i++) {
		for (unsigned int j =0; j < arr[0].size(); j++) {
			if (max <= arr[i][j]){
				max = arr[i][j];
			}
		}
	}
	return 0;
}

void* maxvalLWT1(vector<double> &arr, double &max){
	max = 0;
	for (unsigned int i =0; i < arr.size(); i++) {
		if (max <= arr[i]){
			max = arr[i];
		}

	}
	return 0;
}




RcppExport SEXP liftScheme(SEXP input) {
    BEGIN_RCPP
    Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>
    
    Rcpp::CharacterVector f(input);
    std::string ff(f[0]);
    Rcpp::Rcout << "Loading image:" << ff << std::endl;
    clock_t t0, t1;
    
	Mat img = imread(ff, -1);
	/* if (!img){
		cout << " Can't read Image. Try Different Format." << endl;
		exit(1);
	} */
	int height, width;
	height = img.rows;
	width = img.cols;
	int nc = img.channels();
	//   uchar* ptr2 =(uchar*) img->imageData;
	int pix_depth = img.depth();
    cv::Size size;
	size.width =width;
	size.height=height;
	cout << "depth" << pix_depth <<  "Channels" << nc << endl;


	//cvNamedWindow("Original Image", CV_WINDOW_AUTOSIZE);
	imshow("Original Image", img);
	waitKey(0);
	destroyAllWindows();
	//cvSaveImage("orig.bmp",img);


	int rows =(int) height;
	int cols =(int) width;
	Mat matimg(img);

	vector<double> vec1,copy;


	int k =1;
	for (int i=0; i < rows; i++) {
		for (int j =0; j < cols; j++){
			unsigned char temp;
			temp = ((uchar*) matimg.data + i * matimg.step)[j  * matimg.elemSize() + k ];
			vec1.push_back((double) temp);
		}

	}

	copy=vec1;

    int J=6;
	string nm = "haar";
	liftscheme lift(nm);
    t0 = clock();
	lwt2<double> wt(vec1,rows,cols,lift,J);
    t1 = clock();
    double time_elapsed = (double)((t1 - t0) / CLOCKS_PER_SEC);
    Rcpp::Rcout << "LIFTING." << " took " <<  time_elapsed << " seconds." << std::endl;
    
    

	// unsigned int lf=l1.size();
	//  int rows_n =(int) (rows+ J*(lf-1));
	//  int cols_n =(int)  (cols + J * ( lf -1));

	// Finding 2D DWT Transform of the image using symetric extension algorithm
	// Extension is set to 3 (eg., int e = 3)
/*
	vector<int> length;
	vector<double> output,flag;
	int J =2;
	time_t start,end;
    time (&start);
	dwt_2d_sym(vec1,J,nm,output,flag,length);
	time (&end);

	double max;
	vector<int> length2;
	// This algorithm computes DWT of image of any given size. Together with convolution and
	// subsampling operations it is clear that subsampled images are of different length than
	// dyadic length images. In order to compute the "effective" size of DWT we do additional
	// calculations.
	dwt_output_dim_sym(length,length2,J);
	// length2 is gives the integer vector that contains the size of subimages that will
	// combine to form the displayed output image. The last two entries of length2 gives the
	// size of DWT ( rows_n by cols_n)

	int siz = length2.size();
	int rows_n=length2[siz-2];
	int cols_n = length2[siz-1];

	vector<vector< double> > dwtdisp(rows_n, vector<double>(cols_n));
	dispDWT(output,dwtdisp, length ,length2, J);

	// dispDWT returns the 2D object dwtdisp which will be displayed using OPENCV's image
	// handling functions

	//Extracting Individual Images at the Jth level

	int jrow = length[0];
	int jcol= length[1];

	cout << "No. of rows" << jrow << "No. of columns" << jcol;


	//Let 4 images at Jth level be A,Dh,Dv,Dd
	vector<vector< double> > A(jrow, vector<double>(jcol));
	vector<vector< double> > Dh(jrow, vector<double>(jcol));
	vector<vector< double> > Dv(jrow, vector<double>(jcol));
	vector<vector< double> > Dd(jrow, vector<double>(jcol));
	for (int i=0;i < jrow;i++){
	    for (int j=0; j < jcol;j++) {
	        A[i][j]=output[i*jcol+j];
	        Dh[i][j]=output[jrow*jcol+i*jcol+j];
	        Dv[i][j]=output[2*jrow*jcol+i*jcol+j];
	        Dd[i][j]=output[3*jrow*jcol+i*jcol+j];

	    }
	}

	vector<vector<double> >  dwt_output= dwtdisp;

	maxval(dwt_output,max);// max value is needed to take care of overflow which happens because
	// of convolution operations performed on unsigned 8 bit images

	//Displaying Scaled Image
	// Creating Image in OPENCV
	IplImage *cvImg; // image used for output
	CvSize imgSize; // size of output image

	imgSize.width = cols_n;
	imgSize.height = rows_n;

	cvImg = cvCreateImage( imgSize, 8, 1 );
	// dwt_hold is created to hold the dwt output as further operations need to be
	// carried out on dwt_output in order to display scaled images.
	vector<vector<double> > dwt_hold(rows_n, vector<double>( cols_n));
	dwt_hold = dwt_output;
	// Setting coefficients of created image to the scaled DWT output values
	for (int i = 0; i < imgSize.height; i++ ) {
		for (int j = 0; j < imgSize.width; j++ ){
			if ( dwt_output[i][j] <= 0.0){
				dwt_output[i][j] = 0.0;
			}
			if ( i <= (length2[0]) && j <= (length2[1]) ) {
				((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
						(char) ( (dwt_output[i][j] / max) * 255.0);
			} else {
				((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
						(char) (dwt_output[i][j]) ;
			}
		}
	}

	cvNamedWindow( "DWT Image", 1 ); // creation of a visualisation window
	cvShowImage( "DWT Image", cvImg ); // image visualisation
	cvWaitKey();
	cvDestroyWindow("DWT Image");
	cvSaveImage("dwt.bmp",cvImg);


*/
	 return(R_NilValue);
    END_RCPP
}
