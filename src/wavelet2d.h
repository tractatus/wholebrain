#ifndef WAVELET2D_H
#define WAVELET2D_H
#include <vector>
#include <complex>
using namespace std;



// 1D Functions
// void* dwt(vector<double> &, int ,string , vector<double> &,  vector<double> &);

 void* dwt1(string, vector<double> &, vector<double> &, vector<double> &);

 void* dyadic_zpad_1d(vector<double> &);

 double convol(vector<double> &, vector<double> &, vector<double> &);

 int filtcoef(string , vector<double> &, vector<double> &, vector<double> &,
                vector<double> &);

 void downsamp(vector<double> &, int , vector<double> &);

 void upsamp(vector<double> &, int, vector<double> &);

 void circshift(vector<double> &, int );

 int sign(int);

 void* idwt1(string wname, vector<double> &, vector<double> &, vector<double> &);

 int vecsum(vector<double> &, vector<double> &, vector<double> &);



// 1D Symmetric Extension DWT Functions



 void* dwt_sym(vector<double> &, int ,string , vector<double> &,vector<double> &,
                vector<int> &);

 void* dwt1_sym(string , vector<double> &, vector<double> &, vector<double> &);

 void* idwt_sym(vector<double> &,vector<double> &, string,vector<double> &, vector<int> &);

 void* symm_ext(vector<double> &, int );

 void* idwt1_sym(string, vector<double> &, vector<double> &, vector<double> &); // Not Tested

// 1D Stationary Wavelet Transform

 void* swt(vector<double> &, int , string , vector<double> &, int &) ;

 void* iswt(vector<double> &,int , string, vector<double> &);

 void* per_ext(vector<double> &, int );




// 2D Functions

 void* branch_lp_dn(string , vector<double> &, vector<double> &);

 void* branch_hp_dn(string , vector<double> &, vector<double> &);

 void* branch_lp_hp_up(string ,vector<double> &, vector<double> &, vector<double> &);

// void* dwt_2d(vector<vector<double> > &, int , string , vector<vector<double> > &
  //              , vector<double> &) ;

// void* idwt_2d(vector<vector<double> > &,vector<double> &, string ,vector<vector<double> > &);

 void* dyadic_zpad_2d(vector<vector<double> > &,vector<vector<double> > &);

 void* dwt_output_dim(vector<vector<double> >&, int &, int & );

 void* zero_remove(vector<vector<double> > &,vector<vector<double> > &) ;

 void* getcoeff2d(vector<vector<double> > &, vector<vector<double> > &,
                vector<vector<double> > &,vector<vector<double> > &,vector<double> &, int &);

 void* idwt2(string ,vector<vector<double> > &, vector<vector<double> >  &,
                vector<vector<double> >  &, vector<vector<double> >  &, vector<vector<double> > &);

 void* dwt2(string ,vector<vector<double> > &, vector<vector<double> >  &,
                vector<vector<double> >  &, vector<vector<double> > &, vector<vector<double> > &);

 void* downsamp2(vector<vector<double> > &,vector<vector<double> > &, int, int);

 void* upsamp2(vector<vector<double> > &,vector<vector<double> > &, int, int);

// 2D DWT (Symmetric Extension) Functions

 void* dwt_2d_sym(vector<vector<double> > &, int , string , vector<double> &, vector<double> & ,
                vector<int> &);

 void* dwt2_sym(string ,vector<vector<double> > &, vector<vector<double> >  &,
                vector<vector<double> >  &, vector<vector<double> >  &, vector<vector<double> > &);

 void* idwt_2d_sym(vector<double>  &,vector<double> &, string ,vector<vector<double> > &,
                vector<int> &);

 void* circshift2d(vector<vector<double> > &, int , int );

 void symm_ext2d(vector<vector<double> > &,vector<vector<double> > &, int );

 void* dispDWT(vector<double> &,vector<vector<double> > &, vector<int> &, vector<int> &, int ) ;

 void* dwt_output_dim_sym(vector<int> &,vector<int> &, int );

//2D Stationary Wavelet Transform

 void* swt_2d(vector<vector<double> > &,int , string , vector<double> &);

 void* per_ext2d(vector<vector<double> > &,vector<vector<double> > &, int );

// FFT functions


 double convfft(vector<double> &, vector<double> &, vector<double> &);

 double convfftm(vector<double> &, vector<double> &, vector<double> &);

 void* fft(vector<complex<double> > &,int ,unsigned int);

 void* bitreverse(vector<complex<double> > &);

 void* freq(vector<double> &, vector<double> &);

//New


 void* dwt1_sym_m(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD);//FFTW3 for 2D

 void* idwt1_sym_m(string wname, vector<double> &X, vector<double> &app, vector<double> &detail);

 void* dwt(vector<double> &sig, int J, string nm, vector<double> &dwt_output
                , vector<double> &flag, vector<int> &length );

 void* idwt(vector<double> &,vector<double> &, string,vector<double> &, vector<int> &);

 void* dwt_2d(vector<vector<double> > &, int , string , vector<double> &, vector<double> & ,
                vector<int> &);
 void* dwt1_m(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD) ;

 void* idwt_2d(vector<double>  &dwtop,vector<double> &flag, string nm,
                vector<vector<double> > &idwt_output, vector<int> &length);

 void* idwt1_m(string wname, vector<double> &X, vector<double> &cA, vector<double> &cD);

 void* dwt_output_dim2(vector<int> &length, vector<int> &length2, int J);


#endif/* WAVELET2S_H */
