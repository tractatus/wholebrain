//============================================================================
// Name        : 1D/2D Wavelet Transform
// Author      : Rafat Hussain
// Version     :
// Copyright   : GNU GPL License
// Description : Wavelet Library
//============================================================================
/*Copyright (C) 2011 Rafat Hussain

 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program; if not, write to the Free Software
 *Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
*/

#include <iostream>
#include <fstream>
#include <complex>
#include "wavelet2d.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
using namespace std;

extern "C" int _get_output_format( void ){ return 0; }
fftw_plan plan_forward_inp,plan_forward_filt, plan_backward;
static unsigned int transient_size_of_fft = 0;



void* per_ext2d(vector<vector<double> > &signal,vector<vector<double> > &temp2, int a) {

        unsigned int rows = signal.size();
                unsigned int cols = signal[0].size();
                int cols2;
                if ((cols % 2) != 0) {
                    cols2 = cols + 1;
                } else {
                    cols2 = cols;
                }
                vector<vector<double> > temp_vec(rows ,vector<double>(cols2 + 2* a));
        //	vector<vector<double> > temp2(rows + 2 * a ,vector<double>(cols + 2* a));

                for (unsigned int i=0; i < rows; i++) {
                        vector<double> sig;
                        for (unsigned int j=0; j< cols; j++) {
                                double temp = signal[i][j];
                                sig.push_back(temp);
                        }
                        per_ext(sig,a);
                        for (unsigned int j=0; j< sig.size(); j++) {
                  temp_vec[i][j] = sig[j];
                        }
                }
                for (unsigned int j=0; j < temp_vec[0].size(); j++) {
                        vector<double> sig;
                        for (unsigned int i=0; i< rows; i++) {
                                double temp = temp_vec[i][j];
                                sig.push_back(temp);
                        }
                        per_ext(sig,a);
                        for (unsigned int i=0; i< sig.size(); i++) {
                  temp2[i][j] = sig[i];
                        }
                }


        return 0;
}

void* swt_2d(vector<vector<double> > &sig,int J, string nm, vector<double> &swt_output) {
    int m_size = sig.size(); // No. of rows
    int n_size = sig[0].size(); //No. of columns


    vector<vector<double> > sig2 =sig;
    int rows_n =m_size;
    int cols_n =n_size;
    vector<double> lp1,hp1,lp2,hp2;
    filtcoef(nm,lp1,hp1,lp2,hp2);

    for (int iter =0; iter < J; iter++) {
    int U = (int) pow(2.0,(double)iter);

    vector<double> low_pass, high_pass;
    if(iter > 0) {
        upsamp(lp1,U,low_pass);
        upsamp(hp1,U,high_pass);
    } else {
        low_pass = lp1;
        high_pass = hp1;
    }

    int lf = low_pass.size();

    if ((sig.size() % 2) == 0 ) {
        rows_n = sig.size();
    } else {
        rows_n = sig.size() + 1;
    }

    if ((sig[0].size() % 2) == 0 ) {
        cols_n = sig[0].size();
    } else {
        cols_n = sig[0].size() + 1;
    }

    vector<vector<double> > signal(rows_n + lf,vector<double>(cols_n+lf));
    //    per_ext2d(sig,signal,lf/2); Edit per_ext if you want to use per_ext2d. Remove
    // the even indexing.

    symm_ext2d(sig,signal,lf/2);
    int len_x = signal.size();
    int len_y = signal[0].size();
    vector<vector<double> > sigL(rows_n + lf,vector<double>(cols_n));
    vector<vector<double> > sigH(rows_n + lf,vector<double>(cols_n));
    vector<vector<double> > cA(rows_n,vector<double>(cols_n));
    vector<vector<double> > cH(rows_n,vector<double>(cols_n));
    vector<vector<double> > cV(rows_n,vector<double>(cols_n));
    vector<vector<double> > cD(rows_n,vector<double>(cols_n));

    for (int i=0; i < len_x; i++) {
        vector<double> temp_row;
        for (int j=0; j < len_y; j++) {
                double temp = signal[i][j];
                temp_row.push_back(temp);

        }

        // ------------------Low Pass Branch--------------------------

                vector<double> oup;
                convfftm(temp_row,low_pass,oup);
        oup.erase(oup.begin(), oup.begin()+lf);
                oup.erase(oup.begin()+cols_n,oup.end());

        // ------------------High Pass Branch--------------------------

                vector<double> oup2;
                convfftm(temp_row,high_pass,oup2);
        oup2.erase(oup2.begin(), oup2.begin()+lf);
                oup2.erase(oup2.begin()+cols_n,oup2.end());

                temp_row.clear();

        for (unsigned int j=0; j < oup.size() ; j++) {
                sigL[i][j] = oup[j];
            sigH[i][j] = oup2[j];
        }

    }

    for (int j=0; j < cols_n; j++) {
        vector<double> temp_row;
        for (int i=0; i < len_x; i++){
                double temp = sigL[i][j];
                temp_row.push_back(temp);
        }

        // ------------------Low Pass Branch--------------------------


                vector<double> oup;
                convfftm(temp_row,low_pass,oup);
                oup.erase(oup.begin(), oup.begin()+lf);
                oup.erase(oup.begin()+rows_n,oup.end());

                // ------------------High Pass Branch--------------------------

                vector<double> oup2;
                convfftm(temp_row,high_pass,oup2);
            oup2.erase(oup2.begin(), oup2.begin()+lf);
                oup2.erase(oup2.begin()+rows_n,oup2.end());

                temp_row.clear();


                for (unsigned int i=0; i < oup.size() ; i++) {
                        cA[i][j] = oup[i];
                }

                for (unsigned int i=0; i < oup2.size() ; i++) {
                cH[i][j] = oup2[i];
                }


    }

    for (int j=0; j < cols_n; j++) {
                vector<double> temp_row;
                for (int i=0; i < len_x; i++){
                        double temp = sigH[i][j];
                        temp_row.push_back(temp);
                }

                // ------------------Low Pass Branch--------------------------


                        vector<double> oup;
                        convfftm(temp_row,low_pass,oup);
                oup.erase(oup.begin(), oup.begin()+lf);
                        oup.erase(oup.begin()+rows_n,oup.end());

                // ------------------High Pass Branch--------------------------

                        vector<double> oup2;
                        convfftm(temp_row,high_pass,oup2);
                oup2.erase(oup2.begin(), oup2.begin()+lf);
                        oup2.erase(oup2.begin()+rows_n,oup2.end());

                        temp_row.clear();


                for (unsigned int i=0; i < oup.size() ; i++) {
                        cV[i][j] = oup[i];
                }

                for (unsigned int i=0; i < oup2.size() ; i++) {
                    cD[i][j] = oup2[i];
                }


        }

       sig = cA;
       vector<double> temp_sig2;

       if (iter == J-1) {
                                for(int i =0; i < rows_n; i++){
                                        for (int j =0; j < cols_n; j++){
                                        double temp=cA[i][j];
                                        temp_sig2.push_back(temp);
                                        }
                                }
                           }
                                for(int i =0; i < rows_n; i++){
                                        for (int j = cols_n; j < cols_n * 2; j++){
                                        double temp =cH[i][j - cols_n];
                                        temp_sig2.push_back(temp);
                                        }
                                }

                                for(int i = rows_n; i < rows_n * 2; i++){
                                        for (int j =0; j < cols_n; j++){
                                        double temp=cV[i - rows_n][j];
                                        temp_sig2.push_back(temp);
                                        }
                                }

                                for(int i = rows_n; i < rows_n * 2; i++){
                                        for (int j = cols_n; j < cols_n * 2; j++){
                                        double temp =cD[i- rows_n][j - cols_n];
                                         temp_sig2.push_back(temp);
                                        }
                                }

                                swt_output.insert(swt_output.begin(),temp_sig2.begin(),temp_sig2.end());



    }


        return 0;
}


void* per_ext(vector<double> &sig, int a) {
        unsigned int len;
    len = sig.size();
    if ((len % 2) != 0 ) {
        double temp = sig[len-1];
        sig.push_back(temp);
        len = sig.size();
    }

        for (int i=0; i < a; i++) {
    double temp1 = sig[2 *i];
    double temp2 = sig[len-1];
    sig.insert(sig.begin(), temp2);
    sig.insert(sig.end(), temp1);

        }

        return 0;

}


void* iswt(vector<double> &swtop,int J, string nm, vector<double> &iswt_output) {
         int N = swtop.size() / (J + 1);

     vector<double> lpd,hpd,lpr,hpr;
         filtcoef(nm,lpd,hpd,lpr,hpr);

         vector<double> appx_sig;

     vector<double> low_pass = lpr;
     vector<double> high_pass = hpr;
     int lf = low_pass.size();

     for (int iter = 0; iter < J; iter++) {
     vector<double> det_sig;
         if (iter ==0) {
         for (int i = 0; i < N; i++) {
                 double temp=swtop[i];
                 appx_sig.push_back(temp);
                 double temp1=swtop[(iter + 1) * N + i];
                 det_sig.push_back(temp1);
         }
         } else {
                 for (int i = 0; i < N; i++) {
                 double temp1=swtop[(iter + 1) * N + i];
                 det_sig.push_back(temp1);

         }
         }


           int value =(int) pow(2.0,double(J -1 -iter));
           iswt_output.assign(N,0.0);

          for (int count = 0; count < value; count++) {
        vector<double> appx1, det1;
        for (int index = count; index < N; index+=value){
                double temp = appx_sig[index];
                appx1.push_back(temp);
                double temp1 = det_sig[index];
                det1.push_back(temp1);

        }
        unsigned int len = appx1.size();

        // Shift = 0

         vector<double> appx2, det2;

         for (unsigned int index_shift =0; index_shift < len; index_shift+=2) {
                double temp = appx1[index_shift];
                appx2.push_back(temp);
                double temp1 = det1[index_shift];
                det2.push_back(temp1);
         }

         int U = 2; // Upsampling Factor

         vector<double> cL0,cH0;
         upsamp(appx2,U,cL0);
         upsamp(det2,U,cH0);
         per_ext(cL0,lf/2);
         per_ext(cH0,lf/2);

         vector<double> oup00L, oup00H, oup00;
         convfft(cL0,low_pass,oup00L);
         convfft(cH0,high_pass,oup00H);

         oup00L.erase(oup00L.begin(),oup00L.begin()+lf - 1);
                 oup00L.erase(oup00L.begin()+len,oup00L.end());
         oup00H.erase(oup00H.begin(),oup00H.begin()+lf - 1);
                 oup00H.erase(oup00H.begin()+len,oup00H.end());

                 vecsum(oup00L,oup00H,oup00);

                 // Shift = 1

                 vector<double> appx3, det3;

                         for (unsigned int index_shift =1; index_shift < len; index_shift+=2) {
                                double temp = appx1[index_shift];
                                appx3.push_back(temp);
                                double temp1 = det1[index_shift];
                                det3.push_back(temp1);
                         }


                         vector<double> cL1,cH1;
                         upsamp(appx3,U,cL1);
                         upsamp(det3,U,cH1);
                         per_ext(cL1,lf/2);
                         per_ext(cH1,lf/2);

                         vector<double> oup01L, oup01H, oup01;
                         convfft(cL1,low_pass,oup01L);
                         convfft(cH1,high_pass,oup01H);

                         oup01L.erase(oup01L.begin(), oup01L.begin()+lf - 1);
                                 oup01L.erase(oup01L.begin()+len,oup01L.end());
                         oup01H.erase(oup01H.begin(), oup01H.begin()+lf - 1);
                                 oup01H.erase(oup01H.begin()+len,oup01H.end());

                                 vecsum(oup01L,oup01H,oup01);
                                 circshift(oup01,-1);

         //   Continue
              int index2 = 0;
                        for (int index = count; index < N; index+=value){
                                        double temp = oup00[index2]+oup01[index2];
                                        iswt_output.at(index) = temp/2;
                                        index2++;

                                }


     }
                appx_sig = iswt_output;


     }
        return 0;
}

void* swt(vector<double> &signal1, int J, string nm, vector<double> &swt_output, int &length) {
        vector<double> lpd, hpd, lpr, hpr;
        vector<double> sig = signal1;

            int N = sig.size();
            length = N;

           filtcoef(nm,lpd,hpd,lpr,hpr);

        for (int iter = 0; iter < J; iter++) {
                vector<double> low_pass;
                vector<double> high_pass;
              if ( iter > 0){

                  int M = (int) pow(2.0,iter);
                  upsamp(lpd,M,low_pass);
                  upsamp(hpd,M,high_pass);


              } else {
                  low_pass = lpd;
                  high_pass = hpd;
              }

              unsigned int len_filt = low_pass.size();
               per_ext(sig,len_filt/2);

                        vector<double> cA;
                        convfft(sig,low_pass,cA);
                        vector<double> cD;
                        convfft(sig,high_pass,cD);
                        // Resize cA and cD
                cA.erase(cA.begin(), cA.begin()+len_filt);
                        cA.erase(cA.begin()+N,cA.end());
                cD.erase(cD.begin(), cD.begin()+len_filt);
                        cD.erase(cD.begin()+N,cD.end());
                        // Reset signal value;

                        sig = cA;

             if (iter == J - 1 ) {
             swt_output.insert(swt_output.begin(),cD.begin(),cD.end());
             swt_output.insert(swt_output.begin(),cA.begin(),cA.end());
             } else {
                 swt_output.insert(swt_output.begin(),cD.begin(),cD.end());
             }

                }

  return 0;
}

void* dwt_output_dim_sym(vector<int> &length,vector<int> &length2, int J) {
        unsigned int sz=length.size();
        int rows = length[sz-2];
        int cols = length[sz-1];
    for (int i =0; i < J; i++) {
     rows =(int) ceil((double) rows/ 2.0);
     cols =(int) ceil((double) cols/ 2.0);
     }
    for (int i =0; i < J + 1; i++) {
        length2.push_back(rows);
        length2.push_back(cols);
        rows = rows * 2;
        cols = cols*2;
    }
        return 0;
}

void* dwt_output_dim2(vector<int> &length, vector<int> &length2, int J) {

   int row = length[0];
   int col = length[1];

   for (int i=0; i < J + 1; i++) {

       length2.push_back(row);
       length2.push_back(col);
       row = row * 2;
       col = col * 2;


   }


    return 0;
}

void* dispDWT(vector<double> &output,vector<vector<double> > &dwtdisp, vector<int> &length , vector<int> &length2, int J) {
          int sum = 0;


        for (int iter =0; iter < J; iter++) {
                int d_rows=length[2*iter]-length2[2*iter];
                int d_cols=length[2*iter+1]-length2[2*iter + 1];


                int rows_n =length[2 * iter];
                int cols_n = length[2 * iter + 1];
                vector<vector<double> >  dwt_output(2 * rows_n, vector<double>(2 * cols_n));
        if (iter == 0) {
                for(int i =0; i < rows_n; i++){
                        for (int j =0; j < cols_n; j++){
                                dwt_output[i][j]=output[i*cols_n + j];
                        }
                }

                for(int i =0; i < rows_n; i++){
                        for (int j = cols_n; j < cols_n * 2; j++){
                                dwt_output[i][j]= output[rows_n * cols_n + i * cols_n + (j - cols_n)];
                        }
                }

                for(int i = rows_n; i < rows_n * 2; i++){
                        for (int j =0; j < cols_n; j++){
                                dwt_output[i][j]=output[2 * rows_n * cols_n+ (i - rows_n) * cols_n + j];
                        }
                }


                for(int i = rows_n; i < rows_n * 2; i++){
                        for (int j = cols_n; j < cols_n * 2; j++){
                                dwt_output[i][j]=output[3 * rows_n * cols_n+ (i -rows_n) * cols_n + (j -cols_n)];
                        }
                }
        } else {
                for(int i =0; i < rows_n; i++){
                                        for (int j = cols_n; j < cols_n * 2; j++){
                                                dwt_output[i][j]= output[sum + i * cols_n + (j - cols_n)];
                                        }
                                }

                                for(int i = rows_n; i < rows_n * 2; i++){
                                        for (int j =0; j < cols_n; j++){
                                                dwt_output[i][j]=output[sum + rows_n * cols_n+ (i - rows_n) * cols_n + j];
                                        }
                                }


                                for(int i = rows_n; i < rows_n * 2; i++){
                                        for (int j = cols_n; j < cols_n * 2; j++){
                                                dwt_output[i][j]=output[sum + 2 * rows_n * cols_n+ (i -rows_n) * cols_n + (j -cols_n)];
                                        }
                                }

        }

                int rows_x = length2[2*iter];
        int cols_x =length2[2*iter +1];

        int d_cols2 = (int) ceil( (double) (d_cols - 1) / 2.0);
        int d_rows2 = (int) ceil( (double) (d_rows - 1) / 2.0);
        if (iter ==0) {
        for(int i =0; i < rows_x; i++){
                                for (int j =0; j < cols_x; j++){
                                        if (i + d_rows -1 < 0){
                                                dwtdisp[i][j]=0;
                                        }
                                        else if (j + d_cols -1 < 0){
                                                dwtdisp[i][j]=0;
                                                } else {
                                        dwtdisp[i][j]=dwt_output[i+d_rows -1][j+d_cols -1];
                                                }
                                }
                        }
        }
                        for(int i =0; i < rows_x; i++){
                                for (int j = cols_x; j < cols_x * 2; j++){
                                        if (i + d_rows2 < 0){
                                                dwtdisp[i][j]=0;
                                        }
                                        else if (j + 2* (d_cols -1) +1 > (signed) dwt_output[0].size() - 1){
                                                dwtdisp[i][j]=0;
                                                                                                                        } else {
                                        dwtdisp[i][j]= dwt_output[i+d_rows2 ][j + 2* (d_cols -1)+1 ];
                                                                                                                        }
                                }
                        }

                        for(int i = rows_x; i < rows_x * 2; i++){
                                for (int j =0; j < cols_x; j++){
                                        if (i + 2* (d_rows -1) + 1 > (signed) dwt_output.size() - 1){
                                                dwtdisp[i][j]=0;
                                        }
                                        else if (j + d_cols2 < 0){
                                                dwtdisp[i][j]=0;
                                            } else {

                                        dwtdisp[i][j]=dwt_output[i+2 * (d_rows - 1) + 1 ][j+d_cols2 ];
                                            }
                                }
                        }

                        for(int i = rows_x; i < rows_x * 2; i++){
                                for (int j = cols_x; j < cols_x * 2; j++){

                                        if (i +  (d_rows -1) + 1 + d_rows2 > (signed) dwt_output.size() - 1){
                                                dwtdisp[i][j]=0;
                                                                                }
                                        else if (j + (d_cols -1) + 1 + d_cols2  > (signed) dwt_output[0].size() - 1){
                                                dwtdisp[i][j]=0;
                                                                                        } else {
                                        dwtdisp[i][j]=dwt_output[i +  (d_rows -1) + 1 + d_rows2 ][j + (d_cols -1) + 1 + d_cols2 ];
                                                                                        }
                                }
                        }
                        if (iter == 0) {
                                sum+= 4*rows_n*cols_n;
                        } else {
                                sum+= 3*rows_n * cols_n;
                        }

        }

        return 0;

}

void symm_ext2d(vector<vector<double> > &signal,vector<vector<double> > &temp2, int a) {
        unsigned int rows = signal.size();
        unsigned int cols = signal[0].size();
        vector<vector<double> > temp_vec(rows ,vector<double>(cols + 2* a));
//	vector<vector<double> > temp2(rows + 2 * a ,vector<double>(cols + 2* a));

        for (unsigned int i=0; i < rows; i++) {
                vector<double> sig;
                for (unsigned int j=0; j< cols; j++) {
                        double temp = signal[i][j];
                        sig.push_back(temp);
                }
                symm_ext(sig,a);
                for (unsigned int j=0; j< sig.size(); j++) {
          temp_vec[i][j] = sig[j];
                }
        }
        for (unsigned int j=0; j < temp_vec[0].size(); j++) {
                vector<double> sig;
                for (unsigned int i=0; i< rows; i++) {
                        double temp = temp_vec[i][j];
                        sig.push_back(temp);
                }
                symm_ext(sig,a);
                for (unsigned int i=0; i< sig.size(); i++) {
          temp2[i][j] = sig[i];
                }
        }

}

void* circshift2d(vector<vector<double> > &signal, int x, int y) {
        unsigned int rows = signal.size();
        unsigned int cols = signal[0].size();
        vector<vector<double> > temp_vec(rows,vector<double>(cols));

        for (unsigned int i=0; i < rows; i++) {
                vector<double> sig;
                for (unsigned int j=0; j< cols; j++) {
                        double temp = signal[i][j];
                        sig.push_back(temp);
                }
                circshift(sig,x);
                for (unsigned int j=0; j< cols; j++) {
          temp_vec[i][j] = sig[j];
                }
        }

        for (unsigned int j=0; j < cols; j++) {
                vector<double> sig;
                for (unsigned int i=0; i< rows; i++) {
                        double temp = temp_vec[i][j];
                        sig.push_back(temp);
                }
                circshift(sig,y);
                for (unsigned int i=0; i< rows; i++) {
          signal[i][j] = sig[i];
                }
        }
  return 0;
}

void* idwt_2d_sym(vector<double>  &dwtop,vector<double> &flag, string nm,
                vector<vector<double> > &idwt_output, vector<int> &length){
    int J =(int) flag[0];
    int rows =length[0];
    int cols =length[1];

    int sum_coef =0;
    vector<double> lp1,hp1,lp2,hp2;
    filtcoef(nm,lp1,hp1,lp2,hp2);
    unsigned int lf = lp1.size();
        vector<vector<double> >  cLL(rows, vector<double>(cols));


    for (int iter=0; iter < J; iter++) {

        int rows_n = length[2*iter];
        int cols_n = length[2*iter + 1];

        vector<vector<double> >  cLH(rows_n, vector<double>(cols_n));
        vector<vector<double> >  cHL(rows_n, vector<double>(cols_n));
        vector<vector<double> >  cHH(rows_n, vector<double>(cols_n));

        for (int i = 0 ; i < rows_n; i++ ){
                for (int j = 0; j < cols_n; j++){
                        if (iter == 0) {
                        cLL[i][j] = dwtop[sum_coef+ i * cols_n + j];

                cLH[i][j] = dwtop[sum_coef+ rows_n * cols_n+ i * cols_n + j];

                cHL[i][j] = dwtop[sum_coef+ 2 * rows_n * cols_n + i * cols_n + j];

                cHH[i][j] = dwtop[sum_coef+ 3* rows_n * cols_n + i * cols_n + j];
                        } else {

                        cLH[i][j] = dwtop[sum_coef+  i * cols_n + j];

                        cHL[i][j] = dwtop[sum_coef+ rows_n * cols_n + i * cols_n + j];

                    cHH[i][j] = dwtop[sum_coef+ 2* rows_n * cols_n + i * cols_n + j];

                        }
                }
        }


  //      temp_A = cLL;
  //  	idwt2_sym(nm,idwt_output2, cA, cH,cV,cD);

        unsigned int len_x = cLH.size();
        unsigned int len_y = cLH[0].size();

        // Row Upsampling and Column Filtering at the first LP Stage
        vector<vector<double> > cL(2 *len_x - lf + 2,vector<double>(len_y ));
        vector<vector<double> > cH(2 * len_x - lf +2,vector<double>(len_y ));

        if (iter ==0) {
        for (unsigned int j =0; j < len_y; j++) {

                        vector<double> sigLL,sigLH,oup;

                        for (unsigned int i=0;i <  len_x;i++) {

                                double temp1 = cLL[i][j];
                                double temp2 = cLH[i][j];
                                sigLL.push_back(temp1);
                                sigLH.push_back(temp2);
                        }
                        idwt1_sym_m(nm,oup,sigLL,sigLH);

                                for (int i=0;i < (int) oup.size();i++) {
                                cL[i][j] = oup[i];
                                }

                }
        } else{
                unsigned int rows1 =cLH.size();
                unsigned int cols1 =cLH[0].size();

                for (unsigned int j =0; j < cols1;j++){
                        vector<double> temp_L1,temp_L2,oup;
                        for (unsigned int i =0; i < rows1; i++){
                                double temp = cLL[i][j];
                                temp_L1.push_back(temp);

                                double temp2 = cLH[i][j];
                                temp_L2.push_back(temp2);
                        }
                        idwt1_sym_m(nm,oup,temp_L1,temp_L2);

                                for (unsigned int i =0; i < oup.size(); i++){
                                        cL[i][j]=oup[i];
                                }

               }
            }


        for (unsigned int j =0; j < len_y; j++) {
                        vector<double> sigHL,sigHH,oup2;
                        for (unsigned int i=0;i <  len_x;i++) {
                                double temp3 = cHL[i][j];
                                double temp4 = cHH[i][j];
                                sigHL.push_back(temp3);
                                sigHH.push_back(temp4);
                        }
                        idwt1_sym_m(nm,oup2,sigHL,sigHH);

                                for (int i=0;i < (int) oup2.size();i++) {
                                cH[i][j] = oup2[i];
                                }

                }

       vector<vector<double> > signal(2*len_x-lf +2,vector<double>(2 *len_y - lf +2 ));
        for (unsigned int i =0; i < 2 * len_x - lf +2; i++) {
                        vector<double> sigL,sigH,oup;
                        for (unsigned int j=0;j <  len_y;j++) {
                                double temp5 = cL[i][j];
                                double temp6 = cH[i][j];
                                sigL.push_back(temp5);
                                sigH.push_back(temp6);
                        }
                        idwt1_sym_m(nm,oup,sigL,sigH);

                                for (int j=0;j < (int) oup.size();j++) {
                                signal[i][j] = oup[j];
                                }

                }


        idwt_output = signal;



       if (iter ==0) {
           sum_coef+= 4 *rows_n * cols_n;
       } else {
           sum_coef+= 3 *rows_n * cols_n;
       }
   cLL = signal;


    }


        return 0;
}


void* dwt2_sym(string name,vector<vector<double> > &signal, vector<vector<double> >  &cLL,
                vector<vector<double> >  &cLH, vector<vector<double> >  &cHL, vector<vector<double> > &cHH){
//Analysis
    int rows = signal.size();
    int cols = signal[0].size();
    int cols_lp1 = cLL[0].size();
    int cols_hp1 = cLL[0].size();
    vector<double> lp1,hp1,lp2,hp2;
    filtcoef(name, lp1,hp1,lp2,hp2);
    vector<vector<double> > lp_dn1(rows, vector<double>( cols_lp1));
    vector<vector<double> > hp_dn1(rows, vector<double>( cols_hp1));

    // Implementing row filtering and column downsampling in each branch.
    for (int i =0; i < rows; i++) {
                vector<double> temp_row,oup_lp,oup_hp;
                for (int j=0;j <  cols;j++) {
                        double temp = signal[i][j];
                        temp_row.push_back(temp);
                }
                dwt1_sym_m(name,temp_row,oup_lp,oup_hp);

                        for (int j=0;j < (int) oup_lp.size();j++) {
                        lp_dn1[i][j] = oup_lp[j];
                        hp_dn1[i][j] = oup_hp[j];

                        }

        }


    cols =cols_lp1;
    // Implementing column filtering and row downsampling in Low Pass branch.

    for (int j =0; j < cols; j++) {
                vector<double> temp_row3,oup_lp,oup_hp;
                for (int i=0;i <  rows;i++) {
                        double temp = lp_dn1[i][j];
                        temp_row3.push_back(temp);
                }
                dwt1_sym_m(name,temp_row3,oup_lp,oup_hp);

                        for (int i=0;i < (int) oup_lp.size();i++) {
                        cLL[i][j] = oup_lp[i];
                        cLH[i][j] = oup_hp[i];

                        }


        }



    // Implementing column filtering and row downsampling in High Pass branch.

    for (int j =0; j < cols; j++) {
                vector<double> temp_row5,oup_lp,oup_hp;
                for (int i=0;i <  rows;i++) {
                        double temp = hp_dn1[i][j];
                        temp_row5.push_back(temp);
                }
                dwt1_sym_m(name,temp_row5,oup_lp,oup_hp);

                        for (int i=0;i < (int) oup_lp.size();i++) {
                        cHL[i][j] = oup_lp[i];
                        cHH[i][j] = oup_hp[i];

                        }


        }
    return 0;
}


void* dwt_2d_sym(vector<vector<double> > &origsig, int J, string nm, vector<double> &dwt_output
                , vector<double> &flag , vector<int> &length) {
// flag will contain

   vector<vector<double> >  sig = origsig;
        int rows_n = sig.size(); // No. of rows
    int cols_n = sig[0].size(); //No. of columns
    vector<vector<double> > original_copy(rows_n,vector<double>(cols_n));

        original_copy = sig;
        int Max_Iter;
                    Max_Iter = min((int) ceil(log( double(sig.size()))/log (2.0)),(int) ceil(log( double(sig[0].size()))/log (2.0)));
                    if ( Max_Iter < J) {
                        cout << J << " Iterations are not possible with signals of this dimension "  << endl;
                        exit(1);
                    }
                    vector<double> lp1,hp1,lp2,hp2;

                    flag.push_back(double(J));


                     length.insert(length.begin(),cols_n);
                     length.insert(length.begin(),rows_n);


   // Flag Values
                     /*
                        double temp = (double) (sig2.size() - sig.size()); // Number of zeropad rows
                        flag.push_back(temp);
                        double temp2 = (double) (sig2[0].size() - sig[0].size());// Number of zpad cols
                        flag.push_back(temp2);
                        flag.push_back((double) J); // Number of Iterations
                        */
                int sum_coef = 0;
                 for (int iter = 0; iter < J; iter++) {
                 filtcoef(nm,lp1,hp1,lp2,hp2);
                 unsigned int lf = lp1.size();

                        rows_n =(int) floor((double)(rows_n + lf -1)/2);
                        cols_n =(int) floor((double) (cols_n + lf -1)/2);
                        length.insert(length.begin(),cols_n);
                        length.insert(length.begin(),rows_n);

                                vector<vector<double> >  cA(rows_n, vector<double>(cols_n));
                                vector<vector<double> >  cH(rows_n, vector<double>(cols_n));
                                vector<vector<double> >  cV(rows_n, vector<double>(cols_n));
                                vector<vector<double> >  cD(rows_n, vector<double>(cols_n));

                                if (iter == 0) {
                                dwt2_sym(nm,original_copy,cA,cH,cV,cD);
                                } else {
                                        dwt2_sym(nm,original_copy,cA,cH,cV,cD);

                                }
                                vector<double>   temp_sig2;

                    original_copy = cA;
                    if (iter == J-1) {
                                for(int i =0; i < rows_n; i++){
                                        for (int j =0; j < cols_n; j++){
                                        double temp=cA[i][j];
                                        temp_sig2.push_back(temp);
                                        }
                                }
                    }
                                for(int i =0; i < rows_n; i++){
                                        for (int j = cols_n; j < cols_n * 2; j++){
                                        double temp =cH[i][j - cols_n];
                                        temp_sig2.push_back(temp);
                                        }
                                }

                                for(int i = rows_n; i < rows_n * 2; i++){
                                        for (int j =0; j < cols_n; j++){
                                        double temp=cV[i - rows_n][j];
                                        temp_sig2.push_back(temp);
                                        }
                                }

                                for(int i = rows_n; i < rows_n * 2; i++){
                                        for (int j = cols_n; j < cols_n * 2; j++){
                                        double temp =cD[i- rows_n][j - cols_n];
                                         temp_sig2.push_back(temp);
                                        }
                                }

                                dwt_output.insert(dwt_output.begin(),temp_sig2.begin(),temp_sig2.end());
                                         sum_coef += 4 * rows_n * cols_n;



                        }
/*
                        ofstream dwt2out("dwt2out.dat");
                        for (unsigned int i= 0; i < dwt_output.size(); i++){
                   dwt2out << dwt_output[i] <<endl;
                        }
                        */

        return 0;

}


void* idwt1_sym(string wname, vector<double> &X, vector<double> &app, vector<double> &detail) {

        // Not Tested. Use dwt_sym and idwt_sym for any and all computations
        vector<double> dwtop;
        vector<double> flag;
        vector<double> idwt_output;
    vector<int> length;
    length[0] = app.size();
    length[1] = detail.size();
    dwtop = app;
    dwtop.insert(dwtop.end(),detail.begin(),detail.end());
    flag.push_back(1);
    flag.push_back(0);
    idwt_sym(dwtop,flag,wname,idwt_output,length);
    X = idwt_output;

   return 0;
}

void* idwt1_sym_m(string wname, vector<double> &idwt_output, vector<double> &app, vector<double> &detail) {

    int U = 2; // Upsampling Factor
    vector<double> lpd1,hpd1, lpr1, hpr1;

    filtcoef(wname,lpd1,hpd1,lpr1,hpr1);
     int lf = lpr1.size();


            // Operations in the Low Frequency branch of the Synthesis Filter Bank
            vector<double> X_lp;
            vector<double> cA_up;
            upsamp(app, U,cA_up );
            cA_up.pop_back();
            convfftm(cA_up, lpr1, X_lp);



            // Operations in the High Frequency branch of the Synthesis Filter Bank

            vector<double> X_hp;
            vector<double> cD_up;
            upsamp(detail, U, cD_up);
            cD_up.pop_back();
            convfftm(cD_up, hpr1, X_hp);


    vecsum(X_lp,X_hp,idwt_output);

    idwt_output.erase(idwt_output.begin(),idwt_output.begin()+lf-2);
    idwt_output.erase(idwt_output.end()-(lf - 2),idwt_output.end());

   return 0;
}


void* symm_ext(vector<double> &sig, int a) {
        unsigned int len = sig.size();
        for (int i =0; i < a; i++) {
                double temp1= sig[i * 2];
                double temp2= sig[len - 1];
                sig.insert(sig.begin(),temp1);
                sig.insert(sig.end(),temp2);
        }

        return 0;

}

void* idwt_sym(vector<double> &dwtop,vector<double> &flag, string nm,
                vector<double> &idwt_output, vector<int> &length) {

        int J =(int) flag[1];
        unsigned int lf;

            vector<double> app;
            vector<double> detail;
            unsigned int app_len = length[0];
            unsigned int det_len = length[1];

        vector<double>::iterator dwt;
        dwt = dwtop.begin();
        app.assign(dwt,dwtop.begin()+app_len);
        detail.assign(dwtop.begin()+app_len, dwtop.begin()+ 2* app_len);

            for (int i = 0; i < J; i++) {

                int U = 2; // Upsampling Factor
                vector<double> lpd1,hpd1, lpr1, hpr1;

                filtcoef(nm,lpd1,hpd1,lpr1,hpr1);
                 lf = lpr1.size();


                        // Operations in the Low Frequency branch of the Synthesis Filter Bank
                        vector<double> X_lp;
                        vector<double> cA_up;
                        upsamp(app, U,cA_up );
                        cA_up.pop_back();
                        convfft(cA_up, lpr1, X_lp);



                        // Operations in the High Frequency branch of the Synthesis Filter Bank

                        vector<double> X_hp;
                        vector<double> cD_up;
                        upsamp(detail, U, cD_up);
                        cD_up.pop_back();
                        convfft(cD_up, hpr1, X_hp);


                app_len += det_len;
                vecsum(X_lp,X_hp,idwt_output);

                idwt_output.erase(idwt_output.begin(),idwt_output.begin()+lf-2);
                idwt_output.erase(idwt_output.end()-(lf - 2),idwt_output.end());

                app.clear();
                detail.clear();
                        if ( i < J - 1 ) {
                                det_len = length[i+2];
        //        detail.assign(dwtop.begin()+app_len, dwtop.begin()+ det_len);

            for (unsigned int l = 0; l < det_len;l++) {
                double temp = dwtop[app_len + l];
                detail.push_back(temp);
            }

                        }
            app = idwt_output;

            for (int iter1 = 0; iter1 < (int) (app.size() - det_len);iter1++) {
                app.pop_back();
            }

            }


            // Remove ZeroPadding

            int zerop =(int) flag[0];
            idwt_output.erase(idwt_output.end()- zerop,idwt_output.end());
            return 0;
}

void* dwt1_sym(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD) {

        vector<double> lp1, hp1, lp2, hp2;

                filtcoef(wname,lp1,hp1,lp2,hp2);
                int D = 2; // Downsampling Factor is 2
                int lf = lp1.size();
                symm_ext(signal,lf-1);

                 vector<double> cA_undec;
                     //sig value
                     convfft(signal,lp1,cA_undec);
                     cA_undec.erase(cA_undec.begin(),cA_undec.begin()+lf);
                     cA_undec.erase(cA_undec.end()-lf+1,cA_undec.end());
                     downsamp(cA_undec, D, cA);
                   //  cA.erase(cA.begin(),cA.begin()+(int) ceil(((double)lf-1.0)/2.0));
                   //  cA.erase(cA.end()-(int) ceil(((double)lf-1.0)/2.0),cA.end());


                     //High Pass Branch Computation

                     vector<double> cD_undec;
                     convfft(signal,hp1,cD_undec);
                     cD_undec.erase(cD_undec.begin(),cD_undec.begin()+lf);
                     cD_undec.erase(cD_undec.end()-lf+1,cD_undec.end());
                     downsamp(cD_undec,D,cD);
                  //   cD.erase(cD.begin(),cD.begin()+(int) ceil(((double)lf-1.0)/2.0));
                  //   cD.erase(cD.end()-(int) ceil(((double)lf-1.0)/2.0),cD.end());

                filtcoef(wname,lp1,hp1,lp2,hp2);

  return 0;
}

void* dwt1_sym_m(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD) {

        vector<double> lp1, hp1, lp2, hp2;

                filtcoef(wname,lp1,hp1,lp2,hp2);
                int D = 2; // Downsampling Factor is 2
                int lf = lp1.size();
                symm_ext(signal,lf-1);

                 vector<double> cA_undec;
                     //sig value
                     convfftm(signal,lp1,cA_undec);
                     cA_undec.erase(cA_undec.begin(),cA_undec.begin()+lf);
                     cA_undec.erase(cA_undec.end()-lf+1,cA_undec.end());
                     downsamp(cA_undec, D, cA);
                   //  cA.erase(cA.begin(),cA.begin()+(int) ceil(((double)lf-1.0)/2.0));
                   //  cA.erase(cA.end()-(int) ceil(((double)lf-1.0)/2.0),cA.end());


                     //High Pass Branch Computation

                     vector<double> cD_undec;
                     convfftm(signal,hp1,cD_undec);
                     cD_undec.erase(cD_undec.begin(),cD_undec.begin()+lf);
                     cD_undec.erase(cD_undec.end()-lf+1,cD_undec.end());
                     downsamp(cD_undec,D,cD);
                  //   cD.erase(cD.begin(),cD.begin()+(int) ceil(((double)lf-1.0)/2.0));
                  //   cD.erase(cD.end()-(int) ceil(((double)lf-1.0)/2.0),cD.end());

                filtcoef(wname,lp1,hp1,lp2,hp2);

  return 0;
}

void* dwt_sym(vector<double> &signal, int J,string nm, vector<double> &dwt_output,
                vector<double> &flag, vector<int> &length){

    unsigned int temp_len = signal.size();
        if ( (temp_len % 2) != 0) {
                double temp =signal[temp_len - 1];
                signal.push_back(temp);
                flag.push_back(1);
                temp_len++;
        } else {
                flag.push_back(0);
        }
        length.push_back(temp_len);
        flag.push_back(double(J));
        // flag[2] contains symmetric extension length


    vector<double> original_copy, appx_sig, det_sig;
    original_copy = signal;


        //  Storing Filter Values for GnuPlot
             vector<double> lp1,hp1,lp2,hp2;
             filtcoef(nm,lp1,hp1,lp2,hp2);
        for (int iter = 0; iter < J; iter++) {
            dwt1_sym(nm,signal, appx_sig, det_sig);
                dwt_output.insert(dwt_output.begin(),det_sig.begin(),det_sig.end());
            int l_temp = det_sig.size();
            length.insert(length.begin(),l_temp);

            if (iter == J-1 ) {
                dwt_output.insert(dwt_output.begin(),appx_sig.begin(),appx_sig.end());
                int l_temp = appx_sig.size();
                length.insert(length.begin(),l_temp);

            }

            signal.clear();
            signal = appx_sig;
            appx_sig.clear();
            det_sig.clear();

        }
         signal = original_copy;

return 0;
}


void* freq(vector<double> &sig, vector<double> &freq_resp) {
         unsigned int K = sig.size();
            unsigned int N = (unsigned int) pow(2.0,ceil(log10 (static_cast<double>(K))/log10(2.0)));
            vector<complex<double> > fft_oup;
            for (unsigned int i =0; i < sig.size(); i++) {
                 double temp = sig[i];
                 fft_oup.push_back(complex<double>(temp,0));
            }
            fft(fft_oup,1,N);

            for (unsigned int i = 0; i < N; i++){
                double temp = abs(fft_oup[i]);
                freq_resp.push_back(temp);
            }
            circshift(freq_resp, N/2);
            return 0;
}

double convfft(vector<double> &a, vector<double> &b, vector<double> &c) {
    fftw_complex *inp_data, *filt_data, *inp_fft, *filt_fft, *temp_data, *temp_ifft;
    fftw_plan plan_forward_inp,plan_forward_filt, plan_backward;

    unsigned int sz = a.size() + b.size() - 1;
    inp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    filt_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    inp_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    filt_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    temp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    temp_ifft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    plan_forward_inp  = fftw_plan_dft_1d( sz, inp_data, inp_fft, FFTW_FORWARD, FFTW_ESTIMATE );
    plan_forward_filt  = fftw_plan_dft_1d( sz, filt_data, filt_fft, FFTW_FORWARD, FFTW_ESTIMATE );
    plan_backward = fftw_plan_dft_1d( sz, temp_data, temp_ifft, FFTW_BACKWARD, FFTW_ESTIMATE );


    for (unsigned int i =0; i < sz; i++) {
         if (i < a.size()) {
        inp_data[i][0] = a[i];
         } else {
            inp_data[i][0] = 0.0;

         }
        inp_data[i][1] = 0.0;
        if (i < b.size()) {
                filt_data[i][0] = b[i];
           } else {
                filt_data[i][0] = 0.0;

           }
        filt_data[i][1] = 0.0;

    }



    fftw_execute(plan_forward_inp);

    fftw_execute(plan_forward_filt);

    for (unsigned int i =0; i < sz; i++){
         temp_data[i][0] = inp_fft[i][0]*filt_fft[i][0] - inp_fft[i][1]*filt_fft[i][1];

         temp_data[i][1] = inp_fft[i][0]*filt_fft[i][1] + inp_fft[i][1]*filt_fft[i][0];



    }


    fftw_execute(plan_backward);

    for (unsigned int i = 0; i < sz; i++) {
        double temp1;
        temp1 = temp_ifft[i][0] / (double) sz;
        c.push_back(temp1);

    }
    fftw_free(inp_data);
    fftw_free(filt_data);
    fftw_free(inp_fft);
    fftw_free(filt_fft);
    fftw_free(temp_data);
    fftw_free(temp_ifft);
    fftw_destroy_plan(plan_forward_inp);
    fftw_destroy_plan(plan_forward_filt);
    fftw_destroy_plan(plan_backward);

    return 0;
}

double convfftm(vector<double> &a, vector<double> &b, vector<double> &c) {
    fftw_complex *inp_data, *filt_data, *inp_fft, *filt_fft, *temp_data, *temp_ifft;

    unsigned int sz = a.size() + b.size() - 1;
    inp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    filt_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    inp_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    filt_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    temp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    temp_ifft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    if (sz != transient_size_of_fft) {

        if (transient_size_of_fft != 0) {
            fftw_destroy_plan(plan_forward_inp);
            fftw_destroy_plan(plan_forward_filt);
            fftw_destroy_plan(plan_backward);
        }

    plan_forward_inp  = fftw_plan_dft_1d( sz, inp_data, inp_fft, FFTW_FORWARD, FFTW_MEASURE );
    plan_forward_filt  = fftw_plan_dft_1d( sz, filt_data, filt_fft, FFTW_FORWARD, FFTW_MEASURE );
    plan_backward = fftw_plan_dft_1d( sz, temp_data, temp_ifft, FFTW_BACKWARD, FFTW_MEASURE );
    transient_size_of_fft = sz;

}

    for (unsigned int i =0; i < sz; i++) {
         if (i < a.size()) {
        inp_data[i][0] = a[i];
         } else {
            inp_data[i][0] = 0.0;

         }
        inp_data[i][1] = 0.0;
        if (i < b.size()) {
                filt_data[i][0] = b[i];
           } else {
                filt_data[i][0] = 0.0;

           }
        filt_data[i][1] = 0.0;

    }


    fftw_execute_dft( plan_forward_inp,inp_data, inp_fft);
    fftw_execute_dft( plan_forward_filt,filt_data, filt_fft);



    for (unsigned int i =0; i < sz; i++){
         temp_data[i][0] = inp_fft[i][0]*filt_fft[i][0] - inp_fft[i][1]*filt_fft[i][1];

         temp_data[i][1] = inp_fft[i][0]*filt_fft[i][1] + inp_fft[i][1]*filt_fft[i][0];



    }

    fftw_execute_dft( plan_backward, temp_data, temp_ifft);


    for (unsigned int i = 0; i < sz; i++) {
        double temp1;
        temp1 = temp_ifft[i][0] / (double) sz;
        c.push_back(temp1);

    }
    fftw_free(inp_data);
    fftw_free(filt_data);
    fftw_free(inp_fft);
    fftw_free(filt_fft);
    fftw_free(temp_data);
    fftw_free(temp_ifft);

    return 0;
}

void* fft(vector<complex<double> > &data, int sign,unsigned  int N){
        double pi = - 3.14159265358979;
        if ( sign == 1 || sign == -1) {
        pi = sign * pi;
        } else {
                cout << "Format fft(data, num), num = +1(fft) and num = -1 (Ifft)"  << endl;
            exit(1);
        }
        unsigned int len = data.size();
        vector<complex<double> >::iterator it;
        it = data.end();
        if ( len != N) {
                unsigned int al = N - len;
        data.insert(it,al,complex<double>(0,0));
        }

        unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
        vector<complex<double> >::iterator it1;
        it1 = data.end();
        if ( N < K) {
                unsigned int al = K - N;
                data.insert(it1,al,complex<double>(0,0));
                N = K;
        }

      bitreverse(data);

//	 radix2(data);
         for (unsigned int iter = 1; iter < N; iter <<=1)
            {
               const unsigned int step = iter << 1;

               const double theta =  pi / double(iter);

               double wtemp = sin(theta * .5);
               //   Multipliers
               double wreal = -2 * wtemp * wtemp;
               double wimag = sin(theta);

               //   Factors
               double wr = 1.0;
               double wi = 0.0;
               //   Iteration through two loops

               for (unsigned int m = 0; m < iter; m++)
               {
                  //   Iteration within m
                  for (unsigned int i = m; i < N; i += step)
                  {
                     //   jth position
                     const unsigned int j = i + iter;

                     double tempr= wr * real(data[j]) - wi * imag(data[j]);
                     double tempi= wr * imag(data[j]) + wi * real(data[j]);

                                 complex<double> temp(tempr,tempi);
                                 data[j]= data[i]- temp;
                                 data[i] += temp;

                  }
                  //   Twiddle Factors updated
              wtemp = wr;
                  wr += wr * wreal - wi * wimag;
                  wi += wi * wreal + wtemp * wimag ;
               }

            }

         if ( sign == -1) {
                 double scale = 1.0/double(N);
                 for (unsigned int i = 0; i < N; i++){
                         data[i]*=scale;
                 }
         }



        // Place holder
        return 0;
}


void* bitreverse(vector<complex<double> > &sig) {
        unsigned int len = sig.size();
        unsigned int N = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(len))/log10(2.0)));
         unsigned int rev = 0;
           //   Processing Input Data
           for (unsigned int iter = 0; iter < N; ++iter)
           {
              if (rev > iter)
              {
                 //   Replacing current values with reversed values

                 double tempr = real(sig[rev]);
                         double tempi = imag(sig[rev]);
                         complex<double> temp(tempr,tempi);
                         sig[rev] = sig[iter];
                         sig[iter] = temp;

              }
              //   Using filter "filt" such that the value of reverse changes with each iteration
              unsigned int filt = N;
              while (rev & (filt >>= 1)) {
                 rev &= ~filt;
              }
              rev |= filt;
           }
           return 0;

}


void* dwt(vector<double> &sig, int J, string nm, vector<double> &dwt_output
                , vector<double> &flag, vector<int> &length ) {

        int Max_Iter;
                    Max_Iter = (int) ceil(log( double(sig.size()))/log (2.0)) - 2;

                    if ( Max_Iter < J) {
                      J = Max_Iter;

                    }

    vector<double> original_copy,orig, appx_sig, det_sig;
    original_copy = sig;

    // Zero Pad the Signal to nearest 2^ M value ,where M is an integer.
    unsigned int temp_len = sig.size();
        if ( (temp_len % 2) != 0) {
                double temp =sig[temp_len - 1];
                sig.push_back(temp);
                flag.push_back(1);
                temp_len++;
        } else {
                flag.push_back(0);
        }
        length.push_back(temp_len);
        flag.push_back(double(J));

    orig = sig;


        //  Storing Filter Values for GnuPlot
             vector<double> lp1,hp1,lp2,hp2;
             filtcoef(nm,lp1,hp1,lp2,hp2);


    for (int iter = 0; iter < J; iter++) {
        dwt1(nm,orig, appx_sig, det_sig);
        dwt_output.insert(dwt_output.begin(),det_sig.begin(),det_sig.end());

        int l_temp = det_sig.size();
        length.insert(length.begin(),l_temp);

        if (iter == J-1 ) {
                dwt_output.insert(dwt_output.begin(),appx_sig.begin(),appx_sig.end());
                int l_temp2 = appx_sig.size();
                length.insert(length.begin(),l_temp2);

        }

        orig = appx_sig;
        appx_sig.clear();
        det_sig.clear();

    }

     sig = original_copy;
        return 0;
}


void circshift(vector<double> &sig_cir, int L){
        if ( abs(L) >(signed int) sig_cir.size()) {
                L = sign(L) * (abs(L) % sig_cir.size());
        }

        if ( L < 0 ){
                L = (sig_cir.size() + L) % sig_cir.size();
        //	cout << "L" << L << endl;
        }
                for (int i = 0; i < L; i++){
                        sig_cir.push_back(sig_cir[0]);
                        sig_cir.erase(sig_cir.begin());
                }

}

double convol(vector<double> &a1, vector<double> &b1, vector<double> &c) {
     unsigned int len_c = a1.size() + b1.size() - 1;
     vector<double> a = a1;
     vector<double> b = b1;
     double*  oup= NULL;

     oup = new double[len_c];
     vector<double>::iterator a_it;
     a_it = a.end();
     signed int al = len_c - a.size();
     a.insert(a_it,al,0);


     vector<double>::iterator b_it;
     b_it = b.end();
     signed int bl = len_c - b.size();
     b.insert(b_it,bl, 0);


     for (unsigned int ini = 0; ini < len_c ; ini++){
         double ou1 = 0;
         oup[ini] = 0;
         double temp = 0;
         for (unsigned int jni = 0; jni <= ini; jni++) {
                 ou1 = a[jni] * b[ini - jni];
                 oup[ini]+= ou1;
         }
         temp = oup[ini];
         c.push_back(temp);
     }
     delete [] oup;
     oup = NULL;
     return 0;
}

void downsamp(vector<double> &sig, int M, vector<double> &sig_d){
        int len = sig.size();
        double len_n = ceil( (double) len / (double) M);
        for (int i = 0; i < (int) len_n; i++) {
                double temp = sig[i*M];
                sig_d.push_back(temp);
        }
}





void* dwt1(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD) {

        vector<double> lpd, hpd, lpr, hpr;

                filtcoef(wname,lpd,hpd,lpr,hpr);

                int len_lpfilt = lpd.size();
                int len_hpfilt = hpd.size();
                int len_avg = (len_lpfilt + len_hpfilt) / 2;
                int len_sig = 2 * (int) ceil((double) signal.size() / 2.0);

                // cout << len_lpfilt << "Filter" << endl;
                per_ext(signal,len_avg / 2); // Periodic Extension
                // computations designed to deal with boundary distortions

                // Low Pass Filtering Operations in the Analysis Filter Bank Section
//		int len_cA =(int)  floor(double (len_sig + len_lpfilt -1) / double (2));
                vector<double> cA_undec;
                // convolving signal with lpd, Low Pass Filter, and O/P is stored in cA_undec
                convfft(signal,lpd,cA_undec);
                int D = 2; // Downsampling Factor is 2
                cA_undec.erase(cA_undec.begin(),cA_undec.begin()+len_avg-1);
                cA_undec.erase(cA_undec.end()-len_avg+1,cA_undec.end());
                cA_undec.erase(cA_undec.begin()+len_sig,cA_undec.end());
                cA_undec.erase(cA_undec.begin());


                // Downsampling by 2 gives cA
                downsamp(cA_undec, D, cA);

          //     cA.erase(cA.begin(),cA.begin()+len_avg/2);
          //      cA.erase(cA.end()-len_avg/2,cA.end());

                // High Pass Filtering Operations in the Analysis Filter Bank Section
//		int len_cA =(int)  floor(double (len_sig + len_lpfilt -1) / double (2));

                vector<double> cD_undec;
                // convolving signal with lpd, Low Pass Filter, and O/P is stored in cA_undec
                convfft(signal,hpd,cD_undec);
                cD_undec.erase(cD_undec.begin(),cD_undec.begin()+len_avg-1);
                cD_undec.erase(cD_undec.end()-len_avg+1,cD_undec.end());
                cD_undec.erase(cD_undec.begin()+len_sig,cD_undec.end());
                cD_undec.erase(cD_undec.begin());

                 // Downsampling Factor is 2

                // Downsampling by 2 gives cA
                downsamp(cD_undec, D, cD);

    //            cD.erase(cD.begin(),cD.begin()+len_avg/2);
      //          cD.erase(cD.end()-len_avg/2,cD.end());

                filtcoef(wname,lpd,hpd,lpr,hpr);

  return 0;
}

void* dwt1_m(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD) {

        vector<double> lpd, hpd, lpr, hpr;

                filtcoef(wname,lpd,hpd,lpr,hpr);

                int len_lpfilt = lpd.size();
                int len_hpfilt = hpd.size();
                int len_avg = (len_lpfilt + len_hpfilt) / 2;
                int len_sig = 2 * (int) ceil((double) signal.size() / 2.0);

                // cout << len_lpfilt << "Filter" << endl;
                per_ext(signal,len_avg / 2); // Periodic Extension
                // computations designed to deal with boundary distortions

                // Low Pass Filtering Operations in the Analysis Filter Bank Section
//		int len_cA =(int)  floor(double (len_sig + len_lpfilt -1) / double (2));
                vector<double> cA_undec;
                // convolving signal with lpd, Low Pass Filter, and O/P is stored in cA_undec
                convfftm(signal,lpd,cA_undec);
                int D = 2; // Downsampling Factor is 2
                cA_undec.erase(cA_undec.begin(),cA_undec.begin()+len_avg-1);
                cA_undec.erase(cA_undec.end()-len_avg+1,cA_undec.end());
                cA_undec.erase(cA_undec.begin()+len_sig,cA_undec.end());
                cA_undec.erase(cA_undec.begin());


                // Downsampling by 2 gives cA
                downsamp(cA_undec, D, cA);

          //     cA.erase(cA.begin(),cA.begin()+len_avg/2);
          //      cA.erase(cA.end()-len_avg/2,cA.end());

                // High Pass Filtering Operations in the Analysis Filter Bank Section
//		int len_cA =(int)  floor(double (len_sig + len_lpfilt -1) / double (2));

                vector<double> cD_undec;
                // convolving signal with lpd, Low Pass Filter, and O/P is stored in cA_undec
                convfftm(signal,hpd,cD_undec);
                cD_undec.erase(cD_undec.begin(),cD_undec.begin()+len_avg-1);
                cD_undec.erase(cD_undec.end()-len_avg+1,cD_undec.end());
                cD_undec.erase(cD_undec.begin()+len_sig,cD_undec.end());
                cD_undec.erase(cD_undec.begin());

                 // Downsampling Factor is 2

                // Downsampling by 2 gives cA
                downsamp(cD_undec, D, cD);

    //            cD.erase(cD.begin(),cD.begin()+len_avg/2);
      //          cD.erase(cD.end()-len_avg/2,cD.end());

                filtcoef(wname,lpd,hpd,lpr,hpr);

  return 0;
}


void* dyadic_zpad_1d(vector<double> &signal) {
        unsigned int N = signal.size();
        double M = log10 (static_cast<double> (N)) / log10(2.0);
        int D = (int) ceil(M);
        double int_val = pow(2.0, double(D)) - pow(2.0, M);

        int z = (int) int_val;
    vector<double>::iterator a_it;
    a_it = signal.end();
    double val = signal[N-1];
  //  double val = 0;
    signal.insert(a_it,z,val);
    return 0;

}


void* idwt(vector<double> &dwtop,vector<double> &flag, string nm,
                vector<double> &idwt_output,vector<int> &length) {

        int J =(int) flag[1];
 //       int zpad =(int) flag[0];


            vector<double> app;
            vector<double> detail;
            unsigned int app_len = length[0];
            unsigned int det_len = length[1];

        vector<double>::iterator dwt;
        dwt = dwtop.begin();
        app.assign(dwt,dwtop.begin()+app_len);
        detail.assign(dwtop.begin()+app_len, dwtop.begin()+ 2* app_len);

            for (int i = 0; i < J; i++) {

                idwt1(nm,idwt_output, app,detail);
                app_len +=det_len;
                app.clear();
                detail.clear();
                        if ( i < J - 1 ) {
                            det_len = length[i+2];
                            for (unsigned int l = 0; l < det_len;l++) {

                                double temp = dwtop[app_len + l];
                                detail.push_back(temp);
                            }
                            app = idwt_output;

                            if (app.size() >= detail.size()){
                                int t = app.size() - detail.size();
                               int lent = (int) floor((double)t/2.0);
                               app.erase(app.begin()+detail.size()+lent,app.end());
                               app.erase(app.begin(),app.begin()+lent);
                            }
                        }




//            int value1 = (int) ceil((double)(app.size() - det_len)/2.0);
//           int value2 = (int) floor((double)(app.size() - det_len)/2.0);

//            app.erase(app.end() -value2,app.end());
//           app.erase(app.begin(),app.begin()+value1);

            }


            // Remove ZeroPadding

            int zerop =(int) flag[0];
            idwt_output.erase(idwt_output.end()- zerop,idwt_output.end());

            return 0;
}

void* idwt1_m(string wname, vector<double> &X, vector<double> &cA, vector<double> &cD) {
        vector<double> lpd1,hpd1, lpr1, hpr1;

        filtcoef(wname,lpd1,hpd1,lpr1,hpr1);
        int len_lpfilt = lpr1.size();
        int len_hpfilt = hpr1.size();
        int len_avg = (len_lpfilt + len_hpfilt) / 2;
        unsigned int N = 2 * cD.size();
        int U = 2; // Upsampling Factor

        // Operations in the Low Frequency branch of the Synthesis Filter Bank

        vector<double> cA_up;
        vector<double> X_lp;
       // int len1 = cA_up.size();
        upsamp(cA, U, cA_up);

        per_ext(cA_up,len_avg/2);


        convfftm(cA_up, lpr1, X_lp);


        // Operations in the High Frequency branch of the Synthesis Filter Bank

        vector<double> cD_up;
        vector<double> X_hp;
        upsamp(cD, U, cD_up);
        per_ext(cD_up,len_avg/2);


        convfftm(cD_up, hpr1, X_hp);

   // Remove periodic extension

     //   X.erase(X.begin(),X.begin()+len_avg+len_avg/2-1);
     //   X.erase(X.end()-len_avg-len_avg/2,X.end());


        X_lp.erase(X_lp.begin()+N+len_avg-1,X_lp.end());
        X_lp.erase(X_lp.begin(),X_lp.begin()+len_avg-1);

        X_hp.erase(X_hp.begin()+N+len_avg-1,X_hp.end());
        X_hp.erase(X_hp.begin(),X_hp.begin()+len_avg-1);


        vecsum(X_lp,X_hp,X);


    return 0;
}

void* idwt1(string wname, vector<double> &X, vector<double> &cA, vector<double> &cD) {
        vector<double> lpd1,hpd1, lpr1, hpr1;

        filtcoef(wname,lpd1,hpd1,lpr1,hpr1);
        int len_lpfilt = lpr1.size();
        int len_hpfilt = hpr1.size();
        int len_avg = (len_lpfilt + len_hpfilt) / 2;
        unsigned int N = 2 * cD.size();
        int U = 2; // Upsampling Factor

        // Operations in the Low Frequency branch of the Synthesis Filter Bank

        vector<double> cA_up;
        vector<double> X_lp;
       // int len1 = cA_up.size();
        upsamp(cA, U, cA_up);

        per_ext(cA_up,len_avg/2);


        convfft(cA_up, lpr1, X_lp);


        // Operations in the High Frequency branch of the Synthesis Filter Bank

        vector<double> cD_up;
        vector<double> X_hp;
        upsamp(cD, U, cD_up);
        per_ext(cD_up,len_avg/2);


        convfft(cD_up, hpr1, X_hp);

   // Remove periodic extension

     //   X.erase(X.begin(),X.begin()+len_avg+len_avg/2-1);
     //   X.erase(X.end()-len_avg-len_avg/2,X.end());


        X_lp.erase(X_lp.begin()+N+len_avg-1,X_lp.end());
        X_lp.erase(X_lp.begin(),X_lp.begin()+len_avg-1);

        X_hp.erase(X_hp.begin()+N+len_avg-1,X_hp.end());
        X_hp.erase(X_hp.begin(),X_hp.begin()+len_avg-1);


        vecsum(X_lp,X_hp,X);


    return 0;
}

int sign(int X) {
        if (X >= 0)
                return 1;
        else
                return -1;
}

void upsamp(vector<double> &sig, int M, vector<double> &sig_u) {
        int len = sig.size();
        double len_n = ceil( (double) len * (double) M);

        for (int i = 0; i < (int) len_n; i++) {
                if ( i % M == 0) {
                                double temp = sig[i / M];
                                        sig_u.push_back(temp);

                                }
                                else
                                {
                                         sig_u.push_back(0);
                                }

        }



}

double op_sum(double i, double j) {
        return (i+j);
}

int vecsum(vector<double> &a, vector<double> &b, vector<double> &c){


    c.resize(a.size());
        transform (a.begin(), a.end(), b.begin(), b.begin(), op_sum);
        c = b;
                return 0;
}

void* getcoeff2d(vector<vector<double> > &dwtoutput, vector<vector<double> > &cH,
                vector<vector<double> > &cV,vector<vector<double> > &cD,vector<double> &flag, int &N) {
        if (N > flag[2]) {
                cout << "Signal is decimated only up to " << flag[2] << " levels" << endl;
                exit(1);
        }
        int rows = dwtoutput.size();
        int cols = dwtoutput[0].size();
    // Getting Horizontal Coefficients
        int r = (int) ceil((double) rows /pow(2.0,N)) ;
        int c = (int) ceil((double) cols /pow(2.0,N)) ;

        for (int i =0; i < (int) ceil ((double) rows /pow(2.0,N)); i++){
                for (int j =0; j < (int) ceil ((double) cols /pow(2.0,N)); j++) {
                        cH[i][j]=dwtoutput[i][c+ j];
                }
        }


        for (int i =0; i < (int) ceil ((double) rows /pow(2.0,N)); i++){
                for (int j =0; j < (int) ceil ((double) cols /pow(2.0,N)); j++) {
                        cV[i][j]=dwtoutput[i + r][j];
                }
        }

        for (int i =0; i < (int) ceil ((double) rows /pow(2.0,N)); i++){
                for (int j =0; j < (int) ceil ((double) cols /pow(2.0,N)); j++) {
                        cD[i][j]=dwtoutput[i + r][c+ j];
                }
        }

  return 0;
}

void* zero_remove(vector<vector<double> > &input,vector<vector<double> > &output) {
        int zero_rows = output.size()-input.size();
        int zero_cols = output[0].size()-input[0].size();

        vector<vector<double> >::iterator row = output.end()-zero_rows;


   unsigned int ousize = output.size();
        for (unsigned int i = input.size(); i < ousize; i++){
    output.erase(row);
    row++;

        }

//	unsigned int ousize2 = output[0].size();


        for (unsigned int i = 0; i < ousize; i++){
    vector<double> ::iterator col = output[i].end()-zero_cols;

    output[i].erase(col, output[i].end());

        }
        return 0;
}

void* dwt_output_dim(vector<vector<double> >&signal, int &r, int &c ){
           int rows =signal.size();
           int cols = signal[0].size();

                double Mr = log10 (static_cast<double> (rows)) / log10(2.0);
                int Dr = (int) ceil(Mr);
                double int_val_row = pow(2.0, double(Dr));
                int r1 = (int) int_val_row;

                double Mc = log10 (static_cast<double> (cols)) / log10(2.0);
                int Dc = (int) ceil(Mc);
                double int_val_cols = pow(2.0, double(Dc));
        int c1 = (int) int_val_cols;
        r=max(r1,c1);
        c=max(r1,c1);

        return 0;

}

void* dyadic_zpad_2d(vector<vector<double> > &signal,vector<vector<double> > &mod){
   int rows =signal.size();
   int cols = signal[0].size();

        for (int i=0; i < rows; i++) {
                for (int j = 0; j < cols; j++){
                        mod[i][j] = signal[i][j];
                }

        }
   // Zeropadding the columns

        double Mr = log10 (static_cast<double> (rows)) / log10(2.0);
        int Dr = (int) ceil(Mr);
        double int_val_row = pow(2.0, double(Dr)) - pow(2.0, Mr);

        int zeros_row = (int) int_val_row;

        double Mc = log10 (static_cast<double> (cols)) / log10(2.0);
        int Dc = (int) ceil(Mc);
        double int_val_cols = pow(2.0, double(Dc)) - pow(2.0, Mc);

        int zeros_cols = (int) int_val_cols;

        for (int i=0; i < rows + zeros_row; i++) {
                for (int j = cols; j < cols+zeros_cols; j++){

                        mod[i][j] = 0;
                }

        }

        for (int i= rows; i < rows + zeros_row; i++) {
                for (int j = 0; j < cols+zeros_cols; j++){
                        mod[i][j] = 0;
                }

        }

        return 0;

}

void* idwt_2d(vector<double>  &dwtop,vector<double> &flag, string nm,
                vector<vector<double> > &idwt_output, vector<int> &length){
    int J =(int) flag[0];
    int rows =length[0];
    int cols =length[1];

    int sum_coef =0;
    vector<double> lp1,hp1,lp2,hp2;
    filtcoef(nm,lp1,hp1,lp2,hp2);
        vector<vector<double> >  cLL(rows, vector<double>(cols));


    for (int iter=0; iter < J; iter++) {

        int rows_n = length[2*iter];
        int cols_n = length[2*iter + 1];

        vector<vector<double> >  cLH(rows_n, vector<double>(cols_n));
        vector<vector<double> >  cHL(rows_n, vector<double>(cols_n));
        vector<vector<double> >  cHH(rows_n, vector<double>(cols_n));

        for (int i = 0 ; i < rows_n; i++ ){
                for (int j = 0; j < cols_n; j++){
                        if (iter == 0) {
                        cLL[i][j] = dwtop[sum_coef+ i * cols_n + j];

                cLH[i][j] = dwtop[sum_coef+ rows_n * cols_n+ i * cols_n + j];

                cHL[i][j] = dwtop[sum_coef+ 2 * rows_n * cols_n + i * cols_n + j];

                cHH[i][j] = dwtop[sum_coef+ 3* rows_n * cols_n + i * cols_n + j];
                        } else {

                        cLH[i][j] = dwtop[sum_coef+  i * cols_n + j];

                        cHL[i][j] = dwtop[sum_coef+ rows_n * cols_n + i * cols_n + j];

                    cHH[i][j] = dwtop[sum_coef+ 2* rows_n * cols_n + i * cols_n + j];

                        }
                }
        }


  //      temp_A = cLL;
  //  	idwt2_sym(nm,idwt_output2, cA, cH,cV,cD);

        unsigned int len_x = cLH.size();
        unsigned int len_y = cLH[0].size();

        // Row Upsampling and Column Filtering at the first LP Stage
        vector<vector<double> > cL(2 *len_x,vector<double>(len_y ));
        vector<vector<double> > cH(2 * len_x ,vector<double>(len_y ));

        if (iter ==0) {
        for (unsigned int j =0; j < len_y; j++) {

                        vector<double> sigLL,sigLH,oup;

                        for (unsigned int i=0;i <  len_x;i++) {

                                double temp1 = cLL[i][j];
                                double temp2 = cLH[i][j];
                                sigLL.push_back(temp1);
                                sigLH.push_back(temp2);
                        }
                        idwt1_m(nm,oup,sigLL,sigLH);

                                for (int i=0;i < (int) oup.size();i++) {
                                cL[i][j] = oup[i];
                                }

                }
        } else{
                unsigned int rows1 =cLH.size();
                unsigned int cols1 =cLH[0].size();

                for (unsigned int j =0; j < cols1;j++){
                        vector<double> temp_L1,temp_L2,oup;
                        for (unsigned int i =0; i < rows1; i++){
                                double temp = cLL[i][j];
                                temp_L1.push_back(temp);

                                double temp2 = cLH[i][j];
                                temp_L2.push_back(temp2);
                        }
                        idwt1_m(nm,oup,temp_L1,temp_L2);

                                for (unsigned int i =0; i < oup.size(); i++){
                                        cL[i][j]=oup[i];
                                }

               }
            }


        for (unsigned int j =0; j < len_y; j++) {
                        vector<double> sigHL,sigHH,oup2;
                        for (unsigned int i=0;i <  len_x;i++) {
                                double temp3 = cHL[i][j];
                                double temp4 = cHH[i][j];
                                sigHL.push_back(temp3);
                                sigHH.push_back(temp4);
                        }
                        idwt1_m(nm,oup2,sigHL,sigHH);

                                for (int i=0;i < (int) oup2.size();i++) {
                                cH[i][j] = oup2[i];
                                }

                }

       vector<vector<double> > signal(2*len_x,vector<double>(2 *len_y ));
        for (unsigned int i =0; i < 2 * len_x ; i++) {
                        vector<double> sigL,sigH,oup;
                        for (unsigned int j=0;j <  len_y;j++) {
                                double temp5 = cL[i][j];
                                double temp6 = cH[i][j];
                                sigL.push_back(temp5);
                                sigH.push_back(temp6);
                        }
                        idwt1_m(nm,oup,sigL,sigH);

                                for (int j=0;j < (int) oup.size();j++) {
                                signal[i][j] = oup[j];
                                }

                }

        idwt_output = signal;



       if (iter ==0) {
           sum_coef+= 4 *rows_n * cols_n;
       } else {
           sum_coef+= 3 *rows_n * cols_n;
       }
   cLL = signal;


    }


        return 0;
}



void* dwt_2d(vector<vector<double> > &origsig, int J, string nm, vector<double> &dwt_output
                , vector<double> &flag , vector<int> &length) {
// flag will contain

   vector<vector<double> >  sig = origsig;
        int rows_n = sig.size(); // No. of rows
    int cols_n = sig[0].size(); //No. of columns
    vector<vector<double> > original_copy(rows_n,vector<double>(cols_n));

        original_copy = sig;
        int Max_Iter;
                    Max_Iter = min((int) ceil(log( double(sig.size()))/log (2.0)),(int) ceil(log( double(sig[0].size()))/log (2.0)));
                    if ( Max_Iter < J) {
                        cout << J << " Iterations are not possible with signals of this dimension "  << endl;
                        exit(1);
                    }
                    vector<double> lp1,hp1,lp2,hp2;

                    flag.push_back(double(J));
                    flag.push_back(0);


                     length.insert(length.begin(),cols_n);
                     length.insert(length.begin(),rows_n);


                int sum_coef = 0;
                 for (int iter = 0; iter < J; iter++) {
                 filtcoef(nm,lp1,hp1,lp2,hp2);

                        rows_n =(int) ceil((double)rows_n /2.0);
                        cols_n =(int) ceil((double) cols_n/2.0);
                        length.insert(length.begin(),cols_n);
                        length.insert(length.begin(),rows_n);

                                vector<vector<double> >  cA(rows_n, vector<double>(cols_n));
                                vector<vector<double> >  cH(rows_n, vector<double>(cols_n));
                                vector<vector<double> >  cV(rows_n, vector<double>(cols_n));
                                vector<vector<double> >  cD(rows_n, vector<double>(cols_n));

                                if (iter == 0) {
                                dwt2(nm,original_copy,cA,cH,cV,cD);
                                } else {
                                        dwt2(nm,original_copy,cA,cH,cV,cD);

                                }
                                vector<double>   temp_sig2;

                    original_copy = cA;
                    if (iter == J-1) {
                                for(int i =0; i < rows_n; i++){
                                        for (int j =0; j < cols_n; j++){
                                        double temp=cA[i][j];
                                        temp_sig2.push_back(temp);
                                        }
                                }
                    }
                                for(int i =0; i < rows_n; i++){
                                        for (int j = cols_n; j < cols_n * 2; j++){
                                        double temp =cH[i][j - cols_n];
                                        temp_sig2.push_back(temp);
                                        }
                                }

                                for(int i = rows_n; i < rows_n * 2; i++){
                                        for (int j =0; j < cols_n; j++){
                                        double temp=cV[i - rows_n][j];
                                        temp_sig2.push_back(temp);
                                        }
                                }

                                for(int i = rows_n; i < rows_n * 2; i++){
                                        for (int j = cols_n; j < cols_n * 2; j++){
                                        double temp =cD[i- rows_n][j - cols_n];
                                         temp_sig2.push_back(temp);
                                        }
                                }

                                dwt_output.insert(dwt_output.begin(),temp_sig2.begin(),temp_sig2.end());
                                         sum_coef += 4 * rows_n * cols_n;



                        }
/*
                        ofstream dwt2out("dwt2out.dat");
                        for (unsigned int i= 0; i < dwt_output.size(); i++){
                   dwt2out << dwt_output[i] <<endl;
                        }
                        */

        return 0;

}



void* branch_lp_hp_up(string wname,vector<double> &cA, vector<double> &cD, vector<double> &X) {
        vector<double> lpd1,hpd1, lpr1, hpr1;

                filtcoef(wname,lpd1,hpd1,lpr1,hpr1);
                int len_lpfilt = lpr1.size();
                int len_hpfilt = hpr1.size();
                int len_avg = (len_lpfilt + len_hpfilt) / 2;
                //unsigned int N = 2 * cA.size();
                int U = 2; // Upsampling Factor

                // Operations in the Low Frequency branch of the Synthesis Filter Bank

                vector<double> cA_up;
                vector<double> X_lp;
                per_ext(cA,len_avg/2);

                upsamp(cA, U, cA_up);
                convfftm(cA_up, lpr1, X_lp);


                // Operations in the High Frequency branch of the Synthesis Filter Bank

                vector<double> cD_up;
                vector<double> X_hp;
                per_ext(cD,len_avg/2);
                upsamp(cD, U, cD_up);
                convfftm(cD_up, hpr1, X_hp);




                vecsum(X_lp,X_hp,X);
           // Remove periodic extension

                X.erase(X.begin(),X.begin()+len_avg+len_avg/2-1);
                X.erase(X.end()-len_avg-len_avg/2,X.end());

            return 0;
}

void* branch_hp_dn(string wname, vector<double> &signal, vector<double> &sigop) {

        vector<double> lpd, hpd, lpr, hpr;

                        filtcoef(wname,lpd,hpd,lpr,hpr);
        //		for (unsigned int i = 0;  i < signal.size(); i++) {
        //		cout << signal[i] <<  endl;
        //		out2 << signal[i] <<endl;
        //		}

                        unsigned int temp_len = signal.size();

                        if ( (temp_len % 2) != 0) {
                                double temp =signal[temp_len - 1];
                                signal.push_back(temp);
                        }

                        int len_lpfilt = lpd.size();
                        int len_hpfilt = hpd.size();
                        int len_avg = (len_lpfilt + len_hpfilt) / 2;

                        // cout << len_lpfilt << "Filter" << endl;
                        per_ext(signal,len_avg / 2); // Periodic Extension
                        // computations designed to deal with boundary distortions

                        // Low Pass Filtering Operations in the Analysis Filter Bank Section
        //		int len_cA =(int)  floor(double (len_sig + len_lpfilt -1) / double (2));
                        vector<double> cA_undec;
                        // convolving signal with lpd, Low Pass Filter, and O/P is stored in cA_undec
                        convfftm(signal,hpd,cA_undec);
                        int D = 2; // Downsampling Factor is 2

                        // Downsampling by 2 gives cA
                        downsamp(cA_undec, D, sigop);

                        sigop.erase(sigop.begin(),sigop.begin()+len_avg/2);
                        sigop.erase(sigop.end()-len_avg/2,sigop.end());
return 0;



}
void* branch_lp_dn(string wname, vector<double> &signal, vector<double> &sigop){

        vector<double> lpd, hpd, lpr, hpr;

                        filtcoef(wname,lpd,hpd,lpr,hpr);
        //		for (unsigned int i = 0;  i < signal.size(); i++) {
        //		cout << signal[i] <<  endl;
        //		out2 << signal[i] <<endl;
        //		}

                        unsigned int temp_len = signal.size();

                        if ( (temp_len % 2) != 0) {
                                double temp =signal[temp_len - 1];
                                signal.push_back(temp);
                        }

                        int len_lpfilt = lpd.size();
                        int len_hpfilt = hpd.size();
                        int len_avg = (len_lpfilt + len_hpfilt) / 2;

                        // cout << len_lpfilt << "Filter" << endl;
                        per_ext(signal,len_avg / 2); // Periodic Extension
                        // computations designed to deal with boundary distortions

                        // Low Pass Filtering Operations in the Analysis Filter Bank Section
        //		int len_cA =(int)  floor(double (len_sig + len_lpfilt -1) / double (2));
                        vector<double> cA_undec;
                        // convolving signal with lpd, Low Pass Filter, and O/P is stored in cA_undec
                        convfftm(signal,lpd,cA_undec);
                        int D = 2; // Downsampling Factor is 2

                        // Downsampling by 2 gives cA
                        downsamp(cA_undec, D, sigop);

                        sigop.erase(sigop.begin(),sigop.begin()+len_avg/2);
                        sigop.erase(sigop.end()-len_avg/2,sigop.end());


return 0;

}

void* idwt2(string name,vector<vector<double> > &signal, vector<vector<double> >  &cLL,
                vector<vector<double> >  &cLH, vector<vector<double> >  &cHL, vector<vector<double> > &cHH) {
// Synthesis
        int rows= cLL.size();
        int cols= cLL[0].size();
        int rows_n = 2 * rows;
        // Row Upsampling and Column Filtering at the first LP Stage
        vector<vector<double> > cL(rows_n,vector<double>(cols));
        vector<vector<double> > cH(rows_n,vector<double>(cols));

    for (int j =0; j < cols; j++) {

                vector<double> sigLL;
                vector<double> sigLH;
                for (int i=0;i <  rows;i++) {

                        double temp1 = cLL[i][j];
                        double temp2 = cLH[i][j];
                        sigLL.push_back(temp1);
                        sigLH.push_back(temp2);
                }
                        vector<double> oup;

                        branch_lp_hp_up(name,sigLL,sigLH,oup);
                        sigLL.clear();
                        sigLH.clear();
                        for (int i=0;i < (int) oup.size();i++) {
                        cL[i][j] = oup[i];
                        }

        }

    for (int j =0; j < cols; j++) {
                vector<double> sigHL;
                vector<double> sigHH;
                for (int i=0;i <  rows;i++) {
                        double temp3 = cHL[i][j];
                        double temp4 = cHH[i][j];
                        sigHL.push_back(temp3);
                        sigHH.push_back(temp4);
                }
                        vector<double> oup2;
                        branch_lp_hp_up(name,sigHL,sigHH,oup2);
                        sigHL.clear();
                        sigHH.clear();

                        for (int i=0;i < (int) oup2.size();i++) {
                        cH[i][j] = oup2[i];
                        }

        }

    for (int i =0; i < rows_n; i++) {
                vector<double> sigL;
                vector<double> sigH;
                for (int j=0;j <  cols;j++) {
                        double temp5 = cL[i][j];\
                        double temp6 = cH[i][j];
                        sigL.push_back(temp5);
                        sigH.push_back(temp6);
                }
                        vector<double> oup3;
                        branch_lp_hp_up(name,sigL,sigH,oup3);
                        sigL.clear();
                        sigH.clear();

                        for (int j=0;j < (int) oup3.size();j++) {
                        signal[i][j] = oup3[j];
                        }

        }
        return 0;
}

void* dwt2(string name,vector<vector<double> > &signal, vector<vector<double> >  &cLL,
                vector<vector<double> >  &cLH, vector<vector<double> >  &cHL, vector<vector<double> > &cHH){
//Analysis
    int rows = signal.size();
    int cols = signal[0].size();
    int cols_lp1 = cLL[0].size();
    int cols_hp1 = cLL[0].size();
    vector<double> lp1,hp1,lp2,hp2;
    filtcoef(name, lp1,hp1,lp2,hp2);
    vector<vector<double> > lp_dn1(rows, vector<double>( cols_lp1));
    vector<vector<double> > hp_dn1(rows, vector<double>( cols_hp1));

    // Implementing row filtering and column downsampling in each branch.
    for (int i =0; i < rows; i++) {
                vector<double> temp_row,oup_lp,oup_hp;
                for (int j=0;j <  cols;j++) {
                        double temp = signal[i][j];
                        temp_row.push_back(temp);
                }
                dwt1_m(name,temp_row,oup_lp,oup_hp);

                        for (int j=0;j < (int) oup_lp.size();j++) {
                        lp_dn1[i][j] = oup_lp[j];
                        hp_dn1[i][j] = oup_hp[j];

                        }

        }


    cols =cols_lp1;
    // Implementing column filtering and row downsampling in Low Pass branch.

    for (int j =0; j < cols; j++) {
                vector<double> temp_row3,oup_lp,oup_hp;
                for (int i=0;i <  rows;i++) {
                        double temp = lp_dn1[i][j];
                        temp_row3.push_back(temp);
                }
                dwt1_m(name,temp_row3,oup_lp,oup_hp);

                        for (int i=0;i < (int) oup_lp.size();i++) {
                        cLL[i][j] = oup_lp[i];
                        cLH[i][j] = oup_hp[i];

                        }


        }



    // Implementing column filtering and row downsampling in High Pass branch.

    for (int j =0; j < cols; j++) {
                vector<double> temp_row5,oup_lp,oup_hp;
                for (int i=0;i <  rows;i++) {
                        double temp = hp_dn1[i][j];
                        temp_row5.push_back(temp);
                }
                dwt1_m(name,temp_row5,oup_lp,oup_hp);

                        for (int i=0;i < (int) oup_lp.size();i++) {
                        cHL[i][j] = oup_lp[i];
                        cHH[i][j] = oup_hp[i];

                        }


        }
    return 0;
}


void* downsamp2(vector<vector<double> > & vec1,vector<vector<double> > & vec2, int rows_dn, int cols_dn) {

        int rows = vec1.size();
        int cols = vec1[0].size();
        double rows_n = ceil( (double) rows / (double) rows_dn);
        double cols_n = ceil( (double) cols / (double) cols_dn);
        for (int i =0; i < (int)rows_n; i++){
                for (int j = 0; j< (int) cols_n; j++){

  vec2[i][j] = vec1[i * rows_dn][j*cols_dn];
                }
        }

        return 0;
}

void* upsamp2(vector<vector<double> > & vec1,vector<vector<double> > & vec2, int rows_up, int cols_up){

        int rows = vec1.size();
        int cols = vec1[0].size();
        int rows_n = rows *  rows_up;
        int cols_n =  cols *  cols_up;
        for (int i = 0; i < rows_n; i++){
                for (int j = 0; j < cols_n; j++){
                        if ( i % rows_up == 0 && j % cols_up == 0){
         vec2[i][j]=vec1[(int) (i/rows_up)][(int) (j/cols_up)];
        } else {
         vec2[i][j] = 0;
        }
}
}
 return 0;
}


int filtcoef(string name, vector<double> &lp1, vector<double> &hp1, vector<double> &lp2,
                vector<double> &hp2){
    if (name == "haar" || name == "db1" ) {
        lp1.push_back(0.7071);lp1.push_back(0.7071);
        hp1.push_back(-0.7071);hp1.push_back(0.7071);
        lp2.push_back(0.7071);lp2.push_back(0.7071);
        hp2.push_back(0.7071);hp2.push_back(-0.7071);
  //  	cout << lp2[1] << endl;
//    	hpd = {-0.7071, 0.7071};
//    	lpr = {0.7071, 0.7071};
//    	hpr = {0.7071, -0.7071};
     return 0;
    }
    else if ( name == "db2"){
        double lp1_a[] = {-0.12940952255092145, 0.22414386804185735, 0.83651630373746899,
                        0.48296291314469025};
        lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        double hp1_a[] = {-0.48296291314469025, 0.83651630373746899, -0.22414386804185735,
                        -0.12940952255092145};
         hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        double lp2_a[] = {0.48296291314469025, 0.83651630373746899, 0.22414386804185735,
                        -0.12940952255092145};
         lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

         double hp2_a[] = {-0.12940952255092145, -0.22414386804185735, 0.83651630373746899,
                         -0.48296291314469025};
         hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
         return 0;
    }

    else if ( name == "db3"){
        double lp1_a[] = {0.035226291882100656, -0.085441273882241486, -0.13501102001039084,
                        0.45987750211933132, 0.80689150931333875, 0.33267055295095688};
        lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        double hp1_a[] = {-0.33267055295095688, 0.80689150931333875, -0.45987750211933132,
                        -0.13501102001039084, 0.085441273882241486, 0.035226291882100656 };
         hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        double lp2_a[] = {0.33267055295095688, 0.80689150931333875, 0.45987750211933132,
                        -0.13501102001039084, -0.085441273882241486, 0.035226291882100656 };
         lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

         double hp2_a[] = {0.035226291882100656, 0.085441273882241486, -0.13501102001039084,
                         -0.45987750211933132, 0.80689150931333875, -0.33267055295095688 };
         hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
         return 0;
    }

    else if ( name == "db4"){
        double lp1_a[] = {-0.010597401784997278, 0.032883011666982945, 0.030841381835986965,
                        -0.18703481171888114, -0.027983769416983849, 0.63088076792959036,
                        0.71484657055254153, 0.23037781330885523 };
        lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        double hp1_a[] = {-0.23037781330885523, 0.71484657055254153, -0.63088076792959036,
                        -0.027983769416983849, 0.18703481171888114, 0.030841381835986965,
                        -0.032883011666982945, -0.010597401784997278 };
         hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        double lp2_a[] = {0.23037781330885523, 0.71484657055254153, 0.63088076792959036,
                        -0.027983769416983849, -0.18703481171888114, 0.030841381835986965,
                        0.032883011666982945, -0.010597401784997278 };
         lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

         double hp2_a[] = {-0.010597401784997278, -0.032883011666982945, 0.030841381835986965,
                         0.18703481171888114, -0.027983769416983849, -0.63088076792959036,
                         0.71484657055254153, -0.23037781330885523 };
         hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
         return 0;
    }

    else if ( name == "db5"){
        double lp1_a[] = {0.0033357252850015492, -0.012580751999015526, -0.0062414902130117052,
                        0.077571493840065148, -0.03224486958502952, -0.24229488706619015,
                        0.13842814590110342, 0.72430852843857441, 0.60382926979747287,
                0.16010239797412501 };
        lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        double hp1_a[] = {-0.16010239797412501, 0.60382926979747287, -0.72430852843857441,
                        0.13842814590110342, 0.24229488706619015, -0.03224486958502952,
                        -0.077571493840065148, -0.0062414902130117052, 0.012580751999015526,
                        0.0033357252850015492 };
         hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        double lp2_a[] = {0.16010239797412501, 0.60382926979747287, 0.72430852843857441,
                        0.13842814590110342, -0.24229488706619015, -0.03224486958502952,
                        0.077571493840065148, -0.0062414902130117052, -0.012580751999015526,
                        0.0033357252850015492 };
         lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

         double hp2_a[] = {0.0033357252850015492, 0.012580751999015526, -0.0062414902130117052,
                         -0.077571493840065148, -0.03224486958502952, 0.24229488706619015,
                         0.13842814590110342, -0.72430852843857441, 0.60382926979747287,
                         -0.16010239797412501 };
         hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
         return 0;
    }

    else if ( name == "db6"){
        double lp1_a[] = {-0.0010773010849955799,
                        0.0047772575110106514,
                        0.0005538422009938016,
                        -0.031582039318031156,
                        0.027522865530016288,
                        0.097501605587079362,
                        -0.12976686756709563,
                        -0.22626469396516913,
                        0.3152503517092432,
                        0.75113390802157753,
                        0.49462389039838539,
                        0.11154074335008017
};
        lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        double hp1_a[] = {-0.11154074335008017,
                        0.49462389039838539,
                        -0.75113390802157753,
                        0.3152503517092432,
                        0.22626469396516913,
                        -0.12976686756709563,
                        -0.097501605587079362,
                        0.027522865530016288,
                        0.031582039318031156,
                        0.0005538422009938016,
                        -0.0047772575110106514,
                        -0.0010773010849955799
};
         hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        double lp2_a[] = {0.11154074335008017,
                        0.49462389039838539,
                        0.75113390802157753,
                        0.3152503517092432,
                        -0.22626469396516913,
                        -0.12976686756709563,
                        0.097501605587079362,
                        0.027522865530016288,
                        -0.031582039318031156,
                        0.0005538422009938016,
                        0.0047772575110106514,
                        -0.0010773010849955799
};
         lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

         double hp2_a[] = {-0.0010773010849955799,
                         -0.0047772575110106514,
                         0.0005538422009938016,
                         0.031582039318031156,
                         0.027522865530016288,
                         -0.097501605587079362,
                         -0.12976686756709563,
                         0.22626469396516913,
                         0.3152503517092432,
                         -0.75113390802157753,
                         0.49462389039838539,
                         -0.11154074335008017
};
         hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
         return 0;
    }

    else if ( name == "db7"){
                double lp1_a[] = {0.00035371380000103988,
                                -0.0018016407039998328,
                                0.00042957797300470274,
                                0.012550998556013784,
                                -0.01657454163101562,
                                -0.038029936935034633,
                                0.080612609151065898,
                                0.071309219267050042,
                                -0.22403618499416572,
                                -0.14390600392910627,
                                0.4697822874053586,
                                0.72913209084655506,
                                0.39653931948230575,
                                0.077852054085062364
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.077852054085062364,
                                0.39653931948230575,
                                -0.72913209084655506,
                                0.4697822874053586,
                                0.14390600392910627,
                                -0.22403618499416572,
                                -0.071309219267050042,
                                0.080612609151065898,
                                0.038029936935034633,
                                -0.01657454163101562,
                                -0.012550998556013784,
                                0.0004295779730047027,
                                0.0018016407039998328,
                                0.00035371380000103988
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.077852054085062364,
                                0.39653931948230575,
                                0.72913209084655506,
                                0.4697822874053586,
                                -0.14390600392910627,
                                -0.22403618499416572,
                                0.071309219267050042,
                                0.080612609151065898,
                                -0.038029936935034633,
                                -0.01657454163101562,
                                0.012550998556013784,
                                0.00042957797300470274,
                                -0.0018016407039998328,
                                0.00035371380000103988
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.00035371380000103988,
                                 0.0018016407039998328,
                                 0.00042957797300470274,
                                 -0.01255099855601378,
                                 -0.01657454163101562,
                                 0.038029936935034633,
                                 0.080612609151065898,
                                 -0.071309219267050042,
                                 -0.22403618499416572,
                                 0.14390600392910627,
                                 0.4697822874053586,
                                 -0.72913209084655506,
                                 0.39653931948230575,
                                 -0.077852054085062364
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "db8"){
                double lp1_a[] = {-0.00011747678400228192,
                                0.00067544940599855677,
                                -0.00039174037299597711,
                                -0.0048703529930106603,
                                0.0087460940470156547,
                                0.013981027917015516,
                                -0.044088253931064719,
                                -0.017369301002022108,
                                0.12874742662018601,
                                0.00047248457399797254,
                                -0.28401554296242809,
                                -0.015829105256023893,
                                0.58535468365486909,
                                0.67563073629801285,
                                0.31287159091446592,
                                0.054415842243081609
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.054415842243081609,
                                0.31287159091446592,
                                -0.67563073629801285,
                                0.58535468365486909,
                                0.015829105256023893,
                                -0.28401554296242809,
                                -0.00047248457399797254,
                                0.12874742662018601,
                                0.017369301002022108,
                                -0.044088253931064719,
                                -0.013981027917015516,
                                0.0087460940470156547,
                                0.0048703529930106603,
                                -0.00039174037299597711,
                                -0.00067544940599855677,
                                -0.00011747678400228192
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.054415842243081609,
                                0.31287159091446592,
                                0.67563073629801285,
                                0.58535468365486909,
                                -0.015829105256023893,
                                -0.28401554296242809,
                                0.00047248457399797254,
                                0.12874742662018601,
                                -0.017369301002022108,
                                -0.044088253931064719,
                                0.013981027917015516,
                                0.0087460940470156547,
                                -0.0048703529930106603,
                                -0.00039174037299597711,
                                0.00067544940599855677,
                                -0.00011747678400228192
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.00011747678400228192,
                                 -0.00067544940599855677,
                                 -0.00039174037299597711,
                                 0.0048703529930106603,
                                 0.0087460940470156547,
                                 -0.013981027917015516,
                                 -0.044088253931064719,
                                 0.017369301002022108,
                                 0.12874742662018601,
                                 -0.00047248457399797254,
                                 -0.28401554296242809,
                                 0.015829105256023893,
                                 0.58535468365486909,
                                 -0.67563073629801285,
                                 0.31287159091446592,
                                 -0.054415842243081609
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "db9"){
                double lp1_a[] = {3.9347319995026124e-05,
                                -0.00025196318899817888,
                                0.00023038576399541288,
                                0.0018476468829611268,
                                -0.0042815036819047227,
                                -0.004723204757894831,
                                0.022361662123515244,
                                0.00025094711499193845,
                                -0.067632829059523988,
                                0.030725681478322865,
                                0.14854074933476008,
                                -0.096840783220879037,
                                -0.29327378327258685,
                                0.13319738582208895,
                                0.65728807803663891,
                                0.6048231236767786,
                                0.24383467463766728,
                                0.038077947363167282
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.038077947363167282,
                                0.24383467463766728,
                                -0.6048231236767786,
                                0.65728807803663891,
                                -0.13319738582208895,
                                -0.29327378327258685,
                                0.096840783220879037,
                                0.14854074933476008,
                                -0.030725681478322865,
                                -0.067632829059523988,
                                -0.00025094711499193845,
                                0.022361662123515244,
                                0.004723204757894831,
                                -0.0042815036819047227,
                                -0.0018476468829611268,
                                0.00023038576399541288,
                                0.00025196318899817888,
                                3.9347319995026124e-05
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.038077947363167282,
                                0.24383467463766728,
                                0.6048231236767786,
                                0.65728807803663891,
                                0.13319738582208895,
                                -0.29327378327258685,
                                -0.096840783220879037,
                                0.14854074933476008,
                                0.030725681478322865,
                                -0.067632829059523988,
                                0.00025094711499193845,
                                0.022361662123515244,
                                -0.004723204757894831,
                                -0.0042815036819047227,
                                0.0018476468829611268,
                                0.00023038576399541288,
                                -0.00025196318899817888,
                                3.9347319995026124e-05
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {3.9347319995026124e-05,
                                 0.00025196318899817888,
                                 0.00023038576399541288,
                                 -0.0018476468829611268,
                                 -0.0042815036819047227,
                                 0.004723204757894831,
                                 0.022361662123515244,
                                 -0.00025094711499193845,
                                 -0.067632829059523988,
                                 -0.030725681478322865,
                                 0.14854074933476008,
                                 0.096840783220879037,
                                 -0.29327378327258685,
                                 -0.13319738582208895,
                                 0.65728807803663891,
                                 -0.6048231236767786,
                                 0.24383467463766728,
                                 -0.038077947363167282
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "db10"){
                double lp1_a[] = {-1.3264203002354869e-05,
                                9.3588670001089845e-05,
                                -0.0001164668549943862,
                                -0.00068585669500468248,
                                0.0019924052949908499,
                                0.0013953517469940798,
                                -0.010733175482979604,
                                0.0036065535669883944,
                                0.033212674058933238,
                                -0.029457536821945671,
                                -0.071394147165860775,
                                0.093057364603806592,
                                0.12736934033574265,
                                -0.19594627437659665,
                                -0.24984642432648865,
                                0.28117234366042648,
                                0.68845903945259213,
                                0.52720118893091983,
                                0.18817680007762133,
                                0.026670057900950818
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.026670057900950818,
                                0.18817680007762133,
                                -0.52720118893091983,
                                0.68845903945259213,
                                -0.28117234366042648,
                                -0.24984642432648865,
                                0.19594627437659665,
                                0.12736934033574265,
                                -0.093057364603806592,
                                -0.071394147165860775,
                                0.029457536821945671,
                                0.033212674058933238,
                                -0.0036065535669883944,
                                -0.010733175482979604,
                                -0.0013953517469940798,
                                0.0019924052949908499,
                                0.00068585669500468248,
                                -0.0001164668549943862,
                                -9.3588670001089845e-05,
                                -1.3264203002354869e-05
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.026670057900950818,
                                0.18817680007762133,
                                0.52720118893091983,
                                0.68845903945259213,
                                0.28117234366042648,
                                -0.24984642432648865,
                                -0.19594627437659665,
                                0.12736934033574265,
                                0.093057364603806592,
                                -0.071394147165860775,
                                -0.029457536821945671,
                                0.033212674058933238,
                                0.0036065535669883944,
                                -0.010733175482979604,
                                0.0013953517469940798,
                                0.0019924052949908499,
                                -0.00068585669500468248,
                                -0.0001164668549943862,
                                9.3588670001089845e-05,
                                -1.3264203002354869e-05
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-1.3264203002354869e-05,
                                 -9.3588670001089845e-05,
                                 -0.0001164668549943862,
                                 0.00068585669500468248,
                                 0.0019924052949908499,
                                 -0.0013953517469940798,
                                 -0.010733175482979604,
                                 -0.0036065535669883944,
                                 0.033212674058933238,
                                 0.029457536821945671,
                                 -0.071394147165860775,
                                 -0.093057364603806592,
                                 0.12736934033574265,
                                 0.19594627437659665,
                                 -0.24984642432648865,
                                 -0.28117234366042648,
                                 0.68845903945259213,
                                 -0.52720118893091983,
                                 0.18817680007762133,
                                 -0.026670057900950818
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "db12"){
                double lp1_a[] = {-1.5290717580684923e-06,
                                1.2776952219379579e-05,
                                -2.4241545757030318e-05,
                                -8.8504109208203182e-05,
                                0.00038865306282092672,
                                6.5451282125215034e-06,
                                -0.0021795036186277044,
                                0.0022486072409952287,
                                0.0067114990087955486,
                                -0.012840825198299882,
                                -0.01221864906974642,
                                0.041546277495087637,
                                0.010849130255828966,
                                -0.09643212009649671,
                                0.0053595696743599965,
                                0.18247860592758275,
                                -0.023779257256064865,
                                -0.31617845375277914,
                                -0.044763885653777619,
                                0.51588647842780067,
                                0.65719872257929113,
                                0.37735513521420411,
                                0.10956627282118277,
                                0.013112257957229239
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.013112257957229239,
                                0.10956627282118277,
                                -0.37735513521420411,
                                0.65719872257929113,
                                -0.51588647842780067,
                                -0.044763885653777619,
                                0.31617845375277914,
                                -0.023779257256064865,
                                -0.18247860592758275,
                                0.0053595696743599965,
                                0.09643212009649671,
                                0.010849130255828966,
                                -0.041546277495087637,
                                -0.01221864906974642,
                                0.012840825198299882,
                                0.0067114990087955486,
                                -0.0022486072409952287,
                                -0.0021795036186277044,
                                -6.5451282125215034e-06,
                                0.00038865306282092672,
                                8.8504109208203182e-05,
                                -2.4241545757030318e-05,
                                -1.2776952219379579e-05,
                                -1.5290717580684923e-06
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.013112257957229239,
                                0.10956627282118277,
                                0.37735513521420411,
                                0.65719872257929113,
                                0.51588647842780067,
                                -0.044763885653777619,
                                -0.31617845375277914,
                                -0.023779257256064865,
                                0.18247860592758275,
                                0.0053595696743599965,
                                -0.09643212009649671,
                                0.010849130255828966,
                                0.041546277495087637,
                                -0.01221864906974642,
                                -0.012840825198299882,
                                0.0067114990087955486,
                                0.0022486072409952287,
                                -0.0021795036186277044,
                                6.5451282125215034e-06,
                                0.00038865306282092672,
                                -8.8504109208203182e-05,
                                -2.4241545757030318e-05,
                                1.2776952219379579e-05,
                                -1.5290717580684923e-06
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-1.5290717580684923e-06,
                                 -1.2776952219379579e-05,
                                 -2.4241545757030318e-05,
                                 8.8504109208203182e-05,
                                 0.00038865306282092672,
                                 -6.5451282125215034e-06,
                                 -0.0021795036186277044,
                                 -0.0022486072409952287,
                                 0.0067114990087955486,
                                 0.012840825198299882,
                                 -0.01221864906974642,
                                 -0.041546277495087637,
                                 0.010849130255828966,
                                 0.09643212009649671,
                                 0.0053595696743599965,
                                 -0.18247860592758275,
                                 -0.023779257256064865,
                                 0.31617845375277914,
                                 -0.044763885653777619,
                                 -0.51588647842780067,
                                 0.65719872257929113,
                                 -0.37735513521420411,
                                 0.10956627282118277,
                                 -0.013112257957229239
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "db13"){
                double lp1_a[] = {5.2200350984547998e-07,
                                -4.7004164793608082e-06,
                                1.0441930571407941e-05,
                                3.0678537579324358e-05,
                                -0.00016512898855650571,
                                4.9251525126285676e-05,
                                0.00093232613086724904,
                                -0.0013156739118922766,
                                -0.002761911234656831,
                                0.0072555894016171187,
                                0.0039239414487955773,
                                -0.023831420710327809,
                                0.0023799722540522269,
                                0.056139477100276156,
                                -0.026488406475345658,
                                -0.10580761818792761,
                                0.072948933656788742,
                                0.17947607942935084,
                                -0.12457673075080665,
                                -0.31497290771138414,
                                0.086985726179645007,
                                0.58888957043121193,
                                0.61105585115878114,
                                0.31199632216043488,
                                0.082861243872901946,
                                0.0092021335389622788
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.0092021335389622788,
                                0.082861243872901946,
                                -0.31199632216043488,
                                0.61105585115878114,
                                -0.58888957043121193,
                                0.086985726179645007,
                                0.31497290771138414,
                                -0.12457673075080665,
                                -0.17947607942935084,
                                0.072948933656788742,
                                0.10580761818792761,
                                -0.026488406475345658,
                                -0.056139477100276156,
                                0.0023799722540522269,
                                0.023831420710327809,
                                0.0039239414487955773,
                                -0.0072555894016171187,
                                -0.002761911234656831,
                                0.0013156739118922766,
                                0.00093232613086724904,
                                -4.9251525126285676e-05,
                                -0.00016512898855650571,
                                -3.0678537579324358e-05,
                                1.0441930571407941e-05,
                                4.7004164793608082e-06,
                                5.2200350984547998e-07
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0092021335389622788,
                                0.082861243872901946,
                                0.31199632216043488,
                                0.61105585115878114,
                                0.58888957043121193,
                                0.086985726179645007,
                                -0.31497290771138414,
                                -0.12457673075080665,
                                0.17947607942935084,
                                0.072948933656788742,
                                -0.10580761818792761,
                                -0.026488406475345658,
                                0.056139477100276156,
                                0.0023799722540522269,
                                -0.023831420710327809,
                                0.0039239414487955773,
                                0.0072555894016171187,
                                -0.002761911234656831,
                                -0.0013156739118922766,
                                0.00093232613086724904,
                                4.9251525126285676e-05,
                                -0.00016512898855650571,
                                3.0678537579324358e-05,
                                1.0441930571407941e-05,
                                -4.7004164793608082e-06,
                                5.2200350984547998e-07
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {5.2200350984547998e-07,
                                 4.7004164793608082e-06,
                                 1.0441930571407941e-05,
                                 -3.0678537579324358e-05,
                                 -0.00016512898855650571,
                                 -4.9251525126285676e-05,
                                 0.00093232613086724904,
                                 0.0013156739118922766,
                                 -0.002761911234656831,
                                 -0.0072555894016171187,
                                 0.0039239414487955773,
                                 0.023831420710327809,
                                 0.0023799722540522269,
                                 -0.056139477100276156,
                                 -0.026488406475345658,
                                 0.10580761818792761,
                                 0.072948933656788742,
                                 -0.17947607942935084,
                                 -0.12457673075080665,
                                 0.31497290771138414,
                                 0.086985726179645007,
                                 -0.58888957043121193,
                                 0.61105585115878114,
                                 -0.31199632216043488,
                                 0.082861243872901946,
                                 -0.0092021335389622788
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }

    else if ( name == "db11"){
                double lp1_a[] = {4.4942742772363519e-06,
                                -3.4634984186983789e-05,
                                5.4439074699366381e-05,
                                0.00024915252355281426,
                                -0.00089302325066623663,
                                -0.00030859285881515924,
                                0.0049284176560587777,
                                -0.0033408588730145018,
                                -0.015364820906201324,
                                0.020840904360180039,
                                0.031335090219045313,
                                -0.066438785695020222,
                                -0.04647995511667613,
                                0.14981201246638268,
                                0.066043588196690886,
                                -0.27423084681792875,
                                -0.16227524502747828,
                                0.41196436894789695,
                                0.68568677491617847,
                                0.44989976435603013,
                                0.14406702115061959,
                                0.018694297761470441
  };
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.018694297761470441,
                                0.14406702115061959,
                                -0.44989976435603013,
                                0.68568677491617847,
                                -0.41196436894789695,
                                -0.16227524502747828,
                                0.27423084681792875,
                                0.066043588196690886,
                                -0.14981201246638268,
                                -0.04647995511667613,
                                0.066438785695020222,
                                0.031335090219045313,
                                -0.020840904360180039,
                                -0.015364820906201324,
                                0.0033408588730145018,
                                0.0049284176560587777,
                                0.00030859285881515924,
                                -0.00089302325066623663,
                                -0.00024915252355281426,
                                5.4439074699366381e-05,
                                3.4634984186983789e-05,
                                4.4942742772363519e-06
  };
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.018694297761470441,
                                0.14406702115061959,
                                0.44989976435603013,
                                0.68568677491617847,
                                0.41196436894789695,
                                -0.16227524502747828,
                                -0.27423084681792875,
                                0.066043588196690886,
                                0.14981201246638268,
                                -0.04647995511667613,
                                -0.066438785695020222,
                                0.031335090219045313,
                                0.020840904360180039,
                                -0.015364820906201324,
                                -0.0033408588730145018,
                                0.0049284176560587777,
                                -0.00030859285881515924,
                                -0.00089302325066623663,
                                0.00024915252355281426,
                                5.4439074699366381e-05,
                                -3.4634984186983789e-05,
                                4.4942742772363519e-06
  };
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {4.4942742772363519e-06,
                                 3.4634984186983789e-05,
                                 5.4439074699366381e-05,
                                 -0.00024915252355281426,
                                 -0.00089302325066623663,
                                 0.00030859285881515924,
                                 0.0049284176560587777,
                                 0.0033408588730145018,
                                 -0.015364820906201324,
                                 -0.020840904360180039,
                                 0.031335090219045313,
                                 0.066438785695020222,
                                 -0.04647995511667613,
                                 -0.14981201246638268,
                                 0.066043588196690886,
                                 0.27423084681792875,
                                 -0.16227524502747828,
                                 -0.41196436894789695,
                                 0.68568677491617847,
                                 -0.44989976435603013,
                                 0.14406702115061959,
                                 -0.018694297761470441
  };
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
              }

    else if ( name == "db14"){
                        double lp1_a[] = {-1.7871399683109222e-07,
                                        1.7249946753674012e-06,
                                        -4.3897049017804176e-06,
                                        -1.0337209184568496e-05,
                                        6.875504252695734e-05,
                                        -4.1777245770370672e-05,
                                        -0.00038683194731287514,
                                        0.00070802115423540481,
                                        0.001061691085606874,
                                        -0.003849638868019787,
                                        -0.00074621898926387534,
                                        0.012789493266340071,
                                        -0.0056150495303375755,
                                        -0.030185351540353976,
                                        0.026981408307947971,
                                        0.05523712625925082,
                                        -0.071548955503983505,
                                        -0.086748411568110598,
                                        0.13998901658445695,
                                        0.13839521386479153,
                                        -0.21803352999321651,
                                        -0.27168855227867705,
                                        0.21867068775886594,
                                        0.63118784910471981,
                                        0.55430561794077093,
                                        0.25485026779256437,
                                        0.062364758849384874,
                                        0.0064611534600864905
};
                        lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                        double hp1_a[] = {-0.0064611534600864905,
                                        0.062364758849384874,
                                        -0.25485026779256437,
                                        0.55430561794077093,
                                        -0.63118784910471981,
                                        0.21867068775886594,
                                        0.27168855227867705,
                                        -0.21803352999321651,
                                        -0.13839521386479153,
                                        0.13998901658445695,
                                        0.086748411568110598,
                                        -0.071548955503983505,
                                        -0.05523712625925082,
                                        0.026981408307947971,
                                        0.030185351540353976,
                                        -0.0056150495303375755,
                                        -0.012789493266340071,
                                        -0.00074621898926387534,
                                        0.003849638868019787,
                                        0.001061691085606874,
                                        -0.00070802115423540481,
                                        -0.00038683194731287514,
                                        4.1777245770370672e-05,
                                        6.875504252695734e-05,
                                        1.0337209184568496e-05,
                                        -4.3897049017804176e-06,
                                        -1.7249946753674012e-06,
                                        -1.7871399683109222e-07
};
                         hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                        double lp2_a[] = {0.0064611534600864905,
                                        0.062364758849384874,
                                        0.25485026779256437,
                                        0.55430561794077093,
                                        0.63118784910471981,
                                        0.21867068775886594,
                                        -0.27168855227867705,
                                        -0.21803352999321651,
                                        0.13839521386479153,
                                        0.13998901658445695,
                                        -0.086748411568110598,
                                        -0.071548955503983505,
                                        0.05523712625925082,
                                        0.026981408307947971,
                                        -0.030185351540353976,
                                        -0.0056150495303375755,
                                        0.012789493266340071,
                                        -0.00074621898926387534,
                                        -0.003849638868019787,
                                        0.001061691085606874,
                                        0.00070802115423540481,
                                        -0.00038683194731287514,
                                        -4.1777245770370672e-05,
                                        6.875504252695734e-05,
                                        -1.0337209184568496e-05,
                                        -4.3897049017804176e-06,
                                        1.7249946753674012e-06,
                                        -1.7871399683109222e-07
};
                         lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                         double hp2_a[] = {-1.7871399683109222e-07,
                                        -1.7249946753674012e-06,
                                        -4.3897049017804176e-06,
                                        1.0337209184568496e-05,
                                        6.875504252695734e-05,
                                        4.1777245770370672e-05,
                                        -0.00038683194731287514,
                                        -0.00070802115423540481,
                                        0.001061691085606874,
                                        0.003849638868019787,
                                        -0.00074621898926387534,
                                        -0.012789493266340071,
                                        -0.0056150495303375755,
                                        0.030185351540353976,
                                        0.026981408307947971,
                                        -0.05523712625925082,
                                        -0.071548955503983505,
                                        0.086748411568110598,
                                        0.13998901658445695,
                                        -0.13839521386479153,
                                        -0.21803352999321651,
                                        0.27168855227867705,
                                        0.21867068775886594,
                                        -0.63118784910471981,
                                        0.55430561794077093,
                                        -0.25485026779256437,
                                        0.062364758849384874,
                                        -0.0064611534600864905
};
                         hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                         return 0;
                  }
    else if ( name == "db15"){
                double lp1_a[] = {6.1333599133037138e-08,
                                -6.3168823258794506e-07,
                                1.8112704079399406e-06,
                                3.3629871817363823e-06,
                                -2.8133296266037558e-05,
                                2.579269915531323e-05,
                                0.00015589648992055726,
                                -0.00035956524436229364,
                                -0.00037348235413726472,
                                0.0019433239803823459,
                                -0.00024175649075894543,
                                -0.0064877345603061454,
                                0.0051010003604228726,
                                0.015083918027862582,
                                -0.020810050169636805,
                                -0.025767007328366939,
                                0.054780550584559995,
                                0.033877143923563204,
                                -0.11112093603713753,
                                -0.039666176555733602,
                                0.19014671400708816,
                                0.065282952848765688,
                                -0.28888259656686216,
                                -0.19320413960907623,
                                0.33900253545462167,
                                0.64581314035721027,
                                0.49263177170797529,
                                0.20602386398692688,
                                0.046743394892750617,
                                0.0045385373615773762
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.0045385373615773762,
                                0.046743394892750617,
                                -0.20602386398692688,
                                0.49263177170797529,
                                -0.64581314035721027,
                                0.33900253545462167,
                                0.19320413960907623,
                                -0.28888259656686216,
                                -0.065282952848765688,
                                0.19014671400708816,
                                0.039666176555733602,
                                -0.11112093603713753,
                                -0.033877143923563204,
                                0.054780550584559995,
                                0.025767007328366939,
                                -0.020810050169636805,
                                -0.015083918027862582,
                                0.0051010003604228726,
                                0.0064877345603061454,
                                -0.00024175649075894543,
                                -0.0019433239803823459,
                                -0.00037348235413726472,
                                0.00035956524436229364,
                                0.00015589648992055726,
                                -2.579269915531323e-05,
                                -2.8133296266037558e-05,
                                -3.3629871817363823e-06,
                                1.8112704079399406e-06,
                                6.3168823258794506e-07,
                                6.1333599133037138e-08
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0045385373615773762,
                                0.046743394892750617,
                                0.20602386398692688,
                                0.49263177170797529,
                                0.64581314035721027,
                                0.33900253545462167,
                                -0.19320413960907623,
                                -0.28888259656686216,
                                0.065282952848765688,
                                0.19014671400708816,
                                -0.039666176555733602,
                                -0.11112093603713753,
                                0.033877143923563204,
                                0.054780550584559995,
                                -0.025767007328366939,
                                -0.020810050169636805,
                                0.015083918027862582,
                                0.0051010003604228726,
                                -0.0064877345603061454,
                                -0.00024175649075894543,
                                0.0019433239803823459,
                                -0.00037348235413726472,
                                -0.00035956524436229364,
                                0.00015589648992055726,
                                2.579269915531323e-05,
                                -2.8133296266037558e-05,
                                3.3629871817363823e-06,
                                1.8112704079399406e-06,
                                -6.3168823258794506e-07,
                                6.1333599133037138e-08
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {6.1333599133037138e-08,
                                 6.3168823258794506e-07,
                                 1.8112704079399406e-06,
                                 -3.3629871817363823e-06,
                                 -2.8133296266037558e-05,
                                 -2.579269915531323e-05,
                                 0.00015589648992055726,
                                 0.00035956524436229364,
                                 -0.00037348235413726472,
                                 -0.0019433239803823459,
                                 -0.00024175649075894543,
                                 0.0064877345603061454,
                                 0.0051010003604228726,
                                 -0.015083918027862582,
                                 -0.020810050169636805,
                                 0.025767007328366939,
                                 0.054780550584559995,
                                 -0.033877143923563204,
                                 -0.11112093603713753,
                                 0.039666176555733602,
                                 0.19014671400708816,
                                 -0.065282952848765688,
                                 -0.28888259656686216,
                                 0.19320413960907623,
                                 0.33900253545462167,
                                 -0.64581314035721027,
                                 0.49263177170797529,
                                 -0.20602386398692688,
                                 0.046743394892750617,
                                 -0.0045385373615773762
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "bior1.1"){
                double lp1_a[] = {0.70710678118654757,
                                0.70710678118654757
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.70710678118654757,
                                0.70710678118654757
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.70710678118654757,
                                0.70710678118654757
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.70710678118654757,
                                 -0.70710678118654757
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "bior1.3"){
                double lp1_a[] = {-0.088388347648318447,
                                0.088388347648318447,
                                0.70710678118654757,
                                0.70710678118654757,
                                0.088388347648318447,
                                -0.088388347648318447,
    };
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                -0.70710678118654757,
                                0.70710678118654757,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                 0.0,
                                 0.70710678118654757,
                                 0.70710678118654757,
                                 0.0,
                                 0.0
    };
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.088388347648318447,
                                 -0.088388347648318447,
                                 0.70710678118654757,
                                 -0.70710678118654757,
                                 0.088388347648318447,
                                 0.088388347648318447
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }

    else if ( name == "bior1.5"){
                double lp1_a[] = {0.01657281518405971,
                                -0.01657281518405971,
                                -0.12153397801643787,
                                0.12153397801643787,
                                0.70710678118654757,
                                0.70710678118654757,
                                0.12153397801643787,
                                -0.12153397801643787,
                                -0.01657281518405971,
                                0.01657281518405971
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                -0.70710678118654757,
                                0.70710678118654757,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.70710678118654757,
                                0.70710678118654757,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.01657281518405971,
                                 0.01657281518405971,
                                 -0.12153397801643787,
                                 -0.12153397801643787,
                                 0.70710678118654757,
                                 -0.70710678118654757,
                                 0.12153397801643787,
                                 0.12153397801643787,
                                 -0.01657281518405971,
                                 -0.01657281518405971
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "bior2.2"){
                double lp1_a[] = {0.0,
                                -0.17677669529663689,
                                0.35355339059327379,
                                1.0606601717798214,
                                0.35355339059327379,
                                -0.17677669529663689
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.35355339059327379,
                                -0.70710678118654757,
                                0.35355339059327379,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.35355339059327379,
                                0.70710678118654757,
                                0.35355339059327379,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0,
                                 0.17677669529663689,
                                 0.35355339059327379,
                                 -1.0606601717798214,
                                 0.35355339059327379,
                                 0.17677669529663689

};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "bior2.4"){
                double lp1_a[] = {0.0,
                                0.033145630368119419,
                                -0.066291260736238838,
                                -0.17677669529663689,
                                0.4198446513295126,
                                0.99436891104358249,
                                0.4198446513295126,
                                -0.17677669529663689,
                                -0.066291260736238838,
                                0.033145630368119419
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                0.0,
                                0.35355339059327379,
                                -0.70710678118654757,
                                0.35355339059327379,
                                0.0,
                                0.0,
                                0.0,
                                0.0

};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.0,
                                0.35355339059327379,
                                0.70710678118654757,
                                0.35355339059327379,
                                0.0,
                                0.0,
                                0.0,
                                0.0

};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0,
                                 -0.033145630368119419,
                                 -0.066291260736238838,
                                 0.17677669529663689,
                                 0.4198446513295126,
                                 -0.99436891104358249,
                                 0.4198446513295126,
                                 0.17677669529663689,
                                 -0.066291260736238838,
                                 -0.033145630368119419
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "bior2.6"){
                double lp1_a[] = {0.0,
                                -0.0069053396600248784,
                                0.013810679320049757,
                                0.046956309688169176,
                                -0.10772329869638811,
                                -0.16987135563661201,
                                0.44746600996961211,
                                0.96674755240348298,
                                0.44746600996961211,
                                -0.16987135563661201,
                                -0.10772329869638811,
                                0.046956309688169176,
                                0.013810679320049757,
                                -0.0069053396600248784
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.35355339059327379,
                                -0.70710678118654757,
                                0.35355339059327379,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.35355339059327379,
                                0.70710678118654757,
                                0.35355339059327379,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0,
                                 0.0069053396600248784,
                                 0.013810679320049757,
                                 -0.046956309688169176,
                                 -0.10772329869638811,
                                 0.16987135563661201,
                                 0.44746600996961211,
                                 -0.96674755240348298,
                                 0.44746600996961211,
                                 0.16987135563661201,
                                 -0.10772329869638811,
                                 -0.046956309688169176,
                                 0.013810679320049757,
                                 0.0069053396600248784
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "bior2.8"){
                double lp1_a[] = {0.0,
                                0.0015105430506304422,
                                -0.0030210861012608843,
                                -0.012947511862546647,
                                0.028916109826354178,
                                0.052998481890690945,
                                -0.13491307360773608,
                                -0.16382918343409025,
                                0.46257144047591658,
                                0.95164212189717856,
                                0.46257144047591658,
                                -0.16382918343409025,
                                -0.13491307360773608,
                                0.052998481890690945,
                                0.028916109826354178,
                                -0.012947511862546647,
                                -0.0030210861012608843,
                                0.0015105430506304422
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.35355339059327379,
                                -0.70710678118654757,
                                0.35355339059327379,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.35355339059327379,
                                0.70710678118654757,
                                0.35355339059327379,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0,
                                 -0.0015105430506304422,
                                 -0.0030210861012608843,
                                 0.012947511862546647,
                                 0.028916109826354178,
                                 -0.052998481890690945,
                                 -0.13491307360773608,
                                 0.16382918343409025,
                                 0.46257144047591658,
                                 -0.95164212189717856,
                                 0.46257144047591658,
                                 0.16382918343409025,
                                 -0.13491307360773608,
                                 -0.052998481890690945,
                                 0.028916109826354178,
                                 0.012947511862546647,
                                 -0.0030210861012608843,
                                 -0.0015105430506304422
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }

    else if ( name == "bior3.1"){
                double lp1_a[] = {-0.35355339059327379,
                                1.0606601717798214,
                                1.0606601717798214,
                                -0.35355339059327379
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.17677669529663689,
                                0.53033008588991071,
                                -0.53033008588991071,
                                0.17677669529663689
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.17677669529663689,
                                0.53033008588991071,
                                0.53033008588991071,
                                0.17677669529663689
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.35355339059327379,
                                 -1.0606601717798214,
                                 1.0606601717798214,
                                 0.35355339059327379
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "bior3.3"){
                double lp1_a[] = {0.066291260736238838,
                                -0.19887378220871652,
                                -0.15467960838455727,
                                0.99436891104358249,
                                0.99436891104358249,
                                -0.15467960838455727,
                                -0.19887378220871652,
                                0.066291260736238838
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                -0.17677669529663689,
                                0.53033008588991071,
                                -0.53033008588991071,
                                0.17677669529663689,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.17677669529663689,
                                0.53033008588991071,
                                0.53033008588991071,
                                0.17677669529663689,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.066291260736238838,
                                 0.19887378220871652,
                                 -0.15467960838455727,
                                 -0.99436891104358249,
                                 0.99436891104358249,
                                 0.15467960838455727,
                                 -0.19887378220871652,
                                 -0.066291260736238838
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "bior3.5"){
                double lp1_a[] = {-0.013810679320049757,
                                0.041432037960149271,
                                0.052480581416189075,
                                -0.26792717880896527,
                                -0.071815532464258744,
                                0.96674755240348298,
                                0.96674755240348298,
                                -0.071815532464258744,
                                -0.26792717880896527,
                                0.052480581416189075,
                                0.041432037960149271,
                                -0.013810679320049757
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                -0.17677669529663689,
                                0.53033008588991071,
                                -0.53033008588991071,
                                0.17677669529663689,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.17677669529663689,
                                0.53033008588991071,
                                0.53033008588991071,
                                0.17677669529663689,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.013810679320049757,
                                 -0.041432037960149271,
                                 0.052480581416189075,
                                 0.26792717880896527,
                                 -0.071815532464258744,
                                 -0.96674755240348298,
                                 0.96674755240348298,
                                 0.071815532464258744,
                                 -0.26792717880896527,
                                 -0.052480581416189075,
                                 0.041432037960149271,
                                 0.013810679320049757
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }

    else if ( name == "bior3.7"){
                double lp1_a[] = {0.0030210861012608843,
                                -0.0090632583037826529,
                                -0.016831765421310641,
                                0.074663985074019001,
                                0.031332978707362888,
                                -0.301159125922835,
                                -0.026499240945345472,
                                0.95164212189717856,
                                0.95164212189717856,
                                -0.026499240945345472,
                                -0.301159125922835,
                                0.031332978707362888,
                                0.074663985074019001,
                                -0.016831765421310641,
                                -0.0090632583037826529,
                                0.0030210861012608843
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                -0.17677669529663689,
                                0.53033008588991071,
                                -0.53033008588991071,
                                0.17677669529663689,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.17677669529663689,
                                0.53033008588991071,
                                0.53033008588991071,
                                0.17677669529663689,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0030210861012608843,
                                 0.0090632583037826529,
                                 -0.016831765421310641,
                                 -0.074663985074019001,
                                 0.031332978707362888,
                                 0.301159125922835,
                                 -0.026499240945345472,
                                 -0.95164212189717856,
                                 0.95164212189717856,
                                 0.026499240945345472,
                                 -0.301159125922835,
                                 -0.031332978707362888,
                                 0.074663985074019001,
                                 0.016831765421310641,
                                 -0.0090632583037826529,
                                 -0.0030210861012608843
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "bior3.9"){
                double lp1_a[] = {-0.00067974437278369901,
                                0.0020392331183510968,
                                0.0050603192196119811,
                                -0.020618912641105536,
                                -0.014112787930175846,
                                0.09913478249423216,
                                0.012300136269419315,
                                -0.32019196836077857,
                                0.0020500227115698858,
                                0.94212570067820678,
                                0.94212570067820678,
                                0.0020500227115698858,
                                -0.32019196836077857,
                                0.012300136269419315,
                                0.09913478249423216,
                                -0.014112787930175846,
                                -0.020618912641105536,
                                0.0050603192196119811,
                                0.0020392331183510968,
                                -0.00067974437278369901
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                -0.17677669529663689,
                                0.53033008588991071,
                                -0.53033008588991071,
                                0.17677669529663689,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.17677669529663689,
                                0.53033008588991071,
                                0.53033008588991071,
                                0.17677669529663689,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.00067974437278369901,
                                 -0.0020392331183510968,
                                 0.0050603192196119811,
                                 0.020618912641105536,
                                 -0.014112787930175846,
                                 -0.09913478249423216,
                                 0.012300136269419315,
                                 0.32019196836077857,
                                 0.0020500227115698858,
                                 -0.94212570067820678,
                                 0.94212570067820678,
                                 -0.0020500227115698858,
                                 -0.32019196836077857,
                                 -0.012300136269419315,
                                 0.09913478249423216,
                                 0.014112787930175846,
                                 -0.020618912641105536,
                                 -0.0050603192196119811,
                                 0.0020392331183510968,
                                 0.00067974437278369901
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "bior4.4"){
                double lp1_a[] = {0.0,
                                0.03782845550726404,
                                -0.023849465019556843,
                                -0.11062440441843718,
                                0.37740285561283066,
                                0.85269867900889385,
                                0.37740285561283066,
                                -0.11062440441843718,
                                -0.023849465019556843,
                                0.03782845550726404
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                -0.064538882628697058,
                                0.040689417609164058,
                                0.41809227322161724,
                                -0.7884856164055829,
                                0.41809227322161724,
                                0.040689417609164058,
                                -0.064538882628697058,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                -0.064538882628697058,
                                -0.040689417609164058,
                                0.41809227322161724,
                                0.7884856164055829,
                                0.41809227322161724,
                                -0.040689417609164058,
                                -0.064538882628697058,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0,
                                 -0.03782845550726404,
                                 -0.023849465019556843,
                                 0.11062440441843718,
                                 0.37740285561283066,
                                 -0.85269867900889385,
                                 0.37740285561283066,
                                 0.11062440441843718,
                                 -0.023849465019556843,
                                 -0.03782845550726404
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "bior5.5"){
                double lp1_a[] = {0.0,
                                0.0,
                                0.03968708834740544,
                                0.0079481086372403219,
                                -0.054463788468236907,
                                0.34560528195603346,
                                0.73666018142821055,
                                0.34560528195603346,
                                -0.054463788468236907,
                                0.0079481086372403219,
                                0.03968708834740544,
                                0.0
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.013456709459118716,
                                -0.0026949668801115071,
                                0.13670658466432914,
                                -0.093504697400938863,
                                -0.47680326579848425,
                                0.89950610974864842,
                                -0.47680326579848425,
                                -0.093504697400938863,
                                0.13670658466432914,
                                -0.0026949668801115071,
                                -0.013456709459118716,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.013456709459118716,
                                -0.0026949668801115071,
                                -0.13670658466432914,
                                -0.093504697400938863,
                                0.47680326579848425,
                                0.89950610974864842,
                                0.47680326579848425,
                                -0.093504697400938863,
                                -0.13670658466432914,
                                -0.0026949668801115071,
                                0.013456709459118716,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0,
                                 0.0,
                                 0.03968708834740544,
                                 -0.0079481086372403219,
                                 -0.054463788468236907,
                                 -0.34560528195603346,
                                 0.73666018142821055,
                                 -0.34560528195603346,
                                 -0.054463788468236907,
                                 -0.0079481086372403219,
                                 0.03968708834740544,
                                 0.0
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "bior6.8"){
                double lp1_a[] = {0.0,
                                0.0019088317364812906,
                                -0.0019142861290887667,
                                -0.016990639867602342,
                                0.01193456527972926,
                                0.04973290349094079,
                                -0.077263173167204144,
                                -0.09405920349573646,
                                0.42079628460982682,
                                0.82592299745840225,
                                0.42079628460982682,
                                -0.09405920349573646,
                                -0.077263173167204144,
                                0.04973290349094079,
                                0.01193456527972926,
                                -0.016990639867602342,
                                -0.0019142861290887667,
                                0.0019088317364812906
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0,
                                0.0,
                                0.0,
                                0.014426282505624435,
                                -0.014467504896790148,
                                -0.078722001062628819,
                                0.040367979030339923,
                                0.41784910915027457,
                                -0.75890772945365415,
                                0.41784910915027457,
                                0.040367979030339923,
                                -0.078722001062628819,
                                -0.014467504896790148,
                                0.014426282505624435,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0,
                                0.0,
                                0.0,
                                0.014426282505624435,
                                0.014467504896790148,
                                -0.078722001062628819,
                                -0.040367979030339923,
                                0.41784910915027457,
                                0.75890772945365415,
                                0.41784910915027457,
                                -0.040367979030339923,
                                -0.078722001062628819,
                                0.014467504896790148,
                                0.014426282505624435,
                                0.0,
                                0.0,
                                0.0,
                                0.0
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0,
                                 -0.0019088317364812906,
                                 -0.0019142861290887667,
                                 0.016990639867602342,
                                 0.01193456527972926,
                                 -0.04973290349094079,
                                 -0.077263173167204144,
                                 0.09405920349573646,
                                 0.42079628460982682,
                                 -0.82592299745840225,
                                 0.42079628460982682,
                                 0.09405920349573646,
                                 -0.077263173167204144,
                                 -0.04973290349094079,
                                 0.01193456527972926,
                                 0.016990639867602342,
                                 -0.0019142861290887667,
                                 -0.0019088317364812906
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }

    else if ( name == "coif1"){
                double lp1_a[] = {-0.01565572813546454,
                                -0.072732619512853897,
                                0.38486484686420286,
                                0.85257202021225542,
                                0.33789766245780922,
                                -0.072732619512853897
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.072732619512853897,
                                0.33789766245780922,
                                -0.85257202021225542,
                                0.38486484686420286,
                                0.072732619512853897,
                                -0.01565572813546454
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {-0.072732619512853897,
                                0.33789766245780922,
                                0.85257202021225542,
                                0.38486484686420286,
                                -0.072732619512853897,
                                -0.01565572813546454
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.01565572813546454,
                                 0.072732619512853897,
                                 0.38486484686420286,
                                 -0.85257202021225542,
                                 0.33789766245780922,
                                 0.072732619512853897
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "coif2"){
                double lp1_a[] = {-0.00072054944536451221,
                                -0.0018232088707029932,
                                0.0056114348193944995,
                                0.023680171946334084,
                                -0.059434418646456898,
                                -0.076488599078306393,
                                0.41700518442169254,
                                0.81272363544554227,
                                0.38611006682116222,
                                -0.067372554721963018,
                                -0.041464936781759151,
                                0.016387336463522112
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.016387336463522112,
                                -0.041464936781759151,
                                0.067372554721963018,
                                0.38611006682116222,
                                -0.81272363544554227,
                                0.41700518442169254,
                                0.076488599078306393,
                                -0.059434418646456898,
                                -0.023680171946334084,
                                0.0056114348193944995,
                                0.0018232088707029932,
                                -0.00072054944536451221
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.016387336463522112,
                                -0.041464936781759151,
                                -0.067372554721963018,
                                0.38611006682116222,
                                0.81272363544554227,
                                0.41700518442169254,
                                -0.076488599078306393,
                                -0.059434418646456898,
                                0.023680171946334084,
                                0.0056114348193944995,
                                -0.0018232088707029932,
                                -0.00072054944536451221
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.00072054944536451221,
                                 0.0018232088707029932,
                                 0.0056114348193944995,
                                 -0.023680171946334084,
                                 -0.059434418646456898,
                                 0.076488599078306393,
                                 0.41700518442169254,
                                 -0.81272363544554227,
                                 0.38611006682116222,
                                 0.067372554721963018,
                                 -0.041464936781759151,
                                 -0.016387336463522112
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "coif3"){
                double lp1_a[] = {-3.4599772836212559e-05,
                                -7.0983303138141252e-05,
                                0.00046621696011288631,
                                0.0011175187708906016,
                                -0.0025745176887502236,
                                -0.0090079761366615805,
                                0.015880544863615904,
                                0.034555027573061628,
                                -0.082301927106885983,
                                -0.071799821619312018,
                                0.42848347637761874,
                                0.79377722262562056,
                                0.4051769024096169,
                                -0.061123390002672869,
                                -0.0657719112818555,
                                0.023452696141836267,
                                0.0077825964273254182,
                                -0.0037935128644910141
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.0037935128644910141,
                                0.0077825964273254182,
                                -0.023452696141836267,
                                -0.0657719112818555,
                                0.061123390002672869,
                                0.4051769024096169,
                                -0.79377722262562056,
                                0.42848347637761874,
                                0.071799821619312018,
                                -0.082301927106885983,
                                -0.034555027573061628,
                                0.015880544863615904,
                                0.0090079761366615805,
                                -0.0025745176887502236,
                                -0.0011175187708906016,
                                0.00046621696011288631,
                                7.0983303138141252e-05,
                                -3.4599772836212559e-05
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {-0.0037935128644910141,
                                0.0077825964273254182,
                                0.023452696141836267,
                                -0.0657719112818555,
                                -0.061123390002672869,
                                0.4051769024096169,
                                0.79377722262562056,
                                0.42848347637761874,
                                -0.071799821619312018,
                                -0.082301927106885983,
                                0.034555027573061628,
                                0.015880544863615904,
                                -0.0090079761366615805,
                                -0.0025745176887502236,
                                0.0011175187708906016,
                                0.00046621696011288631,
                                -7.0983303138141252e-05,
                                -3.4599772836212559e-05
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-3.4599772836212559e-05,
                                 7.0983303138141252e-05,
                                 0.00046621696011288631,
                                 -0.0011175187708906016,
                                 -0.0025745176887502236,
                                 0.0090079761366615805,
                                 0.015880544863615904,
                                 -0.034555027573061628,
                                 -0.082301927106885983,
                                 0.071799821619312018,
                                 0.42848347637761874,
                                 -0.79377722262562056,
                                 0.4051769024096169,
                                 0.061123390002672869,
                                 -0.0657719112818555,
                                 -0.023452696141836267,
                                 0.0077825964273254182,
                                 0.0037935128644910141
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "coif4"){
                double lp1_a[] = {-1.7849850030882614e-06,
                                -3.2596802368833675e-06,
                                3.1229875865345646e-05,
                                6.2339034461007128e-05,
                                -0.00025997455248771324,
                                -0.00058902075624433831,
                                0.0012665619292989445,
                                0.0037514361572784571,
                                -0.0056582866866107199,
                                -0.015211731527946259,
                                0.025082261844864097,
                                0.039334427123337491,
                                -0.096220442033987982,
                                -0.066627474263425038,
                                0.4343860564914685,
                                0.78223893092049901,
                                0.41530840703043026,
                                -0.056077313316754807,
                                -0.081266699680878754,
                                0.026682300156053072,
                                0.016068943964776348,
                                -0.0073461663276420935,
                                -0.0016294920126017326,
                                0.00089231366858231456
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.00089231366858231456,
                                -0.0016294920126017326,
                                0.0073461663276420935,
                                0.016068943964776348,
                                -0.026682300156053072,
                                -0.081266699680878754,
                                0.056077313316754807,
                                0.41530840703043026,
                                -0.78223893092049901,
                                0.4343860564914685,
                                0.066627474263425038,
                                -0.096220442033987982,
                                -0.039334427123337491,
                                0.025082261844864097,
                                0.015211731527946259,
                                -0.0056582866866107199,
                                -0.0037514361572784571,
                                0.0012665619292989445,
                                0.00058902075624433831,
                                -0.00025997455248771324,
                                -6.2339034461007128e-05,
                                3.1229875865345646e-05,
                                3.2596802368833675e-06,
                                -1.7849850030882614e-06
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.00089231366858231456,
                                -0.0016294920126017326,
                                -0.0073461663276420935,
                                0.016068943964776348,
                                0.026682300156053072,
                                -0.081266699680878754,
                                -0.056077313316754807,
                                0.41530840703043026,
                                0.78223893092049901,
                                0.4343860564914685,
                                -0.066627474263425038,
                                -0.096220442033987982,
                                0.039334427123337491,
                                0.025082261844864097,
                                -0.015211731527946259,
                                -0.0056582866866107199,
                                0.0037514361572784571,
                                0.0012665619292989445,
                                -0.00058902075624433831,
                                -0.00025997455248771324,
                                6.2339034461007128e-05,
                                3.1229875865345646e-05,
                                -3.2596802368833675e-06,
                                -1.7849850030882614e-06
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-1.7849850030882614e-06,
                                 3.2596802368833675e-06,
                                 3.1229875865345646e-05,
                                 -6.2339034461007128e-05,
                                 -0.00025997455248771324,
                                 0.00058902075624433831,
                                 0.0012665619292989445,

                                 -0.0037514361572784571,
                                 -0.0056582866866107199,
                                 0.015211731527946259,
                                 0.025082261844864097,
                                 -0.039334427123337491,
                                 -0.096220442033987982,
                                 0.066627474263425038,
                                 0.4343860564914685,
                                 -0.78223893092049901,
                                 0.41530840703043026,
                                 0.056077313316754807,
                                 -0.081266699680878754,
                                 -0.026682300156053072,
                                 0.016068943964776348,
                                 0.0073461663276420935,
                                 -0.0016294920126017326,
                                 -0.00089231366858231456
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }
    else if ( name == "coif5"){
                double lp1_a[] = {-9.517657273819165e-08,
                                -1.6744288576823017e-07,
                                2.0637618513646814e-06,
                                3.7346551751414047e-06,
                                -2.1315026809955787e-05,
                                -4.1340432272512511e-05,
                                0.00014054114970203437,
                                0.00030225958181306315,
                                -0.00063813134304511142,
                                -0.0016628637020130838,
                                0.0024333732126576722,
                                0.0067641854480530832,
                                -0.0091642311624818458,
                                -0.019761778942572639,
                                0.032683574267111833,
                                0.041289208750181702,
                                -0.10557420870333893,
                                -0.062035963962903569,
                                0.43799162617183712,
                                0.77428960365295618,
                                0.42156620669085149,
                                -0.052043163176243773,
                                -0.091920010559696244,
                                0.02816802897093635,
                                0.023408156785839195,
                                -0.010131117519849788,
                                -0.004159358781386048,
                                0.0021782363581090178,
                                0.00035858968789573785,
                                -0.00021208083980379827
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.00021208083980379827,
                                0.00035858968789573785,
                                -0.0021782363581090178,
                                -0.004159358781386048,
                                0.010131117519849788,
                                0.023408156785839195,
                                -0.02816802897093635,
                                -0.091920010559696244,
                                0.052043163176243773,
                                0.42156620669085149,
                                -0.77428960365295618,
                                0.43799162617183712,
                                0.062035963962903569,
                                -0.10557420870333893,
                                -0.041289208750181702,
                                0.032683574267111833,
                                0.019761778942572639,
                                -0.0091642311624818458,
                                -0.0067641854480530832,
                                0.0024333732126576722,
                                0.0016628637020130838,
                                -0.00063813134304511142,
                                -0.00030225958181306315,
                                0.00014054114970203437,
                                4.1340432272512511e-05,
                                -2.1315026809955787e-05,
                                -3.7346551751414047e-06,
                                2.0637618513646814e-06,
                                1.6744288576823017e-07,
                                -9.517657273819165e-08
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {-0.00021208083980379827,
                                0.00035858968789573785,
                                0.0021782363581090178,
                                -0.004159358781386048,
                                -0.010131117519849788,
                                0.023408156785839195,
                                0.02816802897093635,
                                -0.091920010559696244,
                                -0.052043163176243773,
                                0.42156620669085149,
                                0.77428960365295618,
                                0.43799162617183712,
                                -0.062035963962903569,
                                -0.10557420870333893,
                                0.041289208750181702,
                                0.032683574267111833,
                                -0.019761778942572639,
                                -0.0091642311624818458,
                                0.0067641854480530832,
                                0.0024333732126576722,
                                -0.0016628637020130838,
                                -0.00063813134304511142,
                                0.00030225958181306315,
                                0.00014054114970203437,
                                -4.1340432272512511e-05,
                                -2.1315026809955787e-05,
                                3.7346551751414047e-06,
                                2.0637618513646814e-06,
                                -1.6744288576823017e-07,
                                -9.517657273819165e-08
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-9.517657273819165e-08,
                                 1.6744288576823017e-07,
                                 2.0637618513646814e-06,
                                 -3.7346551751414047e-06,
                                 -2.1315026809955787e-05,
                                 4.1340432272512511e-05,
                                 0.00014054114970203437,
                                 -0.00030225958181306315,
                                 -0.00063813134304511142,
                                 0.0016628637020130838,
                                 0.0024333732126576722,
                                 -0.0067641854480530832,
                                 -0.0091642311624818458,
                                 0.019761778942572639,
                                 0.032683574267111833,
                                 -0.041289208750181702,
                                 -0.10557420870333893,
                                 0.062035963962903569,
                                 0.43799162617183712,
                                 -0.77428960365295618,
                                 0.42156620669085149,
                                 0.052043163176243773,
                                 -0.091920010559696244,
                                 -0.02816802897093635,
                                 0.023408156785839195,
                                 0.010131117519849788,
                                 -0.004159358781386048,
                                 -0.0021782363581090178,
                                 0.00035858968789573785,
                                 0.00021208083980379827
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
            }

    else if ( name == "sym2"){
                double lp1_a[] = {-0.12940952255092145,
                                  0.22414386804185735,
                                  0.83651630373746899,
                                  0.48296291314469025
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.48296291314469025,
                                  0.83651630373746899,
                                  -0.22414386804185735,
                                  -0.12940952255092145
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.48296291314469025,
                                  0.83651630373746899,
                                  0.22414386804185735,
                                  -0.12940952255092145
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.12940952255092145,
                                   -0.22414386804185735,
                                   0.83651630373746899,
                                   -0.48296291314469025

};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "sym3"){
                double lp1_a[] = {0.035226291882100656,
                                  -0.085441273882241486,
                                  -0.13501102001039084,
                                  0.45987750211933132,
                                  0.80689150931333875,
                                  0.33267055295095688

};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.33267055295095688,
                                  0.80689150931333875,
                                  -0.45987750211933132,
                                  -0.13501102001039084,
                                  0.085441273882241486,
                                  0.035226291882100656
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.33267055295095688,
                                  0.80689150931333875,
                                  0.45987750211933132,
                                  -0.13501102001039084,
                                  -0.085441273882241486,
                                  0.035226291882100656
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.035226291882100656,
                                   0.085441273882241486,
                                   -0.13501102001039084,
                                   -0.45987750211933132,
                                   0.80689150931333875,
                                   -0.33267055295095688
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "sym4"){
                double lp1_a[] = {-0.075765714789273325,
                                  -0.02963552764599851,
                                  0.49761866763201545,
                                  0.80373875180591614,
                                  0.29785779560527736,
                                  -0.099219543576847216,
                                  -0.012603967262037833,
                                  0.032223100604042702
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.032223100604042702,
                                  -0.012603967262037833,
                                  0.099219543576847216,
                                  0.29785779560527736,
                                  -0.80373875180591614,
                                  0.49761866763201545,
                                  0.02963552764599851,
                                  -0.075765714789273325
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.032223100604042702,
                                  -0.012603967262037833,
                                  -0.099219543576847216,
                                  0.29785779560527736,
                                  0.80373875180591614,
                                  0.49761866763201545,
                                  -0.02963552764599851,
                                  -0.075765714789273325
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.075765714789273325,
                                   0.02963552764599851,
                                   0.49761866763201545,
                                   -0.80373875180591614,
                                   0.29785779560527736,
                                   0.099219543576847216,
                                   -0.012603967262037833,
                                   -0.032223100604042702
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "sym5"){
                double lp1_a[] = {0.027333068345077982,
                                  0.029519490925774643,
                                  -0.039134249302383094,
                                  0.1993975339773936,
                                  0.72340769040242059,
                                  0.63397896345821192,
                                  0.016602105764522319,
                                  -0.17532808990845047,
                                  -0.021101834024758855,
                                  0.019538882735286728
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.019538882735286728,
                                  -0.021101834024758855,
                                  0.17532808990845047,
                                  0.016602105764522319,
                                  -0.63397896345821192,
                                  0.72340769040242059,
                                  -0.1993975339773936,
                                  -0.039134249302383094,
                                  -0.029519490925774643,
                                  0.027333068345077982
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.019538882735286728,
                                  -0.021101834024758855,
                                  -0.17532808990845047,
                                  0.016602105764522319,
                                  0.63397896345821192,
                                  0.72340769040242059,
                                  0.1993975339773936,
                                  -0.039134249302383094,
                                  0.029519490925774643,
                                  0.027333068345077982

};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.027333068345077982,
                                   -0.029519490925774643,
                                   -0.039134249302383094,
                                   -0.1993975339773936,
                                   0.72340769040242059,
                                   -0.63397896345821192,
                                   0.016602105764522319,
                                   0.17532808990845047,
                                   -0.021101834024758855,
                                   -0.019538882735286728
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "sym6"){
                double lp1_a[] = {0.015404109327027373,
                                  0.0034907120842174702,
                                  -0.11799011114819057,
                                  -0.048311742585632998,
                                  0.49105594192674662,
                                  0.787641141030194,
                                  0.3379294217276218,
                                  -0.072637522786462516,
                                  -0.021060292512300564,
                                  0.044724901770665779,
                                  0.0017677118642428036,
                                  -0.007800708325034148
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.007800708325034148,
                                  0.0017677118642428036,
                                  -0.044724901770665779,
                                  -0.021060292512300564,
                                  0.072637522786462516,
                                  0.3379294217276218,
                                  -0.787641141030194,
                                  0.49105594192674662,
                                  0.048311742585632998,
                                  -0.11799011114819057,
                                  -0.0034907120842174702,
                                  0.015404109327027373
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {-0.007800708325034148,
                                  0.0017677118642428036,
                                  0.044724901770665779,
                                  -0.021060292512300564,
                                  -0.072637522786462516,
                                  0.3379294217276218,
                                  0.787641141030194,
                                  0.49105594192674662,
                                  -0.048311742585632998,
                                  -0.11799011114819057,
                                  0.0034907120842174702,
                                  0.015404109327027373
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.015404109327027373,
                                   -0.0034907120842174702,
                                   -0.11799011114819057,
                                   0.048311742585632998,
                                   0.49105594192674662,
                                   -0.787641141030194,
                                   0.3379294217276218,
                                   0.072637522786462516,
                                   -0.021060292512300564,
                                   -0.044724901770665779,
                                   0.0017677118642428036,
                                   0.007800708325034148
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "sym7"){
                double lp1_a[] = {0.0026818145682578781,
                                  -0.0010473848886829163,
                                  -0.01263630340325193,
                                  0.03051551316596357,
                                  0.067892693501372697,
                                  -0.049552834937127255,
                                  0.017441255086855827,
                                  0.5361019170917628,
                                  0.76776431700316405,
                                  0.28862963175151463,
                                  -0.14004724044296152,
                                  -0.10780823770381774,
                                  0.0040102448715336634,
                                  0.010268176708511255
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.010268176708511255,
                                  0.0040102448715336634,
                                  0.10780823770381774,
                                  -0.14004724044296152,
                                  -0.28862963175151463,
                                  0.76776431700316405,
                                  -0.5361019170917628,
                                  0.017441255086855827,
                                  0.049552834937127255,
                                  0.067892693501372697,
                                  -0.03051551316596357,
                                  -0.01263630340325193,
                                  0.0010473848886829163,
                                  0.0026818145682578781

};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.010268176708511255,
                                  0.0040102448715336634,
                                  -0.10780823770381774,
                                  -0.14004724044296152,
                                  0.28862963175151463,
                                  0.76776431700316405,
                                  0.5361019170917628,
                                  0.017441255086855827,
                                  -0.049552834937127255,
                                  0.067892693501372697,
                                  0.03051551316596357,
                                  -0.01263630340325193,
                                  -0.0010473848886829163,
                                  0.0026818145682578781
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0026818145682578781,
                                   0.0010473848886829163,
                                   -0.01263630340325193,
                                   -0.03051551316596357,
                                   0.067892693501372697,
                                   0.049552834937127255,
                                   0.017441255086855827,
                                   -0.5361019170917628,
                                   0.76776431700316405,
                                   -0.28862963175151463,
                                   -0.14004724044296152,
                                   0.10780823770381774,
                                   0.0040102448715336634,
                                   -0.010268176708511255
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "sym8"){
                double lp1_a[] = {-0.0033824159510061256,
                                  -0.00054213233179114812,
                                  0.031695087811492981,
                                  0.0076074873249176054,
                                  -0.14329423835080971,
                                  -0.061273359067658524,
                                  0.48135965125837221,
                                  0.77718575170052351,
                                  0.3644418948353314,
                                  -0.051945838107709037,
                                  -0.027219029917056003,
                                  0.049137179673607506,
                                  0.0038087520138906151,
                                  -0.014952258337048231,
                                  -0.0003029205147213668,
                                  0.0018899503327594609
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.0018899503327594609,
                                  -0.0003029205147213668,
                                  0.014952258337048231,
                                  0.0038087520138906151,
                                  -0.049137179673607506,
                                  -0.027219029917056003,
                                  0.051945838107709037,
                                  0.3644418948353314,
                                  -0.77718575170052351,
                                  0.48135965125837221,
                                  0.061273359067658524,
                                  -0.14329423835080971,
                                  -0.0076074873249176054,
                                  0.031695087811492981,
                                  0.00054213233179114812,
                                  -0.0033824159510061256
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0018899503327594609,
                                  -0.0003029205147213668,
                                  -0.014952258337048231,
                                  0.0038087520138906151,
                                  0.049137179673607506,
                                  -0.027219029917056003,
                                  -0.051945838107709037,
                                  0.3644418948353314,
                                  0.77718575170052351,
                                  0.48135965125837221,
                                  -0.061273359067658524,
                                  -0.14329423835080971,
                                  0.0076074873249176054,
                                  0.031695087811492981,
                                  -0.00054213233179114812,
                                  -0.0033824159510061256
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {-0.0033824159510061256,
                                   0.00054213233179114812,
                                   0.031695087811492981,
                                   -0.0076074873249176054,
                                   -0.14329423835080971,
                                   0.061273359067658524,
                                   0.48135965125837221,
                                   -0.77718575170052351,
                                   0.3644418948353314,
                                   0.051945838107709037,
                                   -0.027219029917056003,
                                   -0.049137179673607506,
                                   0.0038087520138906151,
                                   0.014952258337048231,
                                   -0.0003029205147213668,
                                   -0.0018899503327594609
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "sym9"){
                double lp1_a[] = {0.0014009155259146807,
                                  0.00061978088898558676,
                                  -0.013271967781817119,
                                  -0.01152821020767923,
                                  0.03022487885827568,
                                  0.00058346274612580684,
                                  -0.054568958430834071,
                                  0.238760914607303,
                                  0.717897082764412,
                                  0.61733844914093583,
                                  0.035272488035271894,
                                  -0.19155083129728512,
                                  -0.018233770779395985,
                                  0.06207778930288603,
                                  0.0088592674934004842,
                                  -0.010264064027633142,
                                  -0.00047315449868008311,
                                  0.0010694900329086053
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {-0.0010694900329086053,
                                  -0.00047315449868008311,
                                  0.010264064027633142,
                                  0.0088592674934004842,
                                  -0.06207778930288603,
                                  -0.018233770779395985,
                                  0.19155083129728512,
                                  0.035272488035271894,
                                  -0.61733844914093583,
                                  0.717897082764412,
                                  -0.238760914607303,
                                  -0.054568958430834071,
                                  -0.00058346274612580684,
                                  0.03022487885827568,
                                  0.01152821020767923,
                                  -0.013271967781817119,
                                  -0.00061978088898558676,
                                  0.0014009155259146807
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {0.0010694900329086053,
                                  -0.00047315449868008311,
                                  -0.010264064027633142,
                                  0.0088592674934004842,
                                  0.06207778930288603,
                                  -0.018233770779395985,
                                  -0.19155083129728512,
                                  0.035272488035271894,
                                  0.61733844914093583,
                                  0.717897082764412,
                                  0.238760914607303,
                                  -0.054568958430834071,
                                  0.00058346274612580684,
                                  0.03022487885827568,
                                  -0.01152821020767923,
                                  -0.013271967781817119,
                                  0.00061978088898558676,
                                  0.0014009155259146807
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.0014009155259146807,
                                   -0.00061978088898558676,
                                   -0.013271967781817119,
                                   0.01152821020767923,
                                   0.03022487885827568,
                                   -0.00058346274612580684,
                                   -0.054568958430834071,
                                   -0.238760914607303,
                                   0.717897082764412,
                                   -0.61733844914093583,
                                   0.035272488035271894,
                                   0.19155083129728512,
                                   -0.018233770779395985,
                                   -0.06207778930288603,
                                   0.0088592674934004842,
                                   0.010264064027633142,
                                   -0.00047315449868008311,
                                   -0.0010694900329086053
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }

    else if ( name == "sym10"){
                double lp1_a[] = {0.00077015980911449011,
                                  9.5632670722894754e-05,
                                  -0.0086412992770224222,
                                  -0.0014653825813050513,
                                  0.045927239231092203,
                                  0.011609893903711381,
                                  -0.15949427888491757,
                                  -0.070880535783243853,
                                  0.47169066693843925,
                                  0.7695100370211071,
                                  0.38382676106708546,
                                  -0.035536740473817552,
                                  -0.0319900568824278,
                                  0.049994972077376687,
                                  0.0057649120335819086,
                                  -0.02035493981231129,
                                  -0.00080435893201654491,
                                  0.0045931735853118284,
                                  5.7036083618494284e-05,
                                  -0.00045932942100465878
};
                lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                double hp1_a[] = {0.00045932942100465878,
                                  5.7036083618494284e-05,
                                  -0.0045931735853118284,
                                  -0.00080435893201654491,
                                  0.02035493981231129,
                                  0.0057649120335819086,
                                  -0.049994972077376687,
                                  -0.0319900568824278,
                                  0.035536740473817552,
                                  0.38382676106708546,
                                  -0.7695100370211071,
                                  0.47169066693843925,
                                  0.070880535783243853,
                                  -0.15949427888491757,
                                  -0.011609893903711381,
                                  0.045927239231092203,
                                  0.0014653825813050513,
                                  -0.0086412992770224222,
                                  -9.5632670722894754e-05,
                                  0.00077015980911449011
};
                 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                double lp2_a[] = {-0.00045932942100465878,
                                  5.7036083618494284e-05,
                                  0.0045931735853118284,
                                  -0.00080435893201654491,
                                  -0.02035493981231129,
                                  0.0057649120335819086,
                                  0.049994972077376687,
                                  -0.0319900568824278,
                                  -0.035536740473817552,
                                  0.38382676106708546,
                                  0.7695100370211071,
                                  0.47169066693843925,
                                  -0.070880535783243853,
                                  -0.15949427888491757,
                                  0.011609893903711381,
                                  0.045927239231092203,
                                  -0.0014653825813050513,
                                  -0.0086412992770224222,
                                  9.5632670722894754e-05,
                                  0.00077015980911449011
};
                 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                 double hp2_a[] = {0.00077015980911449011,
                                   -9.5632670722894754e-05,
                                   -0.0086412992770224222,
                                   0.0014653825813050513,
                                   0.045927239231092203,
                                   -0.011609893903711381,
                                   -0.15949427888491757,
                                   0.070880535783243853,
                                   0.47169066693843925,
                                   -0.7695100370211071,
                                   0.38382676106708546,
                                   0.035536740473817552,
                                   -0.0319900568824278,
                                   -0.049994972077376687,
                                   0.0057649120335819086,
                                   0.02035493981231129,
                                   -0.00080435893201654491,
                                   -0.0045931735853118284,
                                   5.7036083618494284e-05,
                                   0.00045932942100465878
};
                 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                 return 0;
        }


    else {
        cout << "Filter Not in Database" << endl;
        return -1;
    }

}




