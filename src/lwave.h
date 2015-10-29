//============================================================================
// Name : 1D/2D Wavelet Transform
// Author : Rafat Hussain
// Version :
// Copyright : GNU GPL License
// Description : LiftWave Wavelet Library Component
//============================================================================
/*
* Copyright (c) 2012 Rafat Hussain
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*/

#ifndef LWAVE_H
#define LWAVE_H
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "lift.h"
#include "alg.h"

using namespace std;

template <class T>
struct IsInt 
{
	static const bool value = false;
};

template <>
struct IsInt<int> 
{
	static const bool value = true;
};

template <>
struct IsInt<short> 
{
	static const bool value = true;
};

template <>
struct IsInt<long> 
{
	static const bool value = true;
};

template <typename T>
class lwt {
	vector<T> cA,cD;
	vector<int> cD_length;
	int level;

public:
    lwt(vector<T> &signal, liftscheme &lft){
	level=1;	
	vector<double> coeff;
	vector<int> lenv;
	string lat;
	double K;
	lft.getScheme(coeff,lenv,lat,K);
	
	// Number Of Liftin Stages N
	int N = lat.size();
	//vector<T> sl,dl;
	//split(signal,sl,dl);
    vector<T> sl=signal;
    vector<T> dl=signal;
	int cume_coeff=0;
	
	for (int i=0; i < N ; i++) {
		char lft_type = lat.at(i);
		vector<double> filt;
		int len_filt = lenv[2*i];
		int max_pow = lenv[2*i+1];
		
		for (int j=0; j < len_filt; j++) {
			filt.push_back(coeff[cume_coeff+j]);
		}
		cume_coeff=cume_coeff+len_filt;
		
		if (lft_type == 'd') {
			
			for (int len_dl = 0; len_dl < (int) dl.size();len_dl++) {
				double temp = 0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_dl+max_pow-lf) >= 0 && (len_dl+max_pow-lf) < (int) sl.size()) {
						temp=temp+filt[lf]*sl[len_dl+max_pow-lf];
					
					}
				}
				dl[len_dl]=dl[len_dl]-(T) temp;
				
			}
		
			
		} else if (lft_type == 'p') {
			
			for (int len_sl = 0; len_sl < (int) dl.size();len_sl++) {
				double temp = 0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_sl+max_pow-lf) >= 0 && (len_sl+max_pow-lf) < (int) dl.size()) {
						temp=temp+filt[lf]*dl[len_sl+max_pow-lf];
					
					}
				}
				sl[len_sl]=sl[len_sl]-(T) temp;
				
			}
			
		}
		
	}
	double K1 = 1.0/K;
	if ( !IsInt<T>::value) {
		vecmult(sl,K1);
		vecmult(dl,K);
	}

	cA=sl;
	cD=dl;
	cD_length.clear();
	cD_length.push_back((int) cD.size());
		
	}
	
	lwt(vector<T> &signal, string &name){
	level=1;	
	liftscheme lft(name);
	lwt<T> wavelift(signal,lft);
	vector<T> sx,dx;
	wavelift.getCoeff(sx,dx);
	cA=sx;
	cD=dx;
	cD_length.clear();
	cD_length.push_back((int) cD.size());
	}
	
	lwt(vector<T> &signal, liftscheme &lft, int &J) {
	/*	int Max_Iter;
		Max_Iter = (int) ceil(log( double(signal.size()))/log (2.0)) - 1;

		if ( Max_Iter < J) {
			J = Max_Iter;

		}*/
		
		vector<T> temp=signal;
		vector<T> det,temp2;
		vector<int> len_det;
		
		for (int iter=0; iter < J; iter++) {
			lwt jlevel(temp,lft);
			jlevel.getCoeff(temp,temp2);
			int len_d = temp2.size();
			det.insert(det.begin(),temp2.begin(),temp2.end());
			len_det.insert(len_det.begin(),len_d);
		}
		cA=temp;
		cD=det;
		cD_length=len_det;
		level = J;
		
	}
	
	lwt(vector<T> &signal, string &name, int &J){
	liftscheme lft(name);
	lwt<T> wavelift(signal,lft,J);
	vector<T> sx,dx;
	wavelift.getCoeff(sx,dx);
	cA=sx;
	cD=dx;
	vector<int> cdlen;
	wavelift.getDetailVec(cdlen);
	cD_length=cdlen;
	level=J;
	}
	
void getCoeff(vector<T> &appx, vector<T> &det) {
	appx = cA;
	det = cD;
}

void getDetailVec(vector<int> &detvec) {
	detvec=cD_length;
}

int getLevels() {
	return level;
}
	
	
	virtual ~lwt()
	{
	}

};

template <typename T>
class ilwt {
	vector<T> signal;

public:
	ilwt(vector<T> &sl, vector<T> &dl, liftscheme &lft){
	vector<double> coeff;
	vector<int> lenv;
	string lat;
	double K;
	lft.getScheme(coeff,lenv,lat,K);
	//vector<T> sl,dl;
	//sl=cA;
	//dl=cD;
	double K1 = 1.0/K;
	if ( !IsInt<T>::value) {
		vecmult(sl,K);
		vecmult(dl,K1);
	}
	
	int N = lat.size();
	
	int cume_coeff=coeff.size();
	
	for (int i=N-1; i >= 0 ; i--) {
		char lft_type = lat.at(i);
		vector<double> filt;
		int len_filt = lenv[2*i];
		int max_pow = lenv[2*i+1];
		
		cume_coeff=cume_coeff-len_filt;
		
		for (int j=0; j < len_filt; j++) {
			filt.push_back(coeff[cume_coeff+j]);
		}
		
		
		if (lft_type == 'd') {
			
			for (int len_dl = 0; len_dl < (int) dl.size();len_dl++) {
				double temp =  0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_dl+max_pow-lf) >= 0 && (len_dl+max_pow-lf) < (int) sl.size()) {
						temp=temp+filt[lf]*sl[len_dl+max_pow-lf];
					
					}
				}
				dl[len_dl]=dl[len_dl]+(T) temp;
				
			}
			
		} else if (lft_type == 'p') {
			
			for (int len_sl = 0; len_sl < (int) dl.size();len_sl++) {
				double temp =  0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_sl+max_pow-lf) >= 0 && (len_sl+max_pow-lf) < (int) dl.size()) {
						temp=temp+filt[lf]*dl[len_sl+max_pow-lf];
					
					}
				}
				sl[len_sl]=sl[len_sl]+(T) temp;
				
			}
			
		}
		
	}
	vector<T> idwt_oup;
	merge(idwt_oup,sl,dl);
	
	signal=idwt_oup;
	
	}
	
	ilwt(vector<T> &sl,vector<T> &dl, string &name){
	liftscheme lft(name);
	ilwt<T> wavelift(sl,dl,lft);
	vector<T> sigx;
	wavelift.getSignal(sigx);
	signal=sigx;
	}
	
	ilwt(lwt<T> &wt,liftscheme &lft) {
	int J=wt.getLevels();
	vector<T> sl,dl;
	wt.getCoeff(sl,dl);
	vector<int> detv;
	
	wt.getDetailVec(detv);
	int total=0;
	
	for (int i=0; i < J; i++) {
		vector<T> temp,temp2;
		for (int j=0; j < (int) detv[i]; j++) {
			temp.push_back(dl[total+j]);
		}
		total=total+(int) detv[i];
		ilwt<T> iwt(sl,temp,lft);
		iwt.getSignal(temp2);
		sl=temp2;
		
	}
	signal=sl;
	}
	
	ilwt(lwt<T> &wt, string &name){
	liftscheme lft(name);
	ilwt<T> wavelift(wt,lft);
	vector<T> sigx;
	wavelift.getSignal(sigx);
	signal=sigx;
	}
	
void getSignal(vector<T> &sig) {
	sig=signal;
}
	
	virtual ~ilwt()
	{
	}
};

template <typename T>
class lwt2 {
	vector<T> cLL,cLH,cHL,cHH;
	int rowLL,colLL;
	int rowLH,colLH;
	int rowHL,colHL;
	int rowHH,colHH;
	int level;
	vector<int> coef_lengths;
	
public:
    lwt2(vector<T> &signal,int rows, int cols, liftscheme &lft) {
		vector<T> L,H;
		int rows_L,cols_L,rows_H,cols_H;
		rows_L=rows;
		rows_H=rows;
		for (int i=0; i < rows; i++) {
			vector<T> temp;
			temp.assign(signal.begin()+i*cols,signal.begin()+(i+1)*cols);
			lwt<T> lwt1(temp,lft);
			vector<T> a,d;
	        lwt1.getCoeff(a,d);
			L.insert(L.end(),a.begin(),a.end());
			H.insert(H.end(),d.begin(),d.end());
			if (i==0) {
				cols_L=a.size();
				cols_H=d.size();
			}
			
		}
		
		vector<T> LT,HT;
		transpose(L,rows_L,cols_L,LT);
		transpose(H,rows_H,cols_H,HT);
		int rows_ll,cols_ll,rows_lh,cols_lh;
		
		// Remove cout
		
		
		vector<T> LL,LH;
		
		// Low Pass Stage
		cols_ll=cols_L;
		cols_lh=cols_L;

		for (int i=0; i < cols_L; i++) {
			vector<T> temp;
			temp.assign(LT.begin()+i*rows_L,LT.begin()+(i+1)*rows_L);
			lwt<T> lwt1(temp,lft);
			vector<T> a,d;
	        lwt1.getCoeff(a,d);
			LL.insert(LL.end(),a.begin(),a.end());
			LH.insert(LH.end(),d.begin(),d.end());
			if (i==0) {
				rows_ll=a.size();
				rows_lh=d.size();
			}
			
		}
		
		
		int rows_hl,cols_hl,rows_hh,cols_hh;
		
		vector<T> HL,HH;
		
		// High Pass Stage
		cols_hl=cols_H;
		cols_hh=cols_H;
		for (int i=0; i < cols_H; i++) {
			vector<T> temp;
			temp.assign(HT.begin()+i*rows_H,HT.begin()+(i+1)*rows_H);
			lwt<T> lwt1(temp,lft);
			vector<T> a,d;
	        lwt1.getCoeff(a,d);
			HL.insert(HL.end(),a.begin(),a.end());
			HH.insert(HH.end(),d.begin(),d.end());
			if (i==0) {
				rows_hl=a.size();
				rows_hh=d.size();
			}
			
		}
		
		
		
		//cLL=LL;cLH=LH;cHL=HL;cHH=HH;
		transpose(LL,cols_ll,rows_ll,cLL);
		transpose(LH,cols_lh,rows_lh,cLH);
		transpose(HL,cols_hl,rows_hl,cHL);
		transpose(HH,cols_hh,rows_hh,cHH);
		rowLL=rows_ll;rowLH=rows_lh;
		rowHL=rows_hl;rowHH=rows_hh;
		colLL=cols_ll;colLH=cols_lh;
		colHL=cols_hl;colHH=cols_hh;
		level=1;
		coef_lengths.push_back(rowLL);
		coef_lengths.push_back(colLL);
		coef_lengths.push_back(rowLH);
		coef_lengths.push_back(colLH);
		coef_lengths.push_back(rowHL);
		coef_lengths.push_back(colHL);
		coef_lengths.push_back(rowHH);
		coef_lengths.push_back(colHH);
		
		
	}	
	
	lwt2(vector<T> &signal,int rows, int cols, liftscheme &lft, int J) {
		vector<T> A1,B1,C1,D1;
		vector<int> siglen;
		
		for (int i=0; i < J; i++) {
			lwt2<T> wt2(signal, rows, cols,lft);
			vector<T> tempA,tempB,tempC,tempD;
			vector<int> temp_siglen;
			wt2.getCoef(tempA,tempB,tempC,tempD);
			signal=tempA;
			A1=tempA;
			B1.insert(B1.begin(),tempB.begin(),tempB.end());
			C1.insert(C1.begin(),tempC.begin(),tempC.end());
			D1.insert(D1.begin(),tempD.begin(),tempD.end());
			
			wt2.getDim(temp_siglen);
			rows=temp_siglen[0];
			cols=temp_siglen[1];
			if (i==J-1) {
				siglen.insert(siglen.begin(),temp_siglen.begin(),temp_siglen.end());
			} else {
				siglen.insert(siglen.begin(),temp_siglen.begin()+2,temp_siglen.end());
			}
			
			
			
		}
		cLL=A1;
		cLH=B1;
		cHL=C1;
		cHH=D1;
		level = J;
		coef_lengths=siglen;
		
		
		
	}

    //entered by Daniel
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
	
void getCoef(vector<T> &aLL, vector<T> &aLH, vector<T> &aHL, vector<T> &aHH) {
	aLL=cLL;
	aLH=cLH;
	aHL=cHL;
	aHH=cHH;
	
}	

void getDim(vector<int> &dimvec) {
	
	dimvec = coef_lengths;
}

int getLevels() {
	return level;
}

void getDetails(string align, int slevel, vector<T> &det_vec, vector<int> &det_len) {
	int J = level;
	int lev;
	if (slevel > J) {
		cout << " Decomposition has only " << J << " levels" << endl;
		exit(1);
	} else {
		lev=J-slevel;
	}
	
	vector<int> sig_vec = coef_lengths;
	vector<T> A1,B1,C1,D1;
	A1=cLL;
	B1=cLH;
	C1=cHL;
	D1=cHH;
	int total=0;
	
	if (align == "LH" || align == "lh") {
		
		for (int i=0; i < lev; i++) {
			total=total+sig_vec[2+i*6]*sig_vec[3+i*6];			
		}
		det_vec.assign(B1.begin()+total,B1.begin()+total+sig_vec[2+lev*6]*sig_vec[3+lev*6]);
		det_len.push_back(sig_vec[2+lev*6]);
		det_len.push_back(sig_vec[3+lev*6]);
		
	} else if (align == "HL" || align == "hl") {
		
		for (int i=0; i < lev; i++) {
			total=total+sig_vec[4+i*6]*sig_vec[5+i*6];			
		}
		det_vec.assign(C1.begin()+total,C1.begin()+total+sig_vec[4+lev*6]*sig_vec[5+lev*6]);
		det_len.push_back(sig_vec[4+lev*6]);
		det_len.push_back(sig_vec[5+lev*6]);
		
	} else if (align == "HH" || align == "hh") {
		
		for (int i=0; i < lev; i++) {
			total=total+sig_vec[6+i*6]*sig_vec[7+i*6];			
		}
		det_vec.assign(D1.begin()+total,D1.begin()+total+sig_vec[6+lev*6]*sig_vec[7+lev*6]);
		det_len.push_back(sig_vec[6+lev*6]);
		det_len.push_back(sig_vec[7+lev*6]);
		
	} else {
		cout << "Accepted filter stages are LH or lh, HL or hl and HH or hh" << endl;
		exit(1);
	}
	
}	
	
	virtual ~lwt2() {
		
	}
};

template <typename T>
class ilwt2 {
	vector<T> signal;
	int oup_row, oup_col;
	
public:
ilwt2(vector<T> &A,vector<T> &H,vector<T> &V,vector<T> &D,vector<int> &length,liftscheme &lft) {
		
		int cols_L=length[1];
		
		int rows_LL=length[0];
		int rows_LH=length[2];
		
		vector<T> AT,HT,VT,DT;
		
		transpose(A,length[0],length[1],AT);
		transpose(H,length[2],length[3],HT);
		transpose(V,length[4],length[5],VT);
		transpose(D,length[6],length[7],DT);
		vector<T> L;
		int rows_L;
		
		for (int i=0; i < cols_L; i++) {
			vector<T> temp1,temp2;
			temp1.assign(AT.begin()+i*rows_LL,AT.begin()+(i+1)*rows_LL);
			temp2.assign(HT.begin()+i*rows_LH,HT.begin()+(i+1)*rows_LH);
			ilwt<T> iwt(temp1,temp2,lft);
			vector<T> sig;
			iwt.getSignal(sig);
			L.insert(L.end(),sig.begin(),sig.end());
			if (i==0) {
				rows_L=(int) sig.size();
			}
			
			
		}
		
		vector<T> L2;
		transpose(L,cols_L,rows_L,L2);
		
		int cols_H=length[5];
		
		int rows_HL=length[4];
		int rows_HH=length[6];
		vector<T> H1;
		int rows_H;
		
		for (int i=0; i < cols_H; i++) {
			vector<T> temp1,temp2;
			temp1.assign(VT.begin()+i*rows_HL,VT.begin()+(i+1)*rows_HL);
			temp2.assign(DT.begin()+i*rows_HH,DT.begin()+(i+1)*rows_HH);
			ilwt<T> iwt(temp1,temp2,lft);
			vector<T> sig;
			iwt.getSignal(sig);
			H1.insert(H1.end(),sig.begin(),sig.end());
			if (i==0) {
				rows_H=(int) sig.size();
			}
			
			
		}
		
		vector<T> H2;
		transpose(H1,cols_H,rows_H,H2);
		
		vector<T> oup;
		int cx;
		
		for (int i=0; i < rows_L; i++) {
			vector<T> temp1,temp2;
			temp1.assign(L2.begin()+i*cols_L,L2.begin()+(i+1)*cols_L);
			temp2.assign(H2.begin()+i*cols_H,H2.begin()+(i+1)*cols_H);
			ilwt<T> iwt(temp1,temp2,lft);
			vector<T> sig;
			iwt.getSignal(sig);
			oup.insert(oup.end(),sig.begin(),sig.end());
			if (i==0) {
				cx=(int) sig.size();
			}
			
			
		}
		
		signal=oup;
		oup_row=rows_L;
		oup_col=cx;
		
		
	}	
    
   /* ilwt2(lwt2<T> &wt,liftscheme &lft) {
		vector<T> A,H,V,D;
		wt.getCoef(A,H,V,D);
		vector<int> length;
		wt.getDim(length);
		
		int cols_L=length[1];
		
		int rows_LL=length[0];
		int rows_LH=length[2];
		
		vector<T> AT,HT,VT,DT;
		
		transpose(A,length[0],length[1],AT);
		transpose(H,length[2],length[3],HT);
		transpose(V,length[4],length[5],VT);
		transpose(D,length[6],length[7],DT);
		vector<T> L;
		int rows_L;
		
		for (int i=0; i < cols_L; i++) {
			vector<T> temp1,temp2;
			temp1.assign(AT.begin()+i*rows_LL,AT.begin()+(i+1)*rows_LL);
			temp2.assign(HT.begin()+i*rows_LH,HT.begin()+(i+1)*rows_LH);
			ilwt<T> iwt(temp1,temp2,lft);
			vector<T> sig;
			iwt.getSignal(sig);
			L.insert(L.end(),sig.begin(),sig.end());
			if (i==0) {
				rows_L=(int) sig.size();
			}
			
			
		}
		
		vector<T> L2;
		transpose(L,cols_L,rows_L,L2);
		
		int cols_H=length[5];
		
		int rows_HL=length[4];
		int rows_HH=length[6];
		vector<T> H1;
		int rows_H;
		
		for (int i=0; i < cols_H; i++) {
			vector<T> temp1,temp2;
			temp1.assign(VT.begin()+i*rows_HL,VT.begin()+(i+1)*rows_HL);
			temp2.assign(DT.begin()+i*rows_HH,DT.begin()+(i+1)*rows_HH);
			ilwt<T> iwt(temp1,temp2,lft);
			vector<T> sig;
			iwt.getSignal(sig);
			H1.insert(H1.end(),sig.begin(),sig.end());
			if (i==0) {
				rows_H=(int) sig.size();
			}
			
			
		}
		
		vector<T> H2;
		transpose(H1,cols_H,rows_H,H2);
		
		vector<T> oup;
		int cx;
		
		for (int i=0; i < rows_L; i++) {
			vector<T> temp1,temp2;
			temp1.assign(L2.begin()+i*cols_L,L2.begin()+(i+1)*cols_L);
			temp2.assign(H2.begin()+i*cols_H,H2.begin()+(i+1)*cols_H);
			ilwt<T> iwt(temp1,temp2,lft);
			vector<T> sig;
			iwt.getSignal(sig);
			oup.insert(oup.end(),sig.begin(),sig.end());
			if (i==0) {
				cx=(int) sig.size();
			}
			
			
		}
		
		signal=oup;
		oup_row=rows_L;
		oup_col=cx;
		
		
	}	*/
	
	ilwt2(lwt2<T> &wt,liftscheme &lft) {
		int J = wt.getLevels();
		vector<T> A1,B1,C1,D1;
		wt.getCoef(A1,B1,C1,D1);
		vector<int> len_coef;
		wt.getDim(len_coef);
		int slevel;
		
		int count = 0;
		vector<T> A,H,V,D;
		A=A1;
		vector<int> rxcx;
		
		for (int i=0; i < J; i++) {
			vector<int> temp_coef,detlenH,detlenV,detlenD;
			slevel = J-i;
			if (i==0) {
				temp_coef.assign(len_coef.begin(),len_coef.begin()+8);
				count=count+8;
			} else {
				temp_coef.assign(len_coef.begin()+count,len_coef.begin()+count+6);
				temp_coef.insert(temp_coef.begin(),rxcx.begin(),rxcx.end());
				count=count+6;
				
			}
			// Get LH,HL and HH coefficients
			getDetails(wt,"LH",slevel,H,detlenH);
			getDetails(wt,"HL",slevel,V,detlenV);
			getDetails(wt,"HH",slevel,D,detlenD);

			ilwt2<T> iwt2(A,H,V,D,temp_coef,lft);
			vector<T> temp;
			iwt2.getSignal(temp);
			rxcx.clear();
			iwt2.getDim(rxcx);
			oup_row=rxcx[0];
			oup_col=rxcx[1];
			A=temp;
		}
		    signal=A;
			
			
	}

void getDetails(lwt2<T> &wt, string align, int slevel, vector<T> &det_vec, vector<int> det_len) {
	int J = wt.getLevels();
	int lev;
	if (slevel > J) {
		cout << " Decomposition has only " << J << " levels" << endl;
		exit(1);
	} else {
		lev=J-slevel;
	}
	
	vector<int> sig_vec;
	wt.getDim(sig_vec);
	vector<T> A1,B1,C1,D1;
	wt.getCoef(A1,B1,C1,D1);
	int total=0;
	
	if (align == "LH" || align == "lh") {
		
		for (int i=0; i < lev; i++) {
			total=total+sig_vec[2+i*6]*sig_vec[3+i*6];			
		}
		det_vec.assign(B1.begin()+total,B1.begin()+total+sig_vec[2+lev*6]*sig_vec[3+lev*6]);
		det_len.push_back(sig_vec[2+lev*6]);
		det_len.push_back(sig_vec[3+lev*6]);
		
	} else if (align == "HL" || align == "hl") {
		
		for (int i=0; i < lev; i++) {
			total=total+sig_vec[4+i*6]*sig_vec[5+i*6];			
		}
		det_vec.assign(C1.begin()+total,C1.begin()+total+sig_vec[4+lev*6]*sig_vec[5+lev*6]);
		det_len.push_back(sig_vec[4+lev*6]);
		det_len.push_back(sig_vec[5+lev*6]);
		
	} else if (align == "HH" || align == "hh") {
		
		for (int i=0; i < lev; i++) {
			total=total+sig_vec[6+i*6]*sig_vec[7+i*6];			
		}
		det_vec.assign(D1.begin()+total,D1.begin()+total+sig_vec[6+lev*6]*sig_vec[7+lev*6]);
		det_len.push_back(sig_vec[6+lev*6]);
		det_len.push_back(sig_vec[7+lev*6]);
		
	} else {
		cout << "Accepted filter stages are LH or lh, HL or hl and HH or hh" << endl;
		exit(1);
	}
	
}	

void getSignal(vector<T> &ilwt_oup) {
	ilwt_oup=signal;
}	

void getDim(vector<int> &sigdim) {
	sigdim.push_back(oup_row);
	sigdim.push_back(oup_col);
}
	
	virtual ~ilwt2() {
		
	}
};

#endif // LWAVE_H
