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

#ifndef LAURENT_H_
#define LAURENT_H_
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

template <typename T>
class Laurent {
	int deg,highest;
	vector<T> poly;

public:
	Laurent() {
	int degree=0;
	vector<T> poly;
	int highest = 0;
	// TODO Auto-generated constructor stub

}
	void setPoly(const vector<T> &coef_i, int hdeg) {
	vector<T> coef;
	coef=coef_i;
	
	int i=0;
	
	while  ( (coef[i] == 0 || abs(coef[i]) < 1e-05) && coef.size() > 1 ) {
		coef.erase(coef.begin());
		hdeg--;
		
	}
	
	if (coef.empty()) {
		coef.push_back((T) 0.0);
		hdeg=0;
		deg=0;
	}
	
	i=coef.size();
	
	while ((coef[i-1] == 0 || abs(coef[i-1]) < 1e-05) && coef.size() > 1 ) {
		coef.erase(coef.end()-1);
		i=coef.size();
	}
	
	if (coef.empty()) {
		coef.push_back((T) 0.0);
		hdeg=0;
		deg=0;
	}
	
	deg = coef.size()-1;
	poly = coef;
	highest = hdeg;
	

}
    void getPoly(vector<T> &vec) {
        vec=poly;
    }

	void dispPoly() {
	int sz=poly.size();

	for (int i = 0;i < sz; i++) {
		if (i == sz-1) {
			cout << poly[i] << "*z^" << "(" << highest-i << ")";
		}
		else if (poly[i+1] >= 0) {
			cout << poly[i] << "*z^" << "(" << highest-i << ")" << "+";
		}
		else if (poly[i+1] < 0) {
			cout << poly[i] << "*z^" << "(" << highest-i << ")" ;
		}
	}
	cout << endl;

}

	void printPoly() {
	int sz=poly.size();
	ofstream pwrite("poly.txt",ios::app);
	
    pwrite << "{";
	for (int i = 0;i < sz; i++) {
		if (i == sz-1) {
			pwrite << poly[i] << ",";
		}
		else if (poly[i+1] >= 0) {
			pwrite << poly[i] << ",";
		}
		else if (poly[i+1] < 0) {
			pwrite << poly[i] << ",";
		}
	}
	pwrite << "}" << endl;
	pwrite.close();

}
	int degree()  {
	return deg;
}
	int highdeg() {
	return highest;
}
	void LaurentAdd(Laurent &A,Laurent &B){
	int ha=A.highdeg();
	int hb=B.highdeg();

	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	vector<T> coefA=A.poly;
	vector<T> coefB=B.poly;

	int lc,hc;

	if (ha > hb) {
		hc=ha;
		} else {
		hc=hb;
		}

	if (ha - A.deg < hb - B.deg ) {
		lc=ha-A.deg;
	} else {
		lc=hb-B.deg;
	}

	vector<T> coef_c;

	for (int i=0; i < hc-lc+1; i++) {
		coef_c.push_back(0);
	}

	for (int i=0; i < lenA; i++) {
		coef_c[hc-ha+i]=coef_c[hc-ha+i]+coefA[i];
	}

	for (int i=0; i < lenB; i++) {
		coef_c[hc-hb+i]=coef_c[hc-hb+i]+coefB[i];
	}



    setPoly(coef_c,hc);

}
	void LaurentSub(Laurent &A,Laurent &B){
	int ha=A.highdeg();
	int hb=B.highdeg();

	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	vector<T> coefA=A.poly;
	vector<T> coefB=B.poly;

	int lc,hc;

	if (ha > hb) {
		hc=ha;
		} else {
		hc=hb;
		}

	if (ha - A.deg < hb - B.deg ) {
		lc=ha-A.deg;
	} else {
		lc=hb-B.deg;
	}

	vector<T> coef_c;

	for (int i=0; i < hc-lc+1; i++) {
		coef_c.push_back(0);
	}

	for (int i=0; i < lenA; i++) {
		coef_c[hc-ha+i]=coef_c[hc-ha+i]+coefA[i];
	}

	for (int i=0; i < lenB; i++) {
		coef_c[hc-hb+i]=coef_c[hc-hb+i]-coefB[i];
	}



    setPoly(coef_c,hc);

}
   void LaurentMult(Laurent &A, Laurent &B) {
    int la=A.highdeg();
	int lb=B.highdeg();

	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	vector<T> coefA=A.poly;
	vector<T> coefB=B.poly;
	int degA,degB;
	vector<int> degAB;
	vector<T> coefAB;

	for (int i=0; i < lenA; i++) {
	    degA=i+la;
	    for (int j=0; j < lenB; j++) {
	        degB=j+lb;
	        degAB.push_back(degA+degB);
	        coefAB.push_back(coefA[i]*coefB[j]);

	    }

	}

	int min_deg= *min_element(degAB.begin(),degAB.end());
	int max_deg= *max_element(degAB.begin(),degAB.end());

	vector<T> coef_c;

	for (int i=0; i < abs(max_deg - min_deg)+1; i++) {
		coef_c.push_back(0);
	}

	for (int i=min_deg; i < max_deg+1;i++) {
	    for (int j=0; j < lenA*lenB;j++) {
	        if (i==degAB[j]) {
	            coef_c[i-min_deg]+=coefAB[j];
             }
        }
    }

	setPoly(coef_c,min_deg);


   }
   
   void One() {
	   T temp=(T) 1.0;
	   vector<T> coef_one;
	   coef_one.push_back(temp);
	   int pow =0;
	   setPoly(coef_one,pow);
   }
   
   void Zero() {
	   T temp=(T) 0.0;
	   vector<T> coef_one;
	   coef_one.push_back(temp);
	   int pow =0;
	   setPoly(coef_one,pow);
   }
   
   int isZero() {

	vector<T> coefA=poly;
	
	int lenA= coefA.size();
	int count = 0;
	
	for (int i=0; i < lenA; i++) {
		if (coefA[i] != 0 && abs(coefA[i]) > 1e-05) {
			count++;
		}
		
	}
	
	if (count == 0) {
		return 1;
	} else {
		return 0;
	}
	   
   }
   
   void scale(T s) {
	   vector<T> coeff_s;
	   coeff_s.push_back(s);
	   Laurent<T> sc;
	   sc.setPoly(coeff_s,0);
	   LaurentMult(*this,sc);
	   
   }

   void zinv(Laurent &A) {
    int la=A.highdeg();

	//int lenA=A.degree()+1;

	vector<T> coefA=A.poly;

	la=(la-A.degree())*-1;

	reverse(coefA.begin(),coefA.end());

    setPoly(coefA,la);

   }
   
void nzinv(Laurent &A) {
    int la=A.highdeg();

	//int lenA=A.degree()+1;

	vector<T> coefA=A.poly;

	la=(la-A.degree())*-1;

	reverse(coefA.begin(),coefA.end());
	
	int N=coefA.size();
	
	for (int i=0; i< N ; i++) {
		coefA[i]= coefA[i]* (int) pow(-1.0,la-i);
	}

    setPoly(coefA,la);

   }
   
   void onzinv(Laurent &A) {
    int la=A.highdeg();

	//int lenA=A.degree()+1;

	vector<T> coefA=A.poly;

	la=(la-A.degree())*-1;

	reverse(coefA.begin(),coefA.end());
	
	int N=coefA.size();
	
	for (int i=0; i< N ; i++) {
		coefA[i]= coefA[i]* (int) pow(-1.0,i);
	}

    setPoly(coefA,la);

   }
   
   void nz(Laurent &A)  {
	   
	int la=A.highdeg();

	vector<T> coefA=A.poly;
	int N=coefA.size();
	
	for (int i=0; i<N;i++) {
		coefA[i]=coefA[i]* (int) pow(-1.0,la-i);
	}
	setPoly(coefA,la);
   }
   

   bool isMono() {
    int la=highest;

	int lenA=deg+1;

	vector<T> coefA=poly;
	bool mono;

	if (lenA == 0) {
	    mono = false;
	}

	if (lenA == 1) {
	    if (coefA[0] != 0 || abs(coefA[0]) <= 1e-05) {
	        mono = true;

	    } else {
	        mono = false;

	    }

	}

	if (lenA > 1) {
	    int j=0;
	    for (int i=0; i < lenA; i++) {
	        if (coefA[i] == 0 || abs(coefA[i]) <= 1e-05) {
	            j++;
	        }

	    }

	    if (lenA==j+1) {
	        mono = true;

	    } else {
	        mono = false;
	    }

	}

    return mono;

}

int monoDeg() {
	int mdeg;
	int la=highest;

	int lenA=deg+1;
	vector<T> coefA=poly;
	
	if (lenA == 1) {
		mdeg = la;
	} else {
		int temp=highest;
		for (int i=0; i < lenA; i++) {
			bool val = coefA[i] == 0 || abs(coefA[i]) < 1e-05;
			if (!val) {
	            mdeg=temp;
	        }
			temp--;
			
		}
	}
	
	return mdeg;
	
}

T monoCoef(int md) {
	int la=highest;

	vector<T> coefA=poly;
	T temp = coefA[la];
	return temp;
	
}

void merge(Laurent &X, Laurent &Y) {
	int lx=X.highest;
	int ly=Y.highest;
	
	vector<T> xvec=X.poly;
	vector<T> yvec=Y.poly;
	vector<T> mvec;
	int mhdeg;
	
	int lenA=X.deg+1;
	
	int xhdeg= 2*lx;
	int yhdeg=2*ly-1;
	
	if (xhdeg > yhdeg) {
	
	    for (int i=0; i < lenA; i++) {
		mvec.push_back(xvec[i]);
		mvec.push_back(yvec[i]);
		
	}
	
	mhdeg=xhdeg;
	} else {
		for (int i=0; i < lenA; i++) {
		mvec.push_back(yvec[i]);
		mvec.push_back(xvec[i]);
		
	}
	mhdeg=yhdeg;
	
	}
	
	setPoly(mvec,mhdeg);
}
      

    virtual ~Laurent(){
    }
};

template <typename T>
class LaurentMat {
    Laurent<T> A,B,C,D;


public:
	LaurentMat() {

	}

	void setMat(Laurent<T> &AA, Laurent<T> &BB, Laurent<T> &CC, Laurent<T> &DD) {
	A=AA;
	B=BB;
	C=CC;
	D=DD;
	}

	void Det(Laurent<T> &oup) {
	    Laurent<T> tempMat1,tempMat2;
	    tempMat1.LaurentMult(A,D);
	    tempMat2.LaurentMult(C,B);
	    tempMat1.LaurentSub(tempMat1,tempMat2);
	    oup=tempMat1;

	}
	
	void MatAdd(LaurentMat &X, LaurentMat &Y) {
		Laurent<T> oupA,oupB,oupC,oupD;
		
		oupA.LaurentAdd(X.A,Y.A);
		oupB.LaurentAdd(X.B,Y.B);
		oupC.LaurentAdd(X.C,Y.C);
		oupD.LaurentAdd(X.D,Y.D);
		
		setMat(oupA,oupB,oupC,oupD);
		
	}
	
	void MatSub(LaurentMat &X, LaurentMat &Y) {
		Laurent<T> oupA,oupB,oupC,oupD;
		
		oupA.LaurentSub(X.A,Y.A);
		oupB.LaurentSub(X.B,Y.B);
		oupC.LaurentSub(X.C,Y.C);
		oupD.LaurentSub(X.D,Y.D);
		
		setMat(oupA,oupB,oupC,oupD);
		
	}
	
	void MatMult(LaurentMat &X, LaurentMat &Y) {
		Laurent<T> tempA,tempB;
		Laurent<T> oupA,oupB,oupC,oupD;
		
		tempA.LaurentMult(X.A,Y.A);
		tempB.LaurentMult(X.B,Y.C);
		oupA.LaurentAdd(tempA,tempB);
		
		tempA.LaurentMult(X.A,Y.B);
		tempB.LaurentMult(X.B,Y.D);
		oupB.LaurentAdd(tempA,tempB);
		
		tempA.LaurentMult(X.C,Y.A);
		tempB.LaurentMult(X.D,Y.C);
		oupC.LaurentAdd(tempA,tempB);
		

		tempA.LaurentMult(X.C,Y.B);
		tempB.LaurentMult(X.D,Y.D);
		oupD.LaurentAdd(tempA,tempB);
		
		setMat(oupA,oupB,oupC,oupD);
		
	}
	
	void TZ(Laurent<T> &t) {
		Laurent<T> t1,t2,t3,t4;
		t3=t;
		t1.One();
		t2.Zero();
		t4.One();
		setMat(t1,t2,t3,t4);
		
	}
	
	void SZ(Laurent<T> &s) {
		Laurent<T> t1,t2,t3,t4;
		t2=s;
		t1.One();
		t3.Zero();
		t4.One();
		setMat(t1,t2,t3,t4);
		
	}
	
	void dispMat() {
		cout << " Laurent Matrix" << endl;
		cout << "(1,1) : " ;
        A.dispPoly();
		cout << "(1,2) : " ;
        B.dispPoly();
		cout << "(2,1) : " ;
        C.dispPoly();
		cout << "(2,2) : " ;
        D.dispPoly();
		
	}
	
	void scale(T s) {
	   A.scale(s);
	   B.scale(s);
	   C.scale(s);
	   D.scale(s);
	   
   }
   
   void MatInv(LaurentMat& Inv) {
	   Laurent<T> det_o;
	   Det(det_o);
	   if (det_o.isMono()) {
		   T coeff_d;
		   Laurent<T> A2,B2,C2,D2;
		   A2=A;
		   B2=B;
		   C2=C;
		   D2=D;
		   int mc;
		   mc = det_o.monoDeg();
		   mc = -1 *mc;
		   vector<T> cfd;
		   det_o.getPoly(cfd);
		   coeff_d=cfd[0];
		   coeff_d=((T) 1.0)/coeff_d;
		   cfd.pop_back();
		   cfd.push_back(coeff_d);
		   
		   Laurent<T> new_sc;
		   new_sc.setPoly(cfd,mc);
		   
		   B2.scale(-1.0);
		   C2.scale(-1.0);
		   
		   A2.LaurentMult(A2,new_sc);
		   B2.LaurentMult(B2,new_sc);
		   C2.LaurentMult(C2,new_sc);
		   D2.LaurentMult(D2,new_sc);
		   
		   Inv.setMat(D2,B2,C2,A2);
		   
	   } else {
		   Laurent<T> NA;
		   NA.Zero();
		   Inv.setMat(NA,NA,NA,NA);
	   }
   }
   
   void getLpoly(Laurent<T> &P,int N) {
	   switch(N) {
		   case 1:
				P=A;
				break;
		   case 2:
				P=B;
				break;
		   case 3:
                P=C;
				break;
           case 4:
                P=D;
				break;
           default:
				cout << "getLPoly only takes values from 1 to 4 where 1 -(1,1) 2(1,2) 3 (2,1) 4 (2,2)" << endl;
                break; 
		   
	   }
   }


    virtual ~LaurentMat() {

    }

};

template <typename T>
void Div(Laurent<T> &A, Laurent<T> &B, vector<Laurent<T> > &lcont) {
    int ha=A.highdeg();
	int hb=B.highdeg();
     
	int lenA=A.degree()+1;
	int lenB=B.degree()+1;

	int la=ha-lenA+1;
	int lb=hb-lenB+1;

	int hc = ha - hb;
	int lc = la -lb;

	int lenC = lenA - lenB;

	vector<T> coefA,coefB;
	A.getPoly(coefA);
	B.getPoly(coefB);
	/*
	for (int i=0; i < (int)coefA.size(); i++) {
		if (abs(coefA[i]) <= 1e-07) {
			coefA[i]=0.0;
		} 
		
	}
	
	for (int i=0; i < (int)coefB.size(); i++) {
		if (abs(coefB[i]) <= 1e-07) {
			coefB[i]=0.0;
		} 
		
	} */
	
	/*
	int i=0;
	
	while  (abs(coefA[i]) <= 1e-07 || coefA[i] == 0) {
		coefA.erase(coefA.begin());
		ha--;
		
	}
	
	i=coefA.size();
	
	while (abs(coefA[i-1]) <=1e-07 || coefA[i-1] == 0) {
		coefA.erase(coefA.end()-1);
		i=coefA.size();
	}
	
	i=0;
	
	while  (abs(coefB[i]) <= 1e-07 || coefB[i] == 0) {
		coefB.erase(coefB.begin());
		hb--;
		
	}
	
	i=coefB.size();
	
	while (abs(coefB[i-1]) <=1e-07 || coefB[i-1] == 0) {
		coefB.erase(coefB.end()-1);
		i=coefB.size();
	}
	 */
	
	//A.setPoly(coefA,ha);
	//B.setPoly(coefB,hb);
	
	
	if (lenC > 0) {
		vector<T> coef_q1,coef_q2,coef_r1,coef_r2;
		//vector<T> temp1,temp2;
		T t1,t2;
		coef_r1=coefA;
		
		for (int i=0; i < lenC+1; i++) {
			t1=coef_r1[i]/coefB[0];
			if (i == lenC) {
				t2=coef_r1[lenA-1]/coefB[lenB-1];
				coef_r2=coef_r1;
				int k=0;
			for (int j=lenA-lenB; j < lenA; j++) {
				coef_r2[j]=coef_r2[j]-t2*coefB[k];
				k++;
			}
			     coef_q2=coef_q1;
				 coef_q2.push_back(t2);
				
			}
			coef_q1.push_back(t1);
			int k=0;
			for (int j=i; j < i+lenB; j++) {
				coef_r1[j]=coef_r1[j]-t1*coefB[k];
				k++;
			}
			
		}
		    Laurent<T> q1,q2,r1,r2;
		    q1.setPoly(coef_q1,ha-hb);
		    q2.setPoly(coef_q2,ha-hb);
		    r1.setPoly(coef_r1,ha);
		    r2.setPoly(coef_r2,ha);

		//vector<Laurent<T> > lcont;
		    lcont.push_back(q1);
		    lcont.push_back(r1);
		    lcont.push_back(q2);
			lcont.push_back(r2);
			
			coef_q1.clear(),coef_q2.clear(),coef_r1.clear(),coef_r2.clear();
			coef_r1=coefA;
			
			for (int i=0; i < lenC+1; i++) {
			t1=coef_r1[lenA-i-1]/coefB[lenB-1];
			if (i == lenC) {
				t2=coef_r1[0]/coefB[0];
				coef_r2=coef_r1;
				int k=0;
			for (int j=0; j < lenB; j++) {
				coef_r2[j]=coef_r2[j]-t2*coefB[k];
				k++;
			}
			     coef_q2=coef_q1;
				 coef_q2.push_back(t2);
				
			}
			coef_q1.push_back(t1);
			int k=0;
			for (int j=lenC-i; j < lenA-i; j++) {
				coef_r1[j]=coef_r1[j]-t1*coefB[k];
				k++;
			}
			
		}
		    reverse(coef_q1.begin(),coef_q1.end());
			reverse(coef_q2.begin(),coef_q2.end());
			
			Laurent<T> q3,q4,r3,r4;
		    q3.setPoly(coef_q1,ha-hb);
		    q4.setPoly(coef_q2,(int)la-lb+coef_q1.size()-1);
		    r3.setPoly(coef_r1,ha);
		    r4.setPoly(coef_r2,ha);

		//vector<Laurent<T> > lcont;
		    lcont.push_back(q3);
		    lcont.push_back(r3);
		    lcont.push_back(q4);
			lcont.push_back(r4);
		
		
		
	} else if (lenC == 0) {
		vector<T> coef_q1,coef_q2,coef_r1,coef_r2;
	    T temp1,temp2;
	    temp1=coefA[0]/coefB[0];
	    temp2=coefA[lenA-1]/coefB[lenB-1];

	    coef_q1.push_back(temp1);
	    coef_q2.push_back(temp2);

	    for (int i=0; i < lenA; i++) {
	        coef_r1.push_back(coefA[i]- temp1*coefB[i]);
	        coef_r2.push_back(coefA[i]- temp2*coefB[i]);

	    }
		Laurent<T> q1,q2,r1,r2;
		q1.setPoly(coef_q1,ha-hb);
		q2.setPoly(coef_q2,la-lb);
		r1.setPoly(coef_r1,ha);
		r2.setPoly(coef_r2,ha);

		//vector<Laurent<T> > lcont;
		lcont.push_back(q1);
		lcont.push_back(r1);
		lcont.push_back(q2);
		lcont.push_back(r2);

	} else {
		vector<T> coef_q;
		coef_q.push_back(0);
		Laurent<T> q;
		q.setPoly(coef_q,0);

		//vector<Laurent<T> > lcont;
		lcont.push_back(q);
		lcont.push_back(A);
		
	}

    }

template <typename T>
void EvenOdd(Laurent<T> &A,Laurent<T> &even,Laurent<T> &odd) {
	vector<T> coefA;
	//la = A.highdeg();
	A.getPoly(coefA);
	
	Laurent<T> B;
	B.nz(A);
	
	Laurent<T> ad,su;
	
	ad.LaurentAdd(A,B);
	su.LaurentSub(A,B);
	
	vector<T> coefad,coefsu;
	
	ad.getPoly(coefad);
	su.getPoly(coefsu);
	
	int maxad,maxsu;
	
	maxad=ad.highdeg();
	maxsu=su.highdeg();
	
	int i=0;
	
	while  (coefad[i] == 0.0 || coefad[i] == 0) {
		coefad.erase(coefad.begin());
		maxad--;
		
	}
	
	i=coefad.size();
	
	while (coefad[i-1] == 0.0 || coefad[i-1] == 0) {
		coefad.erase(coefad.end()-1);
		i=coefad.size();
	}
	
	i=0;
	
	while  (coefsu[i] == 0.0 || coefsu[i] == 0) {
		coefsu.erase(coefsu.begin());
		maxsu--;
		
	}
	
	i=coefsu.size();
	
	while (coefsu[i-1] == 0.0 || coefsu[i-1] == 0) {
		coefsu.erase(coefsu.end()-1);
		i=coefsu.size();
	}
	
	vector<T> coefad_n,coefsu_n;
	
	for (int i=0; i < (int) coefad.size(); i+=2) {
		coefad_n.push_back(coefad[i]/2);
	}
	
	for (int i=0; i < (int) coefsu.size(); i+=2) {
		coefsu_n.push_back(coefsu[i]/2);
	}
	
	maxad=maxad/2;
	maxsu=(maxsu+1)/2;
	
	
	even.setPoly(coefad_n,maxad);
	odd.setPoly(coefsu_n,maxsu);
	
	
}
	
	


	

#endif /* LAURENT_H_ */

