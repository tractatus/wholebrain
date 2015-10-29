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
#ifndef LIFT_H
#define LIFT_H
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "Laurent.h"

//template <typename T>
class liftscheme{
	int stages;
	string ltype;
	vector<double> lcoeff;
	vector<int> plen;
    double Kconst;
	string wname;
	
public:
	
	liftscheme(string &name) {
		wname=name;
//		vector<double> coeffs;
        if (name == "lazy" ) {
            ltype="";
			stages=0;
			Kconst=1.0000;
			
			//Stage 0
			
			
			//int pl[]={};
			//plen.assign (pl,pl + sizeof(pl)/sizeof(int));
			
		} else if (name == "haar" || name == "db1" ) {
            ltype="dp";
			stages=2;
			Kconst=0.7071;
			
			//Stage 1,2
			double d1[]={1.0000};
			double p1[]={-0.5000};
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,1,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));
			
		} else if (name == "db2") {

		    ltype="pdp";
			stages=3;
			Kconst=1.93185;
			
			//Stage 1,2,3
			double p1[]={-1.73205};
			double d1[]={0.433013,-0.0669873};
			double p2[]={1.0000};
			
			lcoeff.insert(lcoeff.begin(),p2,p2+1);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,1,1};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "db3") {

		    ltype="dpdp";
			stages=4;
			Kconst=1.9182;
			
			//Stage 1,2,3,4
			
			double d1[]={-0.412287};
			double p1[]={0.352388,-1.56514};
			double d2[]={0.492152,0.0284591};
			double p2[]={-0.38962};
			
			lcoeff.insert(lcoeff.begin(),p2,p2+1);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,0,2,1,1,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "db4") {

		    ltype="pdpdp";
			stages=5;
			Kconst=2.61312;
			
			//Stage 1,2,3,4,5
			
			double p1[]={-3.10293};
			double d1[]={0.291953,-0.0763001};
			double p2[]={-1.66253,5.19949};
			double d2[]={0.0378927,-0.00672237};
			double p3[]={0.314106};
			
			lcoeff.insert(lcoeff.begin(),p3,p3+1);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,2,2,2,-2,1,3};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "db5") {

		    ltype="dpdpdp";
			stages=6;
			Kconst=1.23144;
			
			//Stage 1,2,3,4,5,6
			
			double d1[]={-0.265145};
			double p1[]={0.247729,-0.878163};
			double d2[]={0.534125,0.241421};
			double p2[]={0.1985336258386243,-0.6332784120192370};
			double d3[]={-0.0877884834474499,0.0137333394082371};
			double p3[]={-0.0315951369981596};
			
			lcoeff.insert(lcoeff.begin(),p3,p3+1);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,0,2,1,2,1,2,-1,1,2};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "db6") {

		    ltype="pdpdpdp";
			stages=7;
			Kconst=3.0899;
			
			//Stage 1,2,3,4,5,6,7
			
			double p1[]={-4.43447};
			double d1[]={0.214593,-0.0633132};
			double p2[]={-4.49311,9.97002};
			double d2[]={0.17965,0.400695};
			double p3[]={-0.0430454,0.00550295};
			double d3[]={-0.122882,-0.428777};
			double p4[]={0.0922131,-0.668385,2.33222};
			
			lcoeff.insert(lcoeff.begin(),p4,p4+3);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p3,p3+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,2,2,2,-2,2,2,2,-2,3,5};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "db7") {
			// db7 Factorization needs to be checked. Not sure if it is correct.
		    /*ltype="dpdpdpdp";
			stages=8;
			Kconst=0.0474802;
			
			//Stage 1,2,3,4,5,6,7,8
			
			double d1[]={-0.196329};
			double p1[]={0.189042,-0.622608};
			double d2[]={0.473542,0.549384};
			double p2[]={-2.41574,-0.655465};
			double d3[]={-0.931487,0.366063};
			double p3[]={-4.03893,1.05229};
			double d4[]={-2.01533,0.247251};
			double p4[]={0.496185};
			
			lcoeff.insert(lcoeff.begin(),p4,p4+1);
			lcoeff.insert(lcoeff.begin(),d4,d4+2);
			lcoeff.insert(lcoeff.begin(),p3,p3+2);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,0,2,1,2,1,2,0,2,1,2,0,1,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));*/
			
			ltype="dpdpdpdp";
			stages=8;
			Kconst=0.299011;
			
			//Stage 1,2,3,4,5,6,7,8
			
			double d1[]={5.0935};
			double p1[]={0.0573987,-0.189042};
			double d2[]={-12.2854,5.95921};
			double p2[]={0.0291355,-0.0604279};
			double d3[]={-3.97071,1.56044};
			double p3[]={0.00330657,-0.0126914};
			double d4[]={-0.414198,0.0508159};
			double p4[]={-0.000406214};
			
			lcoeff.insert(lcoeff.begin(),p4,p4+1);
			lcoeff.insert(lcoeff.begin(),d4,d4+2);
			lcoeff.insert(lcoeff.begin(),p3,p3+2);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,1,2,-1,2,3,2,-3,2,5,2,-5,1,6};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "db8") {

		    ltype="pdpdpdpdp";
			stages=9;
			Kconst=3.54936;
			
			//Stage 1,2,3,4,5,6,7
			
			double p1[]={-5.74964};
			double d1[]={0.168817,-0.0522692};
			double p2[]={-7.40211,14.5428};
			double d2[]={0.0609093,-0.0324021};
			double p3[]={-2.7557,5.81872};
			double d3[]={0.242022,0.945295};
			double p4[]={-0.00180382,0.00018884};
			double d4[]={-0.224138,-0.952614};
			double p5[]={0.0271974,-0.246992,1.04974};
			
			lcoeff.insert(lcoeff.begin(),p5,p5+3);
			lcoeff.insert(lcoeff.begin(),d4,d4+2);
			lcoeff.insert(lcoeff.begin(),p4,p4+2);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p3,p3+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,2,2,2,-2,2,4,2,-4,2,4,2,-4,3,7};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior2.2") {

		    ltype="dp";
			stages=2;
			Kconst=0.707107;
			
			//Stage 1,2
			
			double d1[]={0.5,0.5};
			double p1[]={-0.25,-0.25};
			
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,2,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior2.4") {

		    ltype="dp";
			stages=2;
			Kconst=0.707107;
			
			//Stage 1,2
			
			double d1[]={0.5,0.5};
			double p1[]={0.046875,-0.296875,-0.296875,0.046875};
			
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+4);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,4,1};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior2.6") {

		    ltype="dp";
			stages=2;
			Kconst=0.707107;
			
			//Stage 1,2
			
			double d1[]={0.5,0.5};
			double p1[]={-.00976563,0.0761719,-0.316406,-0.316406,0.0761719,-.00976563};
			
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+6);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,6,2};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior2.8") {

		    ltype="dp";
			stages=2;
			Kconst=0.707107;
			
			//Stage 1,2
			
			double d1[]={0.5,0.5};
			double p1[]={0.00213623,-0.0204468,0.0953979,-0.327087,-0.327087,0.0953979,
							-0.0204468,0.00213623};
			
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+8);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,8,3};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior3.1") {

		    ltype="pdp";
			stages=3;
			Kconst=0.471405;
			
			//Stage 1,2,3
			
			double p1[]={0.333333};
			double d1[]={0.375,1.125};
			double p2[]={-0.444444};
			
			
			
			lcoeff.insert(lcoeff.begin(),p2,p2+1);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,-1,2,1,1,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior3.3") {

		    ltype="pdp";
			stages=3;
			Kconst=0.471405;
			
			//Stage 1,2,3
			
			double p1[]={0.333333};
			double d1[]={0.375,1.125};
			double p2[]={0.083333,-0.444444,-0.083333};
			
			
			
			lcoeff.insert(lcoeff.begin(),p2,p2+3);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,-1,2,1,3,1};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior3.5") {

		    ltype="pdp";
			stages=3;
			Kconst=0.471405;
			
			//Stage 1,2,3
			
			double p1[]={0.333333};
			double d1[]={0.375,1.125};
			double p2[]={-0.01736,0.11805,-0.444444,-0.11805,0.01736};
			
			
			
			lcoeff.insert(lcoeff.begin(),p2,p2+5);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,-1,2,1,5,2};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior3.7") {

		    ltype="pdd";
			stages=3;
			Kconst=0.471405;
			
			//Stage 1,2,3
			
			double p1[]={0.333333};
			double d1[]={0.375,1.125};
			double d2[]={0.00379,-0.032552,0.13704,0.444444,-0.13704,0.032552,-0.00379};
			
			
			
			lcoeff.insert(lcoeff.begin(),d2,d2+7);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,-1,2,1,7,4};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior3.9") {

		    ltype="pdd";
			stages=3;
			Kconst=0.471405;
			
			//Stage 1,2,3
			
			double p1[]={0.333333};
			double d1[]={0.375,1.125};
			double d2[]={-0.000854,0.008924,-0.044514,0.149007,0.444444,-0.149007,0.044514,-0.008924,0.000854};
			
			
			
			lcoeff.insert(lcoeff.begin(),d2,d2+9);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,-1,2,1,9,5};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior4.4") {

		    ltype="dpdp";
			stages=4;
			Kconst=-0.869864;
			
			//Stage 1,2,3,4
			
			double d1[]={1.58613,1.58613};
			double p1[]={-1.07964,0.0529801};
			double d2[]={0.882911,0.882911};
			double p2[]={-0.443507,-1.57612};
			
			
			
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,2,0,2,0,2,2};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "sym2") {

		    ltype="pdp";
			stages=3;
			Kconst=1.93185;
			
			//Stage 1,2,3
			double p1[]={-1.73205};
			double d1[]={0.433013,-0.0669873};
			double p2[]={1.0000};
			
			lcoeff.insert(lcoeff.begin(),p2,p2+1);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,1,1};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	}  else if (name == "sym3") {

		    ltype="dpdp";
			stages=4;
			Kconst=1.9182;
			
			//Stage 1,2,3,4
			
			double d1[]={-0.412287};
			double p1[]={0.352388,-1.56514};
			double d2[]={0.492152,0.0284591};
			double p2[]={-0.38962};
			
			lcoeff.insert(lcoeff.begin(),p2,p2+1);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,0,2,1,1,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "sym4") {

		    ltype="pdpdp";
			stages=5;
			Kconst=1.5707;
			
			//Stage 1,2,3,4,5
			
			double p1[]={0.391447};
			double d1[]={-0.339244,-0.12439};
			double p2[]={0.162031,-1.41951};
			double d2[]={0.145983,0.431283};
			double p3[]={-1.04926};
			
			lcoeff.insert(lcoeff.begin(),p3,p3+1);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,2,1,2,1,1,-1};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "sym5") {

		    ltype="dpdpdp";
			stages=6;
			Kconst=0.491434;
			
			//Stage 1,2,3,4,5,6
			
			double d1[]={0.925933};
			double p1[]={-0.498523,-0.131923};
			double d2[]={0.429326,1.45212};
			double p2[]={0.09483,-0.280402};
			double d3[]={-1.95892,-0.768066};
			double p3[]={0.17264};
			
			lcoeff.insert(lcoeff.begin(),p3,p3+1);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,0,2,1,2,1,2,0,1,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "sym6") {

		    ltype="pdpdpdp";
			stages=7;
			Kconst=-1.67071;
			
			//Stage 1,2,3,4,5,6,7
			
			double p1[]={0.226609};
			double d1[]={-0.215541,1.26707};
			double p2[]={4.25516,-0.504776};
			double d2[]={-0.23316,-0.044746};
			double p3[]={-6.62446,18.389};
			double d3[]={0.0567685,-0.144395};
			double p4[]={5.51193};
			
			lcoeff.insert(lcoeff.begin(),p4,p4+1);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p3,p3+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,2,2,2,-2,2,4,2,-4,1,5};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "sym7") {
			// db7 Factorization needs to be checked. Not sure if it is correct.
		    ltype="dpdpdpdp";
			stages=8;
			Kconst=2.14238;
			
			//Stage 1,2,3,4,5,6,7,8
			
			double d1[]={0.390551};
			double p1[]={7.18082,-0.338864};
			double d2[]={-0.137256,-0.0139115};
			double p2[]={0.13389,29.6887};
			double d3[]={-0.00010688,0.128463};
			double p3[]={2.31081,-7.4252};
			double d4[]={0.288609,0.0532701};
			double p4[]={-1.19875};
			
			lcoeff.insert(lcoeff.begin(),p4,p4+1);
			lcoeff.insert(lcoeff.begin(),d4,d4+2);
			lcoeff.insert(lcoeff.begin(),p3,p3+2);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,1,2,-1,2,3,2,-3,2,5,2,-5,1,6};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "sym8") {

		    ltype="pdpdpdpdp";
			stages=9;
			Kconst=0.101801;
			
			//Stage 1,2,3,4,5,6,7
			
			double p1[]={0.16028};
			double d1[]={-0.156265,0.710259};
			double p2[]={1.80785,-0.48815};
			double d2[]={-0.486332,1.73992};
			double p3[]={-0.256576,-0.568637};
			double d3[]={4.77714,0.864854};
			double p4[]={-0.63485,0.495251};
			double d4[]={-3.38375,-11.182};
			double p5[]={0.0185471,-0.027062,0.0894294};
			
			lcoeff.insert(lcoeff.begin(),p5,p5+3);
			lcoeff.insert(lcoeff.begin(),d4,d4+2);
			lcoeff.insert(lcoeff.begin(),p4,p4+2);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p3,p3+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,2,2,2,-2,2,4,2,-4,2,4,2,-4,3,7};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	}
	
	}

int nlifts() {
		return stages;
	}

double K() {
	return Kconst;
}

string getName() {
	return wname;
}	

void getScheme(vector<double> &coeff, vector<int> &lenvec, string &lattice,double &Kc) {
	Kc=Kconst;
	lattice=ltype;
	coeff=lcoeff;
	lenvec=plen;
	
}
	
void disp() {
	cout << "Total Number of Stages : " << stages << endl;
	cout << "--------------------------" << endl;
	int total=0;
	for (int i=0; i < stages; i++) {
			cout << "Stage : " << i+1 << endl;
			if (ltype.at(i) == 'p') {
				cout << "Primal Lifting" << endl;
			} else if (ltype.at(i) == 'd') {
				cout << "Dual Lift" << endl;
			}
			cout << "Coefficients : ";
			int t2=0;
			vector<double> poly_coeff;
			for (int j=0; j < plen[2*i];j++) {
				cout << lcoeff[total+j] << " ";
				poly_coeff.push_back(lcoeff[total+j]);
				t2++;
			}
			total=total+t2;
			cout << endl;
			cout << "Laurent Polynomial : ";
			Laurent<double> polydisp;
			polydisp.setPoly(poly_coeff,plen[2*i+1]);
			polydisp.dispPoly();
			cout << endl;
	}
	cout << "--------------------------" << endl;
	cout << " K : " << Kconst <<endl;
}	

void addLift(string &c,vector<double> &addcoeff, int mp) {
	ltype=ltype+c;
	stages=ltype.size();
	
	int len_add=addcoeff.size();
	plen.push_back(len_add);
	plen.push_back(mp);
	
	lcoeff.insert(lcoeff.end(),addcoeff.begin(),addcoeff.end());
}

void addLift(string &c,vector<double> &addcoeff, int mp, string &pos) {
	if (pos == "end") {
		ltype=ltype+c;
		stages=ltype.size();
	
		int len_add=addcoeff.size();
		plen.push_back(len_add);
		plen.push_back(mp);
	
		lcoeff.insert(lcoeff.end(),addcoeff.begin(),addcoeff.end());
	} else if (pos == "begin") {
		ltype=c+ltype;
		stages=ltype.size();
		int len_add=addcoeff.size();
		
		plen.insert(plen.begin(),mp);
		plen.insert(plen.begin(),len_add);
		
		lcoeff.insert(lcoeff.begin(),addcoeff.begin(),addcoeff.end());
		
		
		
	} else {
		cout << " The program accepts only two positions for the new lifting stage - end / begin" << endl;
	}
	
}

virtual ~liftscheme() {
	
}	
	
};

#endif // LIFT_H
