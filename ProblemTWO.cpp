//		Question # 2, Assignment # 3, PHYS-699
//		Submitted by Aatif Imtiaz Butt
//		This piece of code works fine in ROOT session
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <TRandom>
#include "TSystem.h"
#include "TStopwatch.h"

using namespace std;

double pi=3.14159265;

//****************************************************************************************************************
//***********************************Definitions of Functions Called In Main Program******************************
//****************************************************************************************************************
double SolveForTheta_REJECTION( double _NC, double _a, double _theta )
{
return _NC*(1.0/( sin(_theta)*sin(_theta) + _a*cos(_theta)*cos(_theta) ) );
}


double SolveForTheta_INVERSION( double _NC, double _a, double _x )
{
return atan(sqrt(_a)*tan( (sqrt(_a)/_NC)*( _x + (_NC/sqrt(_a))*(-pi/2) /* atan(tan(-pi)/sqrt(_a)) */ ) ) ); //NUMERIC INCONSISTANCY ADDRESSED ANALYTICALLY
}

Double_t fitf1(Double_t *x, Double_t *par)
{
Double_t fitval1 = 0.0;
fitval1 = par[0]/(par[1]*( sin(x[0])*sin(x[0]) + 0.001*cos(x[0])*cos(x[0]) ));
return fitval1;
}

Double_t fitf2(Double_t *x, Double_t *par)
{
Double_t fitval2 = 0.0;
fitval2 = par[0]/(par[1]*( sin(x[0])*sin(x[0]) + 0.5*cos(x[0])*cos(x[0]) ));
return fitval2;
}
 
//****************************************************************************************************************
//***********************************Main Program Starts Here*****************************************************
//****************************************************************************************************************
int ProblemTWO()
{

//****************************************************************************************************************
//***********************************Declaring Histograms and Initialization of Variables*************************
//****************************************************************************************************************
	gStyle -> SetOptStat(1111111);
	gStyle -> SetOptFit(0011);
	TH1F *ThetaDistribution1_INVERSION = new TH1F("Theta Distribution INVERSION a=0.001", "INVERSION TECHNIQUE::Theta Distribution over [-pi/2 , pi/2]", 320,-1.6,1.6);
	TH1F *ThetaDistribution2_INVERSION = new TH1F("Theta Distribution INVERSION a=0.5", "INVERSION TECHNIQUE::Theta Distribution over [-pi/2 , pi/2]", 320,-1.6,1.6);
	TH1F *ThetaDistribution1_REJECTION = new TH1F("Theta Distribution REJECTION a=0.001", "REJECTION TECHNIQUE::Theta Distribution over [-pi/2 , pi/2]", 320,-1.6,1.6);
	TH1F *ThetaDistribution2_REJECTION = new TH1F("Theta Distribution REJECTION a=0.5", "REJECTION TECHNIQUE::Theta Distribution over [-pi/2 , pi/2]", 320,-1.6,1.6);
	
	TH1F *RandomNumber1_INVERSION = new TH1F("TRandom3 Statistics INVERSION a=0.001", "RANDOM NUMBER GENERATOR for INVERSION TECHNIQUE Theta_Distribution a=0.001", 100, 0, 1);
	TH1F *RandomNumber2_INVERSION = new TH1F("TRandom3 Statistics INVERSION a=0.5", "RANDOM NUMBER GENERATOR for INVERSION TECHNIQUE Theta_Distribution a=0.5", 100, 0, 1);
	TH1F *RandomNumber1_REJECTION = new TH1F("TRandom3 Statistics REJECTION a=0.001", "RANDOM NUMBER GENERATOR for REJECTION TECHNIQUE Theta_Distribution a=0.001", 100, 0, 1);
	TH1F *RandomNumber2_REJECTION = new TH1F("TRandom3 Statistics REJECTION a=0.5", "RANDOM NUMBER GENERATOR for REJECTION TECHNIQUE Theta_Distribution a=0.5", 100, 0, 1);

	double a, theta, NC;
	TRandom3 r(0);
	double x;
	TStopwatch timer;

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	FOR a=0.001	$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
	a=0.001;
	NC = 1.0/99.34588; // NC is the normalization constant for a=0.001
	
//****************************************************************************************************************
//***********************************REJECTION TECHNIQUE ALGORITHM************************************************
//****************************************************************************************************************	
	int i=1;
	timer.Start();
		do
		{
		x = r.Rndm(0);
		RandomNumber1_REJECTION -> Fill(x);
		double theta_trial =0.0;
		theta_trial = -pi/2 + ( pi/2 - (-pi/2) )*x;
		x = r.Rndm(0);
		RandomNumber1_REJECTION -> Fill(x);
			if ( SolveForTheta_REJECTION(NC,a,theta_trial) > x*NC/a )
			{
			ThetaDistribution1_REJECTION -> Fill(theta_trial);
			i++;
			}
		} while (i<=10000);
	timer.Stop();
	double CPUtime_ThetaDistribution1_REJECTION = timer.CpuTime();

//****************************************************************************************************************
//***********************************INVERSION TECHNIQUE ALGORITHM************************************************
//****************************************************************************************************************	
	
	timer.Start();
		for (int j=1; j<=10000; j++)
		{
		x = r.Rndm(0);
		RandomNumber1_INVERSION -> Fill(x);
		theta = SolveForTheta_INVERSION(NC,a,x);
		ThetaDistribution1_INVERSION -> Fill(theta);
		}
	timer.Stop();
	double CPUtime_ThetaDistribution1_INVERSION = timer.CpuTime();
	

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	FOR a=0.5	$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$		
	
	a=0.5;
	NC = 1.0/4.44288; // NC is the normalization constant for a=0.5

//****************************************************************************************************************
//***********************************REJECTION TECHNIQUE ALGORITHM************************************************
//****************************************************************************************************************	
	int i=1;
	timer.Start();
		do
		{
		x = r.Rndm(0);
		RandomNumber2_REJECTION -> Fill(x);
		double theta_trial =0.0;
		theta_trial = -pi/2 + ( pi/2 - (-pi/2) )*x;
		x = r.Rndm(0);
		RandomNumber2_REJECTION -> Fill(x);
			if ( SolveForTheta_REJECTION(NC,a,theta_trial) > x*NC/a )
			{
			ThetaDistribution2_REJECTION -> Fill(theta_trial);
			++i;
			}
		} while (i<=10000);
	timer.Stop();
	double CPUtime_ThetaDistribution2_REJECTION = timer.CpuTime();

//****************************************************************************************************************
//***********************************INVERSION TECHNIQUE ALGORITHM************************************************
//****************************************************************************************************************	
	
	timer.Start();
		for (int j=1; j<=10000; j++)
		{
		x = r.Rndm(0);
		RandomNumber2_INVERSION -> Fill(x);
		theta = SolveForTheta_INVERSION(NC,a,x);
		ThetaDistribution2_INVERSION -> Fill(theta);
		}
	timer.Stop();
	double CPUtime_ThetaDistribution2_INVERSION = timer.CpuTime();

//****************************************************************************************************************
//***********************************	CPU TIME COMSUMED	**************************************************
//****************************************************************************************************************
	cout<<"The CPU Time consumed for four cases will be saved into 'OutPutQuestion2.txt' file."<<endl;
	ofstream myOutPut;
	myOutPut.open("OutPutQuestion2.txt", ios::out);
	myOutPut<<setw(50)<<"CPU Time for INVERSION TECHNIQUE::a=0.001:"<<setw(15)<<CPUtime_ThetaDistribution1_INVERSION<<endl;
	myOutPut<<setw(50)<<"CPU Time for INVERSION TECHNIQUE::a=0.500:"<<setw(15)<<CPUtime_ThetaDistribution2_INVERSION<<endl;
	myOutPut<<setw(50)<<"CPU Time for REJECTION TECHNIQUE::a=0.001:"<<setw(15)<<CPUtime_ThetaDistribution1_REJECTION<<endl;
	myOutPut<<setw(50)<<"CPU Time for REJECTION TECHNIQUE::a=0.500:"<<setw(15)<<CPUtime_ThetaDistribution2_REJECTION<<endl;

//****************************************************************************************************************
//***********************************	Normalization of Histograms	******************************************
//*******************************And Fitting A Smooth Distribution Curve******************************************
//****************************************************************************************************************	

	TF1 *func1 = new TF1( "fit1", fitf1, -pi/2, pi/2, 2 );
	func1->SetLineWidth(1.0);
  	func1->SetLineColor(2);
	func1->SetParameters(1.0, 99.34588);
	
	TF1 *func11 = new TF1( "fit11", fitf1, -pi/2, pi/2, 2 );
	func11->SetLineWidth(1.0);
  	func11->SetLineColor(3);
	func11->SetParameters(1.0, 99.34588);
	func11->FixParameter(0, 1.0);
	func11->FixParameter(1, 99.34588);
	
	TF1 *func2 = new TF1( "fit2", fitf2, -pi/2, pi/2, 2 );
	func2->SetLineWidth(1.0);
  	func2->SetLineColor(2);
	func2->SetParameters(1.0, 4.44288);
	
	TF1 *func22 = new TF1( "fit22", fitf2, -pi/2, pi/2, 2 );
	func22->SetLineWidth(1.0);
  	func22->SetLineColor(3);
	func22->SetParameters(1.0, 4.44288);
	func22->FixParameter(0, 1.0);
	func22->FixParameter(1, 4.44288);

	for (int BN=0;BN<=320;BN++)
	{
	ThetaDistribution1_INVERSION ->SetBinContent(BN, (ThetaDistribution1_INVERSION -> GetBinContent(BN))/( (ThetaDistribution1_INVERSION ->GetBinWidth(BN))*10000) );
	}
	ThetaDistribution1_INVERSION -> Fit("fit1");
	ThetaDistribution1_INVERSION -> Fit("fit11", "+");
	
	for (int BN=0;BN<=320;BN++)
	{
	ThetaDistribution2_INVERSION ->SetBinContent(BN, (ThetaDistribution2_INVERSION -> GetBinContent(BN))/( (ThetaDistribution2_INVERSION ->GetBinWidth(BN))*10000) );
	}
	ThetaDistribution2_INVERSION -> Fit("fit2");
	ThetaDistribution2_INVERSION -> Fit("fit22", "+");
	
	for (int BN=0;BN<=320;BN++)
	{
	ThetaDistribution1_REJECTION ->SetBinContent(BN, (ThetaDistribution1_REJECTION -> GetBinContent(BN))/( (ThetaDistribution1_REJECTION ->GetBinWidth(BN))*10000) );
	}
	ThetaDistribution1_REJECTION -> Fit("fit1");
	ThetaDistribution1_REJECTION -> Fit("fit11", "+");
	
		
	for (int BN=0;BN<=320;BN++)
	{
	ThetaDistribution2_REJECTION ->SetBinContent(BN, (ThetaDistribution2_REJECTION -> GetBinContent(BN))/( (ThetaDistribution2_REJECTION ->GetBinWidth(BN))*10000) );
	}
	ThetaDistribution2_REJECTION -> Fit("fit2");
	ThetaDistribution2_REJECTION -> Fit("fit22", "+");
		
	
//****************************************************************************************************************
//***********************************	Writing Hisograms Into Root Directory	**********************************
//****************************************************************************************************************	
	
	TFile* outputfile = new TFile("ProblemTWO.root", "RECREATE");
	TDirectory *top_dir = gDirectory;
	
	ThetaDistribution1_INVERSION -> Write();
	ThetaDistribution2_INVERSION -> Write();
	ThetaDistribution1_REJECTION -> Write();
	ThetaDistribution2_REJECTION -> Write();
	
	RandomNumber1_INVERSION -> Write();
	RandomNumber2_INVERSION -> Write();
	RandomNumber1_REJECTION -> Write();
	RandomNumber2_REJECTION -> Write();
	

}
