//		Question # 3, Assignment # 3, PHYS-699
//		Submitted by Aatif Imtiaz Butt
//		This piece of code works fine in ROOT session

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <TRandom>

using namespace std;

double pi=3.14159265, m=0.510999E+6, alpha=1.0/137.035999, hBAR=1.05457E-34, m1=9.1093821545E-31, c=2.99792458E+8;

//****************************************************************************************************************
//***********************************Definitions of Functions Called In Main Program******************************
//****************************************************************************************************************
void bad(void)
{
	cerr << "The Maximum Weight Calculated during the exicution is greater than 2.01" << endl;
	cerr << "There is a need to debug the algorithm !!!" << endl;
	exit(1);
}

double SolveForU ( double _k, double _x )
{
return ( (m/_k)*( pow(( 1.0 + 2.0*( _k/ m  ) ),_x) - 1.0 ) ) ;
}

double k_prime_BY_k ( double _k, double _u)
{
return 1.0/( 1.0 + (( _k/m )*_u) );
}

double Klein_Nishina(double *x, double *par)	// For Scaling... Its not a fit !!!
{
double fitval1, fitval2, fitval3, fitval;

fitval1 = (alpha*alpha)/(2.0*m*m);
fitval2 = 1.0/( 1.0 + ( (par[0]/m)*(1-cos(x[0])) ) );
fitval3 = sin(x[0]);

fitval = par[1]*fitval1 * ( pow(fitval2,3) + fitval2 - (pow(fitval2,2)*pow(fitval3,2)) ) * fitval3; 

return fitval;
}

//****************************************************************************************************************
//***********************************Main Program Starts Here*****************************************************
//****************************************************************************************************************

int ComptonScattering()
{

//****************************************************************************************************************
//***********************************Declaring Histograms*********************************************************
//****************************************************************************************************************
	
	gStyle -> SetOptStat(1111111);
	TH1F *Theta_Distribution = new TH1F("Theta Distribution", "Theta Distribution In Radians", 320,0,pi);
	TH1F *Theta_Distribution_Degree = new TH1F("Theta Distribution in Degree", "Theta Distribution In Degree", 360,0,180);
	TH1F *U_Distribution = new TH1F("U Distribution", "U Distribution", 200,0,2);
	TH1F *WEIGHT_Distribution = new TH1F("Weight Distribution", "Weight Distribution", 200, 0, 2);

//****************************************************************************************************************
//***********************************Initialization of Variables**************************************************
//****************************************************************************************************************

double k, u, theta, weight, w_max=-1.0E-35, w_min=1.0E+35, sum=0.0;
int unit, TotalEvents;
cout<<endl<<"Please enter total number of events to be generated"<<endl;
cin>>TotalEvents;
cout<<"Please choose from options below the appropriate units for Photon Energy"<<endl;
cout<<setw(20)<<"Press 1 if KeV"<<endl<<setw(20)<<"Press 2 if MeV"<<endl<<setw(20)<<"Press 3 if GeV"<<endl<<setw(20)<<"Press 4 if TeV"<<endl;
cout<<"As an example, if you want to produce Compton Scattered Events for 5KeV Photon, Choose 1 from above set of options and then press 5"<<endl;
cin>>unit;
cout<<"Enter Initial Energy of Incident Photon"<<endl;
cin>>k;

//****************************************************************************************************************
//***********************************Conversion of Energy in electronVolt*****************************************
//****************************************************************************************************************


switch(unit)
{
	case 1:
	{k = k*1.0E+3; break;}
	case 2:
	{k = k*1.0E+6; break;}
	case 3:
	{k = k*1.0E+9; break;}
	case 4:
	{k = k*1.0E+12; break;}
}
cout.precision(3); cout.setf(std::ios::scientific);
cout<<"Initial Energy of Incident Photon is "<< k <<"eV"<<endl<<endl;


//****************************************************************************************************************
//***********************************Core Calculations & Filling Histograms***************************************
//****************************************************************************************************************

	TRandom3 r(0);
	double x;

		int j=0;
		do
		{
			x = r.Rndm(0);
			u = SolveForU(k,x);
			theta = acos(1.0 - u);
			
			weight = pow( k_prime_BY_k(k,u) , 2 ) + ( k_prime_BY_k(k,u) * pow(u,2) ) - ( 2.0 * u * k_prime_BY_k(k,u) ) + 1.0;
			
			if ( weight < w_min ) w_min = weight;
			if ( weight > w_max ) w_max = weight;
			
			x = r.Rndm(0);
			
			if (weight > x*2.0) 
			{
			Theta_Distribution -> Fill(theta);
			Theta_Distribution_Degree -> Fill(theta*360/(2*pi));
			U_Distribution -> Fill(u);
			WEIGHT_Distribution -> Fill(weight);
			
			sum+=pow(weight,2);	// Adding square of weights to be used in calculating uncertainty
			
			j++;
			}
			
		} while (j<TotalEvents);
		
	cout<<"Maximum Value of Weight is "<<w_max<<endl<<"Minimum Value of Weight is "<<w_min<<endl;
	if (w_max >= 2.01) bad();
	
//****************************************************************************************************************
//***********************************SCALING (NOT FITTING) Klein- Nishina Formula*********************************
//****************************************************************************************************************	
/*
	for (int BN=1;BN<=320;BN++)
	{
	Theta_Distribution ->SetBinContent(BN, (Theta_Distribution -> GetBinContent(BN))/( (Theta_Distribution ->GetBinWidth(BN))*10000) );
	}
*/		
	TF1 *func1 = new TF1( "fit1", Klein_Nishina, 0.0, pi, 2 );
	func1->SetLineWidth(2.0);		
  	func1->SetLineColor(2);			//*************************************************************************************************************************//
	func1->SetParameters(k,1.0E20);		//******This is NOT a fit. We use Parameter[1] to rescale Klein-Nishina Formula with our THETA distribution in radians*****//
	func1->FixParameter(0, k);		//*************************************************************************************************************************//
	Theta_Distribution -> Fit("fit1"); 

//****************************************************************************************************************
//***********************************Saving Histograms in ROOT File***********************************************
//****************************************************************************************************************

	TFile* outputfile = new TFile("ComptonScattering.root", "RECREATE");
	TDirectory *top_dir = gDirectory;

	Theta_Distribution -> Write();
	Theta_Distribution_Degree -> Write();
	U_Distribution -> Write();
	WEIGHT_Distribution -> Write();

//****************************************************************************************************************
//***********************************Calculating Total Cross Section**********************************************
//****************************************************************************************************************
	
	double kk = k*1.602E-19; //converting photon energy from eV to Joules
	double MeanWeight = WEIGHT_Distribution ->GetMean();
	double TotalXSection = MeanWeight * ( 2*pi*pow(alpha,2)*pow(hBAR,2)/(2*kk*m1) )*( log(2*pow(m1,2)*pow(c,2) + 4*m1*kk) - log(2*pow(m1,2)*pow(c,2)) );
	cout<<"Total Cross-Section is "<< TotalXSection <<" Meter-squared"<<endl;

//****************************************************************************************************************
//************BenchMarking Total Cross-Section With WOLFRAM's Demonstration Project::First Three Cases Only*******
//******************http://demonstrations.wolfram.com/KleinNishinaFormulaForPhotonElectronScattering/*************
//****************************************************************************************************************
	
	int WOLFRAM = 0.0;
	if (k==5000.0) WOLFRAM = 8.218;
	else if (k==2.0E+6) WOLFRAM = 1.843;
	else if (k==1.0E+9) WOLFRAM = 0.014;
	cout<<"Total Cross-Section Calculated by Wolfram Demonstration Project is "<< WOLFRAM * 2.818E-15 * 2.818E-15  <<" Meter-squared"<<endl;

//****************************************************************************************************************
//***********************************Calculating Uncertainty In Total Cross Section*******************************
//****************************************************************************************************************	

	double uncertainty = (sum/TotalEvents)- pow(MeanWeight,2);	// Variance Formula On Slide 24 of Lecture#4.
	cout<<"UNCERTAINTY IS "<<uncertainty<<endl;

}
