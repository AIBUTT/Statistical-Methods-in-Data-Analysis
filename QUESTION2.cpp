//		Question # 2:: IT SOLVES ONLY FIRST PART OF THE QUESTION. Assignment # 4, PHYS-699
//		Submitted by Aatif Imtiaz Butt
//		This piece of code works fine in ROOT session

// SOME PORTION OF CODE IS COMMENTED OUT TO MAKE AUTOMATION WORK.
// CODE TAKES TIME WHEN IN AUTOMATION MODE. REASON IS FINE MESHING (10000) OF LOG_LIKELIHOOD FUNCTION AND AT EACH MESS POINT, 200 EVALUATIONS TAKE PLACE.

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <TRandom>
#include <TMath>

using namespace std;

//****************************************************************************************************************
//***********************************USER DEFINED FUNCTIONS AND DECLARATION OF GLOBAL VARIABLES*******************
//****************************************************************************************************************

double t[200], sigma_t=0.5, tau=1.0, sigma_tau_PLUS, sigma_tau_MINUS, weight;
double sigma_Tau, Estimate_Tau; //double sigma_Tau[100], Estimate_Tau[100]; // Uncomment this declaration and comment its counterpart when running Automation.

double Solve_For_t ( double _tau, double _x )
{
//return -1.0*_tau*log( _x*exp(-100.0/_tau) - _x + 1.0 );
return -1.0*_tau*log( _x*exp(-100.0/_tau) - _x*exp(5/_tau) + 1.0*exp(5/_tau) );
}

double pdf_fit ( double *x, double *par)
{
double a = pow(sigma_t,2)/(2*pow(tau,2)) - x[0]/tau;
double b = sigma_t/(sqrt(2.0)*tau) - x[0]/(sqrt(2.0)*sigma_t);
return ( par[0]/(2.0*tau) )*exp(a)*TMath::Erfc(b);
}

double LOG_Likelihood_Function ( double _tau )
{
	double ln_f = 0.0;
	for (int i=0; i<200; i++)
	{
	ln_f += log(1.0/(2.0*_tau)) + (pow(sigma_t,2)/(2*pow(_tau,2))) - (t[i]/_tau) + log(TMath::Erfc( sigma_t/(sqrt(2.0)*_tau) - t[i]/(sqrt(2.0)*sigma_t) ) );
	}	
return ln_f;
}



//****************************************************************************************************************
//***********************************MAIN PROGRAM STARTS HERE*****************************************************
//****************************************************************************************************************

int QUESTION2()
{
//****************************************************************************************************************
//***********************************AUTOMATION*******************************************************************
//****************************************************************************************************************
/*	TCanvas *canvas_2 = new TCanvas("canvas_2", "canvas_2", 0,0,700,500);
	canvas_2 -> cd();
	TH1F *Automate = new TH1F("Automate", "Automate", 40, 0.0, 4.0);

	int Experiment_Number;
	for(Experiment_Number=0; Experiment_Number<100; Experiment_Number++)
	{
*/
		gStyle -> SetOptStat(1111111);

		TCanvas *canvas_1 = new TCanvas("canvas_1", "canvas_1", 0,0,700,500);
		canvas_1 -> Divide(2);
	
	
		TH1F *PDF = new TH1F("PDF", "Given PDF", 150, -5.0, 10.0);
	
		
//****************************************************************************************************************
//***********************************INVERSION TECHNIQUE BASIC ALGORITHM******************************************
//****************************************************************************************************************	

		TRandom3 r(0);
		double x_random;
	
		int i=0;
		do
		{
		x_random = r.Rndm(0);
	
		t[i] = Solve_For_t( tau, x_random );
	
		weight = TMath::Erfc( ( sigma_t/(sqrt(2.0)*tau) ) - ( t[i]/(sqrt(2.0)*sigma_t) ) );
	
		x_random = r.Rndm(0);
		
			
			if (weight > x_random*2.1) 
			{

			PDF -> Fill(t[i]);
			i++;
			}
	
		} while(i<200);
	


//****************************************************************************************************************
//***********************************NORMALIZING PDF AND PLOTTING EXPECTED CURVE**********************************
//****************************************************************************************************************

		canvas_1 -> cd(1);

		PDF -> Scale(1.0/(0.1*200));	//	0.1 is the BinWidth and 200 are total events generated.
		PDF -> Draw();
		PDF -> SetTitle("Probability Distribution Function For Life Time Of Nuclie");
		PDF -> GetXaxis() -> SetTitle("Life Time Of Each Nucleus");
		PDF -> GetYaxis() -> SetTitle("PDF");
	
		TF1 *func = new TF1( "fit", pdf_fit, 0.0, 10.0, 1 );
		func->SetLineWidth(2.0);
  		func->SetLineColor(2);
		func->SetParameter(0, 1.0);	func->FixParameter(0, 1.0);
		PDF -> Fit(fit, "SAME");
	

//****************************************************************************************************************
//***********************************EVALUATION OF LIKELIHOOD FUNCTION FOR PLOTTING PURPOSE***********************
//****************************************************************************************************************
	
		canvas_1 -> cd(2);
		
		int n = 10000, index_max;
		double x[10000], y[10000], y_max=-1.0E30, x_at_y_max;
		for (int k=100; k<n; k++) // We started with k=100 to avoid ln(0) and thus skipped drawing curve near zero for likelihood function
		{
		x[k] = 0.001*k ;
		y[k] = LOG_Likelihood_Function ( x[k] );
			if ( y[k] > y_max ) 
			{
			y_max = y[k];
			x_at_y_max = x[k];
			index_max = k;
			}		
		}
	
		TGraph *log_likelihood_function1 = new TGraph (n, x, y);
		log_likelihood_function1 -> SetTitle("Log[Likelihood] vs Estimate TAU");	
		log_likelihood_function1 -> GetXaxis() -> SetTitle("Estimate TAU");
		log_likelihood_function1 -> GetYaxis() -> SetTitle("Log[L]");
		log_likelihood_function1 -> Draw("AC");
		
		Estimate_Tau = x[index_max];
		
//****************************************************************************************************************
//***********************************SIGMA_TAU+ BEING EVALUATED***************************************************
//****************************************************************************************************************
		int jj=index_max;
		do
		{
		jj++;
		}while( (y[index_max] - y[jj]) - 0.499 < 1.0E-15 );
	
		sigma_tau_PLUS = x[jj] - x[index_max];
	
//****************************************************************************************************************
//***********************************SIGMA_TAU- BEING EVALUATED***************************************************
//****************************************************************************************************************	
		int kk=index_max;
		do
		{
		kk--;
		}while( (y[index_max] - y[kk]) - 0.499 < 1.0E-15 );
	
		sigma_tau_MINUS = x[index_max] - x[kk];
	
		sigma_Tau = (sigma_tau_PLUS + sigma_tau_MINUS)/2.0;
	
		
		cout<<endl<<"The Estimate TAU is "<<Estimate_Tau<<endl;
		cout<<"The Sigma_TAU+ is "<<sigma_tau_PLUS<<endl;
		cout<<"The Sigma_TAU- is "<<sigma_tau_MINUS<<endl;
		cout<<"The Sigma_TAU  is "<<sigma_Tau<<endl;
		

//****************************************************************************************************************
//***********************************AUTOMATION OF LIKELIHOOD FUNCTION********************************************
//****************************************************************************************************************
/*	
		int n_auto = 10000, index_max_auto;
		double x_auto[10000], y_auto[10000], y_max_auto=-1.0E30, x_at_y_max_auto;
		
		for (int k=1; k<n_auto; k++)
		{
		x_auto[k] = 0.001*k ;
		y_auto[k] = LOG_Likelihood_Function ( x_auto[k] );
			if ( y_auto[k] > y_max_auto ) 
			{
			y_max_auto = y_auto[k];
			x_at_y_max_auto = x_auto[k];
			index_max_auto = k;
			}		
		}
	
		Estimate_Tau[Experiment_Number] = x_auto[index_max_auto];

//****************************************************************************************************************
//***********************************SIGMA_TAU+ BEING EVALUATED***************************************************
//****************************************************************************************************************
		int jj=index_max_auto;
		do
		{
		jj++;
		}while( (y_auto[index_max_auto] - y_auto[jj]) - 0.499 < 1.0E-15 );
	
		sigma_tau_PLUS = x_auto[jj] - x_auto[index_max_auto];
	
//****************************************************************************************************************
//***********************************SIGMA_TAU- BEING EVALUATED***************************************************
//****************************************************************************************************************	
		int kk=index_max_auto;
		do
		{
		kk--;
		}while( (y_auto[index_max_auto] - y_auto[kk]) - 0.499 < 1.0E-15 );
	
		sigma_tau_MINUS = x_auto[index_max_auto] - x_auto[kk];
	
		sigma_Tau[Experiment_Number] = (sigma_tau_PLUS + sigma_tau_MINUS)/2.0;
	
		if (Experiment_Number == 0)
		{
		cout<<"The Estimate TAU is "<<Estimate_Tau[0]<<endl;
		cout<<"The Sigma_TAU is "<<sigma_Tau[0]<<endl;
		}
		
		Automate -> Fill( (Estimate_Tau[Experiment_Number] - tau)/sigma_Tau[Experiment_Number] );	
	
	}

//****************************************************************************************************************
//***********************************FITTING GAUSSIAN CURVE*******************************************************
//****************************************************************************************************************
	Automate -> Draw();
	Automate -> Fit("gaus");
*/


return 0;
}
