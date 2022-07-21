//		Question # 1::Assignment # 5, PHYS-699
//		Submitted by Aatif Imtiaz Butt
//		This piece of code works fine in ROOT session

// THIS PIECE OF CODE SIMULATES 100 EXPERIMENTS AND CALCULATES ln_LIKELIHOOD_H0 AND ln_LIKELIHOOD_H1
// AND WRITES THEM TO A FILE QUESTION1.txt. IT TAKES HUGE TIME AS CURRENTLY IT IS CALCULATING 1000 VALUES FOR ln_LIKELIHOOD_H0 &
// 1000x1000 VALUES FOR ln_LIKELIHOOD_H1 FOR EACH EXPERIMENT.

// QUESTION1_DISTRIBUTION.cpp PLOTS THE DISTRIBUTION OF -2ln(LAMBDA), COMPARES WITH CHI_SQUARE DISTRIBUTION, EVALUATES REJECTION
// REGION BY PLOTTING CUMULATIVE CHI_SQUARE CURVE AND FINALLY REPORTS NUMBER OF TRIALS REJECTED.

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

double t[200], Sigma_t=0.5, Tau=1.0, weight;
double Sigma_t_H0, Estimate_Tau_H1, Sigma_t_H1;


double Solve_For_t ( double _Tau, double _x )
{
return -1.0*_Tau*log( _x*exp(-100.0/_Tau) - _x*exp(5/_Tau) + 1.0*exp(5/_Tau) );
}


double LOG_Likelihood_Function_H0 ( double _Sigma_t )
{
	double ln_f_H0 = 0.0;
	for (int i=0; i<200; i++)
	{
	ln_f_H0 += log(1.0/(2.0*Tau)) + (pow(_Sigma_t,2)/(2*pow(Tau,2))) - (t[i]/Tau) + log(TMath::Erfc( _Sigma_t/(sqrt(2.0)*Tau) - t[i]/(sqrt(2.0)*_Sigma_t) ) );
	}	
return ln_f_H0;
}


double LOG_Likelihood_Function_H1 ( double _Tau, double _Sigma_t )
{
	double ln_f_H1 = 0.0;
	for (int i=0; i<200; i++)
	{
	ln_f_H1 += log(1.0/(2.0*_Tau)) + (pow(_Sigma_t,2)/(2*pow(_Tau,2))) - (t[i]/_Tau) + log(TMath::Erfc( _Sigma_t/(sqrt(2.0)*_Tau) - t[i]/(sqrt(2.0)*_Sigma_t) ) );
	}	
return ln_f_H1;
}



//****************************************************************************************************************
//***********************************MAIN PROGRAM STARTS HERE*****************************************************
//****************************************************************************************************************

int QUESTION1()
{
	
	ofstream OUTPUT;
	OUTPUT.open("QUESTION1.txt");
//	OUTPUT<<setw(25)<<"Maximum log(L_H0)"<<setw(25)<<"Maximum log(L_H0)"<<setw(25)<<"log(LAMBDA)"<<endl;
	
	double logLAMBDA;
	for (int ExperimentNumber=0; ExperimentNumber<=99; ExperimentNumber++)
	{
 

//****************************************************************************************************************
//***********************************INVERSION TECHNIQUE BASIC ALGORITHM******************************************
//****************************************************************************************************************	

		TRandom3 r(0);
		double x_random;
	
		int i=0;
		do
		{
		x_random = r.Rndm(0);
	
		t[i] = Solve_For_t( Tau, x_random );
	
		weight = TMath::Erfc( ( Sigma_t/(sqrt(2.0)*Tau) ) - ( t[i]/(sqrt(2.0)*Sigma_t) ) );
	
		x_random = r.Rndm(0);
		
			
			if (weight > x_random*2.1) 
			{
			i++;
			}
	
		} while(i<200);
	
//****************************************************************************************************************
//***********************************EVALUATION OF LIKELIHOOD FUNCTION FOR H0 HYPOTHESIS**************************
//****************************************************************************************************************
	
		int n_H0 = 1000, index_max_H0;
		double x_H0, y_H0, y_max_H0=-1.0E30, x_at_y_max_H0;
		for (int k=1; k<=n_H0; k++)
		{
		x_H0 = 0.01*k ;
		y_H0 = LOG_Likelihood_Function_H0 ( x_H0 );
		
			if ( y_H0 > y_max_H0 ) 
			{
			y_max_H0 = y_H0;
			x_at_y_max_H0 = x_H0;
			index_max_H0 = k;
			}		
		}
		Sigma_t_H0 = x_at_y_max_H0;
		
		//cout<<"Sigma_t for Hypothesis H0 is "<<Sigma_t_H0<<endl;
		cout<<"The Maximum Value of log(Likelihood_H0) is "<<y_max_H0<<endl;
		
//****************************************************************************************************************
//***********************************EVALUATION OF LIKELIHOOD FUNCTION FOR H1 HYPOTHESIS**************************
//****************************************************************************************************************
	
		int n_H1 = 1000, index_max_H1;

		double x_H1, y_H1, z_H1,z_max_H1=-1.0E30;

	
		
		for (int k=1; k<=n_H1; k++)
		{
		x_H1 = 0.01*k ;
			for (int j=1; j<=n_H1; j++)
			{
			y_H1 = 0.01*j ;
			z_H1 = LOG_Likelihood_Function_H1 ( x_H1, y_H1 );
		
			if ( z_H1 > z_max_H1 )	z_max_H1 = z_H1;
			}			
		}
		
		cout<<"The Maximum Value of log(Likelihood_H1) is "<<z_max_H1<<endl;
		
		
//****************************************************************************************************************
//******************** WRITING QUESTION1.txt TO BE USED IN BY QUESTION1_DISTRIBUTION *****************************
//****************************************************************************************************************
		
		logLAMBDA = y_max_H0 - z_max_H1;
		
		OUTPUT.precision(15);
		OUTPUT<<"ln_LIKELIHOOD_H0["<<ExperimentNumber<<"]="<<setw(25)<<y_max_H0<<"; "<<"ln_LIKELIHOOD_H1["<<ExperimentNumber<<"]="<<setw(25)<<z_max_H1<<"; "<<"ln_LAMBDA["<<ExperimentNumber<<"]="<<setw(25)<<logLAMBDA<<"; "<<endl;
	}
return 0;
}
