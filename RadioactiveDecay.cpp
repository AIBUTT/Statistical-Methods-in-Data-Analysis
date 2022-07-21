//		Question # 1, Assignment # 3, PHYS-699
//		Submitted by Aatif Imtiaz Butt
//		This piece of code works fine in ROOT session
#include <iostream>
#include <iomanip>
#include <cmath>
#include <TMath>
#include <TRandom>

using namespace std;

//****************************************************************************************************************
//***********************************Definitions of Functions Called In Main Program******************************
//****************************************************************************************************************

double poissonf(double *x,double *par)                                         
{
	return par[0]*TMath::Poisson(x[0],par[1]);
}

double expo_s(double *x, double *par)
{
	return par[0]*TMath::Exp(par[1]*x[0]); 
}

double Pol_1(double *x, double *par)
{
	return par[0] + ( par[1]*x[0] );
}

//****************************************************************************************************************
//***********************************Main Program Starts Here*****************************************************
//****************************************************************************************************************

int RadioactiveDecay()
{
//cout<<"To simulate Part 1 of Question, Press 1"<<endl<<"To simulate Part 2 of Question, Press 2"<<endl<<endl;
//At the moment, this piece of code solves Question#1 in one go. However, it can be easily generalized by omitting do-while loop and some if statements!!!
int QUESTIONPART=1;
TRandom3 r(0);
double x;

TFile* outputfile = new TFile("RadioactiveDecay.root", "RECREATE");
TDirectory *top_dir = gDirectory;
outputfile -> mkdir("Part_A");
outputfile -> mkdir("Part_B");
gStyle -> SetOptStat(1111111);
gStyle -> SetOptFit(0011);

do{
 
	switch (QUESTIONPART)
	{
		
//****************************************************************************************************************
//***********************************This solves Part A of Question 1*********************************************
//****************************************************************************************************************
		case 1:
		{
		
		outputfile -> cd("Part_A");
	
		int N_initial = 100, N_pre;
		double alpha = 0.01;
	
		int jj=1;	
		do {
		cout<<endl<<"Processing Question 1-A :: Case "<<jj<<endl;
		double t, p, t_total, delta_t;
		t_total = 300.0; delta_t = 1.0; p = delta_t * alpha; 
	
		if (N_initial==100 && alpha==0.01)	//For Part A, Case 1 :: Defining Histograms
		{
		TH1F *N_vs_time_PartA_Case1 = new TH1F("N_vs_time_PartA_Case1", "N(t) vs t PartA Case1", 301,0,300);
		TH1F *lnN_vs_time_PartA_Case1 = new TH1F("lnN_vs_time_PartA_Case1", "ln[N(t)] vs t PartA Case1", 301,0,300);
		}
		else if (N_initial==5000 && alpha==0.03)	//For Part A, Case 2 :: Defining Histograms
		{
		TH1F *N_vs_time_PartA_Case2 = new TH1F("N_vs_time_PartA_Case2", "N(t) vs t PartA Case2", 301,0,300);
		TH1F *lnN_vs_time_PartA_Case2 = new TH1F("lnN_vs_time PartA Case2", "ln[N(t)] vs t PartA Case2", 301,0,300);
		}
			
			//Core Programming For Part-A is this for-loop !!!
			
			N_pre = N_initial;
			for ( t = 0.0; t <= t_total; t += delta_t)
			{
				if (N_initial==100 && alpha==0.01)	//For Part A, Case 1 :: Filling Histograms
				{
					N_vs_time_PartA_Case1 -> SetBinContent(static_cast<int>(t+1), N_pre);
					if ( N_pre != 0 )	lnN_vs_time_PartA_Case1 -> SetBinContent(static_cast<int>(t+1),log ( static_cast<double>(N_pre) ));
				}
				
				else if (N_initial==5000 && alpha==0.03)	//For Part A, Case 2 :: Filling Histograms
				{
					N_vs_time_PartA_Case2 -> SetBinContent(static_cast<int>(t+1), N_pre);
					if ( N_pre != 0 )	lnN_vs_time_PartA_Case2 -> SetBinContent(static_cast<int>(t+1),log ( static_cast<double>(N_pre) ));
				}
			
				int decayed = 0;	// We are calculating the total decayed nuclie at the end of each step !!!
		
					for (int i=1; i<=N_pre; i++)
					{
					x = r.Rndm(0);
					if ( x < p ) decayed = decayed + 1;	
					}
			
				N_pre -= decayed;
	
			}	

			if (N_initial==100 && alpha==0.01)	//For Part A, Case 1 :: Writing Histograms, Fitting Them In RED And Drawing Theoratical Curve in GREEN
			{
			
		//****************************************************************************************************************
		//***********************************Cosmetics for Case-1::Part A of Question 1***********************************
		//****************************************************************************************************************			
				N_vs_time_PartA_Case1 -> Write();
				N_vs_time_PartA_Case1 -> GetXaxis() -> SetTitle("TIME (sec)");
				N_vs_time_PartA_Case1 -> GetYaxis() -> SetTitle("NUMBER OF PARENT NUCLEI");
				
				TF1 *func1 = new TF1( "fit1", expo_s, 0.0, 300.0, 2 );
				func1->SetLineWidth(1.0);
  				func1->SetLineColor(2);
				func1->SetParameter(0, 100.0);
				func1->SetParameter(1, -0.01);
				func1->SetParNames ("No", "-alpha");
				N_vs_time_PartA_Case1->Fit(fit1, "RM", "SAME");
				
				TF1 *func11 = new TF1( "fit11", expo_s, 0.0, 300.0, 2 );
				func11->SetLineWidth(1.0);
  				func11->SetLineColor(3);
				func11->SetParameter(0, 100.0);	func11->FixParameter(0, 100.0);
				func11->SetParameter(1, -0.01);	func11->FixParameter(1, -0.01);
				N_vs_time_PartA_Case1->Fit(fit11 ,"+", "SAME");
				
		//****************************************************************************************************************
		//***********************************Cosmetics for LOGRITHMIC Case-1::Part A of Question 1************************
		//****************************************************************************************************************
				
				lnN_vs_time_PartA_Case1 -> Write();
				lnN_vs_time_PartA_Case1 -> GetXaxis() -> SetTitle("TIME (sec)");
				lnN_vs_time_PartA_Case1 -> GetYaxis() -> SetTitle("ln( NUMBER OF PARENT NUCLEI )");
				
				TF1 *func2 = new TF1( "fit2", Pol_1, 0.0, 300.0, 2 );
				func2->SetLineWidth(1.0);
  				func2->SetLineColor(2);
				func2->SetParameter(0, 5.0);
				func2->SetParameter(1, -0.01);
				func2->SetParNames ("ln(No)", "-alpha");
				lnN_vs_time_PartA_Case1->Fit(fit2, "RM", "SAME");
				
				TF1 *func22 = new TF1( "fit22", Pol_1, 0.0, 300.0, 2 );
				func22->SetLineWidth(1.0);
  				func22->SetLineColor(3);
				func22->SetParameter(0, 5.0);	func22->FixParameter(0, 5.0);
				func22->SetParameter(1, -0.01);	func22->FixParameter(1, -0.01);
				func22->SetParNames ("ln(No)", "-alpha");
				lnN_vs_time_PartA_Case1->Fit(fit22, "+", "SAME");
				
			// This declaration is made just when it finishes Case 1 and starts executing Case 2
				
				N_initial = 5000;
				alpha = 0.03;
			
			}
			else if (N_initial==5000 && alpha==0.03)	//For Part A, Case 2 :: Writing Histograms, Fitting Them In RED And Drawing Theoratical Curve in GREEN 
			{
			
		//****************************************************************************************************************
		//***********************************Cosmetics for Case-2::Part A of Question 1***********************************
		//****************************************************************************************************************
			
				N_vs_time_PartA_Case2 -> Write();
				N_vs_time_PartA_Case2 -> GetXaxis() -> SetTitle("TIME (sec)");
				N_vs_time_PartA_Case2 -> GetYaxis() -> SetTitle("NUMBER OF PARENT NUCLEI");
				
				TF1 *func3 = new TF1( "fit3", expo_s, 0.0, 300.0, 2 );
				func3->SetLineWidth(1.0);
  				func3->SetLineColor(2);
				func3->SetParameter(0, 5000.0);
				func3->SetParameter(1, -0.03);
				func3->SetParNames ("No", "-alpha");
				N_vs_time_PartA_Case2->Fit(fit3, "RM", "SAME");
				
				TF1 *func33 = new TF1( "fit33", expo_s, 0.0, 300.0, 2 );
				func33->SetLineWidth(1.0);
  				func33->SetLineColor(3);
				func33->SetParameter(0, 5000.0);	func33->FixParameter(0, 5000.0);
				func33->SetParameter(1, -0.03);		func33->FixParameter(1, -0.03);
				func33->SetParNames ("No", "-alpha");
				N_vs_time_PartA_Case2->Fit(fit33, "+", "SAME");
		//****************************************************************************************************************
		//***********************************Cosmetics for LOGRITHMIC Case-2::Part A of Question 1************************
		//****************************************************************************************************************
				lnN_vs_time_PartA_Case2 -> Write();
				lnN_vs_time_PartA_Case2 -> GetXaxis() -> SetTitle("TIME (sec)");
				lnN_vs_time_PartA_Case2 -> GetYaxis() -> SetTitle("ln( NUMBER OF PARENT NUCLEI )");
				
				TF1 *func4 = new TF1( "fit4", Pol_1, 0.0, 300.0, 2 );
				func4->SetLineWidth(1.0);
  				func4->SetLineColor(2);
				func4->SetParameter(0, 8.6);
				func4->SetParameter(1, -0.03);
				func4->SetParNames ("ln(No)", "-alpha");
				lnN_vs_time_PartA_Case2->Fit(fit4, "RM", "SAME");

				TF1 *func44 = new TF1( "fit44", Pol_1, 0.0, 300.0, 2 );
				func44->SetLineWidth(1.0);
  				func44->SetLineColor(3);
				func44->SetParameter(0, 8.6);	func44->FixParameter(0, 8.6);
				func44->SetParameter(1, -0.03);	func44->FixParameter(1, -0.03);
				func44->SetParNames ("ln(No)", "-alpha");
				lnN_vs_time_PartA_Case2->Fit(fit44, "+", "SAME");

			}
		
		cout<<endl<<"Processing Question 1-A :: Case "<<jj<<" ...Comleted..."<<endl;
		jj++;	
		
		} while (jj<=2); //As we have only two cases in Question 1, Part A.
	
		break;
		}	// case-1 terminates here.


//****************************************************************************************************************
//***********************************This solves Part B of Question 1*********************************************
//****************************************************************************************************************	
		case 2:
		{
		outputfile -> cd("Part_B");
		
		//Defining Histograms for both cases in Part-B
		
		TH1F *TotalDecayed_PartB_Case1 = new TH1F("TotalDecayed_PartB_Case1", "Total Nuclei Decayed PartB Case1", 30,0,30);
		TH1F *TotalDecayed_PartB_Case2 = new TH1F("TotalDecayed_PartB_Case2", "Total Nuclei Decayed PartB Case2", 30,0,30);
		
		int N_initial = 500, N_pre;
		double alpha = 4.0E-5;
		
		int jj=1;	
		do {
		cout<<endl<<"Processing Question 1-B Case "<<jj<<endl;
		double t, p, t_total, delta_t;
		
		//Core Programming For Part-B is this for-loop !!!
		
		for (int k=1; k<=1000; k++)
		{
			int TotalDecayed=0;
			t_total = 100.0; delta_t = 10.0; p = delta_t * alpha; 	

			N_pre = N_initial;
		
				for ( t = 0.0; t <= t_total; t += delta_t)
				{
				int decayed = 0;
					for (int i=1; i<=N_pre; i++)
					{
					x = r.Rndm(0);
					if ( x < p ) decayed = decayed + 1;	
					}
				N_pre -= decayed;
				TotalDecayed += decayed;				
				}
				
			//For Part B, Case 1 :: Filling Histogram
			if (N_initial==500 && alpha==4.0E-5)		TotalDecayed_PartB_Case1 -> Fill(TotalDecayed);
			//For Part B, Case 2 :: Filling Histogram
			else if (N_initial==500 && alpha==2.0E-4)	TotalDecayed_PartB_Case2 -> Fill(TotalDecayed);
			
		}
			

			if (N_initial==500 && alpha==4.0E-5)	//For Part B, Case 1 :: Writing Histogram and Fitting It
			{
				TotalDecayed_PartB_Case1 -> Write();
				TotalDecayed_PartB_Case1 -> GetXaxis() -> SetTitle("NUMBER OF DECAYS IN 100 SECONDS");
				TotalDecayed_PartB_Case1 -> GetYaxis() -> SetTitle("NUMBER OF EXPERIMENTS");
				TF1 *pois1 = new TF1 ("pois_Case1",poissonf,0,30,2);
				pois1->SetParameter(0,1000.0);
				pois1->SetParameter(1,3.0);
				pois1->SetLineWidth(1.0);
  				pois1->SetLineColor(2);
				TotalDecayed_PartB_Case1->Fit(pois_Case1, "RM", "SAME");
			
			// This declaration is made just when it finishes Case 1 and starts executing Case 2	
				N_initial = 500;
				alpha = 2.0E-4;
			
			}
			else if (N_initial==500 && alpha==2.0E-4)	//For Part B, Case 1 :: Writing Histogram and Fitting It
			{
				TotalDecayed_PartB_Case2 -> Write();
				TotalDecayed_PartB_Case2 -> GetXaxis() -> SetTitle("NUMBER OF DECAYS IN 100 SECONDS");
				TotalDecayed_PartB_Case2 -> GetYaxis() -> SetTitle("NUMBER OF EXPERIMENTS");
				TF1 *pois2 = new TF1 ("pois_Case2",poissonf,0,30,2);
				pois2->SetParameter(0,1000.0);
				pois2->SetParameter(1,10.0);
				pois2->SetLineWidth(1.0);
  				pois2->SetLineColor(2);
				TotalDecayed_PartB_Case2->Fit(pois_Case2, "RM", "SAME");
			}
		
		cout<<endl<<"Processing Question 1-B :: Case "<<jj<<" ...Comleted..."<<endl;
		jj++;	
		
		} while (jj<3);	//As we have only two cases in Question 1, Part B.
		
		break;
		
		} // case-2 ends here.

	} //switch terminates here.

++QUESTIONPART;

}while (QUESTIONPART<=2);

cout<<endl<<"All four cases are written to RadioactiveDecay.root"<<endl;

return 0;

}
