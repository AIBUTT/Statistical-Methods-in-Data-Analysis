//Assignment # 1. Submitted By Aatif Imtiaz Butt.
//This Program runs fine on Linux Machine.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <new>
#include <cmath>

using namespace std;

int Lower_Triangular_Solve(double *L, double B[], double x[], int n)
{
   int i, k;

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.                                      
   
   for (k = 0; k < n; L += n, k++) {
      if (*(L + k) == 0.0) return -1;           // The matrix L is singular
      x[k] = B[k];
      for (i = 0; i < k; i++) x[k] -= x[i] * *(L + i);
      x[k] /= *(L + k);
   }

   return 0;
}
void Unit_Upper_Triangular_Solve(double *U, double B[], double x[], int n)
{
   int i, k;

//         Solve the linear equation Ux = B for x, where U is an upper
//         triangular matrix.                                      
   x[n-1] = B[n-1]; 
   for (k = n-2, U += n * (n - 2); k >= 0; U -= n, k--) {
      x[k] = B[k];
      for (i = k + 1; i < n; i++) x[k] -= x[i] * *(U + i);
   }
}
int Crout_LU_Decomposition(double *A, int n)
{
   int row, i, j, k, p;
   double *p_k, *p_row, *p_col;

//         For each row and column, k = 0, ..., n-1,
//            find the lower triangular matrix elements for column k
//            and if the matrix is non-singular (nonzero diagonal element).
//            find the upper triangular matrix elements for row k. 
 
   for (k = 0, p_k = A; k < n; p_k += n, k++) {
      for (i = k, p_row = p_k; i < n; p_row += n, i++) {
         for (p = 0, p_col = A; p < k; p_col += n, p++)
            *(p_row + k) -= *(p_row + p) * *(p_col + k);
      }  
      if ( *(p_k + k) == 0.0 ) return -1;
      for (j = k+1; j < n; j++) {
         for (p = 0, p_col = A; p < k; p_col += n,  p++)
            *(p_k + j) -= *(p_k + p) * *(p_col + j);
         *(p_k + j) /= *(p_k + k);
      }
   }
   return 0;
}
int Crout_LU_Solve(double *LU, double B[], double x[], int n)
{

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.                                      
   
   if ( Lower_Triangular_Solve(LU, B, x, n) < 0 ) return -1;

//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.
//         The diagonal part of the upper triangular part of the matrix is
//         assumed to be 1.0.

   Unit_Upper_Triangular_Solve(LU, x, x, n);
  
   return 0;
}

double current (int n)
{
	double A[n-2][n-2], RESISTOR[n][n], ROH[n][n];
	double B[n-2], V[n-2];
	double V0=1.0, Vn=0.0;
	double I=0.0;
/*
	// Opening a text file which contains a symmetric matrix containing values of resistances between various nodes.
	// We adopted this method of reading resistances from text file to give our code verstality.
	ifstream MyInput;
	MyInput.open("DataFile5.txt");
	// Reading from file the values of resistors and making a ROH matrix.
	while ( !MyInput.eof() )
	{
		for ( int i = 0; i<n; i++ )
			for ( int j = 0; j<n; j++ )
			{
				MyInput>>RESISTOR[i][j];
				if ( i == j )
					ROH[i][j] = 0.0;
				else if ( i != j && RESISTOR[i][j] != 0.0 )
					ROH[i][j]=1.0/RESISTOR[i][j];
				else
				{cout<<"One of the Resistors is 0 OHM, and ROH is infinity"<<endl; return 0;}
			}
	}
*/
//**********************************************************************************************************************************
	//Printing Resistor Matirx to the terminal for varification.
	//The IF-ELSE structure is determining the RESISTOR matrix and ROH matrix based on formula given in Assignment.
	//If the code is made to read resistor values from text file, we need to switch off this IF-ELSE structure.

//	cout<<endl<<"The values read (in OHM) for RESISTOR Matrix are: "<<endl<<endl;
	for ( int i = 0; i<n; i++ )
	{
		for ( int j = 0; j<n; j++ )
		{
			if (i==j)
			{
				RESISTOR[i][j]=0.0;	ROH[i][j] = 0.0;
			}
			else
			{
				RESISTOR[i][j]=fabs( (double)i -(double)j );	ROH[i][j]=1.0/RESISTOR[i][j];
			}

	//		cout<<setw(15)<<RESISTOR[i][j];
		}
	//	cout<<endl;
	}
//**********************************************************************************************************************************	
	//Printing ROH Matirx to the terminal for varification.
/*
	cout<<endl<<"The values calculated (in Inverse OHM) for ROH Matrix are: "<<endl<<endl;
	for ( int i = 0; i<n; i++ )
	{
		for ( int j = 0; j<n; j++ )
		{
			cout<<setw(15)<<ROH[i][j];
		}
		cout<<endl;
	}
*/
//**********************************************************************************************************************************
	//Making matrix A which will serve as an input to Crout LU Decomposition Function. 
	for ( int i = 0; i<n-2; i++ )
		for ( int j = 0; j<n-2; j++ )
		{
			if ( i == j )
			{ 
				A[i][j] = 0.0;
				for (int k = 0; k<n; k++)
					A[i][j] += -1*ROH[i+1][k];
			}
			else if ( i != j )
				A[i][j]=ROH[i+1][j+1];
		}
//**********************************************************************************************************************************		
/*
cout<<"The Matrix A before Crout's LU Decomposition is"<<endl<<endl;
		for ( int i = 0; i<n-2; i++ )
	{
		for ( int j = 0; j<n-2; j++ )
		{
			cout<<setw(15)<<A[i][j];
		}
		cout<<endl;
	}
*/
//**********************************************************************************************************************************	
	//Making matrix B which will serve as an input to Crout LU Decomposition Function. 
	for ( int i = 0; i<n-2; i++ )
	{
		B[i]=-1*ROH[0][i+1]*V0;
	}
//**********************************************************************************************************************************
/*
cout<<"The Matrix B before Crout's Decomposition is"<<endl<<endl;
	for ( int i = 0; i<n-2; i++ )
	{
		cout<<B[i]<<endl;
	}
*/	
//**********************************************************************************************************************************
	//Using Crout LU Decomposition Function to make Lower Triangular Matrix and Unit Upper Triangular Matrix of Matrix A.
	//It will check if Matrix A is singular. If not, it will proceed further to use Crout_LU_Solve Function and yield solution
	//to the system of linear equations.
int err; 
err = Crout_LU_Decomposition(&A[0][0], n-2);
if (err < 0) {printf(" Matrix A is singular\n"); return 0;}
	else err = Crout_LU_Solve(&A[0][0], B, V, n-2);
	
//**********************************************************************************************************************************
/*
	//We are printing LU Decomposed Matrix A.
	cout<<endl<<"The Matrix A after Crout's LU Decomposition is"<<endl<<endl;
	for ( int i = 0; i<n-2; i++ )
	{
		for ( int j = 0; j<n-2; j++ )
		{
			cout<<setw(15)<<A[i][j];
		}
		cout<<endl;
	}
*/	
//***********************************************************************************************************************************
/*
	//We are printing solution set to the terminal. It is the value of potential difference at each node.
	cout<<endl<<"The Voltage solution (in Volts) to Krichoff's Circuit with "<<n<<" nodes is"<<endl<<endl;
	cout<<setw(15)<<V0<<endl;
	for ( int i = 0; i<n-2; i++ )
	{
		cout<<setw(15)<<V[i]<<endl;
	}
	cout<<setw(15)<<0.0<<endl;
*/	
//***********************************************************************************************************************************
	//We calculate here the current drawn using formula given in Assignment#1.
	I = I + (V0 - Vn)*ROH[0][n-1];
	for ( int i = 0; i<n-2; i++ )
	{ I = I + (V0 - V[i])*ROH[0][i+1];}
//	cout<<endl<<"The CURRENT in Krichoff's Circuit with "<<n<<" nodes is "<<I<<" Amperes."<<endl<<endl;
	
	return I;
}
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
int main()
{
	cout<<endl<<endl<<setw(100)<<"Assignment # 1. Submitted By Aatif Imtiaz Butt"<<endl<<endl;
	
	/*
	int n;
	cout<<"Please Enter Total Number of Nodes in Krichoof's Circuit"<<endl;
	cin>>n;
	*/
	
	cout<<setw(25)<<"Number of Nodes"<<setw(35)<<"Total Current Drawn (Amperes)"<<endl;
	for (int T=3; T<=20; T++)
	{
	cout<<setw(25)<<T<<setw(35)<<current(T)<<endl;
	}
	return 0;
}