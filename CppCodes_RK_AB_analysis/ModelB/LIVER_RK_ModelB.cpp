/*
Code to perform data analysis on the Liver data set in Davies, Coolen and Galla (2021) for retarded kernel model B.
Written by A. Davies (2021).
20 sets of training and test data (obtained from the R code Liver_Split_Average.R) are read in.
Retarded Kernel Model B is fitted to the training data using maximum likelihood.
Minimisation of the negative log likelihood is performed using Powell's method via the Numerical Recipes 3 (NR3) function powell.minimize()
The NR3 source code can be dowloaded here http://numerical.recipes

First we perform the fixed base time analysis with base time t = 3 years:
The prediction time u is varied in steps of 0.2 years from 3 to 10 years
Prediction error calculated using the fitted model and test data for each combination of t and u is stored at each iteration.

Then we perform the fixed prediction window analysis for three windows:
w1=1 year, w2=2 years, w3=3 years
At w1 we use base times t=0-9 years in steps of 0.2 years
At w2 we use base times t=0-8 years in steps of 0.2 years
At w3 we use base times t=0-7 years in steps of 0.2 years
Again, prediction error for each scenario is stored for each iteration.
An average of PE over the 20 iterations is performed subsequently by C++ code Average_Split.cpp

The model fitted here specifies a fixed assocation parameter for individuals with s=0 (n=1)
For the decaying parameter model, the relevant changes to the code are commented out and labelled 's=0 Decaying Association Model'
The changes correpsond to one line in three functions: mu_sum1, mus_sum2 and EXP_B 

Description of retarded kernel model B can be found in the manuscriupt 'Retarded kernels for longitudinal survival analysis and dynamic predictions' Davies, Coolen and Galla (2021). Detailed derivations of the functions evaluated by this code can be found in the corresponding Supplementary Material.

This code took ~ 45 mins to run
*/
#include <iostream>
#include <iomanip>
#include<random>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string> 
#include <algorithm>
#include <cstdlib>
#include<ctime>
#include<vector>
#include "nr3.h"
#include "mins.h"
#include "mins_ndim.h"
#include "time.h"

using namespace std;


//Define global variables----------------------------------------------------------------------------------------

const int p {1}; //number of longitudinal covariates
const int q {1}; //number of fixed covariates

//variables and vectors for training data:
const int N_train {244}; 	//no. of individuals in the training data set
VecInt n(N_train);		//vector containing the no. of observations of each individual in the training data
vector<double> X(N_train);	//vector containing event times in the training data
VecInt Censor(N_train);		//vector of the censoring indicator for the training data 
std::vector<double> T;		//vector of 'switch times' for the training data
std::vector<double> Z;		//vector of longitudinal covariate observations in the training data
std::vector<double> t_obs;	//vector of observation times of individuals in the training data
VecDoub Z_fix(N_train*q);	//vector of fixed covariates in the training data

//variables and vectors for test data:
const int N_test {244};		//no. of individuals in the test data set
VecInt n_test(N_test);		//vector containing the no. of observations of each individual in the test data
vector<double> X_test(N_test);	//vector containing event times in the test data
VecInt Censor_test(N_test);	//vector of the censoring indicator for the test data
std::vector<double> T_test;	//vector of 'switch times' for the test data
std::vector<double> Z_test;	//vector of longitudinal covariate observations in the test data
std::vector<double> t_obs_test;	//vector of observation times of individuals in the test data
VecDoub Z_fix_test(N_test*q);	//vector of fixed covariates in the test data

//---------------------------------------------------------------------------------------------------------------


//Maximum likelihood functions-----------------------------------------------------------------------------------

/*
To perform maximum likelihood inference of retarded kernel model B we minimise the negative log likelihood for the training data.
The function func() outputs the negative log likelihood for model B calling to functions mu_sum1(), j_sum() and mu_sum2().
The equation that correpsonds to this function is Eq (3.22) in Davies, Coolen and Galla (2021) 
We make use of Eq (S16) in the correpsonding supplementary material which gives the result of the integration in (3.22) when using the nearest neighbour interpolation procedure for Model B.
*/

Doub mu_sum1(VecDoub V, vector<double> T, vector<double> Z, VecDoub Z_fix, vector<double> X, int i, int j, int p) {
	double mu_sum = 0.0;
	
	//work out the sum over n[] up to j-1 for use in Z and T indicators
	int sum_nj = 0;
	for(int k=0; k<j; k++){
		sum_nj = sum_nj + n[k];
	}

	double s_j = T[sum_nj + j + n[j]]; //the final observation time = final switch time (when a=n[j])
	double min_XS = min(s_j, X[i]);

	//sum over p longitudinal covariates
	for (int mu = 0; mu < p; mu++) {
		int tau = p + mu;
		//sum over n[j] longitudinal observations
		for(int a = 0; a < n[j]; a++){	
			if(s_j==0){
				mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu];				//s=0 Fixed Association Model
				//mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu]*exp(-X[i]/V[tau]);	//s=0 Decaying Association Model
			}
			else if(X[i]<=T[sum_nj + j + a]){
				mu_sum = mu_sum;
			}
			else{
				double min_XT = min(X[i], T[sum_nj + j + a + 1]);
				
				//first term		
				mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu]*(min_XT - T[sum_nj + j + a])/min_XS;
								
				//second term
				//term 2a:				
				double diffa = X[i] - T[sum_nj+j+a];
				double diffa1 = X[i] - min_XT;
				double term2a = exp(-diffa1/V[tau]) - exp(-diffa/V[tau]);
				
				//term 2b:
				double diffs = X[i] - min_XS; 
				double term2b = (min_XT-T[sum_nj+j+a])*(exp(-diffs/V[tau])-exp(-X[i]/V[tau]))/min_XS;
				
				mu_sum = mu_sum + V[mu]*Z[sum_nj*p + a*p + mu]*(term2a - term2b);
			}
		} //end sum over a
	} //end sum over mu
	//sum over q fixed covariates
	for (int nu=0; nu<q; nu++){
		mu_sum = mu_sum + V[2*p+nu]*Z_fix[j*q+nu];
	}
	return mu_sum;
}

Doub j_sum(VecDoub V, vector<double> T, vector<double> Z, VecDoub Z_fix, vector<double> X, int i, int N_train, int p) {
	double j_sum = 0.0;
	//sum over j (individuals in the training data)
	for (int j = 0; j < N_train; j++) {
		//work out the sum over n[] up to j-1 for use in Z and T indicators
		int sum_nj = 0;
		for(int k=0; k<j; k++){
			sum_nj = sum_nj + n[k];
		}
		double s_j = T[sum_nj + j + n[j]]; //the final observation time = final switch time (when a=n[j])
		
		if (X[j] - X[i] < 0.0 || X[i]< 0.0) {
			j_sum = j_sum;
			
		}
		else {
			double e = exp(mu_sum1(V, T, Z, Z_fix, X, i, j, p));
			j_sum = j_sum + e;
		}
	}
	return j_sum;
}

Doub mu_sum2(VecDoub V, vector<double> T, vector<double> Z, VecDoub Z_fix, vector<double> X, int i, int p) {
	double mu_sum = 0.0;

	//work out the sum over n[] up to i-1 for use in Z and T indicators
	int sum_ni = 0;
	for(int k=0; k<i; k++){
		sum_ni = sum_ni + n[k];
	}

	double s_i = T[sum_ni + i + n[i]]; //the final observation time = final switch time (when a=n[i])
	double min_XS = min(s_i,X[i]); //should always be s_i

	//sum over p longitudinal covariates
	for (int mu = 0; mu < p; mu++) {
		int tau = p + mu;
		//sum over n[i] longitudinal observations
		for(int a = 0; a < n[i]; a++){	
			if(s_i==0){
				mu_sum = mu_sum + V[mu]*Z[sum_ni*p + a*p + mu];				//s=0 Fixed Association Model
				//mu_sum = mu_sum + V[mu]*Z[sum_ni*p + a*p + mu]*exp(-X[i]/V[tau]);	//s=0 Decaying Association Model
			}
			else if(X[i]<=T[sum_ni + i + a]){ //this shouldn't ever be the case
				cout << "oh no mu sum2" << endl;
				mu_sum = mu_sum;
			}
			else{
				double min_XT = min(X[i], T[sum_ni + i + a + 1]); //should always be T
				
				//first term		
				mu_sum = mu_sum + V[mu]*Z[sum_ni*p + a*p + mu]*(min_XT - T[sum_ni + i + a])/min_XS;
								
				//second term
				//term 2a:			
				double diffa = X[i] - T[sum_ni+i+a];
				double diffa1 = X[i] - min_XT;
				double term2a = exp(-diffa1/V[tau]) - exp(-diffa/V[tau]);
				
				//term 2b:
				double diffs = X[i] - min_XS; 
				double term2b = (min_XT-T[sum_ni+i+a])*(exp(-diffs/V[tau])-exp(-X[i]/V[tau]))/min_XS;
				
				mu_sum = mu_sum + V[mu]*Z[sum_ni*p + a*p + mu]*(term2a - term2b);
			}
		} //end sum over a
	} //end sum over mu
	//sum over q fixed covariates
	for (int nu=0; nu<q; nu++){
		mu_sum = mu_sum + V[2*p+nu]*Z_fix[i*q+nu];
	}
	return mu_sum;
}

Doub func(VecDoub_I &V) {

	Doub function = 0.0;
	//penalise for negative tau (NB: for Liver data V[1] represents tau)
	if(V[1]<0){
		function = function + 10e15;
	}
	else{
		//sum over i (individuals in training data)
		for (int i = 0; i < N_train; i++) {
			if(Censor[i]==0){
				function = function; 
			}
			else{
				function = function + log(j_sum(V, T, Z, Z_fix, X, i, N_train, p)) - mu_sum2(V, T, Z, Z_fix, X, i, p);
			}
		} //end sum over i
	}
	return function;
}

//Prediction error functions-------------------------------------------------------------------------------------
/*
Functions to calculate the prediction error for retarded kernel model B fitted to training data evaluated on test data.
Function PE_SQ() returns the prediction error for given base time, prediction time, and training and test data vectors.
PE_SQ() calls to the function S_t() which calculates survival probability.
Function S_t() calls to function EXP_B().
Within function PE_SQ() 'new data' vectors are created from the test data vectors restricting observations to be <= t (base time).
The Eq for prediction error is given in Eq (4.26) of Davies, Coolen and Galla (2021). 
The Eq for survival probability (of the RK model) for a given base time and prediction time is given in Eq (4.27). 
For the integrals we again we make use of Eq (S16) in the correpsonding Supplementary Material.
*/

double EXP_B(double ti, VecDoub kappa, VecDoub tau, VecDoub gamma, vector<double> Z_t, vector<double> T_t, VecDoub Zf_t, VecInt n_t, int j){
	double sum = 0.0;

	//work out the sum over n[] up to j-1 for use in Z and T indicators
	int sum_nj = 0;
	for(int k=0; k<j; k++){
		sum_nj = sum_nj + n_t[k];
	}

	double s_j = T_t[sum_nj + j + n_t[j]]; //the final obs time (or final obs time before basetime) = final switch time (when a=n[k])		
	double min_XS = min(s_j,ti);

	//sum over p longitudinal covariates
	for (int mu = 0; mu < p; mu++) {
		//sum over n[j] longitudinal observations
		for(int a = 0; a < n_t[j]; a++){	
			if(s_j==0){
				sum = sum + kappa[mu]*Z_t[sum_nj*p + a*p + mu];				//s=0 Fixed Association Model
				//sum = sum + kappa[mu]*Z_t[sum_nj*p + a*p + mu]*exp(-ti/tau[mu]);	//s=0 Decaying Association Model
			}
			else if(ti<=T_t[sum_nj + j + a]){
				sum = sum;
			}
			else{
				double min_XT = min(ti, T_t[sum_nj + j + a + 1]);
				
				//first term		
				sum = sum + kappa[mu]*Z_t[sum_nj*p + a*p + mu]*(min_XT - T_t[sum_nj + j + a])/min_XS;
				
				//second term
				//term 2a:			
				double diffa = ti - T_t[sum_nj+j+a];
				double diffa1 = ti - min_XT;
				double term2a = exp(-diffa1/tau[mu]) - exp(-diffa/tau[mu]);
				
				//term 2b:
				double diffs = ti - min_XS; 
				double term2b = (min_XT-T_t[sum_nj+j+a])*(exp(-diffs/tau[mu])-exp(-ti/tau[mu]))/min_XS;
				
				sum = sum + kappa[mu]*Z_t[sum_nj*p + a*p + mu]*(term2a - term2b);
			}
		} //end of sum over a
	} //end of sum over mu
	//sum over q fixed covariates
	for(int nu = 0; nu<q; nu++){
		sum = sum + gamma[nu]*Zf_t[j*q+nu];
	}
	double e = exp(sum); 
	return e;

}


double S_t(VecDoub kappa, VecDoub tau, VecDoub gamma, double predTime, double baseTime, vector<double> X, VecDoub Z_fix, VecDoub Z_fix_test, vector<double> Z, vector<double> Z_new, vector<double> T, vector<double> T_new, VecInt n, VecInt n_new, VecInt Censor, int j, int N_train){ 

	//j refers to the indiviudal in the test data
	//i and k label individuals in the training data

	//sum over individuals (event times) in the training data
	double sum_i=0.0;
	for(int i=0; i<N_train; i++){ 
		if(Censor[i]==0 || X[i] < baseTime || X[i] > predTime){ 
			sum_i = sum_i;  
		}
		else{
			//base hazard - sum over individuals in the training data:
			double denom = 0.0;
			for(int k=0; k<N_train; k++){
				int sum_nk = 0;
				for(int l=0; l<k; l++){
					sum_nk = sum_nk + n[l];
				}
				double s_k = t_obs[sum_nk + n[k] - 1]; //the final observation time = final switch time (when a=n[k]-1)

				if(X[i] < 0.0 || X[i] > X[k]){ denom = denom; }
				else{
					denom = denom + EXP_B(X[i], kappa, tau, gamma, Z, T, Z_fix, n, k);
				}
			}
			//for individual j (test) we only have data up to baseTime hence 'new' vectors (from test data)
			double numer = EXP_B(X[i], kappa, tau, gamma, Z_new, T_new, Z_fix_test, n_new, j);
			sum_i = sum_i + numer/denom;
		}
	}
	double S = exp(-sum_i);
	return S;
}


double PE_SQ(double predTime, double baseTime, VecDoub kappa, VecDoub tau, VecDoub gamma, vector<double> X, vector<double> X_test, VecInt Censor, VecInt Censor_test, vector<double> Z, vector<double> Z_test, VecDoub Z_fix, VecDoub Z_fix_test, vector<double> t_obs, vector<double> t_obs_test, VecInt n, VecInt n_test, int N_train, int N_test){
	
	//create Z_new, n_new, t_obs_new and T_new from baseTime and TEST vectors
	//keep only observations >= baseTime

	//No. of observations
	VecInt n_new(N_test);
	int sum_new = 0;
	for(int j=0; j<N_test; j++){ 
		int sum_nj = 0;
		for(int i=0; i<j; i++){
			sum_nj = sum_nj + n_test[i];
		} 	
		
		int count_n = 0;
		for(int a=0; a<n_test[j]; a++){	
			if(t_obs_test[sum_nj+a]<= baseTime){
				count_n = count_n + 1;
			}
			else{ count_n = count_n; }
		}
		n_new[j] = count_n; //new n defined 
		sum_new = sum_new + n_new[j];
	}

	vector<double> Z_new(p*sum_new);	//longitudinal covariates
	vector<double> t_obs_new(sum_new);	//observation times
	vector<double> T_new(N_test+sum_new);	//switch times
	
	for(int j=0; j<N_test; j++){ 
	
		int new_nj = 0;
		for(int i=0; i<j; i++){
			new_nj = new_nj + n_new[i];
		} 
		int sum_nj = 0;
		for(int i=0; i<j; i++){
			sum_nj = sum_nj + n_test[i];
		} 

		//t_obs new************************************
		for(int a=0; a<n_new[j]; a++){
			t_obs_new[new_nj+a] = t_obs_test[sum_nj+a]; 	
		}

		//Z_new****************************************
		for(int mu=0; mu<p; mu++){
			for(int a=0; a<n_new[j]; a++){
				Z_new[p*new_nj + p*a + mu] = Z_test[p*sum_nj + p*a + mu];
			}
		}

		//T_new****************************************
		//define vector of observation times for individual j (length n_new[j])
		VecDoub t_o(n_new[j]);		
		for(int i=0; i<n_new[j]; i++){
			t_o[i] = t_obs_test[sum_nj+i]; 
		}
		
		T_new[new_nj + j] = 0.0; 
		T_new[new_nj + j + n_new[j]] = t_o[n_new[j]-1];
		for(int b=1; b<n_new[j]; b++){
			T_new[new_nj+j+b] = (t_o[b-1] + t_o[b])/2.0;
		}
	}


	double r = 0.0; //no. of subjects in test data at risk at baseTime
	double sum = 0.0;
	//sum over individuals in test data still at risk at baseTime (Xj>=baseTime)
	for(int j=0; j<N_test; j++){
				
		if(X_test[j] < baseTime){ 
			sum = sum; 
		}
		else{
			r = r + 1.0;

			//Prob of j's survival to predTime given covariates up to baseTime (new vectors) and survival to baseTime
			double SprobTau = S_t(kappa, tau, gamma, predTime, baseTime, X, Z_fix, Z_fix_test, Z, Z_new, T, T_new, n, n_new, Censor, j, N_train);
			//term 1: still alive at predTime
			if(X_test[j] >= predTime){
				sum = sum + pow((1.0-SprobTau),2.0);
			}
			//term 2: experienced an event by predTime (and after baseTime)
			else if(Censor_test[j]==1){
				sum = sum + pow(SprobTau,2.0);
			}
			//term 3: censored by predTime (and after baseTime)
			else{
				//Prob of j's survival to predTime given covariates up to baseTime (new vectors) and survival to X_test[j]
				double SprobTj = S_t(kappa, tau, gamma, predTime, X_test[j], X, Z_fix, Z_fix_test, Z, Z_new, T, T_new, n, n_new, Censor, j, N_train);
				sum = sum + pow((1.0-SprobTau),2.0)*SprobTj + pow(SprobTau,2.0)*(1.0-SprobTj);
			}
		}
	}	
	double M_Y_sq = sum/r;
	return M_Y_sq;

}



int main() {

	time_t start, end; 
	time(&start);

	//base folder where training and test data are stored for Liver data
	string fn_base = "~.../LIVER_DATA/";

	//iterate over 20 splits of training/test data
	int rep = 20;
	for(int loop=0; loop<rep; loop++){

		//Folder label for iteration 'loop'
		cout << "**************loop = " << loop+1 << endl;
		string label = "LOOP"+to_string(loop+1);
	
		//---------------------------------------------------------------------------------------------------------------------------
		//READ IN TRAIN DATA LIVER---------------------------------------------------------------------------------------------------
		//---------------------------------------------------------------------------------------------------------------------------
		
		int ev = 0;
		double value1;
		ifstream file;
		file.open(fn_base+label+"/TRAIN/Events.txt");
		while (file >> value1)
		{
			X[ev] = value1;   
			ev++;
		}
		file.close();		

		int no = 0;
		int value2;
		ifstream file2;
		file2.open(fn_base+label+"/TRAIN/n_obs.txt");
		while (file2 >> value2)
		{
			n[no] = value2;   
			no++;
		}
		file2.close();
		int sum_n = 0;
		//count number of observations 
		for(int j=0; j<N_train; j++){ 
			sum_n = sum_n + n[j];

		}
				
		Z.resize(p*sum_n);
		int zl = 0;
		double value3;
		ifstream file3;
		file3.open(fn_base+label+"/TRAIN/Zlong_proth.txt");
		while (file3 >> value3)
		{
			Z[zl] = value3;   
			zl++;
		}
		file3.close();

		int zf = 0;
		int value4; //will need to change if fixed covariates aren't all int
		ifstream file4;
		file4.open(fn_base+label+"/TRAIN/Zfix_drug.txt");
		while (file4 >> value4)
		{
			Z_fix[zf] = (double)value4;   
			zf++;
		}
		file4.close();

		t_obs.resize(sum_n);
		int to = 0;
		double value5;
		ifstream file5;
		file5.open(fn_base+label+"/TRAIN/t_obs.txt");
		while (file5 >> value5)
		{
			t_obs[to] = value5;   
			to++;
		}
		file5.close();
		
		int c = 0;
		int count = 0;
		int value6;
		ifstream file6;
		file6.open(fn_base+label+"/TRAIN/Censor.txt");
		while (file6 >> value6)
		{
			Censor[c] = value6;   
			if(Censor[c]==0){ count = count + 1; }
			c++;
		}
		file6.close();
		
		//Create vector of switch times
		T.resize(N_train+sum_n);
		for(int j=0; j<N_train; j++){ 
			int sum_nj = 0;
			for(int i=0; i<j; i++){
				sum_nj = sum_nj + n[i];
			} 
			//define vector of observation times of individual i (length n[i])
			VecDoub t_o(n[j]);
			for(int i=0; i<n[j]; i++){
				t_o[i] = t_obs[sum_nj+i];
			}
			//switch times
			T[sum_nj + j] = 0.0; 
			T[sum_nj + j + n[j]] = t_o[n[j]-1];
			for(int b=1; b<n[j]; b++){
				T[sum_nj+j+b] = (t_o[b-1] + t_o[b])/2.0;
			}			
		}

		
		//---------------------------------------------------------------------------------------------------------------------------
		//READ IN TEST DATA LIVER----------------------------------------------------------------------------------------------------
		//---------------------------------------------------------------------------------------------------------------------------

		int eva = 0;
		double value1a;
		ifstream filea;
		filea.open(fn_base+label+"/TEST/Events.txt");
		while (filea >> value1a)
		{
			X_test[eva] = value1a;   
			eva++;
		}
		filea.close();
			
		int noa = 0;
		int value2a;
		ifstream file2a;
		file2a.open(fn_base+label+"/TEST/n_obs.txt");
		while (file2a >> value2a)
		{
			n_test[noa] = value2a;   
			noa++;
		}
		file2a.close();
		int sum_n_test = 0;
		//count number of observations 
		for(int j=0; j<N_test; j++){ 
			sum_n_test = sum_n_test + n_test[j];

		}
		
		Z_test.resize(p*sum_n_test);
		int zla = 0;
		double value3a;
		ifstream file3a;
		file3a.open(fn_base+label+"/TEST/Zlong_proth.txt");
		while (file3a >> value3a)
		{
			Z_test[zla] = value3a;  
			zla++;
		}
		file3a.close();
		
		int zfa = 0;
		int value4a; //will need to change if fixed covariates aren't all int
		ifstream file4a;
		file4a.open(fn_base+label+"/TEST/Zfix_drug.txt");
		while (file4a >> value4a)
		{
			Z_fix_test[zfa] = (double)value4a;   
			zfa++;
		}
		file4a.close();

		t_obs_test.resize(sum_n_test);
		int ta = 0;
		double value5a;
		ifstream file5a;
		file5a.open(fn_base+label+"/TEST/t_obs.txt");
		while (file5a >> value5a)
		{
			t_obs_test[ta] = value5a;  
			ta++;
		}
		file5a.close();

		int ca = 0;
		int counta = 0;
		int value6a;
		ifstream file6a;
		file6a.open(fn_base+label+"/TEST/Censor.txt");
		while (file6a >> value6a)
		{
			Censor_test[ca] = value6a;   
			if(Censor_test[ca]==0){ counta = counta + 1; }
			ca++;
		}
		file6a.close();

		//Create switch times of test data
		T_test.resize(N_test+sum_n_test);
		for(int j=0; j<N_test; j++){ 
			int sum_nj = 0;
			for(int i=0; i<j; i++){
				sum_nj = sum_nj + n_test[i];
			} 
			//define observation times of individual j length n_test[j]
			VecDoub t_o(n_test[j]);
			for(int i=0; i<n_test[j]; i++){
				t_o[i] = t_obs_test[sum_nj+i];
			}
			//switch times
			T_test[sum_nj + j] = 0.0; 
			T_test[sum_nj + j + n_test[j]] = t_o[n_test[j]-1];
			for(int b=1; b<n_test[j]; b++){
				T_test[sum_nj+j+b] = (t_o[b-1] + t_o[b])/2.0;
			}
		}

	
		cout << "Data read in " << endl;


		//---------------------------------------------------------------------------------------------------------------------------
		// FIT MODEL B --------------------------------------------------------------------------------------------------------------
		//---------------------------------------------------------------------------------------------------------------------------

		VecDoub kappa(p);
		VecDoub tau(p);
		VecDoub gamma(q);
		try{

			//dimensions = number of long. covariates*number of parameters in beta parameterisation + no. of fixed covariates
			int DIM = 2*p + q; 

			//Initialise parameter estimates 
			//order is kappa(proth), tau(proth), gamma(drug)
			Doub initx[] = { 
				0.0, 0.1, 0.0
			};

			VecDoub v(DIM, initx); // Initialise the vector with the initial guess
			VecDoub vv;

			// Feed the function to the constructor
			Powell <Doub(VecDoub_I &)> powell(func);

			// Perform the minimization
			vv = powell.minimize(v); 

			//final (minimum) value of function = powell.fret
			cout << "-F(end) = "  << powell.fret << endl;
			
			//maximum likelihood estimates of model parameters
			cout << "estimated k " << endl;
			for(int i = 0; i<p; i++){
				cout << setprecision(15) << vv[i] << ",  ";
				kappa[i] = vv[i];
			}
			
			cout << endl << "estimated t " << endl;
			for(int i = 0; i<p; i++){
				cout << setprecision(15) << vv[p+i] << ",  ";
				tau[i] = vv[p+i];
			}
			
			cout << endl << "estimated g " << endl;
			for(int i = 0; i<q; i++){
				cout << setprecision(15) << vv[2*p+i] << ",  ";
				gamma[i] = vv[2*p+i];
			}
			cout << endl << endl;


			//-------------------------------------------------------------------------------------------------------------------
			// FIX T ------------------------------------------------------------------------------------------------------------
			//-------------------------------------------------------------------------------------------------------------------

			//for fixed t=3 years, iterate over prediction time u = 3-10 years in steps of 0.2 years
  			//calculate PE for Model B for each combination of (t,u) using newdata = test data

			cout << "begin fix t " << endl;

			ofstream filePE;
			filePE.open(fn_base+label+"/ModelB_fixt.txt");

			double baseTime = 3.0;
			for(int l=0; l<36; l++){
				double predTime = baseTime+(double)l*0.2;
				double predErr = PE_SQ(predTime, baseTime, kappa, tau, gamma, X, X_test, Censor, Censor_test, Z, Z_test, Z_fix, Z_fix_test, t_obs, t_obs_test, n, n_test, N_train, N_test);
				filePE << setprecision(15) << predErr << "\n";
			}
			filePE.close();

			cout << "fix t done" << endl;

			//-------------------------------------------------------------------------------------------------------------------
			// Pred Window ------------------------------------------------------------------------------------------------------
			//-------------------------------------------------------------------------------------------------------------------
			double w1 = 1.0;
			double w2 = 2.0;
			double w3 = 3.0;

			//For fixed window w1=1 year work out PE for base time t=0-9 years (steps of 0.2)

			cout << "begin window 1 " << endl;
			
			ofstream filePE1;
			filePE1.open(fn_base+label+"/ModelB_w1.txt");

			for(int l=0; l<46; l++){
				double t = (double)l*0.2;
				double predErr = PE_SQ(t + w1, t, kappa, tau, gamma, X, X_test, Censor, Censor_test, Z, Z_test, Z_fix, Z_fix_test, t_obs, t_obs_test, n, n_test, N_train, N_test);
				filePE1 << setprecision(15) << predErr << "\n";
			}
			filePE1.close();

			cout << "done" << endl;

			//For fixed window w2=2 years work out PE for base time t=0-8 years (steps of 0.2)

			cout << "begin window 2 " << endl;

			ofstream filePE2;
			filePE2.open(fn_base+label+"/ModelB_w2.txt");

			for(int l=0; l<41; l++){
				double t = (double)l*0.2;
				double predErr = PE_SQ(t + w2, t, kappa, tau, gamma, X, X_test, Censor, Censor_test, Z, Z_test, Z_fix, Z_fix_test, t_obs, t_obs_test, n, n_test, N_train, N_test);
				filePE2 << setprecision(15) << predErr << "\n";
			}
			filePE2.close();

			cout << "done" << endl;

			//For fixed window w3=3 years work out PE for base time t=0-7 years (steps of 0.2)

			cout << "begin window 3" << endl;

			ofstream filePE3;
			filePE3.open(fn_base+label+"/ModelB_w3.txt");

			for(int l=0; l<36; l++){
				double t = (double)l*0.2;
				double predErr = PE_SQ(t + w3, t, kappa, tau, gamma, X, X_test, Censor, Censor_test, Z, Z_test, Z_fix, Z_fix_test, t_obs, t_obs_test, n, n_test, N_train, N_test);
				filePE3 << setprecision(15) << predErr << "\n";
			}
			filePE3.close();	

			cout << "done" << endl;

			
		}catch(...){
			
			cout << "can't minimise carry on" << endl;
		}	
		
	}	

	time(&end);

	double time_taken = double(end-start); 
    	cout << "Time taken by program is : " << fixed << time_taken << setprecision(5); 
    	cout << " sec " << endl; 
	

	return 0;	
	

}



