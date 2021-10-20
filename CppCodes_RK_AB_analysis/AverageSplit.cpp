/*
Code to read in prediction error vector files for all models (joint, landmarking and retarded kernel) for each loop and average over the loops
Each loop represents a random split of the data into training and test data sets
Same code can be used for all data sets - AIDS, Liver, PBC - folder and model names should be changed accordingly
*/
#include <iostream>
#include <iomanip>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string> 
#include <algorithm>
#include <cstdlib>
#include<ctime>
#include<vector>

using namespace std;

//function to read in the results from each loop
//model is either JM, LM, ModelA or ModelB
//type is either "fixt" "w1" "w2" or "w3" 
vector<double> read(string fn_base, string label, string model, string type, int L){

	vector<double> values(L);
	int a = 0;
	string val;
	ifstream file;
	file.open(fn_base + label + "/" + model + "_" + type + ".txt");
	cout << "filename " << label + "/" + model + "_" + type + ".txt" << endl;
	while (file >> val)
	{
		 //we should only get NA for JM t=u (Tstart=Thoriz)
		//i.e. JM, fixt, a=0
		//we can treat this as 0. 
		if(val == "NA"){
			val = 0.0;
			cout << "NA at a = " << a << endl; //to check we don't get NA anywhere else
		}
		else{
			values[a] = stod(val);   
		}
		a++;
	}
	file.close();
	return values;

}

//Function to perform average and create/write to files
//L = length of vector of PE values (fixt, w1, w2 or w3)
//rep = no. of loops over which to average 
void Average(string fn_base, string model, string type, int L, int rep){
	ofstream file_av;
	file_av.open("Av_" + model + "_" + type + ".txt");
	
	vector<double> Av_vals(L);
	for(int loop=0; loop<rep; loop++){

		cout << "**************loop = " << loop+1 << endl;

		string label = "LOOP"+to_string(loop+1);
		vector<double> Loop_vals = read(fn_base, label, model, type, L);

		for(int i=0; i<L; i++){
			Av_vals[i] = Av_vals[i] + Loop_vals[i];
		}

	}

	for(int i=0; i<L; i++){
		Av_vals[i] = Av_vals[i]/(double)rep;
		file_av << setprecision(15) << Av_vals[i] << "\n";
	}
	file_av.close();

}

int main(){

	int rep = 20; //number of loops (splits of data) over which to average

	//length of fixt, window 1, 2, and 3 vectors (change as appropriate)
	int Lf = 26;
	int Lw1 = 46;
	int Lw2 = 41;
	int Lw3 = 36;

	//base folder where results are stored (replace 'DATA_RESULTS' with appropriate folder name for chosen data)
	string fn_base = "~.../DATA_RESULTS/";
		
	//Perform average for the different models
	//if different versions of the models are required - edit the model name accordingly (e.g. JM_linear_cen or ModelA_n1Decay etc)

	//JM 
	Average(fn_base, "JM", "fixt", Lf, rep);
	Average(fn_base, "JM", "w1", Lw1, rep);
	Average(fn_base, "JM", "w2", Lw2, rep);
	Average(fn_base, "JM", "w3", Lw3, rep);
	
	//LM:
	Average(fn_base, "LM", "fixt", Lf, rep);
	Average(fn_base, "LM", "w1", Lw1, rep);
	Average(fn_base, "LM", "w2", Lw2, rep);
	Average(fn_base, "LM", "w3", Lw3, rep);
	
	//ModelA: 
	Average(fn_base, "ModelA", "fixt", Lf, rep);
	Average(fn_base, "ModelA", "w1", Lw1, rep);
	Average(fn_base, "ModelA", "w2", Lw2, rep);
	Average(fn_base, "ModelA", "w3", Lw3, rep);

	//ModelB: 
	Average(fn_base, "ModelB", "fixt", Lf, rep);
	Average(fn_base, "ModelB", "w1", Lw1, rep);
	Average(fn_base, "ModelB", "w2", Lw2, rep);
	Average(fn_base, "ModelB", "w3", Lw3, rep);

	
	return 0;

}
