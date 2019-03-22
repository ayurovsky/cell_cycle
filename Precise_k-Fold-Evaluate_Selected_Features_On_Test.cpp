#include <array>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <string>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include<numeric>  // inner_product, iota
#include<cmath>    // sqrt
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
using boost::property_tree::ptree;
using boost::property_tree::read_json;
using boost::property_tree::write_json;
using namespace std;


int main(int argc, char *argv[])
{
	vector<int> available_train_gene_indeces;
	string train_file_prefix = argv[1];	
	int k = stoi(argv[2]);
	for (int i = 3; i < argc; ++i){
		int gene_index = stoi(argv[i]);
		available_train_gene_indeces.push_back(gene_index);
	}

	// read in train gene names
	vector<string> canonical_train_gene_names;
	ifstream file(train_file_prefix + "_genes.csv");
	string line = "";
	getline(file, line);
	boost::algorithm::split(canonical_train_gene_names, line, boost::is_any_of(","));
	file.close();


	// read in the train angles 
	vector<double> train_cell_angles;
	vector<string> read;
	ifstream file2(train_file_prefix + "_angles.csv");
	getline(file2, line);
	stringstream tokenStream2(line);
	string token2;
	while (getline(tokenStream2, token2, ',')){
		train_cell_angles.push_back(stod(token2)); 
	}
	file2.close();


	// read in train cell names
	vector<string> train_cell_names;
	ifstream file3(train_file_prefix + "_cells.csv");
	//ifstream file3("SC_df_z_train_cells.csv");
	getline(file3, line);
	boost::algorithm::split(train_cell_names, line, boost::is_any_of(","));
	file3.close();
	

	// read in the train data matrix
	vector<vector<double>> trainDataMatrix;
	ifstream file4(train_file_prefix + "_array.csv");		
	while (getline(file4, line))
	{
		vector<double> row;
		stringstream tokenStream(line);
		string token;
		int counter = 0;
		while (getline(tokenStream, token, ' ')){
			if (token == "nan") {
				row.push_back(0.0);
				//cout<<"removing gene from can use set: "<<canonical_test_gene_names[counter]<<endl;
				available_train_gene_indeces.erase(remove(available_train_gene_indeces.begin(),\
				available_train_gene_indeces.end(), counter), available_train_gene_indeces.end());
			} else {	
				row.push_back(stod(token));
			}
			counter += 1;
		}
		trainDataMatrix.push_back(row);
	}
	file4.close();


	int gene_indeces_size = available_train_gene_indeces.size();
	cout<<"size of gene intersection is "<<gene_indeces_size<<endl;


	srand(time(NULL));
	// go through all the test cells one by one 	
	double total_error = 0.0;
	double total_monkey_error = 0.0;
	double total_g1_error = 0.0;
	double total_g2_error = 0.0;
	double total_s_error = 0.0;
	int num_g1 = 0;
	int num_g2 = 0;
	int num_s = 0;
	double total_pvalue = 0.0;
	vector<string> predicted_angles;
	for(int i = 0; i < train_cell_names.size(); i++) {
		cout<<endl<<"starting for cell "<<i<<endl;
		vector<tuple<double, double> > dist_names_list;
		vector<double> test_cell = trainDataMatrix[i];
		/*string phase = "";
		if (strstr(train_cell_names[i].c_str(), "_G1"))
			phase = "g1";	
			total_range = 5.9 - 4.8;
		if (strstr(train_cell_names[i].c_str(), "_G2"))
			phase = "g2";
			total_range = 4.8 - 2;
		if (strstr(train_cell_names[i].c_str(), "_S"))
			phase = "s";
			total_range = 2.0 + (6.283 - 5.9);*/
		// go through all the train cells one by one
		// train cells are all cells not identical to the test cell
		cout<<"_"<<train_cell_names[i]<<endl;
		for(int j = 0; j < train_cell_names.size(); j++){
			if (i == j)
				continue;
			/*if (not ( (strstr(train_cell_names[i].c_str(), "_G1") and (strstr(train_cell_names[j].c_str(), "_G1")))
				or   (strstr(train_cell_names[i].c_str(), "_G2") and (strstr(train_cell_names[j].c_str(), "_G2")))
				or   (strstr(train_cell_names[i].c_str(), "_S") and (strstr(train_cell_names[j].c_str(), "_S")))))
				continue;*/
			vector<double> train_cell = trainDataMatrix[j];
			//cout<<"Train cell is "<<cell_names[j]<<endl;

			// select two vectors with distance
			vector<double> test_vector(gene_indeces_size, 0);	
			vector<double> train_vector(gene_indeces_size, 0);	

			transform(available_train_gene_indeces.begin(), available_train_gene_indeces.end(), \
				test_vector.begin(), [test_cell](size_t pos) {return test_cell[pos];});
			transform(available_train_gene_indeces.begin(), available_train_gene_indeces.end(),  \
				train_vector.begin(), [train_cell](size_t pos) {return train_cell[pos];});

			// euclidean distance	
			vector<double>	auxiliary;
			transform(test_vector.begin(), test_vector.end(), train_vector.begin(), back_inserter(auxiliary),//
			[](double element1, double element2) {return pow((element1-element2),2);});
			double distance = sqrt(accumulate(auxiliary.begin(), auxiliary.end(), 0.0));

			// record distance 
			dist_names_list.push_back(make_tuple(distance, train_cell_angles[j]));
				
			//break;
		}
		cout<<"got this many comparisons: "<<dist_names_list.size()<<endl;
		// sort the training cells based on distances, get the best one for now
		sort(dist_names_list.begin(), dist_names_list.end());
	
		double sine_sum = 0.0;
		double cos_sum = 0.0;
		//cerr<<cell_names[i]<<endl;
		for (int idx = 0; idx <k; idx++) {
			cout<<get<1>(dist_names_list[idx])<<endl;
			sine_sum += sin(get<1>(dist_names_list[idx]));
			cos_sum += cos(get<1>(dist_names_list[idx]));
		}
		double avg_angle = atan2(sine_sum, cos_sum);
		if (avg_angle < 0)
			avg_angle = 2 * M_PI + avg_angle;
		cout<<train_cell_angles[i]<<" "<<avg_angle<<endl;
		double test_error = (M_PI - abs(abs(train_cell_angles[i] - avg_angle) - M_PI))/M_PI;
		cout<<"Error is: "<<test_error<<endl;
		total_error += test_error;
		
		// add to predicted labels
		predicted_angles.push_back(to_string(avg_angle));
	
		// Now get the monkey error
		double sum_random_error = 0.0;
		int num_repeat = 1000;
		//random.seed();
		for (int p=0; p<num_repeat; p++) { 

			avg_angle =  ((float)rand() / RAND_MAX) * 2 * M_PI;  
			double rand_error = (M_PI - abs(abs(train_cell_angles[i] - avg_angle) - M_PI))/M_PI;
			sum_random_error += rand_error;
		}
		double cell_random_error = sum_random_error/num_repeat;
		cout<<"Random Monkey error is: "<<cell_random_error;
		total_monkey_error += cell_random_error;

		// now get the p-value for the error
		/*auto rng = default_random_engine {};
		double counter = 0.0;
		int num_permute = 10000;
		for (int p=0; p<num_permute; p++) { 
			shuffle(dist_names_list.begin(), dist_names_list.end(), rng);

			sine_sum = 0.0;
			cos_sum = 0.0;
			//cerr<<cell_names[i]<<endl;
			for (int idx = 0; idx <k; idx++) {
				sine_sum += sin(get<1>(dist_names_list[idx]));
				cos_sum += cos(get<1>(dist_names_list[idx]));
			}
			avg_angle = atan2(sine_sum, cos_sum);
			if (avg_angle < 0)
				avg_angle = 2 * M_PI + avg_angle;
			double rand_error = abs(train_cell_angles[i] - avg_angle)/total_range;
			if (rand_error < test_error)
				counter += 1;
		}
		total_pvalue += counter/num_permute;
		cout<<"p-value for "<<i<<" is: "<<counter/num_permute<<endl;*/
		/*if (phase == "g1") {
			num_g1 += 1;
			total_g1_error += test_error;
		} else if (phase == "g2") {
			num_g2 += 1;
			total_g2_error += test_error;
		} else if (phase == "s") {
			num_s += 1;
			total_s_error += test_error;
		}*/
			
		//break;
	}
	cout<<"Done"<<endl;
	cout<<"Average error: "<<total_error/(train_cell_names.size()-1)<<endl;
	cout<<"Average Monkey error: "<<total_monkey_error/(train_cell_names.size()-1)<<endl;
	//cout<<"Average pvalue: "<<total_pvalue/(train_cell_names.size()-1)<<endl;
	/*cout<<"G1 Average Error: "<<total_g1_error/num_g1<<endl;
	cout<<"G2 Average Error: "<<total_g2_error/num_g2<<endl;
	cout<<"S Average Error: "<<total_s_error/num_s<<endl;*/
	ofstream outfile(train_file_prefix + "_predicted_angles.csv");
	string joined_output = boost::algorithm::join(predicted_angles, ",");
	outfile<<joined_output<<"\n";
	outfile.close();
}


