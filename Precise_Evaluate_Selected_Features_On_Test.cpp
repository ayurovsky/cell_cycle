#include <array>
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
	string test_file_prefix = argv[1];	
	string train_file_prefix = argv[2];	
	string phase = argv[3];	
	double total_range = stod(argv[4]);	
	int k = stoi(argv[5]);
	for (int i = 6; i < argc; ++i){
		int gene_index = stoi(argv[i]);
		available_train_gene_indeces.push_back(gene_index);
	}
	vector<int> can_use_train_gene_indeces;

	// read in train gene names
	vector<string> canonical_train_gene_names;
	ifstream file(test_file_prefix + "_genes.csv");
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

	// read in the test angles 
	vector<double> test_cell_angles;
	vector<string> read2;
	ifstream file6(test_file_prefix + "_angles.csv");
	getline(file6, line);
	stringstream tokenStream3(line);
	string token3;
	while (getline(tokenStream3, token3, ',')){
		test_cell_angles.push_back(stod(token3)); 
	}
	file6.close();

	// read in test gene names
	vector<string> canonical_test_gene_names;
	vector<string> can_use_test_gene_names;
	ifstream file0(test_file_prefix + "_genes.csv");
	getline(file0, line);
	boost::algorithm::split(canonical_test_gene_names, line, boost::is_any_of(","));
	file0.close();
	can_use_test_gene_names = canonical_test_gene_names;
	vector<int> can_use_test_gene_indeces;


	// read in test cell names
	vector<string> test_cell_names;
	ifstream file3(test_file_prefix + "_cells.csv");
	//ifstream file3("SC_df_z_train_cells.csv");
	getline(file3, line);
	boost::algorithm::split(test_cell_names, line, boost::is_any_of(","));
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

	// read in the test data matrix - remove any suspect genes - ones with nan for example
	vector<vector<double>> testDataMatrix;
	ifstream file5(test_file_prefix + "_array.csv");		
	while(getline(file5, line))
	{
		vector<double> row;
		stringstream tokenStream(line);
		string token;
		int counter = 0;
		while (getline(tokenStream, token, ' ')){
			if (token == "nan") {
				row.push_back(0.0);
				//cout<<"removing gene from can use set: "<<canonical_test_gene_names[counter]<<endl;
				can_use_test_gene_names.erase(remove(can_use_test_gene_names.begin(),\
				can_use_test_gene_names.end(), canonical_test_gene_names[counter]), can_use_test_gene_names.end());
			} else {	
				row.push_back(stod(token));
			}
			counter += 1;
		}
		testDataMatrix.push_back(row);
	}
	file5.close();

	// get intersection between genes in train and test sets
	for (int j = 0; j < available_train_gene_indeces.size(); j++) {
		int i = available_train_gene_indeces[j];
		//cout<<j<<" ";
		//cout<<i<<" ";
		string gene_name = canonical_train_gene_names[i];
		//cout<<gene_name<<" ";
		// see if this name is in can_use_test_gene_names
		int pos = find(can_use_test_gene_names.begin(), can_use_test_gene_names.end(), gene_name) - can_use_test_gene_names.begin();
		if (pos < can_use_test_gene_names.size()) {
			// add this index to train indeces we can use
			can_use_train_gene_indeces.push_back(i);
			// find this gene name in the canonical test names list and get the index 
			int pos_test = find(canonical_test_gene_names.begin(), canonical_test_gene_names.end(), gene_name) - canonical_test_gene_names.begin();
			if (pos_test >= canonical_test_gene_names.size())
				cout<<"Major bug"<<endl;
			can_use_test_gene_indeces.push_back(pos_test);
		}	
	}
	
	/*cout<<"train gene names we can use"<<endl;
	for (auto i = can_use_train_gene_indeces.begin(); i != can_use_train_gene_indeces.end(); ++i)
    	cout << canonical_train_gene_names[*i] << ' ';
	cout<<endl<<"test gene names we can use"<<endl;
	for (auto i = can_use_test_gene_indeces.begin(); i != can_use_test_gene_indeces.end(); ++i)
    	cout << canonical_test_gene_names[*i] << ' ';
	cout<<endl;	*/
	if (can_use_train_gene_indeces.size() != can_use_test_gene_indeces.size()) {
		cout<< "Major bug!, sizes of indeces for test and train do not match!!!"<<endl;
	}
	int gene_indeces_size = can_use_train_gene_indeces.size();
	cout<<"size of gene intersection is "<<gene_indeces_size<<endl;

	//for (int m = 0; m <= gene_indeces_size; m++)
	//	cout<<testDataMatrix[0][m]<<" "<<trainDataMatrix[0][m]<<endl;

		// go through all the test cells one by one 	
	for(int i = 0; i < test_cell_names.size(); i++) { //30; i++) {
		vector<tuple<double, double> > dist_names_list;
		vector<double> test_cell = testDataMatrix[i];
		if (phase == "g1" and not strstr(test_cell_names[i].c_str(), "_G1"))
			continue;
		if (phase == "g2" and not strstr(test_cell_names[i].c_str(), "_G2"))
			continue;
		if (phase == "s" and not strstr(test_cell_names[i].c_str(), "_S"))
			continue;
		// go through all the train cells one by one
		for(int j = 0; j < train_cell_angles.size(); j++){
			vector<double> train_cell = trainDataMatrix[j];
			//cout<<"Train cell is "<<cell_names[j]<<endl;

			// select two vectors with distance
			vector<double> test_vector(gene_indeces_size, 0);	
			vector<double> train_vector(gene_indeces_size, 0);	

			transform(can_use_test_gene_indeces.begin(), can_use_test_gene_indeces.end(), \
				test_vector.begin(), [test_cell](size_t pos) {return test_cell[pos];});
			transform(can_use_train_gene_indeces.begin(), can_use_train_gene_indeces.end(),  \
				train_vector.begin(), [train_cell](size_t pos) {return train_cell[pos];});

			// euclidean distance	
			vector<double>	auxiliary;
			transform(test_vector.begin(), test_vector.end(), train_vector.begin(), back_inserter(auxiliary),//
			[](double element1, double element2) {return pow((element1-element2),2);});
			double distance = sqrt(accumulate(auxiliary.begin(), auxiliary.end(), 0.0));
		/*	if (i == j) {
				cout<<distance<<endl;
			}*/

			// record distance 
			dist_names_list.push_back(make_tuple(distance, train_cell_angles[j]));
				
			//break;
		}
		// sort the training cells based on distances0, get the best one for now
		cout<<"got this many comparisons: "<<dist_names_list.size()<<endl;
		sort(dist_names_list.begin(), dist_names_list.end());
		for (int l=0; l<dist_names_list.size(); l++)
        	cout<<get<1>(dist_names_list[l])<<" ";
		cout<<endl;	
	
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
		cout<<test_cell_angles[i]<<" "<<avg_angle<<endl;
		double test_error = abs(test_cell_angles[i] - avg_angle)/total_range;
		cout<<"Error is: "<<test_error<<endl;

		// now get the p-value for the error
		auto rng = default_random_engine {};
		double counter = 0.0;
		for (int p=0; p<100; p++) { 
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
			double rand_error = abs(test_cell_angles[i] - avg_angle)/total_range;
			cout<<"random error: "<<rand_error<<endl;
			if (rand_error < test_error)
				counter += 1;
		}
		cout<<"p-value: "<<counter/100<<endl;
	}
}


