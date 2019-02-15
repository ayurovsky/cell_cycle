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
	int k = stoi(argv[2]);
	for (int i = 3; i < argc; ++i){
		int gene_index = stoi(argv[i]);
		available_train_gene_indeces.push_back(gene_index);
	}
	vector<int> can_use_train_gene_indeces;

	// read in train gene names
	vector<string> canonical_train_gene_names;
	ifstream file("SC_df_z_train_genes.csv");
	string line = "";
	getline(file, line);
	boost::algorithm::split(canonical_train_gene_names, line, boost::is_any_of(","));
	file.close();

	// read in train cell names
	vector<string> train_cell_names;
	ifstream file2("SC_df_z_train_cells.csv");
	getline(file2, line);
	boost::algorithm::split(train_cell_names, line, boost::is_any_of(","));
	file2.close();

	// read in test gene names
	vector<string> canonical_test_gene_names;
	vector<string> can_use_test_gene_names;
	ifstream file0(test_file_prefix + "_test_genes.csv");
	getline(file0, line);
	boost::algorithm::split(canonical_test_gene_names, line, boost::is_any_of(","));
	file0.close();
	can_use_test_gene_names = canonical_test_gene_names;
	vector<int> can_use_test_gene_indeces;

	// read in test cell names
	vector<string> test_cell_names;
	ifstream file3(test_file_prefix + "_test_cells.csv");
	//ifstream file3("SC_df_z_train_cells.csv");
	getline(file3, line);
	boost::algorithm::split(test_cell_names, line, boost::is_any_of(","));
	file3.close();
	

	// read in the train data matrix
	vector<vector<double>> trainDataMatrix;
	ifstream file4("SC_df_z_train_array.csv");		
	while (getline(file4, line))
	{
		vector<double> row;
		stringstream ss(line);
		double d = 0.0;
    	while (ss >> d)
        	row.push_back(d);
		trainDataMatrix.push_back(row);
	}
	file4.close();

	// read in the test data matrix - remove any suspect genes - ones with nan for example
	vector<vector<double>> testDataMatrix;
	ifstream file5(test_file_prefix + "_test_array.csv");		
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
	
	cout<<"train gene names we can use"<<endl;
	for (auto i = can_use_train_gene_indeces.begin(); i != can_use_train_gene_indeces.end(); ++i)
    	cout << canonical_train_gene_names[*i] << ' ';
	cout<<endl<<"test gene names we can use"<<endl;
	for (auto i = can_use_test_gene_indeces.begin(); i != can_use_test_gene_indeces.end(); ++i)
    	cout << canonical_test_gene_names[*i] << ' ';
	cout<<endl;	
	if (can_use_train_gene_indeces.size() != can_use_test_gene_indeces.size()) {
		cout<< "Major bug!, sizes of indeces for test and train do not match!!!"<<endl;
	}

	int gene_indeces_size = can_use_train_gene_indeces.size();
	cout<<"size of gene intersection is "<<gene_indeces_size<<endl;
	double total_correct = 0;
		// go through all the test cells one by one 	
	for(int i = 0; i < test_cell_names.size(); i++) { //30; i++) {
		vector<tuple<double, string> > dist_names_list;
		vector<double> test_cell = testDataMatrix[i];
		cout<<"Test cell is "<<test_cell_names[i]<<endl;
		// go through all the train cells one by one
		for(int j = 0; j < 500; j++){
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
			//cout<<distance<<endl;

			// record distance 
			dist_names_list.push_back(make_tuple(distance, train_cell_names[j])); //cell_type_names[j]));
				
			//break;
		}
		// sort the training cells based on distances, get the best one for now
		sort(dist_names_list.begin(), dist_names_list.end());
	
		double sum_q = 0.0;
		for (int q = 0; q < 10; q++)
			sum_q += get<0>(dist_names_list[0]);
		if (sum_q == 0.0) {
			total_correct += 0.33;
			continue;
		}
		double num_G2M = 0.0;
		double num_G1 = 0.0;
		double num_S = 0.0;
		//cerr<<cell_names[i]<<endl;
		for (int idx = 0; idx <k; idx++) {
			double dist = get<0>(dist_names_list[idx]);
			//cerr<<get<1>(dist_names_list[idx])<<" "<<get<0>(dist_names_list[idx])<<endl;
			if (strstr(get<1>(dist_names_list[idx]).c_str(), "_G2"))
				num_G2M += 1; //(dist) ? 1.0/dist : 0;
			else if (strstr(get<1>(dist_names_list[idx]).c_str(),"_G1"))
				num_G1 += 1; //(dist) ? 1.0/dist : 0; //get<0>(dist_names_list[idx]);
			else if (strstr(get<1>(dist_names_list[idx]).c_str(),"_S")) 
				num_S += 1; //(dist) ? 1.0/dist : 0; //get<0>(dist_names_list[idx]);
		}

		if ( ((num_G2M >= num_G1) and (num_G2M >= num_S) and strstr(test_cell_names[i].c_str(),"_G2")) or
			 ((num_G1 >= num_G2M) and (num_G1 >=  num_S) and strstr(test_cell_names[i].c_str(),"_G1")) or
			 ((num_S >= num_G2M) and (num_S >= num_G1) and  strstr(test_cell_names[i].c_str(),"_S"))) { 
			total_correct += 1;

		}
	}
	cout<<round(total_correct)<<endl;
	cout<<round(total_correct)/test_cell_names.size()<<endl;
}


