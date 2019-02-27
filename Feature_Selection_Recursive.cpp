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
	vector<int> gene_indeces;
	int k = stoi(argv[1]);
	string prefix = argv[2];
	//cout<<prefix<<endl;
	for (int i = 3; i < argc; ++i){
		int gene_index = stoi(argv[i]);
		gene_indeces.push_back(gene_index);
	}
	// read in gene names
	vector<string> gene_names;
	ifstream file(prefix + "_genes.csv");
	string line = "";
	getline(file, line);
	boost::algorithm::split(gene_names, line, boost::is_any_of(","));
	file.close();

	// read in cell names
	vector<string> cell_names;
	ifstream file2(prefix + "_cells.csv");
	getline(file2, line);
	boost::algorithm::split(cell_names, line, boost::is_any_of(","));
	file2.close();
	
	// read in the data matrix
	vector<vector<double>> dataMatrix;
	ifstream file3(prefix + "_array.csv");		
	while (getline(file3, line))
	{
		vector<double> row;
		stringstream tokenStream(line);
		string token;
		int counter = 0;
		while (getline(tokenStream, token, ' ')){
			if (token == "nan") {
				row.push_back(0.0);
				//cout<<"removing gene from can use set: "<<canonical_test_gene_names[counter]<<endl;
				gene_indeces.erase(remove(gene_indeces.begin(),\
				gene_indeces.end(), counter), gene_indeces.end());
			} else {	
				row.push_back(stod(token));
			}
			counter += 1;
		}
		dataMatrix.push_back(row);
	}
	file3.close();

	//cerr<<"gene indeces size: "<<gene_indeces.size()<<endl;
	int total_correct = 0;
	if (gene_indeces.size() >  0) {
		// go through all the test cells one by one 	
		for(int i = 0; i < cell_names.size(); i++) {
			vector<tuple<double, string> > dist_names_list;
			vector<double> test_cell = dataMatrix[i];
			// go through all the train cells one by one
			for(int j = 0; j < cell_names.size(); j++){
				vector<double> train_cell = dataMatrix[j];
				//cout<<"Train cell is "<<cell_names[j]<<endl;
				if (i == j) 
					continue;

				// select two vectors with distance
				vector<double> test_vector(gene_indeces.size(), 0);	
				vector<double> train_vector(gene_indeces.size(), 0);	
				transform(gene_indeces.begin(), gene_indeces.end(), test_vector.begin(), [test_cell](size_t pos) {return test_cell[pos];});
				transform(gene_indeces.begin(), gene_indeces.end(), train_vector.begin(), [train_cell](size_t pos) {return train_cell[pos];});
				// euclidean distance	
				vector<double>	auxiliary;
				transform(test_vector.begin(), test_vector.end(), train_vector.begin(), back_inserter(auxiliary),//
				[](double element1, double element2) {return pow((element1-element2),2);});
				double distance = sqrt(accumulate(auxiliary.begin(), auxiliary.end(), 0.0));
				//cout<<distance<<endl;

				// record distance 
				dist_names_list.push_back(make_tuple(distance, cell_names[j])); //cell_type_names[j]));
					
				//break;
			}
			// sort the training cells based on distances, get the best one for now
			sort(dist_names_list.begin(), dist_names_list.end());
		
			/*string best_train_name = get<1>(dist_names_list[0]);
			if ((strstr(best_train_name.c_str(),"_G2")	and strstr(cell_names[i].c_str(),"_G2")) or
				(strstr(best_train_name.c_str(),"_G1") and strstr(cell_names[i].c_str(),"_G1")) or
				(strstr(best_train_name.c_str(),"_S") and strstr(cell_names[i].c_str(),"_S"))) {
					total_correct += 1;
			}*/
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
			//cerr<<num_G2M<<"_"<<num_G1<<"_"<<num_S<<endl;
			if ( !num_G2M and !num_G1 and !num_S and gene_indeces.size() > 1) {
				cerr<<"Zero for num genes "<<gene_indeces.size()<<endl; 
				cerr<<get<1>(dist_names_list[0]).c_str()<<" "<<get<1>(dist_names_list[1]).c_str()<<" "<<get<1>(dist_names_list[2]).c_str()<<endl;
			}
			//	continue;

			if ( ((num_G2M >= num_G1) and (num_G2M >= num_S) and strstr(cell_names[i].c_str(),"_G2")) or
				 ((num_G1 >= num_G2M) and (num_G1 >=  num_S) and strstr(cell_names[i].c_str(),"_G1")) or
				 ((num_S >= num_G2M) and (num_S >= num_G1) and  strstr(cell_names[i].c_str(),"_S"))) { 
				total_correct += 1;

			}
		}
	}
	cout<<round(total_correct);
}


