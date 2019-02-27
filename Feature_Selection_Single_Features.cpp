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


int main()
{
//#if __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1900)
//    cout << "C++11 is supported" << "\n";
//#else
//    cout << "C++11 is not supported" << \n";
//#endif

	// read in gene names
	vector<string> gene_names;
	ifstream file("SC_df_z_train_genes.csv");
	string line = "";
	getline(file, line);
	boost::algorithm::split(gene_names, line, boost::is_any_of(","));
	/*for(string data : gene_names)
	{
		cout<<data << " , ";
	}
	cout<<endl;
	cout<<gene_names.size()<<endl;
	cout<<gene_names[0]<<endl;	*/
	file.close();

	// read in cell names
	vector<string> cell_names;
	ifstream file2("SC_df_z_train_cells.csv");
	getline(file2, line);
	boost::algorithm::split(cell_names, line, boost::is_any_of(","));
	/*for(string data : cell_names)
	{
		cout<<data << " , ";
	}
	cout<<endl;
	cout<<cell_names.size()<<endl;
	cout<<cell_names[0]<<endl;	*/
	file2.close();

	// read in the data matrix
	vector<vector<double>> dataMatrix;
	ifstream file3("SC_df_z_train_array.csv");		
	while (getline(file3, line))
	{
		vector<double> row;
		stringstream ss(line);
		double d = 0.0;
    	while (ss >> d)
        	row.push_back(d);
		dataMatrix.push_back(row);
	}
	/*cout<<dataMatrix.size()<<endl;
	cout<<dataMatrix[0].size()<<endl;
	cout<<dataMatrix[499][251]<<endl;
	cout<<dataMatrix[0][0]<<endl;*/
	file3.close();

	map<string,double> final_accuracies;

	for(int g = 0; g < gene_names.size(); g++) {

		vector<int> gene_indeces = {g};
		double total_correct = 0;
		// go through all the test cells one by one 	
		for(int i = 0; i < 500; i++) {
			vector<tuple<double, string> > dist_names_list;
			vector<double> test_cell = dataMatrix[i];
			// go through all the train cells one by one
			for(int j = 0; j < 500; j++){
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
				dist_names_list.push_back(make_tuple(distance, cell_names[j]));
					
				//break;
			}
			// sort the training cells based on distances, get the best one for now
			sort(dist_names_list.begin(), dist_names_list.end());

			// if the first 10 distances are all 0, don't count this cell
			double sum_q = 0.0;
			for (int q = 0; q < 10; q++)
				sum_q += get<0>(dist_names_list[0]);
			if (sum_q == 0.0) {
				total_correct += 0.33;
			} else {
				string best_train_name = get<1>(dist_names_list[0]);
				if ((strstr(best_train_name.c_str(),"_G2")	and strstr(cell_names[i].c_str(),"_G2")) or
					(strstr(best_train_name.c_str(),"_G1") and strstr(cell_names[i].c_str(),"_G1")) or
					(strstr(best_train_name.c_str(),"_S") and strstr(cell_names[i].c_str(),"_S"))) {
						total_correct += 1;
				}
			}
			//cout<<"Total_correct "<<total_correct<<endl;
		}
		final_accuracies.insert(pair<string,double>(gene_names[g],round(total_correct)));
		//cout<<"Accuracy for gene "<<gene_names[g]<<" is "<<total_correct<<endl;
	}
	ptree pt; 
	for (auto& entry: final_accuracies) 
      pt.put (entry.first, entry.second);
	ofstream file4("cpp_single_features_accuracy_on_training_set.json");
  	write_json(file4, pt, false); 
	file4.close();
}


