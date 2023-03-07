/* ===============================================================================================
   MatInOut.h

   Small functions that read in a matrix of values.
  
   Authors:  Yossi Bokor
   Date:    March 7, 2022
   Version: 1

   Methods:
 
   =============================================================================================== */

#ifndef _MATINOUT_H_
#define _MATINOUT_H_

#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

class MATINOUT {

	public:
		vector<vector<double>> read(string infile);

};

vector<vector<double>> read(string infile) {
   vector<vector<double>> X;
   ifstream in(infile);
   string line;
   double val;

   while (getline(in, line)){
       istringstream iss(line);
       vector<double> x_i;
       while (iss >> val){
          x_i.push_back(val);
       };
       X.push_back(x_i);
   };
   return X;
};

#endif
