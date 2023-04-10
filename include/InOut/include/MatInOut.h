/*
   MatInOut.h

 	Authors: Patrice Koehl, Department of Computer Science, University of California, Davis
				Yossi Bokor Bleile, Department of Mathematical Sciences, University of Aalborg, Aalborg
 	Date: April 2023
	Version: 1
*/


#ifndef _MATINOUT_H_
#define _MATINOUT_H_

#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;
namespace correa {


class MATINOUT {

	public:
		vector<vector<double> > read(string infile);

};

vector<vector<double> > read(string infile) {
   vector<vector<double> > X;
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
} //end namespace correa
#endif
