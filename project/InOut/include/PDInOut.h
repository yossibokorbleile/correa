/*
   ReadPersistenceDiagram.h

   Small functions that read and write a persistence diagram. This is based on InOut.h by Patrice Koehl.
  
   Authors:  Yossi Bokor
   Date:    January 17, 2022
   Version: 0.1

   Methods:
 
  */

#ifndef _PDINOUT_H_
#define _PDINOUT_H_


#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <fstream>

/*! \brief readPD reads in an existing persistence diagram
*/
void readPD(std::string inFile, int ndim, int *npoints, double *coord){
	
	//Create input file stream "input" which is associated with the file name inFile
	std::ifstream input;
	input.open(inFile);
	
	//The file is read line by line
	std::string line;
	
	int nval = 0;
	double val;
	
	//For each line, we neet to extract the coordinates of the represented point
	while(getline(input, line)){
		//Create a string stream based on the line
		std::istringstream iss(line);
		//Then loop over the ndim coordinates
		for(int i = 0; i < ndim; i++){
			iss >> val;
			coord[ndim*nval+1]=val;
		};
		nval++;
	};
	*npoints = nval;
};

#endif
