/*!
* @file	InOut.h
* @brief read and write data
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
*/

#ifndef _INOUT_H_
#define _INOUT_H_

#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

// Check if Python bindings are being built
#ifdef CORREA_PYTHON_BINDINGS
namespace correa {
	extern bool g_verbose;
}
#define CORREA_VERBOSE_PRINT(x) if (correa::g_verbose) { x; }
#else
#define CORREA_VERBOSE_PRINT(x) x;
#endif

namespace correa {
  /* ===============================================================================================
   class
   =============================================================================================== */

   class INOUT {

	public:
		void read(std::string infile, int *ndim, int *npoint, double **Coord);

		void write(std::string outfile, int ndim, int npoint, double *Coord);

	private:

		bool exists(const std::string& name);

		void countPoints(std::ifstream & inFile, int *ndim, int *npoint);

		void readPoints(std::ifstream & inFile, int ndim, int npoint, double *coord);

		std::string replaceAll(std::string str, const std::string &from, const std::string &to);

   };

  /* ===============================================================================================
   Checks if a file is accessible
   =============================================================================================== */

 bool INOUT::exists(const std::string& name)
 {
	struct stat buffer;   
	return (stat (name.c_str(), &buffer) == 0); 
 }

  /* ===============================================================================================
   Replace all instances of a substring by another string
   =============================================================================================== */

  std::string INOUT::replaceAll(std::string str, const std::string &from, const std::string &to)
  {
	size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos)
	{
		str.replace(start_pos, from.length(), to);
		start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
	}
	return str;
  }

  /* ===============================================================================================
   Read input file
   =============================================================================================== */

  void INOUT::read(std::string infile, int *ndim, int *npoint, double **Coord) 
  {

	if(!exists(infile)) {
		std::cout << "File " << infile << " does not exists...Please check command" << std::endl;
		exit(1);
	}

	std::ifstream input;
	input.open(infile);

	int n1, n2;
	countPoints(input, &n1, &n2);
	*ndim = n1;
	*npoint = n2;

	input.clear();
	input.seekg(0);

	*Coord = new double[n1*n2];
	readPoints(input, n1, n2, *Coord);
	CORREA_VERBOSE_PRINT(std::cout << "Number of points loaded: " << n2 << std::endl)
	input.close();

  }

/* ===============================================================================================
   Count # of data points in input file
   =============================================================================================== */

void INOUT::countPoints(std::ifstream & inFile, int *ndim, int *npoint)
{
	std::string line, line0;

	*npoint = 0;
	int nval=0;
	double val;

/* 	==========================================================================================
	Go over each line in the file
   	========================================================================================== */

	while (getline(inFile, line0)) // until reach the end of file 
	{
/* 		==================================================================================
		Only keep lines that do not contain //
   		================================================================================== */

		if (line0.substr(0,2) != "//")
		{
			*npoint = *npoint + 1;

/* 			==========================================================================
			Count number of records on this line: transform the line into a stringstream
			and go over the number of record
   			========================================================================== */

			line = replaceAll(line0, ",", " ");
			std::istringstream iss(line);
			nval = 0;
			while ( iss >> val)
			{
				nval = nval + 1;
			}

		}
	}
	*ndim = nval;

}
	
/* ===============================================================================================
   Read points from input file
	- coordinates x,y (and maybe z)
   =============================================================================================== */

void INOUT::readPoints(std::ifstream & inFile, int ndim, int npoint, double *coord)
{
	std::string line0, line;

	int nat=0;
	double val;

	while (getline(inFile, line0))
	{
		if (line0.substr(0,2) != "//")
		{
			line = replaceAll(line0, ",", " ");
			std::istringstream iss(line);
			int ncoord = 0;
			while (iss >> val)
			{
				coord[ndim*nat + ncoord] = val;
				ncoord++;
			}
			nat++;
		}
	}
}
  /* ===============================================================================================
   Write result file
   =============================================================================================== */

  void INOUT::write(std::string outfile, int ndim, int npoint, double *Coord) 
  {

	std::string out = outfile+"_r.crd";
	std::ofstream output;
	output.open(out);

	for(int i = 0; i < npoint; i++) {
		for(int k = 0; k < ndim; k++) {
			output << Coord[k+i*ndim] << " ";
		}
		output << std::endl;
	}

	output.close();
  }
} //end namespace corera
#endif
