/* ===============================================================================================
   Author:  Yossi Bokor
   Date:    February 2 2022
   Version: 1
   =============================================================================================== */

/* ===============================================================================================
   Includes
   =============================================================================================== */

#include "Clustering.h"

/* ===============================================================================================
   Main program
   =============================================================================================== */

int main(int argc, char **argv)
{

/*	==========================================================================================
	Show usage if needed
	========================================================================================== */

	if( argc < 2 )
	{
		usage(argv);
		return -1;
	}

	std::string input = argv[1];
	if( input == "-h" || input == "-help" )
	{
		usage(argv);
		return -1;
	}

/*	==========================================================================================
	Read in all inputs (some values may have been preset)
	========================================================================================== */

	std::string INfile;
	
	int nclus;
	int nobj;

        if (!parse_args(argc, argv, &INfile, &nobj, &nclus)) return 1;

/*	==========================================================================================
	Read in the distances from input file
	========================================================================================== */

	vector<vector<double>> X;

   	ifstream in(INfile); //open the file which has the matrix of coordinates
   	string line;
   	double val;
	/* Create a vector of vectors to store the points in, note that we are treating each line as a point, rather than each column */
   	while (getline(in, line)){ 
       istringstream iss(line);
       std::vector<double> x_i;
       while (iss >> val){
		   cout << "val is " << val << endl;
          x_i.push_back(val);
       };
       X.push_back(x_i);
   	};
	
	std::cout << "The " << nobj << " points are:" << std::endl;
	for (int i=0; i < nobj; i++){
		std::cout << " point " << i << " which has length " << X[i].size() << std::endl;
		for (auto j: X[i]){
			std::cout << j << " ";
		};
		std::cout << std::endl;
	};
	kMeans clus(X, nclus); //initialise the clustering with random seed centroids

	clus.Cluster(); //perform the clustering until there are no changes in the assignments of points to clusters

	for (int i=0; i < nclus; i++){
		cout << "Cluster " << i << " has centroid (" << endl;
		int dim = clus.centroids[i].coords.size();
		for (int j = 0; j < dim-1; j++){
			cout << clus.centroids[i].coords[j]<<", " ; //print coordinate wise
		};
		cout << clus.centroids[i].coords[dim-1] << ")." << endl; //print last coordinate sperately so that we have a nice format
	};

	return 0;

};

/* ===============================================================================================
   Usage
   =============================================================================================== */

static void usage(char** argv)
{
    std::cout << "\n\n" << "     " << 
	"================================================================================================" <<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=                                       Clustering                                             ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     Usage is:                                                                                ="<<std::endl;
    std::cout << "     " << "=          Clustering -i INFILE -n NOBJ -k NCLUS                                               ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     where:                                                                                   ="<<std::endl;
    std::cout << "     " << "=                 -i  INFILE --> Input file (Distances; ascii or csv file with 1 point / line) ="<<std::endl;
    std::cout << "     " << "=                 -n  NOBJ   --> number of objects                                             ="<<std::endl;
    std::cout << "     " << "=                 -k  NCLUS  --> number of clusters to use                                     ="<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "\n\n" <<std::endl;
};

/* ===============================================================================================
   Parse Argument from command line:

   =============================================================================================== */

bool parse_args(int argc, char **argv, std::string *infile, int *nobj, int *nclus)
{
//
// Make sure we have at least two parameters....
//
	std::string param;
	if (argc == 1)
	{
		return false;
	}
	else
	{
		for (int i = 1; i < argc - 1; i = i + 2)
		{
			param = argv[i];

			if (param == "-i") {
				*infile = argv[i + 1];
			}
			if (param == "-n") {
				*nobj = std::atoi(argv[i + 1]);
			}
			if (param == "-k") {
				*nclus = std::atoi(argv[i + 1]);
			}
		}
  	}
	return true;
};


