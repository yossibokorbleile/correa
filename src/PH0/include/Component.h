/*
 Components.hpp

 Create a filtration of a polygon

 Author: Yossi Bokor
 Date: January 28, 2022
 Version 1
 */



#ifndef _COMPONENT_H
#define _COMPONENT_H


#include <iostream>
#include <tuple>

/* Comp class */
class Comp{

public:

	int num_points;

	int start;

	int end;

	double birth;

	double death;

	int parent;

	bool Alive;

	Comp();

	Comp(int i, double height, int n_points);

	Comp(std::pair<int, double> x, int n_points);

	int Parent();

	int AddNode(std::pair<int, double> x);
};



/* Constructor from a node */
Comp::Comp(std::pair<int, double> x, int n_points){
	parent = x.first;
	birth = x.second;
	start = x.first;
	end = x.first;
	death = x.second;
	num_points = n_points;
	Alive =true;
};

//Constructor from an integer and a double
Comp::Comp(int i, double height, int n_points){
	parent = i;
	birth = height;
	start = i;
	end = i;
	num_points = n_points;
	Alive = true;
};

/* Find the parent */
int Comp::Parent(){
	return parent;
};

/* Add a node to a component */
int Comp::AddNode(std::pair<int, double> x){
	if (x.first != (start % num_points)-1 && x.first != (end % num_points)+1){
		std::cout << "The node you are trying is not adjacent to this component." << std::endl;
		return 0;
	} else if (x.first == start-1){
		start = x.first;
		death = x.second;
		return 1;
	}else if (x.first == end+1){
		end = x.first;
		death = x.second;
		return 1;
	} else {
		return 0;
	};
};

#endif
