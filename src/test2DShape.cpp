#include <iostream>
#include <tuple>
#include <vector>
#include "Polygon.h"
#include "Ellipse.h"
int main(){
	std::cout << "Let us create a tuple" << std::endl;
	std::tuple<int, int> test{0,1};
	std::vector<std::tuple<int, double, int> > unsorted_nodes;
	std::cout << "Tuple is " << std::get<0>(test) << " and " << std::get<1>(test) << std::endl;
	correa::Polygon polygon();
	std::cout << "Hello world!" << std::endl;
	correa::Ellipse ellipse;
	return 1;
}