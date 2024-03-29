/*!
* @file PH0.h
* @brief functions to compute dimension 0 persistent homology of polygon.
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
*/

#ifndef _PH0_H
#define _PH0_H


#include <iostream>
#include <vector>
#include <tuple>
#include "Vertex.h"
#include "Polygon.h"
//#include "Component.h"
#include <utility>
using namespace std;

/* Lets do the the persistent homology in dimension 0 now */
namespace correa {
/* Class to compute the persistent homology in dimension 0. */
class PH0{
	using PersistenceDiagram = std::vector<std::pair<double,double>>;

	private:
		vector<std::tuple<int, double, int> > unsorted_nodes;
		//vector<Comp> comps;
		PersistenceDiagram pd;

	public:

		int n_pts; // Number of nodes in the polygon

		PH0(vector<pair<int, double>> x); // Construct from a vector pairs, where the pair is of the form <int index, double height>. 
		PH0(vector<Vertex> points);
		void AddNode(pair<int, double>& y);
		void Persistence();
		PersistenceDiagram persistence_diagram();
		void printPD();

};

/* Constructor from a vector of nodes */
PH0::PH0(vector<pair<int, double> > x){
	n_pts = x.size();
	for (int i = 0; i < n_pts; i++){
		tuple<int, double, int> pt;
		get<0>(pt) = x[i].first;
		get<1>(pt) = x[i].second;
		get<2>(pt) = x[i].first;
		unsorted_nodes.push_back(pt);
	};
};

/* Constructor from a vector of Vertex */
PH0::PH0(vector<Vertex> points){
	n_pts = points.size();
	for (int i = 0; i < n_pts; i++){
		Vertex vert = points[i];
		double height;
		height = vert.height();
		tuple<int, double, int> vertex_tuple(i, height, i);
		unsorted_nodes.push_back(vertex_tuple);
	};
};

/* Add nodes to the filtation */
void PH0::AddNode(pair<int, double>& y){
	tuple<int, double, int> pt;
	get<0>(pt) = y.first;
	get<1>(pt) = y.second;
	get<2>(pt) = y.first;
	unsorted_nodes.push_back(pt);
};

/* Compare the heights of two nodes */
bool HeightComparison(tuple<int, double>& x, tuple<int, double>& y){
	//if (get<1>(x) == get<1>(y)) {
	//	return true;
	//} else {
		return (get<1>(x)< get<1>(y));
	//};
};


/* Obtain the components and their birth/death times */
void PH0::Persistence(){
	vector<tuple<int, double> > sorted_nodes;
	for (int i = 0; i < n_pts; i++){
		tuple<int, double> x;
		get<0>(x) = get<0>(unsorted_nodes[i]);
		get<1>(x) = get<1>(unsorted_nodes[i]);
		sorted_nodes.push_back(x);
	};
	sort(sorted_nodes.begin(), sorted_nodes.end(), HeightComparison);
	for (int i = 0; i < n_pts; i ++){
		int looking_at;
		looking_at = get<0>(sorted_nodes[i]);
		if (looking_at == 0){ //If the point we are looking at is the first point in the polygon, this requires special treament
			if ((get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[1])) && (get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[n_pts-1]))){ //If point 0 appears before point 1 and the last point, then we make it it's own parent.
				get<2>(unsorted_nodes[looking_at]) = looking_at;
			} else if ( (get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[1])) && (get<1>(sorted_nodes[i]) > get<1>(unsorted_nodes[n_pts-1]))){ //If point 0 appears before point 1 and after the last point, then we make the last points parent the parent of point 0.
				get<2>(unsorted_nodes[looking_at]) = get<2>(unsorted_nodes[n_pts-1]);
			} else if ( (get<1>(unsorted_nodes[looking_at]) > get<1>(unsorted_nodes[1])) && (get<1>(sorted_nodes[i]) < get<1>(unsorted_nodes[n_pts-1]))){ //If point 0 appears after point 1 and before the last point, then we make the point 1s parent to as the parent of point 0.
				get<2>(unsorted_nodes[looking_at]) = get<2>(unsorted_nodes[1]);
			} else if ( (get<1>(unsorted_nodes[looking_at]) >= get<1>(unsorted_nodes[1])) && (get<1>(sorted_nodes[i]) >= get<1>(unsorted_nodes[n_pts-1]))){ //If point 0 appears after point 1 and the last point, then we need to follow the elder rule.
				if (get<1>(unsorted_nodes[get<2>(unsorted_nodes[n_pts-1])]) <= get<1>(unsorted_nodes[get<2>(unsorted_nodes[1])])){//We make the point npts-1 parent to as the parent of point 0 and all the points connected to point 1, and add a point to the persistence diagram.
					int parent, old_parent;
					parent = get<2>(unsorted_nodes[n_pts-1]);
					old_parent = get<2>(unsorted_nodes[1]);
					get<2>(unsorted_nodes[looking_at]) = parent;
					std::pair<double,double> pt{get<1>(unsorted_nodes[get<2>(unsorted_nodes[old_parent])]), get<1>(unsorted_nodes[looking_at])}; //Add a point to the persistence diagram.
					pd.push_back(pt);
					for (int j =0; j < unsorted_nodes.size(); j++){
						if (get<2>(unsorted_nodes[j]) == old_parent){
							get<2>(unsorted_nodes[j]) = parent;
						};
					};
				} else if (get<1>(unsorted_nodes[get<2>(unsorted_nodes[n_pts-1])]) > get<1>(unsorted_nodes[get<2>(unsorted_nodes[1])])){//We make the point 1 parent to as the parent of point 0 and all the points connected to point n_pts-1, and add a point to the persistence diagram.
					int parent, old_parent;
					parent = get<2>(unsorted_nodes[1]);
					old_parent = get<2>(unsorted_nodes[n_pts-1]);
					get<2>(unsorted_nodes[looking_at]) = parent;
					std::pair<double,double> pt{get<1>(unsorted_nodes[old_parent]), get<1>(unsorted_nodes[looking_at])}; //Add a point to the persistence diagram.
					pd.push_back(pt);
					for (int j =0; j < unsorted_nodes.size(); j++){
						if (get<2>(unsorted_nodes[j]) == old_parent){
							get<2>(unsorted_nodes[j]) = parent;
						};
					};
				};
			};
		} else if (looking_at == n_pts-1){ //If the point we are looking at is the last point in the polygon, this requires special treament
			if ((get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[0])) && (get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[n_pts-2]))){ //If point looking_at appears before point looking_at+1 and point looking_at-1, then we make it it's own parent.
				get<2>(unsorted_nodes[looking_at]) = looking_at;
			} else if ( (get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[0])) && (get<1>(sorted_nodes[i]) > get<1>(unsorted_nodes[n_pts-2]))){ //If point looking_at appears before point looking_at+1 and after point looking_at-1, then we make it it's own parent.
				get<2>(unsorted_nodes[looking_at]) = get<2>(unsorted_nodes[n_pts-2]);
			} else if ( (get<1>(unsorted_nodes[looking_at]) > get<1>(unsorted_nodes[0])) && (get<1>(sorted_nodes[i]) < get<1>(unsorted_nodes[n_pts-2]))){ //If point looking_at appears after point looking_at+1 and before point looking_at-1, then we make it it's own parent.
				get<2>(unsorted_nodes[looking_at]) = get<2>(unsorted_nodes[0]);
			} else if ( (get<1>(unsorted_nodes[looking_at]) >= get<1>(unsorted_nodes[0])) && (get<1>(sorted_nodes[i]) >= get<1>(unsorted_nodes[n_pts-2]))){ //If point looking_at appears after point looking_at+1 and point looking_at-1, then we follow the elder rule.
				if (get<1>(unsorted_nodes[get<2>(unsorted_nodes[looking_at-1])]) <= get<1>(unsorted_nodes[get<2>(unsorted_nodes[0])])){//Check which neighbouring point has the older parent, in this case it is looking_at-1, and so we take its parents.
					int parent;
					parent = get<2>(unsorted_nodes[looking_at-1]);
					get<2>(unsorted_nodes[looking_at]) = parent;
					std::pair<double,double> pt{get<1>(unsorted_nodes[get<2>(unsorted_nodes[0])]), get<1>(unsorted_nodes[looking_at])}; //Add a point to the persistence diagram.
					pd.push_back(pt);
					for (int j = 0; j < unsorted_nodes.size(); j++){
						if (get<2>(unsorted_nodes[j]) ==get<2>(unsorted_nodes[0])){
							get<2>(unsorted_nodes[j]) = parent;
						};
					};
				} else if (get<1>(unsorted_nodes[get<2>(unsorted_nodes[looking_at-1])]) > get<1>(unsorted_nodes[get<2>(unsorted_nodes[0])])){//Check which neighbouring point has the older parent, in this case it is looking_at+1, and so we take its parents.
					int parent, old_parent;
					parent = get<2>(unsorted_nodes[0]);
					old_parent = get<2>(unsorted_nodes[looking_at-1]);
					get<2>(unsorted_nodes[looking_at]) = parent;
					std::pair<double,double> pt{get<1>(unsorted_nodes[old_parent]), get<1>(unsorted_nodes[looking_at])}; //Add a point to the persistence diagram.
					pd.push_back(pt);
					for (int j =0; j < unsorted_nodes.size(); j++){
						if (get<2>(unsorted_nodes[j]) == old_parent){
							get<2>(unsorted_nodes[j]) = parent;
						};
					};
				};
			};
		} else { //Any other point does not require special treatment
			if ((get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[looking_at+1])) && (get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[looking_at-1]))){ //If point looking_at appears before point 1 and the last point, then we make it it's own parent.
				get<2>(unsorted_nodes[looking_at]) = looking_at;
			} else if ( (get<1>(unsorted_nodes[looking_at]) < get<1>(unsorted_nodes[looking_at+1])) && (get<1>(sorted_nodes[i]) > get<1>(unsorted_nodes[looking_at-1]))){ //If point looking_at appears before point looking_at+1 and after point looking_at-1, then we make point looking_at+1s parent the parent of point looking_at.
				get<2>(unsorted_nodes[looking_at]) = get<2>(unsorted_nodes[looking_at-1]);
			} else if ( (get<1>(unsorted_nodes[looking_at]) > get<1>(unsorted_nodes[looking_at+1])) && (get<1>(sorted_nodes[i]) < get<1>(unsorted_nodes[looking_at-1]))){ //If point looking_at appears after point looking_at+1 and before point looking_at-, then we make the point looking_at-1s parent to as the parent of point looking_at.
				get<2>(unsorted_nodes[looking_at]) = get<2>(unsorted_nodes[looking_at+1]);
			} else if ( (get<1>(unsorted_nodes[looking_at]) >= get<1>(unsorted_nodes[looking_at+1])) && (get<1>(sorted_nodes[i]) >= get<1>(unsorted_nodes[looking_at-1]))){ //If point looking_at appears after point looking_at+1 and point looking_at-1, then we check the elder rule.
				if (get<1>(unsorted_nodes[get<2>(unsorted_nodes[looking_at-1])]) <= get<1>(unsorted_nodes[get<2>(unsorted_nodes[looking_at+1])])) {//Check which neighbouring point has the older parent, in this case it is looking_at-1, and so we take its parents.
					int parent, old_parent;
					parent = get<2>(unsorted_nodes[looking_at-1]);
					old_parent = get<2>(unsorted_nodes[looking_at+1]);
					get<2>(unsorted_nodes[looking_at]) = parent;
					std::pair<double,double> pt{get<1>(unsorted_nodes[get<2>(unsorted_nodes[looking_at+1])]), get<1>(unsorted_nodes[looking_at])}; //Add a point to the persistence diagram.
					pd.push_back(pt);
					for (int j = 0; j < unsorted_nodes.size(); j++){
						if (get<2>(unsorted_nodes[j]) == old_parent){
							get<2>(unsorted_nodes[j]) = parent;
						};
					};
				} else if (get<1>(unsorted_nodes[get<2>(unsorted_nodes[looking_at-1])]) > get<1>(unsorted_nodes[get<2>(unsorted_nodes[looking_at+1])])) {//Check which neighbouring point has the older parent, in this case it is looking_at+1, and so we take its parents.
					int parent, old_parent;
					parent = get<2>(unsorted_nodes[looking_at+1]);
					old_parent = get<2>(unsorted_nodes[looking_at-1]);
					get<2>(unsorted_nodes[looking_at]) = parent;
					std::pair<double,double> pt{get<1>(unsorted_nodes[get<2>(unsorted_nodes[looking_at-1])]), get<1>(unsorted_nodes[looking_at])}; //Add a point to the persistence diagram.
					pd.push_back(pt);
					for (int j =0; j < unsorted_nodes.size(); j++){
						if (get<2>(unsorted_nodes[j]) == old_parent){
							get<2>(unsorted_nodes[j]) = parent;
						};
					};
				};
			};
		};
	};
};	

// Return the persistence diagram.
vector<std::pair<double, double>> PH0::persistence_diagram() {
	return pd;
};

// Print the persistence diagram
void PH0::printPD() {
	for (int i = 0; i < pd.size(); i++) {
		std::cout << "(" << get<0>(pd[i]) << ", " << get<1>(pd[i]) << ")" << std::endl;		
		}	
};

} //end namespace correa
#endif
