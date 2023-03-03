/*
    kMeans.h
    
    kMeans clustering, following https://reasonabledeviations.com/2019/10/02/k-means-in-cpp/.
    
    Author: Yossi Bokor
    Date: January 14, 2022
    Version 0.1
*/

#ifndef _KMEANS_H
#define _KMEANS_H

#include <ctime>     // for a random seed
#include <fstream>   // for file-reading
#include <iostream>  // for file-reading
#include <sstream>   // for file-reading
#include <vector>
#include <math.h>
#include <algorithm>
#include <set>

using namespace std;

/*Define a class to store each data point in*/
class Point{
	
	public:
		
		vector<double> coords;
		
		int clust;
		
		double minDist;
		
		Point(vector<double> coords);

		/*Point(double coords);*/
		
};

/*Constructor*/
Point::Point(vector<double> Coords){
	coords = Coords;
	clust = -1;
};

/*Constructor from array */
/*Point::Point(double Coords){
	for (double i: Coords){
		coords.push_back(i);
	};
	clust = -1;
};
*/

/*Define a class for the centroids*/
class Centroid{
	
	public:
		
		vector<double> coords;
		
		Centroid(vector<Point>& points);
		
		Centroid(vector<Point>& points, vector<int> assigned);
		
		Centroid(Point& pt);
		
		void PrintCentroid();
};

/*Constructor* from a vector of points*/
Centroid::Centroid(vector<Point>& points){
	int n_points;
	int n_dim;
	n_points = points.size();
	n_dim = points[0].coords.size();
	
	for (int i=0; i < n_dim; i++){
		double i_coord;
		i_coord = 0;
		for(int j=0; j < n_points; j++){
			i_coord += points[j].coords[i];
		};
		coords.push_back(i_coord/n_points);
	};
};

/*Constructor from a vector of points and the ones assigned*/
Centroid::Centroid(vector<Point>& points, vector<int> assigned){
	int n_points;
	int n_dim;
	n_points = assigned.size();
	n_dim = points[0].coords.size();
	
	for (int i=0; i < n_dim; i++){
		double i_coord;
		i_coord = 0;
		for(int j=0; j < n_points; j++){
			i_coord += points[assigned[j]].coords[i];
		};
		coords.push_back(i_coord/n_points);
	};
};

/*Cosntructor from a coordinate vector*/
Centroid::Centroid(Point& pt){
	coords = pt.coords;
};

void Centroid::PrintCentroid(){
	int dim = coords.size();
	cout << "(";
	for (int i = 0; i < dim-1; i++){
		cout << coords[i] <<", ";
	};
	cout << coords[dim] << ")" << endl;
};
/* Distance between a point and a centroid */
double distance(Point p, Centroid c){
	int n_dim;
	n_dim = c.coords.size();
	if(n_dim != p.coords.size()){
		cout << "The point and centroid are in different dimensions: " << n_dim << " vs " << p.coords.size() << endl;
		assert(n_dim == p.coords.size());
	};
	double dist_sq;
	dist_sq = 0;
	for (int i=0; i < n_dim; i++){
		dist_sq += pow(p.coords[i] - c.coords[i], 2);
	};
	
	return sqrt(dist_sq);
};

/*Class to do kMeans in*/
class kMeans {
	
	public:
	
		int n_pts;
	
		int k;
		
		bool Done;
	
		vector<Point> points;
	
		vector<Centroid> centroids;
	
		vector<int> assignment;
	
		kMeans(vector<vector<double>> points, int n_clus);
		
		void updateCentroids();
		
		void updateAssignment();
		
		void Cluster();
	
};

/*kMeans Constructor*/
kMeans::kMeans(vector<vector<double>> pts, int n_clus){
	Done = false;
	n_pts = pts.size();
	for (int i=0; i < n_pts; i++){ //Add the points.
		points.push_back(Point(pts[i]));
	};
	k = n_clus;
	cout << "We have "<< n_pts << " points." << endl;
	srand(time(0));
	set<int> Chosen;
	for (int i = 0; i < k; i++){ //Initialise the centroids as random points.
		int j = rand() & (n_pts - 1);
		bool good = false;
		while (good == false){
			if (Chosen.count(j) == 0){
				good = true;
			} else {
				j = rand() & (n_pts - 1);
			};
		};
		cout << "random choice is " << j << endl;
		Point& chosenRef = points[j];
		cout << "chosen point is in dimension " << points[j].coords.size() << endl;
		Centroid c(chosenRef);
		cout << "Centroid is in dimension " << c.coords.size() << endl;
		centroids.push_back(Centroid(chosenRef));
		Chosen.insert(j);
	};
	for (int i=0; i< n_pts; i++){
		vector<double> i_dists;
		for (int j=0; j < k; j++){
			i_dists.push_back(distance(points[i], centroids[j]));
		};
		int closestCentroid = std::min_element(i_dists.begin(),i_dists.end()) - i_dists.begin();
		assignment.push_back(closestCentroid);
	};
};


/*Update the centroids*/
void kMeans::updateCentroids(){
	vector<vector<int>> clusters; //for each centroid, a vector of points assigned to it
	for (int i =0; i < k; i++){
		vector<int> i_vec;
		clusters.push_back(i_vec);
	};
	for (int i = 0; i < n_pts; i++){
		cout<< "Point number " << i << " is assigned to " << assignment[i] << endl;
			int a = assignment[i];
			clusters[a].push_back(i);
	};
	vector<Point>& pointsRef = points;
	for(int i =0; i< k; i++){
		cout << "The following points are in cluster " << i << endl;
		vector<int> i_vec = clusters[i];
		for(int j=0; j < i_vec.size(); j++){
			cout << i_vec[j] << endl;
		};
	};
	for (int i = 0; i <k; i++){
		centroids[i] = Centroid(pointsRef, clusters[i]);//update the centroid
	};
};


/*Update the assignments*/
void kMeans::updateAssignment(){
	bool changed = false;
	for (int i =0; i < n_pts; i++){
		vector<double> i_dists;
		for (int j=0; j < k; j++){
			i_dists.push_back(distance(points[i], centroids[j]));
		};
		int closestCentroid = std::min_element(i_dists.begin(),i_dists.end()) - i_dists.begin();
		if (assignment[i] != closestCentroid){
			assignment[i] = closestCentroid;
			changed = true;
		};
	};
	if (changed == false){
		Done = true;
	};
};

/*Perform clusering until stablised*/
void kMeans::Cluster(){
	while(Done == false){
		updateCentroids();
		for(int i = 0; i < k; i++){
			Centroid c = centroids[i];
			c.PrintCentroid();
		};
		updateAssignment();
		cout << "we have the following assignemnts:" << endl;
		for (int i =0; i < n_pts; i++){
			cout << i << " assigned to " << assignment[i] << endl;
		};
	};
};

#endif /* _KMEANS_H */
