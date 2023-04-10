/*
   PersistenceDiagram.h
    
    Calculate the 2-Wasserstein distance between two finite stes of points in the plane. Note, in this case we assume there are no points at infinity.
    
    Author: Yossi Bokor
    Date: January 28, 2022
    Version 1
*/
#ifndef _PERSISTENCEDIAGRAM_H_
#define _PERSISTENCEDIAGRAM_H_

#include <cmath>
#include "Hungarian.h"
//#include <hera/wasserstein.h> 



namespace correa {
	//template<typename Real>
	class PersistencePoint{
		
		private:
			double birth;
			double death;
		
		public:
			double Birth(){
				return birth;
			}

			double Death(){
				return death;
			}

			double Lifetime(){
				return death - birth;
			}

	};

	//template<typename Real>
	class PersistenceDiagram{
		using Point = std::vector<double>;//PersistencePoint;

		private:
			

		public:

			std::vector<Point> points;

			int np;

			PersistenceDiagram(){
				np = 0;
			}
			
			PersistenceDiagram(vector<Point> points_) :
				points (points_) {
					np = points.size();
			}

			void addPoint(Point pt);

			void printPD() {

			};

		
	};


	
	/* Constructor */
	

	/* Add a point to the persistence diagram*/
	void PersistenceDiagram::addPoint(vector<double> pt){
		points.push_back(pt);
		np += 1;
	};


	/*Wasserstein distance between two persistence diagrams */
	/*double WassersteinDistance(PersistenceDiagram D1, PersistenceDiagram D2){
		vector<vector<double>> points1 = D1.points;
		int n1 = D1.np;
		
		vector<vector<double>> points2 = D2.points;
		int n2 = D2.np;
		
		for(int i=0; i<n1; i++){ //For each point in D1, add the projection to the diagnoal to D2
			vector<double> diag((points1[i][0]+points1[i][1])/2, (points1[i][0]+points1[i][1])/2);
			points2.push_back(diag);
		};
		for(int i=0; i<n2; i++){ //For each point in D2, add the projection to the diagnoal to D1
			vector<double> diag((points2[i][0]+points2[i][1])/2, (points2[i][0]+points2[i][1])/2);
			points2.push_back(diag);
		};
		
		vector<vector<double>> CostMatrix;
		
		for(int i=0; i < n1+n2; i++){
			vector<double> i_cost;
			for(int j=0; j<n1+n2; j++){
				if(i > n1 && j > n2+1){
					i_cost.push_back(0);//If both points are on the diagonal then the cost is 0. std::numeric_limits<double>::max());
				} else {
					double cost;
					cost = pow(points1[i][0]-points2[j][0],2)+pow(points1[i][1]-points2[j][1],2); //As we are taking the p=q=2 Wasserstein distance, we don't need to take the square root here.
					i_cost.push_back(cost);
				};
			};
			CostMatrix.push_back(i_cost);
		};
		vector<int> Assignment;
		for(int i=0; i <n1+n2; i++){
			Assignment.push_back(-1); //Assign -1 to each point
		};

		double cost;
		D_HungarianAlgorithm Hungarian;
		cost = sqrt(Hungarian.Solve(CostMatrix, Assignment)); //Hungarian algorithm
		
		return cost;
	};*/


	/*void PersistenceDiagram::printPD(){
		for (int i = 0; i < np; i++){
			cout<< "\n(" << points[i][0] << ", " << points[i][1] << ")";
		};
		cout << "." << endl;
	};*/

	/*
	template<typename Real>
	std::ostream& operator<<(std::ostream& out, const PersistenceDiagram<Real>& PD){
	*/
} //end namespace correa
#endif
