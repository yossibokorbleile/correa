/*!
* @file PersistenceDiagram.h
* @brief define classes needed for peristence diagrams, including persistence points
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date January 2022
* @version 1
*/

#ifndef _PERSISTENCEDIAGRAM_H_
#define _PERSISTENCEDIAGRAM_H_

#include <cmath>
#include "Hungarian.h"



namespace correa {
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

	
	/* Add a point to the persistence diagram*/
	void PersistenceDiagram::addPoint(vector<double> pt){
		points.push_back(pt);
		np += 1;
	};

} //end namespace correa
#endif
