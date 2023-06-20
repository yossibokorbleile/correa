"""Correa 

Documentation for the Python functions in Correa.
"""

import _correa

from matplotlib import pyplot as plt

def create_polygon(poly_path : str):
	"""! 	Read a polygon from file. This will recenter the polygon to the center of mass of the vertices, so be careful.
			@param poly_path 	path to the file you want to read in.
			@return 			PyPolygon object.	
	"""
	return _correa.PyPolygon(poly_path)

def create_polygon_focal_point(poly_path : str, focal_point):
	"""! 	Read a polygon from file, and specify a focal point.
			@param poly_path 	path to the file containing the polygon.
			@param focal_point		either a path to the file containing the focal point, or a list with the coordinates.
			@return 	PyPolygon object
	"""
	return _correa.PyPolygon(poly_path, focal_point)


def print_polygon(poly) :
	"""!	Print information about the polygon.
			@param poly		the polygon you want information about.
	"""
	_correa.print_polygon(poly)

def plot_polygon(poly : _correa.PyPolygon):
	"""!	Plot a polygon.
			@param poly the polygon you want to plot.
			@return 	a matplotlib.pyplot figure.
	"""
	verts = poly.vertices()
	x = []
	y = []
	for i in range(len(verts)):
		x.append(verts[i][0])
		y.append(verts[i][1])
	x.append(verts[0][0])
	y.append(verts[0][1])
	fig = plt.figure(1, figsize=(5,5), dpi=90)  
	ax = fig.add_subplot(111)
	ax.plot(x, y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
	ax.set_title('Polygon')
	fig.show()
	return fig

def compare_polygons(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon, q=2, verbose=False):
	"""!	Given two polygons (poly1, poly2), compare them.
			@param poly1 	first polygon to compare.
			@param poly2	second polygon to compare.
			@param q		q for the Wasserstein distance,.
			@param verbose	setting for verbose output.
			@return 		dWasserstein, dFrechet, dMax, dMin, dLSQ, dWillmore, dCurvOT.
	"""
	return _correa.compare_polygons(poly1, poly2, q, verbose)

def curv_ot_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	"""!	Given two polygons (poly1, poly2), calculate the curv_ot distance between them.
			@param poly1 	first polygon to compare.
			@param poly2	second polygon to compare.
			@return 		distance.
	"""
	return _correa.curv_ot_distance(poly1, poly2)

def frechet_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	"""!	Given two polygons (poly1, poly2), calculate the frechet distance between them.
			@param poly1 	first polygon to compare.
			@param poly2	second polygon to compare.
			@return 		distance.
	"""
	return _correa.frechet_distance(poly1, poly2)

def max_ellipse_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	"""!	Given two polygons (poly1, poly2), calculate the max ellipse distance between them.
			@param poly1 	first polygon to compare.
			@param poly2	second polygon to compare.
			@return 		distance.
	"""
	return _correa.max_ellipse_distance(poly1 , poly2)

def min_ellipse_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	"""!	Given two polygons (poly1, poly2), calculate the min ellipse distance between them.
			@param poly1 	first polygon to compare.
			@param poly2	second polygon to compare.
			@return 		distance.
	"""
	return _correa.min_ellipse_distance(poly1, poly2)

def lsq_ellipse_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	"""!	Given two polygons (poly1, poly2), calculate the least square ellipse distance between them.
			@param poly1 	first polygon to compare.
			@param poly2	second polygon to compare.
			@return 		distance.
	"""
	return _correa.lsq_ellipse_distance(poly1, poly2)

def wasserstein_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon, q=2):
	"""!	Given two polygons (poly1, poly2), calculate the wasserstein distance between them via their persistence diagrams.
			@param poly1 	first polygon to compare.
			@param poly2	second polygon to compare.
			@return 		distance.
	"""
	return _correa.wasserstein_distance(poly1, poly2, q)

def willmore_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	"""!	Given two polygons (poly1, poly2), calculate the willmore distance between them.
			@param poly1 	first polygon to compare.
			@param poly2	second polygon to compare.
			@return 		distance.
	"""
	return _correa.willmore_distance(poly1, poly2)
