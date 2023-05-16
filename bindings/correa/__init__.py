##
#@file Correa 

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../build', '_correa'))
import _correa

from matplotlib import pyplot as plt


def create_polygon(path_file : str):
	return _correa.PyPolygon(path_file)

def create_polygon_focal_point(polygon_path : str, focal_point):
	return _correa.PyPolygon(polygon_path, focal_point)


def print_polygon(p) :
	_correa.print_polygon(p)

def plot_polygon(p : _correa.PyPolygon):
	verts = p.vertices()
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
	
def polygon_properties(p : _correa.PyPolygon):
	print_polygon(p)

def compare_polygons(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon, q=2, verbose=False):
	return _correa.compare_polygons(poly1, poly2, q, verbose)

def curv_ot_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	return _correa.curv_ot_distance(poly1, poly2)

def frechet_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	return _correa.frechet_distance(poly1, poly2)

def max_ellipse_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	return _correa.max_ellipse_distance(poly1 , poly2)

def min_ellipse_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	return _correa.min_ellipse_distance(poly1, poly2)

def lsq_ellipse_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	return _correa.lsq_ellipse_distance(poly1, poly2)

def wasserstein_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon, q=2):
	return _correa.wasserstein_distance(poly1, poly2, q)

def willmore_distance(poly1 : _correa.PyPolygon, poly2 : _correa.PyPolygon):
	return _correa.willmore_distance(poly1, poly2)



