##
#@file Correa 

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../build', '_correa'))
import _correa

from matplotlib import pyplot as plt


def create_polygon(path_file : str):
	return _correa.PyPolygon(path_file)


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

def ComparePolygons():
	return _correa.ComparePolygons()

def AllDistances(cp, poly1, poly2, q=2, verbose=False):
    return cp.AllDistances(poly1, poly2, q, verbose)