# Correa 
[![Build bindings](https://github.com/yossibokor/correa/actions/workflows/bindings.yml/badge.svg?branch=master)](https://github.com/yossibokor/correa/actions/workflows/bindings.yml)
[![Executables](https://github.com/yossibokor/correa/actions/workflows/executables.yml/badge.svg?branch=master)](https://github.com/yossibokor/correa/actions/workflows/executables.yml)
[![Documentation](https://github.com/yossibokor/correa/actions/workflows/documentation.yml/badge.svg?branch=master&event=push)](https://github.com/yossibokor/correa/actions/workflows/documentation.yml)
[![DOI](https://zenodo.org/badge/609015232.svg)](https://zenodo.org/badge/latestdoi/609015232)

## UNDER DEVELOPMENT: This is a public beta. 

Correa is a C++ library, you can choose to build Python bindings and/or executables. It is a package for analysing simple closed curves in $\mathbb{R}^2$ using a variety of tools from geometry and topology. In particular, it calculates the persistent homology of the radial growth function from a focal point in the interior of the curve. This focal point can be determined automatically or be specified by the user. As the curves are simple and closed, we avoid any infinite life persistence points and retain information about the *size* of the curve by pairing the birth of the essential $0$-cycle with the birth of the essential $1$-cycle in the $0^{th}$ persistent diagram. We assume the curves are represented by ordered lists of points, with a straight line connecting any two consecutive points. We also assume that the first point is the last point, and so this should not be repeated in the file.

If you use this software in your research, please cite it as:
Bokor Bleile, Y., Koehl, P., 'Correa', https://zenodo.org/badge/latestdoi/609015232, 2023.

### License

Correa is released under the [BSD 3-Clause License](LICENSE.md), which you should have received a copy of.


## Installation

First clone the [repository](https://github.com/yossibokor/correa) with recursive submodules, and ensure you have [CMake](https://cmake.org/), [BLAS](https://www.netlib.org/blas/) and [LAPACK](https://www.netlib.org/lapack/) installed.

Then, use [cmake](https://cmake.org/) from the directory to build Correa with the options you would like. The options are
- `CORREA_BUILD_EXECUTABLES` (default OFF)
- `CORREA_BUILD_PYTHON_BINDINGS` (default ON)


To build the Python bindings:
```bash
git clone --recurse-submodules https://github.com/yossibokor/correa.git
cd correa
cmake -S . -B build
cd build
make
```

and then ether launch python from `.../correa/build` or add `.../correa/build` to your `PYTHONPATH`.


## Usage
Once installed, you can load import a polygon from file using either `correa.create_polygon` or `correa.create_polygon_focal_point`. To create a polygon and use the center of mass of the vertices as the focal point, use `correa.create_polygon`:
```python
poly = correa.create_polygon(".../correa/exmaples/contour1.csv")
```
and if you want to specify the focal point, use `correa.reate_polygon_focal_point`:
```python
poly = correa.create_polygon_focal_point(".../correa/exmaples/contour1.csv", focal_point)
```
where `focal_point` is either the path a file containing the two coordinates of the focal point, or a list containing the two coordinates. The focal point (either automatic or user specified) forms the base point for the persistent homology computations. The persistence diagram for a polygon is the $0^{th}$ persistent homology of the distance from the focal point function, where the birth of the essential $0$-cycle is matched with the birth of the essential $1$-cycle. This is possible as we assumed the curves the polygons represent are simple closed curves, hence have a unique essential $0$-cycle and unique essential $1$-cycle.

You can then print the properties of such a polygon (including the persistence diagram via the radial distance function from the focal point) using `print_polygon:
```python
correa.print_polygon(poly)
```

For a pair of polygons `poly1` and `poly2`, we can compare them using a varitey of metrics:
- Wassterstein distance between the persistence diagrams (`correa.wasserstein_distance`)
- Frechet distance (`correa.frechet_distance`)
- max ellipse distance (`correa.max_ellipse_distance`)
- min ellipse distance (`correa.min_ellipse_distance`)
- least square ellipse distance (`correa.lsq_ellipse_distance`)
- Willmore distance (`correa.willmore_distance`)

For more information about these distances, see [distances](DISTANCES.md)
  
which can be done individually or all at the same time. To calculate all the distances at once, use `compare_polygons`, otherwise use the corresponding function. 

```python
correa.compare_polygons(poly1, poly2, q=2, verbose=False)
```

Using `correa.compare_polygons` returns the distances 
```python
dWasserstein, dFrechet, dMax, dMin, dLSQ, dWillmore, dCurvOT
```


## UNDER DEVELOPMENT: This is a public beta. 
Please use [issues](https://github.com/yossibokor/correa/issues) to report any bugs.
