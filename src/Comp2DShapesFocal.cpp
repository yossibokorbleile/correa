/*!
* @file	Comp2DShapesFocal.cpp
* @brief compare polygons with specific focal points
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date December 2025
* @version 1.1
* @copyright BSD 3-Clause License.
*/

/*!
 * @page distances Distance Metrics for 2D Shape Comparison
 *
 * @section distances_overview Overview
 *
 * This page describes the various distance metrics implemented in the Correa shape comparison
 * library. These metrics quantify the dissimilarity between two 2D polygonal shapes from different
 * geometric and topological perspectives.
 *
 * @section distances_types Available Distance Metrics
 *
 * @subsection frechet_distance Fréchet Distance (Type 0)
 *
 * The Fréchet distance measures the similarity between two curves by considering their trajectories.
 * It can be intuitively understood as the minimum length of a leash required for a person walking
 * along one curve and their dog walking along another curve, where both can control their speed
 * but cannot go backwards.
 *
 * **Mathematical Definition:**
 * Given two curves P and Q, the Fréchet distance is:
 * \f[
 * \delta_F(P, Q) = \inf_{\alpha, \beta} \max_{t \in [0,1]} d(P(\alpha(t)), Q(\beta(t)))
 * \f]
 * where $\alpha$ and $\beta$ are continuous monotone reparameterizations.
 *
 * **Properties:**
 * - Captures both location and ordering of points along curves
 * - Sensitive to the shape of the entire trajectory
 * - Metric: satisfies triangle inequality, symmetry, and non-negativity
 *
 * **Use Cases:**
 * - Comparing shapes where the ordering of points matters
 * - Path similarity analysis
 * - Shape matching where global structure is important
 *
 * @subsection ellipse_distances Ellipse-Based Distances (Type 1)
 *
 * These metrics compare shapes based on the properties of fitted ellipses. Three types of
 * ellipse fits are computed for each polygon:
 *
 * **1. Maximum Volume Inscribed Ellipse:**
 * - The largest ellipse that fits completely inside the polygon
 * - Provides a measure of the polygon's interior capacity
 *
 * **2. Minimum Volume Inscribing Ellipse:**
 * - The smallest ellipse that completely contains the polygon
 * - Provides an outer bound on the shape
 *
 * **3. Least-Squares Ellipse:**
 * - Ellipse fitted by minimizing squared distances to polygon vertices
 * - Balances interior and exterior fit
 *
 * **Distance Computation:**
 * For each ellipse type, we compute the aspect ratio r = a/b (where a and b are semi-major
 * and semi-minor axes). The distance between two polygons is the absolute difference in
 * their aspect ratios:
 * \f[
 * d_E(P_1, P_2) = |r_1 - r_2|
 * \f]
 *
 * **Properties:**
 * - Captures overall elongation and orientation
 * - Robust to small perturbations in polygon vertices
 * - Scale-invariant when comparing aspect ratios
 *
 * **Use Cases:**
 * - Comparing overall shape elongation
 * - Orientation-independent shape classification
 * - Quality control in manufacturing (elliptical objects)
 *
 * @subsection curvature_distances Curvature-Based Distances (Type 2)
 *
 * These metrics analyze shapes based on their curvature properties, capturing local geometric
 * features and bending characteristics.
 *
 * **1. Willmore Energy:**
 * The Willmore energy measures the total squared curvature of a curve:
 * \f[
 * W = \int \kappa^2 \, ds
 * \f]
 * where $\kappa$ is the curvature and s is arc length.
 *
 * The distance between two polygons is:
 * \f[
 * d_W(P_1, P_2) = |W_1 - W_2|
 * \f]
 *
 * **Properties:**
 * - Measures the "bending energy" of the curve
 * - Minimized by circles (constant curvature)
 * - Sensitive to sharp corners and irregularities
 *
 * **2. Wasserstein Distance on Curvature Distributions:**
 * Treats the curvature distribution along each polygon as a probability measure and
 * computes the optimal transport distance between them.
 *
 * \f[
 * d_{OT}(\mu_1, \mu_2) = \inf_{\gamma \in \Gamma(\mu_1, \mu_2)} \int |x - y| \, d\gamma(x, y)
 * \f]
 *
 * **Properties:**
 * - Captures differences in local curvature patterns
 * - More informative than global curvature measures
 * - Robust metric structure (satisfies triangle inequality)
 *
 * **Use Cases:**
 * - Detecting shape deformations and irregularities
 * - Smoothness comparison
 * - Medical imaging (comparing organ boundaries)
 *
 * @subsection persistence_distance Persistence Diagram Distance (Type 3)
 *
 * This metric uses topological data analysis to compare shapes based on their persistent
 * homology features, specifically 0-dimensional persistence (connected components).
 *
 * **Persistence Diagrams:**
 * A persistence diagram is a multiset of points (b, d) in the plane, where:
 * - b (birth): the scale at which a topological feature appears
 * - d (death): the scale at which the feature disappears
 * - Persistence = d - b measures feature significance
 *
 * **2-Wasserstein Distance:**
 * The distance between two persistence diagrams $D_1$ and $D_2$ is:
 * \f[
 * W_2(D_1, D_2) = \left( \inf_{\gamma: D_1 \to D_2} \sum_{p \in D_1} \|p - \gamma(p)\|_2^2 \right)^{1/2}
 * \f]
 * where the infimum is over all bijections $\gamma$ (including the diagonal).
 *
 * **Properties:**
 * - Topologically robust: insensitive to small perturbations
 * - Captures multi-scale shape features
 * - Mathematically well-founded metric
 * - Stable under continuous deformations
 *
 * **Use Cases:**
 * - Shape classification with topological invariance
 * - Robust feature detection in noisy data
 * - Multi-scale shape analysis
 * - Pattern recognition in scientific data
 *
 * @section distances_sphericity Sphericity Measure
 *
 * All comparison programs also compute the sphericity difference:
 * \f[
 * d_S(P_1, P_2) = 4\pi \left| \frac{A_1}{L_1^2} - \frac{A_2}{L_2^2} \right|
 * \f]
 * where A is area and L is perimeter length.
 *
 * **Properties:**
 * - Measures how close a shape is to a circle (circle has sphericity = 1)
 * - Scale-invariant (dimensionless quantity)
 * - Simple and computationally efficient
 *
 * @section distances_choosing Choosing a Distance Metric
 *
 * Different metrics are appropriate for different applications:
 *
 * | Metric | Best For | Computational Cost |
 * |--------|----------|-------------------|
 * | Fréchet | Path similarity, ordered matching | Medium |
 * | Ellipse-based | Overall shape, elongation | Low |
 * | Curvature-based | Local features, smoothness | Medium |
 * | Persistence | Topological features, noisy data | High |
 * | Sphericity | Quick comparison, roundness | Very Low |
 *
 * **Recommendations:**
 * - Use **Type 4 (All)** for comprehensive analysis when computational cost is not a concern
 * - Use **Fréchet** for ordered point sequences and path matching
 * - Use **Ellipse-based** for fast, orientation-independent comparison
 * - Use **Curvature-based** when local shape details matter
 * - Use **Persistence** for robust topological comparison resistant to noise
 *
 * @section distances_implementation Implementation Notes
 *
 * All distance metrics are implemented in the following programs:
 * - Comp2DShapesFocal: User-specified focal points for alignment
 * - Comp2DShapes: Automatic centroid-based alignment
 *
 * The distance type is selected via the `-d` command-line flag.
 *
 * @see Comp2DShapesFocal.cpp
 * @see Comp2DShapes.cpp
 * @see parse_args()
 */


/* ===============================================================================================
   Includes
   =============================================================================================== */

#include "Comp2DShapesFocal.h"


/* ===============================================================================================
   Helper function to parse focal point from either file or inline coordinates
   =============================================================================================== */

/*!
 * @brief Parse focal point from either a file path or inline coordinates
 *
 * This function accepts either a file path containing focal point coordinates
 * or a string with inline coordinates in the format "x,y". If a valid file
 * exists at the path, it reads the coordinates from the file. Otherwise, it
 * parses the string as inline coordinates.
 *
 * @param focal_input String containing either a file path or inline coordinates (e.g., "698.678,816.662")
 * @param verbose If true, prints debug information during parsing (default: false)
 * @return Vector2D containing the parsed focal point coordinates
 * @throws std::exit(1) if the focal point cannot be parsed or has fewer than 2 coordinates
 */
Vector2D parse_focal_point(const std::string& focal_input, bool verbose = false) {
	std::vector<double> focal_coords;

	// Try to open as a file first
	std::ifstream infile(focal_input);
	if (infile.good()) {
		// Read from file
		std::string line;
		double val;
		while (getline(infile, line)) {
			size_t start_pos = 0;
			while ((start_pos = line.find(",", start_pos)) != std::string::npos) {
				line.replace(start_pos, 1, " ");
				start_pos += 1;
			}
			std::istringstream iss(line);
			while (iss >> val) {
				if (verbose) std::cout << "val is " << val << std::endl;
				focal_coords.push_back(val);
			}
		}
		infile.close();
	} else {
		// Parse as inline coordinates (e.g., "698.678,816.662")
		std::string coords_str = focal_input;
		size_t start_pos = 0;
		while ((start_pos = coords_str.find(",", start_pos)) != std::string::npos) {
			coords_str.replace(start_pos, 1, " ");
			start_pos += 1;
		}
		std::istringstream iss(coords_str);
		double val;
		while (iss >> val) {
			if (verbose) std::cout << "val is " << val << std::endl;
			focal_coords.push_back(val);
		}
	}

	if (focal_coords.size() < 2) {
		std::cerr << "Error: Could not parse focal point from '" << focal_input << "'" << std::endl;
		std::cerr << "Expected either a file path or coordinates in format 'x,y'" << std::endl;
		exit(1);
	}

	if (verbose) std::cout << "focal is (" << focal_coords[0] << ", " << focal_coords[1] << ")." << std::endl;
	return Vector2D(focal_coords[0], focal_coords[1]);
}

/* ===============================================================================================
   Main program
   =============================================================================================== */

/*!
 * @brief Main function for comparing 2D shapes with specified focal points
 *
 * This program compares two 2D polygons using various distance metrics, with user-specified
 * focal points for each polygon. The polygons are shifted so that their focal points are at
 * the origin before comparison, allowing for focal-point-relative shape analysis.
 *
 * The program supports multiple distance metrics:
 * - Fréchet distance: Measures curve similarity (minimum leash length)
 * - Ellipse-based: Compares fitted ellipse properties (inscribed, inscribing, least-squares)
 * - Curvature-based: Willmore energy and Wasserstein distance on curvature distributions
 * - Persistence diagrams: Topological comparison using persistent homology
 *
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return 0 on success, 1 on failure
 *
 * @see parse_args() for command-line argument details
 * @see parse_focal_point() for focal point specification format
 */
int main(int argc, char **argv)
{

/*	==========================================================================================
	Show usage if needed
	========================================================================================== */

	if( argc < 2 )
	{
		usage(argv);
		return -1;
	}

	std::string input = argv[1];
	if( input == "-h" || input == "-help" )
	{
		usage(argv);
		return -1;
	}

/*	==========================================================================================
	Read in all inputs (some values may have been preset)
	========================================================================================== */

	std::string INfile1;
	std::string INfile2;
	std::string INfocal1;
	std::string INfocal2;
	int disttype = 0;
	double microns_per_pixel1 = 1.0; // default: no conversion (1 pixel = 1 micron)
	double microns_per_pixel2 = 1.0; // default: no conversion (1 pixel = 1 micron)
	bool verbose = false; // default: no verbose output

        if (!parse_args(argc, argv, &INfile1, &INfocal1, &INfile2, &INfocal2, &disttype, &microns_per_pixel1, &microns_per_pixel2, &verbose)) return 1;

/*	==========================================================================================
	Read in the polygon1 from input file1
	========================================================================================== */

	int ndim;
	int npoint1;

	double *X1;
	X1 = nullptr;
	inout.read(INfile1, &ndim, &npoint1, &X1);

	Vector2D focal1 = parse_focal_point(INfocal1, verbose);

/*	==========================================================================================
	Store polygon1 as a polygon
	========================================================================================== */

	Polygon polygon1;
	pbuilder.clean_points(&npoint1, X1);
	pbuilder.buildPolygon(npoint1, X1, polygon1);

	// Check if focal point 1 is inside the polygon
	if (!polygon1.isPointInside(focal1)) {
		std::cerr << "WARNING: Focal point 1 (" << focal1.x << ", " << focal1.y
		          << ") is NOT inside polygon 1!" << std::endl;
	}

	polygon1.shift(focal1);
	// Apply pixel-to-micron conversion to polygon 1 and its focal point
	pbuilder.convertPixelsToMicrometers(polygon1, microns_per_pixel1);
	// focal1.x *= microns_per_pixel1;
	// focal1.y *= microns_per_pixel1;

	

/*	==========================================================================================
	Read in the polygon2 from input file2
	========================================================================================== */

	int npoint2;

	double *X2;
	X2 = nullptr;
	inout.read(INfile2, &ndim, &npoint2, &X2);

	Vector2D focal2 = parse_focal_point(INfocal2, verbose);

/*	==========================================================================================
	Store polygon2 as a polygon
	========================================================================================== */

	Polygon polygon2;
	pbuilder.clean_points(&npoint2, X2);
	pbuilder.buildPolygon(npoint2, X2, polygon2);

	// Check if focal point 2 is inside the polygon
	if (!polygon2.isPointInside(focal2)) {
		std::cerr << "WARNING: Focal point 2 (" << focal2.x << ", " << focal2.y
		          << ") is NOT inside polygon 2!" << std::endl;
	}

	polygon2.shift(focal2);
	// Apply pixel-to-micron conversion to polygon 2 and its focal point
	pbuilder.convertPixelsToMicrometers(polygon2, microns_per_pixel2);
	// focal2.x *= microns_per_pixel2;
	// focal2.y *= microns_per_pixel2;

	

/*	==========================================================================================
	Compute distances
	========================================================================================== */

	double dFrechet, dE_M, dE_m, dE_l, dW, dcOT1, w_q_dist;
	double a1_M, b1_M, r1_M, a2_M, b2_M, r2_M;
	double a1_m, b1_m, r1_m, a2_m, b2_m, r2_m;
	double a1_l, b1_l, r1_l, a2_l, b2_l, r2_l;
	double willmore1, willmore2;
	correa::PH0 f1(polygon1.vertices);
	correa::PH0 f2(polygon2.vertices);

	if(disttype == 0) {
		dFrechet = frechet.dFD(polygon1, polygon2);
	} else if(disttype == 1) {
		dE_m = ellipse.dEllipseMin(polygon1, polygon2, &a1_m, &b1_m, &a2_m, &b2_m);
		r1_m = a1_m/b1_m; r2_m = a2_m/b2_m;
		dE_M = ellipse.dEllipseMax(polygon1, polygon2, &a1_M, &b1_M, &a2_M, &b2_M);
		r1_M = a1_M/b1_M; r2_M = a2_M/b2_M;
		dE_l = ellipse.dEllipseLSQ(polygon1, polygon2, &a1_l, &b1_l, &a2_l, &b2_l);
		r1_l = a1_l/b1_l; r2_l = a2_l/b2_l;
	} else if(disttype == 2) {
		willmore1 = curv.Willmore(polygon1);
		willmore2 = curv.Willmore(polygon2);
		dW = std::abs(willmore1-willmore2);
		dcOT1 = curv.curvOT(polygon1, polygon2);
	} else if(disttype == 3){
		f1.Persistence();
		f2.Persistence();
		hera::AuctionParams<double> hera_params;
    	hera_params.max_num_phases = 800;
		hera_params.wasserstein_power = 2;
		hera_params.internal_p = 2;
		const std::vector<std::pair<double,double>> pd1 = f1.persistence_diagram();
		const std::vector<std::pair<double,double>> pd2 = f2.persistence_diagram();
		const hera::AuctionParams<double> params = hera_params;
		w_q_dist = hera::wasserstein_dist<std::vector<std::pair<double,double>>>(pd1, pd2, params);
	} else {
		dFrechet = frechet.dFD(polygon1, polygon2);
		dE_m = ellipse.dEllipseMin(polygon1, polygon2, &a1_m, &b1_m, &a2_m, &b2_m);
		r1_m = a1_m/b1_m; r2_m = a2_m/b2_m;
		dE_M = ellipse.dEllipseMax(polygon1, polygon2, &a1_M, &b1_M, &a2_M, &b2_M);
		r1_M = a1_M/b1_M; r2_M = a2_M/b2_M;
		dE_l = ellipse.dEllipseLSQ(polygon1, polygon2, &a1_l, &b1_l, &a2_l, &b2_l);
		r1_l = a1_l/b1_l; r2_l = a2_l/b2_l;
		willmore1 = curv.Willmore(polygon1);
		willmore2 = curv.Willmore(polygon2);
		dW = std::abs(willmore1-willmore2);
		dcOT1 = curv.curvOT(polygon1, polygon2);
		f1.Persistence();
		f2.Persistence();
		hera::AuctionParams<double> hera_params;
    	hera_params.max_num_phases = 800;
		hera_params.wasserstein_power = 2;
		hera_params.internal_p = 2;
		const std::vector<std::pair<double,double>> pd1 = f1.persistence_diagram();
		const std::vector<std::pair<double,double>> pd2 = f2.persistence_diagram();
		const hera::AuctionParams<double> params = hera_params;
		w_q_dist = hera::wasserstein_dist<std::vector<std::pair<double,double>>>(pd1, pd2, params);
	}

/*	==========================================================================================
	Print results
	========================================================================================== */

	// Info on polygon 1:

	double l1 = polygon1.length();
	double A1 = polygon1.area();

	if (verbose) {
		std::cout << " " << std::endl;
		std::cout << "Polygon 1 : " << std::endl;
		std::cout << "===========" << std::endl;
		std::cout << "Number of points in polygon       : " << npoint1 << std::endl;
		std::cout << "Length of polygon                 : " << l1 << std::endl;
		std::cout << "Area of polygon                   : " << A1 << std::endl;
		std::cout << "Sphericity (4*Pi*Area/L^2)        : " << 4*M_PI*A1/(l1*l1) << std::endl;
		if(disttype==1 || disttype==4) {
			std::cout << "Maximum volume inscribed ellipse  : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a1_M << " b : " << b1_M << " Aspect ratio: " << r1_M << std::endl;
			std::cout << "Minimum volume inscribing ellipse : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a1_m << " b : " << b1_m << " Aspect ratio: " << r1_m << std::endl;
			std::cout << "Least square ellipse              : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a1_l << " b : " << b1_l << " Aspect ratio: " << r1_l << std::endl;
		}
		if(disttype==2 || disttype==4) {
			std::cout << "Willmore energy of polygon        : " << willmore1 << std::endl;
		}
		if(disttype==3 || disttype==4) {
			std::cout << "Persistence diagram of polygon 1  : ";
			std::cout<< "There are " << f1.persistence_diagram().size() << " points in the persistence diagram." << std::endl;
			for (int i = 0; i < f1.persistence_diagram().size(); i++){
				std::cout << "(" << get<0>(f1.persistence_diagram()[i]) << ", " << get<1>(f1.persistence_diagram()[i]) << ")" << std::endl;
			};
		}
	}

	// Info on polygon 2:

	double l2 = polygon2.length();
	double A2 = polygon2.area();

	if (verbose) {
		std::cout << " " << std::endl;
		std::cout << "Polygon 2 : " << std::endl;
		std::cout << "===========" << std::endl;
		std::cout << "Number of points in polygon       : " << npoint2 << std::endl;
		std::cout << "Length of polygon                 : " << l2 << std::endl;
		std::cout << "Area of polygon                   : " << A2 << std::endl;
		std::cout << "Sphericity (4*Pi*Area/L^2)        : " << 4*M_PI*A2/(l2*l2) << std::endl;
		if(disttype==1 || disttype==4) {
			std::cout << "Maximum volume inscribed ellipse  : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a2_M << " b : " << b2_M << " Aspect ratio: " << r2_M << std::endl;
			std::cout << "Minimum volume inscribing ellipse : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a2_m << " b : " << b2_m << " Aspect ratio: " << r2_m << std::endl;
			std::cout << "Least square ellipse              : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a2_l << " b : " << b2_l << " Aspect ratio: " << r2_l << std::endl;
		}
		if(disttype==2 || disttype==4) {
			std::cout << "Willmore energy of polygon        : " << willmore2 << std::endl;
		}
		if(disttype==3 || disttype==4) {
			std::cout << "Persistence diagram of polygon 2  : ";
			std::cout<< "There are " << f2.persistence_diagram().size() << " points in the persistence diagram." << std::endl;
			for (int i = 0; i < f2.persistence_diagram().size(); i++){
				std::cout << "(" << get<1>(f2.persistence_diagram()[i]) << ", " << get<1>(f2.persistence_diagram()[i]) << ")" << std::endl;
			};
		}
	}

	if (verbose) {
		std::cout << " " << std::endl;

		std::cout << std::setw(10) << std::fixed << std::setprecision(5);

		// Distances between the two polygons
		std::cout << " " << std::endl;
		std::cout << "Distances between the 2 polygons : " << std::endl;
		std::cout << "===================================" << std::endl;
		std::cout << "Sphericity 4*Pi*(A1/l1^2 - A2/l2^2) : " << 4*M_PI*std::abs(A1/(l1*l1)-A2/(l2*l2)) << std::endl;
		if(disttype==0 || disttype == 4) {
			std::cout << "Frechet distance                                    : " << dFrechet << std::endl;
		}
		if(disttype==1 || disttype == 4) {
			std::cout << "Distance (inscribed ellipse)                        : " << dE_M << std::endl;
			std::cout << "Distance (inscribing ellipse)                       : " << dE_m << std::endl;
			std::cout << "Distance (least square ellipse)                     : " << dE_l << std::endl;
		}
		if(disttype==2 || disttype == 4) {
			std::cout << "Distance (Willmore)                                 : " << dW << std::endl;
		}
		if(disttype == 4) {
			std::cout << "Distance (Wasserstein-curvature)                    : " << std::scientific << dcOT1 << std::endl;
		}
		if(disttype==3 || disttype == 4) {
			std::cout << "Distance (2-Wasserstein between persitence diagrams): " << w_q_dist << std::endl;
		}
	}
	
	return 0;

}

/* ===============================================================================================
   Usage
   =============================================================================================== */

static void usage(char** argv)
{
    std::cout << "\n\n" <<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=                                Comp2DShapesFocal                                             ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     Usage is:                                                                                ="<<std::endl;
    std::cout << "     " << "=          Comp2DShapesFocal -i1 FILE1 -f1 FOCAL1 -i2 FILE2 -f2 FOCAL2 -d disttype [options]  ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     Required arguments:                                                                      ="<<std::endl;
    std::cout << "     " << "=               -i1 FILE1     --> Input file 1 (Curve; ascii or csv file with 1 point / line)  ="<<std::endl;
    std::cout << "     " << "=               -f1 FOCAL1    --> Focal point 1 (file path or inline coords like \"x,y\")        ="<<std::endl;
	std::cout << "     " << "=               -i2 FILE2     --> Input file 2 (Curve; ascii or csv file with 1 point / line)  ="<<std::endl;
	std::cout << "     " << "=               -f2 FOCAL2    --> Focal point 2 (file path or inline coords like \"x,y\")        ="<<std::endl;
    std::cout << "     " << "=               -d disttype   --> Distance type flag:                                          ="<<std::endl;
    std::cout << "     " << "=                                   (0) Frechet distance                                       ="<<std::endl;
    std::cout << "     " << "=                                   (1) Aspect ratio distances (based on ellipses)             ="<<std::endl;
    std::cout << "     " << "=                                   (2) Curvature-based distances                              ="<<std::endl;
    std::cout << "     " << "=                                   (3) 2-Wasserstein distance between persistence diagrams    ="<<std::endl;
    std::cout << "     " << "=                                   (4) All                                                    ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     Optional arguments:                                                                      ="<<std::endl;
    std::cout << "     " << "=               -mpp VALUE    --> Microns per pixel for both contours (default: 1.0)           ="<<std::endl;
    std::cout << "     " << "=               -mpp1 VALUE   --> Microns per pixel for contour 1 (default: 1.0)               ="<<std::endl;
    std::cout << "     " << "=               -mpp2 VALUE   --> Microns per pixel for contour 2 (default: 1.0)               ="<<std::endl;
    std::cout << "     " << "=               -v, --verbose --> Enable verbose output (print detailed polygon info)           ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     Note: Use -mpp1 and -mpp2 to specify different scales for each contour                   ="<<std::endl;
    std::cout << "     " << "=           (useful when contours are from different imaging sessions/magnifications)          ="<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "\n\n" <<std::endl;
}

/* ===============================================================================================
   Parse Argument from command line:

   =============================================================================================== */

/*!
 * @brief Parse command-line arguments for Comp2DShapesFocal
 *
 * Parses command-line arguments for the Comp2DShapesFocal program, extracting
 * input files, focal points, distance type, microns per pixel settings, and
 * verbose flag.
 *
 * @param argc Argument count from main()
 * @param argv Argument vector from main()
 * @param[out] file1 Path to first input polygon file (set by -i1 flag)
 * @param[out] focal1 Focal point for first polygon - file path or "x,y" format (set by -f1 flag)
 * @param[out] file2 Path to second input polygon file (set by -i2 flag)
 * @param[out] focal2 Focal point for second polygon - file path or "x,y" format (set by -f2 flag)
 * @param[out] disttype Distance type flag (set by -d flag):
 *                      - 0: Fréchet distance - measures similarity between curves, representing the minimum
 *                           leash length needed for a person walking a dog along each curve
 *                      - 1: Ellipse-based distances - compares aspect ratios of inscribed, inscribing,
 *                           and least-squares fitted ellipses
 *                      - 2: Curvature-based distances - includes Willmore energy (bending energy) and
 *                           Wasserstein distance between curvature distributions
 *                      - 3: 2-Wasserstein distance between persistence diagrams - topological shape descriptor
 *                           based on persistent homology
 *                      - 4: All distances - computes all of the above metrics
 * @param[out] microns_per_pixel1 Microns per pixel conversion for polygon 1 (set by -mpp1 or -mpp flag, default: 1.0)
 * @param[out] microns_per_pixel2 Microns per pixel conversion for polygon 2 (set by -mpp2 or -mpp flag, default: 1.0)
 * @param[out] verbose Enable verbose output (set by -v or --verbose flag, default: false)
 * @return true if arguments were parsed successfully, false otherwise
 */
bool parse_args(int argc, char **argv, std::string *file1, std::string *focal1, std::string *file2, std::string *focal2, int *disttype, double *microns_per_pixel1, double *microns_per_pixel2, bool *verbose)
{
//
// Make sure we have at least two parameters....
//
	std::string param;
	if (argc == 1)
	{
		return false;
	}
	else
	{
		for (int i = 1; i < argc; i++)
		{
			param = argv[i];

			// Handle flag-only arguments (no value)
			if (param == "-v" || param == "--verbose") {
				*verbose = true;
				continue;
			}

			// Skip if this is the last argument and it's not a flag
			if (i >= argc - 1) continue;

			if (param == "-i1") {
				*file1 = argv[i + 1];
				i++;
			}
			else if (param == "-i2") {
				*file2 = argv[i + 1];
				i++;
			}
			else if (param == "-d") {
				*disttype = std::atoi(argv[i + 1]);
				i++;
			}
			else if (param == "-f1") {
				*focal1 = argv[i + 1];
				i++;
			}
			else if (param == "-f2") {
				*focal2 = argv[i + 1];
				i++;
			}
			else if (param == "-mpp1") {
				*microns_per_pixel1 = std::atof(argv[i + 1]);
				i++;
			}
			else if (param == "-mpp2") {
				*microns_per_pixel2 = std::atof(argv[i + 1]);
				i++;
			}
			// Backward compatibility: if -mpp is used, apply to both
			else if (param == "-mpp") {
				double mpp = std::atof(argv[i + 1]);
				*microns_per_pixel1 = mpp;
				*microns_per_pixel2 = mpp;
				i++;
			}
		}
  	}
	return true;
}
