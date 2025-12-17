#!/usr/bin/env python3
"""
Basic functionality tests for Correa library.

Tests basic polygon creation, comparison, and metric calculations.
"""

import sys
import os
import pytest

# Add the build directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bindings'))

import correa


class TestPolygonCreation:
    """Test basic polygon creation functionality."""

    @pytest.fixture
    def example_contour1(self):
        """Path to example contour 1."""
        return os.path.join(os.path.dirname(__file__), '..', 'examples', 'contour1.csv')

    @pytest.fixture
    def example_contour2(self):
        """Path to example contour 2."""
        return os.path.join(os.path.dirname(__file__), '..', 'examples', 'contour2.csv')

    def test_create_polygon_basic(self, example_contour1):
        """Test basic polygon creation."""
        p1 = correa.create_polygon(example_contour1)
        assert p1.size() > 0
        assert p1.area() > 0
        assert p1.length() > 0

    def test_create_polygon_with_cleaning(self, example_contour1):
        """Test polygon creation with point cleaning."""
        p1 = correa.create_polygon(example_contour1, clean_points=True)
        assert p1.size() > 0

    def test_create_polygon_without_cleaning(self, example_contour1):
        """Test polygon creation without point cleaning."""
        p1 = correa.create_polygon(example_contour1, clean_points=False)
        assert p1.size() > 0

    def test_create_polygon_with_area_scaling(self, example_contour1):
        """Test polygon creation with area scaling."""
        p1 = correa.create_polygon(example_contour1, scale_by_area=True)
        # After area scaling, area should be 100
        assert abs(p1.area() - 100.0) < 1e-6

    def test_create_polygon_with_micron_conversion(self, example_contour1):
        """Test polygon creation with pixel-to-micron conversion."""
        scale_factor = 0.5
        p1 = correa.create_polygon(example_contour1, convert_to_microns_factor=scale_factor)
        assert p1.size() > 0


class TestPolygonProperties:
    """Test polygon property calculations."""

    @pytest.fixture
    def polygon(self):
        """Create a test polygon."""
        path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'contour1.csv')
        return correa.create_polygon(path)

    def test_polygon_size(self, polygon):
        """Test polygon size (number of vertices)."""
        assert polygon.size() > 0
        assert isinstance(polygon.size(), int)

    def test_polygon_area(self, polygon):
        """Test polygon area calculation."""
        area = polygon.area()
        assert area > 0
        assert isinstance(area, float)

    def test_polygon_length(self, polygon):
        """Test polygon perimeter length."""
        length = polygon.length()
        assert length > 0
        assert isinstance(length, float)

    def test_polygon_vertices(self, polygon):
        """Test polygon vertices retrieval."""
        vertices = polygon.vertices()
        assert len(vertices) == polygon.size()
        # Each vertex should be a pair of coordinates
        for vertex in vertices:
            assert len(vertex) == 2
            assert isinstance(vertex[0], float)
            assert isinstance(vertex[1], float)

    def test_ellipse_properties(self, polygon):
        """Test ellipse fitting properties."""
        # Maximum inscribed ellipse
        a_max = polygon.ellipse_max_a()
        b_max = polygon.ellipse_max_b()
        ratio_max = polygon.ellipse_max_ratio()
        assert a_max > 0
        assert b_max > 0
        assert abs(ratio_max - a_max/b_max) < 1e-6

        # Minimum inscribing ellipse
        a_min = polygon.ellipse_min_a()
        b_min = polygon.ellipse_min_b()
        ratio_min = polygon.ellipse_min_ratio()
        assert a_min > 0
        assert b_min > 0
        assert abs(ratio_min - a_min/b_min) < 1e-6

        # Least squares ellipse
        a_lsq = polygon.ellipse_lsq_a()
        b_lsq = polygon.ellipse_lsq_b()
        ratio_lsq = polygon.ellipse_lsq_ratio()
        assert a_lsq > 0
        assert b_lsq > 0
        assert abs(ratio_lsq - a_lsq/b_lsq) < 1e-6

    def test_willmore_energy(self, polygon):
        """Test Willmore energy calculation."""
        willmore = polygon.willmore()
        assert isinstance(willmore, float)
        assert willmore >= 0  # Willmore energy should be non-negative

    def test_persistence_diagram(self, polygon):
        """Test persistence diagram calculation."""
        pd = polygon.persistence_diagram()
        assert isinstance(pd, list)
        # Each point in persistence diagram should have birth and death values
        for point in pd:
            assert len(point) == 2
            birth, death = point
            assert isinstance(birth, float)
            assert isinstance(death, float)
            # Death should be >= birth
            assert death >= birth


class TestPolygonComparison:
    """Test polygon comparison metrics."""

    @pytest.fixture
    def polygons(self):
        """Create two test polygons."""
        path1 = os.path.join(os.path.dirname(__file__), '..', 'examples', 'contour1.csv')
        path2 = os.path.join(os.path.dirname(__file__), '..', 'examples', 'contour2.csv')
        p1 = correa.create_polygon(path1)
        p2 = correa.create_polygon(path2)
        return p1, p2

    def test_compare_polygons_basic(self, polygons):
        """Test basic polygon comparison."""
        p1, p2 = polygons
        result = correa.compare_polygons(p1, p2, q=2, verbose=False)
        assert isinstance(result, list)
        assert len(result) == 7  # Should return 7 distance metrics

    def test_frechet_distance(self, polygons):
        """Test Fréchet distance calculation."""
        p1, p2 = polygons
        distance = correa.frechet_distance(p1, p2)
        assert isinstance(distance, float)
        assert distance >= 0

    def test_wasserstein_distance(self, polygons):
        """Test Wasserstein distance calculation."""
        p1, p2 = polygons
        distance = correa.wasserstein_distance(p1, p2, q=2)
        assert isinstance(distance, float)
        assert distance >= 0

    def test_ellipse_distances(self, polygons):
        """Test ellipse-based distance metrics."""
        p1, p2 = polygons

        max_dist = correa.max_ellipse_distance(p1, p2)
        assert isinstance(max_dist, float)
        assert max_dist >= 0

        min_dist = correa.min_ellipse_distance(p1, p2)
        assert isinstance(min_dist, float)
        assert min_dist >= 0

        lsq_dist = correa.lsq_ellipse_distance(p1, p2)
        assert isinstance(lsq_dist, float)
        assert lsq_dist >= 0

    def test_willmore_distance(self, polygons):
        """Test Willmore distance calculation."""
        p1, p2 = polygons
        distance = correa.willmore_distance(p1, p2)
        assert isinstance(distance, float)
        assert distance >= 0

    def test_curv_ot_distance(self, polygons):
        """Test curvature optimal transport distance."""
        p1, p2 = polygons
        distance = correa.curv_ot_distance(p1, p2)
        assert isinstance(distance, float)
        assert distance >= 0

    def test_compare_polygons_verbose(self, polygons, capfd):
        """Test polygon comparison with verbose output."""
        p1, p2 = polygons
        result = correa.compare_polygons(p1, p2, q=2, verbose=True)
        captured = capfd.readouterr()
        # Verbose mode should print comparison results
        assert len(captured.err) > 0 or len(captured.out) > 0


class TestVerbosityControl:
    """Test verbosity control functionality."""

    def test_verbose_default_false(self):
        """Test that verbose is False by default."""
        assert correa.get_verbose() == False

    def test_set_verbose_true(self):
        """Test setting verbose to True."""
        correa.set_verbose(True)
        assert correa.get_verbose() == True
        correa.set_verbose(False)  # Clean up

    def test_set_verbose_false(self):
        """Test setting verbose to False."""
        correa.set_verbose(True)
        correa.set_verbose(False)
        assert correa.get_verbose() == False


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "--tb=short"])
