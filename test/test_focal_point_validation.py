#!/usr/bin/env python3
"""
Test suite for focal point validation functionality.

Tests the ray casting algorithm (Jordan curve theorem) implementation
and verifies that focal points outside polygons are properly detected
and cause appropriate errors.
"""

import sys
import os
import pytest

# Add the build directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bindings'))

import _correa
import correa


class TestFocalPointValidation:
    """Test focal point validation for polygons."""

    @pytest.fixture
    def example_contour1(self):
        """Path to example contour 1."""
        return os.path.join(os.path.dirname(__file__), '..', 'examples', 'contour1.csv')

    @pytest.fixture
    def example_contour2(self):
        """Path to example contour 2."""
        return os.path.join(os.path.dirname(__file__), '..', 'examples', 'contour2.csv')

    @pytest.fixture
    def simple_square(self, tmp_path):
        """Create a simple square polygon for testing (0,0) to (10,10)."""
        square_file = tmp_path / "square.csv"
        with open(square_file, 'w') as f:
            # Create a square from (0,0) to (10,10)
            for i in range(11):
                f.write(f"{i},0\n")  # Bottom edge
            for i in range(1, 11):
                f.write(f"10,{i}\n")  # Right edge
            for i in range(9, -1, -1):
                f.write(f"{i},10\n")  # Top edge
            for i in range(9, 0, -1):
                f.write(f"0,{i}\n")  # Left edge
        return str(square_file)

    def test_focal_point_inside_square(self, simple_square):
        """Test that a focal point clearly inside the polygon is accepted."""
        # Center of the square (5, 5) should be inside
        polygon = _correa.PyPolygon(simple_square, [5.0, 5.0], True, False, 1.0)
        assert polygon.size() > 0

    def test_focal_point_outside_square_left(self, simple_square, capfd):
        """Test that a focal point to the left of the polygon prints a warning."""
        # Should not raise an exception, but should print a warning
        polygon = _correa.PyPolygon(simple_square, [-5.0, 5.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err
        assert "NOT inside" in captured.err

    def test_focal_point_outside_square_right(self, simple_square, capfd):
        """Test that a focal point to the right of the polygon prints a warning."""
        polygon = _correa.PyPolygon(simple_square, [15.0, 5.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err
        assert "NOT inside" in captured.err

    def test_focal_point_outside_square_above(self, simple_square, capfd):
        """Test that a focal point above the polygon prints a warning."""
        polygon = _correa.PyPolygon(simple_square, [5.0, 15.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err
        assert "NOT inside" in captured.err

    def test_focal_point_outside_square_below(self, simple_square, capfd):
        """Test that a focal point below the polygon prints a warning."""
        polygon = _correa.PyPolygon(simple_square, [5.0, -5.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err
        assert "NOT inside" in captured.err

    def test_focal_point_at_origin_outside(self, simple_square, capfd):
        """Test that a clearly outside point prints a warning."""
        # (0, 0) is on the boundary/vertex, which should be treated as outside
        # or at minimum, the ray casting might have edge cases
        # Let's use a clearly outside point
        polygon = _correa.PyPolygon(simple_square, [-1.0, -1.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err

    def test_focal_point_near_edge_inside(self, simple_square):
        """Test that a focal point very close to the edge but inside is accepted."""
        # Point at (0.5, 5.0) should be inside
        polygon = _correa.PyPolygon(simple_square, [0.5, 5.0], True, False, 1.0)
        assert polygon.size() > 0

    def test_focal_point_example_contour_inside(self, example_contour1):
        """Test with example contour 1 using a focal point likely inside."""
        # Based on the contour data (around 83-122 x, 52-56 y), use center
        polygon = _correa.PyPolygon(example_contour1, [105.0, 54.0], True, False, 1.0)
        assert polygon.size() > 0

    def test_focal_point_example_contour_outside(self, example_contour1, capfd):
        """Test with example contour 1 using a focal point clearly outside."""
        polygon = _correa.PyPolygon(example_contour1, [0.0, 0.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err
        assert "NOT inside" in captured.err

    def test_focal_point_far_outside(self, example_contour1, capfd):
        """Test with a focal point very far from the polygon."""
        polygon = _correa.PyPolygon(example_contour1, [1000.0, 1000.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err

    def test_high_level_api_focal_point_inside(self, example_contour1):
        """Test the high-level API with a focal point inside."""
        polygon = correa.create_polygon_focal_point(
            example_contour1,
            [105.0, 54.0],
            clean_points=True,
            scale_by_area=False
        )
        assert polygon.size() > 0

    def test_high_level_api_focal_point_outside(self, example_contour1, capfd):
        """Test the high-level API with a focal point outside."""
        polygon = correa.create_polygon_focal_point(
            example_contour1,
            [0.0, 0.0],
            clean_points=True,
            scale_by_area=False
        )
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err

    def test_verbose_mode_shows_check(self, simple_square, capfd):
        """Test that verbose mode shows the focal point check."""
        _correa.set_verbose(True)
        try:
            polygon = _correa.PyPolygon(simple_square, [5.0, 5.0], True, False, 1.0)
            captured = capfd.readouterr()
            assert "Checking if focal point" in captured.out
        finally:
            _correa.set_verbose(False)

    def test_warning_message_contains_coordinates(self, simple_square, capfd):
        """Test that warning messages include the focal point coordinates."""
        polygon = _correa.PyPolygon(simple_square, [-5.0, -5.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "-5" in captured.err
        assert "focal point" in captured.err.lower() or "Focal point" in captured.err


class TestEdgeCases:
    """Test edge cases for the point-in-polygon algorithm."""

    @pytest.fixture
    def triangle(self, tmp_path):
        """Create a simple triangle for testing."""
        triangle_file = tmp_path / "triangle.csv"
        with open(triangle_file, 'w') as f:
            # Triangle with vertices at (0,0), (10,0), (5,10)
            for i in range(11):
                f.write(f"{i},0\n")  # Base
            for i in range(1, 11):
                x = 10 - i/2
                y = i
                f.write(f"{x},{y}\n")  # Right edge
            for i in range(9, -1, -1):
                x = i/2
                y = 10 - (10-i)
                f.write(f"{x},{y}\n")  # Left edge
        return str(triangle_file)

    @pytest.fixture
    def concave_polygon(self, tmp_path):
        """Create a concave polygon (L-shape) for testing."""
        l_shape_file = tmp_path / "l_shape.csv"
        with open(l_shape_file, 'w') as f:
            # L-shaped polygon
            points = [
                (0, 0), (5, 0), (5, 5), (10, 5), (10, 10),
                (0, 10), (0, 0)
            ]
            for x, y in points:
                f.write(f"{x},{y}\n")
        return str(l_shape_file)

    def test_concave_polygon_inside_main_area(self, concave_polygon):
        """Test focal point inside the main area of a concave polygon."""
        # Point (2, 2) should be inside the L-shape
        polygon = _correa.PyPolygon(concave_polygon, [2.0, 2.0], True, False, 1.0)
        assert polygon.size() > 0

    def test_concave_polygon_inside_extension(self, concave_polygon):
        """Test focal point inside the extension of a concave polygon."""
        # Point (7, 7) should be inside the top-right part of the L
        polygon = _correa.PyPolygon(concave_polygon, [7.0, 7.0], True, False, 1.0)
        assert polygon.size() > 0

    def test_concave_polygon_outside_concave_part(self, concave_polygon, capfd):
        """Test focal point in the concave (cut-out) area."""
        # Point (7, 2) should be outside (in the cut-out area)
        polygon = _correa.PyPolygon(concave_polygon, [7.0, 2.0], True, False, 1.0)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err


class TestScaleAndConversion:
    """Test that focal point validation works with scaling and unit conversion."""

    @pytest.fixture
    def simple_square(self, tmp_path):
        """Create a simple square polygon for testing."""
        square_file = tmp_path / "square.csv"
        with open(square_file, 'w') as f:
            for i in range(11):
                f.write(f"{i},0\n")
            for i in range(1, 11):
                f.write(f"10,{i}\n")
            for i in range(9, -1, -1):
                f.write(f"{i},10\n")
            for i in range(9, 0, -1):
                f.write(f"0,{i}\n")
        return str(square_file)

    def test_focal_point_with_micron_conversion(self, simple_square):
        """Test that focal point validation works before unit conversion."""
        # The validation should happen BEFORE unit conversion
        # So the focal point should be in the original pixel coordinates
        polygon = _correa.PyPolygon(simple_square, [5.0, 5.0], True, False, 0.5)
        assert polygon.size() > 0

    def test_focal_point_outside_with_micron_conversion(self, simple_square, capfd):
        """Test that invalid focal points are detected before unit conversion."""
        polygon = _correa.PyPolygon(simple_square, [-5.0, 5.0], True, False, 0.5)
        assert polygon.size() > 0
        captured = capfd.readouterr()
        assert "WARNING" in captured.err

    def test_focal_point_with_area_scaling(self, simple_square):
        """Test that focal point validation works with area scaling."""
        polygon = _correa.PyPolygon(simple_square, [5.0, 5.0], True, True, 1.0)
        assert polygon.size() > 0


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "--tb=short"])
