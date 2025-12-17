#!/usr/bin/env python3
"""
Demonstration of focal point validation warnings.

This script shows how the library warns when a focal point
is outside the polygon contour, but continues computation.
"""

import sys
import os

# Add the build directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bindings'))

import _correa
import correa

def main():
    print("=" * 70)
    print("Focal Point Validation Demo")
    print("=" * 70)
    print()

    example_contour = os.path.join(os.path.dirname(__file__), '..', 'examples', 'contour1.csv')

    print("Example contour coordinates range: ~(83-122, 52-56)")
    print()

    # Test 1: Focal point inside
    print("-" * 70)
    print("Test 1: Focal point INSIDE the contour")
    print("-" * 70)
    print("Creating polygon with focal point (105, 54)...")
    poly1 = _correa.PyPolygon(example_contour, [105.0, 54.0], True, False, 1.0)
    print(f"✓ Polygon created successfully with {poly1.size()} vertices")
    print(f"  Area: {poly1.area():.2f}")
    print(f"  No warnings (focal point is inside)")
    print()

    # Test 2: Focal point outside (origin)
    print("-" * 70)
    print("Test 2: Focal point OUTSIDE the contour (at origin)")
    print("-" * 70)
    print("Creating polygon with focal point (0, 0)...")
    print("Expected: WARNING message (but computation continues)")
    print()
    poly2 = _correa.PyPolygon(example_contour, [0.0, 0.0], True, False, 1.0)
    print()
    print(f"✓ Polygon created successfully with {poly2.size()} vertices")
    print(f"  Area: {poly2.area():.2f}")
    print(f"  (Computation continued despite warning)")
    print()

    # Test 3: Focal point far outside
    print("-" * 70)
    print("Test 3: Focal point FAR OUTSIDE the contour")
    print("-" * 70)
    print("Creating polygon with focal point (1000, 1000)...")
    print("Expected: WARNING message (but computation continues)")
    print()
    poly3 = _correa.PyPolygon(example_contour, [1000.0, 1000.0], True, False, 1.0)
    print()
    print(f"✓ Polygon created successfully with {poly3.size()} vertices")
    print(f"  Area: {poly3.area():.2f}")
    print(f"  (Computation continued despite warning)")
    print()

    # Test 4: With verbose mode
    print("-" * 70)
    print("Test 4: With VERBOSE mode enabled")
    print("-" * 70)
    print("Enabling verbose mode...")
    _correa.set_verbose(True)
    print("Creating polygon with focal point (200, 100) - outside...")
    print()
    poly4 = _correa.PyPolygon(example_contour, [200.0, 100.0], True, False, 1.0)
    _correa.set_verbose(False)
    print()
    print(f"✓ Polygon created with {poly4.size()} vertices")
    print(f"  (Verbose mode shows detailed validation steps)")
    print()

    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print("• Focal points inside the contour: ✓ No warnings")
    print("• Focal points outside the contour: ⚠ Warning printed")
    print("• Computation continues in both cases")
    print("• Verbose mode shows detailed validation process")
    print("=" * 70)

if __name__ == "__main__":
    main()
