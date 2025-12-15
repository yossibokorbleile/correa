#!/usr/bin/env python3
"""
Test script to verify that the updated create_polygon_focal_point function works correctly.
"""

import sys
import os

# Add the bindings directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'bindings'))

try:
    import correa
    
    # Test the function call with the parameters from the notebook
    # Using a sample contour file if it exists
    test_contour_path = "/Users/yossi/PDaMSoC-data/X2/cell/contours/001_contour.csv"
    test_focal_point = [100.0, 100.0]  # Sample focal point coordinates
    test_scale_factor = 0.3155
    
    if os.path.exists(test_contour_path):
        print("Testing create_polygon_focal_point function...")
        
        # Test with convert_to_microns_factor
        polygon_scaled = correa.create_polygon_focal_point(
            test_contour_path, 
            test_focal_point, 
            scale_by_area=True, 
            convert_to_microns_factor=test_scale_factor
        )
        print("✓ Function call with convert_to_microns_factor succeeded")
        print(f"  Polygon has {polygon_scaled.size()} vertices")
        
        # Test with default convert_to_microns_factor (should be 1.0)
        polygon_default = correa.create_polygon_focal_point(
            test_contour_path, 
            test_focal_point, 
            scale_by_area=False
        )
        print("✓ Function call with default convert_to_microns_factor succeeded")
        print(f"  Polygon has {polygon_default.size()} vertices")
        
        print("\n✓ All tests passed! The function call should now work in the notebook.")
        
    else:
        print(f"Test contour file not found at {test_contour_path}")
        print("Please ensure the PDaMSoC-data directory is accessible.")
        
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure the correa module is properly built and accessible.")
except Exception as e:
    print(f"Error during testing: {e}")
    import traceback
    traceback.print_exc()
