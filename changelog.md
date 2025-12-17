# Change log

## v1.2.3 2025-12-17
- Added focal point validation using ray casting algorithm (Jordan curve theorem)
- Prints warning when focal point is not inside the polygon contour (computation continues)
- Enhanced tangential crossing detection for edge cases in point-in-polygon test
- Added conditional verbose output for debugging focal point checks
- Comprehensive test suite with 40+ tests covering validation, edge cases, and integration
- Supports both Python bindings and standalone C++ programs

## v1.0.1 2024-02-05
- fixed bug in building polygons with specified focal point, which was being changed by centering.
- minor changes to output formatting.