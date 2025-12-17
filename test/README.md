# Correa Tests

This directory contains test suites for the Correa library.

## Test Files

- **test_focal_point_validation.py**: Comprehensive tests for focal point validation functionality
  - Tests the ray casting algorithm (Jordan curve theorem) implementation
  - Verifies that focal points outside polygons are properly detected with warnings
  - Tests edge cases including concave polygons, boundary conditions, and tangential crossings
  - Tests integration with scaling and unit conversion features
  - Validates that computation continues even when focal point is outside (with warning)

- **tests.py**: Basic integration test for polygon comparison

- **demo_focal_warning.py**: Interactive demonstration of focal point validation warnings
  - Shows warnings when focal points are outside contours
  - Demonstrates that computation continues despite warnings
  - Examples with verbose mode enabled

## Requirements

Install test dependencies:

```bash
pip install pytest pytest-cov
```

## Running Tests

### Run all tests with verbose output:
```bash
pytest test/ -v
```

### Run specific test file:
```bash
pytest test/test_focal_point_validation.py -v
```

### Run specific test class:
```bash
pytest test/test_focal_point_validation.py::TestFocalPointValidation -v
```

### Run specific test:
```bash
pytest test/test_focal_point_validation.py::TestFocalPointValidation::test_focal_point_inside_square -v
```

### Run with coverage report:
```bash
pytest test/ --cov=correa --cov-report=html
```

### Run with detailed output on failure:
```bash
pytest test/ -v --tb=long
```

## Test Structure

### TestFocalPointValidation
Tests the main focal point validation functionality:
- Focal points inside polygons (should succeed)
- Focal points outside polygons (should raise errors)
- High-level API compatibility
- Verbose mode output
- Error message formatting

### TestEdgeCases
Tests edge cases for the point-in-polygon algorithm:
- Concave polygons
- Points in cut-out areas
- Boundary conditions
- Tangential crossings

### TestScaleAndConversion
Tests that validation works correctly with:
- Pixel-to-micron conversion
- Area-based scaling
- Ensures validation happens before transformations

## Adding New Tests

When adding new tests:

1. Follow the existing test structure using pytest
2. Use descriptive test names starting with `test_`
3. Add docstrings explaining what each test validates
4. Use fixtures for reusable test data
5. Group related tests in classes
6. Test both success and failure cases

## Expected Test Data

Tests use:
- **Simple geometric shapes**: Created dynamically in fixtures (squares, triangles, L-shapes)
- **Example contours**: From `../examples/contour*.csv`
- **Temporary files**: Created in pytest's temporary directory

## Continuous Integration

These tests should be run:
- Before committing changes
- In CI/CD pipelines
- After building the C++ extensions

## Known Issues

None currently. If you encounter issues:
1. Ensure the `build/` directory contains the compiled `_correa` module
2. Check that example files exist in `examples/` directory
3. Verify pytest is installed: `pytest --version`
