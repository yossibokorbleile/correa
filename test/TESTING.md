# Testing Guide for Correa v1.2.3

## Overview

This guide explains the testing infrastructure for Correa, with a focus on the new focal point validation feature added in v1.2.3.

## Test Suite Structure

```
test/
├── conftest.py                        # Pytest configuration and shared fixtures
├── pytest.ini                         # Pytest settings
├── requirements.txt                   # Test dependencies
├── run_tests.sh                       # Convenient test runner script
├── test_focal_point_validation.py     # Tests for focal point validation (NEW in v1.2.3)
├── test_basic_functionality.py        # Tests for core functionality
├── tests.py                           # Original basic test (legacy)
├── data/                              # Test data directory
└── README.md                          # Testing documentation
```

## New in v1.2.3: Focal Point Validation Tests

The `test_focal_point_validation.py` file contains comprehensive tests for the new focal point validation feature:

### Test Coverage

1. **Basic Validation**
   - Points clearly inside polygons (should succeed without warnings)
   - Points clearly outside polygons (should print warnings but continue)
   - Points near boundaries
   - Points at various positions (left, right, above, below)

2. **Edge Cases**
   - Concave polygons (L-shapes, etc.)
   - Points in cut-out areas of concave polygons
   - Tangential ray crossings
   - Horizontal edge handling

3. **Integration**
   - Works with pixel-to-micron conversion
   - Works with area-based scaling
   - Validation happens before transformations
   - Both low-level and high-level API compatibility

4. **Warning Handling**
   - Warning messages printed to stderr for outside focal points
   - Warning messages contain focal point coordinates
   - Verbose mode shows validation steps
   - Computation continues even when warnings are issued

## Running Tests

### Quick Start

```bash
# From the project root
cd test
./run_tests.sh
```

### Specific Test Suites

```bash
# Test only focal point validation
pytest test_focal_point_validation.py -v

# Test only basic functionality
pytest test_basic_functionality.py -v

# Run all tests
pytest -v
```

### With Coverage

```bash
./run_tests.sh --coverage
```

This generates an HTML coverage report in `htmlcov/index.html`.

### Fast Tests Only

```bash
./run_tests.sh --fast
```

Skips tests marked as slow.

## Test Fixtures

### Shared Fixtures (in conftest.py)

- `examples_dir`: Path to example contour files
- `test_data_dir`: Path to test data directory
- `reset_verbose`: Automatically resets verbose mode after each test

### Test-Specific Fixtures

- `simple_square`: Creates a 10x10 square polygon for testing
- `triangle`: Creates a triangular polygon
- `concave_polygon`: Creates an L-shaped polygon for concave testing
- `example_contour1`, `example_contour2`: Paths to example files

## Writing New Tests

### Template

```python
import pytest
import _correa

class TestMyFeature:
    """Test description."""

    @pytest.fixture
    def my_fixture(self):
        """Fixture description."""
        # Setup code
        return test_data

    def test_something(self, my_fixture):
        """Test that something works."""
        result = my_function(my_fixture)
        assert result == expected_value

    def test_warning_case(self, capsys):
        """Test that warnings are properly printed."""
        result = my_function(invalid_input)
        # Function should still work but print a warning
        assert result is not None
        captured = capsys.readouterr()
        assert "WARNING" in captured.err
```

### Best Practices

1. **Use descriptive names**: Test names should clearly describe what they test
2. **One assertion per test**: Focus each test on a single behavior
3. **Use fixtures**: Share setup code via fixtures
4. **Test both success and failure**: Test expected behavior AND error cases
5. **Add docstrings**: Explain what each test validates
6. **Group related tests**: Use classes to organize related tests

## Continuous Integration

These tests should be integrated into your CI/CD pipeline:

```yaml
# Example GitHub Actions
- name: Install dependencies
  run: pip install -r test/requirements.txt

- name: Build project
  run: |
    mkdir build && cd build
    cmake ..
    cmake --build .

- name: Run tests
  run: cd test && ./run_tests.sh --coverage

- name: Upload coverage
  uses: codecov/codecov-action@v3
```

## Debugging Failed Tests

### Verbose Output

```bash
pytest test_focal_point_validation.py -v -s
```

The `-s` flag shows print statements.

### Debug Specific Test

```bash
pytest test_focal_point_validation.py::TestFocalPointValidation::test_focal_point_inside_square -vv
```

### Enable Correa Verbose Mode

```python
import _correa
_correa.set_verbose(True)
```

This shows internal validation steps.

### Use pdb Debugger

```bash
pytest --pdb
```

Drops into debugger on first failure.

## Known Issues and Limitations

None currently. If you discover issues:

1. Check that the `build/` directory contains `_correa*.so`
2. Verify example files exist in `examples/`
3. Ensure pytest >= 7.0.0 is installed
4. Check Python path includes `build/` and `bindings/`

## Performance Considerations

- **Test Duration**: Full suite runs in < 5 seconds on typical hardware
- **Slow Tests**: Mark computationally expensive tests with `@pytest.mark.slow`
- **Fixtures**: Use `scope="session"` for expensive setup that can be reused

## Reporting Issues

If tests fail unexpectedly:

1. Note the test name and error message
2. Check the environment (Python version, OS, build configuration)
3. Run with `-vv` for detailed output
4. Enable verbose mode for Correa-specific debugging
5. Report issues with full error trace

## Future Test Additions

Potential areas for expanded test coverage:

- Performance benchmarks for point-in-polygon algorithm
- Stress tests with very large polygons (10,000+ vertices)
- Randomized property-based testing (using hypothesis)
- Cross-platform validation (Linux, macOS, Windows)
- Memory leak detection (using valgrind or similar)
