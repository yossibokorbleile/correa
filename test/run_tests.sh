#!/bin/bash
# Test runner script for Correa

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

echo "========================================"
echo "Correa Test Suite"
echo "========================================"
echo ""

# Check if pytest is installed
if ! command -v pytest &> /dev/null; then
    echo "Error: pytest is not installed"
    echo "Install with: pip install pytest pytest-cov"
    exit 1
fi

# Check if the build directory exists
if [ ! -d "$PROJECT_ROOT/build" ]; then
    echo "Error: build directory not found"
    echo "Please build the project first:"
    echo "  mkdir -p build && cd build && cmake .. && cmake --build ."
    exit 1
fi

# Check if _correa module exists
if [ ! -f "$PROJECT_ROOT/build/_correa"*.so ]; then
    echo "Error: _correa module not found in build directory"
    echo "Please build the project first"
    exit 1
fi

echo "Running tests..."
echo ""

# Parse command line arguments
if [ "$1" == "--coverage" ]; then
    echo "Running with coverage report..."
    pytest "$SCRIPT_DIR" -v --cov=correa --cov-report=term --cov-report=html
    echo ""
    echo "Coverage report generated in htmlcov/index.html"
elif [ "$1" == "--fast" ]; then
    echo "Running fast tests only..."
    pytest "$SCRIPT_DIR" -v -m "not slow"
elif [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --coverage    Run tests with coverage report"
    echo "  --fast        Run only fast tests (skip slow tests)"
    echo "  --help, -h    Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                    # Run all tests"
    echo "  $0 --coverage         # Run with coverage"
    echo "  $0 --fast             # Run fast tests only"
    exit 0
else
    pytest "$SCRIPT_DIR" -v "$@"
fi

echo ""
echo "========================================"
echo "Tests completed successfully!"
echo "========================================"
