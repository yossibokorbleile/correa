"""
Pytest configuration and shared fixtures for Correa tests.
"""

import sys
import os
import pytest

# Add the build directory and bindings to the path for all tests
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bindings'))


@pytest.fixture(scope="session")
def examples_dir():
    """Return the path to the examples directory."""
    return os.path.join(os.path.dirname(__file__), '..', 'examples')


@pytest.fixture(scope="session")
def test_data_dir():
    """Return the path to the test data directory."""
    return os.path.join(os.path.dirname(__file__), 'data')


@pytest.fixture(autouse=True)
def reset_verbose():
    """Reset verbose mode after each test."""
    yield
    # Clean up: ensure verbose is turned off after each test
    try:
        import _correa
        _correa.set_verbose(False)
    except ImportError:
        pass


def pytest_configure(config):
    """Pytest configuration hook."""
    # Add custom markers
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "requires_examples: marks tests that require example data files"
    )
