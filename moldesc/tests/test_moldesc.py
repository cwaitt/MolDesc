"""
Unit and regression test for the moldesc package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import moldesc


def test_moldesc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "moldesc" in sys.modules
