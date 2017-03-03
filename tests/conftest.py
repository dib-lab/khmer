import pytest


@pytest.fixture()
def seed():
    """Base seed used by all tests that require random numbers"""
    return 1
