import pytest

def pytest_addoption(parser):
    parser.addoption("--run-config", default=None,
                     help="Path to integration test config JSON")