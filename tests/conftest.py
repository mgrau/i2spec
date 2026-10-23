"""Shared fixtures."""
import pytest


def pytest_addoption(parser):
    parser.addoption("--slow", action="store_true", help="also run the slow refits of the fit pipeline")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: a full refit (minutes); run with --slow")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--slow"):
        return
    skip = pytest.mark.skip(reason="slow: run with --slow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)


@pytest.fixture
def hand_made(monkeypatch):
    """Hand-made line lists taken as they are: the model's level tables neither move their lines nor give them
    uncertainties, so no table has to be built and the lines stay where the test put them."""
    from i2spec import lookup
    monkeypatch.setattr(lookup, "preferred_positions", lambda master: master)
    monkeypatch.setattr(lookup, "_held", lambda line: None)
