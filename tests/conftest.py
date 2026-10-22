"""Shared fixtures."""
import pytest


@pytest.fixture
def hand_made(monkeypatch):
    """Hand-made line lists taken as they are: the model's level tables neither move their lines nor give them
    uncertainties, so no table has to be built and the lines stay where the test put them."""
    from i2spec import lookup
    monkeypatch.setattr(lookup, "preferred_positions", lambda master: master)
    monkeypatch.setattr(lookup, "_held", lambda line: None)
