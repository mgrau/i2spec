"""Shared fixtures."""
import pytest


@pytest.fixture
def fallback_model(monkeypatch):
    """The apps on i2spec2026o (I2SPEC_MODEL), for tests of its own rules or with hand-made line lists."""
    monkeypatch.setenv("I2SPEC_MODEL", "i2spec2026o")
