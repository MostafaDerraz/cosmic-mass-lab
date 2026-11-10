"""Package skeleton (task P1): the engine package imports and reports its version."""

from importlib.metadata import version

import cosmicmass


def test_package_imports():
    assert cosmicmass.__name__ == "cosmicmass"


def test_version():
    assert cosmicmass.__version__ == "0.1.0"


def test_installed_version_matches_source():
    # pyproject.toml reads the version from cosmicmass.__version__, so the
    # installed metadata and the source can never disagree.
    assert version("cosmicmass") == cosmicmass.__version__
