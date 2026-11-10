"""P2: the pinned environment is written in two files; they must agree (decision D13)."""

import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def pins(lines):
    return sorted(line.strip() for line in lines if line.strip() and not line.strip().startswith("#"))


def test_pyproject_dependencies_equal_requirements_txt():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    requirements = (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
    assert pins(project["dependencies"]) == pins(requirements)


def test_every_dependency_is_pinned_exactly():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    extras = [d for group in project["optional-dependencies"].values() for d in group]
    assert all("==" in d for d in project["dependencies"] + extras)
