"""Tests that need a tool skip when it is missing, unless CI says it must be there.

Mark a test with @pytest.mark.needs("cadquery", "gmsh", ...). Tool names are
Python modules, or "ccx" for the CalculiX executable. In CI each job sets
CADFEA_REQUIRE to the tools it installed, so a missing one fails instead of
silently skipping; scripts/check_no_skips.py then checks that every test ran in
at least one job.
"""
import importlib.util
import os
import shutil

import pytest

REQUIRED = {t for t in os.environ.get("CADFEA_REQUIRE", "").split(",") if t}


def have(tool: str) -> bool:
    if tool == "ccx":
        return shutil.which("ccx") is not None
    return importlib.util.find_spec(tool) is not None


def pytest_configure(config):
    config.addinivalue_line("markers", "needs(*tools): skip unless these tools are installed")


def pytest_runtest_setup(item):
    for mark in item.iter_markers(name="needs"):
        missing = [t for t in mark.args if not have(t)]
        if not missing:
            continue
        if REQUIRED & set(missing):
            pytest.fail(f"required tool missing: {', '.join(missing)}")
        pytest.skip(f"not installed here: {', '.join(missing)}")
