"""Comparisons against the private legacy oracle.

These read the original model and its archived outputs from the directory named by
the NMC_ORACLE_DIR environment variable. That material is never distributed with
this repository, so every test here is skipped when the variable is unset.
"""
import os
from pathlib import Path

import pytest

ORACLE_DIR = os.environ.get("NMC_ORACLE_DIR")


def pytest_collection_modifyitems(config, items):
    if ORACLE_DIR and Path(ORACLE_DIR).is_dir():
        return
    skip = pytest.mark.skip(reason="NMC_ORACLE_DIR not set; private oracle comparisons skipped")
    for item in items:
        if "private_checks" in str(item.fspath):
            item.add_marker(skip)


@pytest.fixture
def oracle_dir() -> Path:
    return Path(ORACLE_DIR)


