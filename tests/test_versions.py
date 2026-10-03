"""Every version string agrees: pyproject.toml, the package, CITATION.cff, CMakeLists.txt and the newest CHANGELOG heading."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _find(path, pattern):
    m = re.search(pattern, (ROOT / path).read_text(), re.M)
    assert m, f"no version in {path}"
    return m.group(1).strip()


def test_versions_agree():
    import nmc_model
    versions = {
        "pyproject.toml": _find("pyproject.toml", r'^version = "(.+)"'),
        "__version__": nmc_model.__version__,
        "CITATION.cff": _find("CITATION.cff", r"^version: (.+)$"),
        "CMakeLists.txt": _find("CMakeLists.txt", r"^project\(\w+ VERSION ([\d.]+)"),
        "CHANGELOG.md": _find("CHANGELOG.md", r"^## ([\d.]+)"),
    }
    assert len(set(versions.values())) == 1, versions
