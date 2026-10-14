"""The golden fixture trees are content-addressed.

Every golden tree ships a ``CHECKSUMS.sha256`` listing the SHA-256 of each
fixture, so a fixture cannot change without a test failing:

* ``golden_extended`` stores the large ``solver_adj_path.tsv`` matrices gzipped
  (they were ~80 % of the 57 MB tree); the recorded digest is that of the
  **decompressed** content, so it certifies that the container changed and the
  data did not;
* ``golden``, ``golden_app`` and ``golden_stats`` are stored uncompressed and
  their digest is that of the file on disk.

Note on scope: a manifest proves that nothing changed **from the moment the
manifest was created onwards**. It says nothing about the history of a file
before that (for the history, use ``git log`` / ``git diff`` against a tag).
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pytest

HERE = Path(__file__).parent
TREES = ["golden", "golden_app", "golden_stats", "golden_extended"]
MANIFEST = "CHECKSUMS.sha256"


def _root(tree: str) -> Path:
    return HERE / tree


def _logical_path(root: Path, rel: str) -> Path:
    """Path on disk holding ``rel`` (``rel`` may live in a ``.gz`` companion)."""
    direct = root / rel
    return direct if direct.is_file() else direct.with_name(direct.name + ".gz")


def _sha256(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as fh:
            data = fh.read()
    else:
        data = path.read_bytes()
    return hashlib.sha256(data).hexdigest()


def _entries(root: Path) -> list[tuple[str, str]]:
    entries = []
    for line in (root / MANIFEST).read_text().splitlines():
        if line.strip():
            digest, _, rel = line.partition("  ")
            entries.append((digest, rel))
    return entries


@pytest.mark.parametrize("tree", TREES)
def test_manifest_exists_and_is_not_empty(tree: str) -> None:
    root = _root(tree)
    assert (root / MANIFEST).is_file(), f"{tree}/{MANIFEST} is missing"
    assert _entries(root), f"{tree}/{MANIFEST} has no entries"


@pytest.mark.parametrize("tree", TREES)
def test_checksums_cover_every_file(tree: str) -> None:
    root = _root(tree)
    listed = {rel for _, rel in _entries(root)}
    on_disk = {
        (
            path.relative_to(root).with_name(path.name[: -len(".gz")])
            if path.suffix == ".gz"
            else path.relative_to(root)
        ).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path.name != MANIFEST
    }
    assert listed == on_disk


@pytest.mark.parametrize("tree", TREES)
def test_golden_content_is_unchanged(tree: str) -> None:
    """Every fixture still hashes to the value recorded for its logical content."""
    root = _root(tree)
    mismatched = []
    missing = []
    for digest, rel in _entries(root):
        path = _logical_path(root, rel)
        if not path.is_file():
            missing.append(rel)
        elif _sha256(path) != digest:
            mismatched.append(rel)
    assert not missing, f"{tree}: missing golden files: {missing}"
    assert not mismatched, f"{tree}: golden files whose content changed: {mismatched}"
