"""The golden fixture tree is content-addressed.

``golden_extended`` stores the large ``solver_adj_path.tsv`` matrices gzipped
(they were ~80 % of the 57 MB tree).  ``CHECKSUMS.sha256`` records the SHA-256 of
the **decompressed** content, so it certifies that the container changed and the
data did not.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

ROOT = Path(__file__).parent / "golden_extended"
CHECKSUMS = ROOT / "CHECKSUMS.sha256"


def _logical_path(rel: str) -> Path:
    """Path on disk holding ``rel`` (``rel`` may live in a ``.gz`` companion)."""
    direct = ROOT / rel
    return direct if direct.is_file() else direct.with_name(direct.name + ".gz")


def _sha256(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as fh:
            data = fh.read()
    else:
        data = path.read_bytes()
    return hashlib.sha256(data).hexdigest()


def _entries() -> list[tuple[str, str]]:
    entries = []
    for line in CHECKSUMS.read_text().splitlines():
        if line.strip():
            digest, _, rel = line.partition("  ")
            entries.append((digest, rel))
    return entries


def test_checksums_cover_every_file() -> None:
    listed = {rel for _, rel in _entries()}
    on_disk = {
        (
            path.relative_to(ROOT).with_name(path.name[: -len(".gz")])
            if path.suffix == ".gz"
            else path.relative_to(ROOT)
        ).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.name != "CHECKSUMS.sha256"
    }
    assert listed == on_disk


def test_golden_content_is_unchanged() -> None:
    """Every fixture still hashes to the value recorded for its logical content."""
    mismatched = []
    missing = []
    for digest, rel in _entries():
        path = _logical_path(rel)
        if not path.is_file():
            missing.append(rel)
        elif _sha256(path) != digest:
            mismatched.append(rel)
    assert not missing, f"missing golden files: {missing}"
    assert not mismatched, f"golden files whose content changed: {mismatched}"
