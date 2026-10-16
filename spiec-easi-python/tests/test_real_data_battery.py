"""Slow checks over the real-data battery results (reads the CSVs only, no R).

These verify the *declared* reading rules of docs/VERIFICA_DATI_REALI.md against
the raw results produced by ``tests/real_data/run_battery.py``. They do NOT
re-run R and do NOT touch golden data, tolerances or production code.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

RESULTS = Path(__file__).parent / "real_data" / "results"
CSV_PATH = RESULTS / "battery_results.csv"
VALIDITY_L2 = RESULTS / "validity_L2.csv"

pytestmark = pytest.mark.skipif(
    not CSV_PATH.is_file(),
    reason="battery results not found (run tests/real_data/run_battery.py first)",
)

# block (a) cases -> expected number of seeds completed
EXPECTED = {
    "R1": 10, "R2": 10, "R3": 10, "R4": 10, "R5": 10,
    "L1": 10, "L2": 10, "L3": 5, "L4": 10,
    "S1": 10, "S2": 10, "S3": 10, "S4": 10,
    "P1": 10, "P2": 10, "P3": 10, "R6": 5,
}


def _rows() -> list[dict]:
    with CSV_PATH.open() as fh:
        return list(csv.DictReader(fh))


def test_expected_cases_and_seeds_present() -> None:
    rows = _rows()
    seen: dict[str, set[int]] = {}
    for r in rows:
        seen.setdefault(r["case"], set()).add(int(r["seed"]))
    for case, n in EXPECTED.items():
        assert case in seen, f"{case} missing from results"
        assert len(seen[case]) == n, f"{case}: expected {n} seeds, got {len(seen[case])}"


def test_classification_is_internally_consistent() -> None:
    for r in _rows():
        cls = r["classification"]
        assert cls in {"identico", "equivalente", "diverso"}, r
        if cls == "identico":
            assert int(r["solver_diff_edges"]) == 0, r
            assert r["opt_index_match"] == "1", r
            fd = int(r["final_diff_edges"])
            if fd:
                # only the documented lambda.max exception: at most ONE edge
                # (the argmax|cor| pair) at the least-regularised lambda.
                assert r["opt_index_py"] == "1" and fd <= 1, r
        if cls == "diverso":
            assert r["first_div_step"] != "", r


def test_lambda_selection_and_net_agree_everywhere() -> None:
    # pre-registered observation: the chosen lambda always agrees; the final net
    # agrees except the documented lambda.max exception (at most ONE edge).
    for r in _rows():
        if r["status"] != "ok":
            continue
        assert r["opt_index_match"] == "1", r
        fd = int(r["final_diff_edges"])
        if fd:
            assert r["opt_index_py"] == "1" and fd <= 1, r


def test_l2_is_the_only_solver_divergence() -> None:
    diverging = {r["case"] for r in _rows() if int(r["solver_diff_edges"] or 0) > 0}
    assert diverging <= {"L2"}, diverging


@pytest.mark.skipif(not VALIDITY_L2.is_file(), reason="L2 validity CSV not found")
def test_l2_validity_thresholds() -> None:
    with VALIDITY_L2.open() as fh:
        rows = list(csv.DictReader(fh))
    assert rows, "empty L2 validity results"
    # KKT is satisfied by BOTH solutions on ALL reported divergences
    for r in rows:
        assert float(r["kkt_R"]) <= 1e-3, r
        assert float(r["kkt_Py"]) <= 1e-3, r
    # objective-gap criterion (1e-5) is met by the large majority but not all
    ok = sum(float(r["obj_gap"]) <= 1e-5 for r in rows)
    assert ok / len(rows) > 0.9
