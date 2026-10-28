#!/usr/bin/env bash
# Regenerate R references, run all three Python-vs-R suites, and seal the output
# with checksums. This writes only to a new, caller-selected output directory.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-}"
if [[ -z "$OUT" ]]; then
  echo "Usage: bash scripts/reproduce_equivalence.sh <new-output-directory>" >&2
  exit 2
fi
if [[ -e "$OUT" && -n "$(find "$OUT" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Refusing to overwrite non-empty output directory: $OUT" >&2
  exit 2
fi
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

PY_CCLASSO="${PY_CCLASSO:-python3}"
PY_SPIEC="${PY_SPIEC:-python3}"
PY_SPRING="${PY_SPRING:-python3}"

run_component() {
  local name="$1" python_bin="$2" package="$3"
  local out_dir="$OUT/$name"
  mkdir -p "$out_dir"
  echo "== $name: R golden generation =="
  case "$name" in
    cclasso)
      (cd "$ROOT/$package" && Rscript tests/r_reference/generate_golden.R "$out_dir")
      ;;
    spieceasi)
      (cd "$ROOT/$package" && GOLDEN_OUT="$out_dir" \
        Rscript tests/r_reference/generate_golden.R)
      ;;
    spring)
      (cd "$ROOT/$package" && Rscript tests/r_reference/generate_golden.R "$out_dir")
      ;;
  esac
  echo "== $name: environment manifest =="
  Rscript "$ROOT/scripts/write_r_environment.R" "$out_dir/r_environment.json"
  "$python_bin" "$ROOT/scripts/write_python_environment.py" \
    "$out_dir/python_environment.json" --component "$name"
  echo "== $name: Python equivalence suite =="
  case "$name" in
    cclasso)
      (cd "$ROOT/cclasso-python" && CCLASSO_GOLDEN_DIR="$out_dir" \
        "$python_bin" -m pytest -q)
      ;;
    spieceasi)
      (cd "$ROOT/spiec-easi-python" && SPIEC_EASI_GOLDEN_DIR="$out_dir" \
        PYTHONPATH=src "$python_bin" -m pytest -q tests/test_equivalence.py)
      ;;
    spring)
      (cd "$ROOT/spring-python" && SPRING_GOLDEN_DIR="$out_dir" \
        "$python_bin" -m pytest -q)
      ;;
  esac
  (cd "$out_dir" && find . -type f ! -name MANIFEST.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > MANIFEST.sha256)
}

run_component cclasso "$PY_CCLASSO" cclasso-python
run_component spieceasi "$PY_SPIEC" spiec-easi-python
run_component spring "$PY_SPRING" spring-python

echo "Reproduction outputs and manifests: $OUT"
