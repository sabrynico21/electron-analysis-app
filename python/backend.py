#!/usr/bin/env python3
"""Single entry point for the analysis backend.

The Electron main process talks to the backend through this file, whether it
runs it with the system `python3` (development / source checkout) or as the
frozen executable produced by PyInstaller (packaged installs).

    backend analyze --params <params.json> --output <result.json>
    backend cluster --input <clustering-input.json>

Both commands live behind a single entry point on purpose: the frozen build
then ships one copy of numpy/scipy/networkx instead of one per script, which
roughly halves the size of the installer.
"""

from __future__ import annotations

import os
import sys

USAGE = "usage: analysis-backend {analyze|cluster} ..."


def _configure_stdio() -> None:
    """Force UTF-8 on stdout/stderr.

    On Windows a *piped* stdout uses the legacy ANSI codepage, so a log line
    containing a character outside cp1252 — the arrows and comparison signs used
    in the analysis log, for instance — raises UnicodeEncodeError and aborts the
    run. Reconfiguring the streams keeps the output identical on every platform.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError, OSError):
            # Streams without a buffer (already wrapped, or closed) are fine.
            pass


def _base_dir() -> str:
    """Directory holding the Python modules, frozen or not."""
    if getattr(sys, "frozen", False):
        # PyInstaller unpacks the bundled modules into _MEIPASS.
        return getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.dirname(os.path.abspath(__file__))


def _prepare_import_path() -> None:
    base = _base_dir()
    for candidate in (base, os.path.join(base, "clustering")):
        if candidate not in sys.path:
            sys.path.insert(0, candidate)


def _run_analyze(argv: list[str]) -> int:
    import analysis

    sys.argv = ["analysis"] + list(argv)
    analysis.main()
    return 0


def _run_cluster(argv: list[str]) -> int:
    from run_clustering import main as cluster_main

    sys.argv = ["run_clustering"] + list(argv)
    return int(cluster_main() or 0)


def main() -> int:
    _configure_stdio()
    argv = sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help"):
        print(USAGE, file=sys.stderr)
        return 2

    command, rest = argv[0], argv[1:]
    if command == "analyze":
        _prepare_import_path()
        return _run_analyze(rest)
    if command == "cluster":
        _prepare_import_path()
        return _run_cluster(rest)

    print(f"{USAGE}\nunknown command: {command}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
