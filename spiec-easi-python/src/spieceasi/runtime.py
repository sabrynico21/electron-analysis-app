"""Runtime helpers: worker-count defaults based on *physical* cores.

The StARS subsamples are CPU-bound and independent, and the compiled MB core
releases the GIL, so a **thread** pool is the right backend: no pickling, no
re-execution of the interpreter (which is what breaks a frozen build) and
bit-identical results for any worker count.

How many workers to use by default is decided from the number of *physical*
cores, never the logical ones: with SMT the logical count oversubscribes the
floating-point units and slows the pipeline down.
"""

from __future__ import annotations

import os
import platform
import subprocess


def _linux_physical_cores() -> int | None:
    """Count distinct (physical id, core id) pairs in ``/proc/cpuinfo``."""
    try:
        pairs: set[tuple[str, str]] = set()
        current: dict[str, str] = {}
        with open("/proc/cpuinfo", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip():
                    if "physical id" in current and "core id" in current:
                        pairs.add((current["physical id"], current["core id"]))
                    current = {}
                    continue
                key, _, value = line.partition(":")
                key = key.strip()
                if key in ("physical id", "core id"):
                    current[key] = value.strip()
        if "physical id" in current and "core id" in current:
            pairs.add((current["physical id"], current["core id"]))
        return len(pairs) or None
    except OSError:
        return None


def _macos_physical_cores() -> int | None:
    try:
        done = subprocess.run(
            ["sysctl", "-n", "hw.physicalcpu"],
            capture_output=True,
            text=True,
            timeout=5,
            check=True,
        )
        value = int(done.stdout.strip())
        return value if value > 0 else None
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def _windows_physical_cores() -> int | None:
    commands = (
        ["wmic", "cpu", "get", "NumberOfCores", "/value"],
        [
            "powershell",
            "-NoProfile",
            "-Command",
            "(Get-CimInstance Win32_Processor | "
            "Measure-Object -Property NumberOfCores -Sum).Sum",
        ],
    )
    for command in commands:
        try:
            done = subprocess.run(
                command, capture_output=True, text=True, timeout=10, check=False
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if done.returncode != 0:
            continue
        total = 0
        for line in done.stdout.splitlines():
            line = line.strip()
            if line.startswith("NumberOfCores="):
                total += int(line.split("=", 1)[1])
            elif line.isdigit():
                total += int(line)
        if total > 0:
            return total
    return None


def physical_core_count() -> int | None:
    """Number of physical CPU cores, or ``None`` if it cannot be determined."""
    try:  # optional dependency, used only if the app happens to ship it
        import psutil
    except ImportError:
        pass
    else:
        try:
            value = psutil.cpu_count(logical=False)
        except Exception:  # noqa: BLE001 - any psutil failure means "unknown"
            value = None
        if value:
            return int(value)

    system = platform.system()
    if system == "Linux":
        return _linux_physical_cores()
    if system == "Darwin":
        return _macos_physical_cores()
    if system == "Windows":
        return _windows_physical_cores()
    return None


def default_n_jobs() -> int:
    """Prudent default worker count for the StARS subsamples.

    * physical cores when they can be detected;
    * otherwise assume SMT (two hardware threads per core) and use
      ``logical // 2`` — conservative, because over-subscribing is worse than
      under-using;
    * ``1`` (serial) when even the logical count is unknown.

    Never less than 1.  Callers can always override it (the desktop app exposes
    ``spieceasiJobs``).
    """
    physical = physical_core_count()
    if physical:
        return max(1, int(physical))
    logical = os.cpu_count()
    if logical:
        return max(1, int(logical) // 2)
    return 1
