"""Setuptools build for the optional compiled MB core.

``spieceasi._mb_core`` is a near-literal translation of ``huge_core.cpp::mb``
(see ``src/spieceasi/_mb_core.cpp``). It is an *accelerator*: the pure-Python
solver in ``mb_solver_reference.py`` remains the reference implementation and the
automatic fallback when the extension cannot be imported.

Floating-point flags are chosen to match R's own package build (``-O2``, no
fast-math) so that results stay bit-comparable:

* ``-O2``            — same optimisation level as R/huge;
* ``-ffp-contract=off`` — never fuse multiply-add (FMA); an FMA would change the
  rounding of the scalar accumulation loops and make results compiler-dependent;
* no ``-ffast-math``, no ``-march=native`` — these reorder/simplify floating point
  and would break parity.

On MSVC the equivalent is ``/O2 /fp:precise`` (MSVC never contracts by default).
"""

from __future__ import annotations

import sys
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

_POSIX_FLAGS = ["-O2", "-ffp-contract=off", "-fno-lto"]
_MSVC_FLAGS = ["/O2", "/fp:precise"]

ext_modules = [
    Pybind11Extension(
        "spieceasi._mb_core",
        [str(Path("src") / "spieceasi" / "_mb_core.cpp")],
        cxx_std=17,
        extra_compile_args=_MSVC_FLAGS if sys.platform == "win32" else _POSIX_FLAGS,
    )
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
