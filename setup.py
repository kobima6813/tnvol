# SPDX-License-Identifier: AGPL-3.0-only
"""Compile the native codec as a plain shared library (ctypes, no Python.h)."""
from __future__ import annotations

import os
import subprocess
import sys

from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
from setuptools.extension import Extension

WLC_FLAGS = [
    "-O3",
    "-mavx2",
    "-mfma",
    "-ffp-contract=off",
    "-DNPAIR=16",  # required to match the published r10 numbers
    "-fPIC",
]


def compile_libwlc(dest: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(dest)) or ".", exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    src = os.path.join(here, "src", "tnvol", "wlc.c")
    cc = os.environ.get("CC") or "gcc"
    cmd = [cc, *WLC_FLAGS, "-shared", "-o", dest, src, "-lm"]
    print("building libwlc:", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


class BuildPy(build_py):
    def run(self):
        super().run()
        compile_libwlc(os.path.join(self.build_lib, "tnvol", "libwlc.so"))


class BuildExt(build_ext):
    """Also used by editable installs, which do not always run build_py hooks."""

    def run(self):
        if getattr(self, "inplace", False):
            here = os.path.dirname(os.path.abspath(__file__))
            compile_libwlc(os.path.join(here, "src", "tnvol", "libwlc.so"))
        elif self.build_lib:
            compile_libwlc(os.path.join(self.build_lib, "tnvol", "libwlc.so"))


if sys.platform != "linux":
    sys.stderr.write(
        "warning: tnvol's native codec is currently supported on Linux x86-64 with AVX2\n"
    )

# A dummy extension so pip/setuptools run build_ext (needed for editable installs).
# The empty sources list means the default compiler is never invoked.
ext = Extension("tnvol._native_marker", sources=[])

setup(
    ext_modules=[ext],
    cmdclass={"build_py": BuildPy, "build_ext": BuildExt},
)
