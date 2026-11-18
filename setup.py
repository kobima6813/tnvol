# SPDX-License-Identifier: AGPL-3.0-only
"""Compile the native codec as a plain shared library (ctypes, no Python.h)."""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
import sysconfig
import tempfile

from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
from setuptools.extension import Extension

def _is_msvc(cc: str) -> bool:
    base = os.path.basename(cc).lower()
    return base in ("cl", "cl.exe") or base.endswith("cl.exe")


def _detect_cc() -> str:
    env = os.environ.get("CC")
    if env:
        return env
    if sys.platform == "win32":
        for cand in ("cl", "clang-cl", "gcc", "clang", "cc"):
            if shutil.which(cand):
                return cand
        return "cl"
    for cand in ("cc", "gcc", "clang"):
        if shutil.which(cand):
            return cand
    return "cc"


def _target_is_x86_64() -> bool:
    archflags = os.environ.get("ARCHFLAGS", "")
    if "-arch x86_64" in archflags:
        return True
    if "-arch arm64" in archflags or "-arch aarch64" in archflags:
        return False
    plat = (sysconfig.get_platform() or "").lower()
    if "x86_64" in plat or "amd64" in plat or plat.endswith("x64"):
        return True
    if "arm64" in plat or "aarch64" in plat:
        return False
    # cibuildwheel linux cross: AUDITWHEEL_ARCH / _PYTHON_HOST_PLATFORM
    for key in ("AUDITWHEEL_ARCH", "CBUILDWHEEL_ARCH", "_PYTHON_HOST_PLATFORM"):
        val = (os.environ.get(key) or "").lower()
        if "x86_64" in val or "amd64" in val:
            return True
        if "aarch64" in val or "arm64" in val:
            return False
    mach = platform.machine().lower()
    return mach in ("x86_64", "amd64", "x64")


def _lib_names() -> tuple[str, str]:
    if sys.platform == "darwin":
        return "libwlc.dylib", "libwlc_avx2.dylib"
    if sys.platform == "win32":
        return "libwlc.dll", "libwlc_avx2.dll"
    return "libwlc.so", "libwlc_avx2.so"


def _compile_one(cc: str, src: str, dest: str, avx2: bool) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(dest)) or ".", exist_ok=True)
    msvc = _is_msvc(cc)
    if msvc:
        # /fp:precise matches -ffp-contract=off (no mul+add contraction).
        flags = ["/nologo", "/O2", "/fp:precise", "/DNPAIR=16", "/LD"]
        if avx2:
            flags.append("/arch:AVX2")
        # Keep .obj/.exp/.lib out of the package directory.
        with tempfile.TemporaryDirectory() as td:
            dll_name = os.path.basename(dest)
            obj = os.path.join(td, "wlc.obj")
            out = os.path.join(td, dll_name)
            cmd = [cc, *flags, src, f"/Fo:{obj}", f"/Fe:{out}"]
            print("building", dll_name + ":", " ".join(cmd), flush=True)
            subprocess.check_call(cmd)
            shutil.copy2(out, dest)
        return
    flags = ["-O3", "-ffp-contract=off", "-DNPAIR=16", "-fPIC"]
    if avx2:
        flags.extend(["-mavx2", "-mfma"])
    if sys.platform == "darwin":
        archflags = os.environ.get("ARCHFLAGS", "")
        if "-arch x86_64" in archflags:
            flags.extend(["-arch", "x86_64"])
        elif "-arch arm64" in archflags:
            flags.extend(["-arch", "arm64"])
    cmd = [cc, *flags, "-shared", "-o", dest, src]
    if sys.platform != "win32":
        cmd.append("-lm")
    print("building", os.path.basename(dest) + ":", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def compile_libwlc(dest_dir: str) -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    src = os.path.join(here, "src", "tnvol", "wlc.c")
    cc = _detect_cc()
    portable, avx2 = _lib_names()
    _compile_one(cc, src, os.path.join(dest_dir, portable), avx2=False)
    if _target_is_x86_64():
        _compile_one(cc, src, os.path.join(dest_dir, avx2), avx2=True)


class BuildPy(build_py):
    def run(self):
        super().run()
        compile_libwlc(os.path.join(self.build_lib, "tnvol"))


class BuildExt(build_ext):
    """Also used by editable installs, which do not always run build_py hooks."""

    def run(self):
        if getattr(self, "inplace", False):
            here = os.path.dirname(os.path.abspath(__file__))
            compile_libwlc(os.path.join(here, "src", "tnvol"))
        elif self.build_lib:
            compile_libwlc(os.path.join(self.build_lib, "tnvol"))


# A dummy extension so pip/setuptools run build_ext (needed for editable installs).
# The empty sources list means the default compiler is never invoked.
ext = Extension("tnvol._native_marker", sources=[])

setup(
    ext_modules=[ext],
    cmdclass={"build_py": BuildPy, "build_ext": BuildExt},
)
