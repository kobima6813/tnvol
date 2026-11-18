# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Yackov Mozgori
"""tnvol: error-bounded compressor for 3-D scientific volumes."""

from .codec import compress, decompress, encode, decode, info

__version__ = "0.10.0"
__all__ = ["compress", "decompress", "encode", "decode", "info", "__version__"]
