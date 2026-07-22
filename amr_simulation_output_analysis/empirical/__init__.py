#!/usr/bin/env python3
"""
Provenance-controlled comparison overlays for AMR simulation analysis.

This package contains:
- data_loader: loading observed comparisons and optional best-guess placeholders
- provenance: strict observed-versus-placeholder classification
- normalizers: name normalization for matching simulation and overlay data
"""

from .data_loader import load_empirical_calibration_data, get_empirical_data_for_plot
from .normalizers import normalize_name_for_empirical_matching
from .provenance import filter_overlay_rows

__all__ = [
    'load_empirical_calibration_data', 'get_empirical_data_for_plot',
    'normalize_name_for_empirical_matching', 'filter_overlay_rows'
]
