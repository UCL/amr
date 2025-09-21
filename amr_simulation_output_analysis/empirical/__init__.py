#!/usr/bin/env python3
"""
Empirical data handling for AMR simulation output analysis.

This package contains:
- data_loader: Loading and caching empirical surveillance data
- normalizers: Name normalization for matching simulation and empirical data
"""

from .data_loader import load_empirical_calibration_data, get_empirical_data_for_plot
from .normalizers import normalize_name_for_empirical_matching

__all__ = [
    'load_empirical_calibration_data', 'get_empirical_data_for_plot',
    'normalize_name_for_empirical_matching'
]