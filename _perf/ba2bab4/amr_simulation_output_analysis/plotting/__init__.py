#!/usr/bin/env python3
"""
Plotting modules for AMR simulation output analysis.

This package contains:
- base_plots: Base classes and common plotting utilities
- grouped_plots: Grouped summary plots
- detail_plots: Granular analysis plots saved to subfolders
"""

from .base_plots import BasePlot, StandardizedPlot
from .grouped_plots import create_grouped_plots
from .detail_plots import create_detail_plots

__all__ = [
    'BasePlot', 'StandardizedPlot', 
    'create_grouped_plots', 'create_detail_plots'
]
