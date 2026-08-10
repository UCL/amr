#!/usr/bin/env python3
"""
Base plotting classes and utilities for AMR simulation output analysis.

This module provides shared plotting infrastructure for consistent styling,
validation, saving, and cleanup.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, Tuple, List
from abc import ABC, abstractmethod
import logging

from ..config import PlotConfig
from ..utils import (
    safe_plot_creation, ensure_output_directory, save_and_show_plot,
    validate_plot_data, add_statistics_text_box
)

logger = logging.getLogger(__name__)

class BasePlot(ABC):
    """
    Abstract base class for all AMR simulation plots.
    
    Provides standardized infrastructure for:
    - Plot configuration and styling
    - Data validation
    - Provenance-controlled comparison overlay
    - Error handling and cleanup
    - Consistent saving and display
    """
    
    def __init__(self, 
                 data: pd.DataFrame,
                 plot_config: PlotConfig,
                 plotting_config: PlotConfig,
                 empirical_data: Optional[Dict[str, pd.DataFrame]] = None):
        """
        Initialize base plot.
        
        Args:
            data: Simulation data DataFrame
            plot_config: Configuration for this specific plot type
            plotting_config: General plotting configuration
            empirical_data: Optional empirical data for overlay
        """
        self.data = data
        self.plot_config = plot_config
        self.plotting_config = plotting_config
        self.empirical_data = empirical_data or {}
        
        # Set up plotting style
        plt.style.use(self.plotting_config.style)
        
    @abstractmethod
    def get_required_columns(self) -> List[str]:
        """Return list of columns required for this plot."""
        pass
    
    @abstractmethod
    def create_plot_content(self, fig, axes) -> None:
        """Create the actual plot content. Must be implemented by subclasses."""
        pass
    
    @abstractmethod
    def get_plot_title(self) -> str:
        """Return the main title for this plot."""
        pass
    
    def validate_data(self) -> bool:
        """Validate that data is suitable for this plot."""
        required_columns = self.get_required_columns()
        return validate_plot_data(self.data, required_columns, self.get_plot_title())
    
    def setup_figure(self, figsize: Optional[Tuple[int, int]] = None) -> Tuple[plt.Figure, plt.Axes]:
        """Set up matplotlib figure and axes."""
        if figsize is None:
            figsize = self.plotting_config.figure_size_single
            
        fig, ax = plt.subplots(figsize=figsize)
        return fig, ax
    
    def setup_subplots(self, 
                      nrows: int = 2, 
                      ncols: int = 2,
                      figsize: Optional[Tuple[int, int]] = None) -> Tuple[plt.Figure, np.ndarray]:
        """Set up matplotlib figure with subplots."""
        if figsize is None:
            figsize = self.plotting_config.adaptive_figure_size
            
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        if isinstance(axes, plt.Axes):
            axes = np.array([axes])
        return fig, axes.flatten() if axes.ndim > 1 else axes
    
    def add_empirical_overlay(self, ax: plt.Axes, 
                             bacteria: str = None, 
                             drug: str = None,
                             data_source: str = None) -> None:
        """
        Add an eligible comparison overlay to a plot if available.
        
        Args:
            ax: Matplotlib axis to add overlay to
            bacteria: Bacteria name for filtering empirical data
            drug: Drug name for filtering empirical data  
            data_source: Type of empirical data ('drug_failure', 'mic_values', etc.)
        """
        if not self.plot_config.empirical_overlay:
            return
            
        if data_source not in self.empirical_data:
            return
            
        empirical_df = self.empirical_data[data_source]
        if empirical_df is None:
            return
            
        try:
            # Import here to avoid circular imports
            from ..empirical.data_loader import get_empirical_data_for_plot
            
            sim_years, means, p5, p95 = get_empirical_data_for_plot(
                empirical_df, 
                drug=drug, 
                bacteria=bacteria,
                data_source=data_source,
                include_best_guess_placeholders=(
                    self.plot_config.show_best_guess_placeholder_overlays
                ),
            )
            
            if sim_years is not None and means is not None:
                overlay_label = (
                    'Best-guess placeholder'
                    if self.plot_config.show_best_guess_placeholder_overlays
                    else 'Observed comparison'
                )
                ax.scatter(sim_years, means, 
                          color='red', alpha=0.7, s=30, 
                          label=overlay_label, zorder=10)
                
                # Add confidence intervals if available
                if p5 is not None and p95 is not None:
                    ax.fill_between(sim_years, p5, p95, 
                                   color='red', alpha=0.2, 
                                   label=f'{overlay_label} interval')
                
                logger.info(f"Added comparison overlay: {len(sim_years)} points")
                
        except Exception as e:
            logger.warning(f"Could not add comparison overlay: {e}")
    
    def finalize_plot(self, fig: plt.Figure, output_filename: str) -> None:
        """Finalize plot with consistent styling and save."""
        # Add main title
        fig.suptitle(self.get_plot_title(), fontsize=16)
        
        # Tight layout
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for title
        
        # Save plot
        output_path = Path(self.plot_config.output_dir) / output_filename
        save_and_show_plot(
            output_path,
            title=self.get_plot_title(),
            dpi=self.plotting_config.dpi,
            bbox_inches=self.plotting_config.bbox_inches,
            show_plot=self.plot_config.show_plot
        )
    
    @safe_plot_creation
    def create(self, output_filename: str) -> Optional[plt.Figure]:
        """
        Main method to create the plot.
        
        Args:
            output_filename: Name of output file
            
        Returns:
            Figure object if successful, None if failed
        """
        if not self.validate_data():
            return None
            
        # Ensure output directory exists
        ensure_output_directory(Path(self.plot_config.output_dir))
        
        # Create figure and content
        fig, axes = self.setup_figure()
        self.create_plot_content(fig, axes)
        
        # Finalize and save
        self.finalize_plot(fig, output_filename)
        
        return fig

class StandardizedPlot(BasePlot):
    """
    Concrete implementation for standardized single-panel plots.
    
    Useful for simple time series and other common plot types.
    """
    
    def __init__(self, 
                 data: pd.DataFrame,
                 plot_config: PlotConfig,
                 plotting_config: PlotConfig,
                 title: str,
                 required_columns: List[str],
                 plot_function: callable,
                 empirical_data: Optional[Dict[str, pd.DataFrame]] = None):
        """
        Initialize standardized plot.
        
        Args:
            data: Simulation data
            plot_config: Plot configuration
            plotting_config: General plotting config
            title: Plot title
            required_columns: Required data columns
            plot_function: Function that creates plot content (takes fig, ax, data)
            empirical_data: Optional empirical data
        """
        super().__init__(data, plot_config, plotting_config, empirical_data)
        self.title = title
        self.required_columns = required_columns
        self.plot_function = plot_function
    
    def get_required_columns(self) -> List[str]:
        return self.required_columns
    
    def get_plot_title(self) -> str:
        return self.title
    
    def create_plot_content(self, fig: plt.Figure, ax: plt.Axes) -> None:
        """Create plot content using the provided plot function."""
        self.plot_function(fig, ax, self.data)

class GroupedPlot(BasePlot):
    """
    Base class for grouped plots with multiple subplots.
    """
    
    def __init__(self,
                 data: pd.DataFrame,
                 plot_config: PlotConfig, 
                 plotting_config: PlotConfig,
                 empirical_data: Optional[Dict[str, pd.DataFrame]] = None):
        super().__init__(data, plot_config, plotting_config, empirical_data)
    
    def setup_figure(self, figsize: Optional[Tuple[int, int]] = None) -> Tuple[plt.Figure, np.ndarray]:
        """Override to return subplots by default."""
        return self.setup_subplots(figsize=figsize)
    
    @abstractmethod
    def create_subplot_content(self, ax: plt.Axes, subplot_index: int) -> None:
        """Create content for a specific subplot. Must be implemented by subclasses.""" 
        pass
    
    def create_plot_content(self, fig: plt.Figure, axes: np.ndarray) -> None:
        """Create content for all subplots."""
        for i, ax in enumerate(axes):
            try:
                self.create_subplot_content(ax, i)
            except Exception as e:
                logger.warning(f"Failed to create subplot {i}: {e}")
                # Create placeholder for failed subplot
                ax.text(0.5, 0.5, f'Subplot {i+1}\n(Error: {str(e)})', 
                       ha='center', va='center', fontsize=12, color='red')
                ax.set_axis_off()

def create_simple_time_series_plot(fig: plt.Figure, ax: plt.Axes, data: pd.DataFrame,
                                  x_col: str, y_col: str, 
                                  title: str, ylabel: str,
                                  color: str = 'blue',
                                  smoothing_window: int = 1095) -> None:
    """
    Helper function to create a simple time series plot.
    
    Args:
        fig: Figure object
        ax: Axis object  
        data: Data DataFrame
        x_col: X-axis column name
        y_col: Y-axis column name
        title: Plot title
        ylabel: Y-axis label
        color: Line color
        smoothing_window: Rolling window for smoothing
    """
    if y_col not in data.columns:
        ax.text(0.5, 0.5, f'Data not available\n(missing {y_col})', 
               ha='center', va='center', fontsize=12, color='gray')
        ax.set_title(title)
        ax.set_axis_off()
        return
    
    # Apply smoothing
    smoothed_data = data[y_col].rolling(
        window=min(smoothing_window, len(data)), 
        min_periods=1, 
        center=True
    ).mean()
    
    # Create plot
    ax.plot(data[x_col], smoothed_data, color=color, linewidth=2)
    ax.set_title(title)
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    
    # Add statistics text box
    add_statistics_text_box(ax, smoothed_data)
