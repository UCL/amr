#!/usr/bin/env python3
"""
Utility Functions for AMR Simulation Output Analysis

This module provides common utility functions including error handling,
logging setup, and data validation functions.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from pathlib import Path
from typing import Optional, Union, List, Dict, Any, Iterable
from functools import wraps

def setup_logging(log_level: str = "INFO", log_file: Optional[str] = None) -> logging.Logger:
    """
    Set up logging for the analysis pipeline.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional log file path
        
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger('amr_analysis')
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Clear existing handlers
    logger.handlers.clear()
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s'
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    
    return logger

def safe_divide(numerator: Union[np.ndarray, pd.Series, float],
               denominator: Union[np.ndarray, pd.Series, float],
               default: float = np.nan) -> Union[np.ndarray, pd.Series, float]:
    """Safe division avoiding division by zero and suppressing runtime warnings."""
    numerator_arr = np.asarray(numerator, dtype=float)
    denominator_arr = np.asarray(denominator, dtype=float)
    result = np.full_like(numerator_arr, default, dtype=float)

    with np.errstate(divide='ignore', invalid='ignore'):
        np.divide(numerator_arr, denominator_arr, out=result, where=denominator_arr != 0)

    return result

def coerce_policy_identifier(policy_value: Any) -> Optional[int]:
    """Attempt to convert a policy identifier (int/float/str) into an integer label."""
    if policy_value is None:
        return None

    if isinstance(policy_value, (int, np.integer)):
        return int(policy_value)

    try:
        numeric = int(float(policy_value))
        return numeric
    except (ValueError, TypeError):
        pass

    policy_str = str(policy_value).lower()
    digits = ''.join(ch for ch in policy_str if ch.isdigit() or ch == '-')
    if digits not in ('', '-'):
        try:
            return int(digits)
        except ValueError:
            return None

    return None

def normalize_policy_identifier_list(policy_values: Optional[Union[int, float, str, Iterable[Any]]]) -> Optional[List[int]]:
    """Normalize a user-provided list (or scalar) of policy identifiers into distinct ints."""
    if policy_values is None:
        return None

    if isinstance(policy_values, (int, float, str)):
        iterable: Iterable[Any] = [policy_values]
    else:
        try:
            iterable = list(policy_values)
        except TypeError:
            iterable = [policy_values]

    normalized: List[int] = []
    for value in iterable:
        numeric = coerce_policy_identifier(value)
        if numeric is None or numeric in normalized:
            continue
        normalized.append(numeric)

    return normalized or None

def validate_dataframe(df: pd.DataFrame, required_columns: List[str] = None) -> bool:
    """
    Validate that a DataFrame meets basic requirements.
    
    Args:
        df: DataFrame to validate
        required_columns: List of required column names
        
    Returns:
        True if valid, False otherwise
    """
    if df is None:
        return False
        
    if df.empty:
        logging.warning("DataFrame is empty")
        return False
        
    if required_columns:
        missing_cols = [col for col in required_columns if col not in df.columns]
        if missing_cols:
            logging.warning(f"Missing required columns: {missing_cols}")
            return False
    
    return True

def validate_plot_data(df: pd.DataFrame, 
                      required_columns: List[str],
                      plot_name: str = "plot") -> bool:
    """
    Validate data before creating a plot.
    
    Args:
        df: DataFrame with plot data
        required_columns: Columns required for the plot
        plot_name: Name of the plot for error messages
        
    Returns:
        True if data is valid for plotting, False otherwise
    """
    if not validate_dataframe(df, required_columns):
        logging.warning(f"Cannot create {plot_name}: data validation failed")
        return False
        
    # Check for sufficient data
    if len(df) < 2:
        logging.warning(f"Cannot create {plot_name}: insufficient data points ({len(df)})")
        return False
        
    # Check for all-zero or all-NaN columns
    for col in required_columns:
        if col in df.columns:
            values = df[col].dropna()
            if len(values) == 0 or values.sum() == 0:
                logging.warning(f"Cannot create {plot_name}: column '{col}' has no valid data")
                return False
    
    return True

def safe_plot_creation(func):
    """
    Decorator for safe plot creation with error handling and memory management.
    
    Handles common plotting errors gracefully and ensures matplotlib resources
    are properly cleaned up.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        plot_name = func.__name__
        logger = logging.getLogger('amr_analysis')
        
        try:
            logger.info(f"Creating plot: {plot_name}")
            result = func(*args, **kwargs)
            logger.info(f"Successfully created plot: {plot_name}")
            return result
            
        except Exception as e:
            logger.error(f"Error creating plot {plot_name}: {str(e)}")
            print(f"[WARNING] Failed to create {plot_name}: {str(e)}")
            
            # Clean up any open matplotlib figures
            plt.close('all')
            return None
            
    return wrapper

def ensure_output_directory(output_path: Union[str, Path]) -> Path:
    """
    Ensure output directory exists, creating it if necessary.
    
    Args:
        output_path: Path to output file or directory
        
    Returns:
        Path object for the directory
    """
    path = Path(output_path)
    
    # If it's a file path, get the directory
    if path.suffix:
        directory = path.parent
    else:
        directory = path
        
    directory.mkdir(parents=True, exist_ok=True)
    return directory

def get_consistent_color_for_drug(drug_name: str, drug_list: List[str]) -> str:
    """
    Generate a consistent color for a drug based on its position in the sorted drug list.
    This ensures the same drug always gets the same color across all plots.
    """
    try:
        drug_index = drug_list.index(drug_name)
        # Use matplotlib's tab20 colormap for up to 20 drugs, then extend
        if drug_index < 20:
            colors = plt.cm.tab20(np.linspace(0, 1, 20))
            return colors[drug_index]
        elif drug_index < 40:
            colors = plt.cm.tab20b(np.linspace(0, 1, 20))
            return colors[drug_index - 20]
        else:
            colors = plt.cm.tab20c(np.linspace(0, 1, 20))
            return colors[(drug_index - 40) % 20]
    except ValueError:
        # Fallback to default color if drug not in list
        return 'blue'

def get_consistent_color_for_bacteria(bacteria_name: str, bacteria_list: List[str]) -> str:
    """
    Generate a consistent color for bacteria based on its position in the sorted list.
    """
    try:
        bacteria_index = bacteria_list.index(bacteria_name)
        # Use different colormaps to distinguish from drugs
        if bacteria_index < 20:
            colors = plt.cm.Set3(np.linspace(0, 1, 12))  # Set3 has 12 colors
            return colors[bacteria_index % 12]
        else:
            colors = plt.cm.Paired(np.linspace(0, 1, 12))  # Paired has 12 colors  
            return colors[bacteria_index % 12]
    except ValueError:
        return 'red'  # Fallback

def format_large_numbers(number: float, precision: int = 2) -> str:
    """
    Format large numbers in a readable way (e.g., 1.5M, 2.3K).
    
    Args:
        number: Number to format
        precision: Decimal places to show
        
    Returns:
        Formatted string
    """
    if abs(number) >= 1e6:
        return f"{number/1e6:.{precision}f}M"
    elif abs(number) >= 1e3:
        return f"{number/1e3:.{precision}f}K"
    else:
        return f"{number:.{precision}f}"

def calculate_percentage_change(old_value: float, new_value: float) -> float:
    """
    Calculate percentage change between two values.
    
    Args:
        old_value: Original value
        new_value: New value
        
    Returns:
        Percentage change (positive for increase, negative for decrease)
    """
    if old_value == 0:
        return float('inf') if new_value > 0 else 0
    return ((new_value - old_value) / old_value) * 100

def add_statistics_text_box(ax, data: pd.Series, 
                           box_position: tuple = (0.02, 0.98),
                           box_props: dict = None) -> None:
    """
    Add a text box with basic statistics to a plot.
    
    Args:
        ax: Matplotlib axis object
        data: Data series for statistics
        box_position: (x, y) position for text box in axis coordinates
        box_props: Properties for the text box styling
    """
    if box_props is None:
        box_props = dict(boxstyle='round', facecolor='lightblue', alpha=0.7)
    
    # Calculate statistics
    valid_data = data.dropna()
    if len(valid_data) == 0:
        return
        
    mean_val = valid_data.mean()
    max_val = valid_data.max()
    min_val = valid_data.min()
    
    # Format text
    if mean_val < 0.01:
        textstr = f'Mean: {mean_val:.2e}\nMax: {max_val:.2e}\nMin: {min_val:.2e}'
    else:
        textstr = f'Mean: {mean_val:.3f}\nMax: {max_val:.3f}\nMin: {min_val:.3f}'
    
    # Add text box
    ax.text(box_position[0], box_position[1], textstr, 
           transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=box_props)

def setup_plot_style(style: str = 'seaborn-v0_8') -> None:
    """Configure matplotlib plot style."""
    try:
        plt.style.use(style)
    except Exception as e:
        logging.warning(f"Could not set plot style '{style}': {e}")
        # Fallback to default style
        plt.style.use('default')

def save_and_show_plot(filename: Union[str, Path], 
                      title: str = None,
                      dpi: int = 300,
                      bbox_inches: str = 'tight',
                      show_plot: bool = False) -> None:
    """
    Standardized plot saving and optional display.
    
    Args:
        filename: Output filename
        title: Plot title for logging
        dpi: Resolution for saved plot
        bbox_inches: Bounding box setting
        show_plot: Whether to display plot interactively
    """
    try:
        # Ensure output directory exists
        output_path = Path(filename)
        ensure_output_directory(output_path.parent)
        
        plt.tight_layout()
        plt.savefig(filename, dpi=dpi, bbox_inches=bbox_inches)
        
        if show_plot:
            plt.show()
        else:
            plt.close()  # Close to free memory
            
        print(f"[OK] {title or 'Plot'} saved as '{filename}'")
        
    except Exception as e:
        logging.error(f"Error saving plot to {filename}: {e}")
        plt.close()  # Ensure cleanup even on error
        raise

# Data extraction utilities
def extract_bacteria_list_from_csv(df: pd.DataFrame) -> List[str]:
    """
    Dynamically extract the list of bacteria from CSV column headers.
    This replaces hardcoded bacteria lists and automatically adapts to any BACTERIA_LIST configuration.
    """
    bacteria_list = []
    for col in df.columns:
        if col.endswith('_currently_infected'):
            bacteria_name = col.replace('_currently_infected', '')
            bacteria_list.append(bacteria_name)
    
    bacteria_list.sort()  # For consistent ordering
    print(f"Detected {len(bacteria_list)} bacteria from CSV headers:")
    for i, bacteria in enumerate(bacteria_list, 1):
        print(f"   {i}. {bacteria}")
    
    return bacteria_list

def extract_drug_list_from_csv(df: pd.DataFrame) -> List[str]:
    """
    Dynamically extract the list of drugs from CSV column headers.
    """
    drugs = []
    for col in df.columns:
        if col.endswith('_currently_on_drug'):
            drug_name = col.replace('_currently_on_drug', '')
            drugs.append(drug_name)
    
    drugs.sort()  # For consistent ordering
    print(f"Detected {len(drugs)} drugs from CSV headers")
    return drugs

def extract_resistance_mechanisms_from_csv(df: pd.DataFrame) -> List[str]:
    """
    Dynamically extract resistance mechanisms from CSV column headers.
    """
    mechanisms = set()
    for col in df.columns:
        if '_infected_with_' in col:
            # Extract mechanism name from column like "escherichia_coli_infected_with_esbl"
            parts = col.split('_infected_with_')
            if len(parts) == 2:
                mechanism = parts[1]
                mechanisms.add(mechanism)
    
    mechanism_list = sorted(list(mechanisms))
    print(f"Detected {len(mechanism_list)} resistance mechanisms: {', '.join(mechanism_list)}")
    return mechanism_list