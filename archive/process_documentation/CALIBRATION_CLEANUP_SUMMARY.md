# Calibration CSV File Cleanup Summary

## Final Clean Calibration Structure (September 18, 2025)

### ✅ **Essential Calibration Files Remaining (4 files)**
These are the active files used by the current clean implementation:
- **`calibration_drug_usage_empirical.csv`** (4.8 MB) - Drug usage patterns
- **`calibration_resistance_empirical.csv`** (33.9 MB) - Resistance patterns  
- **`calibration_infection_incidence_empirical.csv`** (5.0 MB) - Infection incidence patterns
- **`calibration_deaths_empirical.csv`** (4.3 MB) - Mortality patterns

### 🗑️ **Deleted Development Files (10 files)**
These superseded development artifacts were deleted (not archived):

#### Enhanced Versions (3 files deleted)
- `calibration_deaths_empirical_ENHANCED.csv`
- `calibration_infection_incidence_empirical_ENHANCED.csv`  
- `calibration_resistance_empirical_ENHANCED.csv`

#### Targeted Enhanced Versions (3 files deleted)
- `calibration_deaths_empirical_ENHANCED_TARGETED.csv`
- `calibration_infection_incidence_empirical_ENHANCED_TARGETED.csv`
- `calibration_resistance_empirical_ENHANCED_TARGETED.csv`

#### Comprehensive Versions (2 files deleted)
- `calibration_drug_usage_empirical_COMPREHENSIVE.csv`
- `calibration_resistance_empirical_COMPREHENSIVE.csv`

#### Other Development Versions (2 files deleted)
- `calibration_resistance_empirical_PHASE2_ENHANCED.csv`
- `calibration_drug_usage_empirical_BACKUP.csv`
- `calibration_drug_usage_empirical_ENHANCED.csv`

### 📄 **Other CSV Files Kept (4 files)**
Non-calibration simulation output and logging files:
- `individuals_log.csv` - Individual patient simulation log
- `simulation_run_log.csv` - Simulation run tracking
- `simulation_summary.csv` - Simulation summary output
- `summary_statistics.csv` - Statistical analysis output

### 📦 **Original Data Preserved**
Original synthetic calibration files remain archived in `archive/original_synthetic/`:
- `calibration_deaths_empirical.csv` (original synthetic)
- `calibration_infection_incidence_empirical.csv` (original synthetic)  
- `calibration_resistance_empirical.csv` (original synthetic)

## Verification Results
✅ All 4 essential calibration files verified functional (13 columns each)  
✅ No development artifacts remain in workspace  
✅ Clean implementation uses standard file names without suffixes  
✅ Total CSV files reduced from 19 → 8 files  

## Achievement
- **Eliminated clutter**: Removed 10 superseded development files
- **Clean workspace**: Only essential calibration files + simulation outputs
- **Standard naming**: No confusing suffixes (_ENHANCED, _TARGETED, etc.)
- **Preserved history**: Original synthetic data safely archived
- **Simple deletion**: No archival needed for development artifacts

## Usage
- **Current workflow**: `analyze_simulation.py` automatically uses the 4 clean calibration files
- **Data generation**: `empirical_enhancement.py` generates/updates the 4 clean calibration files
- **Archive reference**: Original synthetic versions available in `archive/original_synthetic/`