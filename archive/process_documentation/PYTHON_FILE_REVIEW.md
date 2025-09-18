# Python File Review and Archival Plan

## 📁 **Current Python Files Analysis**

### **✅ KEEP (Essential for ongoing work)**

**1. `analyze_simulation.py`**
- **Purpose**: Main analysis and visualization script
- **Status**: ✅ **ESSENTIAL** - This is the primary user interface
- **Notes**: Updated with clean empirical data integration

**2. `empirical_enhancement.py`**
- **Purpose**: Single integrated empirical data generation module
- **Status**: ✅ **ESSENTIAL** - This is our consolidated enhancement solution
- **Notes**: Replaces all other enhancement scripts

**3. `individuals_log_view.py`**
- **Purpose**: Individual-level simulation data analysis
- **Status**: ✅ **KEEP** - Specialized analysis tool, not related to empirical enhancement
- **Notes**: Independent functionality for detailed simulation inspection

---

### **📦 ARCHIVE (Superseded by clean implementation)**

**Enhancement Development Files (Superseded by `empirical_enhancement.py`)**

**4. `targeted_empirical_enhancement.py`**
- **Purpose**: Phase 1 enhancement (WHO GLASS core patterns)
- **Status**: 📦 **ARCHIVE** - Superseded by integrated approach
- **Reason**: Functionality integrated into `empirical_enhancement.py`

**5. `phase2_resistance_enhancement.py`**
- **Purpose**: Phase 2 resistance enhancement (ECDC + CDC patterns)
- **Status**: 📦 **ARCHIVE** - Superseded by integrated approach
- **Reason**: Functionality integrated into `empirical_enhancement.py`

**6. `comprehensive_resistance_enhancement.py`**
- **Purpose**: Comprehensive resistance coverage expansion
- **Status**: 📦 **ARCHIVE** - Superseded by integrated approach
- **Reason**: Functionality integrated into `empirical_enhancement.py`

**7. `enhanced_empirical_data_collection.py`**
- **Purpose**: Data collection strategy exploration
- **Status**: 📦 **ARCHIVE** - Development/exploration file
- **Reason**: Concepts integrated into final implementation

**8. `enhanced_empirical_data_config.py`**
- **Purpose**: Configuration for enhanced data collection
- **Status**: 📦 **ARCHIVE** - Development/exploration file
- **Reason**: Configuration now integrated into main modules

**9. `enhanced_empirical_data_strategy.py`**
- **Purpose**: Strategy development for empirical enhancement
- **Status**: 📦 **ARCHIVE** - Development/exploration file
- **Reason**: Strategy implemented in final solution

**10. `generate_empirical_calibration.py`**
- **Purpose**: Original empirical calibration generation
- **Status**: 📦 **ARCHIVE** - Superseded by integrated approach
- **Reason**: Replaced by `empirical_enhancement.py`

**Documentation/Guide Files**

**11. `implementation_guide.py`**
- **Purpose**: Implementation documentation in code form
- **Status**: 📦 **ARCHIVE** - Documentation file
- **Reason**: Information now in markdown documentation

**12. `demo_real_data_integration.py`**
- **Purpose**: Demonstration of real data integration
- **Status**: 📦 **ARCHIVE** - Demo/example file
- **Reason**: Functionality now standard in main workflow

---

### **🤔 REVIEW (Uncertain)**

**13. `empirical_data_parsers.py`**
- **Purpose**: Parsing utilities for empirical data
- **Status**: 🤔 **REVIEW NEEDED** - May contain useful utilities
- **Action**: Check if used by current implementation

---

## 🗂️ **Archival Plan**

### **Step 1: Create Archive Structure**
```powershell
mkdir archive/development_scripts
mkdir archive/exploration_scripts
mkdir archive/documentation_scripts
```

### **Step 2: Archive Enhancement Development Files**
```powershell
# Phase-based enhancement scripts (superseded)
mv targeted_empirical_enhancement.py archive/development_scripts/
mv phase2_resistance_enhancement.py archive/development_scripts/
mv comprehensive_resistance_enhancement.py archive/development_scripts/

# Enhanced data exploration scripts (superseded)
mv enhanced_empirical_data_collection.py archive/exploration_scripts/
mv enhanced_empirical_data_config.py archive/exploration_scripts/
mv enhanced_empirical_data_strategy.py archive/exploration_scripts/

# Original generation script (superseded)
mv generate_empirical_calibration.py archive/development_scripts/
```

### **Step 3: Archive Documentation Files**
```powershell
# Documentation and demo scripts
mv implementation_guide.py archive/documentation_scripts/
mv demo_real_data_integration.py archive/documentation_scripts/
```

### **Step 4: Review Uncertain Files**
```powershell
# Check if empirical_data_parsers.py is used
grep -r "empirical_data_parsers" *.py
```

---

## 📋 **Final Clean State**

### **Working Directory (Post-Cleanup)**
```
📁 Active Python Files:
├── analyze_simulation.py          # Main analysis script ✅
├── empirical_enhancement.py       # Empirical data generation ✅  
├── individuals_log_view.py        # Individual-level analysis ✅
└── [empirical_data_parsers.py]    # TBD after review 🤔

📁 Archive:
├── development_scripts/           # Enhancement development history
├── exploration_scripts/           # Data strategy exploration
└── documentation_scripts/         # Documentation and demos
```

### **Benefits of Cleanup**
- **Clear purpose**: Each remaining file has a distinct role
- **No confusion**: No multiple scripts doing similar things
- **Maintainable**: Easy to understand what each file does
- **Historical record**: Development history preserved in archive
- **Clean workspace**: Focus on essential tools only

---

## ✅ **Recommendation**

**Proceed with archival** - We have a clear consolidation where:
- `empirical_enhancement.py` replaces 6 separate enhancement scripts
- `analyze_simulation.py` is the main user interface
- Development history is preserved in archive
- Clean, purposeful workspace remains

This follows software development best practices of consolidating functionality while preserving development history.