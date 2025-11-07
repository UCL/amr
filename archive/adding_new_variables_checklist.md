# Adding New Variables to AMR Simulation: Complete Implementation Checklist

This document provides a comprehensive checklist for adding new variables to the AMR simulation system, ensuring proper data flow from individual tracking through CSV export.

## 📋 **Pre-Implementation Planning**

### 1. **Variable Definition & Scope**
- [ ] **Variable name and purpose**: Clearly define what the variable tracks
- [ ] **Data type**: Determine appropriate type (bool, usize, f64, Vec<T>, etc.)
- [ ] **Granularity level**: Per-individual, per-bacteria, per-drug, per-region, per-timestep, or combinations
- [ ] **Temporal behavior**: Does it reset each timestep, accumulate over time, or persist?
- [ ] **Default value**: What should the initial/reset value be?
- [ ] **Expected range**: Min/max values for validation

### 2. **Data Collection Requirements**
- [ ] **When to collect**: Which simulation phase (rules application, post-processing, etc.)
- [ ] **Collection frequency**: Every timestep, periodic, or event-triggered
- [ ] **Thread safety**: Will data be collected in parallel operations?
- [ ] **Performance impact**: Consider computational cost of collection
- [ ] **Memory requirements**: Estimate storage needs for full simulation

### 3. **CSV Output Specifications**
- [ ] **Column naming convention**: Follow existing patterns (use underscores, descriptive names)
- [ ] **Column positioning**: Early for key metrics, later for detailed breakdowns
- [ ] **Multiple columns needed**: Per-bacteria, per-drug, per-region expansions
- [ ] **Data format**: Integer counts, percentages, rates, etc.

## 🔧 **Implementation Steps**

### Step 1: Individual-Level Data Storage
**File: `src/simulation/population.rs`**

- [ ] **Add field to Individual struct**
  - Choose appropriate data type and size
  - Add comprehensive documentation comment
  - Consider Vec sizes (num_bacteria, num_drugs, combinations)

- [ ] **Initialize field in Individual::new()**
  - Set appropriate default values
  - Ensure Vec sizing matches constants (BACTERIA_LIST.len(), DRUG_SHORT_NAMES.len())
  - Consider initialization performance for large populations

- [ ] **Add field to individual logging CSV** (if `log_individuals = true`)
  - Update CSV header generation in simulation.rs
  - Update individual data writing loop
  - Consider flattening multi-dimensional data

### Step 2: Parallel Data Collection Infrastructure
**File: `src/simulation/simulation.rs`**

- [ ] **Add field to LocalTotals struct**
  - Choose collection data type (Vec<usize>, HashMap, etc.)
  - Consider thread-local capacity pre-allocation
  - Match aggregation needs (sums, counts, lists)

- [ ] **Initialize field in LocalTotals::new()**
  - Pre-allocate with appropriate capacity
  - Set correct Vec sizes for multi-dimensional data
  - Consider per-thread memory usage

- [ ] **Add aggregation in LocalTotals::merge()**
  - Implement proper combining logic (addition, concatenation, etc.)
  - Ensure thread-safe operations
  - Handle different data types appropriately

- [ ] **Add data collection in fold operation**
  - Locate appropriate collection point in simulation loop
  - Implement efficient collection logic
  - Consider early exit conditions for performance

- [ ] **Add field extraction in destructuring**
  - Add to LocalTotals destructuring pattern
  - Ensure field name matches struct definition

### Step 3: Summary Data Structure
**File: `src/simulation/simulation.rs`**

- [ ] **Add field to TimeStepSummary struct**
  - Use descriptive name and appropriate type
  - Add documentation comment
  - Position logically with related fields

- [ ] **Add field to summary creation**
  - Include in TimeStepSummary initialization
  - Ensure field name matches struct definition
  - Apply any necessary data transformations

### Step 4: CSV Export Implementation
**File: `src/simulation/simulation.rs` - `export_summary_to_csv()` function**

- [ ] **Add header generation**
  - Position columns appropriately (early vs. late in CSV)
  - Follow naming conventions (bacteria_name_variable_name)
  - Generate headers for all dimensions (per-bacteria, per-drug, etc.)
  - Use consistent formatting (replace spaces with underscores)

- [ ] **Add data writing**
  - Write data in same order as headers
  - Handle multi-dimensional data with nested loops
  - Use appropriate formatting for data type
  - Ensure consistent row structure

- [ ] **Consider pre-allocation**
  - Update row String capacity if adding many columns
  - Consider header String capacity for long headers

### Step 5: Data Source Implementation
**Files: Various depending on variable purpose**

- [ ] **Implement data setting logic**
  - **Rules engine** (`src/rules/mod.rs`): For simulation rule outcomes
  - **Population dynamics**: For demographic changes
  - **Drug interactions**: For pharmacological effects
  - **Resistance mechanisms**: For genetic/evolutionary changes

- [ ] **Add reset mechanism** (if needed)
  - Reset per-timestep variables before new timestep
  - Add to existing reset loops or create new ones
  - Ensure proper timing (before vs. after data collection)

### Step 6: Data Validation & Testing
- [ ] **Compilation testing**
  - Verify all struct field additions compile
  - Check all initialization points
  - Validate data type consistency

- [ ] **Simulation testing**
  - Run short test simulation
  - Verify CSV generation
  - Check column headers are correct
  - Validate data ranges are reasonable

- [ ] **Data verification**
  - Check for expected patterns (zeros, non-zeros, trends)
  - Verify column counts match expectations
  - Test edge cases (empty data, maximum values)

## ⚠️ **Common Pitfalls & Considerations**

### Memory & Performance
- [ ] **Vec sizing**: Always initialize Vecs with correct capacity (num_bacteria, num_drugs)
- [ ] **Thread safety**: Ensure parallel collection doesn't create race conditions
- [ ] **Memory usage**: Large Vec<Vec<T>> structures can consume significant memory
- [ ] **Collection frequency**: Per-individual collection in large populations can be expensive

### Data Consistency
- [ ] **Field order**: Maintain consistent field ordering between struct definition and usage
- [ ] **Naming consistency**: Use same naming pattern throughout (snake_case, descriptive)
- [ ] **Type consistency**: Ensure data types match between collection and storage
- [ ] **Reset timing**: Reset flags/counters at appropriate simulation phases

### CSV Export Alignment
- [ ] **Header-data ordering**: Headers and data must be written in identical order
- [ ] **Multi-dimensional expansion**: Nested loops must match header generation logic
- [ ] **Column count**: Total columns must match between header and data rows
- [ ] **Missing data handling**: Handle cases where data might be missing or invalid

### Debugging & Maintenance
- [ ] **Documentation**: Add clear comments explaining variable purpose and behavior
- [ ] **Logging**: Consider adding debug logging for new variable collection
- [ ] **Error handling**: Add appropriate error checking for data validation
- [ ] **Future extensibility**: Design with potential future enhancements in mind

## 📝 **Variable Request Template**

When requesting a new variable, provide:

```markdown
## New Variable Request

**Variable Name**: `variable_name`
**Purpose**: What does this variable track/measure?
**Data Type**: bool/usize/f64/Vec<T>
**Scope**: Per-individual/per-bacteria/per-drug/per-region/combination
**Temporal Behavior**: Resets each timestep/accumulates/persists
**Collection Trigger**: When should data be collected?
**CSV Output**: How many columns? Naming pattern?
**Expected Usage**: How will this data be analyzed?
**Performance Considerations**: Any special computational requirements?
```

## 🔍 **Verification Checklist**

After implementation:

- [ ] Code compiles without errors or warnings
- [ ] Simulation runs to completion
- [ ] CSV file generates successfully  
- [ ] Column headers are present and correctly named
- [ ] Data values are within expected ranges
- [ ] No missing or corrupted data in CSV
- [ ] Column count matches header count
- [ ] Data collection performs adequately
- [ ] Memory usage is reasonable
- [ ] Implementation is properly documented

## 📚 **Reference Examples**

**Simple per-bacteria boolean**: `infection_prevented_by_drug: Vec<bool>`
**Complex multi-dimensional**: `resistance_by_bacteria_drug: Vec<Vec<usize>>`
**Regional breakdown**: `deaths_by_region: Vec<usize>` 
**Event counting**: `newly_infected_count: usize`
**Rate calculation**: `background_all_cause_mortality_rate: f64`

---

*This checklist should be consulted for every new variable addition to ensure comprehensive and correct implementation across the entire simulation pipeline.*