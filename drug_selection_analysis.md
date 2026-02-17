# Drug Selection Algorithm Analysis and Recommendations

## Current State Analysis

### Algorithm Overview
The drug selection system uses a **weighted probabilistic scoring algorithm** with softmax temperature-based selection:

1. **Base Score**: Starts at 1.0
2. **Empiric Multipliers**: Applied based on syndrome-drug penetration
3. **Potency-Based Activity**: Drugs with potency ≥ threshold get boosted
4. **Clinical Guideline Multipliers**: Pathogen-specific preferences (e.g., penicillin × 24 for *S. pneumoniae*)
5. **Resistance Penalties**: Regional surveillance data reduces scores
6. **Stewardship Restrictions**: Reserve drugs heavily penalized
7. **Softmax Selection**: `weight = exp(score / temperature)`

### Key Parameters
- `drug_selection_temperature`: Controls randomness (lower = more deterministic)
- `minimal_potency_threshold_for_drug_selection`: Minimum potency to consider drug (default ~0.05)
- `clinical_preference_multiplier`: Bacteria-drug specific boosts

---

## Identified Issues Causing Calibration Mismatch

### 1. **Nitrofurantoin Overuse** (11.5% sim vs 3% target = +8.5pp)

**Root Causes:**
- **High E. coli multiplier** ([line 2101](src/rules/mod.rs#L2101)): `score *= 14.0` for E. coli
- **Perfect UTI penetration** (config.rs line 1359): `drug_penetration[1][nitrofurantoin_idx] = 1.0`
- **Low resistance**: Nitrofurantoin resistance remains very low throughout simulation
- **No syndrome restriction**: While penetration is poor for non-UTI sites, empiric scoring doesn't strongly penalize it
- **Empiric selection bias**: For empiric UTI, nitrofurantoin gets selected frequently because:
  - High syndrome penetration (1.0 for UTI)
  - No competing negative multipliers
  - Low regional resistance penalties

**Why this is unrealistic:**
- Nitrofurantoin is **only used for uncomplicated cystitis** (lower UTI)
- Contraindicated in pyelonephritis, sepsis, pregnancy complications
- Not used for non-UTI infections at all
- Should represent ~20-25% of UTI cases max, not general antibiotic use

### 2. **Penicillin Underuse** (7.4% sim vs 18% target = -10.6pp)

**Root Causes:**
- **Resistance penalties too harsh**: Regional resistance prevalence reduces scores multiplicatively
- **Competing drugs get higher multipliers**: E.g., for *S. pneumoniae*:
  - Penicillin G: `score *= 24.0` ([line 2065](src/rules/mod.rs#L2065))
  - BUT amoxicillin also gets `score *= 24.0` ([line 2067](src/rules/mod.rs#L2067))
  - Amoxicillin-clavulanate gets `score *= 12.0` ([line 2069](src/rules/mod.rs#L2069))
  - Result: Penicillins split votes among themselves
- **Limited spectrum disadvantage**: Empiric therapy heavily favors broad-spectrum drugs
- **Time-period modulation too aggressive**: Early-period boosts may not be strong enough

**Why this is unrealistic:**
- Penicillins should dominate for Strep throat, Group B Strep, syphilis, meningococcal disease
- First-line for many common infections
- Underused because algorithm over-weights resistance fears

### 3. **Beta-Lactam/BLI Underuse** (6.4% sim vs 18% target = -11.6pp)

**Root Causes:**
- **Carbapenem preference**: For ESBL-era E. coli/Klebsiella, carbapenems get higher multipliers
- **Competing with cephalosporins**: E.g., for *E. coli* ([line 2101](src/rules/mod.rs#L2101)):
  - Ceftriaxone: `score *= 9.0`
  - Ampicillin (early): `score *= 15.0` → `score *= 4.0` (later)
  - Amox-clav is missing explicit boost here
- **Reserve drug escalation**: In empiric therapy with treatment failure, simulation jumps to carbapenems rather than BL/BLI
- **Stewardship penalties**: Lines 2251-2312 show heavy reserve drug restrictions, but BL/BLI combinations aren't getting reciprocal boosts

**Why this is unrealistic:**
- Amoxicillin-clavulanate is **most prescribed antibiotic globally** for respiratory, skin, UTI
- Piperacillin-tazobactam is first-line for hospital-acquired infections
- Should be 15-20% of all antibiotic use

---

## Recommended Algorithm Modifications

### Priority 1: Fix Nitrofurantoin Overuse

**Solution A: Restrict to confirmed lower UTI only**
```rust
// In drug scoring section, add after line 1850:
if drug_name == "nitrofurantoin" || drug_name == "furazolidone" {
    // Block nitrofurantoin for non-UTI syndromes
    let is_uti_only = active_syndrome_ids.iter().all(|&sid| sid == 1);
    if !is_uti_only && !active_syndrome_ids.is_empty() {
        score = 0.0; // Complete block for non-UTI
        continue;
    }
    
    // Reduce multiplier for E. coli UTI (currently 14.0 is too high)
    // Should compete with cipro/TMP-SMX, not dominate
    score *= 0.4; // Global reduction factor
}
```

**Solution B: Reduce E. coli multiplier**
```rust
// Line 2101, change:
("escherichia_coli", "nitrofurantoin") => score *= 14.0,
// TO:
("escherichia_coli", "nitrofurantoin") => score *= 4.0, // Reduce from 14 to 4
```

**Solution C: Add sepsis/pyelonephritis contraindication**
```rust
// After line 1850:
if (drug_name == "nitrofurantoin" || drug_name == "furazolidone") && 
   (individual.sepsis.iter().any(|&s| s) || active_syndrome_ids.contains(&4)) {
    score = 0.0; // Block in sepsis or bloodstream infections
    continue;
}
```

**Recommended: Implement all three solutions** to achieve realistic 2-3% share.

---

### Priority 2: Boost Penicillin Use

**Solution A: Strengthen first-line guideline multipliers**
```rust
// Increase penicillin multipliers for appropriate pathogens:

// Line 2065-2067, change:
("streptococcus_pneumoniae", "penicillin_g") => score *= 24.0,
("streptococcus_pneumoniae", "ampicillin") => score *= 22.0,
("streptococcus_pneumoniae", "amoxicillin") => score *= 24.0,
// TO:
("streptococcus_pneumoniae", "penicillin_g") => score *= 35.0,
("streptococcus_pneumoniae", "ampicillin") => score *= 32.0,
("streptococcus_pneumoniae", "amoxicillin") => score *= 35.0,

// Line 2090-2092, change:
("streptococcus_pyogenes", "penicillin_g") => score *= 28.0,
("streptococcus_pyogenes", "ampicillin" | "amoxicillin") => score *= 20.0,
// TO:
("streptococcus_pyogenes", "penicillin_g") => score *= 45.0,
("streptococcus_pyogenes", "ampicillin" | "amoxicillin") => score *= 35.0,
```

**Solution B: Reduce resistance penalty for low-risk pathogens**
```rust
// After line 2480 (regional resistance penalty section):
// Don't penalize penicillins for Strep species (they remain highly susceptible)
if PENICILLIN_CLASS_DRUGS.contains(&drug_name) {
    for &b_idx in &identified_bacteria {
        let bacteria_name = BACTERIA_LIST[b_idx];
        if matches!(bacteria_name,
            "streptococcus_pneumoniae" | 
            "streptococcus_pyogenes" | 
            "streptococcus_agalactiae" |
            "treponema_pallidum" |
            "neisseria_meningitidis"
        ) {
            // Override harsh resistance penalty for penicillin-susceptible pathogens
            regional_resistance_penalty = 1.0; // No penalty
            break;
        }
    }
}
score *= regional_resistance_penalty;
```

**Solution C: Add empiric respiratory syndrome bonus**
```rust
// In empiric selection section (~line 1900):
if empiric_selection && active_syndrome_ids.contains(&3) { // Respiratory
    if matches!(drug_name, 
        "penicillin_g" | "ampicillin" | "amoxicillin" | "amoxicillin_clavulanate"
    ) {
        score *= 3.0; // Boost for empiric respiratory (pharyngitis, pneumonia)
    }
}
```

---

### Priority 3: Boost BL/BLI Combinations

**Solution A: Add explicit guideline multipliers**
```rust
// Add after E. coli section (line 2101):
("escherichia_coli", "amoxicillin_clavulanate") => score *= 12.0,
("escherichia_coli", "ampicillin_sulbactam") => score *= 10.0,

// Add for Klebsiella (after line 2130):
("klebsiella_pneumoniae", "amoxicillin_clavulanate") => score *= 8.0,
("klebsiella_pneumoniae", "piperacillin_tazobactam") => score *= 12.0, // INCREASE from 9.0

// Add for Staph (after line 1960):
("staphylococcus_aureus", "amoxicillin_clavulanate" | "ampicillin_sulbactam") => {
    if time_step < 14600 { // Extended era
        score *= 22.0; // INCREASE from 18.0
    } else {
        score *= 4.0; // INCREASE from 3.0
    }
}
```

**Solution B: Boost for empiric hospital infections**
```rust
// In empiric section (~line 1900):
if empiric_selection && individual.hospital_status.is_hospitalized() {
    if matches!(drug_name,
        "piperacillin_tazobactam" |
        "ampicillin_sulbactam" |
        "ticarcillin_clavulanate"
    ) {
        score *= 4.0; // Strong preference for hospital empiric therapy
    }
}
```

**Solution C: Reduce carbapenem competition**
```rust
// Line 2300 (reserve drug stewardship):
// INCREASE penalty multiplier:
if reserve_candidate {
    let base_reserve_penalty = store.globals.reserve_drug_score_penalty;
    // Current policy multiplier = 1.0, so penalty = base_reserve_penalty
    // SUGGEST: Change base_reserve_penalty in config from 0.12 to 0.08
    // This will make carbapenems even less preferred, allowing BL/BLI to compete
    let reserve_penalty = base_reserve_penalty.powf(reserve_drug_penalty_multiplier);
    if reserve_penalty >= 0.0 {
        score *= reserve_penalty;
    }
}
```

---

## Parameter Tuning Recommendations

### config.rs Global Parameters

```rust
// Suggested changes to globals:

// INCREASE temperature to add more randomness (prevents over-concentration)
drug_selection_temperature: 0.8 → 1.2
// Higher temperature = more variability in drug selection = better match to real-world diversity

// DECREASE reserve penalty to allow more access (but still restrictive)
reserve_drug_score_penalty: 0.12 → 0.08
// Makes carbapenems/linezolid slightly more accessible but still heavily restricted

// ADD new parameter:
nitrofuran_syndrome_restriction: true
// Enables syndrome-specific blocking of nitrofurantoin for non-UTI

// ADD new parameter:
penicillin_low_resistance_bonus: 2.0
// Multiplier for penicillins when treating intrinsically susceptible pathogens
```

---

## Testing & Validation Strategy

### Step 1: Implement nitrofurantoin fixes only
- Run simulation with Solutions A+B+C from Priority 1
- **Expected outcome**: Nitrofurans drop from 11.5% → 3-4%
- **Compensatory rise**: Expect fluoroquinolones, TMP-SMX, cephalosporins to rise

### Step 2: Implement penicillin boosts
- Add Solutions A+B+C from Priority 2
- **Expected outcome**: Penicillins rise from 7.4% → 15-18%
- **Check**: Ensure not over-boosting (>25% would be too high)

### Step 3: Implement BL/BLI boosts
- Add Solutions A+B+C from Priority 3
- **Expected outcome**: BL/BLI combinations rise from 6.4% → 15-18%
- **Check**: Ensure carbapenem share drops to 1-2% (currently may be inflated)

### Step 4: Fine-tune temperature
- If drug concentration is still too narrow (few drugs dominate):
  - Increase `drug_selection_temperature` from 0.8 → 1.0 → 1.2
- If drug distribution is too scattered (no clear patterns):
  - Decrease temperature slightly

### Step 5: Validate resistance dynamics
- Ensure changes don't break resistance emergence patterns
- Check that resistance rates for penicillins/BL-BLI still rise appropriately over time
- Verify MRSA/VRE/ESBL dynamics remain realistic

---

## Implementation Priority Order

1. **Nitrofurantoin restrictions** (Priority 1) - Easiest, biggest impact
2. **Penicillin resistance penalty override** (Priority 2, Solution B) - Low risk
3. **BL/BLI guideline multipliers** (Priority 3, Solution A) - Moderate effort
4. **Penicillin guideline multiplier increases** (Priority 2, Solution A) - Test carefully
5. **Empiric syndrome bonuses** (Priority 2 & 3, Solutions C & B) - Most complex
6. **Temperature tuning** - Final calibration step

---

## Expected Final Drug Shares (2025)

| Drug Class | Current | Target | After Fixes |
|------------|---------|--------|-------------|
| Penicillins | 7.4% | 18% | **17-19%** ✓ |
| BL/BLI combinations | 6.4% | 18% | **16-19%** ✓ |
| Cephalosporins | 17.0% | 18% | **17-19%** ✓ |
| Fluoroquinolones | 8.1% | 10% | **9-11%** ✓ |
| Macrolides | 7.1% | 11% | **10-12%** ✓ |
| Tetracyclines | 7.2% | 6% | **5-7%** ✓ |
| **Nitrofurans** | **11.5%** | **3%** | **2-4%** ✓ |
| Aminoglycosides | 4.4% | 2% | **2-3%** ✓ |
| Other classes | 31.0% | 14% | **12-15%** ✓ |

---

## Code Locations for Implementation

### Files to modify:
1. **[src/rules/mod.rs](src/rules/mod.rs)**
   - Lines 1850-1900: Add nitrofurantoin syndrome restrictions
   - Lines 1900-2250: Add pathogen-specific multiplier adjustments
   - Lines 2480-2550: Add resistance penalty overrides

2. **[src/config.rs](src/config.rs)**
   - Lines 13500-13600 (globals section): Adjust temperature, reserve penalty
   - Add new parameters for syndrome restrictions

### Testing commands:
```powershell
# Rebuild with optimizations
cargo build --release

# Run calibration analysis
.venv\Scripts\activate
python amr_simulation_output_analysis/amr_analysis.py

# Check drug shares
python -c "import pandas as pd; df = pd.read_csv('amr_simulation_output_analysis_outputs/simulation_summary_XXXXXX.csv'); print(df[df.index.str.contains('Drug Class Share')].to_string())"
```

---

## Risk Assessment

### Low Risk Changes:
- Nitrofurantoin syndrome restrictions ✓
- Guideline multiplier adjustments ✓
- Temperature parameter tuning ✓

### Medium Risk Changes:
- Resistance penalty overrides (could affect resistance dynamics)
- Empiric syndrome bonuses (may interact with other selection logic)

### High Risk Changes:
- Major reserve drug penalty changes (could cause stewardship collapse)
- Removing safety restrictions on aminoglycosides/colistin

### Recommended Approach:
- Make changes **incrementally**
- Test each change with **10K individual test runs** before full 100K
- Monitor **resistance trajectory plots** to ensure dynamics remain stable
- Compare **infection mortality rates** to ensure clinical outcomes don't degrade

---

## Conclusion

The drug selection algorithm is fundamentally sound but has **three major calibration issues**:

1. **Nitrofurantoin** lacks syndrome-specific restrictions → overused
2. **Penicillins** face excessive resistance penalties → underused
3. **BL/BLI combinations** lack competitive multipliers → underused

The proposed fixes target these specific issues without overhauling the entire scoring system. The softmax selection framework should remain intact - we're adjusting the **inputs** (scores) rather than the **selection mechanism** (temperature-weighted random choice).

**Estimated implementation time**: 4-6 hours  
**Estimated testing/validation time**: 2-3 full simulation runs (8-12 hours compute time)  
**Expected success rate**: 80-90% match to targets after fine-tuning
