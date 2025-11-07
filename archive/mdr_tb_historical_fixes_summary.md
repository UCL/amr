# MDR TB Historical Modeling Fixes

## Issues Identified

### 1. **Historical Implausibility (1946 MDR TB Case)**
- **Problem**: MDR TB cases occurring in 1946, before effective TB drugs existed
- **Medical Issue**: MDR TB requires resistance to rifampicin + isoniazid, neither available until 1966+
- **Simulation Issue**: 4-year-old with bone/joint MDR TB cleared by sulfanilamide + penicillin (impossible)

### 2. **Treatment Efficacy Paradox**  
- **Problem**: Ineffective drugs (sulfanilamide 0.06 potency, penicillin 0.0001 potency) clearing MDR TB
- **Medical Reality**: Pre-antibiotic TB mortality was 50%+ in children
- **Simulation Reality**: Near-universal clearance with ineffective treatment

### 3. **Historical Drug Timeline Issues**
- **Streptomycin**: Available from 1944 (day 5,115) - first TB drug
- **Rifampicin**: Available from 1966 (day 13,140) - defines MDR TB
- **Problem**: Model assigns 90% rifampicin resistance before rifampicin exists

## Fixes Implemented

### 1. **Time-Dependent MDR TB Incidence**
```rust
// New parameters in config.rs:
map.insert("mdr_tb_pre_antibiotic_era_multiplier", 0.0001); // Pre-1944: 99.99% reduction
map.insert("mdr_tb_early_antibiotic_era_multiplier", 0.01); // 1944-1965: 99% reduction  
map.insert("mdr_tb_modern_era_multiplier", 1.0);           // 1966+: full rates
```

### 2. **Historical Acquisition Rate Modulation**
```rust
// In rules/mod.rs infection acquisition:
let simulation_year = 1930.0 + (time_step as f64 / 365.0);
if bacteria == "mdr mycobacterium tuberculosis" {
    let mdr_tb_multiplier = match simulation_year {
        y if y < 1944.0 => 0.0001,  // Pre-antibiotic era
        y if y < 1966.0 => 0.01,    // Early antibiotic era  
        _ => 1.0                    // Modern era
    };
    acquisition_probability *= mdr_tb_multiplier;
}
```

### 3. **Time-Dependent Treatment Effectiveness**
```rust
// Historical treatment effectiveness in TB synergy calculation:
if simulation_year < 1944.0 {
    background_effectiveness *= 0.01; // 99% reduction - no effective treatment
} else if simulation_year < 1966.0 {
    background_effectiveness *= 0.3;  // 70% reduction - limited monotherapy
}
```

### 4. **Era-Appropriate Rifampicin Resistance**
```rust
// Only apply guaranteed rifampicin resistance after rifampicin availability:
let guaranteed_rifampicin_resistance = if is_tb && simulation_year >= 1966.0 {
    get_global_param("mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance")
        .unwrap_or(0.90)
} else {
    0.0
};
```

## Expected Outcomes

### **Pre-Antibiotic Era (1930-1944)**
- **MDR TB Incidence**: ~0.01% of normal rates (virtually eliminated)
- **Treatment**: 99% reduced effectiveness (reflecting historical reality)
- **Mortality**: Should approach historical levels (50%+ for children)

### **Early Antibiotic Era (1944-1966)** 
- **MDR TB Incidence**: 1% of modern rates (monotherapy resistance emergence)
- **Treatment**: 70% reduced effectiveness (limited to streptomycin monotherapy)
- **Rifampicin Resistance**: 0% (drug not yet available)

### **Modern Era (1966+)**
- **MDR TB Incidence**: Full model rates
- **Treatment**: Full effectiveness with multi-drug regimens
- **Rifampicin Resistance**: 90% guaranteed (defines MDR TB)

## Clinical Accuracy Improvements

1. **Historical Authenticity**: No MDR TB before drugs to resist
2. **Treatment Realism**: Appropriate mortality in pre-effective treatment eras  
3. **Drug Timeline Accuracy**: Resistance only to available drugs
4. **Epidemiological Accuracy**: MDR TB emergence follows antibiotic introduction patterns

## Testing Notes

- 1946 case should now be extremely rare (0.01% probability)
- Pre-1944 MDR TB cases should be virtually eliminated
- Treatment failures should increase dramatically in early eras
- Post-1966 behavior should remain unchanged

## Files Modified

1. **src/config.rs**: Added time-dependent incidence parameters
2. **src/rules/mod.rs**: Implemented historical logic in:
   - Infection acquisition probability
   - Treatment effectiveness calculation  
   - Rifampicin resistance assignment

## Impact on Simulation

- **Backward Compatibility**: Post-1966 behavior unchanged
- **Historical Accuracy**: Pre-antibiotic era now medically realistic
- **Performance**: Minimal impact (simple multiplier calculations)
- **Validation**: Eliminates impossible 1946 MDR TB scenarios