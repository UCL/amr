# Resistance Mechanism Consistency Audit

## Date: February 13, 2026

## Issue Summary

The resistance mechanism emergence multiplier blocks in `config.rs` are:
1. **Incomplete** - Missing 10 mechanisms that are defined in `population.rs`
2. **Inconsistent** - Mechanisms listed in different orders for each bacteria

## Complete List of Resistance Mechanisms (from population.rs)

Based on `src/simulation/population.rs` lines 72-145, all 25 mechanisms are:

### 1. ESBL Enzymes (Extended-Spectrum Beta-Lactamases)
- `enzyme_esbl_ctx_m` ✅ Present in config
- `enzyme_esbl_tem` ❌ **MISSING from config**
- `enzyme_esbl_shv` ❌ **MISSING from config**

### 2. Carbapenemase Enzymes
- `enzyme_kpc` ✅ Present in config
- `enzyme_ndm_vim` ✅ Present in config
- `enzyme_oxa_48` ✅ Present in config

### 3. AmpC Beta-Lactamases
- `enzyme_ampc_cmy` ✅ Present in config
- `enzyme_ampc_dha` ❌ **MISSING from config**

### 4. Target Site Modifications
- `target_site_pbp2a_meca` ✅ Present in config (MRSA)
- `target_site_van_a` ✅ Present in config (VRE)
- `target_site_van_b` ❌ **MISSING from config** (VRE variant)
- `target_site_erm_b` ✅ Present in config (Macrolide resistance)
- `target_site_cfr` ❌ **MISSING from config** (Linezolid/Chloramphenicol)

### 5. Fluoroquinolone Resistance
- `mutation_gyra_primary` ✅ Present in config
- `mutation_gyra_parc_secondary` ❌ **MISSING from config** ⚠️ USER SPOTTED
- `protection_qnr` ✅ Present in config

### 6. Aminoglycoside Resistance
- `enzyme_16s_rrmt` ✅ Present in config

### 7. Other Enzymatic Resistance
- `enzyme_cat` ✅ Present in config (Chloramphenicol)

### 8. Efflux Pumps
- `efflux_acrab_tolc` ✅ Present in config
- `efflux_mexxy_oprm` ❌ **MISSING from config** (Pseudomonas-specific)
- `global_efflux_pump` ❌ **MISSING from config**

### 9. Porin Loss
- `porin_loss_ompk35_36` ✅ Present in config
- `porin_loss_oprd` ❌ **MISSING from config** (Pseudomonas-specific)
- `global_porin_loss` ❌ **MISSING from config**

### 10. Surface Modifications
- `modification_mcr_1` ✅ Present in config (Colistin)

---

## Current State: 16/25 mechanisms present (64%)

**Missing 10 mechanisms (40%):**
1. enzyme_esbl_tem
2. enzyme_esbl_shv
3. enzyme_ampc_dha
4. target_site_van_b
5. target_site_cfr
6. mutation_gyra_parc_secondary ⭐ User-identified
7. efflux_mexxy_oprm
8. porin_loss_oprd
9. global_efflux_pump
10. global_porin_loss

---

## Recommended Standard Order (Biological Grouping)

All bacteria blocks should use this order for consistency:

```rust
// --- ESBL Enzymes (Beta-lactam resistance) ---
enzyme_esbl_ctx_m
enzyme_esbl_tem
enzyme_esbl_shv

// --- AmpC Enzymes (Cephalosporin resistance) ---
enzyme_ampc_cmy
enzyme_ampc_dha

// --- Carbapenemases (Carbapenem resistance) ---
enzyme_kpc
enzyme_ndm_vim
enzyme_oxa_48

// --- Other Enzymatic ---
enzyme_cat           // Chloramphenicol
enzyme_16s_rrmt      // Aminoglycosides

// --- Target Site Modifications ---
target_site_pbp2a_meca    // Methicillin/MRSA
target_site_van_a         // Vancomycin/VRE
target_site_van_b         // Vancomycin/VRE variant
target_site_erm_b         // Macrolides
target_site_cfr           // Linezolid/Chloramphenicol

// --- Fluoroquinolone Resistance ---
mutation_gyra_primary
mutation_gyra_parc_secondary
protection_qnr

// --- Efflux Pumps ---
efflux_acrab_tolc
efflux_mexxy_oprm
global_efflux_pump

// --- Porin Loss ---
porin_loss_ompk35_36
porin_loss_oprd
global_porin_loss

// --- Surface Modifications ---
modification_mcr_1   // Colistin
```

---

## Impact of Missing Mechanisms

### High Impact:
- **mutation_gyra_parc_secondary**: Affects fluoroquinolone resistance evolution (2nd-step mutations)
- **enzyme_esbl_tem/shv**: Missing major ESBL variants (historical TEM, Klebsiella-specific SHV)
- **target_site_van_b**: VRE variant with different teicoplanin susceptibility

### Medium Impact:
- **enzyme_ampc_dha**: Alternative AmpC variant (less common than CMY)
- **target_site_cfr**: Linezolid resistance (rare but clinically important)
- **efflux_mexxy_oprm**: Pseudomonas-specific multidrug efflux
- **porin_loss_oprd**: Pseudomonas carbapenem resistance

### Lower Impact (Generic mechanisms):
- **global_efflux_pump**: Catch-all efflux
- **global_porin_loss**: Catch-all porin

---

## Implementation Plan

### Step 1: Generate complete mechanism blocks for all 42 bacteria
- Use standard order
- Include all 25 mechanisms
- Set appropriate default multipliers:
  - 1.0 for E. coli (reference organism)
  - 1000.0 for most other bacteria (effectively disables emergence unless specifically configured)
  - Lower values (10-100) for clinically important mechanisms in specific pathogens

### Step 2: Replace existing blocks in config.rs
- Find each bacteria's mechanism block
- Replace with complete, standardized block
- Preserve any custom multiplier values that differ from 1.0/1000.0

### Step 3: Verification
- Compile to ensure parameter keys match
- Run short test to verify mechanisms are loaded
- Check resistance trajectories haven't changed drastically

---

## Example: Complete Block for E. coli

```rust
// E. coli - Complete mechanism emergence multipliers (standard order)
map.insert("bacteria_escherichia_coli_mechanism_enzyme_esbl_ctx_m_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_esbl_tem_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_esbl_shv_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_ampc_cmy_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_ampc_dha_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_kpc_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_ndm_vim_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_oxa_48_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_cat_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_enzyme_16s_rrmt_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_target_site_pbp2a_meca_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_target_site_van_a_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_target_site_van_b_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_target_site_erm_b_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_target_site_cfr_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_mutation_gyra_primary_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_mutation_gyra_parc_secondary_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_protection_qnr_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_efflux_acrab_tolc_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_efflux_mexxy_oprm_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_global_efflux_pump_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_porin_loss_ompk35_36_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_porin_loss_oprd_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_global_porin_loss_emergence_multiplier".to_string(), 1.0);
map.insert("bacteria_escherichia_coli_mechanism_modification_mcr_1_emergence_multiplier".to_string(), 1.0);
```

---

## Bacteria Requiring Updates

All 42 bacteria need their mechanism blocks updated:
1. escherichia_coli
2. klebsiella_pneumoniae
3. enterobacter_cloacae
4. enterobacter_spp.
5. citrobacter_spp.
6. serratia_spp.
7. morganella_spp.
8. proteus_spp.
9. salmonella_enterica_serovar_typhi
10. salmonella_enterica_serovar_paratyphi_a
11. invasive_non-typhoidal_salmonella_spp.
12. shigella_spp.
13. campylobacter_jejuni
14. helicobacter_pylori
15. pseudomonas_aeruginosa
16. acinetobacter_baumannii
17. stenotrophomonas_maltophilia
18. burkholderia_cepacia_complex
19. staphylococcus_aureus
20. staphylococcus_epidermidis
21. enterococcus_faecalis
22. enterococcus_faecium
23. streptococcus_pneumoniae
24. streptococcus_pyogenes
25. streptococcus_agalactiae
26. haemophilus_influenzae
27. neisseria_meningitidis
28. neisseria_gonorrhoeae
29. moraxella_catarrhalis
30. clostridioides_difficile
31. bacteroides_fragilis
32. vibrio_cholerae
33. treponema_pallidum
34. chlamydia_trachomatis
35. mycoplasma_pneumoniae
36. mycobacterium_tuberculosis
37. mdr_mycobacterium_tuberculosis
38. listeria_monocytogenes
39. legionella_pneumophila
40. bordetella_pertussis
41. yersinia_pestis
42. francisella_tularensis

---

## Files to Modify

- `src/config.rs`: Lines ~10535-11500 (mechanism emergence multiplier blocks)

---

## Testing Checklist

After implementation:
- [ ] `cargo build --release` compiles successfully
- [ ] No new warnings about unused parameters
- [ ] Resistance mechanism prevalence cache loads correctly
- [ ] Run short simulation (1K individuals, 1000 timesteps) to verify:
  - [ ] ESBL emergence in E. coli still works
  - [ ] MRSA emergence in S. aureus still works
  - [ ] Fluoroquinolone resistance still evolves
  - [ ] No crashes or panics
- [ ] Check log output for "mechanism not found" errors

---

## Priority

**High Priority** - This affects:
- Resistance evolution dynamics
- Cross-resistance patterns
- Horizontal gene transfer logic
- Drug selection algorithm (resistance-aware prescribing)

The missing mechanisms may cause simulation to:
- Under-estimate resistance prevalence
- Miss important cross-resistance patterns (e.g., 2nd-step FQ resistance)
- Fail to properly model specific pathogens (Pseudomonas, VRE variants)

---

## Next Steps

1. Generate complete standardized blocks for all 42 bacteria
2. Create multi_replace_string_in_file operations to update config.rs
3. Compile and test
4. Document changes in changelog
5. Consider adding automated test to prevent future inconsistencies
