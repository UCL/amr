# Enums and Constants

> **Source Files**: 
> - `src/simulation/population.rs` → BACTERIA_LIST, DRUG_SHORT_NAMES, Individual
> - `src/config.rs` → additional constants

This document provides a reference for all enumerated types, constant lists, and magic numbers in the simulation.

---

## Table of Contents
1. [Bacteria List](#bacteria-list)
2. [Drug List](#drug-list)
3. [Resistance Mechanisms](#resistance-mechanisms)
4. [Syndromes](#syndromes)
5. [Regions](#regions)
6. [Hospital Status](#hospital-status)
7. [Drug Introduction Dates](#drug-introduction-dates)
8. [Array Sizes](#array-sizes)
9. [Threshold Constants](#threshold-constants)
10. [Cross-Reference Tables](#cross-reference-tables)

---

## Bacteria List

### Complete List (39 bacteria)

| Index | Name | Full Name | Gram | Category |
|-------|------|-----------|------|----------|
| 0 | `e_coli` | Escherichia coli | - | Enterobacteriaceae |
| 1 | `k_pneumoniae` | Klebsiella pneumoniae | - | Enterobacteriaceae |
| 2 | `k_oxytoca` | Klebsiella oxytoca | - | Enterobacteriaceae |
| 3 | `e_cloacae` | Enterobacter cloacae | - | Enterobacteriaceae |
| 4 | `e_aerogenes` | Enterobacter aerogenes | - | Enterobacteriaceae |
| 5 | `c_freundii` | Citrobacter freundii | - | Enterobacteriaceae |
| 6 | `c_koseri` | Citrobacter koseri | - | Enterobacteriaceae |
| 7 | `s_marcescens` | Serratia marcescens | - | Enterobacteriaceae |
| 8 | `p_mirabilis` | Proteus mirabilis | - | Enterobacteriaceae |
| 9 | `p_vulgaris` | Proteus vulgaris | - | Enterobacteriaceae |
| 10 | `m_morganii` | Morganella morganii | - | Enterobacteriaceae |
| 11 | `p_stuartii` | Providencia stuartii | - | Enterobacteriaceae |
| 12 | `p_aeruginosa` | Pseudomonas aeruginosa | - | Non-fermenter |
| 13 | `a_baumannii` | Acinetobacter baumannii | - | Non-fermenter |
| 14 | `s_maltophilia` | Stenotrophomonas maltophilia | - | Non-fermenter |
| 15 | `h_influenzae` | Haemophilus influenzae | - | HACEK |
| 16 | `n_meningitidis` | Neisseria meningitidis | - | Neisseria |
| 17 | `n_gonorrhoeae` | Neisseria gonorrhoeae | - | Neisseria |
| 18 | `m_catarrhalis` | Moraxella catarrhalis | - | Moraxella |
| 19 | `b_fragilis` | Bacteroides fragilis | - | Anaerobe |
| 20 | `s_aureus` | Staphylococcus aureus | + | Staphylococcus |
| 21 | `s_epidermidis` | Staphylococcus epidermidis | + | Staphylococcus |
| 22 | `s_lugdunensis` | Staphylococcus lugdunensis | + | Staphylococcus |
| 23 | `s_saprophyticus` | Staphylococcus saprophyticus | + | Staphylococcus |
| 24 | `s_pneumoniae` | Streptococcus pneumoniae | + | Streptococcus |
| 25 | `s_pyogenes` | Streptococcus pyogenes | + | Streptococcus |
| 26 | `s_agalactiae` | Streptococcus agalactiae | + | Streptococcus |
| 27 | `s_anginosus` | Streptococcus anginosus | + | Streptococcus |
| 28 | `e_faecalis` | Enterococcus faecalis | + | Enterococcus |
| 29 | `e_faecium` | Enterococcus faecium | + | Enterococcus |
| 30 | `c_difficile` | Clostridioides difficile | + | Anaerobe |
| 31 | `l_monocytogenes` | Listeria monocytogenes | + | Listeria |
| 32 | `m_tuberculosis` | Mycobacterium tuberculosis | + | Mycobacteria |
| 33 | `m_avium` | Mycobacterium avium | + | Mycobacteria |
| 34 | `c_pneumoniae` | Chlamydia pneumoniae | - | Atypical |
| 35 | `m_pneumoniae` | Mycoplasma pneumoniae | - | Atypical |
| 36 | `l_pneumophila` | Legionella pneumophila | - | Atypical |
| 37 | `h_pylori` | Helicobacter pylori | - | Helicobacter |
| 38 | (reserved) | - | - | - |

### Code Definition

```rust
pub const BACTERIA_LIST: &[&str] = &[
    "e_coli", "k_pneumoniae", "k_oxytoca", "e_cloacae", "e_aerogenes",
    "c_freundii", "c_koseri", "s_marcescens", "p_mirabilis", "p_vulgaris",
    "m_morganii", "p_stuartii",
    "p_aeruginosa", "a_baumannii", "s_maltophilia",
    "h_influenzae", "n_meningitidis", "n_gonorrhoeae", "m_catarrhalis",
    "b_fragilis",
    "s_aureus", "s_epidermidis", "s_lugdunensis", "s_saprophyticus",
    "s_pneumoniae", "s_pyogenes", "s_agalactiae", "s_anginosus",
    "e_faecalis", "e_faecium",
    "c_difficile", "l_monocytogenes",
    "m_tuberculosis", "m_avium",
    "c_pneumoniae", "m_pneumoniae", "l_pneumophila",
    "h_pylori",
];

pub const N_BACTERIA: usize = 39;
```

### Lookup Functions

```rust
pub fn get_bacteria_index(name: &str) -> Option<usize> {
    BACTERIA_LIST.iter().position(|&b| b == name)
}

pub fn get_bacteria_name(index: usize) -> Option<&'static str> {
    BACTERIA_LIST.get(index).copied()
}
```

---

## Drug List

### Complete List (52 drugs)

| Index | Short Name | Full Name | Class |
|-------|------------|-----------|-------|
| 0 | `sulfanilamide` | Sulfanilamide | Sulfonamide |
| 1 | `penicilling` | Penicillin G | Penicillin |
| 2 | `ampicillin` | Ampicillin | Penicillin |
| 3 | `amoxicillin` | Amoxicillin | Penicillin |
| 4 | `piperacillin` | Piperacillin | Penicillin |
| 5 | `ticarcillin` | Ticarcillin | Penicillin |
| 6 | `cephalexin` | Cephalexin | Cephalosporin 1st |
| 7 | `cefazolin` | Cefazolin | Cephalosporin 1st |
| 8 | `cefuroxime` | Cefuroxime | Cephalosporin 2nd |
| 9 | `ceftriaxone` | Ceftriaxone | Cephalosporin 3rd |
| 10 | `ceftazidime` | Ceftazidime | Cephalosporin 3rd |
| 11 | `cefepime` | Cefepime | Cephalosporin 4th |
| 12 | `ceftaroline` | Ceftaroline | Cephalosporin 5th |
| 13 | `meropenem` | Meropenem | Carbapenem |
| 14 | `imipenem_c` | Imipenem-Cilastatin | Carbapenem |
| 15 | `ertapenem` | Ertapenem | Carbapenem |
| 16 | `aztreonam` | Aztreonam | Monobactam |
| 17 | `erythromycin` | Erythromycin | Macrolide |
| 18 | `azithromycin` | Azithromycin | Macrolide |
| 19 | `clarithromycin` | Clarithromycin | Macrolide |
| 20 | `clindamycin` | Clindamycin | Lincosamide |
| 21 | `gentamicin` | Gentamicin | Aminoglycoside |
| 22 | `tobramycin` | Tobramycin | Aminoglycoside |
| 23 | `amikacin` | Amikacin | Aminoglycoside |
| 24 | `ciprofloxacin` | Ciprofloxacin | Fluoroquinolone |
| 25 | `levofloxacin` | Levofloxacin | Fluoroquinolone |
| 26 | `moxifloxacin` | Moxifloxacin | Fluoroquinolone |
| 27 | `ofloxacin` | Ofloxacin | Fluoroquinolone |
| 28 | `tetracycline` | Tetracycline | Tetracycline |
| 29 | `doxycycline` | Doxycycline | Tetracycline |
| 30 | `minocycline` | Minocycline | Tetracycline |
| 31 | `vancomycin` | Vancomycin | Glycopeptide |
| 32 | `teicoplanin` | Teicoplanin | Glycopeptide |
| 33 | `dalbavancin` | Dalbavancin | Glycopeptide |
| 34 | `linezolid` | Linezolid | Oxazolidinone |
| 35 | `tedizolid` | Tedizolid | Oxazolidinone |
| 36 | `quinu_dalfo` | Quinupristin-Dalfopristin | Streptogramin |
| 37 | `trim_sulf` | Trimethoprim-Sulfamethoxazole | Sulfonamide combo |
| 38 | `chlorampheni` | Chloramphenicol | Phenicol |
| 39 | `nitrofurantoin` | Nitrofurantoin | Nitrofuran |
| 40 | `retapamulin` | Retapamulin | Pleuromutilin |
| 41 | `fusidic_a` | Fusidic Acid | Fusidane |
| 42 | `metronidazole` | Metronidazole | Nitroimidazole |
| 43 | `furazolidone` | Furazolidone | Nitrofuran |
| 44 | `rifampicin` | Rifampicin | Rifamycin |
| 45 | `amoxicillin_clavulanate` | Amoxicillin-Clavulanate | BL/BLI |
| 46 | `piperacillin_tazobactam` | Piperacillin-Tazobactam | BL/BLI |
| 47 | `ampicillin_sulbactam` | Ampicillin-Sulbactam | BL/BLI |
| 48 | `ticarcillin_clavulanate` | Ticarcillin-Clavulanate | BL/BLI |
| 49 | `ceftazidime_avibactam` | Ceftazidime-Avibactam | BL/BLI |
| 50 | `meropenem_vaborbactam` | Meropenem-Vaborbactam | BL/BLI |
| 51 | `colistin` | Colistin | Polymyxin |

### Code Definition

```rust
pub const DRUG_SHORT_NAMES: &[&str] = &[
    "sulfanilamide",
    "penicilling", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin",
    "cephalexin", "cefazolin", "cefuroxime",
    "ceftriaxone", "ceftazidime", "cefepime", "ceftaroline",
    "meropenem", "imipenem_c", "ertapenem",
    "aztreonam",
    "erythromycin", "azithromycin", "clarithromycin",
    "clindamycin",
    "gentamicin", "tobramycin", "amikacin",
    "ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin",
    "tetracycline", "doxycycline", "minocycline",
    "vancomycin", "teicoplanin", "dalbavancin",
    "linezolid", "tedizolid",
    "quinu_dalfo", "trim_sulf", "chlorampheni", "nitrofurantoin",
    "retapamulin", "fusidic_a", "metronidazole", "furazolidone", "rifampicin",
    "amoxicillin_clavulanate", "piperacillin_tazobactam",
    "ampicillin_sulbactam", "ticarcillin_clavulanate",
    "ceftazidime_avibactam", "meropenem_vaborbactam",
    "colistin",
];

pub const N_DRUGS: usize = 52;
```

### Drug Classes

```rust
pub const PENICILLINS: &[&str] = &[
    "penicilling", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin",
    "amoxicillin_clavulanate", "ampicillin_sulbactam", 
    "piperacillin_tazobactam", "ticarcillin_clavulanate"
];

pub const CEPHALOSPORINS: &[&str] = &[
    "cephalexin", "cefazolin", "cefuroxime",
    "ceftriaxone", "ceftazidime", "cefepime", "ceftaroline",
    "ceftazidime_avibactam"
];

pub const CARBAPENEMS: &[&str] = &[
    "meropenem", "imipenem_c", "ertapenem", "meropenem_vaborbactam"
];

pub const FLUOROQUINOLONES: &[&str] = &[
    "ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"
];

pub const AMINOGLYCOSIDES: &[&str] = &[
    "gentamicin", "tobramycin", "amikacin"
];

pub const GLYCOPEPTIDES: &[&str] = &[
    "vancomycin", "teicoplanin", "dalbavancin"
];
```

---

## Resistance Mechanisms

### Mechanism List

| Index | Name | Mobile? | Confers Resistance To |
|-------|------|---------|----------------------|
| 0 | `esbl` | Yes | 3rd gen cephalosporins |
| 1 | `ampc` | Some | Cephalosporins |
| 2 | `carbapenemase_kpc` | Yes | Carbapenems |
| 3 | `carbapenemase_ndm` | Yes | Carbapenems |
| 4 | `carbapenemase_oxa` | Yes | Carbapenems |
| 5 | `meca` | No | β-lactams (MRSA) |
| 6 | `vana` | Yes | Vancomycin |
| 7 | `vanb` | Yes | Vancomycin |
| 8 | `gyra_mutation` | No | Fluoroquinolones |
| 9 | `qnr` | Yes | Fluoroquinolones |
| 10 | `mcr` | Yes | Colistin |
| 11 | `aac` | Yes | Aminoglycosides |
| 12 | `erm` | Yes | Macrolides |
| 13 | `tet_efflux` | Yes | Tetracyclines |

### Code Definition

```rust
pub const MECHANISMS: &[&str] = &[
    "esbl", "ampc", "carbapenemase_kpc", "carbapenemase_ndm", "carbapenemase_oxa",
    "meca", "vana", "vanb", "gyra_mutation", "qnr", "mcr", "aac", "erm", "tet_efflux"
];

pub const N_MECHANISMS: usize = 14;
```

---

## Syndromes

### Clinical Syndromes

| Index | Name | Description |
|-------|------|-------------|
| 0 | `UTI` | Urinary tract infection |
| 1 | `Pneumonia` | Lower respiratory infection |
| 2 | `SkinSoftTissue` | Cellulitis, abscess, wound infection |
| 3 | `Bacteremia` | Bloodstream infection |
| 4 | `GI` | Gastrointestinal infection |
| 5 | `Meningitis` | Central nervous system infection |
| 6 | `BoneJoint` | Osteomyelitis, septic arthritis |
| 7 | `Endocarditis` | Heart valve infection |
| 8 | `Other` | Other/unspecified |

### Code Definition

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Syndrome {
    UTI = 0,
    Pneumonia = 1,
    SkinSoftTissue = 2,
    Bacteremia = 3,
    GI = 4,
    Meningitis = 5,
    BoneJoint = 6,
    Endocarditis = 7,
    Other = 8,
}

pub const N_SYNDROMES: usize = 9;
```

---

## Regions

### Geographic Regions

| Index | Name | Description |
|-------|------|-------------|
| 0 | `north_america` | USA, Canada |
| 1 | `europe_west` | Western Europe |
| 2 | `europe_east` | Eastern Europe |
| 3 | `asia_east` | China, Japan, Korea |
| 4 | `asia_south` | India, Pakistan, Bangladesh |
| 5 | `asia_southeast` | Thailand, Vietnam, etc. |
| 6 | `africa_north` | North Africa, Middle East |
| 7 | `africa_sub` | Sub-Saharan Africa |
| 8 | `south_america` | Latin America |
| 9 | `oceania` | Australia, New Zealand |

### Code Definition

```rust
pub const REGIONS: &[&str] = &[
    "north_america", "europe_west", "europe_east",
    "asia_east", "asia_south", "asia_southeast",
    "africa_north", "africa_sub",
    "south_america", "oceania"
];

pub const N_REGIONS: usize = 10;
```

---

## Hospital Status

### Status Values

| Value | Meaning |
|-------|---------|
| 0 | Not hospitalized (community) |
| 1 | Hospitalized - General ward |
| 2 | Hospitalized - Surgical ward |
| 3 | Hospitalized - Medical ward |
| 4 | Hospitalized - ICU |
| 5 | Long-term care facility |

### Code Definition

```rust
// Simple integer representation
// 0 = community, >0 = some form of healthcare
pub type HospitalStatus = u8;

pub const HOSPITAL_COMMUNITY: u8 = 0;
pub const HOSPITAL_GENERAL: u8 = 1;
pub const HOSPITAL_SURGICAL: u8 = 2;
pub const HOSPITAL_MEDICAL: u8 = 3;
pub const HOSPITAL_ICU: u8 = 4;
pub const HOSPITAL_LTCF: u8 = 5;
```

---

## Drug Introduction Dates

### Historical Timeline

```rust
lazy_static! {
    pub static ref DRUG_INTRODUCTION_DATES: HashMap<&'static str, usize> = {
        let mut map = HashMap::new();
        // Format: (drug_name, time_step)
        // time_step 0 = Jan 1, 1930
        // time_step 365 = Jan 1, 1931
        
        // 1930s
        map.insert("sulfanilamide", 2555);     // 1937
        
        // 1940s
        map.insert("penicilling", 4380);       // 1942
        map.insert("chlorampheni", 6570);      // 1948
        
        // 1950s
        map.insert("erythromycin", 8030);      // 1952
        map.insert("vancomycin", 9490);        // 1956
        map.insert("tetracycline", 8395);      // 1953
        
        // 1960s
        map.insert("ampicillin", 11315);       // 1961
        map.insert("gentamicin", 13140);       // 1966
        map.insert("cephalexin", 14235);       // 1969
        
        // 1970s
        map.insert("amoxicillin", 15695);      // 1973
        map.insert("cefazolin", 15330);        // 1972
        map.insert("trim_sulf", 14600);        // 1970
        map.insert("metronidazole", 16060);    // 1974
        
        // 1980s
        map.insert("ceftriaxone", 19345);      // 1983
        map.insert("ceftazidime", 19710);      // 1984
        map.insert("imipenem_c", 20075);       // 1985
        map.insert("ciprofloxacin", 20805);    // 1987
        map.insert("azithromycin", 21535);     // 1989
        map.insert("amoxicillin_clavulanate", 19345);  // 1983
        
        // 1990s
        map.insert("meropenem", 23360);        // 1994
        map.insert("levofloxacin", 24090);     // 1996
        map.insert("piperacillin_tazobactam", 23360);  // 1994
        
        // 2000s
        map.insert("linezolid", 25550);        // 2000
        map.insert("daptomycin", 26645);       // 2003
        map.insert("tigecycline", 27375);      // 2005
        map.insert("ceftaroline", 28835);      // 2010
        
        // 2010s
        map.insert("ceftazidime_avibactam", 31025);   // 2015
        map.insert("meropenem_vaborbactam", 31755);   // 2017
        
        map
    };
}
```

---

## Array Sizes

### Core Constants

```rust
pub const N_BACTERIA: usize = 39;
pub const N_DRUGS: usize = 52;
pub const N_MECHANISMS: usize = 14;
pub const N_SYNDROMES: usize = 9;
pub const N_REGIONS: usize = 10;
```

### Array Dimensions

| Array | Dimensions | Total Elements |
|-------|------------|----------------|
| `level` | `[N_BACTERIA]` | 39 |
| `symptoms` | `[N_BACTERIA]` | 39 |
| `resistances` | `[N_BACTERIA][N_DRUGS]` | 2,028 |
| `microbiome_r` | `[N_BACTERIA][N_DRUGS]` | 2,028 |
| `mechanisms` | `[N_BACTERIA][N_MECHANISMS]` | 546 |
| `cur_use_drug` | `[N_DRUGS]` | 52 |
| `cur_level_drug` | `[N_DRUGS]` | 52 |

---

## Threshold Constants

### Infection Thresholds

```rust
pub const INFECTION_EPS: f64 = 0.01;         // Below = cleared
pub const SYMPTOM_THRESHOLD: f64 = 0.1;      // Level for symptoms
pub const SHEDDING_THRESHOLD: f64 = 0.05;    // Level to transmit
pub const MAX_INFECTION_LEVEL: f64 = 10.0;   // Cap
pub const INITIAL_INFECTION_LEVEL: f64 = 0.1;  // Starting level
```

### Resistance Thresholds

```rust
pub const MAX_RESISTANCE_LEVEL: f64 = 1.0;   // Cap
pub const RESISTANCE_EPS: f64 = 0.001;       // Negligible
pub const MAJORITY_R_THRESHOLD: f64 = 0.5;   // "Majority" definition
```

### Drug Thresholds

```rust
pub const THERAPEUTIC_DRUG_LEVEL: f64 = 3.0;  // Effective
pub const SUB_THERAPEUTIC_LEVEL: f64 = 0.1;   // Below = not effective
pub const STANDARD_DOSE_LEVEL: f64 = 10.0;    // Normal dose
```

### Time Constants

```rust
pub const DAYS_PER_YEAR: usize = 365;
pub const SIMULATION_START_YEAR: usize = 1930;
pub const TREATMENT_FAILURE_ASSESSMENT_DAY: usize = 3;
pub const REINFECTION_PROTECTION_WINDOW: usize = 30;  // days
```

---

## Cross-Reference Tables

### Bacteria to Common Syndromes

| Bacteria | Primary | Secondary | Tertiary |
|----------|---------|-----------|----------|
| e_coli | UTI | Bacteremia | GI |
| k_pneumoniae | Pneumonia | UTI | Bacteremia |
| s_aureus | SkinSoftTissue | Bacteremia | Pneumonia |
| s_pneumoniae | Pneumonia | Meningitis | Bacteremia |
| p_aeruginosa | Pneumonia | UTI | Bacteremia |
| e_faecalis | UTI | Bacteremia | Endocarditis |
| c_difficile | GI | - | - |

### Drug to Spectrum (Simplified)

| Drug | Gram+ | Gram- | Anaerobes | Atypicals |
|------|-------|-------|-----------|-----------|
| Penicillin G | +++ | - | + | - |
| Ampicillin | ++ | ++ | + | - |
| Ceftriaxone | ++ | +++ | - | - |
| Meropenem | +++ | +++ | +++ | - |
| Vancomycin | +++ | - | + | - |
| Ciprofloxacin | + | +++ | - | ++ |
| Azithromycin | ++ | + | - | +++ |
| Metronidazole | - | - | +++ | - |

### Mechanism to Drug Classes

| Mechanism | Affects |
|-----------|---------|
| ESBL | 3rd gen cephalosporins |
| KPC | All β-lactams incl. carbapenems |
| NDM | All β-lactams incl. carbapenems |
| MecA | All β-lactams (MRSA) |
| VanA | Vancomycin, teicoplanin |
| VanB | Vancomycin only |
| GyrA | Fluoroquinolones |
| MCR | Colistin |

---

## Index Lookup Utilities

```rust
// Bacteria
pub fn bacteria_idx(name: &str) -> Option<usize> {
    BACTERIA_LIST.iter().position(|&b| b == name)
}

// Drugs
pub fn drug_idx(name: &str) -> Option<usize> {
    DRUG_SHORT_NAMES.iter().position(|&d| d == name)
}

// Regions
pub fn region_idx(name: &str) -> Option<usize> {
    REGIONS.iter().position(|&r| r == name)
}

// Example usage
let ecoli = bacteria_idx("e_coli").unwrap();  // 0
let cipro = drug_idx("ciprofloxacin").unwrap();  // 24
let resistance = individual.resistances[ecoli][cipro].any_r;
```
