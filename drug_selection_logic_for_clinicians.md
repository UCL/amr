# Antibiotic Selection Logic for Clinicians

## Overview

This document describes how the antimicrobial resistance (AMR) simulation model selects antibiotics for patients. The model attempts to replicate realistic clinical prescribing patterns across different historical eras (1930–present), accounting for pathogen identification, regional resistance patterns, clinical syndromes, and antimicrobial stewardship principles.

**This is a reference document for clinicians.** It explains exactly how the model decides when to start antibiotics and which drugs to choose. All parameter values shown are the actual values used in the simulation.

---

## Two-Stage Selection Process

The model uses a **two-stage probabilistic process** for antibiotic selection:

### Stage 1: Decision to Initiate Antibiotic Therapy

Each day, the model calculates the probability that a patient will start antibiotic therapy. This probability starts at a **baseline rate** and is then modified by several clinical factors.

#### Baseline Probability

- **Daily baseline rate**: 0.1% per day (0.001)
  - This represents the small chance any patient might receive antibiotics (e.g., prophylaxis, suspected infections)

#### Clinical Factors That Increase Treatment Probability

| Factor | Multiplier | Meaning |
|--------|------------|---------|
| **Symptomatic infection present** | **260×** | If the patient has a detectable bacterial infection causing symptoms, they are 260 times more likely to start antibiotics that day |
| **Pathogen identified (culture positive)** | **2.5×** | If laboratory testing has confirmed the organism, an additional 2.5× increase applies (on top of the infection multiplier) |
| **Immunocompromised host** | **8×** | Immunosuppressed patients (e.g., transplant recipients, chemotherapy patients) are 8× more likely to receive prophylactic or empiric antibiotics |
| **Already on antibiotics** | **0.1×** (reduces) | If already receiving one antibiotic, 90% less likely to add another (to prevent poly-pharmacy unless needed) |

#### Syndrome-Specific Multipliers

Different infection syndromes have different treatment urgency. When a patient presents with a clinical syndrome, the following multipliers apply to the probability of starting antibiotics:

| Syndrome ID | Clinical Syndrome | Multiplier | Interpretation |
|-------------|-------------------|------------|----------------|
| 1 | Urinary tract infection (UTI) | 1.0× | Baseline (no change) |
| 2 | Skin / soft tissue infection | 1.0× | Baseline |
| 3 | Respiratory infection | **10.0×** | Strong push to treat—pneumonia requires rapid antibiotic initiation |
| 4 | Bloodstream infection (bacteremia) | 1.0× | Baseline (severity captured via hospitalization) |
| 5 | Intra-abdominal infection | 1.0× | Baseline |
| 6 | Central nervous system (meningitis) | 1.0× | Baseline (urgent but rare) |
| 7 | Gastrointestinal infection | **8.0×** | Often empirically treated (travelers' diarrhea, gastroenteritis) |
| 8 | Genital/pelvic infection | **12.0×** | STIs and pelvic infections often treated empirically |
| 9 | Bone/joint infection | 1.0× | Baseline |
| 10 | Other severe infection | 1.0× | Baseline |

#### Hospitalization Effect

- Hospitalized patients have higher baseline treatment rates due to severity of illness
- Hospital-based clinicians have more rapid access to IV antibiotics
- This is captured through the syndrome and infection presence multipliers

#### Maximum Concurrent Drugs

- **Maximum 3 drugs simultaneously**
- Once a patient is on 3 antibiotics, no new drugs can be added
- This reflects clinical practice for severe/complicated infections (e.g., TB treatment, MDR infections)

---

### Stage 2: Which Antibiotic to Select

Once the model decides to start treatment, it must choose *which* antibiotic to prescribe. This is done through a **weighted random selection** process:

1. Every available antibiotic receives a **score** based on multiple clinical factors
2. Higher-scoring drugs are more likely to be chosen
3. The selection is **probabilistic**, not deterministic—this mimics real clinical variability where different clinicians may make different reasonable choices

#### How Scores Become Probabilities

The model uses a mathematical formula to convert scores into selection probabilities. Here is the intuition:

- **If Drug A has score 10 and Drug B has score 5**, Drug A is roughly twice as likely to be selected as Drug B
- **Very high-scoring drugs dominate** the selection, but lower-scoring drugs may occasionally be chosen
- A "randomness parameter" (set to 1.1) controls how much variability exists—at this setting, preferred drugs are usually selected, but alternatives sometimes are too

**Example**: If three drugs have scores of 100, 50, and 10:
- Drug 1 (score 100) will be chosen ~60% of the time
- Drug 2 (score 50) will be chosen ~30% of the time
- Drug 3 (score 10) will be chosen ~10% of the time

---

## Drug Scoring Algorithm

Each antibiotic receives a score based on multiple factors applied **multiplicatively**. A drug starts with score = 1.0, then various multipliers are applied.

### 1. Pathogen-Specific Clinical Guidelines (Targeted Therapy)

When the causative organism has been **identified** (culture-positive), the model applies pathogen-specific scoring multipliers reflecting real-world clinical guidelines.

#### Gram-Positive Organisms

**Staphylococcus aureus**

| Drug | Score Multiplier | Notes |
|------|------------------|-------|
| Penicillin G | **25×** (before 1950), **2.5×** (after) | Strong preference early; MRSA reduces utility later |
| Amoxicillin/clavulanate | **18×** (before 1960), **3×** (after) | Beta-lactam with beta-lactamase cover |
| Ampicillin/sulbactam | **18×** (before 1960), **3×** (after) | Similar to amox/clav |
| Vancomycin | **1.5×** (before 1950), **18×** (after) | Reserved for MRSA era |
| Linezolid | **0.5×** (before 1960), **12×** (after) | Modern reserve agent |
| Tedizolid | **0.5×** (before 1960), **12×** (after) | Modern reserve agent |
| Clindamycin | **5×** | Consistently useful alternative |

**Staphylococcus epidermidis** (coagulase-negative staph, often device-associated)

| Drug | Score Multiplier |
|------|------------------|
| Vancomycin | **14×** |
| Linezolid | **10×** |
| Tedizolid | **10×** |
| Quinupristin/dalfopristin | **6×** |
| TMP-SMX | **4×** |
| Penicillin G, ampicillin, amoxicillin, cephalexin, cefazolin, ceftriaxone | **0.05×** (95% reduction—these drugs rarely work) |

**Streptococcus pneumoniae**

| Drug | Score Multiplier |
|------|------------------|
| Penicillin G | **24×** |
| Amoxicillin | **24×** |
| Ampicillin | **22×** |
| Amoxicillin/clavulanate | **12×** |
| Ampicillin/sulbactam | **12×** |
| Ceftriaxone | **6×** |
| Azithromycin | **6×** |
| Clarithromycin | **6×** |
| Meropenem, imipenem, colistin, linezolid, tedizolid | **0.15×** (85% reduction—avoid reserve drugs) |

**Streptococcus pyogenes** (Group A Strep)

| Drug | Score Multiplier |
|------|------------------|
| Penicillin G | **28×** |
| Ampicillin | **20×** |
| Amoxicillin | **20×** |
| Amoxicillin/clavulanate | **10×** |
| Meropenem, imipenem, colistin, linezolid, tedizolid | **0.1×** (avoid reserve drugs) |

**Enterococcus faecalis**

| Drug | Score Multiplier | Notes |
|------|------------------|-------|
| Ampicillin | **20×** | First-line treatment |
| Vancomycin | **5×** (before 1960), **12×** (VRE era) | Increases as VRE emerges |
| Linezolid | **2×** (before 1970), **10×** (VRE era) | Reserve for resistant strains |

**Enterococcus faecium** (inherently more resistant)

| Drug | Score Multiplier | Notes |
|------|------------------|-------|
| Ampicillin | **4×** | Less effective than for E. faecalis |
| Vancomycin | **8×** (before 1960), **15×** (VRE era) | |
| Linezolid | **3×** (before 1970), **12×** (VRE era) | |
| Quinupristin/dalfopristin | **10×** (from 1975 onward) | Late introduction |

#### Gram-Negative Organisms

**Pseudomonas aeruginosa**

| Drug | Score Multiplier | Notes |
|------|------------------|-------|
| Piperacillin/tazobactam | **12×** | First-line anti-pseudomonal |
| Ceftazidime | **10×** | Third-gen cephalosporin with activity |
| Cefepime | **10×** | Fourth-gen cephalosporin |
| Ciprofloxacin | **7×** | Fluoroquinolone with activity |
| Tobramycin | **7×** | Aminoglycoside of choice |
| Meropenem | **6×** | Carbapenem (reserve) |
| Imipenem | **5×** | Carbapenem (reserve) |
| Colistin | **4×** | Last resort |
| Penicillin G, ampicillin, amoxicillin, cephalexin, ceftriaxone, vancomycin | **Score = 0** (completely blocked—no intrinsic activity) |

**Escherichia coli**

| Drug | Score Multiplier | Notes |
|------|------------------|-------|
| Nitrofurantoin | **14×** | First-line for UTI |
| Ciprofloxacin | **12×** | Very commonly used |
| TMP-SMX | **10×** | First-line for UTI |
| Ceftriaxone | **9×** | For serious infections |
| Ampicillin | **15×** (before 1950), **4×** (after) | Resistance emerged early |
| Meropenem, imipenem | **0.3×** (before 1970), **6×** (ESBL era) | Reserved for resistant strains |

**Klebsiella pneumoniae**

| Drug | Score Multiplier | Notes |
|------|------------------|-------|
| Ceftriaxone | **10×** (before ESBL), **6×** (ESBL era) | Efficacy decreases as ESBL spreads |
| Meropenem, imipenem | **4×** (early), **12×** (ESBL era) | Becomes more important |
| Ciprofloxacin | **7×** | |
| Piperacillin/tazobactam | **9×** | |

**Acinetobacter baumannii**

| Drug | Score Multiplier | Notes |
|------|------------------|-------|
| Meropenem, imipenem | **12×** (before 2010), **6×** (after) | Resistance emerging |
| Colistin | **5×** (early), **10×** (MDR era) | Last resort |
| Ampicillin/sulbactam | **12×** | Unique activity |

**Stenotrophomonas maltophilia** (inherently resistant to carbapenems)

| Drug | Score Multiplier | Notes |
|------|------------------|-------|
| TMP-SMX | **14×** | Drug of choice |
| Minocycline | **10×** | Alternative |
| Doxycycline | **10×** | Alternative |
| Levofloxacin | **6×** | |
| Ciprofloxacin | **6×** | |
| Piperacillin/tazobactam, ceftazidime, meropenem, imipenem, aminoglycosides | **0.05×** (strongly penalized—intrinsically resistant) |

**Haemophilus influenzae**

| Drug | Score Multiplier |
|------|------------------|
| Amoxicillin/clavulanate | **14×** |
| Ampicillin/sulbactam | **12×** |
| Amoxicillin | **10×** |
| Ceftriaxone | **6×** |
| Cefuroxime | **6×** |
| Meropenem, imipenem, colistin | **0.25×** (avoid unnecessary broad spectrum) |

**Neisseria meningitidis**

| Drug | Score Multiplier |
|------|------------------|
| Penicillin G | **18×** |
| Ampicillin | **18×** |
| Ceftriaxone | **10×** |
| Cefepime | **10×** |
| Meropenem, imipenem, colistin, linezolid | **0.2×** |

#### Other Bacteria in the Model

The simulation includes 39 bacteria in total. Here is the complete list:

**Enterobacterales:**
- *Escherichia coli*
- *Klebsiella pneumoniae*
- *Enterobacter* spp.
- *Enterobacter cloacae*
- *Citrobacter* spp.
- *Serratia* spp.
- *Proteus* spp.
- *Morganella* spp.
- *Providencia stuartii*
- *Salmonella enterica* serovar Typhi
- *Salmonella enterica* serovar Paratyphi A
- Invasive non-typhoidal *Salmonella* spp.
- *Shigella* spp.
- *Yersinia enterocolitica*

**Non-fermenters:**
- *Pseudomonas aeruginosa*
- *Acinetobacter baumannii*
- *Stenotrophomonas maltophilia*

**Gram-positive cocci:**
- *Staphylococcus aureus*
- *Staphylococcus epidermidis*
- *Streptococcus pneumoniae*
- *Streptococcus pyogenes*
- *Streptococcus agalactiae*
- *Enterococcus faecalis*
- *Enterococcus faecium*

**Other Gram-negatives:**
- *Haemophilus influenzae*
- *Neisseria meningitidis*
- *Neisseria gonorrhoeae*
- *Moraxella catarrhalis*
- *Campylobacter jejuni*
- *Vibrio cholerae*
- *Bordetella pertussis*

**Atypicals:**
- *Chlamydia trachomatis*
- *Mycoplasma genitalium*
- *Treponema pallidum*
- *Helicobacter pylori*

**Anaerobes:**
- *Bacteroides fragilis*
- *Clostridioides difficile*

**Other:**
- *Listeria monocytogenes*
- MDR *Mycobacterium tuberculosis*

For bacteria without specific guidelines listed above, the model uses generic drug activity (potency) data to guide selection.

#### Carbapenem Stewardship Penalty

Carbapenems (meropenem, imipenem, ertapenem, meropenem/vaborbactam) are powerful reserve agents. The model applies an **88% score reduction** (multiply by 0.12) when carbapenems are prescribed for organisms where they are not specifically indicated—such as community-acquired E. coli UTI.

**What this means**: If a carbapenem would otherwise have a score of 100, it is reduced to 12 when used for non-indicated organisms. This makes it far less likely to be selected, enforcing antimicrobial stewardship principles.

---

### 2. First-Line vs Off-Guideline Therapy

The model maintains pathogen-specific lists of **first-line and second-line drugs**. If a drug is NOT on these lists, it receives an **85% score reduction** (multiply by 0.15).

**What this means**: An off-guideline drug with score 100 becomes score 15, dramatically reducing its selection probability.

#### Complete First/Second-Line Drug Lists by Pathogen

| Pathogen | First/Second-Line Drugs |
|----------|------------------------|
| **P. aeruginosa** | Piperacillin/tazobactam, meropenem, imipenem, ceftazidime, cefepime, ciprofloxacin, tobramycin |
| **S. aureus** | Penicillin G, amoxicillin/clavulanate, ampicillin/sulbactam, vancomycin, linezolid, tedizolid, clindamycin |
| **S. epidermidis** | Vancomycin, linezolid, tedizolid, quinupristin/dalfopristin, TMP-SMX |
| **S. maltophilia** | TMP-SMX, minocycline, doxycycline, levofloxacin, ciprofloxacin |
| **S. pneumoniae** | Penicillin G, ampicillin, amoxicillin, amoxicillin/clavulanate, ceftriaxone, cefuroxime, azithromycin, clarithromycin |
| **S. pyogenes** | Penicillin G, ampicillin, amoxicillin, amoxicillin/clavulanate, clindamycin, azithromycin |
| **H. influenzae** | Amoxicillin, ampicillin, amoxicillin/clavulanate, ampicillin/sulbactam, cefuroxime, ceftriaxone |
| **N. meningitidis** | Penicillin G, ampicillin, ceftriaxone, cefepime |
| **E. coli** | Ciprofloxacin, nitrofurantoin, TMP-SMX, ceftriaxone, ampicillin, cefuroxime |
| **K. pneumoniae** | Ceftriaxone, ceftazidime, cefepime, piperacillin/tazobactam, ciprofloxacin |
| **E. faecalis** | Ampicillin, vancomycin, linezolid, tedizolid |
| **E. faecium** | Vancomycin, linezolid, tedizolid, quinupristin/dalfopristin |
| **A. baumannii** | Meropenem, imipenem, colistin, ampicillin/sulbactam, minocycline |

For other bacteria, no specific list restriction applies.

---

### 3. Drug Potency (Expected Activity)

For targeted therapy (when the pathogen is identified), drugs receive bonus scores based on their expected **in vitro potency** against the identified organism.

"Potency" is a value from 0.0 to 1.0, where:
- **1.0** = drug completely eliminates the bacteria at standard doses
- **0.5** = drug has moderate activity
- **0.0** = drug has no activity against this organism

| Potency Range | Score Multiplier | Meaning |
|---------------|------------------|---------|
| **≥0.50** (high potency) | **15×** | Excellent activity—strongly favored |
| **0.30 to 0.49** | **10×** | Good activity—favored |
| **0.15 to 0.29** | **6×** | Moderate activity—acceptable |
| **0.15** (minimal threshold) to 0.14 | **2×** | Marginal activity—only if no better options |
| **Below 0.15** | Drug excluded | Drug has insufficient activity and is not considered |

**What this means**: Drugs that work well against the identified pathogen receive dramatically higher scores, making them much more likely to be selected.

---

### 4. Historical Era / Time Period Effects

The model simulates realistic prescribing patterns across multiple antibiotic eras. Many pathogen-specific multipliers (shown above) change based on the simulated calendar year.

#### Key Antibiotic Eras

| Era | Approximate Years | Key Events |
|-----|-------------------|------------|
| **Pre-antibiotic** | Before 1937 | No antibiotics available |
| **Early sulfonamide** | 1937–1942 | Only sulfanilamide available |
| **Early penicillin** | 1942–1961 | Penicillin G dominant; limited drug choice |
| **Golden age** | 1961–1985 | Many new classes introduced; resistance rare |
| **MRSA era** | ~1950 onward | Methicillin-resistant S. aureus emerges; vancomycin use increases |
| **ESBL era** | ~1960 onward | Extended-spectrum beta-lactamases spread; carbapenem use increases |
| **VRE era** | ~1970 onward | Vancomycin-resistant enterococci emerge; linezolid becomes important |
| **MDR/XDR era** | ~2000 onward | Multi-drug resistant organisms; colistin and novel agents required |

#### Time Step Reference

The simulation uses "time steps" where each step = 1 day, starting from January 1, 1930.

| Simulation Year | Time Step (approximate) |
|-----------------|------------------------|
| 1930 | 0 |
| 1940 | 3,650 |
| 1950 | 7,300 |
| 1960 | 10,950 |
| 1970 | 14,600 |
| 1980 | 18,250 |
| 1990 | 21,900 |
| 2000 | 25,550 |
| 2010 | 29,200 |
| 2020 | 32,850 |
| 2035 | 38,325 |

#### Example Era-Specific Drug Scoring Changes

| Pathogen + Drug | Before Era Change | After Era Change | Inflection Year |
|-----------------|-------------------|------------------|-----------------|
| E. coli + ampicillin | **15×** | **4×** | ~1950 |
| E. coli + carbapenems | **0.3×** | **6×** | ~1970 |
| K. pneumoniae + ceftriaxone | **10×** | **6×** | ~1960 |
| K. pneumoniae + carbapenems | **4×** | **12×** | ~1960 |
| S. aureus + penicillin | **25×** | **2.5×** | ~1950 |
| S. aureus + vancomycin | **1.5×** | **18×** | ~1950 |
| S. aureus + linezolid | **0.5×** | **12×** | ~1960 |
| E. faecalis + vancomycin | **5×** | **12×** | ~1960 |
| E. faecalis + linezolid | **2×** | **10×** | ~1970 |
| E. faecium + linezolid | **3×** | **12×** | ~1970 |
| A. baumannii + carbapenems | **12×** | **6×** | ~2010 |
| A. baumannii + colistin | **5×** | **10×** | ~1970 |

---

### 5. Drug Introduction Year

Antibiotics only become available after their historical introduction date. Before that date, the drug **cannot be selected** (score = 0).

#### Complete Drug Introduction Timeline

| Drug | Model Name | Introduction Year |
|------|------------|-------------------|
| Sulfanilamide | sulfanilamide | 1937 |
| Penicillin G | penicilling | 1942 |
| Tetracycline | tetracycline | 1948 |
| Chloramphenicol | chlorampheni | 1949 |
| Erythromycin | erythromycin | 1952 |
| Colistin | colistin | 1952 |
| Nitrofurantoin | nitrofurantoin | 1953 |
| Furazolidone | furazolidone | 1955 |
| Vancomycin | vancomycin | 1958 |
| Metronidazole | metronidazole | 1960 |
| Ampicillin | ampicillin | 1961 |
| Fusidic acid | fusidic_a | 1962 |
| Gentamicin | gentamicin | 1963 |
| Rifampicin | rifampicin | 1966 |
| Doxycycline | doxycycline | 1967 |
| Clindamycin | clindamycin | 1968 |
| TMP-SMX | trim_sulf | 1968 |
| Cephalexin | cephalexin | 1970 |
| Minocycline | minocycline | 1971 |
| Amoxicillin | amoxicillin | 1972 |
| Cefazolin | cefazolin | 1973 |
| Tobramycin | tobramycin | 1975 |
| Amikacin | amikacin | 1976 |
| Ticarcillin | ticarcillin | 1977 |
| Cefuroxime | cefuroxime | 1978 |
| Piperacillin | piperacillin | 1981 |
| Ceftriaxone | ceftriaxone | 1984 |
| Piperacillin/tazobactam | piperacillin_tazobactam | 1984 |
| Imipenem/cilastatin | imipenem_c | 1985 |
| Ceftazidime | ceftazidime | 1985 |
| Amoxicillin/clavulanate | amoxicillin_clavulanate | 1985 |
| Aztreonam | aztreonam | 1986 |
| Ciprofloxacin | ciprofloxacin | 1987 |
| Teicoplanin | teicoplanin | 1988 |
| Ofloxacin | ofloxacin | 1990 |
| Clarithromycin | clarithromycin | 1990 |
| Ampicillin/sulbactam | ampicillin_sulbactam | 1990 |
| Ticarcillin/clavulanate | ticarcillin_clavulanate | 1990 |
| Azithromycin | azithromycin | 1991 |
| Meropenem | meropenem | 1996 |
| Cefepime | cefepime | 1996 |
| Levofloxacin | levofloxacin | 1996 |
| Moxifloxacin | moxifloxacin | 1999 |
| Quinupristin/dalfopristin | quinu_dalfo | 1999 |
| Linezolid | linezolid | 2000 |
| Ertapenem | ertapenem | 2001 |
| Ceftazidime/avibactam | ceftazidime_avibactam | 2006 |
| Retapamulin | retapamulin | 2007 |
| Ceftaroline | ceftaroline | 2010 |
| Tedizolid | tedizolid | 2014 |
| Dalbavancin | dalbavancin | 2014 |
| Meropenem/vaborbactam | meropenem_vaborbactam | 2018 |

---

## Empiric Therapy (Organism Unknown)

When no pathogen has been identified (no positive culture), the model uses **syndrome-based empiric scoring**. Each clinical syndrome has preset drug scores reflecting standard-of-care empiric regimens.

**How it works**: The score shown below is directly proportional to the probability of selection. If Drug A has score 18 and Drug B has score 9, Drug A is about twice as likely to be chosen.

### Complete Syndrome-Specific Empiric Drug Scores

#### Syndrome 1: Urinary Tract Infection (UTI)

| Drug | Score |
|------|-------|
| Nitrofurantoin | **18.0** |
| TMP-SMX | **14.0** |
| Ciprofloxacin | **12.0** |
| Amoxicillin/clavulanate | **12.0** |
| Levofloxacin | **10.0** |
| Amoxicillin | **10.0** |
| Ampicillin | **8.5** |
| Ceftriaxone | **8.0** |
| Cefuroxime | **7.0** |
| Piperacillin/tazobactam | **5.0** |
| Cefepime | **4.0** |
| Ceftazidime | **4.0** |
| Meropenem | **4.0** |
| Imipenem | **4.0** |
| Ertapenem | **4.0** |
| Meropenem/vaborbactam | **3.0** |
| Ceftazidime/avibactam | **3.0** |
| Colistin | **0.2** |
| Vancomycin | **0.1** |
| Linezolid | **0.1** |

#### Syndrome 2: Skin and Soft Tissue Infection

| Drug | Score |
|------|-------|
| Penicillin G | **16.0** |
| Amoxicillin/clavulanate | **14.0** |
| Amoxicillin | **14.0** |
| Ampicillin | **13.0** |
| Cephalexin | **13.0** |
| Cefazolin | **12.0** |
| Clindamycin | **12.0** |
| Vancomycin | **11.0** |
| Linezolid | **10.0** |
| Tedizolid | **9.0** |
| Dalbavancin | **9.0** |
| TMP-SMX | **9.0** |
| Doxycycline | **9.0** |
| Minocycline | **9.0** |
| Quinupristin/dalfopristin | **8.0** |
| Rifampicin | **6.0** |
| Ciprofloxacin | **4.0** |
| Piperacillin/tazobactam | **3.0** |

#### Syndrome 3: Respiratory Infection

| Drug | Score |
|------|-------|
| Amoxicillin/clavulanate | **16.0** |
| Amoxicillin | **15.5** |
| Penicillin G | **14.0** |
| Ampicillin | **13.5** |
| Azithromycin | **10.5** |
| Clarithromycin | **9.5** |
| Ceftriaxone | **9.5** |
| Cefuroxime | **8.5** |
| Piperacillin/tazobactam | **8.0** |
| Levofloxacin | **8.0** |
| Moxifloxacin | **8.0** |
| Cefepime | **7.5** |
| Erythromycin | **7.5** |
| Linezolid | **7.0** |
| Doxycycline | **6.5** |
| Vancomycin | **6.5** |
| Meropenem | **6.0** |
| Imipenem | **6.0** |
| Ofloxacin | **6.0** |
| Minocycline | **5.5** |

#### Syndrome 4: Bloodstream Infection (Bacteremia/Sepsis)

| Drug | Score |
|------|-------|
| Piperacillin/tazobactam | **14.0** |
| Meropenem | **13.0** |
| Imipenem | **13.0** |
| Meropenem/vaborbactam | **13.0** |
| Ceftazidime/avibactam | **12.5** |
| Cefepime | **12.0** |
| Ampicillin/sulbactam | **11.5** |
| Ceftazidime | **11.0** |
| Vancomycin | **11.0** |
| Amoxicillin/clavulanate | **10.5** |
| Ceftriaxone | **10.0** |
| Ampicillin | **10.0** |
| Linezolid | **10.0** |
| Amoxicillin | **9.5** |
| Tedizolid | **9.0** |
| Quinupristin/dalfopristin | **8.5** |
| Dalbavancin | **8.0** |
| Amikacin | **7.0** |
| Gentamicin | **7.0** |
| Tobramycin | **6.5** |
| Penicillin G | **6.5** |
| Ciprofloxacin | **6.0** |
| Colistin | **6.0** |
| Levofloxacin | **5.5** |
| Rifampicin | **4.0** |

#### Syndrome 5: Intra-abdominal Infection

| Drug | Score |
|------|-------|
| Metronidazole | **15.0** |
| Piperacillin/tazobactam | **13.0** |
| Meropenem | **13.0** |
| Ampicillin/sulbactam | **12.5** |
| Imipenem | **12.5** |
| Amoxicillin/clavulanate | **11.5** |
| Ertapenem | **11.0** |
| Ceftazidime/avibactam | **10.0** |
| Ceftriaxone | **9.0** |
| Ceftazidime | **9.0** |
| Cefepime | **9.0** |
| Meropenem/vaborbactam | **10.0** |
| Ciprofloxacin | **7.0** |
| Levofloxacin | **6.5** |
| TMP-SMX | **4.0** |
| Colistin | **3.5** |

#### Syndrome 6: Central Nervous System (Meningitis)

| Drug | Score |
|------|-------|
| Ceftriaxone | **15.0** |
| Vancomycin | **13.0** |
| Ampicillin | **13.0** |
| Ceftazidime | **12.0** |
| Cefepime | **12.0** |
| Penicillin G | **11.0** |
| Meropenem | **11.0** |
| Imipenem | **10.0** |
| Linezolid | **10.0** |
| Chloramphenicol | **9.0** |
| Rifampicin | **7.0** |
| Piperacillin/tazobactam | **6.0** |

#### Syndrome 7: Gastrointestinal Infection

| Drug | Score |
|------|-------|
| Ciprofloxacin | **12.0** |
| Metronidazole | **12.0** |
| Furazolidone | **11.0** |
| Azithromycin | **10.0** |
| Levofloxacin | **10.0** |
| Doxycycline | **8.5** |
| TMP-SMX | **8.5** |
| Amoxicillin/clavulanate | **7.0** |
| Minocycline | **6.5** |
| Amoxicillin | **6.5** |
| Ampicillin | **6.5** |
| Rifampicin | **5.0** |

#### Syndrome 8: Genital/Pelvic Infection

| Drug | Score |
|------|-------|
| Azithromycin | **13.0** |
| Ceftriaxone | **13.0** |
| Penicillin G | **12.0** |
| Doxycycline | **12.0** |
| Metronidazole | **12.0** |
| Cefuroxime | **9.0** |
| Amoxicillin/clavulanate | **9.5** |
| Amoxicillin | **9.0** |
| Clindamycin | **9.0** |
| Ciprofloxacin | **7.0** |
| Levofloxacin | **6.5** |
| TMP-SMX | **5.0** |
| Rifampicin | **4.0** |

#### Syndrome 9: Bone/Joint Infection

| Drug | Score |
|------|-------|
| Cefazolin | **13.0** |
| Vancomycin | **12.0** |
| Cephalexin | **11.0** |
| Ceftriaxone | **11.0** |
| Linezolid | **11.0** |
| Tedizolid | **10.0** |
| Dalbavancin | **10.0** |
| Clindamycin | **10.0** |
| Ciprofloxacin | **9.0** |
| Levofloxacin | **9.0** |
| Rifampicin | **9.0** |
| TMP-SMX | **8.0** |
| Meropenem | **7.0** |
| Piperacillin/tazobactam | **6.5** |

#### Syndrome 10: Other Severe Infection

| Drug | Score |
|------|-------|
| Piperacillin/tazobactam | **8.0** |
| Cefepime | **8.0** |
| Ceftriaxone | **8.0** |
| Meropenem | **8.0** |
| Imipenem | **8.0** |
| Vancomycin | **8.0** |
| Linezolid | **7.0** |
| Ciprofloxacin | **7.0** |
| Azithromycin | **6.0** |

### Drug Spectrum and Empiric Therapy

In addition to the syndrome-specific scores, the model considers **drug spectrum breadth** during empiric therapy:

- **Narrow-spectrum drugs** (spectrum score ≤2.0): Receive a modest penalty (**0.85×**) during empiric therapy
- **Broad-spectrum drugs** (spectrum score ≥3.5): Receive a small bonus (**0.85×** penalty removed)

**Why this matters**: When the organism is unknown, broader-spectrum agents provide better coverage of potential pathogens. However, syndrome-specific scores already capture most of this effect (e.g., carbapenems score higher for bacteremia than for UTI).

---

## Reserve Drug Gating (Antimicrobial Stewardship)

The model enforces **escalation requirements** for reserve antibiotics. These agents are heavily restricted to preserve their effectiveness against resistant organisms.

### Reserve Drug List

The following drugs are classified as "reserve" agents:
- **Carbapenems**: Meropenem, imipenem, ertapenem, meropenem/vaborbactam
- **Colistin** (polymyxin)
- **Oxazolidinones**: Linezolid, tedizolid
- **Quinupristin/dalfopristin**
- **Dalbavancin**

### Stewardship Requirements for Reserve Drugs

**For Targeted Therapy (pathogen identified):**

Even when the causative organism is identified, reserve drugs require **documented prior treatment failure** to be easily accessible.

- **Without prior failure**: Score is multiplied by **0.02** (a 98% reduction)
  - This means a drug that would otherwise have score 100 becomes score 2
  - The drug can still be selected, but it is 50 times less likely than it would otherwise be
- **With documented failure**: No penalty (full score)

**What counts as "documented failure"?** If a patient was on antibiotics and their infection did not improve (bacteria level remained high or increased) within the failure assessment window (typically 5 days), this is recorded as a treatment failure.

**For Empiric Therapy (pathogen unknown):**

Reserve drugs face even stricter requirements when the organism is not yet identified:

1. **Documented prior treatment failure is required** (same as above)
2. **AND** the regional resistance prevalence for that drug must be "high" (≥45%)

If BOTH conditions are not met, the reserve drug receives a score of **zero** and cannot be selected for empiric therapy.

**Clinical rationale**: This reflects real-world antimicrobial stewardship programs where carbapenems and other reserve agents require:
- Prior culture data showing resistance to first-line agents, OR
- Local antibiogram showing high resistance rates justifying empiric use

---

## Regional Resistance Surveillance

The model incorporates **regional resistance prevalence** into drug selection for empiric therapy. This simulates how clinicians adjust empiric regimens based on local antibiogram data.

### How Regional Resistance Affects Drug Selection

When prescribing empiric therapy (organism unknown), the model checks the resistance prevalence for each drug in the patient's region. If resistance is high, the drug receives a score penalty.

| Regional Resistance Level | Threshold | Score Penalty | Result |
|--------------------------|-----------|---------------|--------|
| **Very high** | ≥60% | **0.2×** (80% reduction) | Drug rarely selected |
| **High** | ≥45% | **0.4×** (60% reduction) | Drug discouraged |
| **Moderate** | ≥10% | **0.7×** (30% reduction) | Mild penalty |
| **Low** | <10% | **1.0×** (no penalty) | No change |

### Clinical Example

If ciprofloxacin resistance in E. coli is 50% in a particular region:
- Ciprofloxacin for empiric UTI therapy would receive a **0.4× penalty** (high resistance)
- Its score of 12.0 becomes 4.8
- This shifts prescribing toward alternatives like nitrofurantoin (score 18.0) or TMP-SMX (score 14.0)

---

## Resistance Testing Impact

When antimicrobial susceptibility testing (AST) results are available:

1. **Resistant isolate detected**: Drug is **completely excluded** from selection (score = 0)
2. **Susceptible isolate detected**: Drug remains a candidate (no penalty)
3. **Testing delay**: Results take approximately **3 days** to become available after culture is obtained

### Historical Testing Availability

Laboratory testing capabilities evolved over time:

| Technology | Year First Available |
|-----------|---------------------|
| Bacterial culture/identification | 1945 |
| Antibiotic susceptibility testing | 1955 |

Before 1945, no cultures are performed in the model. Before 1955, resistance testing is not available. This means targeted therapy based on AST results is not possible before these dates.

---

## Special Patient Populations

### Immunocompromised Hosts
- **8× higher baseline antibiotic initiation rate** (prophylaxis common)
- May receive antibiotics even without symptomatic infection
- More aggressive empiric therapy when infection suspected

### Hospitalized Patients
- **Higher testing rates**: 8× more likely to have bacterial cultures; 5× more likely to have AST performed
- Access to broader range of IV antibiotics
- Separate resistance tracking for hospital-acquired vs community-acquired infections

### Drug Allergies
- **Penicillin allergy**: All penicillin-class beta-lactams are excluded from selection
- Cross-reactivity with cephalosporins is modeled (some patients with penicillin allergy also cannot receive cephalosporins)

---

## Summary: What Determines Drug Choice?

The model chooses antibiotics through a systematic process. Here are the key questions, in order:

1. **Is the pathogen identified?**
   - **Yes** → Use targeted therapy with pathogen-specific guidelines and potency-based scoring
   - **No** → Use empiric therapy based on clinical syndrome

2. **What year is it in the simulation?**
   - Only drugs introduced by that year are available
   - Time-dependent scoring reflects evolving resistance patterns (e.g., vancomycin becomes more important as MRSA spreads)

3. **Has prior treatment failed?**
   - Required for easy access to reserve agents (carbapenems, colistin, linezolid, etc.)
   - Without prior failure, reserve drugs receive 98% score penalty

4. **What is the regional resistance profile?**
   - High local resistance (≥45%) penalizes empiric use of affected drugs by 60%
   - Very high resistance (≥60%) penalizes by 80%

5. **What AST results are available?**
   - Resistant isolates completely exclude that drug from selection

6. **What is the patient's clinical status?**
   - Syndrome type (UTI vs pneumonia vs bacteremia) determines base empiric scores
   - Hospitalization and immunosuppression increase treatment urgency

The final result is a **clinically plausible prescribing pattern** that evolves appropriately across the ~105-year simulation period (1930–2035), reflecting both therapeutic advances and the emergence of antimicrobial resistance.
