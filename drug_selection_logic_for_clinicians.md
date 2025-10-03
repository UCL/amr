# Drug Selection Logic in AMR Simulation Model
 

#### C. CLINICAL SCENARIO MODIFIERS

This document explains how the AMR simulation model makes decisions about when to start antibiotics and which drugs to prescribe. 

## Overview: How Treatment Decisions Are Made

The model uses a **two-stage decision process** for prescribing antibiotics:

1. **Stage 1**: Should we start ANY antibiotic at all?
2. **Stage 2**: If yes, WHICH specific antibiotic should we choose?


---

## Stage 1: Decision to Start Any Antibiotic

### Base Probability
- Every patient has a small baseline chance of receiving an antibiotic each day (even without infection)
- This represents prophylactic use, empirical treatment of suspected infections, etc.

### Factors That INCREASE the Likelihood of Starting Antibiotics:

**1. Active Infection Present**
- If bacteria levels are detectable (>0.001), probability increases significantly
- Multiplier applied: typically 5-10x higher chance

**2. Infection Identified by Testing**
- If diagnostic testing has confirmed the infection, probability increases further
- Represents the clinical reality that confirmed infections get treated more aggressively

**3. Patient Already on Antibiotics**
- If patient is already taking other antibiotics, there's reduced chance of adding more

**4. Immunocompromised Patients**
- Patients with immunodeficiency get prophylactic antibiotics much more frequently
- Multiplier: typically 8x higher chance
- Reflects clinical practice of aggressive prophylaxis in high-risk patients

**5. Clinical Syndrome**
- Different infection types (pneumonia, bloodstream infection, etc.) have different treatment thresholds
- Some syndromes prompt more aggressive antibiotic use

### Important Restrictions:
- **Maximum 3 drugs simultaneously**: Patients cannot be started on a 4th antibiotic
- **Must have available drugs**: If no antibiotics are available in the region/time period, no treatment possible

---

## Stage 2: Which Specific Antibiotic to Choose

If the decision is made to start an antibiotic, the model then scores every available drug and selects one using weighted probability (favoring higher-scoring drugs but allowing some randomness).

### Step 1: Filter Available Drugs

**Drugs are EXCLUDED if:**
- Not available in the patient's region at this time period
- Not yet historically introduced (e.g., no ceftriaxone in 1950)
- Resistance testing shows the bacteria is resistant to that drug
- Drug has insufficient activity against the infection (potency below 10% threshold)

### Step 2: Score Each Remaining Drug

Every eligible drug gets a score based on multiple factors:

#### A. INTRINSIC ACTIVITY WHEN NO RESISTANCE PRESENT (Potency)
- Drugs must have meaningful activity against the current infection
- Minimum threshold: usually 10% potency
- Higher potency against the infection = higher score

#### B. PATHOGEN-SPECIFIC CLINICAL GUIDELINES

The model incorporates evidence-based treatment preferences and heavily penalizes (95% score reduction) drugs that aren't first/second-line for the specific pathogen:

**Pseudomonas aeruginosa infections:**
- STRONGLY favors anti-pseudomonal agents:
  - Piperacillin-tazobactam: 25x multiplier
  - Meropenem: 25x multiplier  
  - Cefepime: 22x multiplier
  - Ceftazidime: 20x multiplier
- COMPLETELY BLOCKS inappropriate drugs:
  - Penicillin, ampicillin, cephalexin, vancomycin get 0 score
  - Reflects that these have no activity against Pseudomonas

**Staphylococcus aureus infections (time-dependent):**
- **Early era (pre-1970s - MSSA dominant):**
  - Penicillin: 50x multiplier (highly effective)
  - Methicillin: 40x multiplier
  - Vancomycin: 2x multiplier (rarely needed)
- **Later era (post-1970s - MRSA emergence):**
  - Vancomycin: 35x multiplier (becomes essential)
  - Penicillin: 2x multiplier (often ineffective)
  - Newer agents (linezolid, daptomycin): High multipliers

**MDR Tuberculosis (we only model MDR TB):**
- Model simulates MDR-TB with guaranteed 90% rifampicin resistance
- TB drugs in model: rifampicin, fluoroquinolones (ciprofloxacin, levofloxacin, moxifloxacin)
- *See Multi-Drug Therapy section below for TB-specific synergy details*



#### C. CLINICAL SCENARIO MODIFIERS

**Severe Infections:**
- ICU patients, sepsis, high bacteria levels
- Broad-spectrum agents get preference multipliers (2-5x)

**Age-Based Preferences:**
- **Pediatric patients:** fluoroquinolones generally avoided (safety concerns), beta-lactams preferred
- **Elderly patients:** potentially nephrotoxic drugs penalized (aminoglycosides: gentamicin, tobramycin, amikacin; colistin)
- Age-dependent sepsis risk multipliers applied for drug initiation decisions

**Hospital vs Community:**
- **Hospital-acquired infections favor broad-spectrum agents:** carbapenems (meropenem, imipenem), piperacillin-tazobactam, vancomycin, linezolid, colistin
- **Community infections favor narrower-spectrum agents:** beta-lactams (penicillin, ampicillin), ceftriaxone, fluoroquinolones (ciprofloxacin), trimethoprim-sulfamethoxazole

**Region-Specific Guidelines:**
- **Clinical preference multipliers:** Each bacteria-drug combination can have region-specific preference scores
- **Drug availability:** Time and region-aware availability (e.g., limited carbapenem access in low-resource settings)
- **Resistance surveillance:** Regional resistance data influences empirical therapy choices
- **Examples:** High vancomycin preference in MRSA-endemic regions, limited fluoroquinolone use in high-resistance areas

#### D. RESISTANCE CONSIDERATIONS

**Known Resistance:**
- If resistance testing shows resistance: drug gets 0 score (blocked)
- Recent treatment failures: temporary avoidance of recently failed drugs

**Resistance Risk:**
- High-resistance-risk drugs may be penalized in certain scenarios

### Step 3: Weighted Selection

- All scores are converted to probabilities using a "temperature" parameter
- Lower temperature = more deterministic (usually chooses highest-scoring drug)
- Higher temperature = more random selection
- Default setting favors best drugs but allows some variation (mimics clinical judgment differences)

---

## Multi-Drug Therapy Decisions

### When Multiple Drugs Are Used:
1. **Sequential Addition**: Model starts one drug, may add others on subsequent days using same selection logic
2. **Severe Infections**: Up to 3 drugs allowed simultaneously
3. **MDR TB Treatment**: 
   - Sequential drug selection (TB drugs chosen individually)
   - Automatic synergy when ≥2 TB drugs active: 2.5x effectiveness multiplier
   - Background effectiveness bonus: 0.8 (representing unmodeled TB drugs: bedaquiline, pretomanid, delamanid, cycloserine)
   - Reflects biological requirement for multi-drug TB therapy
4. **Combination Therapy Scenarios**:
   - Immunocompromised patients (8x higher prophylactic antibiotic rates)
   - Severe sepsis
   - Known resistant pathogens (resistance testing guided)
   - ICU settings (broad-spectrum preference)

### Combination Restrictions:
- Maximum 3 drugs to prevent unrealistic polypharmacy
- Must have clinical justification (severe infection, etc.)
- Regional availability must support combination therapy

---

## Treatment Failure and Drug Switching

### Treatment Failure Assessment:
- **Timing**: Evaluated on day 4 of treatment (configurable)
- **Criteria**: Bacteria level still ≥50% of initial level
- **Action**: Switch to alternative drug if failure detected

### Drug Switching Logic:
1. **Identify failure**: Compare current vs initial bacteria levels
2. **Find alternatives**: Score available drugs (excluding recently failed ones)
3. **Memory effect**: Avoid drugs that failed recently (30-day window)
4. **New selection**: Use same scoring system as initial prescription

---

## Special Populations and Scenarios

### Immunocompromised Patients:
- 8x higher probability of receiving prophylactic antibiotics
- Broader spectrum agents preferred
- Earlier initiation thresholds

### ICU Patients:
- Higher initiation rates
- Preference for broad-spectrum agents
- More aggressive combination therapy

### Pediatric Patients:
- Age-appropriate drug selection
- Dosing considerations built into scoring

### Historical Accuracy:
- Drug availability changes over time
- Prescribing patterns evolve with evidence
- Resistance patterns influence choices realistically

---

## Appendix: Comprehensive Pathogen-Specific Drug Selection Guidelines

### Escherichia coli
- **Ceftriaxone**: Strongly preferred for susceptible isolates (20x multiplier)
- **Meropenem/Imipenem**: Reserved for resistant cases, especially ESBL producers (25x multiplier for years 1-3)
- **Ciprofloxacin**: First-line oral agent (35x multiplier, but decreased with resistance)
- **Nitrofurantoin**: UTI-specific agent (30x multiplier)
- **Trimethoprim-sulfamethoxazole**: Alternative oral option (25x multiplier)
- **Piperacillin-Tazobactam**: Broad-spectrum option (18x multiplier)
- **Ampicillin**: Time-dependent use - 25x before 1970s resistance, 3x after
- **Ertapenem**: Carbapenem-sparing alternative (moderate preference)

### Klebsiella pneumoniae
- **Ceftriaxone**: Time-dependent preference - 25x before ESBL era, 8x after resistance emergence
- **Meropenem/Imipenem**: Critical for ESBL and carbapenem-resistant strains (30x in ESBL era, 3x before)
- **Ciprofloxacin**: Moderately preferred option (15x multiplier)
- **Piperacillin-Tazobactam**: Broad-spectrum alternative (18x multiplier)
- **Ertapenem**: Carbapenem option with narrower spectrum

### Enterococcus species
- **E. faecalis**: Ampicillin first-line (40x multiplier), vancomycin (8-30x based on VRE prevalence)
- **E. faecium**: Higher resistance profile - vancomycin (15-35x), linezolid for VRE era
- **Linezolid**: Alternative for VRE, especially in serious infections (18x multiplier)
- **Teicoplanin**: Alternative glycopeptide option
- **Vancomycin**: Primary treatment for VRE and severe infections (20x multiplier)

### Acinetobacter baumannii
- **Carbapenems**: Time-dependent use - 40x before resistance, 15x after emergence
- **Colistin**: Last resort for MDR cases (35x in resistance era)
- **Ampicillin-Sulbactam**: Intrinsic anti-Acinetobacter activity (25x multiplier)
- Limited treatment options reflecting clinical reality of MDR Acinetobacter

### Pseudomonas aeruginosa
- **Meropenem/Imipenem**: Anti-pseudomonal carbapenems for serious infections (25x multiplier for years 1-3)
- **Piperacillin-Tazobactam**: Broad-spectrum anti-pseudomonal option (18x multiplier)
- **Ceftazidime**: Anti-pseudomonal cephalosporin (specific preference)
- **Ciprofloxacin**: Oral option with anti-pseudomonal activity (15x multiplier)
- **Colistin**: Last resort for MDR Pseudomonas

### Staphylococcus aureus
- **Vancomycin**: Gold standard for MRSA infections (20x multiplier)
- **Linezolid**: Alternative for MRSA with good tissue penetration (18x multiplier)
- **Clindamycin**: Option for susceptible strains and soft tissue infections
- **Teicoplanin**: Alternative glycopeptide for MRSA
- **Methicillin/Flucloxacillin**: First-line for MSSA (when susceptible)

### Streptococcus pneumoniae
Infection site-specific preferences (probability distribution):
- **Amoxicillin**: 70% probability - Primary oral choice for pneumonia
- **Cephalexin**: 15% - Alternative oral cephalosporin  
- **Piperacillin**: 8% - IV option for severe infections
- **Sulfanilamide**: 4% - Historical agent, limited modern use
- **Penicillin G**: 2% - IV penicillin for severe infections
- **Ceftriaxone**: 1% - IV cephalosporin for severe cases

### Streptococcus pyogenes (Group A Strep)
Infection site-specific preferences:
- **Penicillin G**: 50% - Remains first-line treatment (still universally susceptible)
- **Amoxicillin**: 25% - Oral alternative
- **Piperacillin**: 15% - IV option for severe infections
- **Ceftazidime**: 5% - Broad-spectrum alternative
- **Ticarcillin**: 3% - Extended-spectrum penicillin
- **Sulfanilamide**: 2% - Limited historical use

### Enterobacter species
- **Meropenem/Imipenem**: Preferred for serious infections due to AmpC resistance risk (25x multiplier)
- **Ciprofloxacin**: Alternative fluoroquinolone (15x multiplier)
- **Piperacillin-Tazobactam**: Beta-lactamase inhibitor combination (18x multiplier)
- Systematic avoidance of cephalosporins due to inducible AmpC production

### Citrobacter species
- **Meropenem/Imipenem**: First-line for serious infections (25x multiplier for years 1-3)
- **Ciprofloxacin**: Alternative quinolone option (15x multiplier)
- **Piperacillin-Tazobactam**: Broad-spectrum choice (18x multiplier)
- Similar resistance patterns and treatment approach to Enterobacter

### Serratia marcescens
- **Meropenem/Imipenem**: Preferred carbapenems (25x multiplier for years 1-3)
- **Ciprofloxacin**: Quinolone alternative (15x multiplier)
- **Piperacillin-Tazobactam**: Beta-lactam option (18x multiplier)
- Intrinsically resistant to ampicillin and first-generation cephalosporins

### Morganella morganii
- **Meropenem/Imipenem**: Primary carbapenem choice (25x multiplier for years 1-3)
- **Ciprofloxacin**: Fluoroquinolone option (15x multiplier)
- **Piperacillin-Tazobactam**: Alternative beta-lactam (18x multiplier)
- Natural resistance to multiple beta-lactams requires careful selection

### Additional Gram-positive Pathogens
**Streptococcus agalactiae (Group B Strep)**: Similar patterns to S. pyogenes
**Coagulase-negative Staphylococci**: Often treated similar to S. aureus with MRSA coverage