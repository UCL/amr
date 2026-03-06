# Antimicrobial Resistance (AMR) Individual-Based Model: A Clinical Overview

## 1. Introduction
This document outlines the core architecture and biological assumptions of our Antimicrobial Resistance (AMR) simulation model. It is completely abstracted from the underlying code, designed specifically for infectious disease physicians, clinical microbiologists, and epidemiological researchers to review the clinical logic and mechanics driving the simulation.

The platform simulates the transmission, resistance dynamics, and clinical trajectories of bacterial infections across a synthetic global population from 1930 to 2035. It dynamically models **42 distinct bacterial pathogens**, **58 antibiotics (across 18 drug classes)**, and **35 unique resistance mechanisms**, observing how resistance emerges from both *de novo* mutation and horizontal gene transfer (HGT), and how it is subsequently propagated across the community.

---

## 2. The "Virtual Patient" (Host Characteristics)
The model tracks each individual patient day-by-day. Every virtual patient has fundamental characteristics that dictate their baseline risk, care-seeking behavior, and clinical outcomes.

### Demographics and Age
Patients age dynamically, moving through distinct clinical risk epochs (Neonatal, Pediatric, Young Adult, and Elderly). Age significantly moderates:
* **Infection site risk:** (e.g., STI peaks in young adults, respiratory infections exhibit a U-shaped risk curve peaking in infants and the elderly).
* **Sepsis vulnerability:** Neonates and the elderly face heavily weighted risks of clinical deterioration.

### Immunocompetence
Patients may suffer from temporary (e.g., chemotherapy, post-viral) or chronic (e.g., HIV, organ transplant) immunosuppression. Immunocompromised hosts experience:
* Higher baseline probabilities of initially acquiring pathogens.
* A drastically accelerated transition from localized infection to systemic sepsis.
* Lower rates of natural, non-pharmacologic immune clearance.

### Geography and Movement
The population mimics the global distribution over 6 major continental regions. Patients probabilistically travel between regions, acquiring the localized resistance profile of their destination, acting as vectors transferring resistant plasmids and strains back to their home communities.

---

## 3. Pathogen Acquisition and The Microbiome
The simulation strictly delineates between asymptomatic colonization (the microbiome) and clinically active disease.

### Bacterial Colonization vs. Active Infection
Patients acquire bacterial strains continually from the community. Usually, these organisms enter the microbiome as colonizing flora.
* **Autoinfection / Endogenous Infection:** A substantial portion of active infections (e.g., *E. coli* causing a UTI, *S. aureus* skin infections) erupt from the patient's existing colonized flora. 
* **Community Transmission:** New strains are transmitted horizontally from patient to patient.
* **Hospital-Acquired (Nosocomial):** Hospitalized patients face drastically elevated acquisition rates of opportunistic flora (e.g., *C. difficile*, *A. baumannii*, *P. aeruginosa*).

### Clinical Syndromes
When an infection becomes active, it manifests as one of 10 distinct clinical syndromes (e.g., UTI, Respiratory, Intra-abdominal, CNS, Bloodstream). 
* The syndrome dictates intrinsic bacterial growth rates (e.g., rapid doubling in bloodstream vs slower in bone/joint).
* The syndrome modifies the patient's likelihood of seeking treatment and the physician's empiric antibiotic selection.

---

## 4. Clinical Progression

### Symptom Onset, Infection Sites, and Sepsis
Unchecked bacterial proliferation rapidly pushes a patient from an asymptomatic to a symptomatic state. The anatomical site of the infection (the "syndrome") dynamically drives both the pathogen’s growth rate and the patient's likelihood of seeking medical care. For instance, a bloodstream infection features rapid bacterial doubling and almost guaranteed healthcare seeking, whereas localized skin infections grow more slowly and are often managed outpatient or ignored.

* **Sepsis Onset:** If untreated or inappropriately treated, high bacterial loads trigger systemic decompensation. The likelihood of a patient deteriorating into sepsis is calculated daily using a massive logistic cascade of risk factors starting from a very low base log-odds of **-14.0**. This risk aggressively stacks based on:
   * **Bacterial Load:** Every unit of bacterial load adds **+0.8** log-odds.
   * **Age Vulnerability:** Neonates face a massive foundational penalty (**+1.10** log-odds, roughly a 3x risk increase), and the elderly face **+0.69** (a 2x increase). Young adults possess the baseline (0.0).
   * **Pathogen Aggression:** Bugs differ wildly in their intrinsic capacity to cause sepsis. Uncomplicated pathogens like *C. trachomatis* or *H. pylori* are deeply penalized (**-19.0** and **-250.0** log-odds respectively, effectively preventing sepsis organically), whereas invasive organisms like *S. pneumoniae* or *N. meningitidis* sit at a much more dangerous **-9.0** base.
   * **Care Status Penalty:** Patients lingering outside of medical care are heavily penalized (**+1.0** log-odds), mimicking the real-world deterioration seen in delayed presentations.

### Hospitalization
Hospital admission is a probabilistic event driven primarily by symptom severity and sepsis onset, mitigated by regional healthcare capacity. The model uses "log-odds" to calculate the daily probability of a patient being admitted to the hospital. A base log-odds of **-10.4** (~0.003% daily chance) is modified heavily by the patient's state:

* **Symptomatic Infection:** Adds **+2.5** to log-odds (~12x increase in admission likelihood).
* **Sepsis Onset:** Adds **+4.4** to log-odds (an extreme ~80x driver for immediate hospitalization).
* **Age Penalty:** Adds **+0.02** for every year of the patient's life (meaning older adults are steadily more likely to be admitted).
* **Regional Modifiers:** North America and Europe apply a bonus (**+0.5, +0.6**) simulating robust healthcare access, while Africa applies a penalty (**-0.5**) simulating restricted access.

Once admitted, the patient is exposed to the local nosocomial environment, severely elevating their risk for multi-drug resistant (MDR) superinfections (e.g., *C. difficile*, *A. baumannii*). Hospitalized patients are also tested much more frequently.

---

## 5. Antimicrobial Stewardship and Therapeutics

### Diagnostic Testing
The model simulates both historical timelines (testing becomes globally available roughly post-WWII) and modern constraints. The base probability of ordering a diagnostic culture is roughly **15%** per day, but this is multiplied by **8.0x** if the patient is hospitalized, and **4.0x** if they are septic.

* **Antibiotic Susceptibility Testing (AST)** results have a built-in delay (typically **3 days**). 
* This delay forces the simulation into a strictly realistic two-stage prescribing model: **Empiric Therapy** followed by **Targeted Therapy**.

### 5.1 The Empiric Therapy Phase (Blind Prescribing)
Before culture results return, physicians must choose a drug blindly. The model explicitly scores every available antibiotic on a real-time point system. The highest scoring drug is prescribed:

1. **Syndrome-Specific Templates:** Each of the 10 anatomical syndromes has a unique hierarchical "guideline" of preferred antibiotics. For example:
   * **UTI (Syndrome 1):** Nitrofurantoin (Score: 15.0), TMP-SMX (14.0), Ciprofloxacin (12.0). Meropenem is deliberately scored terribly low (3.0) to prevent its use as an outpatient empiric UTI drug.
   * **Bloodstream (Syndrome 4):** Piperacillin-Tazobactam (18.0), Meropenem (14.0), Vancomycin (13.0). Mild oral agents score very poorly (e.g., Azithromycin 4.0).
   * **Intra-abdominal (Syndrome 5):** Piperacillin-Tazobactam (14.0), Meropenem (13.0), Metronidazole (11.0).

2. **Broad-Spectrum Bias:** During empiric therapy, the model assumes clinicians are nervous and intentionally applies a **+0.85 broad-spectrum bonus** to agents that cover a wide net of organisms. This artificially elevates drugs like Piperacillin-Tazobactam over narrower agents during the blind phase.

3. **Regional Resistance Penalties (Antibiograms):** Clinician agents review regional surveillance data. If widespread resistance to a drug crosses specific thresholds, the physician explicitly abandons it empirically:
   * Moderate Resistance (>10%): The drug's score is multiplied by **0.80** (a 20% penalty).
   * High Resistance (>45%): The score is halved (**0.50 multiplier**).
   * Very High Resistance (>60%): The score is decimated (**0.30 multiplier**), effectively removing the drug from the empiric formulary (e.g., rendering fluoroquinolones unusable for empiric GN infections).

4. **Tissue Penetration (Pharmacokinetics):** A drug's theoretical efficacy is heavily modified by anatomical penetration (0.0 = no penetration, 1.0 = perfect penetration). For example:
   * **Aminoglycosides (AGs)** penetrate CNS (Syndrome 6) at **0.05**, meaning even if the pathogen is fully susceptible on the AST panel, the clinical agent recognizes it will fail the patient in a meningitis scenario.
   * **Daptomycin** penetrates the lung (Syndrome 3) at **0.0**, rendering it useless for pneumonia.
   * **Fluoroquinolones (FQs)** penetrate the Urinary Tract at **1.0**, driving rapid clearance if the bug is susceptible.

### 5.2 The Targeted Therapy Phase (De-escalation)
Once the AST culture sensitivities return, the clinician agent perfectly understands the pathogen's resistance profile and reorganizes its therapy:
1. **Microbiologic Match:** Drugs demonstrating *in vitro* resistance immediately face an **ineffective drug penalty (0.001 multiplier)**, essentially hard-banning them from selection.
2. **Narrow-Spectrum Bonus:** In a direct reversal of the empiric phase, the model now demands good stewardship. A massive **+5.0 narrow-spectrum bonus** is applied to targeted agents, alongside a **-0.10 penalty** for remaining on a broad-spectrum drug unnecessarily. Given two active drugs, the clinician agent will confidently de-escalate (e.g., stepping down from empiric Cefepime to targeted Ampicillin). 

### 5.3 Microbiome Collateral Damage
Every antibiotic administered carries a "spectrum breadth" rating (e.g., Penicillin G = **2.0**, Azithromycin = **4.0**, Meropenem = **5.0**). While a broad-spectrum drug (rating 5.0) might consistently cure the active infection, it simultaneously inflicts massive collateral damage on the patient’s asymptomatic microbiome. This ecological wipeout vacates host niches, allowing asymptomatic resistant clones to proliferate freely without competition—readying the patient to become a silent vector for MDR transmission once discharged.

---

## 6. Resistance Dynamics (The Core AMR Engine)
The simulation utilizes highly rigorous biological models for how bacteria evolve and maintain resistance under drug pressure, tracking exactly how genes translate into clinical failure.

### 6.1 The 35 Resistance Mechanisms
Rather than simply logging a bug as "resistant," the model explicitly tracks **35 distinct genetic or biochemical resistance mechanisms**. When a patient acquires resistance, they acquire a specific mechanism.
* **Beta-Lactamases:** The model tracks exact enzymatic variants: `ESBL CTX-M`, `ESBL TEM`, `ESBL SHV`, `AmpC CMY/DHA`.
* **Carbapenemases:** `KPC`, `NDM/VIM`, and `OXA-48`.
* **Target Site Alterations:** `PBP2a/MecA` for MRSA, `VanA/VanB` for VRE, `GyrA` and `GyrA+ParC` double mutations for fluoroquinolones.
* **Efflux Pumps & Porins:** `AcrAB-TolC`, `MexXY-OprM`, `OmpK35/36`, and `OprD`.
* **Other Elements:** `MCR-1` (Colistin), `ErmB` (Macrolides), `Cfr` (Linezolid/Multi), `16S rRMT` (Aminoglycosides), etc.

### 6.2 The Resistance Vocabulary (Metrics)
When evaluating the patient or the population, the simulation relies on a specific set of clinical metrics to determine drug behavior:
* **Potency:** The *intrinsic* baseline activity of a drug against a perfect wild-type bug before any resistance is considered. E.g., Meropenem has a potency of 0.95 against *E. coli*, but exactly 0.0 against *E. faecalis* (since it is intrinsically resistant).
* **Microbiome Resistance (microbiome_r):** The percentage of the patient's asymptomatic carriage flora that possesses resistance genes.
* **Any Resistance (any_r):** A binary or percentage check: Does this patient carry *any* resistant clones of this species, either in active infection or silent microbiome carriage?
* **Majority Resistance (majority_r):** The public health metric. In the general population or the specific region, what percentage of *observed* infections express majority resistance to the drug? This is the "Antibiogram" metric that informs clinical empiric penalties (Section 5.1).
* **Activity Resistance (activity_r):** The ultimate clinical metric. This represents "Effective Resistance." By combining the mechanism's genetic strength with the drug's inherent potency and patient tissue penetration, `activity_r` dictates the likelihood that the drug will physically fail to clear the active infection.

### 6.3 Cross-Resistance Groupings
Bacteria do not view drugs independently; they block classes. The model strictly bundles antibiotics into biological "cross-resistance" groups for each species.
* Eradicating susceptible *E. coli* with Ciprofloxacin might accidentally trigger the selection of a `GyrA+ParC` mutant. Because of cross-groupings, that *E. coli* is now simultaneously counted as resistant to Levofloxacin, Moxifloxacin, and Ofloxacin, instantly destroying the entirety of the fluoroquinolone empiric templates.
* For *S. aureus*, acquiring `PBP2a/MecA` instantly provides cross-resistance to *all* penicillins and *all* cephalosporins (excluding novel anti-MRSA cephalosporins if enabled), perfectly mimicking MRSA clinical behavior.

### 6.4 Mutation vs. Horizontal Gene Transfer (HGT)
* **De Novo Emergency:** Spontaneous chromosomal mutation rates are biologically specific and heavily suppressed. Rare chromosomal mutations (like early `GyrA`) occur at ~1x10^-9, whereas plasmid-mediated events occur at ~1x10^-6. We enforce "Incidence Band Multipliers" (e.g., highly prevalent bugs like *E. coli* or *S. pneumoniae* have their mutation rates scaled down x0.1 to avoid mathematically overheating global resistance just by sheer volume).
* **HGT (The Primary Driver):** AMR spreads primarily through the exchange of plasmids and transposons. The model groups bugs into biological "HGT Plasmid Pools" (e.g., *Enteric Gram-Negatives* exchange easily at 1x10^-10; transferring across Gram boundaries is biologically blocked at `0.0`). This HGT exchange happens continually, unseen, in the **asymptomatic microbiome**.

### 6.5 Fitness Cost, Reversion, and Competitive Exclusion
A critical biological assumption is the *fitness cost* of carrying resistance genes. Resistant organisms expend continuous biological energy replicating plasmids or expressing efflux pumps, making them intrinsically less fit than wild-type strains in a drug-free environment.
* **In the Patient (Reversion Rate):** Without active antibiotic pressure killing off susceptible competitors, resistance plasmids face a daily probability of being lost (reversion). Faster-replicating plasmids decay quickly: high-level `VanA` decays at 0.002/day, whereas hardcoded chromosomal mutations like `GyrA` decay glacially slowly at 0.0001/day.
* **In the Community (Dilution Factor):** When a patient sheds bacteria into the community, resistant strains face severe competitive exclusion from the massive environmental pool of susceptible wild-type bacteria. The simulation artificially suppresses the transmission strength of resistant strains (a "dilution factor" of ~50%) to perfectly mirror how resistant mutants fail to organically out-compete wild-type bugs unless public antibiotics are heavily overutilized.

### 6.6 Selection Pressure
Antibiotic administration provides the overriding selective pressure required for resistant strains to flourish. By functionally eliminating the susceptible (fitter) competitors, resistance explodes within the treated host's microbiome (transitioning from a minority `microbiome_r` to a dominant state), perfectly readying the patient to become an active vector for MDR organisms.

## 7. Strategic Utility
By bridging individual-level clinical behavior with population-level evolutionary biology, the simulation allows public health doctors to test interventions such as:
1. Shortening antibiotic durations.
2. Increasing diagnostic testing speed/availability (reducing time spent on broad empiric therapy).
3. Implementing vaccination campaigns (removing vulnerable transmission hubs).
4. Controlling nosocomial outbreaks.