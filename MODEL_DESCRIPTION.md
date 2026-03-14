# AMR Simulation — Model Description



## Contents

1. [Overview](#1-overview)
2. [Population and Demographics](#2-population-and-demographics)
3. [Infection Acquisition](#3-infection-acquisition)
4. [Clinical Progression](#4-clinical-progression)
5. [Diagnostic Testing](#5-diagnostic-testing)
6. [Antibiotic Treatment](#6-antibiotic-treatment)
7. [Resistance Dynamics](#7-resistance-dynamics)
8. [Microbiome and Carriage](#8-microbiome-and-carriage)
9. [Horizontal Gene Transfer](#9-horizontal-gene-transfer)
10. [Mortality](#10-mortality)
11. [Policy Evaluation](#11-policy-evaluation)
12. [Limitations](#12-limitations)
- [Appendix A — Bacteria, Drugs, Mechanisms and Enums](#appendix-a--bacteria-drugs-mechanisms-and-enums)
- [Appendix B — Parameter Reference](#appendix-b--parameter-reference)
- [Appendix C — Output Specification](#appendix-c--output-specification)

---



## 1. Overview


### 1.1 What this model does

Antimicrobial resistance (AMR) — the ability of bacteria to survive antibiotic treatment — is one of the most serious threats to global health. Understanding how resistance emerges, spreads, and responds to policy changes requires a model that captures the interplay between antibiotic use, bacterial biology, and healthcare systems.

This model simulates the emergence and dynamics of AMR across a synthetic human population from **1930 to 2035**. The simulation starts in 1930 because that is before antibiotics were widely available; by beginning at that point, the model can reproduce the entire historical arc of antibiotic introduction, rising consumption, and the gradual accumulation of resistance that followed.

The model tracks **42 bacterial species**, **58 antibiotics** (grouped into **36 drug classes**), and **40 resistance mechanisms**. The population is distributed across **6 world regions** (North America, Europe, Asia, Oceania, South America, Africa), each with distinct epidemiological, travel, hospitalisation, and healthcare profiles.


### 1.2 How the model works — a brief guide for non-modellers

This is an **individual-based model** (sometimes called an agent-based model). Rather than using equations to describe an entire population at once, it creates a virtual population of individual people — typically 100,000 — and simulates what happens to each of them, day by day, over more than 100 years.

**Time steps.** The simulation advances in discrete daily steps. Each simulated day, every living person in the population is processed through a sequence of **21 mechanistic rules**. These rules govern the events that can happen to a person on any given day:

- Ageing and demographic changes (births, deaths from non-infectious causes)
- Acquiring a new bacterial infection (from the community, in hospital, or from bacteria they already carry)
- The infection getting better or worse, potentially progressing to sepsis
- Being tested — first to identify which bacterium is causing the infection, then to check which antibiotics it is resistant to
- Starting, continuing, or stopping antibiotic treatment
- The bacteria developing new resistance (either by spontaneous mutation or by acquiring resistance genes from other bacteria)
- Dying from the infection, from sepsis, or from drug side-effects

**Stochastic (random) processes.** The model does not say "this person *will* get an infection today." Instead, it calculates a *probability* of each event occurring, then uses a random number to decide whether it actually happens. This means that running the same model twice will produce slightly different results — just as two otherwise identical hospitals would see different patients on any given day. This randomness is a feature, not a bug: it lets us see the range of plausible outcomes, not just a single prediction.

**Log-odds — a note on the mathematics.** Many sections of this document describe probabilities using **log-odds** (also called logit values). This is a standard technique in medical statistics. If you are not familiar with it, here is a brief explanation:

- A probability of 50% corresponds to log-odds of **0**.
- Negative log-odds mean the event is unlikely (log-odds of −2 ≈ 12% probability; −4 ≈ 2%).
- Positive log-odds mean the event is likely (log-odds of +2 ≈ 88%; +4 ≈ 98%).
- The model adds together multiple log-odds terms (e.g., a baseline term, an age term, a severity term) and then converts the total into a probability. This is exactly how logistic regression works — the same technique used in many clinical risk scores.

For example, the model might calculate the daily probability of starting antibiotics as: *baseline log-odds (−6.5) + symptomatic infection (+6.5) + sepsis (+6.0) + immunodeficiency (+2.08) = ...*. The sum is then converted to a probability between 0 and 1. A more negative sum means "very unlikely today"; a more positive sum means "almost certain today."

**Calibration.** The model's parameters (the numbers that control how fast infections spread, how often drugs are prescribed, how quickly resistance emerges, etc.) are adjusted — *calibrated* — so that the model's outputs match real-world data. For example, the model is calibrated against observed antibiotic consumption rates, resistance prevalence reported by surveillance networks (such as ECDC and CDC), and infection incidence data. Sections 2–10 describe what the model does; Appendix B lists all the parameter values.


### 1.3 Scope and purpose

The model is specifically designed for reconstructing the historical emergence and growth of AMR over time by mechanistically linking antibiotic consumption, biological mutability, and transmission. It evaluates the potential impact of antibiotic stewardship policies by recreating empirical observations of resistance incidence and separating resistance acquisition across different care settings (e.g., community-acquired versus hospital-acquired).


### 1.4 How to read the rest of this document

The document is structured to follow the journey a person takes through the model:

| Section | What it covers | Clinical analogy |
|---------|---------------|-----------------|
| **2. Population** | Who the simulated people are — age, sex, region, immune status | The patient demographics of a hospital |
| **3. Infection Acquisition** | How people catch bacteria | Epidemiology — incidence, risk factors, hospital vs community |
| **4. Clinical Progression** | What happens once infected — symptoms, syndromes, sepsis | The natural history of infectious disease |
| **5. Diagnostic Testing** | When and how bacteria and resistance are identified | Microbiology lab workflow |
| **6. Antibiotic Treatment** | How drugs are started, chosen, dosed, and stopped | Empiric and targeted prescribing |
| **7. Resistance Dynamics** | How bacteria become resistant and how resistance spreads | The biology of AMR — mechanisms, gene transfer, selection pressure |
| **8. Microbiome & Carriage** | Asymptomatic bacterial colonisation | Carrier states (e.g., MRSA nasal carriage) |
| **9. Horizontal Gene Transfer** | Bacteria sharing resistance genes | Plasmid transfer between species |
| **10. Mortality** | How people die in the model | Case fatality rates, sepsis mortality |
| **11. Policy Evaluation** | Comparing "what if" scenarios | Stewardship interventions, counterfactuals |
| **12. Limitations** | What the model does not capture | Caveats for interpretation |
| **Appendices** | Reference tables of all bacteria, drugs, parameters, and outputs | Data dictionary |

Each section gives the *what* (what the model does), the *why* (what real-world phenomenon it is trying to capture), and the *how* (the specific rules and parameter values). Parameter tables are included for completeness; you do not need to memorise them to understand the model's logic.

---



## 2. Population and Demographics

This section describes the virtual people in the model — who they are, where they live, and the health states they can be in. These characteristics determine each person's risk of infection, their likelihood of receiving treatment, and their chance of dying. The model needs these details because AMR outcomes in the real world differ enormously by age, geography, immune status, and care setting.


### 2.1 Initialisation

The population is created at day 0 (representing the calendar year 1930). Each individual is assigned:

- **Age**: Drawn from a continuous demographic distribution that encodes both living individuals and future births. Negative age values at initialisation represent individuals who have not yet been born; they enter the simulation exactly when their age reaches zero. This is how the model handles births over the 105-year simulation period without needing a separate birth process.
- **Sex**: Male or female, assigned with equal probability.
- **Region**: Sampled from demographic weights reflecting the global population distribution.

The six regions and their approximate population shares determine the starting geographical distribution:

| Region | Population Share |
|--------|------------------|
| Asia | ~55% |
| Europe | ~15% |
| Africa | ~12% |
| North America | ~9% |
| South America | ~6% |
| Oceania | ~3% |

These regions matter because they differ in antibiotic availability (some drugs reach low-income settings decades later), hospital capacity, testing rates, and the prevalence of specific pathogens. A person's region shapes nearly every aspect of their simulated clinical journey.


### 2.2 Ageing and age categories

Each day, every individual's age increments by one day. The model groups people into age categories that determine their risk profiles. This mirrors real clinical practice — a doctor treats a neonate very differently from an elderly patient.

**General age categories** (used for most risk calculations):

| Age Category | Age Range | Clinical relevance |
|--------------|-----------|-------------------|
| Infant | 0–1 year | Immature immune system, high infection susceptibility |
| Preschool | 1–5 years | Frequent respiratory and enteric infections |
| School Age | 5–18 years | Generally lowest infection risk |
| Young Adult | 18–50 years | Reference group for most risk calculations |
| Middle Age | 50–70 years | Increasing comorbidities |
| Elderly | 70+ years | Immunosenescence, highest mortality risk |

**Sepsis/mortality age categories** (a separate, finer grouping):

Neonates (0–28 days) have dramatically different infection risks and case-fatality rates compared to older infants — for example, Group B *Streptococcus* sepsis in a 5-day-old neonate is a very different clinical entity from a respiratory infection in a 10-month-old. To avoid averaging over these differences, the model uses a separate age classification for sepsis onset and infection-related mortality:

| Category | Age Range |
|----------|-----------|
| Neonatal | 0–28 days |
| Paediatric | 28 days–18 years |
| Young Adult | 18–50 years |
| Elderly | 50+ years |


### 2.3 Immunodeficiency

In real life, patients with weakened immune systems — whether from HIV, chemotherapy, organ transplantation, or simply old age — are at much higher risk of infection, harder to treat, and more likely to die. The model captures this through two types of immunosuppression.

**Temporary immunosuppression** represents acute episodes such as a short course of steroids, a viral illness that transiently suppresses immunity, or post-surgical immunosuppression. People enter this state at a rate of `0.00005` per day and recover at `0.01` per day (average duration ~100 days).

**Chronic immunosuppression** represents long-term conditions like HIV/AIDS, solid organ transplant, or autoimmune disease requiring ongoing immunosuppressive therapy. It develops at `0.00006` per day and recovers much more slowly at `0.0012` per day.

At the start of the simulation (and for new births), some individuals begin with pre-existing chronic immunodeficiency. The probability increases with age to reflect the accumulating burden of disease over a lifetime:

| Age group | Probability | Examples |
|-----------|------------|---------|
| 0–1 year | 30% | Congenital immunodeficiencies (e.g., SCID, DiGeorge) |
| 1–18 years | 20% | Childhood leukaemia, congenital conditions |
| 18–65 years | 40% | HIV, transplant recipients, autoimmune disease on biologics |
| 65+ years | 60% | Immunosenescence, multiple comorbidities, haematological malignancies |

**How immunodeficiency affects the clinical journey:**

The table below summarises all the ways immunosuppression changes a person's trajectory through the model. Each effect has a real-world clinical rationale:

| Effect | Parameter | Value | What this means in practice |
|--------|-----------|-------|-----------------------------|
| More likely to receive empiric antibiotics | `antibiotic_initiation_log_odds_immunodeficiency` | +2.08 | ~8× higher odds of being started on antibiotics, reflecting the lower threshold for prescribing in immunocompromised patients |
| More diagnostic testing | `testing_immunosuppressed_multiplier` | ×2.5 | Clinicians investigate more aggressively in immunocompromised hosts |
| Higher sepsis risk | `log_odds_sepsis_onset_immunosuppressed` | +0.7 | ~2× higher daily risk of developing sepsis |
| Harder to recover from sepsis | `sepsis_log_odds_immunosuppressed` | −1.0 | ~2.7× lower odds of daily recovery, reflecting poor immune clearance |
| Higher mortality from sepsis | sepsis death log-odds | +1.5 | ~4.5× higher risk of dying during sepsis |
| Higher mortality from drug toxicity | toxicity death log-odds | +0.9 | ~2.5× higher risk — reflects drug interactions and organ dysfunction |
| Higher background mortality | `log_odds_mortality_immunosuppressed` | +0.916 | ~2.5× overall mortality uplift |


### 2.4 Hospitalisation

Hospital admission matters for AMR because hospitals are where the most resistant organisms are found, where the broadest-spectrum antibiotics are used, and where vulnerable patients are concentrated. The model captures this by simulating daily admission decisions, length of stay, and the elevated risks of hospital-acquired (nosocomial) infection.

**Who gets admitted?** Each day, the model calculates a probability of hospital admission for every person using a logistic model (see Section 1.2 for an explanation of log-odds). The key factors are:

| Factor | Log-odds contribution | What it means |
|--------|----------------------|---------------|
| Baseline (healthy person) | −10.4 | Very low daily risk (~0.003%) — most people are not admitted on any given day |
| Age | +0.02 per year | Older patients are progressively more likely to be admitted |
| Sepsis | +4.4 | Sepsis is a strong driver of admission (~80× multiplier) |
| Symptomatic infection (severity > 3.0) | +2.5 | Moderate-to-severe infections prompt admission (~12× multiplier) |
| Regional healthcare access | varies (see below) | Reflects real-world differences in hospital capacity |

**Length of stay:** Once admitted, patients are discharged at a rate of `0.28` per day (average stay ~3.6 days), with a hard maximum of 30 days. Patients with active sepsis cannot be discharged — they remain in hospital until sepsis resolves or they die.

**Regional healthcare access:**

Not everyone in the world has equal access to hospitals. The model uses regional modifiers that adjust the admission threshold:

| Region | Modifier | Interpretation |
|--------|----------|---------------|
| Europe | +0.6 | Highest access (universal healthcare systems) |
| North America | +0.5 | Good access |
| Oceania | +0.4 | Good access in developed areas |
| Asia | 0.0 | Reference baseline (mixed access) |
| South America | −0.2 | Variable access |
| Africa | −0.5 | Most limited hospital capacity |

Negative values mean patients are *less* likely to be admitted — not because they are less sick, but because hospital beds are less available. This matters for AMR because patients who cannot access hospital care may not receive appropriate antibiotics or diagnostics.

**Nosocomial (hospital-acquired) risks:**

Being in hospital dramatically changes a patient's infection risk profile. Hospital patients are exposed to multi-drug-resistant organisms on surfaces, devices, and other patients. The model captures this with pathogen-specific hospital acquisition modifiers:

| Pathogen | Hospital modifier | Approximate risk multiplier | Clinical context |
|----------|------------------|---------------------------|-----------------|
| *A. baumannii* | +3.4 | ~30× | Ventilator-associated pneumonia, ICU pathogen |
| *E. faecium* (VRE) | +3.3 | ~27× | Line infections, post-surgical |
| *P. aeruginosa* | +3.0 | ~20× | Burns, wounds, ventilators |
| *S. aureus* (MRSA) | +2.3 | ~10× | Surgical site, line-related infections |
| *K. pneumoniae* | +2.0 | ~7× | Carbapenem-resistant strains in ICU |
| Community pathogens (*C. trachomatis*, *T. pallidum*, *Campylobacter*) | −0.6 to −1.5 | Lower in hospital | Sexually transmitted or food-borne — acquired in the community, not in hospital |

Hospital patients also face higher baseline mortality (+0.262 log-odds, ~1.3×) and higher sepsis onset risk (+0.5 log-odds, ~1.6×), but they also have a higher probability of *recovering* from sepsis (+0.8 log-odds) because of access to intensive care.


### 2.5 Travel

International travel is a well-documented driver of AMR spread. Travellers who visit regions with high resistance prevalence can acquire resistant bacteria and bring them home — this is how, for example, ESBL-producing *E. coli* from South and South-East Asia has spread to European populations.

The model simulates this by giving each person a small daily probability of travelling to another region (`0.00005` per day, roughly one trip every 55 years per person). Travel frequency varies by region of origin, reflecting real-world patterns:

| Region | Travel multiplier | Rationale |
|--------|------------------|-----------|
| Europe | ×3.5 | High international travel rates |
| North America | ×3.0 | High travel, large business travel |
| Oceania | ×2.5 | Geographic distance drives air travel |
| Asia | ×1.5 | Rapidly growing travel volumes |
| South America | ×0.8 | Moderate travel rates |
| Africa | ×0.3 | Lowest international travel rates |

When a person travels, they are temporarily exposed to the infection risks and drug availability of the destination region. This can mean acquiring bacteria with resistance patterns typical of that region. Age-specific modifiers capture the higher risk of travel-related enteric diseases in younger adults — for example, young European adults travelling to endemic areas face elevated risk of *Salmonella enterica* serovar Typhi (+0.8 log-odds) and *Shigella* spp. (+0.5 log-odds), while *V. cholerae* risk is suppressed (−1.0) for these demographics unless visiting highly endemic zones.

---



## 3. Infection Acquisition

This section describes how people in the model catch bacterial infections. In the real world, a person can acquire bacteria from three main sources: the community (e.g., food, water, close contacts), the hospital environment (e.g., ventilators, catheters, other patients), or their own body (bacteria they are already carrying asymptomatically can flare into active infection). The model captures all three pathways.


### 3.1 Community acquisition

Each day, every person who does not already have an active infection has a chance of acquiring any of the 42 bacterial species. The model calculates a separate probability for each species using a logistic model (see Section 1.2) that combines several risk factors:

- **Baseline acquisition rate** for the specific bacterium — some bacteria (e.g., *E. coli*) cause infections far more frequently than others (e.g., *L. monocytogenes*)
- **Region** — infection rates vary by geography due to climate, sanitation, and population density
- **Age** — infants and the elderly are more susceptible to most infections; sexually transmitted infections peak in young adults
- **Immune status** — immunosuppressed individuals are at higher risk
- **Season** — respiratory pathogens (e.g., *S. pneumoniae*) follow a sinusoidal seasonal pattern, peaking in winter
- **Calendar era** — some infections have become more or less common over the decades
- **Existing population-level resistance** (`majority_r`) — this shapes how likely it is that a newly acquired bacterium is already resistant to common drugs (see Section 3.4)

| Variable pattern | What it controls |
|------------------|-----------------|
| `bacteria_{name}_acquisition_log_odds` | How common this bacterium is overall |
| `{region}_bacteria_{name}_acquisition_log_odds` | Regional differences for this bacterium |
| `bacteria_{name}_log_odds_{age_category}` | Age-specific risk for this bacterium |
| `{bacteria}_{region}_log_odds_{age_category}` | Interaction between bacterium, region, and age |


#### Age risk templates

Different bacteria infect different age groups. Rather than setting individual age parameters for all 42 species, the model assigns each bacterium a **risk template** — a pattern describing how infection risk varies across six age bands. The multipliers below are applied to the baseline acquisition rate:

| Template | Typical use | 0–1y | 1–5y | 5–18y | 18–50y | 50–70y | 70+y | Clinical rationale |
|----------|------------|------|------|-------|--------|--------|------|--------------------|
| `respiratory` | *S. pneumoniae*, *H. influenzae* | 3.0 | 1.8 | 0.8 | 1.0 | 1.3 | 2.5 | U-shaped: infants and elderly most vulnerable |
| `gastrointestinal` | *Salmonella*, *Shigella* | 2.5 | 2.0 | 1.2 | 1.0 | 1.1 | 1.8 | Young children and elderly via food/water |
| `urogenital` | *E. coli* (UTI) | 1.2 | 0.8 | 0.9 | 1.0 | 1.4 | 2.2 | Rises with age, especially in women |
| `skin_soft_tissue` | *S. aureus* | 1.5 | 1.3 | 1.1 | 1.0 | 1.2 | 1.8 | Moderate age variation |
| `bloodstream` | *P. aeruginosa* | 4.0 | 2.0 | 0.7 | 1.0 | 1.5 | 3.0 | Neonates and elderly at highest risk |
| `sexually_transmitted` | *N. gonorrhoeae*, *C. trachomatis* | 0.1 | 0.2 | 0.8 | 1.0 | 0.8 | 0.3 | Peaks in sexually active adults |
| `flat` | Default | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | Equal risk across all ages |

A multiplier of 3.0 for infants with the `respiratory` template means that an infant is three times as likely to acquire that bacterium compared to a young adult (the reference group at 1.0).


### 3.2 Hospital acquisition

Hospitalised patients are exposed to a different set of pathogens and at different rates than people in the community. Instead of using community acquisition rates, the model uses separate hospital-specific acquisition parameters (`{bacteria}_log_odds_hospital_acquired`) for each species.

This reflects the clinical reality that hospitals concentrate drug-resistant organisms: patients on ventilators are exposed to *Acinetobacter* and *Pseudomonas*, patients with central lines to *Staphylococcus* and *Enterococcus*, and patients on broad-spectrum antibiotics to *C. difficile*. The specific hospital modifiers for each pathogen are listed in Section 2.4.


### 3.3 Carrier-derived infection

People can carry bacteria in their gut, skin, or respiratory tract without being ill — this is called **asymptomatic carriage** (see Section 8 for details). Occasionally, these carried bacteria can cause an active infection in the same person. This is called **endogenous infection** and is extremely important for AMR because:

- The carried bacteria may already be resistant (having been selected by previous antibiotic courses)
- The person's resistance profile passes directly from carriage to infection

This pathway is governed by two parameters:

| Parameter | Value | What it means |
|-----------|-------|---------------|
| `carrier_resistance_inheritance_probability` | 0.50 | 50% chance that the new infection inherits the exact resistance profile of the carried strain |
| `infection_from_microbiome_dampening` | 0.70 | A dampening factor that prevents carriage from converting to clinical infection too frequently — not every colonised patient develops disease |


### 3.4 Resistance at acquisition

When a new infection is acquired from the community, the model needs to decide: is this bacterium resistant to any drugs, and if so, which ones?

Rather than rolling a separate dice for each drug independently (which would produce unrealistic resistance patterns), the model uses a three-step process that reflects how resistance actually exists in bacterial populations:

1. **Population-level resistance prevalence** (`majority_r`): The model continuously tracks a rolling average of resistance across all infected individuals for every bacterium–drug combination. This represents "how much resistance is circulating in the community right now." The averaging window is ~3 years (`majority_r_window_days` = 1,000 days), capturing the inertia of resistance in populations.

2. **Community dilution**: Clinical samples (from hospitals and sick patients) tend to over-represent resistant strains. To account for the fact that community bacteria are less resistant than those seen in clinics, the model applies a dilution factor (`community_resistance_dilution_factor` = 0.50). In other words, if 40% of clinical *E. coli* UTIs are resistant to ciprofloxacin, the model assumes ~20% of community *E. coli* are resistant.

3. **Correlated mechanism profiles**: Instead of assigning resistance genes one at a time (which would miss the real-world phenomenon of multi-drug resistance genes travelling together on the same plasmid), the model samples **complete resistance profiles** from other currently infected individuals via a `MechanismProfileCache`. This ensures that if a person acquires a resistant *E. coli*, its resistance pattern looks biologically realistic — for example, an ESBL-producing strain that is also fluoroquinolone-resistant, mirroring how these resistances co-occur on real plasmids.

---



## 4. Clinical Progression

Once a person has acquired a bacterial infection, the model simulates the clinical course: which body site is affected, how the infection grows, whether it progresses to sepsis, and whether the body can clear it without treatment. This section mirrors the natural history of infectious disease — the journey from first exposure to clinical outcome.


### 4.1 Syndrome assignment

When a person develops an active infection, the model assigns an **anatomical syndrome** — the body site where the infection is located. This is one of the most consequential decisions in the model, because the syndrome determines:

- **Which drugs a doctor would choose** (empiric prescribing guidelines differ by site — see Section 6.2)
- **How well drugs can reach the infection** (drug penetration varies enormously by tissue — see Section 6.4)
- **How fast the bacteria multiply** (some sites, like the bloodstream, support rapid replication)
- **How likely the patient is to develop sepsis or die** (bloodstream infections are far more dangerous than skin infections)

The 10 syndromes in the model correspond to the major infectious disease presentations a doctor encounters:

| Syndrome | Index | Examples in clinical practice |
|----------|-------|------------------------------|
| UTI | 1 | Cystitis, pyelonephritis — the most common bacterial infection |
| Skin/soft tissue | 2 | Cellulitis, wound infections, abscesses |
| Respiratory | 3 | Community-acquired and hospital-acquired pneumonia |
| Bloodstream | 4 | Bacteraemia, line-related infections |
| Intra-abdominal | 5 | Peritonitis, appendicitis, biliary sepsis |
| CNS | 6 | Meningitis, brain abscess — drugs must cross the blood-brain barrier |
| Gastrointestinal | 7 | Gastroenteritis, food poisoning |
| Genital/pelvic | 8 | Sexually transmitted infections, pelvic inflammatory disease |
| Bone/joint | 9 | Osteomyelitis, septic arthritis — slow to resolve, needs prolonged treatment |
| Other | 10 | Device-related infections, undifferentiated febrile illness |


#### How syndromes affect disease behaviour

Each syndrome modifies two key aspects of the infection:

- **Treatment initiation multiplier** — how urgently the patient seeks care. A patient with pneumonia (×10) presents to a doctor far more quickly than one with a mild UTI (×1).
- **Bacterial growth rate multiplier** — how fast the bacteria replicate at that body site. Bacteria in the bloodstream (×1.4) multiply faster than bacteria embedded in bone (×0.85).

| Syndrome | Treatment-seeking multiplier | Growth multiplier | Clinical rationale |
|----------|-----------------------------|--------------------|-------------------|
| UTI | ×1.0 | ×1.0 | Reference group |
| Skin | ×1.0 | ×1.1 | Slightly faster growth in necrotic tissue |
| Respiratory | ×10.0 | ×1.2 | Breathlessness and fever drive rapid presentation |
| Bloodstream | ×1.0 | ×1.4 | Nutrient-rich blood supports rapid replication |
| Intra-abdominal | ×1.0 | ×1.15 | Moderate growth rate |
| CNS | ×1.0 | ×1.3 | Rapid replication in cerebrospinal fluid |
| GI | ×8.0 | ×1.1 | Diarrhoea and vomiting drive rapid presentation |
| Genital | ×12.0 | ×0.9 | High care-seeking for STI symptoms; indolent course |
| Bone/joint | ×1.0 | ×0.85 | Slow, deep-seated infection |


### 4.2 Infection dynamics

The model does not simply label someone as "infected" or "not infected." Instead, it tracks a numerical **infection level** — an abstract measure of bacterial burden — that rises and falls over time. This is conceptually similar to a bacterial colony count increasing on serial blood cultures.

- **Starting level**: When a person first acquires an infection, the bacterial load is low (`initial_infection_level` = 0.01).
- **Growth**: Each day, the bacteria multiply. The growth rate depends on the specific bacterium, the syndrome site (see above), and whether antibiotics are active.
- **Symptom threshold**: When the infection level reaches `3.0` (`symptomatic_infection_level_threshold`), the person develops noticeable symptoms — fever, pain, cough, etc. — and begins seeking medical care. Below this threshold, they are infected but feel well enough that they do not present to a doctor.

This mechanism matters for AMR because there is a window between acquiring an infection and becoming symptomatic during which bacteria are replicating without antibiotic pressure — and during which resistance can emerge or be selected.


### 4.3 Sepsis

**Sepsis** is a life-threatening organ dysfunction caused by a dysregulated host response to infection. In clinical practice, it is a medical emergency with high mortality. In the model, sepsis is a distinct state that dramatically increases both the urgency of treatment and the risk of death.

Each day, the model calculates the probability of a person's infection progressing to sepsis using a logistic model that combines:

| Risk factor | Parameter | Value | What it means |
|-------------|-----------|-------|---------------|
| Bacterial load | `log_odds_sepsis_infection_level` | +0.9 per unit | The sicker the patient (higher bacterial burden), the more likely sepsis becomes — the single strongest driver |
| Duration of infection | `log_odds_sepsis_infection_duration` | +0.005 per day | Untreated infections gradually become more dangerous |
| Neonatal age | `sepsis_age_log_odds_neonatal` | +1.10 | Neonates are ~3× more likely to develop sepsis |
| Elderly age | `sepsis_age_log_odds_elderly` | +0.69 | Over-70s are ~2× more likely |

**Not all bacteria or body sites carry equal sepsis risk.** Clinically, a bloodstream infection with *N. meningitidis* is incomparably more dangerous than a *C. trachomatis* genital infection. The model captures this through:

- **Per-bacterium baseline**: Ranges from very low (*E. coli* UTI: −21.0, making sepsis extremely rare for routine UTIs) to high (*N. meningitidis*: −1.2, reflecting its aggressive clinical course)
- **Per-syndrome modifier**: Bloodstream (+1.5) and CNS (+1.2) infections are far more likely to cause sepsis; genitourinary (−2.0) and skin (−1.0) infections far less so

**Regional factors** also affect sepsis risk, reflecting differences in healthcare access and sanitation:
- Europe: −0.6 (best mitigation)
- North America, Oceania: −0.5
- Asia: −0.1 (reference)
- Africa: +0.1 (least mitigation — patients present later and with fewer resources)


### 4.4 Natural clearance

Not all infections need antibiotics. The body's immune system can clear many infections on its own — this is why mild gastroenteritis usually resolves without treatment. The model captures this with a daily probability of natural clearance:

- **Baseline**: 1% per day chance of the immune system clearing the infection without treatment (`default_microbiome_clearance_probability_per_day` = 0.01)
- **Duration penalty**: The longer an infection has been established, the harder it is to clear. Old infections become up to 7× harder to eliminate (`carriage_duration_log_odds_coefficient` = −0.01 per day, capping at −2.0 log-odds). This reflects the real phenomenon of bacteria forming biofilms and adapting to the host environment over time.
- **Antibiotic-assisted clearance**: When the patient is on an effective antibiotic (i.e., one that the bacteria are susceptible to), the clearance probability jumps to 80% per day (`microbiome_clearance_probability_on_drug_treatment` = 0.80) — reflecting the dramatic difference between treated and untreated infection outcomes.

---



## 5. Diagnostic Testing

Diagnostic testing is the bridge between empiric prescribing (guessing which antibiotic to use) and targeted prescribing (knowing which antibiotic will work). In real clinical practice, a doctor sends a sample (blood, urine, sputum) to the microbiology laboratory, which first identifies which bacterium is causing the infection (culture), and then tests which antibiotics it is susceptible or resistant to (antimicrobial susceptibility testing, or AST). This process takes days — and during that waiting period, the patient is treated empirically based on clinical judgement.

The model simulates this entire workflow: the decision to send a test, the delay in getting results, the possibility of laboratory errors, and the historical availability of testing technology.


### 5.1 Historical introduction

Modern diagnostic microbiology did not exist in 1930. The model introduces testing capabilities at historically appropriate time points:

| Technology | Available from | ~ Calendar year | Clinical context |
|------------|---------------|-----------------|-----------------|
| **Bacterial culture** | Day 5,478 | ~1945 | Basic culture techniques became routine in the mid-20th century |
| **Antimicrobial susceptibility testing (AST)** | Day 9,131 | ~1955 | Standardised AST methods (e.g., disc diffusion) followed about a decade later |

Before these dates, all prescribing in the model is entirely empiric — doctors have no laboratory information to guide drug choice. This accurately represents the early antibiotic era, when penicillin was prescribed without knowing the susceptibility of the infecting organism.


### 5.2 The testing process

Once testing is available and ordered, the model simulates a realistic laboratory workflow:

| Step | Parameter | Value | What happens |
|------|-----------|-------|-------------|
| **Lab turnaround time** | `test_delay_days` | 3 days | Results are not available until 3 days after the sample is sent — the patient is treated empirically during this time |
| **AST completion rate** | `prob_test_r_done` | 95% | If a culture grows a bacterium, there is a 95% chance AST is performed (occasionally omitted for low-priority isolates or technical reasons) |
| **Reporting error rate** | `test_r_error_probability` | 2% | AST results are wrong 2% of the time — the lab reports a resistant organism as susceptible or vice versa. This reflects real-world issues with breakpoint interpretation, contaminated samples, and technical failures |

The 3-day delay is clinically significant: a patient with sepsis will receive 3 days of empiric therapy before any lab results arrive. If the empiric choice was wrong (e.g., the bacterium was resistant), those 3 days of ineffective treatment allow the infection to progress.


### 5.3 Who gets tested?

Not every infected patient gets tested. In practice, a community GP managing a simple UTI may prescribe antibiotics without sending a urine culture, while a hospital patient with sepsis will have blood cultures drawn immediately. The model captures these differences:

| Factor | Parameter | Value | Clinical meaning |
|--------|-----------|-------|-----------------|
| **Baseline culture rate** | `bacterial_testing_base_rate_per_day` | 15% per day | A symptomatic outpatient has a 15% daily chance of having a culture sent |
| **AST reflex rate** | `resistance_testing_base_rate_per_day` | 95% per day | Once a culture is positive, AST is almost always performed |
| **Sepsis** | `testing_sepsis_multiplier` | ×4.0 | Septic patients are tested urgently |
| **Immunosuppressed** | `testing_immunosuppressed_multiplier` | ×2.5 | Clinicians investigate more aggressively |
| **Hospitalised (culture)** | `bacterial_testing_hospital_multiplier` | ×8.0 | Hospital patients have far greater access to microbiology labs |
| **Hospitalised (AST)** | `resistance_testing_hospital_multiplier` | ×5.0 | Hospitals perform AST more routinely |

**Regional differences:** Laboratory capacity varies dramatically around the world. Many hospitals in sub-Saharan Africa lack the microbiological infrastructure that is routine in European hospitals. The model captures this with regional testing multipliers:

| Region | Testing multiplier | Context |
|--------|-------------------|---------|
| Europe | ×1.2 | Highest testing density |
| North America | ×1.1 | High infrastructure |
| Oceania | ×0.8 | Good but geographically dispersed |
| Asia | ×0.7 | Highly variable by country |
| South America | ×0.6 | Variable access |
| Africa | ×0.3 | Very limited lab infrastructure in many settings |

These regional differences have direct consequences for AMR: in settings where testing is rare, patients are more likely to continue on ineffective empiric therapy, creating selection pressure for resistance without the feedback loop of culture results to guide narrower prescribing.

---



## 6. Antibiotic Treatment

This section covers the entire antibiotic prescribing process as the model simulates it — from the decision to start an antibiotic, through drug selection and dosing, to stopping the course. This is the heart of the AMR model: antibiotic use drives the selection pressure that causes resistance to emerge and spread.

The model aims to reproduce how antibiotics are actually prescribed in clinical practice — including imperfect decisions, regional variation in drug access, and the difference between empiric therapy (best-guess prescribing before lab results are available) and targeted therapy (prescribing guided by culture and susceptibility results).


### 6.1 Treatment initiation — deciding to start antibiotics

Each day, the model decides whether to start a new antibiotic course for each person, using a logistic model (see Section 1.2). The probability of starting antibiotics depends on the person's clinical state:

| Factor | Log-odds | Approximate effect | Clinical rationale |
|--------|----------|-------------------|-------------------|
| Baseline (no symptoms) | −5.5 | ~0.4% daily chance | Occasionally antibiotics are prescribed without clear indication — this captures "just in case" prescribing |
| Symptomatic infection | +6.0 | Jumps to ~62% | Once a patient has obvious symptoms (fever, pain, etc.), prescribing becomes likely |
| Sepsis | +6.0 | Near-certain | Sepsis is a medical emergency requiring immediate antibiotics |
| Immunodeficiency | +2.08 | ~8× more likely | Clinicians have a much lower threshold for prescribing in immunocompromised patients |
| No clinical indication | −1.05 | ~3× less likely | A protective factor: if investigation shows no active infection, prescribing is dampened |
| Lab-confirmed infection | +0.92 | ~2.5× more likely | Positive culture results prompt targeted therapy |
| Already on an antibiotic | +0.18 | ~1.2× more likely | Patients in the "pharmacy loop" may accumulate additional agents (combination therapy) |

**Worked example:** Consider a 65-year-old immunosuppressed patient in hospital with a symptomatic *E. coli* UTI, with no lab results yet. Their daily initiation log-odds would be roughly: −5.5 (baseline) + 6.0 (symptomatic) + 2.08 (immunodeficiency) = +2.58, which converts to about an 93% probability of starting antibiotics today.

**Regional variation in antibiotic access:**

Not everyone who needs antibiotics can get them. The model captures large global differences in antibiotic access:

| Region | Log-odds modifier | Effect on prescribing | Rationale |
|--------|------------------|----------------------|-----------|
| North America, Europe, Oceania | 0.0 | Reference | Good pharmaceutical access |
| Asia | −0.5 | ~38% reduction | Variable access across countries |
| South America | −0.8 | ~55% reduction | Limited access in some settings |
| Africa | −1.4 | ~75% reduction | Major access barriers in many countries |

These access barriers have a paradoxical effect on AMR: in settings where antibiotics are hard to obtain, people die of treatable infections, but the selection pressure for resistance is also lower. The model captures both sides of this equation.


### 6.2 Drug selection — choosing which antibiotic to use

Once the model decides to start an antibiotic, it must choose *which* antibiotic. This is one of the most clinically complex parts of the model, because the choice depends on what the doctor knows at the time.

**Two modes of prescribing:**

1. **Empiric therapy** — the doctor has no lab results yet and must choose a drug based on clinical judgement. "The patient has a UTI; guidelines say to use nitrofurantoin or trimethoprim." The model uses syndrome-specific scoring templates (see below) to replicate this guideline-based prescribing.

2. **Targeted therapy** — lab results have identified the bacterium and its susceptibility profile. The doctor can now choose a drug known to work. The model strongly rewards narrow-spectrum choices at this stage (×5.0 bonus for narrow-spectrum drugs) and penalises unnecessary broad-spectrum use (×0.1 penalty), reflecting the principle of antibiotic de-escalation.

**How drug scoring works:**

For each candidate drug, the model calculates a score based on several factors:

| Scoring factor | Empiric phase | Targeted phase | What it captures |
|---------------|---------------|----------------|-----------------|
| Syndrome-specific template score | Primary driver | Secondary | How well this drug matches guidelines for the infection site |
| Spectrum width | Slight bonus (×0.85) for broad-spectrum | Strong penalty (×0.1) for broad-spectrum | Empiric: cast a wide net. Targeted: use the narrowest effective drug |
| Known ineffectiveness | Near-zero score (×0.001) | Near-zero score (×0.001) | Never select a drug that is known to not work |
| Narrow-spectrum bonus | — | ×5.0 | Reward de-escalation to targeted therapy |

**Regional resistance surveillance:** If population-level resistance data shows that a drug class is failing frequently in the region, the model penalises empiric use of that drug — mimicking real-world guideline updates when local resistance rates exceed thresholds:

| Local resistance rate | Empiric score penalty | Clinical parallel |
|----------------------|----------------------|------------------|
| >60% resistant | ×0.3 | Drug dropped from guidelines (e.g., ciprofloxacin for *E. coli* UTI in South-East Asia) |
| >45% resistant | ×0.5 | Drug used cautiously, alternatives preferred |
| >10% resistant | ×0.8 | Drug still used but with awareness of resistance risk |


#### Treatment cessation — stopping antibiotics

Patients stop their antibiotic course based on several factors:

| Scenario | Daily stop probability | Typical course length | Real-world parallel |
|----------|----------------------|----------------------|-------------------|
| Default course | 0.45% per day | ~14 days | Standard course for many infections |
| No active infection found | 15% per day | ~3 days | Antibiotics stopped when investigation shows no infection |
| Cholera / *E. coli* GI | 2.5% per day | 3–5 days | Short courses per guidelines |
| *S. aureus* / *S. pneumoniae* | 1.5% per day | ~7 days | Standard courses |
| MDR-TB | 0.06% per day | 6–24 months | Prolonged anti-TB regimens |


#### Syndrome-specific empiric scoring templates

The tables below show which drugs score highest for empiric prescribing in each syndrome. Higher scores mean the drug is more likely to be selected. These templates are calibrated to match real-world prescribing guidelines — for example, nitrofurantoin and trimethoprim-sulfamethoxazole score highest for UTI, while piperacillin-tazobactam and meropenem score highest for bloodstream infections.

**Syndrome 1 — UTI** *(most common bacterial infection; oral drugs preferred)*

| Drug | Score |
|------|-------|
| nitrofurantoin | 15.0 |
| trim_sulf | 14.0 |
| ciprofloxacin | 12.0 |
| levofloxacin | 11.0 |
| cephalexin | 10.0 |
| amoxicillin_clavulanate | 9.0 |
| cefuroxime | 8.0 |
| ampicillin | 6.0 |
| amoxicillin | 6.0 |
| ceftriaxone | 5.0 |
| gentamicin | 4.0 |
| meropenem | 3.0 |

**Syndrome 2 — Skin/Soft Tissue** *(anti-staphylococcal and streptococcal coverage; penicillins and cephalosporins first-line)*

| Drug | Score |
|------|-------|
| penicillin_g | 13.0 |
| amoxicillin_clavulanate | 12.0 |
| cephalexin | 12.0 |
| clindamycin | 11.0 |
| doxycycline | 10.0 |
| trim_sulf | 9.0 |
| cefazolin | 8.0 |
| vancomycin | 7.0 |
| linezolid | 6.0 |
| ciprofloxacin | 5.0 |
| azithromycin | 5.0 |
| metronidazole | 4.0 |

**Syndrome 3 — Respiratory** *(community-acquired pneumonia guidelines: beta-lactam + atypical cover; co-amoxiclav scores highest)*

| Drug | Score |
|------|-------|
| amoxicillin_clavulanate | 20.0 |
| amoxicillin | 12.0 |
| azithromycin | 12.0 |
| doxycycline | 11.0 |
| levofloxacin | 10.0 |
| ceftriaxone | 9.0 |
| cefuroxime | 8.0 |
| trim_sulf | 7.0 |
| moxifloxacin | 6.0 |
| penicillin_g | 6.0 |
| clarithromycin | 5.0 |
| ampicillin | 4.0 |
| ciprofloxacin | 3.0 |
| erythromycin | 3.0 |

**Syndrome 4 — Bloodstream** *(medical emergency; broad-spectrum IV agents with good bactericidal activity)*

| Drug | Score |
|------|-------|
| piperacillin_tazobactam | 18.0 |
| meropenem | 14.0 |
| vancomycin | 13.0 |
| ceftriaxone | 12.0 |
| cefepime | 11.0 |
| ampicillin_sulbactam | 10.0 |
| ciprofloxacin | 9.0 |
| levofloxacin | 9.0 |
| amoxicillin_clavulanate | 8.0 |
| linezolid | 7.0 |
| metronidazole | 5.0 |
| trim_sulf | 5.0 |
| azithromycin | 4.0 |
| gentamicin | 1.0 |
| tobramycin | 1.0 |
| amikacin | 1.0 |

**Syndrome 5 — Intra-abdominal** *(must cover Gram-negatives and anaerobes; combination therapy common)*

| Drug | Score |
|------|-------|
| piperacillin_tazobactam | 14.0 |
| meropenem | 13.0 |
| ampicillin_sulbactam | 11.0 |
| metronidazole | 11.0 |
| ceftriaxone | 10.0 |
| ciprofloxacin | 9.0 |
| levofloxacin | 9.0 |
| cefepime | 8.0 |
| amoxicillin_clavulanate | 7.0 |
| clindamycin | 6.0 |
| vancomycin | 4.0 |
| azithromycin | 3.0 |

**Syndrome 6 — CNS** *(meningitis; only drugs that cross the blood-brain barrier are useful — see Section 6.4)*

| Drug | Score |
|------|-------|
| ceftriaxone | 15.0 |
| vancomycin | 14.0 |
| ampicillin | 12.0 |
| meropenem | 11.0 |
| penicillin_g | 10.0 |
| linezolid | 8.0 |
| chloramphenicol | 7.0 |
| rifampicin | 6.0 |
| ciprofloxacin | 5.0 |
| trim_sulf | 5.0 |
| metronidazole | 4.0 |
| doxycycline | 3.0 |

**Syndrome 7 — Gastrointestinal** *(fluoroquinolones and azithromycin for bacterial gastroenteritis; metronidazole for anaerobes)*

| Drug | Score |
|------|-------|
| ciprofloxacin | 12.0 |
| furazolidone | 11.0 |
| azithromycin | 11.0 |
| metronidazole | 10.0 |
| doxycycline | 9.0 |
| trim_sulf | 8.0 |
| ceftriaxone | 7.0 |
| amoxicillin_clavulanate | 7.0 |

**Syndrome 8 — Genital/Pelvic** *(STI guidelines: ceftriaxone + azithromycin for gonorrhoea; doxycycline for chlamydia)*

| Drug | Score |
|------|-------|
| azithromycin | 13.0 |
| ceftriaxone | 13.0 |
| doxycycline | 12.0 |
| penicillin_g | 12.0 |
| amoxicillin_clavulanate | 9.5 |
| amoxicillin | 9.0 |
| cefuroxime | 9.0 |
| clindamycin | 9.0 |
| ciprofloxacin | 7.0 |
| levofloxacin | 6.5 |
| trim_sulf | 5.0 |
| rifampicin | 4.0 |
| metronidazole | 2.5 |

**Syndrome 9 — Bone/Joint** *(prolonged courses required; good bone penetration essential — rifampicin, fluoroquinolones, linezolid)*

| Drug | Score |
|------|-------|
| penicillin_g | 14.0 |
| cefazolin | 13.0 |
| ampicillin | 12.0 |
| vancomycin | 12.0 |
| linezolid | 11.0 |
| cephalexin | 11.0 |
| ceftriaxone | 11.0 |
| tedizolid | 10.0 |
| dalbavancin | 10.0 |
| clindamycin | 10.0 |
| ciprofloxacin | 9.0 |
| levofloxacin | 9.0 |
| rifampicin | 9.0 |
| trim_sulf | 8.0 |
| meropenem | 7.0 |
| piperacillin_tazobactam | 6.5 |

**Syndrome 10 — Other/Device-Related** *(broad empiric cover when site is uncertain; even scoring reflects clinical uncertainty)*

| Drug | Score |
|------|-------|
| piperacillin_tazobactam | 8.0 |
| cefepime | 8.0 |
| ceftriaxone | 8.0 |
| meropenem | 8.0 |
| imipenem_c | 8.0 |
| vancomycin | 8.0 |
| linezolid | 7.0 |
| ciprofloxacin | 7.0 |
| azithromycin | 6.0 |



### 6.3 Drug pharmacokinetics

In the real world, antibiotics do not stay at constant levels in the body — they peak after administration and then decay as the body metabolises and excretes them. The model uses a simplified pharmacokinetic model where each drug has a **half-life** (the time for its level to halve) and a **starting level** at administration.

This matters for AMR because sub-therapeutic drug levels — where the drug concentration is too low to kill bacteria but high enough to exert selection pressure — are a key driver of resistance emergence (see "danger zone" in Section 7.3).

| Parameter | Default | What it represents |
|-----------|---------|-------------------|
| `drug_{name}_half_life_days` | Drug-specific | How quickly the drug is cleared from the body |
| `drug_{name}_initial_level` | 10.0 | Drug level immediately after dosing |
| `drug_{name}_double_dose_multiplier` | 2.0 | Level when a double dose is given |
| `drug_{name}_spectrum_breadth` | 3.0 | How broadly the drug disrupts the microbiome (higher = kills more bystander bacteria = more collateral damage) |


#### Selected drug half-lives

Half-lives vary enormously — from penicillin G (cleared within an hour, needing frequent dosing) to dalbavancin (which persists for two weeks, enabling single-dose therapy):

| Drug | Half-life (days) | Clinical note |
|------|-----------------|---------------|
| penicillin_g | 0.042 (~1 hour) | Very short — needs IV infusion or frequent dosing |
| ampicillin | 0.063 (~1.5 hours) | Short-acting penicillin |
| meropenem | 0.042 (~1 hour) | Short — given as IV infusion TDS |
| cefiderocol | 0.10 (~2.4 hours) | Short-acting novel siderophore cephalosporin |
| ciprofloxacin | 0.17 (~4 hours) | Moderate — allows twice-daily oral dosing |
| linezolid | 0.21 (~5 hours) | Moderate |
| vancomycin | 0.25 (~6 hours) | Requires therapeutic drug monitoring |
| sulfanilamide | 0.29 (~7 hours) | Historical agent |
| ceftriaxone | 0.33 (~8 hours) | Long enough for once-daily dosing |
| doxycycline | 0.75 (~18 hours) | Long — convenient once or twice-daily oral |
| azithromycin | 2.92 (~70 hours) | Very long tissue half-life — enables 3–5 day courses |
| dalbavancin | 14.0 (2 weeks) | Ultra-long — allows single-dose outpatient treatment |


#### Spectrum breadth — collateral damage to the microbiome

Broad-spectrum antibiotics kill not only the target pathogen but also many commensal ("friendly") bacteria in the gut, skin, and respiratory tract. This collateral damage creates ecological niches for resistant organisms to fill — one of the key mechanisms driving AMR. The model captures this via a `spectrum_breadth` parameter:

| Drug | Breadth | Meaning |
|------|---------|---------|
| penicillin_g | 2.0 (Narrow) | Minimal disruption to the microbiome |
| linezolid | 2.0 (Narrow) | Targets Gram-positives only |
| vancomycin | 2.5 (Narrow-medium) | Mainly Gram-positive spectrum |
| trim_sulf | 3.5 (Medium-broad) | Moderate disruption |
| azithromycin | 4.0 (Broad) | Significant microbiome disruption |
| ceftriaxone | 4.0 (Broad) | Major disruption; linked to *C. difficile* risk |
| ciprofloxacin | 4.5 (Very broad) | Extensive gut microbiome disruption |
| meropenem | 5.0 (Very broad) | Maximum disruption — the "sledgehammer" antibiotic |


### 6.4 Drug penetration by syndrome

A drug can only work if it reaches the infection site at adequate concentrations. This is particularly important for certain anatomical sites:

- **CNS (meningitis):** The blood-brain barrier blocks most antibiotics. Only a few drugs (ceftriaxone, metronidazole, chloramphenicol, linezolid) achieve therapeutic levels in cerebrospinal fluid.
- **Bone/joint:** Drugs must penetrate dense, poorly vascularised tissue. Rifampicin and fluoroquinolones penetrate well; aminoglycosides do not.
- **Bloodstream:** By definition, any IV drug achieves full levels here (penetration = 1.0 for all drugs).

Penetration values range from 0.0 (no drug reaches the site) to 1.0 (full systemic concentration available):

| Syndrome | Best penetration | Poorest penetration |
|----------|-----------------|---------------------|
| UTI (1) | FQ, TMP-SMX, nitrofurantoin, fosfomycin (1.0) | Macrolides (0.4), clindamycin (0.3), daptomycin (0.1) |
| Skin (2) | Daptomycin (0.95), FQ (0.9), oxazolidinones (0.9) | Nitrofurantoin (0.2) |
| Respiratory (3) | Macrolides (0.95), FQ (0.95), oxazolidinones (0.9) | Daptomycin (0.0), aminoglycosides (0.4) |
| Bloodstream (4) | All 1.0 (reference compartment) | — |
| Intra-abdominal (5) | Metronidazole (0.9), FQ (0.75), carbapenems (0.75) | Aminoglycosides (0.3) |
| CNS (6) | Metronidazole (0.80), oxazolidinones (0.70), chloramphenicol (0.70) | Aminoglycosides (0.05), colistin (0.05), daptomycin (0.05) |
| GI (7) | Fidaxomicin (1.0), metronidazole (0.95), oral vancomycin (0.90) | Glycopeptides IV (0.35) |
| Genital (8) | FQ (0.9), metronidazole (0.8), TMP-SMX (0.8) | Aminoglycosides (0.35) |
| Bone/joint (9) | Rifampicin (0.80), oxazolidinones (0.75), FQ (0.70) | Aminoglycosides (0.25), colistin (0.2) |

These penetration values directly affect treatment outcomes in the model: a drug with 0.05 penetration to the CNS will be nearly ineffective for meningitis even if the bacterium is fully susceptible.


### 6.5 Drug potency matrix

Not all antibiotics work against all bacteria. Penicillin G is highly effective against *Streptococcus pneumoniae* (potency 0.90) but has zero activity against *Pseudomonas aeruginosa* (intrinsically resistant). The model encodes this in a **potency matrix** — a 42×52 table (42 bacteria × 52 drug groups) where each cell represents the intrinsic activity of that drug against that bacterium when no acquired resistance is present.

Values range from 0.0 (no activity — the drug simply does not work against this organism) to 1.0 (maximum activity). These potency values are based on published MIC (minimum inhibitory concentration) data and clinical breakpoints.

Key examples:
- Meropenem vs *E. coli*: 0.95 (very high potency — a carbapenem is one of the most effective drugs against Gram-negatives)
- Vancomycin vs *E. coli*: 0.0 (vancomycin does not work against Gram-negative bacteria)
- Ceftriaxone vs *S. pneumoniae*: 0.90 (standard treatment for pneumococcal meningitis)



### 6.6 Drug availability by region and era

The model simulates antibiotics becoming available at their historical introduction dates. Before penicillin was introduced in 1942, there were no antibiotics in the model. Before ciprofloxacin was introduced in 1987, fluoroquinolones did not exist. This historical layering is essential for reproducing the sequential emergence of resistance over the 20th century.

**Regional availability:** Even after a drug is introduced globally, not all regions have equal access. Newer, more expensive drugs may be unavailable or rarely used in low-income settings:

| Region | Access pattern |
|--------|---------------|
| North America | Full access to all drugs |
| Europe | Full access to all drugs |
| Asia | Most drugs available; limited access to tedizolid, ceftaroline (30%) |
| Oceania | Good access; limited novel agents (50%) |
| South America | Limited newer drugs (tedizolid 10%, linezolid 50%, carbapenems 60–70%) |
| Africa | Basic antibiotics available (80–100%); ceftriaxone 60%; vancomycin 30%; carbapenems 10–20%; most novel drugs 0–10% |

This has major implications for AMR: in Africa, where carbapenems are rarely available, carbapenem resistance may emerge more slowly — but when it does arrive (via travel or HGT), there are no last-resort drugs available to treat it.



#### Drug introduction dates

The 58 antibiotics in the model span 82 years of pharmaceutical development:

| Drug | ~Year | Drug | ~Year |
|------|-------|------|-------|
| sulfanilamide | 1937 | ceftriaxone | 1984 |
| penicillin_g | 1942 | piperacillin_tazobactam | 1984 |
| tetracycline | 1948 | ceftazidime | 1985 |
| chloramphenicol | 1949 | imipenem_c | 1985 |
| colistin | 1952 | amoxicillin_clavulanate | 1985 |
| erythromycin | 1952 | aztreonam | 1986 |
| nitrofurantoin | 1953 | ciprofloxacin | 1987 |
| furazolidone | 1955 | teicoplanin | 1988 |
| vancomycin | 1958 | ampicillin_sulbactam | 1990 |
| fosfomycin | 1959 | clarithromycin | 1990 |
| metronidazole | 1960 | ofloxacin | 1990 |
| ampicillin | 1961 | azithromycin | 1991 |
| fusidic_a | 1962 | cefepime | 1996 |
| gentamicin | 1963 | meropenem | 1996 |
| rifampicin | 1966 | levofloxacin | 1996 |
| doxycycline | 1967 | moxifloxacin | 1999 |
| clindamycin | 1968 | linezolid | 2000 |
| trim_sulf | 1968 | ertapenem | 2001 |
| cephalexin | 1970 | daptomycin | 2005 |
| minocycline | 1971 | ceftazidime_avibactam | 2006 |
| amoxicillin | 1972 | tigecycline | 2007 |
| cefazolin | 1973 | ceftaroline | 2010 |
| tobramycin | 1975 | fidaxomicin | 2011 |
| amikacin | 1976 | tedizolid | 2014 |
| ticarcillin | 1977 | dalbavancin | 2014 |
| cefuroxime | 1978 | ceftolozane_tazobactam | 2014 |
| piperacillin | 1981 | meropenem_vaborbactam | 2018 |
| ticarcillin_clavulanate | 1990 | cefiderocol | 2019 |
| quinu_dalfo | 1999 | retapamulin | 2007 |

**Special case — Colistin:** Colistin was introduced in 1952 but withdrawn from routine use between ~1970 and ~1995 due to severe nephrotoxicity. It was then reintroduced as a last-resort agent for multi-drug-resistant Gram-negative infections. The model reflects this by dropping colistin availability to 5% during the withdrawal window.



### 6.7 Drug toxicity

Antibiotics are not without harm. Some drugs — particularly aminoglycosides (nephrotoxicity, ototoxicity) and colistin (nephrotoxicity) — carry significant toxicity risks. The model simulates drug toxicity as a **reservoir** that accumulates with continued use and decays when the drug is stopped.

Toxicity can cause two outcomes:

**1. Drug discontinuation (sub-lethal toxicity):** When toxicity accumulates, the treating clinician may stop the drug. This is the more common outcome — in real life, a rising creatinine on gentamicin prompts the team to switch to a less nephrotoxic alternative.

| Factor | Parameter | Value | What it means |
|--------|-----------|-------|---------------|
| Baseline discontinuation risk | `toxicity_discontinuation_base_log_odds` | -3.0 | Low baseline — drugs are not stopped without reason |
| Toxicity level | `toxicity_discontinuation_log_odds_per_reservoir_unit` | +1.5 | Higher toxicity → more likely to stop |
| Sepsis protection | `toxicity_discontinuation_log_odds_sepsis` | -1.5 | Clinicians tolerate more toxicity during sepsis because the alternative (no antibiotic) is worse |
| Recent toxicity avoidance | `toxicity_avoidance_penalty_multiplier` | x0.05 | A drug that caused toxicity recently is unlikely to be re-selected |
| Avoidance window | `toxicity_avoidance_window_days` | 14 days | How long the prescriber "remembers" to avoid that drug |

**2. Drug-related death (lethal toxicity):** Rarely, severe drug toxicity can be fatal — for example, acute kidney injury from colistin leading to multiorgan failure. The baseline risk is very low (log-odds -8.0), but it increases with accumulated drug level, age, and immunosuppression.

| Factor | Parameter | Value |
|--------|-----------|-------|
| Baseline | `toxicity_death_base_log_odds` | -8.0 (very rare) |
| Per unit of toxicity | `toxicity_death_log_odds_per_reservoir_unit` | +2.0 |
| Infants | `toxicity_death_log_odds_age_infant` | +0.6 |
| Elderly | `toxicity_death_log_odds_age_elderly` | +0.8 |
| Immunosuppressed | `toxicity_death_log_odds_immunosuppressed` | +0.9 |


### 6.8 Antibiotic infection prevention

Patients who are already receiving an effective antibiotic are partially protected against acquiring new infections — the drug in their system kills susceptible bacteria before they can establish. This mirrors the real-world concept of antibiotic prophylaxis (e.g., surgical prophylaxis with cefazolin prevents wound infections).

The model applies a 70% reduction in new infection risk for susceptible organisms when the patient is already on an active antibiotic (`antibiotic_infection_prevention_efficacy` = 0.7). This does *not* protect against resistant organisms — a crucial point, because it means patients on antibiotics are selectively more likely to acquire resistant infections relative to susceptible ones, creating further selection pressure for resistance.
---


## 7. Resistance Dynamics

This section describes the heart of the model — how bacteria become resistant to antibiotics. For a clinician, this section explains the mechanisms behind the resistance patterns you see in microbiology reports. For example, when your lab reports "ESBL-producing *E. coli*", the model tracks the specific enzyme (CTX-M, TEM, or SHV) that produces that phenotype, which drugs it affects, and how it spreads.

The model tracks resistance at the level of individual **mechanisms** — the specific biological tools bacteria use to evade antibiotics. This matters because the same phenotype (e.g., "carbapenem-resistant *K. pneumoniae*") can arise from very different mechanisms (KPC, NDM, OXA-48), each with different implications for treatment, spread, and even which novel drugs might still work.


### 7.1 Resistance mechanisms

The model explicitly tracks **40** distinct resistance mechanisms. Each mechanism represents a specific biological pathway: an enzyme that destroys the drug, a mutation that changes the drug's target, a pump that ejects the drug from the cell, or a barrier that prevents the drug entering.

The table below lists every mechanism, the drugs it affects, and which bacterial groups can acquire it. You do not need to memorise this table — it is a reference. The key insight is that each mechanism has a defined scope: ESBL enzymes (rows 1–3) destroy penicillins and cephalosporins but not carbapenems, while KPC and NDM (rows 6–7) destroy carbapenems as well.


  | Mechanism | Variable name | Description | Explicit Drugs Affected | Bacterial Classes Affected |
  |-----------|--------------|-------------|-------------------------|----------------------------|
  | ESBL CTX-M | `esbl_ctx_m` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | ESBL TEM | `esbl_tem` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | ESBL SHV | `esbl_shv` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AmpC CMY | `ampc_cmy` | Plasmid-mediated AmpC β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AmpC DHA | `ampc_dha` | Plasmid-mediated AmpC β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | KPC | `kpc` | *K. pneumoniae* carbapenemase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | NDM/VIM | `ndm_vim` | Metallo-β-lactamases | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `ceftazidime_avibactam`, `meropenem_vaborbactam`, `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | OXA-48 | `oxa_48` | Oxacillinase-type carbapenemase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftaroline`, `meropenem`, `imipenem_c`, `ertapenem`, `meropenem_vaborbactam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | OXA-Acinetob. | `oxa_acinetobacter` | OXA-23/40/58 carbapenemases (A. baumannii) | `meropenem`, `imipenem_c`, `ertapenem`, `ceftazidime`, `cefepime`, `ceftazidime_avibactam` | Nonfermenters |
  | blaZ | `blaz` | Staphylococcal penicillinase | `penicillin_g`, `ampicillin`, `amoxicillin` | Gram-Positives |
  | PBP2a/MecA | `pbp2a_meca` | PBP alteration (MRSA) | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `cefiderocol`, `ceftazidime_avibactam`, `meropenem_vaborbactam`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem` | Gram-Positives, Helicobacter |
  | VanA | `vana` | High-level vancomycin resistance | `vancomycin`, `teicoplanin`, `dalbavancin` | Gram-Positives, Helicobacter |
  | VanB | `vanb` | Variable-level vancomycin resistance | `vancomycin` | Gram-Positives, Helicobacter |
  | GyrA (pri.) | `gyra_primary` | DNA gyrase mutation (step 1) | `ciprofloxacin`, `ofloxacin` | All |
  | GyrA + ParC | `gyra_parc` | Additional topoisomerase mutation | `ciprofloxacin`, `ofloxacin`, `levofloxacin`, `moxifloxacin` | All |
  | Qnr | `qnr` | Quinolone resistance protein | `ciprofloxacin`, `ofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | 16S rRMT | `16s_rrmt` | 16S rRNA methyltransferase | `gentamicin`, `tobramycin`, `amikacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AAC/APH/ANT | `aac_aph` | Aminoglycoside-modifying enzymes | `gentamicin`, `tobramycin`, `amikacin`, `streptomycin`, `neomycin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Gram-Positives |
  | ErmB | `ermb` | Erythromycin ribosome methylase | `erythromycin`, `azithromycin`, `clarithromycin`, `clindamycin`, `quinu_dalfo` | Gram-Positives, Anaerobes, Fastidious, Helicobacter |
  | 23S rRNA | `23s_rrna` | 23S rRNA point mutation | `erythromycin`, `azithromycin`, `clarithromycin` | Helicobacter, Enteric Pathogens, Fastidious, Gram-Positives |
  | Cfr | `cfr` | 23S rRNA methyltransferase | `linezolid`, `tedizolid`, `chloramphenicol`, `clindamycin`, `retapamulin` | Gram-Positives, Anaerobes, Fastidious, Helicobacter |
  | CAT | `cat` | Chloramphenicol acetyltransferase | `chloramphenicol` | All |
  | MCR-1 | `mcr_1` | Mobilised colistin resistance | `colistin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AcrAB-TolC | `acrab_tolc` | Gram-negative efflux pump | `tetracycline`, `doxycycline`, `minocycline`, `tigecycline`, `chloramphenicol`, `ciprofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | MexXY-OprM | `mexxy_oprm` | Pseudomonas-specific efflux pump | `tetracycline`, `doxycycline`, `minocycline`, `gentamicin`, `tobramycin`, `amikacin`, `chloramphenicol`, `ciprofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | Global eff. | `global_efflux` | Non-specific efflux upregulation | `tetracycline`, `doxycycline`, `minocycline`, `tigecycline`, `chloramphenicol`, `ciprofloxacin` | All |
  | TetA/B/C | `tet_abc` | Gram-negative tetracycline efflux | `tetracycline`, `doxycycline` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious |
  | TetM/TetO | `tetm` | Ribosomal protection | `tetracycline`, `doxycycline`, `minocycline` | All |
  | OmpK35/36 | `ompk35_36` | Outer membrane porin loss (Klebsiella) | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `ceftaroline`, `cefiderocol`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem`, `ciprofloxacin`, `levofloxacin`, `moxifloxacin`, `ofloxacin`, `gentamicin`, `tobramycin`, `amikacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | OprD | `oprd` | Outer membrane porin loss (Pseudomonas) | `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | Global por. | `global_porin_loss` | Non-specific porin downregulation | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `ceftaroline`, `cefiderocol`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem`, `ciprofloxacin`, `levofloxacin`, `moxifloxacin`, `ofloxacin`, `gentamicin`, `tobramycin`, `amikacin` | All |
  | Folate path | `folate_pathway` | Altered dihydrofolate reductase | `sulfanilamide`, `trim_sulf` | All |
  | Nitroreduct | `nitroreductase` | Nitroreductase loss | `metronidazole`, `nitrofurantoin`, `furazolidone` | Enterobacterales, Enteric Pathogens, Anaerobes, Fastidious, Helicobacter |
  | FosA | `fosa` | Fosfomycin-modifying enzyme | `fosfomycin` | Enterobacterales, Nonfermenters, Enteric Pathogens |
  | MprF | `mprf` | Membrane charge modification | `daptomycin` | Gram-Positives |
  | RpoB | `rpob` | RNA polymerase mutation | `fidaxomicin` | All |
  | FusB | `fusb` | Fusidic acid resistance determinant | `fusidic_a` | Gram-Positives |
  | Unknown 1–3 | `as_yet_unknown_{1..3}` | Placeholder mechanisms | Evaluates `true` dynamically for all applied overrides | All (Calibration Placeholders) |



### 7.2 Mechanism–drug-class enhancement multipliers

When a bacterium possesses a resistance mechanism, it does not simply become immune to every drug. Instead, each mechanism **reduces** drug efficacy by a specific amount. The "enhancement multiplier" (0.0–1.0) represents **how much** of a drug's effectiveness is knocked out:

- **0.0** = the mechanism has no effect on this drug (e.g., a tetracycline efflux pump does nothing against meropenem)
- **0.95** = the mechanism eliminates 95% of the drug's activity (e.g., NDM metallo-β-lactamase virtually destroys carbapenem efficacy)
- **1.0** = complete resistance (the drug is useless)

In clinical terms, an enhancement multiplier of 0.85 for ESBL CTX-M against cephalosporins means: if a patient has an ESBL-producing *E. coli* UTI turned treated with ceftriaxone, the drug retains only 15% of its normal killing power — enough to provide some marginal activity but not enough to reliably cure the infection.

There are 40 mechanisms × 36 drug classes = 1,440 individual values. The table below shows the **global default** multiplier for each mechanism (used when a specific per-class value has not been configured):

| Mechanism | Multiplier | Clinical interpretation |
|-----------|-----------|----------------------|
| NDM/VIM | 0.95 | Near-complete resistance — these metallo-β-lactamases destroy almost all β-lactams |
| VanA | 0.95 | Near-complete vancomycin resistance |
| KPC | 0.90 | Very high — KPC carbapenemases severely compromise carbapenems |
| PBP2a/MecA | 0.90 | Very high — defines MRSA; eliminates nearly all β-lactam activity |
| ESBL CTX-M | 0.85 | High — but β-lactamase inhibitor combinations retain partial activity |
| VanB | 0.85 | High vancomycin resistance (but teicoplanin may still work) |
| GyrA + ParC | 0.85 | High-level fluoroquinolone resistance (double mutation) |
| 16S rRMT | 0.85 | High-level aminoglycoside resistance |
| ESBL TEM | 0.80 | High |
| OXA-48 | 0.80 | High — but with variable carbapenem MICs |
| ErmB | 0.80 | MLS_B resistance (macrolides, lincosamides) |
| RpoB | 0.80 | Rifampicin resistance |
| ESBL SHV | 0.75 | High |
| Cfr | 0.75 | Cross-resistance to oxazolidinones and phenicols |
| AmpC CMY/DHA | 0.70 | Moderate-high — overcomes β-lactamase inhibitors too |
| CAT | 0.70 | Chloramphenicol resistance |
| GyrA primary | 0.70 | First-step fluoroquinolone resistance (partial) |
| Folate pathway | 0.70 | Trimethoprim-sulfamethoxazole resistance |
| FusB | 0.70 | Fusidic acid resistance |
| FosA | 0.65 | Fosfomycin resistance |
| MCR-1 | 0.60 | Colistin resistance — critically important as colistin is the last resort |
| Nitroreductase | 0.60 | Nitrofurantoin resistance |
| OprD | 0.55 | Porin loss — carbapenem resistance (mainly in *Pseudomonas*) |
| MprF | 0.55 | Daptomycin resistance |
| OmpK35/36 | 0.50 | Porin loss — broad resistance in Enterobacterales |
| Qnr | 0.50 | Low-level quinolone resistance (facilitates further mutation) |
| Global porin loss | 0.45 | Broad, non-specific resistance via reduced permeability |
| MexXY-OprM | 0.45 | Efflux pump — aminoglycoside/FQ resistance in *Pseudomonas* |
| AcrAB-TolC | 0.40 | Gram-negative efflux — modest broad-spectrum resistance |
| Global efflux | 0.35 | Non-specific efflux — weakest single mechanism |
| As-yet-unknown 1–3 | 0.50 each | Calibration placeholders |



### 7.3 Resistance emergence

New resistance arises in the model only when a patient is **actively receiving antibiotics**. This is biologically accurate — antibiotic exposure creates the selective pressure that favours resistant mutants. Without antibiotics, a random mutation conferring resistance offers no survival advantage and is unlikely to become established.

**The "danger zone" — why sub-therapeutic dosing drives resistance:**

The probability of resistance emerging is not simply proportional to drug concentration. Instead, it follows a bell-shaped curve that peaks when drug levels are at **roughly half the therapeutic concentration**:

- **Very low drug levels:** No selective pressure — susceptible and resistant bacteria coexist equally. Resistance has no advantage.
- **Sub-therapeutic levels (the danger zone):** The drug kills susceptible bacteria, removing competition, but is too weak to eliminate the resistant mutant. The mutant thrives in the ecological vacuum. *This is the most dangerous window for resistance emergence.*
- **Full therapeutic levels:** Both susceptible and resistant bacteria are suppressed. Even if a mutant arises, the drug concentration is high enough to limit its growth.

This is why incomplete antibiotic courses, poor adherence, and underdosing are such powerful drivers of AMR — they create prolonged periods in the danger zone.

**The emergence formula:**

Each day, for each active infection under antibiotic treatment, the model calculates:

```
emergence_rate = mechanism_rate × bacteria_level_factor × drug_factor × multi_drug_penalty
```

| Factor | What it represents | Clinical analogy |
|--------|-------------------|------------------|
| `mechanism_rate` | How biologically likely this bacterium is to acquire this specific mechanism | Some mutations are common (e.g., *gyrA* point mutations); others are extremely rare (e.g., acquiring NDM by conjugation) |
| `bacteria_level_factor` | Logarithmic scaling by bacterial load | A bloodstream infection with 10⁸ bacteria generates more mutants per day than a colonisation with 10⁴ |
| `drug_factor` | Gaussian curve peaking at 50% of therapeutic concentration | The "danger zone" — resistance is most likely when drug levels are sub-therapeutic |
| `multi_drug_penalty` | Suppression when multiple drugs are used together | Combination therapy (e.g., meropenem + amikacin) makes it much harder for a single mechanism to confer survival |


#### Incidence band multipliers — scaling for population realism

In a simulation of 100,000 people, common bacteria like *E. coli* (carried by nearly everyone) would accumulate resistance mutations unrealistically fast without adjustment. Real-world constraints — clonal competition, immune clearance, fitness costs — prevent this in nature, but the model must explicitly correct for it.

The solution is to assign each bacterium to an "incidence band" that scales its mutation rate:

| Band | Multiplier | Rationale | Bacteria |
|------|-----------|-----------|----------|
| **High incidence** | ×0.1 | Very common organisms — rate scaled down 10-fold | *E. coli*, *S. aureus*, *S. pneumoniae*, *K. pneumoniae*, *H. pylori*, *H. influenzae*, *C. trachomatis*, *M. pneumoniae*, *S. epidermidis* |
| **Moderate incidence** | ×1.0 | Reference rate — no adjustment | *S. pyogenes*, *N. gonorrhoeae*, *C. jejuni*, *Enterobacter spp.*, *E. cloacae*, *Proteus spp.*, *S. agalactiae*, *T. pallidum*, *M. catarrhalis*, *M. genitalium* |
| **Low incidence** | ×3.0 | Uncommon organisms — rate scaled up to ensure resistance emerges during the simulation | *P. aeruginosa*, *C. difficile*, *B. fragilis*, *Citrobacter spp.*, *Serratia spp.*, *Salmonella (non-typhoidal)*, *Shigella spp.*, *E. faecalis*, *B. pertussis*, *L. pneumophila* |
| **Very low incidence** | ×10.0 | Rare pathogens — strong upscaling needed | *A. baumannii*, *E. faecium*, *MDR M. tuberculosis*, *N. meningitidis*, *S. Typhi*, *S. Paratyphi A*, *V. cholerae*, *L. monocytogenes*, *Y. enterocolitica*, *S. maltophilia*, *Morganella spp.*, *P. stuartii*, *B. cepacia complex* |

Most mechanism emergence rates are set to `0.0` for biologically impossible combinations (e.g., *S. pyogenes* acquiring NDM carbapenemase). When a resistance event is triggered, there is an 80% chance the specific mechanism is assigned (`mechanism_assignment_probability` = 0.8); otherwise it becomes a generalised resistance state.




### 7.4 Resistance reversion and fitness costs

Resistance is not free. Maintaining resistance mechanisms costs the bacterium energy and resources — like carrying a heavy suitcase through an airport. In the absence of antibiotics, resistant bacteria grow more slowly than their susceptible competitors and are gradually outcompeted. This is why resistance can decline after antibiotic use is reduced — a key insight for stewardship policy.

The model assigns each mechanism a daily **reversion rate** — the probability of losing resistance per day when no antibiotic pressure is present. Higher rates mean the mechanism is "expensive" and lost quickly; lower rates mean it is nearly cost-free and persists indefinitely.

Key patterns:
- **Most stable:** Single point mutations (e.g., *gyrA* fluoroquinolone resistance, reversion 0.0001/day) — the mutation barely affects the bacterium's fitness, so it persists for years even without ciprofloxacin pressure
- **Least stable:** Complex multi-gene cassettes (e.g., VanA/VanB vancomycin resistance, reversion 0.002/day; *rpoB* rifampicin resistance, 0.002/day) — these impose significant metabolic costs and are lost relatively quickly without glycopeptide or rifampicin exposure
- **Default** for non-mechanism-specific resistance: 0.0004/day

The full reversion rates by mechanism category:

### Enzymatic Inactivation
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| **KPC** (*bla*~KPC~) | `0.001` | Plasmid-mediated carbapenemase; moderate maintenance cost. |
| **NDM / VIM** | `0.0015` | Metallo-$eta$-lactamases, frequently on large, high-burden mobile genetic elements. |
| **OXA-48** | `0.0005` | Class D carbapenemase; comparatively lower fitness burden. |
| **ESBL CTX-M / TEM / SHV** | `0.0006` | Standard extended-spectrum $eta$-lactamases. |
| **AmpC DHA** | `0.0006` | Plasmid-mediated AmpC; typical cost profile. |
| **AmpC CMY** | `0.0001` | Often native gene upregulation; minimal fitness loss to maintain. |
| **FosA** | `0.0005` | Plasmid-mediated fosfomycin resistance; moderate cost. |
| **CAT** | `0.0005` | Chloramphenicol acetyltransferase. |
| **16S rRMTase** | `0.0005` | Ribosomal RNA methyltransferases conferring high-level aminoglycoside resistance. |

### Target Site Alterations
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| **PBP2a / *mecA*** | `0.0009` | High energetic cost associated with maintaining the staphylococcal cassette chromosome *mec* (SCC*mec*). |
| ***erm(B)*** | `0.002` | High reversion rate; target methylation for macrolide-lincosamide-streptogramin B (MLS~B~) resistance. |
| **VanA / VanB** | `0.002` | Highly complex target reprograming (D-Ala-D-Ala to D-Ala-D-Lac); significant energetic drain in the absence of glycopeptide exposure. |
| **CFR** | `0.0005` | RNA methyltransferase (oxazolidinone/phenicol cross-resistance). |

### Structural Mutations
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| ***gyrA* (Primary)** | `0.0001` | Single-step topoisomerase mutation; essentially stable with negligible fitness penalty. |
| ***parC* (Secondary)** | `0.0002` | Secondary topoisomerase IV mutations; compounding structural cost. |
| **Folate Pathway** | `0.0001` | Low cost; largely integron-associated (e.g., *sul* or *dfrA* elements). |
| **Nitroreductase** | `0.0003` | Loss-of-function mutations affecting nitrofurantoin activation; moderate cost. |
| ***mprF*** | `0.001` | Membrane lipid modification (daptomycin resistance); structural cell wall alterations bear a distinctive fitness cost. |
| ***rpoB*** | `0.002` | High fitness cost due to structurally significant alterations in RNA polymerase (rifampicin resistance). |

### Target Protection & Target Modification
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| ***mcr-1*** | `0.0015` | Plasmid-mediated phosphoethanolamine transferase (colistin resistance); substantial lipid A modification burden. |
| **Tet(M)** | `0.0005` | Ribosomal protection protein; moderate cost to maintain on transposons (e.g., Tn*916*). |
| **Qnr** | `0.0001` | Plasmid-mediated quinolone resistance protein; relatively stable. |
| **FusB** | `0.0005` | Target protection mechanism for fusidic acid resistance. |

### Porin Loss & Efflux Pumps
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| **OprD Loss** | `0.0005` | Loss of outer membrane channel (carbapenem resistance); hinders nutrient acquisition. |
| **OmpK35 / OmpK36 Loss** | `0.0005` | Analogous mechanism in Enterobacterales. |
| **Global / Generic Porin Loss** | `0.0005` | Broad permeability phenotype cost. |
| **AcrAB-TolC** | `0.0005` | Overexpression of major RND-family efflux pump complex. |
| **MexXY-OprM** | `0.0005` | Endogenous efflux system upregulation (common in *Pseudomonas aeruginosa*). |
| **Global / Generic Efflux** | `0.0005` | Broad, non-specific transport energy costs. |

*Note: The system reserves three placeholder variables (`as_yet_unknown_1-3`, baseline rate `0.001`) designated for future, non-prescribed clinical trials or empirical calibration against newly emergent paradigms.*

### 7.5 Resistance floors

Some bacteria are **intrinsically resistant** to certain antibiotics — *Stenotrophomonas maltophilia*, for example, produces L1/L2 metallo-β-lactamases that make it naturally resistant to carbapenems. This is not acquired resistance; every isolate has it.

In a simulation of only 100,000 people, rare pathogens like *S. maltophilia* produce so few infections that their resistance levels can randomly drift to zero — a modelling artefact, not real biology. To prevent this, the model enforces **resistance floors**: minimum resistance levels that certain organisms cannot drop below.

| Parameter | Value | What it does |
|-----------|-------|-------------|
| `resistance_floor_feature_enabled` | 1.0 (on) | Master switch for the floor system |
| `bacteria_{name}_resistance_floor_enabled` | Per-organism | Turns floors on for specific species |
| `bacteria_{name}_resistance_floor_ramp_years` | 10.0 | Years to reach full floor level after drug introduction |
| `bacteria_{name}_{drug_class}_resistance_floor` | 0.0–1.0 | The minimum resistance prevalence enforced |

Currently configured:

- ***S. maltophilia***: **Enabled** — preserves known intrinsic resistance to carbapenems and cephalosporins
- ***E. faecium***: **Disabled** (infrastructure exists but not active in current configuration)



### 7.6 Cross-resistance groups

When a bacterium becomes resistant to one antibiotic, it often becomes resistant to related antibiotics at the same time. For instance, if *E. coli* acquires an ESBL enzyme and becomes resistant to ceftriaxone, you would expect it to also be resistant to other cephalosporins (cefazolin, cefuroxime) and penicillins — because the enzyme destroys the same β-lactam ring in all of them.

The model captures this by defining **cross-resistance groups**: sets of drugs for which resistance is always acquired and lost together. When an individual's *E. coli* becomes resistant to any drug in Group 1 (β-lactams), it simultaneously becomes resistant to all drugs in that group.

The table below shows all cross-resistance groups for the major organisms (selected examples — the full list covers all 42 bacteria):

| Bacteria | Group | Drugs sharing resistance |
|----------|-------|------------------------|
| E. coli | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Cephalexin, Cefazolin, Cefuroxime, Ceftriaxone, Amoxicillin Clavulanate, Ampicillin Sulbactam, Piperacillin Tazobactam, Ticarcillin Clavulanate |
| E. coli | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| E. coli | Group 3 | Gentamicin, Tobramycin, Amikacin |
| A. baumannii | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Cephalexin, Cefazolin, Cefuroxime, Amoxicillin Clavulanate, Ampicillin Sulbactam, Piperacillin Tazobactam, Ticarcillin Clavulanate |
| A. baumannii | Group 2 | Meropenem, Imipenem C, Ertapenem, Meropenem Vaborbactam |
| A. baumannii | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin |
| A. baumannii | Group 4 | Gentamicin, Tobramycin, Amikacin |
| K. pneumoniae | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Cephalexin, Cefazolin, Cefuroxime, Ceftriaxone, Amoxicillin Clavulanate, Ampicillin Sulbactam, Piperacillin Tazobactam, Ticarcillin Clavulanate |
| K. pneumoniae | Group 2 | Meropenem, Imipenem C, Ertapenem, Meropenem Vaborbactam |
| K. pneumoniae | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. pneumoniae | Group 1 | Erythromycin, Azithromycin, Clarithromycin |
| S. pneumoniae | Group 2 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| S. aureus | Group 1 | Penicillin G, Ampicillin, Amoxicillin |
| S. aureus | Group 2 | Cephalexin, Cefazolin, Cefuroxime, Ceftriaxone |
| S. aureus | Group 3 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| S. epidermidis | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Cephalexin, Cefazolin, Cefuroxime, Ceftriaxone |
| S. epidermidis | Group 2 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| S. epidermidis | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. maltophilia | Group 1 | Trim Sulf |
| S. maltophilia | Group 2 | Tetracycline, Doxycycline, Minocycline |
| S. maltophilia | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| P. aeruginosa | Group 1 | Piperacillin, Piperacillin Tazobactam, Ceftazidime, Ceftazidime Avibactam, Cefepime |
| P. aeruginosa | Group 2 | Meropenem, Meropenem Vaborbactam, Imipenem C |
| P. aeruginosa | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| P. aeruginosa | Group 4 | Gentamicin, Tobramycin, Amikacin |
| E. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| E. spp. | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| E. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Mdr mycobacterium tuberculosis | Group 1 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Mdr mycobacterium tuberculosis | Group 2 | Gentamicin, Tobramycin, Amikacin |
| Mdr mycobacterium tuberculosis | Group 3 | Rifampicin |
| Mdr mycobacterium tuberculosis | Group 4 | Linezolid |
| E. faecalis | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Piperacillin, Piperacillin Tazobactam, Ticarcillin, Ticarcillin Clavulanate |
| E. faecalis | Group 2 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| E. faecalis | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| E. faecalis | Group 4 | Tetracycline, Doxycycline, Minocycline |
| E. faecalis | Group 5 | Vancomycin, Teicoplanin, Dalbavancin |
| E. faecium | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Piperacillin, Piperacillin Tazobactam |
| E. faecium | Group 2 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| E. faecium | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| E. faecium | Group 4 | Tetracycline, Doxycycline, Minocycline |
| E. faecium | Group 5 | Vancomycin, Teicoplanin, Dalbavancin |
| C. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| C. spp. | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| C. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| C. spp. | Group 4 | Gentamicin, Tobramycin, Amikacin |
| E. cloacae | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| E. cloacae | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| E. cloacae | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| M. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| M. spp. | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| M. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| P. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| P. spp. | Group 2 | Ceftriaxone, Ceftazidime |
| P. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| S. spp. | Group 2 | Ceftriaxone, Ceftazidime, Cefepime |
| S. spp. | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. spp. | Group 4 | Gentamicin, Tobramycin, Amikacin |
| P. stuartii | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin |
| P. stuartii | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Salmonella enterica serovar typhi | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate, Cephalexin, Cefazolin, Cefuroxime |
| Salmonella enterica serovar typhi | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Salmonella enterica serovar typhi | Group 3 | Ceftriaxone, Ceftazidime |
| Salmonella enterica serovar paratyphi a | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| Salmonella enterica serovar paratyphi a | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Salmonella enterica serovar paratyphi a | Group 3 | Ceftriaxone, Ceftazidime |
| Invasive non-typhoidal salmonella spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| Invasive non-typhoidal salmonella spp. | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Invasive non-typhoidal salmonella spp. | Group 3 | Ceftriaxone, Ceftazidime |
| S. spp. | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| S. spp. | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. spp. | Group 3 | Ceftriaxone, Ceftazidime |
| S. spp. | Group 4 | Tetracycline, Doxycycline |
| V. cholerae | Group 1 | Tetracycline, Doxycycline |
| V. cholerae | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| V. cholerae | Group 3 | Erythromycin, Azithromycin, Clarithromycin |
| C. jejuni | Group 1 | Erythromycin, Azithromycin, Clarithromycin |
| C. jejuni | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| C. jejuni | Group 3 | Tetracycline, Doxycycline |
| Y. enterocolitica | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| Y. enterocolitica | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| Y. enterocolitica | Group 3 | Tetracycline, Doxycycline |
| H. pylori | Group 1 | Clarithromycin, Erythromycin, Azithromycin |
| H. pylori | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| H. pylori | Group 3 | Tetracycline, Doxycycline |
| S. pyogenes | Group 1 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| S. pyogenes | Group 2 | Tetracycline, Doxycycline |
| S. pyogenes | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| S. agalactiae | Group 1 | Erythromycin, Azithromycin, Clarithromycin, Clindamycin |
| S. agalactiae | Group 2 | Tetracycline, Doxycycline |
| S. agalactiae | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| H. influenzae | Group 1 | Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| H. influenzae | Group 2 | Erythromycin, Azithromycin, Clarithromycin |
| H. influenzae | Group 3 | Ciprofloxacin, Levofloxacin |
| H. influenzae | Group 4 | Tetracycline, Doxycycline |
| M. catarrhalis | Group 1 | Ampicillin, Amoxicillin, Penicillin G, Piperacillin, Ticarcillin, Amoxicillin Clavulanate, Ampicillin Sulbactam, Piperacillin Tazobactam, Ticarcillin Clavulanate |
| M. catarrhalis | Group 2 | Erythromycin, Azithromycin, Clarithromycin |
| N. gonorrhoeae | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| N. gonorrhoeae | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| N. gonorrhoeae | Group 3 | Tetracycline, Doxycycline |
| N. gonorrhoeae | Group 4 | Erythromycin, Azithromycin |
| N. gonorrhoeae | Group 5 | Ceftriaxone, Cefixime |
| N. meningitidis | Group 1 | Penicillin G, Ampicillin, Amoxicillin, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| N. meningitidis | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| N. meningitidis | Group 3 | Rifampicin |
| C. difficile | Group 1 | Vancomycin |
| C. difficile | Group 2 | Metronidazole |
| B. fragilis | Group 1 | Metronidazole |
| B. fragilis | Group 2 | Clindamycin |
| B. fragilis | Group 3 | Meropenem, Imipenem C |
| L. monocytogenes | Group 1 | Ampicillin, Amoxicillin, Penicillin G, Ampicillin Sulbactam, Amoxicillin Clavulanate |
| L. monocytogenes | Group 2 | Tetracycline, Doxycycline |
| C. trachomatis | Group 1 | Azithromycin, Erythromycin, Clarithromycin |
| C. trachomatis | Group 2 | Doxycycline, Tetracycline |
| C. trachomatis | Group 3 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| M. genitalium | Group 1 | Azithromycin, Erythromycin, Clarithromycin |
| M. genitalium | Group 2 | Ciprofloxacin, Levofloxacin, Moxifloxacin, Ofloxacin |
| M. genitalium | Group 3 | Doxycycline, Tetracycline |
| T. pallidum | Group 1 | Erythromycin, Azithromycin, Clarithromycin |
| T. pallidum | Group 2 | Doxycycline, Tetracycline |
| B. pertussis | Group 1 | Erythromycin, Azithromycin, Clarithromycin |


---

## 8. Microbiome and Carriage

Most bacteria in and on the human body are harmless commensals — they live on skin, in the gut, and in the respiratory tract without causing disease. This is **carriage** (or colonisation), and it is the normal state. Only a small fraction of carried bacteria ever cause active infection.

However, carriage is critically important for AMR because the microbiome is where resistance is **stored and exchanged**. A patient who was treated with ciprofloxacin last month may still carry ciprofloxacin-resistant *E. coli* in their gut microbiome. If they develop a UTI from that same resistant strain, the empiric therapy may fail.


### 8.1 Carriage compartments

Each bacterium in the model has a designated ecological niche — where it naturally lives in (or on) the body:

| Compartment | Example bacteria | Clinical relevance |
|-------------|-----------------|-------------------|
| Gut | *E. coli*, *K. pneumoniae*, *Enterococcus spp.*, *Shigella*, *Salmonella*, *C. difficile* | Largest reservoir; disrupted by broad-spectrum antibiotics |
| Respiratory | *S. pneumoniae*, *H. influenzae*, *P. aeruginosa*, *A. baumannii*, *M. catarrhalis*, *M. tuberculosis* | Carriage often precedes pneumonia |
| Skin/Soft tissue | *S. aureus*, *S. epidermidis* | Nasal/skin MRSA carriage drives surgical wound infections |
| Genitourinary | *N. gonorrhoeae*, *C. trachomatis*, *M. genitalium*, *T. pallidum*, *S. agalactiae* | Asymptomatic STI carriage enables transmission |


### 8.2 Resistance in the microbiome

The microbiome serves as a hidden reservoir of resistance. Each individual carries a resistance tracking matrix for every organism, even when no clinical infection is present.

Key dynamics:

| Process | Parameter | Value | What it means |
|---------|-----------|-------|---------------|
| Established colonies harder to clear | `carriage_duration_log_odds_coefficient` | -0.01/day (caps at -2.0) | The longer a resistant strain has been carried, the harder it is to eradicate — mature colonies are ~7× harder to clear than newly acquired ones |
| Auto-infection dampening | `infection_from_microbiome_dampening` | 0.70 | When a carrier develops an infection from their own gut flora, the starting bacterial load is reduced to 70% (not all commensal bacteria transition to pathogens) |
| Silent mutation rate | `microbiome_resistance_emergence_rate_per_day_baseline` | 1.0e-20 (effectively zero) | Resistance does NOT evolve spontaneously in asymptomatic carriage — antibiotic exposure is required |

The last point is an important design decision: the model says that AMR is driven by antibiotic-treated infections, not by silent mutation in commensal flora.


## 9. Horizontal Gene Transfer (HGT)

Bacteria can share resistance genes directly with each other — even between different species. This is **horizontal gene transfer (HGT)**, and it is the main reason a resistance gene that evolves in one species can rapidly appear in others. Clinically, this is why you see the same ESBL plasmids in *E. coli*, *Klebsiella*, and *Proteus* from the same hospital ward.


### 9.1 Transfer compatibility

Not all bacteria can exchange genes equally. The model uses a compatibility matrix that reflects known biological barriers:

| Transfer type | Probability | Example |
|--------------|-------------|---------|
| Same species | 0.95 | *E. coli* → *E. coli* — virtually unobstructed plasmid exchange |
| Related Gram-negatives | 0.80–0.90 | *E. coli* → *Klebsiella* — highly fluid exchange (shared plasmid types) |
| Unrelated Gram-positives | 0.10–0.20 | *S. aureus* → *Enterococcus* — rare transposon-mediated events |
| Gram-negative → Gram-positive | 0.0 | Biologically prohibited (fundamentally different cell wall architectures) |


### 9.2 The HGT process

Each day, for every individual carrying resistant bacteria in their microbiome, the model evaluates potential gene transfer events:

| Step | Parameter | Value | Clinical parallel |
|------|-----------|-------|-------------------|
| Base transfer rate | `microbiome_resistance_transfer_probability_per_day` | 0.0001 | Background rate — equivalent to a conjugation event occurring every ~27 years per carrier, reflecting how rare HGT is without antibiotic pressure |
| Amplification during antibiotic therapy | `hgt_antibiotic_pressure_multiplier` | 1.50 (×1.5) | Antibiotic stress triggers the bacterial SOS response, which activates mobile genetic elements and increases conjugation rates by 50% — one of the reasons antibiotic use drives resistance even beyond the target pathogen |



## 10. Mortality

The model tracks mortality from both background (non-infection) causes and from active bacterial infections. This matters for AMR projections because resistant infections carry higher mortality — if empiric therapy fails, the patient stays infected longer, and the risk of death accumulates.


### 10.1 Background mortality

Everyone faces a baseline mortality risk that increases with age:

| Factor | Parameter | Value | What it means |
|--------|-----------|-------|---------------|
| Aging penalty | `log_odds_mortality_per_year_of_age` | 0.04 | Each year of age adds ~4% relative increase in daily death risk (exp(0.04) = 1.04) |
| Elderly frailty acceleration | `log_odds_mortality_per_year_of_age_squared` | 0.05 | A quadratic term that makes mortality rise faster above ~70 — capturing how an 85-year-old is much frailer than a 65-year-old |

These parameters operate on a log-odds scale, so they compound multiplicatively over time.


### 10.2 Infection mortality

When a person has an active symptomatic infection (above the `symptomatic_infection_level_threshold`), the model evaluates a daily mortality check. The risk depends on three things:

**1. How dangerous the bacterium is** (base mortality intercept, on a log-odds scale):

| Category | Example | Intercept | Clinical context |
|----------|---------|-----------|-----------------|
| Mild | *C. trachomatis* | -15.0 | Almost never fatal — values this negative make daily death probability near zero |
| Moderate | *E. coli* | -9.0 | Fatal mainly in elderly, immunosuppressed, or when resistance delays treatment |
| Severe | *S. pneumoniae* | -6.0 | Untreated pneumococcal pneumonia historically killed ~30% |
| Critical | *N. meningitidis* (-1.2), *M. tuberculosis* (-0.8) | — | Can kill rapidly without effective therapy |

**2. Which body site is infected** (syndrome multiplier):

These multipliers reflect how anatomically dangerous each infection site is:

| Syndrome | Multiplier | Rationale |
|----------|-----------|-----------|
| Genital | 0.05 | Rarely fatal (localised mucosal infections) |
| Skin / Ear | 0.1 | Low systemic risk unless secondary bacteraemia |
| UTI | 0.5 | Usually self-limiting but can ascend to urosepsis |
| Bone/Joint | 0.8 | Serious but slow-progressing; mortality from surgical complications |
| Intra-abdominal | 1.5 | Peritonitis carries high mortality even with surgery |
| Respiratory | 1.5 | Pneumonia — leading infectious cause of death globally |
| CNS | 3.0 | Meningitis/brain abscess — poor penetration of many antibiotics |
| Bloodstream | 4.0 | Bacteraemia/sepsis — the most immediately life-threatening |

**3. Sepsis override**: If the infection progresses to sepsis (organ dysfunction from uncontrolled infection), normal mortality limits are overridden and the model applies an aggressively escalated daily death risk. This reflects the clinical reality that sepsis without effective antibiotics is rapidly fatal — and resistant organisms that are untreatable with empiric therapy are exactly the scenario where this matters most.



## 11. Policy Evaluation

The purpose of this model is not just to simulate AMR — it is to compare what happens under different policy choices. To do this, the simulation runs a single shared history and then **branches** into parallel futures, each with different antibiotic prescribing rules.


### 11.1 How branching works

At a configurable year (default: **2027**), the simulation saves a complete snapshot of the entire population — every person's age, infections, microbiome resistance, treatment history, everything. It then runs three independent scenarios forward from that identical starting point:

| Branch | What it represents |
|--------|--------------------|
| **Baseline** | Business as usual — no policy changes. This is the "do nothing" future. |
| **Stewardship** | Antibiotic stewardship interventions are introduced: narrower prescribing, more diagnostic testing, stronger disincentives for reserve drugs. |
| **Counterfactual** | A hypothetical world where resistance is eliminated at the branch point. This lets you measure the total burden attributable to AMR by comparing outcomes against this "no resistance" scenario. |

Because all three branches start from an identical population state, any differences in outcomes (deaths, treatment failures, resistance prevalence) are **causally attributable** to the policy differences alone.


### 11.2 Policy parameters

Each branch can adjust the following parameters. A dash (—) means the parameter is left at its default value:

| Parameter | Baseline | Stewardship | Counterfactual | What it does |
|-----------|----------|-------------|----------------|-------------|
| `drug_selection_temperature` | — | x0.65 | — | Makes prescribing more deterministic (less random variation in drug choice) — stewardship guidelines reduce idiosyncratic prescribing |
| `bacterial_testing_rate_multiplier` | — | x1.5 | — | 50% more bacterial cultures ordered — better pathogen identification |
| `resistance_testing_rate_multiplier` | — | x1.5 | — | 50% more susceptibility testing — clinicians know which drugs will work |
| `reserve_drug_penalty_multiplier` | — | x2.0 | — | Doubles the prescribing barrier for last-resort antibiotics — discourages casual use of carbapenems, colistin, etc. |
| `drug_initiation_rate_multiplier` | — | x0.85 | — | 15% fewer antibiotic courses started — reflects "watchful waiting" and avoiding unnecessary prescriptions |
| `drug_cessation_rate_multiplier` | — | x1.2 | — | Courses are 20% shorter on average — reflecting evidence that shorter courses are often equally effective |
| `counterfactual_resistance_multiplier` | — | — | 0.0 | Multiplies all resistance levels by zero — instantly creates a resistance-free world |
| `clear_all_resistance_on_branch_start` | false | false | true | Wipes all microbiome and infection resistance at the branch point |


### 11.3 Key simulation constants

| Constant | Value | Purpose |
|----------|-------|---------|
| `SIMULATION_START_YEAR` | 1930.0 | Calendar year at day 0 — early enough to capture the pre-antibiotic era |
| `POLICY_BRANCH_YEAR` | 2027.0 | Year when the three policy branches diverge |
| `INFECTION_EPS` | 0.001 | Minimum meaningful infection level (below this, the infection is treated as cleared) |
| `MICROBIOME_MAJORITY_THRESHOLD` | 0.5 | If >50% of a species in the microbiome carries a resistance mechanism, it is classified as "majority resistant" |
| `MAX_MECHANISM_PROFILES` | 200 | Reservoir sample size per bacteria for mechanism profile caching (performance optimisation) |

---



## 12. Known Limitations

Every model is a simplification of reality. These are the main areas where this model knowingly trades accuracy for tractability:

1. **Abstract drug levels**: Antibiotic concentrations are modelled as dimensionless units rather than true pharmacokinetic concentrations (mg/L). This means we capture the *relative* dynamics of drug accumulation and clearance, but cannot directly compare model values to MIC breakpoints from a clinical microbiology report.

2. **No explicit strain competition**: Within the microbiome, resistant and susceptible strains do not explicitly compete for resources. Instead, resistance is promoted by antibiotic pressure and decays through reversion rates. This means the model cannot capture scenarios where a fitness-cost-free resistant strain permanently outcompetes susceptible strains in the absence of antibiotics.

3. **No within-host spatial structure**: Infections are treated as homogeneous within a body compartment. Biofilm formation, abscess walling-off, and planktonic-vs-sessile distinctions are not modelled. In reality, biofilm-embedded bacteria can survive antibiotic concentrations 100–1000× higher than planktonic cells.

4. **Static vaccine model**: Vaccinated individuals have a fixed proportional reduction in infection risk. Vaccine effects do not depend on background prevalence (no herd immunity dynamics), and vaccine-driven serotype replacement is not captured.

5. **Broad regional groupings**: The model uses continental-level regions (e.g., "Europe", "Africa") rather than country-level or hospital-level variation. Antibiotic consumption patterns and resistance rates can vary dramatically between countries within the same region.

---



## Appendix A — Bacteria, Drugs, Mechanisms and Enums

This appendix lists every entity in the model. Use it as a lookup reference when you encounter a specific bacterium, drug, or mechanism code in the main text.



### A.1 Bacteria (42 species)

| Index | Species | Group | Carriage compartment |
|-------|---------|-------|---------------------|
| 0 | Acinetobacter baumannii | NonFermenter | Respiratory |
| 1 | Citrobacter spp. | Enterobacterales | Gut |
| 2 | Enterobacter spp. | Enterobacterales | Gut |
| 3 | Enterococcus faecalis | GramPositive | Gut |
| 4 | Enterococcus faecium | GramPositive | Gut |
| 5 | Escherichia coli | Enterobacterales | Gut |
| 6 | Klebsiella pneumoniae | Enterobacterales | Gut |
| 7 | Morganella spp. | Enterobacterales | Gut |
| 8 | Proteus spp. | Enterobacterales | Gut |
| 9 | Serratia spp. | Enterobacterales | Gut |
| 10 | Providencia stuartii | Enterobacterales | Genitourinary |
| 11 | Pseudomonas aeruginosa | NonFermenter | Respiratory |
| 12 | Stenotrophomonas maltophilia | NonFermenter | Respiratory |
| 13 | Staphylococcus aureus | GramPositive | Skin/Soft Tissue |
| 14 | Staphylococcus epidermidis | GramPositive | Skin/Soft Tissue |
| 15 | Streptococcus pneumoniae | GramPositive | Respiratory |
| 16 | Salmonella enterica serovar Typhi | Enterobacterales | Gut |
| 17 | Salmonella enterica serovar Paratyphi A | Enterobacterales | Gut |
| 18 | Invasive non-typhoidal Salmonella spp. | Enterobacterales | Gut |
| 19 | Shigella spp. | Enterobacterales | Gut |
| 20 | Neisseria gonorrhoeae | Fastidious | Genitourinary |
| 21 | Streptococcus pyogenes | GramPositive | Respiratory |
| 22 | Streptococcus agalactiae | GramPositive | Genitourinary |
| 23 | Haemophilus influenzae | Fastidious | Respiratory |
| 24 | Chlamydia trachomatis | Fastidious | Genitourinary |
| 25 | Mycoplasma genitalium | Fastidious | Genitourinary |
| 26 | Vibrio cholerae | EntericPathogen | Gut |
| 27 | Neisseria meningitidis | Fastidious | Respiratory |
| 28 | Listeria monocytogenes | GramPositive | Gut |
| 29 | Clostridioides difficile | Anaerobe | Gut |
| 30 | Bacteroides fragilis | Anaerobe | Gut |
| 31 | Campylobacter jejuni | Helicobacter | Gut |
| 32 | Enterobacter cloacae | Enterobacterales | Gut |
| 33 | Yersinia enterocolitica | Enterobacterales | Gut |
| 34 | Moraxella catarrhalis | Fastidious | Respiratory |
| 35 | Treponema pallidum | Spirochete | Genitourinary |
| 36 | Bordetella pertussis | Fastidious | Respiratory |
| 37 | Helicobacter pylori | Helicobacter | Gut |
| 38 | MDR Mycobacterium tuberculosis | Mycobacteria | Respiratory |
| 39 | Mycoplasma pneumoniae | Fastidious | Respiratory |
| 40 | Legionella pneumophila | Fastidious | Respiratory |
| 41 | Burkholderia cepacia complex | NonFermenter | Respiratory |



### A.2 Antibiotics (58 drugs in 36 classes)

| Drug | Class |
|------|-------|
| sulfanilamide | Sulfonamides |
| penicillin_g | Penicillins |
| ampicillin | Penicillins |
| amoxicillin | Penicillins |
| piperacillin | Penicillins |
| ticarcillin | Penicillins |
| cephalexin | Cephalosporins 1–2G |
| cefazolin | Cephalosporins 1–2G |
| cefuroxime | Cephalosporins 1–2G |
| ceftriaxone | Cephalosporins 3G |
| ceftazidime | Cephalosporins 3G |
| cefepime | Cephalosporins 4G |
| ceftaroline | Anti-MRSA Cephalosporins (5G) |
| ceftolozane_tazobactam | Cephalosporins 3G/BLI |
| cefiderocol | Siderophore Cephalosporins |
| meropenem | Carbapenems |
| imipenem_c | Carbapenems |
| ertapenem | Carbapenems |
| aztreonam | Monobactams |
| erythromycin | Macrolides |
| azithromycin | Macrolides |
| clarithromycin | Macrolides |
| clindamycin | Macrolides |
| gentamicin | Aminoglycosides |
| tobramycin | Aminoglycosides |
| amikacin | Aminoglycosides |
| ciprofloxacin | Fluoroquinolones |
| levofloxacin | Fluoroquinolones |
| moxifloxacin | Fluoroquinolones |
| ofloxacin | Fluoroquinolones |
| tetracycline | Tetracyclines |
| doxycycline | Tetracyclines |
| minocycline | Tetracyclines |
| tigecycline | Tetracyclines |
| vancomycin | Glycopeptides |
| teicoplanin | Glycopeptides |
| dalbavancin | Glycopeptides |
| linezolid | Oxazolidinones |
| tedizolid | Oxazolidinones |
| daptomycin | Other |
| quinu_dalfo | Other |
| trim_sulf | Sulfonamides |
| chloramphenicol | Chloramphenicol |
| nitrofurantoin | Other |
| fosfomycin | Other |
| retapamulin | Other |
| fusidic_a | Other |
| metronidazole | Other |
| fidaxomicin | Other |
| furazolidone | Other |
| rifampicin | Other |
| amoxicillin_clavulanate | BLI Combinations |
| piperacillin_tazobactam | BLI Combinations |
| ampicillin_sulbactam | BLI Combinations |
| ticarcillin_clavulanate | BLI Combinations |
| ceftazidime_avibactam | Novel BL/BLI |
| meropenem_vaborbactam | Novel BL/BLI |
| colistin | Polymyxins |



### A.3 Drug Classes (36)

| Code | Full name | Drugs |
|------|-----------|-------|
| `pen` | Penicillins | penicillin G, ampicillin, amoxicillin, piperacillin, ticarcillin |
| `bli` | BLI Combinations | amoxicillin-clavulanate, piperacillin-tazobactam, ampicillin-sulbactam, ticarcillin-clavulanate |
| `c1_2g` | Cephalosporins 1–2G | cephalexin, cefazolin, cefuroxime |
| `c3g` | Cephalosporins 3G | ceftriaxone, ceftazidime |
| `c3g_bli` | Cephalosporins 3G/BLI | ceftolozane-tazobactam |
| `c4g` | Cephalosporins 4G | cefepime |
| `c5g` | Anti-MRSA Cephalosporins (5G) | ceftaroline |
| `sid_ceph` | Siderophore Cephalosporins | cefiderocol |
| `bl_ni` | Novel BL/BLI | ceftazidime-avibactam, meropenem-vaborbactam |
| `carb` | Carbapenems | meropenem, imipenem, ertapenem |
| `mono` | Monobactams | aztreonam |
| `fq` | Fluoroquinolones | ciprofloxacin, levofloxacin, moxifloxacin, ofloxacin |
| `ag` | Aminoglycosides | gentamicin, tobramycin, amikacin |
| `mls` | Macrolides/Lincosamides | erythromycin, azithromycin, clarithromycin, clindamycin |
| `glyc` | Glycopeptides | vancomycin, teicoplanin, dalbavancin |
| `tet` | Tetracyclines | tetracycline, doxycycline, minocycline, tigecycline |
| `poly` | Polymyxins | colistin |
| `oxa` | Oxazolidinones | linezolid, tedizolid |
| `chl` | Chloramphenicol | chloramphenicol |
| `sulf` | Sulfonamides | sulfanilamide, trimethoprim-sulfamethoxazole |
| `other` | Other | daptomycin, quinupristin-dalfopristin, nitrofurantoin, fosfomycin, retapamulin, fusidic acid, metronidazole, fidaxomicin, furazolidone, rifampicin |



### A.4 Resistance Mechanisms (40)

See [Section 7.1](#71-resistance-mechanisms) for the full table.



### A.5 Enumerations



#### BacteriaGroup (9 groups)

| Group | Description |
|-------|-------------|
| `Enterobacterales` | Gram-negative enteric rods |
| `NonFermenter` | Non-fermenting Gram-negatives |
| `GramPositive` | Gram-positive cocci and rods |
| `Fastidious` | Fastidious Gram-negatives and atypicals |
| `EntericPathogen` | Specific enteric pathogens (V. cholerae) |
| `Anaerobe` | Obligate anaerobes |
| `Spirochete` | Spirochetes (T. pallidum) |
| `Helicobacter` | Helicobacter/Campylobacter |
| `Mycobacteria` | Mycobacteria (M. tuberculosis) |



#### CarriageCompartment (5)

`Gut`, `Respiratory`, `SkinSoftTissue`, `Genitourinary`, `Systemic`



#### ResistanceAcquisitionType (5)

| Type | Description |
|------|-------------|
| `AtInfection` | Acquired at community infection |
| `AtInfectionHosp` | Acquired at hospital infection |
| `AtInfectionTB` | Acquired at carrier-to-infection (treated-by) conversion |
| `DuringInfection` | De novo emergence during treatment |
| `HGT` | Horizontal gene transfer |



#### InfectionResolutionType (6)

`DrugTreatment`, `NaturalClearance`, `Death`, `SepsisDeath`, `ToxicityDeath`, `BackgroundDeath`



#### ImmunodeficiencyType (2)

`Temporary`, `Chronic`



#### AgeCategory (7)

`Infant`, `Preschool`, `SchoolAge`, `YoungAdult`, `MiddleAge`, `Elderly`, `NotYetBorn`



#### HospitalStatus (2)

`Community`, `Hospital`



#### Region (7)

`NorthAmerica`, `Europe`, `Asia`, `Oceania`, `SouthAmerica`, `Africa`, `Home` (fallback)



#### MicrobiomeResistanceLevel (4)

`None`, `Low`, `Medium`, `High`

---



## Appendix B — Parameter Reference

This appendix lists all ~120 top-level parameters that control the model's behaviour. Each parameter is a single number stored in the configuration. If you want to understand *why* the model produces a particular result, these are the numbers driving it. Parameters on a "log-odds" scale are described in the log-odds primer in Section 1.3.



### B.1 Global Scalars

These are the ~120 top-level parameters stored in the `GlobalScalars` struct:



#### Infection acquisition

| Parameter | Baseline Value | Description |
|----------|---------|-------------|
| `infection_growth_rate_per_day` | 0.1 | Daily bacterial growth increment |
| `infection_initial_level` | 1.0 | Starting bacterial load |
| `infection_clearance_threshold` | 0.5 | Level below which infection resolves |
| `infection_death_threshold` | 50.0 | Level at which death may occur |
| `symptom_onset_threshold` | 3.0 | Level for symptom development |
| `symptom_recheck_interval_days` | 7.0 | Re-evaluation interval |
| `symptom_onset_rate_per_day` | 0.1 | Base symptom development rate |
| `not_under_care_fraction` | 0.05 | Fraction not seeking medical care |



#### Antibiotic treatment initiation

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `antibiotic_initiation_base_log_odds` | −6.5 |
| `antibiotic_initiation_log_odds_symptomatic_infection` | 6.5 |
| `antibiotic_initiation_log_odds_sepsis` | 6.0 |
| `antibiotic_initiation_log_odds_test_identified` | 0.92 |
| `antibiotic_initiation_log_odds_already_on_drug` | 0.18 |
| `antibiotic_initiation_log_odds_immunodeficiency` | 2.08 |
| `antibiotic_initiation_log_odds_no_indication` | −1.05 |



#### Antibiotic treatment cessation

| Parameter | Baseline Value |
|----------|---------|
| `treatment_stop_improvement_threshold` | 2.0 |
| `treatment_stop_rate_per_day` | 0.03 |
| `treatment_duration_base_days` | 7.0 |



#### Drug efficacy

| Parameter | Baseline Value |
|----------|---------|
| `drug_effect_on_bacteria_per_day` | 0.5 |
| `drug_minimum_effective_level` | 0.1 |



#### Drug selection

| Parameter | Baseline Value |
|----------|---------|
| `empiric_therapy_broad_spectrum_bonus` | 0.85 |
| `empiric_therapy_ineffective_drug_penalty` | 0.001 |
| `targeted_therapy_narrow_spectrum_bonus` | 5.0 |
| `targeted_therapy_broad_spectrum_penalty` | 0.1 |
| `targeted_therapy_ineffective_drug_penalty` | 0.001 |
| `regional_resistance_penalty_very_high` | 0.3 |
| `regional_resistance_penalty_high` | 0.5 |
| `regional_resistance_penalty_moderate` | 0.8 |
| `regional_resistance_threshold_very_high` | 0.6 |
| `regional_resistance_threshold_high` | 0.45 |
| `regional_resistance_threshold_moderate` | 0.1 |



#### Testing

| Parameter | Baseline Value |
|----------|---------|
| `bacterial_testing_available_from_day` | 5,478 |
| `resistance_testing_available_from_day` | 9,131 |
| `test_delay_days` | 3.0 |
| `test_rate_per_day` | 0.2 |
| `prob_test_r_done` | 0.95 |
| `test_r_error_probability` | 0.02 |
| `bacterial_testing_base_rate_per_day` | 0.15 |
| `resistance_testing_base_rate_per_day` | 0.95 |
| `bacterial_testing_hospital_multiplier` | 8.0 |
| `resistance_testing_hospital_multiplier` | 5.0 |
| `testing_immunosuppressed_multiplier` | 2.5 |
| `testing_sepsis_multiplier` | 4.0 |
| `bacterial_testing_initial_adoption_rate` | 0.1 |
| `bacterial_testing_max_temporal_multiplier` | 1.0 |
| `resistance_testing_initial_adoption_rate` | 0.05 |
| `resistance_testing_max_temporal_multiplier` | 1.0 |



#### Sepsis onset

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `sepsis_baseline_log_odds` | −14.0 |
| `log_odds_sepsis_infection_level` | 0.8 |
| `log_odds_sepsis_infection_duration` | 0.005 |
| `log_odds_sepsis_onset_immunosuppressed` | 0.7 |
| `log_odds_sepsis_onset_hospitalized` | 0.5 |
| `log_odds_sepsis_onset_not_under_care` | 1.0 |
| `sepsis_age_log_odds_neonatal` | 1.10 |
| `sepsis_age_log_odds_pediatric` | 0.18 |
| `sepsis_age_log_odds_young_adult` | 0.0 |
| `sepsis_age_log_odds_elderly` | 0.69 |



#### Sepsis death

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `sepsis_death_base_log_odds` | −5.0 |
| `sepsis_death_log_odds_age_infant` | 1.1 |
| `sepsis_death_log_odds_age_child` | −0.7 |
| `sepsis_death_log_odds_age_adult` | 0.0 |
| `sepsis_death_log_odds_age_elderly` | 0.9 |
| `sepsis_death_log_odds_immunosuppressed` | 1.5 |
| `sepsis_death_log_odds_bacteria_level` | 0.35 |
| `sepsis_death_log_odds_duration` | 0.04 |
| `sepsis_death_log_odds_early_phase` | 0.8 |
| `sepsis_death_early_phase_days` | 3.0 |
| `sepsis_death_log_odds_not_under_care` | 1.4 |



#### Sepsis recovery

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `sepsis_base_log_odds_of_recovery_per_day` | 0.0 |
| `sepsis_log_odds_bacteria_level` | −0.3 |
| `sepsis_log_odds_in_hospital` | 0.8 |
| `sepsis_log_odds_age_infant` | −0.5 |
| `sepsis_log_odds_age_child` | 0.4 |
| `sepsis_log_odds_age_adult` | 0.0 |
| `sepsis_log_odds_age_elderly` | −0.7 |
| `sepsis_log_odds_immunosuppressed` | −1.0 |
| `sepsis_minimum_duration_days` | 1.0 |



#### Background mortality

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `background_mortality_baseline_log_odds` | −14.0 |
| `log_odds_mortality_per_year_of_age` | 0.04 |
| `log_odds_mortality_per_year_of_age_squared` | 0.05 |
| `mortality_baseline_1930_multiplier` | 3.0 |
| `mortality_baseline_2035_multiplier` | 1.0 |
| `mortality_improvement_half_life_years` | 35.0 |
| `log_odds_mortality_immunosuppressed` | 0.916 |
| `log_odds_mortality_hospitalized` | 0.262 |



#### Hospitalisation

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `hospitalization_base_log_odds` | −10.4 |
| `hospitalization_log_odds_per_age_year` | 0.02 |
| `hospitalization_log_odds_sepsis` | 4.4 |
| `hospitalization_log_odds_symptomatic_infection` | 2.5 |
| `hospitalization_symptomatic_infection_level_threshold` | 3.0 |
| `hospitalization_recovery_rate_per_day` | 0.28 |
| `hospitalization_max_days` | 30.0 |
| `hospitalization_prevent_discharge_with_sepsis` | 1.0 |



#### Immunodeficiency

| Parameter | Baseline Value |
|----------|---------|
| `temporary_immunosuppression_onset_rate_per_day` | 0.00005 |
| `temporary_immunosuppression_recovery_rate_per_day` | 0.01 |
| `chronic_immunosuppression_onset_rate_per_day` | 0.00006 |
| `chronic_immunosuppression_recovery_rate_per_day` | 0.0012 |
| `chronic_immunodeficiency_probability_age_0_1` | 0.3 |
| `chronic_immunodeficiency_probability_age_1_18` | 0.2 |
| `chronic_immunodeficiency_probability_age_18_65` | 0.4 |
| `chronic_immunodeficiency_probability_age_65_plus` | 0.6 |
| `antibiotic_infection_prevention_efficacy` | 0.7 |



#### Resistance

| Parameter | Baseline Value |
|----------|---------|
| `mechanism_assignment_probability` | 0.8 |



#### Microbiome and carriage

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `default_microbiome_clearance_probability_per_day` | 0.01 |
| `microbiome_clearance_probability_on_drug_treatment` | 0.8 |
| `default_microbiome_disruption_log_odds` | 0.3 |
| `microbiome_resistance_multiplier_on_acquisition` | 0.50 |
| `infection_from_microbiome_dampening` | 0.70 |
| `antibiotic_disruption_decay_half_life_days` | 30.0 |
| `carriage_duration_log_odds_coefficient` | −0.01 |
| `carriage_duration_max_log_odds_effect` | −2.0 |
| `antibiotic_clearance_log_odds_per_unit_activity` | 0.5 |
| `carrier_resistance_inheritance_probability` | 0.50 |
| `community_resistance_dilution_factor` | 0.50 |
| `majority_r_window_days` | 100 |
| `majority_r_min_total_samples` | 10 |
| `majority_r_freeze_at_last_positive` | 0.0 |
| `microbiome_majority_promotion_rate_per_day` | 0.02 |



#### Toxicity

| Parameter | Baseline (Log-odds ratio) |
|----------|---------|
| `default_drug_toxicity_death_hazard_per_unit_level` | 0.0 |
| `default_toxicity_reservoir_half_life_days` | 1.5 |
| `toxicity_death_base_log_odds` | −8.0 |
| `toxicity_death_log_odds_per_reservoir_unit` | 2.0 |
| `toxicity_death_log_odds_age_infant` | 0.6 |
| `toxicity_death_log_odds_age_child` | 0.2 |
| `toxicity_death_log_odds_age_adult` | 0.0 |
| `toxicity_death_log_odds_age_elderly` | 0.8 |
| `toxicity_death_log_odds_immunosuppressed` | 0.9 |
| `toxicity_death_log_odds_hospitalized` | 0.25 |
| `toxicity_discontinuation_base_log_odds` | −3.0 |
| `toxicity_discontinuation_log_odds_per_reservoir_unit` | 1.5 |
| `toxicity_discontinuation_log_odds_sepsis` | −1.5 |
| `toxicity_avoidance_penalty_multiplier` | 0.05 |
| `toxicity_avoidance_window_days` | 14.0 |



#### HGT

| Parameter | Baseline Value |
|----------|---------|
| `hgt_base_probability` | 1×10⁻⁵ |
| `hgt_co_infection_multiplier` | 10.0 |
| `hgt_hospital_multiplier` | 5.0 |
| `hgt_microbiome_multiplier` | 2.0 |
| `hgt_gut_compartment_multiplier` | 2.0 |



#### Travel

| Parameter | Baseline Value |
|----------|---------|
| `travel_probability_per_day` | 0.00005 |
| `north_america_travel_multiplier` | 3.0 |
| `europe_travel_multiplier` | 3.5 |
| `oceania_travel_multiplier` | 2.5 |
| `asia_travel_multiplier` | 1.5 |
| `south_america_travel_multiplier` | 0.8 |
| `africa_travel_multiplier` | 0.3 |



### B.2 Per-Bacteria Parameters (42 bacteria × N parameters)

Generated with key pattern `bacteria_{name}_{param}`:

| Parameter suffix | Description |
|------------------|-------------|
| `acquisition_log_odds` | Baseline acquisition rate |
| `hospital_acquisition_log_odds` | Hospital acquisition rate |
| `microbiome_vs_infection_log_odds` | Likelihood of carriage vs infection |
| `clearance_probability_per_day_no_treatment` | Natural clearance rate |
| `growth_rate_multiplier` | Bacteria-specific growth modifier |
| `symptom_onset_override` | Override for symptom onset threshold |
| `hgt_donor_probability` | Probability of being an HGT donor |
| `resistance_floor_enabled` | Whether resistance floors apply |
| `resistance_floor_ramp_years` | Ramp-up period for floors |
| `{drug_class}_resistance_floor` | Per-class floor target |



### B.3 Per-Drug Parameters (58 drugs × N parameters)

Generated with key pattern `drug_{name}_{param}`:

| Parameter suffix | Default | Description |
|------------------|---------|-------------|
| `half_life_days` | Drug-specific | PK half-life |
| `initial_level` | 10.0 | Administration level |
| `double_dose_multiplier` | 2.0 | Double-dose level |
| `spectrum_breadth` | 3.0 | Microbiome disruption breadth |
| `toxicity_death_hazard_per_unit_level` | 0.0 | Per-drug toxicity hazard |
| `toxicity_reservoir_half_life_days` | 1.5 | Per-drug toxicity decay |



### B.4 Per-Region Parameters (6 regions × N parameters)

Generated with key pattern `{region}_{param}`:

| Parameter suffix | Description |
|------------------|-------------|
| `hospitalization_log_odds` | Regional hospitalisation modifier |
| `sepsis_onset_log_odds` | Regional sepsis onset modifier |
| `sepsis_mortality_multiplier` | Regional sepsis death modifier |
| `sepsis_recovery_log_odds` | Regional sepsis recovery modifier |
| `mortality_log_odds` | Regional background mortality modifier |
| `testing_multiplier` | Regional testing rate modifier |
| `antibiotic_initiation_log_odds` | Regional prescribing modifier |
| `drug_{drug}_availability` | Per-drug availability (0.0–1.0) |



### B.5 Per-Mechanism Parameters (40 mechanisms × N parameters)

Generated with key pattern for mechanism-specific parameters:

| Parameter pattern | Description |
|-------------------|-------------|
| `mech_{mechanism}_emergence_rate` | Base emergence rate |
| `mech_{mechanism}_reversion_rate` | Fitness-cost reversion rate |
| `mech_{mechanism}_enhancement_{drug_class}` | Enhancement multiplier (how much resistance the mechanism confers against each drug class) |
| `mech_{mechanism}_global_enhancement_multiplier` | Legacy global multiplier |



### B.6 Cross-Indexed Parameters



#### Potency matrix (42 × 52)

Key: `drug_{drug}_for_bacteria_{bacteria}_potency_when_no_r`

Values: 0.0 (no activity) to 1.0 (maximum activity). See Section 6.5.



#### HGT probability matrix (42 × 42)

Key: `hgt_prob_{source}_to_{target}`

Values: see Section 9.2 for the probability rules.



#### Resistance emergence rates (42 × 40)

Key: `bacteria_{bacteria}_mechanism_{mechanism}_emergence_rate`

Values: base rate × incidence band multiplier. See Section 7.3.



#### Demographic distribution (6 × 18)

Key: `demo_{region}_age_{start}_{end}`

Values: probability weight for each region-age combination.



### B.7 Syndrome Scoring Templates (10 syndromes)

See Section 6.2 for the full empiric scoring tables for all 10 syndromes.



### B.8 Drug Penetration Table (10 syndromes × 36 drug classes)

See Section 6.4 for the penetration values by syndrome and drug class.



### B.9 Cross-Resistance Groups

See Section 7.6 for per-bacteria cross-resistance groupings.



### B.10 Drug Introduction Dates

See Section 6.6 for the full table of 58 drug introduction time steps.



### B.11 Age-Specific Risk Templates

See Section 3.1 for the 7 risk templates and their multiplier vectors.



### B.12 Regional Drug Availability

See Section 6.6 for the per-region availability tables.

---



## Appendix C — Output Specification

The simulation produces a single large CSV file per run. This appendix describes the column structure so you can interpret the output data.



### C.1 Output File

Each simulation run produces a single CSV file:

```
amr_simulation_output_analysis_outputs/simulation_summary_NNNNNN.csv
```

where `NNNNNN` is a zero-padded run identifier.



### C.2 Row Structure

Each row represents one simulated day. The number of rows equals the total number of time steps (default 38,325).



### C.3 Column Categories



#### Scalar columns (per-timestep)

| Column | Type | Description |
|--------|------|-------------|
| `day` | int | Simulation day (0-indexed) |
| `year` | float | Calendar year (1930.0 + day/365.25) |
| `total_alive` | int | Living individuals |
| `total_infected` | int | Individuals with ≥1 active infection |
| `total_on_treatment` | int | Individuals receiving ≥1 antibiotic |
| `total_in_hospital` | int | Hospitalised individuals |
| `total_sepsis` | int | Individuals with sepsis |
| `total_died_infection` | int | Cumulative infection deaths |
| `total_died_sepsis` | int | Cumulative sepsis deaths |
| `total_died_background` | int | Cumulative background deaths |
| `total_died_toxicity` | int | Cumulative drug toxicity deaths |
| `total_new_infections` | int | New infections this day |
| `drug_stops_due_to_toxicity` | int | Drug courses stopped for toxicity this day |
| `total_carriers` | int | Individuals carrying ≥1 bacteria |
| `total_immunosuppressed` | int | Immunosuppressed individuals |
| `policy_name` | string | Active policy branch name |



#### Per-bacteria columns (~42 each)

| Pattern | Description |
|---------|-------------|
| `{bacteria}_infected` | Currently infected count |
| `{bacteria}_carriers` | Current carrier count |
| `{bacteria}_deaths` | Cumulative deaths from this species |
| `{bacteria}_new_infections` | New infections this day |
| `{bacteria}_new_infections_community` | From community acquisition |
| `{bacteria}_new_infections_hospital` | From hospital acquisition |
| `{bacteria}_new_infections_carrier` | From carrier-to-infection |
| `{bacteria}_sepsis` | Current sepsis count |



#### Per-drug columns (~58 each)

| Pattern | Description |
|---------|-------------|
| `{drug}_prescribed` | Courses initiated this day |
| `{drug}_active_treatments` | Currently on this drug |



#### Per-bacteria × per-drug columns (~2,436 each)

| Pattern | Description |
|---------|-------------|
| `{bacteria}_{drug}_activity_r` | Mean resistance (activity_r) |
| `{bacteria}_{drug}_majority_r` | Population-level resistance prevalence |



#### Per-region columns (~6 each)

| Pattern | Description |
|---------|-------------|
| `{region}_infected` | Regional infection count |
| `{region}_hospitalized` | Regional hospital count |
| `{region}_deaths` | Regional death count |



### C.4 Total Column Count

With 42 bacteria, 58 drugs, and 6 regions, the CSV contains approximately:

- ~16 scalar columns
- ~336 per-bacteria columns (42 × 8)
- ~116 per-drug columns (58 × 2)
- ~4,872 per-bacteria-per-drug columns (42 × 58 × 2)
- ~18 per-region columns (6 × 3)
- **Total: ~5,358 columns**



### C.5 Infection Journey Logs

When enabled, individual infection journeys are logged to the `infection_journeys/` directory as CSV files, capturing:

- Infection acquisition details
- Resistance profile at acquisition and over time
- Treatment episodes
- Clinical outcome (clearance, death, ongoing)
- Mechanism gains and losses

---

*This document describes the model as implemented in the Rust codebase. All variable names correspond to parameter keys used in `src/config.rs`.*
