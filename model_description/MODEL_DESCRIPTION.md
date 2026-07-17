# AMR Simulation — Technical Model Description



## Contents

1. [Overview](#1-overview)
2. [Population and Demographics](#2-population-and-demographics)
3. [Infection Acquisition](#3-infection-acquisition)
4. [Clinical Progression](#4-clinical-progression)
5. [Diagnostic Testing](#5-diagnostic-testing)
6. [Antibiotic Treatment](#6-antibiotic-treatment)
7. [Resistance Dynamics](#7-resistance-dynamics)
8. [Microbiome and Carriage](#8-microbiome-and-carriage)
9. [Horizontal Gene Transfer](#9-horizontal-gene-transfer-hgt)
10. [Mortality](#10-mortality)
11. [Potential future model uses](#11-potential-future-model-uses)
12. [Limitations](#12-limitations)
- [Appendix A — Bacteria, Drugs, Mechanisms and Enums](#appendix-a-bacteria-drugs-mechanisms-and-enums)
- [Appendix B — Parameter Reference](#appendix-b-parameter-reference)
  - [B.1 Global Parameters](#b1-global-parameters)
    - [B.2 Drug Properties](#b2-drug-properties)
    - [B.3 Bacteria Properties](#b3-bacteria-properties)
    - [B.4 Drug–Bacteria Potency Matrix](#b4-drugbacteria-potency-matrix)
    - [B.5 Regional Parameters](#b5-regional-parameters)
    - [B.6 Age-Dependent Parameters](#b6-age-dependent-parameters)
    - [B.7 Syndrome Parameters](#b7-syndrome-parameters)
    - [B.8 Clearance Parameters](#b8-clearance-parameters)
    - [B.9 Immunodeficiency, Sex, and Vaccination Parameters](#b9-immunodeficiency-sex-and-vaccination-parameters)
    - [B.10 Resistance Mechanisms](#b10-resistance-mechanisms)
    - [B.11 Horizontal Gene Transfer Matrix](#b11-horizontal-gene-transfer-matrix)
- [Appendix C — Output Specification](#appendix-c-output-specification)

---



## 1. Overview

### 1.1 Framework structure

We present a **framework** in the form of an **individual-based model** which simulates infection incidence, antibacterial usage, resistance emergence and risk of sepsis and death.  In the current configuration / set-up, we **simulate the lives of a repreentative sample the global population from 1930 to the present**, before the use of antibiotics, to the present day in **daily** time steps.  We typically simulate 3 million people alive at some time during this period.  The framework code is open-source and we encourage colleagues to use, and perhaps further develop, the framework to enable them to answer questions around how resistance has emerged over time (including through use of counterfactuals) and, critically, to predict the effect of various potential policies  on antibacterial drug resistance and infection related mortality.  Possible policies could relate to those aimed at minimizing drug resistance emergence and growth but also those that for which future drug resistance is a concern as a negative effect of a policy, such as those involving wider antibiotic use, and the framework should enable these negative effects to be quantified and compared with the positive effects.

As currently structured, the framework tracks 42 bacterial species, 62 antibiotics (grouped into 39 internal drug classes), and 46 resistance mechanisms. The population is distributed across 6 world regions (North America, Europe, Asia, Oceania, South America, Africa), each with distinct epidemiological, travel, hospitalisation, and healthcare profiles.

The framework aims to capture how bacterial resistance in one individual can affect risk in others in the same region (and accounts for movement between regions).  It does not, however, dynamically model infection spread per se; infection incidence in a given region in a given age, immunodeficiency and hospitalisation status is directly specified. 

The simulation advances in discrete daily steps. Each simulated day, every living person in the population is processed through a sequence of mechanistic rules. These rules govern the events that can happen to a person on any given day and include:

- Demographic processes (ageing, births, background mortality)
- Infection acquisition from community, hospital, or endogenous sources
- Infection progression, including potential development of sepsis
- Diagnostic testing — bacterial identification followed by antimicrobial susceptibility testing
- Antibiotic initiation, continuation, and cessation
- Resistance emergence via de novo mutation or horizontal gene transfer
- Mortality from infection, sepsis, or drug toxicity

**Implementation order within a simulated day.** In the live `apply_rules()` implementation, the ordering is: age/living-state checks and birth-cohort vaccination; reset daily microbiome flags; immunodeficiency onset/recovery; hospital admission or discharge; travel; sepsis onset or clearance for currently active infections; antibiotic stopping, pharmacokinetic decay, and new drug initiation/selection; toxicity and treatment-failure/restart checks; mortality and sepsis recovery; then the per-bacterium loop for new acquisition, microbiome acquisition/clearance, resistance emergence/reversion/floors, infection growth, symptoms, clearance, testing, and finally HGT. This order matters: sepsis onset is evaluated before new infections acquired later that same day, and new acquisition is evaluated separately for each bacterium rather than skipped because the person already has another active infection.

**Stochastic processes.** The model does not deterministically assign events such as infection, testing, treatment, or death. Instead, it calculates a *probability* for each event and then samples whether that event occurs. Repeated runs therefore produce slightly different trajectories.

**Calibration.** As laid out below, the framework contains thousands of parameters.  In the current configuration we have just a single value for each parameter which we argue provides a reasonable calibration to observed antibiotic consumption rates, resistance prevalence reported by named surveillance systems and burden studies (including WHO GLASS, ECDC EARS-Net, CDC AR Threats, and GRAM/GBD), and infection incidence data.  Nevertheless, we recognise that there is uncertainty over most of the parameter values, sometimes large uncertainty. Future users of the framework are likely to want to try to identify multiple sets of parameter values that produce an acceptable calibration in order to be able to express uncertainty over parameter values when comparing future policy options.


### 1.2 Document structure

The document is organised to follow the progression of an individual through the simulation:

| Section | Content |
|---------|---------------|
| **2. Population** | Who the simulated people are — age, sex, region, immune status |
| **3. Infection Acquisition** | How people catch bacteria (epidemiology — incidence, risk factors, hospital vs community) |
| **4. Clinical Progression** | What happens once infected — symptoms, syndromes, sepsis |
| **5. Diagnostic Testing** | When and how bacteria and resistance are identified |
| **6. Antibiotic Treatment** | How drugs are started, chosen, dosed, and stopped (empiric and targeted prescribing) |
| **7. Resistance Dynamics** | How bacteria become resistant and how resistance spreads (biology of AMR — mechanisms, selection pressure) |
| **8. Microbiome & Carriage** | Asymptomatic bacterial colonisation |
| **9. Horizontal Gene Transfer** | Bacteria sharing resistance genes (plasmid transfer between species) |
| **10. Mortality** | Case fatality rates, sepsis mortality |
| **11. Counterfactual Design and AMR-Attributable Burden** | How the resistance-suppressed counterfactual branch is constructed and used to estimate AMR-attributable deaths |
| **12. Limitations** | What the framework does not capture / caveats for interpretation |
| **Appendices** | Reference tables of all bacteria, drugs, parameters, and outputs |


Each section describes the modelling choices, their rationale, and the specific rules and parameter values. 

---


## 2. Population and Demographics

This section describes the virtual people in the model — who they are, where they live, and the health states they can be in. These characteristics determine each individual's risk of infection, treatment probability, and mortality. Since AMR outcomes differ substantially by age, geography, immune status, and care setting, these host attributes are required for realistic policy evaluation. The host layer is deliberately parsimonious: it represents the host differences most likely to matter for policy questions, rather than a full comorbidity-level clinical phenotyping framework.  Future users may wish to add further details (variables) for each individual.

### 2.1 Initialisation

The population is created at day 0 (representing the calendar year 1930). Each individual is assigned:

- **Age**: Drawn from a continuous demographic distribution that encodes both living individuals and future births. Negative age values at initialisation represent individuals who have not yet been born; they enter the simulation exactly when their age reaches zero. This is how the model handles births over the 105-year simulation period without needing a separate birth process.
- **Birth sex**: Male or female, assigned with equal probability.
- **Region**: Sampled from demographic weights reflecting the global population distribution.

The six regions and their approximate population shares determine the starting geographical distribution:

Where a table in this document includes a **Citation / source** column, that citation should usually be read as support for the presence, direction, or broad ordering of the modeled effect rather than as a claim that the exact tuned numeric value is taken directly from a single empirical estimate.

| Region | Population Share | Citation / source |
|--------|------------------|-------------------|
| Asia | ~55% | UN DESA Population Division, 2024 |
| Europe | ~15% | UN DESA Population Division, 2024 |
| Africa | ~12% | UN DESA Population Division, 2024 |
| North America | ~9% | UN DESA Population Division, 2024 |
| South America | ~6% | UN DESA Population Division, 2024 |
| Oceania | ~3% | UN DESA Population Division, 2024 |

These shares are intended as a coarse world-population partition for simulation purposes rather than a literal census reconstruction of any single year. Their ordering and approximate magnitudes are consistent with the United Nations World Population Prospects 2024, which provides official demographic estimates and projections across global regions and countries (UN DESA Population Division, 2024).

The regions differ in antibiotic availability, hospital capacity, testing rates, and the prevalence of specific pathogens. A person's region shapes nearly every aspect of their simulated clinical journey.

### 2.2 Ageing and age categories

Each day, every individual's age increments by one day. The model groups people into age categories that determine their risk profiles, reflecting the reality that risk of infection, presentation, and outcome differ substantially across the age spectrum.

**General age categories** (used for most risk calculations):

| Age Category | Age Range | Clinical relevance |
|--------------|-----------|-------------------|
| Infant | 0–1 year | Immature immune system, high infection susceptibility |
| Preschool | 1–5 years | Frequent respiratory and enteric infections |
| School Age | 5–18 years | Generally lowest infection risk |
| Young Adult | 18–50 years | Reference group for most risk calculations |
| Middle Age | 50–70 years | Increasing comorbidities |
| Elderly | 70+ years | Immunosenescence, highest mortality risk |



Within the broad 0–1 year category, the neonatal period (0–28 days) carries especially high infection risk and infection-attributable mortality. We keep neonates inside the general infant category for most risk calculations, but treat them separately where that distinction is most clinically important, namely sepsis onset and infection-related mortality (GBD 2019 Lower Respiratory Infections Collaborators, 2022).

Likewise, the broad 70+ category compresses clinically important heterogeneity within later life: the risk of poor infection outcomes rises further with increasing age, especially above age 80. We retain a single elderly category, while continuous age effects and the separate sepsis/mortality classification capture part of that additional late-life risk.


**Sepsis/mortality age categories** :

Neonates (0–28 days) have dramatically different case-fatality rates compared to older infants and the model uses a separate age classification for sepsis onset and infection-related mortality:

| Category | Age Range |
|----------|-----------|
| Neonatal | 0–28 days |
| Paediatric | 28 days–18 years |
| Young Adult | 18–50 years |
| Elderly | 50+ years |



### 2.3 Immunodeficiency

Since immunocompromised hosts — from HIV, chemotherapy, transplantation, advanced frailty — face substantially higher infection risk, treatment difficulty, and mortality (Fishman JA, 2007; Taplitz RA et al., 2018), the model captures this through two types of immunosuppression.  

**Temporary immunosuppression** represents medium-duration higher-risk episodes more compatible with prolonged corticosteroid exposure, chemotherapy or radiotherapy-related suppression, or other treatment-associated immunosuppression lasting weeks to months than with a brief viral illness or only a few post-operative days.   In the initial set-up we are presenting here people enter this state at a rate of `0.00005` per day and recover at `0.01` per day (average duration ~100 days).

Pregnancy is not currently represented as a separate maternal immunologic state within this immunosuppression framework, but future users of the framework may wish to incorporate this. 

**Chronic immunosuppression** represents long-term conditions like HIV/AIDS with low CD4 count, solid organ transplant, or autoimmune disease requiring ongoing immunosuppressive therapy. It develops at `0.00006` per day and recovers much more slowly at `0.0012` per day.

At simulation start, a configurable fraction (currently 5%) of the population are allocated to having immmunodeficiency (e.g. Martinson ML et al., 2024).

When a new immunodeficiency episode occurs in the model, the following age-band probabilities determine whether that episode is typed as **chronic** rather than **temporary**. They are therefore best read as a structural mapping from age to chronic-vs-temporary assignment, not as literal age-specific prevalence estimates of diagnosed immunodeficiency in the underlying population:

These immunodeficiency probabilities should be read as part of a **composite infection-vulnerability state**, not as literal prevalence estimates of formal immunodeficiency diagnoses. The model therefore aggregates classic immunodeficiency, transplant medicine, chemotherapy-related neutropenia, advanced HIV, frailty, and other clinically important causes of impaired host defence into one tractable state variable (Fishman JA, 2007; Taplitz RA et al., 2018). 

**How immunodeficiency affects the clinical risks:**

The table below summarises all the ways immunosuppression changes a person's trajectory through the model. Each effect has a real-world clinical rationale:

| Effect | Parameter | Value | Clinical effect |
|--------|-----------|-------|-----------------------------|
| Weaker direct empiric-start trigger in the absence of symptoms | `antibiotic_initiation_log_odds_immunodeficiency` | +0.2 | Immunodeficiency alone does not usually trigger treatment in the current model; symptoms, sepsis, and test results drive most starts |
| More diagnostic testing | `testing_immunosuppressed_multiplier` | ×2.5 | Clinicians investigate more aggressively in immunocompromised hosts |
| Higher sepsis risk | `log_odds_sepsis_onset_immunosuppressed` | +0.7 | ~2× higher daily risk of developing sepsis for a given bacteria level |
| Harder to recover from sepsis | `sepsis_recovery_log_odds_immunosuppressed` | −1.0 | ~2.7× lower odds of daily recovery, reflecting poor immune clearance |
| Higher mortality from sepsis | sepsis death log-odds | +1.5 | ~4.5× higher risk of dying during sepsis |
| Higher mortality from drug toxicity | toxicity death log-odds | +0.9 | ~2.5× higher risk — reflects drug interactions and organ dysfunction |
| Higher background mortality | `log_odds_mortality_immunosuppressed` | +0.916 | ~2.5× overall mortality uplift |


Clinically, severely immunocompromised patients often also receive broader empiric antibiotic cover, sometimes extending to agents such as carbapenems or aminoglycosides because of opportunistic pathogens, resistant organisms, repeated prior antibiotic exposure, and heavier healthcare contact. The current model does not encode a separate immunodeficiency-specific bonus for broad-spectrum drug selection. Instead, that real-world tendency is captured indirectly through the model's general empiric preference for broader-spectrum therapy, increased testing, higher sepsis risk, higher hospitalization exposure, and a small constrained prophylaxis pool for some immunocompromised hosts.

***Log-odds**.  Many sections of this document describe probabilities using log-odds (also called logit values), which is standard in medical statistics.  With our parameter names we refer to effects of covaraiates as "log_odds" when strictly these are log odds ratios:

- A probability of 50% corresponds to log-odds of 0.
- Negative log-odds mean the event is unlikely (log-odds of −2 ≈ 12% probability; −4 ≈ 2%).
- Positive log-odds mean the event is likely (log-odds of +2 ≈ 88%; +4 ≈ 98%).
- The model adds together multiple log-odds terms (e.g., a baseline term, an age term, a severity term) and then converts the total into a probability using the standard logistic function.

For example, the daily probability of starting antibiotics for a hospitalised septic immunocompromised patient with symptoms would include: *baseline log-odds (-5.5) + symptomatic infection (+6.0) + sepsis (+6.0) + hospitalised (+0.7) + immunodeficiency (+0.2) = +7.4*, yielding near-certainty because the acute clinical syndrome dominates the decision.

### 2.4 Hospitalisation

Given the concentration of resistant organisms, broad-spectrum antibiotic use, and vulnerable patients in hospital settings (Magill SS et al., 2018), the model simulates daily admission decisions, length of stay, and the elevated risks of nosocomial infection.

**Admission criteria.** Each day, the model calculates a probability of hospital admission for every individual using a logistic model. The key factors are:

| Factor | Log-odds contribution | Interpretation |
|--------|----------------------|---------------|
| Baseline (healthy person) | −10.4 | Very low daily risk (~0.003%) |
| Age | +0.02 per year | Older patients are progressively more likely to be admitted |
| Symptomatic infection (severity > 3.0) | +9.5 | Severe symptomatic infection materially increases admission probability even without sepsis |
| Sepsis | +13.0 | Sepsis is a strong driver of admission, producing near-immediate inpatient escalation in most cases |
| Regional healthcare access | varies | Reflects real-world differences in hospital capacity |


Independent of this baseline logistic admission process, starting a **hospital-managed antibiotic** also triggers inpatient management. In the current model this includes a broad set of parenteral hospital drugs plus a narrow oral reserve subset (`linezolid`, `tedizolid`) used as a proxy for infections that would usually be managed in hospital. This is a simplification: in real practice, some prolonged IV courses are delivered through outpatient parenteral antimicrobial therapy (OPAT), especially in higher-income settings and particularly for infections such as bone and joint disease that may require 4-6 weeks of IV treatment.


**Length of stay:** Once admitted, patients face a baseline discharge hazard of `0.28` per day (average stay ~3.6 days), with a hard maximum of 30 days. This baseline applies only to relatively uncomplicated admissions. Patients with active sepsis, any still-active infection above the model threshold, or a current **hospital-managed antibiotic** cannot be discharged; in the current model, septic patients therefore remain admitted until the sepsis episode has resolved, the infection has cleared below the discharge threshold, and any hospital-managed treatment course has finished. Real systems are less rigid because some patients complete part of a prolonged IV course via OPAT rather than remaining continuously hospitalised. The `0.28` figure should therefore be interpreted as an **effective all-cause discharge hazard for clinically stable inpatients**, not as a claim that sepsis admissions average only 3.6 days or that every real-world admission has the same geometric length-of-stay distribution.

**Regional healthcare access:**

Hospital access varies substantially across regions. The model uses regional modifiers that adjust the admission threshold:

| Region | Modifier | Interpretation | Citation / source |
|--------|----------|---------------|-------------------|
| Europe | +0.6 | Highest access (universal healthcare systems) | WHO UHC fact sheet, 2025; World Bank, `SH.MED.BEDS.ZS` |
| North America | +0.5 | Good access | WHO UHC fact sheet, 2025; World Bank, `SH.MED.BEDS.ZS` |
| Oceania | +0.4 | Good access in developed areas | WHO UHC fact sheet, 2025; World Bank, `SH.MED.BEDS.ZS` |
| Asia | 0.0 | Reference baseline (mixed access) | WHO UHC fact sheet, 2025; World Bank, `SH.MED.BEDS.ZS` |
| South America | −0.2 | Variable access | WHO UHC fact sheet, 2025; World Bank, `SH.MED.BEDS.ZS` |
| Africa | −0.5 | Most limited hospital capacity | WHO UHC fact sheet, 2025; World Bank, `SH.MED.BEDS.ZS` |



These modifiers should be read as a qualitative ordering of effective hospital access rather than literal estimates of admission probabilities. The ranking is consistent with broad cross-country differences in service coverage and infrastructure documented by WHO's universal health coverage monitoring framework and the World Bank's hospital-bed indicator, which show persistent between-country variation in effective access to care and inpatient capacity even as global service coverage has improved over time (World Health Organization UHC fact sheet, 2025; World Bank, `SH.MED.BEDS.ZS`).

Negative values mean patients are *less* likely to be admitted — not because they are less sick, but because hospital bed capacity is limited. This matters for AMR because patients who cannot access hospital care may not receive appropriate antibiotics or diagnostics, whereas international sepsis-care programmes have associated better structured in-hospital and ICU bundle delivery with lower hospital mortality (Evans L et al., 2021; Levy MM et al., 2010).




**Nosocomial (hospital-acquired) risks:**

The model captures pathogen-specific increased risk of infection acquisition for hospitalised people:

| Pathogen group | Current pattern in the live configuration | Clinical context |
|----------|-----------------------------------------|-----------------|
| Classic nosocomial opportunists | Strongly positive bacterium-specific hospital-acquisition terms | *A. baumannii*, *P. aeruginosa*, *S. maltophilia*, and related device-associated pathogens remain heavily hospital-enriched |
| Hospital-enriched Enterobacterales and enterococci | Moderate-to-strong positive bacterium-specific hospital-acquisition terms | Reflects line infections, postoperative infections, ICU outbreaks, and ward-level amplification |
| Mixed hospital/community organisms | Small positive or near-neutral tuned values depending on calibration | Captures organisms such as *S. aureus*, *E. coli*, and respiratory pathogens that remain important in both settings |
| Primarily community or STI pathogens | Neutral or only modestly positive tuned values | These organisms are still more often acquired in community transmission networks than from ward ecology |



Hospital patients also face higher baseline mortality (+0.262 log-odds, ~1.3×) and higher sepsis onset risk (+0.5 log-odds, ~1.6×), but they also have a higher probability of *recovering* from sepsis (+0.8 log-odds) because of access to intensive care. The background-mortality term here should be read as a residual inpatient case-mix / frailty adjustment, not as a hospital-acquired-infection term; HCAI pressure is modelled separately through the hospital-acquisition modifiers above.


### 2.5 Travel

Since international travel is a well-established vector for AMR importation — as illustrated by ESBL-producing *E. coli* acquired by European travellers in South and South-East Asia (Arcilla MS et al., 2017) — the framework needs a cross-region mixing mechanism.

The framework simulates this by giving each person a daily probability of travelling to another region (`0.0002` per day).  Travel frequency varies by region of origin, reflecting real-world patterns:

| Region | Travel multiplier | Rationale | 
|--------|------------------|-----------|
| Europe | ×3.5 | High international travel rates | 
| North America | ×3.0 | High travel, large business travel | 
| Oceania | ×2.5 | Geographic distance drives air travel | 
| Asia | ×1.5 | Rapidly growing travel volumes | 
| South America | ×0.8 | Moderate travel rates | 
| Africa | ×0.3 | Lowest international travel rates | 

Source:  UN Tourism, 2025; World Bank, `ST.INT.DPRT`; World Bank, `IS.AIR.PSGR` 


These multipliers are intended as a semi-quantitative ranking of cross-region mixing intensity, not as literal estimates of per-capita trip counts. The ordering is supported by broad regional patterns in UN Tourism's global and regional tourism dashboard and World Bank indicators for international tourism departures and air passenger volumes, which collectively show very high international mobility in Europe and North America, strong air-travel dependence in Oceania, rapid growth but substantial heterogeneity across Asia, intermediate volumes in South America, and lower outbound tourism and aviation intensity across much of Africa (UN Tourism, 2025; World Bank, `ST.INT.DPRT`; World Bank, `IS.AIR.PSGR`).

When a person travels, they are temporarily exposed to the infection risks and drug availability of the destination region. This can mean acquiring bacteria with resistance patterns typical of that region. In the current implementation, travel can only start while the person is at home and not hospitalised. A trip lasts 30 days, then the person returns to their home region. Destinations are sampled from a fixed origin-specific matrix and the home region is excluded from the candidate destinations for that trip.

**Age-specific travel risk.** Age-specific modifiers capture the higher risk of travel-related enteric diseases in younger adults — for example, young European adults travelling to endemic areas face elevated risk of *Salmonella enterica* serovar Typhi (+0.8 log-odds) and *Shigella* spp. (+0.5 log-odds), while *V. cholerae* risk is suppressed (−1.0) for these demographics unless visiting highly endemic zones. These modifiers reflect both behavioural factors (young adults engaging in higher-risk food/water exposure, street food consumption, water sports) and biological factors (previously unexposed immune systems). Older travellers and children show different risk profiles, with older adults generally having lower risk for enteric infections (likely due to prior exposure immunity, more conservative food/water behaviour, and sometimes reduced travel frequency) and very young children having elevated risk for multiple enteric pathogens but no risk for sexually transmitted infections acquired during travel.

---



## 3. Infection Acquisition

This section describes how people in the model acquire bacterial infections. A person can acquire bacteria from three main sources: the community (e.g., food, water, close contacts), the hospital environment (e.g., ventilators, catheters, other patients), or their own body (bacteria they are already carrying asymptomatically can flare into active infection). The model captures all three pathways, but it does so through a deliberately compressed acquisition architecture that preserves the main epidemiological distinctions needed for long-run AMR policy analysis rather than every route-specific exposure mechanism.


### 3.1 Community acquisition

Each day, the model evaluates each of the 42 bacterial species separately. A person can acquire a bacterium if they are not already actively infected with that same bacterium; having another active infection does not block acquisition. Co-infections are therefore possible. The model calculates a separate probability for each species using a logistic model (see Section 1.2) that combines several risk factors:

- **Baseline acquisition rate** for the specific bacterium — some bacteria (e.g., *E. coli*) cause infections far more frequently than others (e.g., *L. monocytogenes*)
- **Region** — infection rates vary by geography due to climate, sanitation, and population density
- **Age** — infants and the elderly are more susceptible to most infections; sexually transmitted infections peak in young adults
- **Immune status** — immunosuppressed individuals are at higher risk
- **Calendar era** — some infections have become more or less common over the decades


| Variable pattern | Function |
|------------------|-----------------|
| `bacteria_{name}_acquisition_log_odds` | How common this bacterium is overall |
| `{region}_bacteria_{name}_acquisition_log_odds` | Regional differences for this bacterium |
| `bacteria_{name}_log_odds_{age_category}` | Age-specific risk for this bacterium |
| `{bacteria}_{region}_log_odds_{age_category}` | Interaction between bacterium, region, and age |

#### Vaccination

Vaccination is implemented as a per-bacterium prevention layer that acts before infection or carriage is acquired. Each person therefore has a yes/no vaccination record (`vaccination_status` in the code) for every bacterium. Vaccination is assigned once, at cohort entry: at the start of the first active day, while a newborn's age is zero and before the daily age increment, the code checks the historically available vaccines and vaccinates that birth cohort with a probability determined by the vaccine's rollout progress at that calendar year. Once that record is set to `true`, it remains on for the rest of the simulation; there is currently no waning, revaccination, booster logic, or catch-up campaign.

The vaccine layer currently supports four bacterial vaccines:

| Vaccine | Target bacterium | Availability year |
| --- | --- | ---: |
| Pneumococcal | *Streptococcus pneumoniae* | 2000 |
| Meningococcal | *Neisseria meningitidis* | 1981 |
| Hib | *Haemophilus influenzae* | 1985 |
| Pertussis | *Bordetella pertussis* | 1948 |

Vaccination affects acquisition in two places:

- **Infection acquisition**: if an individual is vaccinated against bacterium *b*, the model adds `log_odds_vaccinated` for that bacterium to the infection-acquisition log-odds.
- **Microbiome / carriage acquisition**: the same log-odds adjustment is applied when modelling asymptomatic carriage acquisition.

The default fallback is `log_odds_vaccinated = -2.0`, corresponding to an odds multiplier of approximately $e^{-2} \approx 0.135$, so vaccination reduces acquisition odds by about 86.5% for bacteria that use the default. Vaccination does **not** directly modify bacterial growth after infection has started, symptom onset, sepsis progression, mortality, treatment choice, or transmission. It is therefore best interpreted as a static reduction in susceptibility rather than a full immune-history or herd-immunity model.

Under the current default parameter map, vaccination is active rather than dormant: each vaccine has a non-zero target birth-cohort coverage and a rollout duration. 



#### Age risk templates

Since age-specific infection risk varies by organism and syndrome, the model assigns each bacterium a **risk template** — a pattern describing how acquisition probability varies across six age bands. The multipliers below are applied to the baseline acquisition rate:

| Template | Typical use | 0–1y | 1–5y | 5–18y | 18–50y | 50–70y | 70+y | Clinical rationale |
|----------|------------|------|------|-------|--------|--------|------|--------------------|
| `respiratory` | *S. pneumoniae*, *H. influenzae* | 3.0 | 1.8 | 0.8 | 1.0 | 1.3 | 2.5 | U-shaped: infants and elderly most vulnerable |
| `gastrointestinal` | *Salmonella*, *Shigella* | 2.5 | 2.0 | 1.2 | 1.0 | 1.1 | 1.8 | Young children and elderly via food/water |
| `urogenital` | *E. coli* (UTI) | 1.2 | 0.8 | 0.9 | 1.0 | 1.4 | 2.2 | Female risk is bimodal: sexually active adult women and postmenopausal ages; coarse bands mainly capture the later-life rise |
| `skin_soft_tissue` | *S. aureus* | 1.5 | 1.3 | 1.1 | 1.0 | 1.2 | 1.8 | Moderate age variation |
| `bloodstream` | *P. aeruginosa* | 4.0 | 2.0 | 0.7 | 1.0 | 1.5 | 3.0 | Neonates and elderly at highest risk |
| `largely sexually_transmitted` | *N. gonorrhoeae*, *C. trachomatis* | 0.1 | 0.2 | 0.8 | 1.0 | 0.8 | 0.3 | Peaks in sexually active adults |
| `flat` | Default | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | Equal risk across all ages |



A multiplier of 3.0 for infants with the `respiratory` template means that an infant is three times as likely to acquire that bacterium compared to a young adult (the reference group at 1.0).

These community-acquisition templates should be read as structured relative-risk shapes rather than literal incidence-rate estimates for each age-region-organism cell. They preserve broad, globally observed patterns — such as the concentration of enteric disease in children (Troeger C et al., 2018), respiratory vulnerability at the extremes of age (GBD 2019 Lower Respiratory Infections Collaborators, 2022), and young-adult concentration of sexually transmitted infections (Rowley J et al., 2019) — while leaving exact organism-level burden to calibration of the bacterium-specific baseline and interaction terms against the model's target outputs.

### 3.2 Hospital acquisition

Since hospitals concentrate nosocomial pathogens — *Acinetobacter* and *Pseudomonas* on ventilators, *Staphylococcus* (including coagulase-negative staphylococci) and *Enterococcus* on central lines, *C. difficile* in antibiotic-exposed patients — the model uses separate hospital-specific acquisition parameters (`{bacteria}_log_odds_hospital_acquired`) for each species.

These hospital-acquisition terms are best interpreted as semi-quantitative rankings of nosocomial exposure pressure rather than directly sourced empirical estimates. That is consistent with the way global AMR surveillance systems aggregate routine clinical microbiology data: they show that healthcare-associated pathogen mixes differ systematically from community mixes, but with large between-country differences in sampling intensity, bed capacity, case mix, and laboratory coverage (World Health Organization GLASS dashboard, accessed 2026; World Health Organization Global antibiotic resistance surveillance report, 2025).

For pathogens whose transmission is overwhelmingly sexual or foodborne, the current configuration assumes effectively zero hospital acquisition risk.  So while organisms such as *C. trachomatis*, *N. gonorrhoeae*, *M. genitalium*, *T. pallidum*, and *Campylobacter* can still be diagnosed while a patient is in hospital, the model treats those episodes as infections acquired in the community rather than true nosocomial transmission.


### 3.3 Carrier-derived infection

Asymptomatic carriage (see Section 8) can give rise to endogenous infection when commensal organisms transition to an active infection site. This pathway is important for AMR because:

- The carried bacteria may already be resistant (having been selected by previous antibiotic courses)
- The person's resistance profile passes directly from carriage to infection and each resistance mechanism present in the microbiome compartment (`mechanism_microbiome`) is independently considered for transfer to the infection compartment (`mechanism_any`)

This pathway is governed by two parameters:

| Parameter | Value | Interpretation |
|-----------|-------|---------------|
| `carrier_resistance_inheritance_probability` | 0.50 | 50% chance that the carrier-inheritance pathway fires at all — when it does, individual mechanisms are copied from the microbiome to the infection compartment |
| `infection_from_microbiome_dampening` | 0.70 | Per-mechanism transfer probability: each mechanism in the microbiome has a 70% chance of being copied to the infection site, reflecting that not all colonising lineages successfully transition to the infection site |

The pathway-level inheritance probability is also scaled by the run-level `run_pathway_carrier_inheritance_multiplier`. If the pathway fires, only mechanisms already present in the person's carriage compartment are eligible, and those mechanisms are copied into `mechanism_any` before infection resistance is re-derived. In practice, this bridge is checked after a new infection has been established, provided the person is already carrying the same organism in `presence_microbiome`. It is distinct from sampling from the regional stored profile library during acquisition (Section 3.4), because it uses the person's own carriage profile rather than importing a profile from the surrounding community or hospital pool.



### 3.4 Resistance at acquisition

When a new active infection is established, the modelling framework determines whether (and, if so, what) resistance is present. The framework uses a set of linked mechanism-level pathways that reflect how resistance co-occurs in bacterial populations and how the simulation's policy branches interact with that process:

**Step 1 — Stored library of circulating resistance profiles and prevalence tracking**

After every simulated day, the model refreshes a **stored library of circulating resistance profiles** (`MechanismCache` in the code) containing up to 1000 complete resistance profiles for each combination of region × care setting (community / hospital) × bacterium. Each eligible active infection contributes its current majority-mechanism profile to that day's collection, including a zero mask for a fully susceptible majority strain. Every eligible active infection-day has an equal chance of entering its daily region × setting × bacterium reservoir, but a persistent infection can contribute again on successive days. The retained library is therefore infection-day weighted rather than a sample of unique people, incident infections, or unique bacterial lineages. Profiles are assigned according to the person's **current location** (region and hospital/community status) at the time of collection. Thread-local samples are combined using their full observation counts: an exact hypergeometric draw determines how many profiles to take from each local reservoir, followed by uniform sampling within each reservoir. Unequal processing-chunk sizes therefore do not bias which eligible infection-days are represented.

The stored profile library blends retained profiles with newly collected profiles each day. Every old profile has the same marginal probability of retention, implemented by stochastic rounding of the retained count followed by uniform removal; survival is therefore independent of profile position and resistance status. In the live configuration, both `community_profile_cache_retention` and `hospital_profile_cache_retention` remain 0.999. For an existing cache entry this corresponds to an expected marginal half-life of roughly 693 daily refreshes. It is a memory parameter for retained cache entries, not a biological half-life for host carriage, a bacterial clone, or a resistance mechanism. Hospital profile groups have one additional guard: if random retention and refill would remove the last resistant profile from a group that previously contained one, the refresh step restores one resistant example. That explicit hospital exception prevents chance extinction at the boundary and therefore does not follow neutral 0.999 survival in that particular circumstance.

Because the stored profile library includes both resistant and susceptible profiles, the **prevalence** of resistance to any given drug can be estimated directly from it: for each stored profile, the model checks whether at least one recorded mechanism would affect that drug, and the resistant fraction across all stored profiles gives the current prevalence estimate. The same stored profile library therefore supplies both the profiles sampled at acquisition and the prevalence estimates used later by prescribing logic (Section 6).

**Hospital enrichment levers.** Hospital-acquired active infections can be made more likely to draw already-observed resistant profiles in two ways. `hospital_resistance_prune_susceptible_percent` temporarily removes a chosen fraction of fully susceptible hospital profiles before the draw; resistant profiles are never removed. `hospital_resistance_concentration_factor` can then give extra weight to profiles carrying more resistance mechanisms. In the current live configuration, the concentration factors are neutral (`f = 1.0`) and hospital enrichment is mainly carried by temporary removal of some fully susceptible hospital profiles (default 50%, with higher organism-specific values for classic nosocomial pathogens). These levers change the probability of selecting already-observed profiles; they do not create artificial new combinations of mechanisms.

**Step 2 — Community dilution**

Clinical samples tend to over-represent resistant strains because susceptible infections are more likely to resolve quickly, generate less urgent microbiology, and be under-sampled in surveillance systems. To account for the fact that many community acquisitions are not drawn from the same resistant pool seen in clinical sampling, the model applies a per-bacterium **community resistance dilution factor** (`community_resistance_dilution_factor`). The term "dilution" is therefore a source-mixture shorthand, not a direct multiplier applied to resistance after a strain has been sampled.

Operationally, the factor is the probability that a community-acquired infection is treated as coming from the local human circulating profile library. If the random check fails, sampling from that library is skipped and the acquisition enters the exogenous/no-profile path. Static environmental and dynamic ratchet probabilities may add mechanisms on that path; the separate universal acquisition guard and later pathways such as same-person carriage inheritance can also add mechanisms independently. The factor therefore limits how much local human AMR is inherited directly from the profile library. It is not intended to force every non-human source to be susceptible: organisms whose food-chain, wastewater, or human-faecal source already carries substantial resistance are assigned higher human-reservoir values, supported by explicit exogenous probabilities, or both. This dilution step is only applied to community-acquired infections; hospital-acquired infections are always treated as drawing from a human-associated circulating pool.

The dilution factor is assigned by ecological category, reflecting the strength of each organism's link to the circulating human reservoir and to resistant exogenous reservoirs:

| Category | Dilution range | Example bacteria | Rationale |
|----------|---------------:|------------------|-----------|
| Environmental / waterborne | 0.03 | *A. baumannii*, *Pseudomonas*, *Stenotrophomonas*, *Burkholderia*, *Legionella*, *V. cholerae* | Community acquisition is dominated by environmental exposure rather than the local human circulating pool |
| Foodborne / mixed exogenous reservoir | 0.07–0.95 | *Campylobacter*, iNTS, *Yersinia*, *Listeria*, *S. Typhi*, *S. Paratyphi*, *Shigella* | The non-human or broad exogenous fraction remains important, but some organisms in this group are still largely human-maintained |
| Healthcare-associated | 0.05–0.40 | *C. difficile*, *Enterobacter*, *Citrobacter*, *Serratia*, *Morganella*, *Proteus*, *P. stuartii*, *S. epidermidis*, *K. pneumoniae*, *E. faecium*, *E. faecalis* | Resistance is amplified in hospitals, so the community pool is materially more susceptible |
| Endogenous flora / commensal | 0.60–0.80 | *E. coli*, *S. aureus*, *S. pneumoniae*, *B. fragilis*, *H. influenzae*, *H. pylori* | Community strains substantially reflect recent human ecology and carriage |
| Obligate human pathogen / STI | 1.00 | *N. gonorrhoeae*, *Chlamydia*, *Mycoplasma*, *Treponema*, MDR-TB, *Bordetella* | Human-only transmission means the community pool is the human circulating pool |

Per-bacteria values are calibrated within these ecological bands.

**Step 3 — Correlated mechanism profiles — weighted profile sampling**

Rather than deciding each mechanism independently from a separate prevalence estimate, the framework draws a **complete resistance profile** from the stored profile library described in Step 1, assigning all mechanisms recorded in that profile to the newly infected individual at the same time. The result is that newly acquired *E. coli*, for example, arrives with a resistance profile that mirrors an actual circulating strain — e.g., ESBL CTX-M together with fluoroquinolone resistance, as these co-occur on real plasmids (Partridge SR et al., 2018).

**Community-acquired infections** that pass the Step 2 dilution gate — meaning the model treats them as arising from the human circulating reservoir rather than from a wholly exogenous environmental exposure — draw a profile uniformly at random from the community profile group for the relevant region and bacterium. Because that stored group already includes susceptible as well as resistant profiles, this uniform draw reproduces the resistance prevalence already present in the community library after the ecological dilution step has done its work. In the very early warm-up phase, if the requested community profile group is still empty, the model temporarily falls back to the other stratum rather than inventing a synthetic profile.

**Hospital-acquired enriched sampling.** Hospital-acquired active infections can, when configured, draw from a hospital-weighted version of the stored local profile group. Early in a run, before the hospital profile group for that region and bacterium has accumulated a resistant example of its own, the model temporarily samples from the combined local community-plus-hospital set, so early hospital acquisitions can bootstrap from resistance already circulating nearby. Once resistant hospital profiles exist, sampling reverts to the hospital group alone. The model can also temporarily remove a chosen fraction of fully susceptible hospital profiles before the draw, while never removing resistant profiles. It can then apply the $f^k$ weighting described above, so profiles carrying more resistance mechanisms are more likely to be selected. If removing susceptible profiles would leave nothing to sample from, the unpruned set is used instead. The chosen profile is therefore always one that already exists in the stored library; this step changes the probability of selecting previously observed profiles, but does not create artificial new resistance combinations.

Carriage acquisition (Section 8.2) is related but not identical. Hospitalized carriage acquisitions can use the same $f^k$ weighted sampling when `hospital_resistance_concentration_factor > 1`, but they do not use the additional susceptible-profile removal step reserved for active hospital infections. Carriage also has its own acquisition gate, `microbiome_resistance_multiplier_on_acquisition`, plus hospital and community dilution modifiers. A sampled carriage profile is recorded only in `mechanism_microbiome` and projected to `microbiome_r`; it does not create active-infection `mechanism_any` or `any_r` state. Transfer into a later same-organism infection occurs through the separate carrier-inheritance pathway.

If no profile exists in either relevant care-setting stratum for a given region and bacterium (early in the simulation, before enough infections have accumulated), no profile is assigned and the individual remains susceptible through the profile-sampling pathway.

**Counterfactual gating.** Profile sampling is gated by the `counterfactual_resistance_multiplier`. Before writing the sampled resistance profile into the individual, the model makes a random probability check against that multiplier (default 1.0; set to 0.0 in the counterfactual branch). If the check fails, none of the sampled profile's mechanisms are written, so the newly infected individual remains susceptible through this pathway. The same multiplier is also used in the de novo, microbiome-acquisition, and HGT pathways elsewhere in the rules engine. It is not, however, a single blanket switch around every acquisition-time rule: the exogenous floor path and the MDR-TB rifampicin rule are handled separately below.


**Step 4 — Exogenous reseeding and extinction-guard floors**

After profile sampling, the model applies two distinct floor-style preservation mechanisms.

First, when a community acquisition falls into the **non-human / exogenous fraction** (`!from_human_reservoir`), the code rolls each mechanism against an effective floor equal to:

$$\text{effective floor} = \max(\text{environmental floor}, \text{ratchet floor})$$

The **environmental floor** is configured only for selected bacteria-mechanism pairs and represents resistance maintained outside the local human circulating pool, for example in agricultural, food-chain, wastewater, or other exogenous reservoirs (see Section 7.7). The **ratchet floor** is dynamic and is derived from stepped historical peak values in the stored profile library (see Section 7.8). Despite the word “floor,” both values are Bernoulli mechanism-assignment probabilities on each exogenous acquisition. They provide reseeding pressure but do not mathematically guarantee that resistance prevalence across all infections remains above the same percentage.

Second, independently of the exogenous floor path, the model evaluates `calculate_resistance_floor(bacteria, drug, current_day)`. In the current code this is a **universal 1% acquisition reseeding guard**, not a per-organism ramped or prevalence-triggered floor. For every new infection and each eligible drug, if floors are enabled and the drug class has been introduced, the function returns `resistance_floor_default_level` (currently 0.01). A successful 1% check enters mechanism assignment; an applicable non-placeholder mechanism is then assigned with probability `mechanism_assignment_probability_on_any_r_gain` (currently 0.8), and `any_r` is derived from that mechanism. The rule operates whether current prevalence is low or high. It can reduce the chance that sparse resistance disappears while it remains represented, but it does not directly clamp prevalence to 1%.

**Causal correctness guard — the mechanism must be represented in the retained global cache.** The universal guard can assign a mechanism only if that mechanism is currently present in at least one stored profile for the bacterium across all regions and both community and hospital strata. Although the helper is named `mechanism_has_ever_emerged_globally`, it delegates to current retained-cache membership; there is no separate permanent registry of every historical emergence. The rule therefore cannot introduce a mechanism before it appears in the circulating cache, and it cannot resurrect one after it has disappeared from every retained profile.

The universal guard is active for all organisms once the relevant drug class exists, subject to current cache representation and mechanism applicability. Because it is attempted separately for eligible drugs and one mechanism can affect several drugs, its effect is not equivalent to setting a 1% organism-level resistance prevalence. On top of this guard, selected bacteria-mechanism pairs also receive explicit environmental reseeding probabilities, and mechanisms below the configured reversion-rate threshold can receive ratchet support through the exogenous pathway.

**S. aureus lineage enrichment.** One additional completion step is applied only to human-reservoir *Staphylococcus aureus* acquisitions. If the sampled profile already contains `TargetSitePbp2aMecA` (that is, the newly acquired infection is already MRSA-like), the model may add a small linked set of co-travelling mechanisms: `EnzymeBlaZ`, `TargetSiteErmB`, `EnzymeAacAph`, `MutationGyrAPrimary`, `MutationGyrAParCSecondary` (only if primary gyrA resistance is present), `ProtectionTetM`, and `ProtectionFusB`. Each addition is controlled by an explicit probability parameter and is allowed only if that mechanism has already appeared in the simulated circulating *S. aureus* population recorded in the global stored profile library (`MechanismCache` in the code).

This rule is currently restricted to *S. aureus* because it represents a specific MRSA-lineage phenomenon rather than a generic co-resistance shortcut. In *S. aureus*, `mecA`/PBP2a marks a well-defined transmitted lineage context (SCCmec/MRSA) where additional penicillinase, macrolide-lincosamide, aminoglycoside, fluoroquinolone, tetracycline, and fusidic-acid resistance determinants often co-travel in successful hospital and community clones. Other organisms still acquire correlated resistance through the sampled genotype profile itself, HGT, floors, and de novo selection, but the code does not currently contain an equally explicit lineage-completion rule for them.

The purpose of this step is not to claim that one newly infected person has spontaneously generated several new resistance determinants at once. Instead, it approximates the epidemiological reality that successful MRSA lineages often circulate as partially linked co-resistance packages, while a sparse genotype cache draw may capture only part of that lineage's resistance content. In other words, the enrichment step is a transmitted-lineage completion rule, not a multi-mutation de novo rule. Hospital-acquired *S. aureus* infections apply a modest multiplier to these conditional probabilities, reflecting the greater concentration of hospital-adapted MRSA lineages in wards and device-associated care settings.

**Step 5 — Converting recorded mechanisms into the summary resistance measure `any_r`**

After the resistance profile has been recorded (by profile sampling, floor enforcement, or both), the model calls `propagate_mechanism_resistance()` to translate the mechanisms recorded as present into the summary resistance measure `any_r` reported in outputs. For every drug, this function computes the **multiplicative susceptibility** across all active mechanisms:

$$\text{any\_r} = 1 - \prod_{m : \text{mechanism\_any}_m = \text{true}} (1 - e_m)$$

where $e_m$ is the enhancement multiplier for mechanism $m$ against that drug (see Section 7.2). At this stage the function recalculates resistance from the full recorded mechanism profile (`raise_only = false` in the code), so the derived value replaces whatever was there before. If later steps add additional mechanisms, the same function is called again in an upward-only update mode (`raise_only = true`), so `any_r` can rise when new mechanisms are added but is not reduced by that partial update.

**Step 6 — MDR-TB guaranteed rifampicin resistance**

*M. tuberculosis* is defined as multi-drug-resistant (MDR) by its rifampicin resistance — any isolate labelled MDR-TB will already carry a rifampicin-resistance mutation at the time of diagnosis (World Health Organization drug-resistant TB treatment guidelines, 2020). The model reflects this definitional requirement with a dedicated post-sampling step. For all MDR-TB acquisitions occurring in or after 1966 (when rifampicin was introduced), the current code records every resistance mechanism applicable to rifampicin for that bacterium — specifically the `rpoB` RNA-polymerase mutation — in both `mechanism_any` and `mechanism_majority`, regardless of what the earlier profile-sampling step assigned. `propagate_mechanism_resistance()` is then called again in upward-only mode so that the rifampicin `any_r` is updated to reflect this guaranteed mechanism. The acquisition provenance for rifampicin is stamped with `ResistanceAcquisitionType::AtInfectionTB` rather than `AtInfectionCommunity`, which allows the model to track this pathway separately in output statistics. The parameter `mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance` is currently used as an on/off gate: any positive value enables the hard-set rule, so the live default 0.90 behaves as "enabled" rather than as a 90% random probability check.

**Note**

After these six steps, the carrier-inheritance bridge described in Section 3.3 may apply additional mechanisms from `mechanism_microbiome` if the individual is already a carrier of the same organism. That step is separate from acquisition from the regional stored profile library: it is evaluated for any newly established same-organism infection with pre-existing carriage and is governed by `carrier_resistance_inheritance_probability`, `run_pathway_carrier_inheritance_multiplier`, and the per-mechanism `infection_from_microbiome_dampening` transfer probability.

#### Resistance appearance and expression summary

The model distinguishes three related but different records of resistance: mechanisms present anywhere in the active infection, mechanisms present in the dominant infecting strain, and mechanisms present in asymptomatic carriage. In the code these are `mechanism_any`, `mechanism_majority`, and `mechanism_microbiome`, respectively. `mechanism_any` is the source of `any_r`. `mechanism_majority` is what active infections contribute back into the stored profile library and it also affects HGT donor strength. `mechanism_microbiome` is the carriage compartment and is the source of `microbiome_r`. Mechanism acquisition, dominant-strain evolution, and drug-level resistance calculation are therefore separate operations.

| Pathway | When it runs | What changes in the model | Interpretation |
|---------|--------------|---------------|----------------|
| Human-reservoir profile sampling | New infection passes the human-reservoir draw | Records the sampled mechanism profile in `mechanism_any` and `mechanism_majority` | Acquisition of an already-circulating strain genotype |
| Hospital-enriched profile sampling | Hospital-acquired active infection with hospital enrichment configured | Records a hospital-weighted sampled profile in `mechanism_any` and `mechanism_majority` | Makes already-observed resistant ward profiles more likely to be selected without creating new combinations |
| Community dilution | Community-acquired infection before profile sampling | Determines whether the acquisition draws from the human circulating profile library at all | A failed draw does not itself add resistance; it sends the acquisition to the exogenous/no-profile path |
| Static environmental floors | Community acquisition not sampled from the human reservoir | Records one or more mechanisms as present in `mechanism_any` and `mechanism_majority` | Explicit exogenous reseeding for configured bacteria-mechanism pairs |
| Dynamic ratchet floors | Same exogenous path, for mechanisms below the configured reversion-rate threshold | Records one or more mechanisms as present in `mechanism_any` and `mechanism_majority` | Uses stepped historical cache peaks as exogenous reseeding probabilities; does not impose a hard prevalence bound |
| Universal acquisition reseeding guard | Every new infection, per introduced drug class | Can record one applicable mechanism currently represented in the global cache in `mechanism_any` and `mechanism_majority` | Repeated weak 1% entry check followed by mechanism assignment; not a direct 1% prevalence clamp |
| MDR-TB rifampicin rule | MDR-TB acquisition from 1966 onward when the gate is enabled | Sets rifampicin-applicable mechanisms in `mechanism_any` and `mechanism_majority` | Enforces the definitional rifampicin resistance of MDR-TB |
| Same-person microbiome inheritance | New infection when the person already carries the same organism | Copies eligible mechanisms from `mechanism_microbiome` into `mechanism_any` only | Carrier-to-infection bridge; copied mechanisms start as active-infection minority mechanisms |
| De novo infection emergence | Existing active infection under drug pressure | Adds newly emerged mechanisms to `mechanism_any` | Within-host emergence under selection; majority status can follow later under drug pressure |
| De novo microbiome emergence | Existing carriage with at least one positive-level active drug applicable to an absent mechanism | Gives each eligible mechanism one daily emergence roll and records successful transitions in `mechanism_microbiome` | Carriage-compartment emergence under binary drug pressure without multiplying the daily rate by the number of selecting drugs |
| Microbiome mechanism reversion | Existing carriage, evaluated independently for each present mechanism when no positive-level active drug is clinically applicable to that bacterium-mechanism pair | Removes successful transitions from `mechanism_microbiome` and recalculates `microbiome_r` | Effective loss of an unselected carriage mechanism; pressure selecting another mechanism does not block the transition |
| Majority-strain evolution | Existing active infection with at least one active drug applicable to a minority mechanism | Gives each eligible mechanism one daily promotion roll and records successful transitions in `mechanism_majority` | Converts a minority resistant subpopulation into the dominant/circulating strain without multiplying the daily rate by the number of selecting drugs |
| Infection-microbiome transfer | Same organism present in both compartments and they contain different mechanisms | Copies mechanisms so `mechanism_any` and `mechanism_microbiome` contain the combined profile | Within-host spillover between active infection and carriage |
| HGT | After the per-bacterium loop, between co-present bacteria | Adds transferable donor mechanisms to each recipient compartment where that organism is present: `mechanism_any` for active infection and `mechanism_microbiome` for carriage | Inter-species mechanism transfer, scaled by context and donor majority/minority status |
| Mechanism-to-drug propagation | After mechanism changes | Recalculates `any_r` and/or `microbiome_r` from the mechanisms recorded as present | Drug-level resistance is recalculated according to which mechanisms affect which drugs |
| Cross-resistance propagation | End of each per-bacterium update | Equalises `any_r` within configured drug groups | Phenotype bundling across related drugs; no mechanism bits are added or removed |

---



## 4. Clinical Progression

Once a person has acquired a bacterial infection, the model simulates the clinical course: which body site is affected, how the infection grows, whether it progresses to sepsis, and whether the body can clear it without treatment. The level of syndromic and host detail is chosen to support policy comparisons around antibiotic resistance rather than attempting to model the clinical course of each infection in high level detail.


### 4.1 Syndrome assignment

When a person develops an active infection, the model assigns an **anatomical syndrome**. This assignment is consequential because syndrome determines:

- **Empiric drug choice** (prescribing guidelines differ by site — see Section 6.2)
- **Drug penetration** (varies by tissue — see Section 6.4)
- **Replication rate** (e.g. bloodstream supports rapid growth; bone does not)
- **Sepsis and mortality risk** (bloodstream infections are far more dangerous than skin infections)

The 10 syndromes correspond to the major infectious disease presentations encountered in clinical microbiology:

| Syndrome | Index | Examples in clinical practice |
|----------|-------|------------------------------|
| UTI | 1 | Simple UTI, pyelonephritis — the most common bacterial infection |
| Skin/soft tissue | 2 | Cellulitis, wound infections, abscesses |
| Respiratory | 3 | Community-acquired and hospital-acquired pneumonia |
| Bloodstream | 4 | Bacteraemia, line-related infections |
| Intra-abdominal | 5 | Peritonitis, appendicitis, diverticulitis, cholangitis |
| CNS | 6 | Meningitis, brain abscess — drugs must cross the blood-brain barrier |
| Gastrointestinal | 7 | Gastroenteritis, food poisoning |
| Genital/pelvic | 8 | Sexually transmitted infections, pelvic inflammatory disease |
| Bone/joint | 9 | Osteomyelitis, septic arthritis — slow to resolve, needs prolonged treatment |
| Other | 10 | Device-related infections, undifferentiated febrile illness |



#### How syndromes affect disease behaviour

Each syndrome modifies two key aspects of the infection:

- **Syndrome antibiotic-initiation odds ratio** — how strongly symptoms at that site increase the odds of antibiotic initiation once the infection has become symptomatic. A patient with bacteraemia (×16) or meningitis (×14) is much more likely to start antibiotics than a patient with a simple UTI (×6).
- **Bacterial growth rate multiplier** — how fast the bacteria replicate at that body site. Bacteria in the bloodstream (×1.4) multiply faster than bacteria embedded in bone (×0.85).

| Syndrome | Antibiotic-initiation odds ratio | Growth multiplier | Clinical rationale |
|----------|-----------------------------|--------------------|-------------------|
| UTI | ×6.0 | ×1.0 | Common symptomatic outpatient presentation; dysuria and fever usually prompt treatment once symptoms emerge |
| Skin | ×6.0 | ×1.1 | Painful cellulitis, wound infection, and abscesses are commonly brought for treatment |
| Respiratory | ×10.0 | ×1.2 | Dyspnoea and fever drive rapid presentation |
| Bloodstream | ×16.0 | ×1.4 | Systemic toxicity and rapid incapacitation create a near-immediate treatment imperative |
| Intra-abdominal | ×10.0 | ×1.15 | Severe pain and systemic upset usually drive prompt review even if progression is not as explosive as bacteraemia |
| CNS | ×14.0 | ×1.3 | Severe headache, altered mental status, meningism, and neurological compromise should trigger urgent assessment |
| GI | ×8.0 | ×1.1 | Diarrhoea, vomiting, abdominal pain, and dehydration drive presentation |
| Genital | ×12.0 | ×0.9 | Symptomatic urethritis, cervicitis, PID, and genital ulcer disease often prompt treatment, but asymptomatic infections must still first cross the symptom threshold |
| Bone/joint | ×4.0 | ×0.85 | Important to treat once recognized, but many bone and joint infections are slower and less immediately explosive |
| Other | ×4.0 | ×1.0 | Catch-all for clinically recognized infections that usually merit treatment but lack a more urgent syndrome-specific cue |



### 4.2 Infection dynamics

To reflect the continuum from initial infection to sepsis, the model tracks a numerical **infection level** — an abstract measure of bacterial burden — that rises and falls over time, rather than using a binary infected/uninfected state.

- **Starting level**: When a person first acquires an infection, the bacterial load is low (`initial_infection_level` = 0.01).
- **Growth**: Each day, the bacteria multiply. The growth rate depends on the specific bacterium, the syndrome site (see above), and whether antibiotics are present and, if so, their activity level.
- **Symptom threshold**: When the infection level reaches `3.0` (`symptomatic_infection_level_threshold`), the person develops noticeable symptoms — fever, pain, cough, etc. — and begins seeking medical care. Below this threshold, they are infected but feel well enough that they do not present for assessment.

This mechanism matters for AMR because there is a window between acquiring an infection and becoming symptomatic during which bacteria can replicate before treatment is initiated.


### 4.3 Sepsis

Sepsis — the dysregulated host response to infection carrying high mortality (Singer M et al., 2016; Evans L et al., 2021) — is modelled as a distinct state that substantially increases both treatment urgency and death risk.

Each day, for each active infection that is not already septic, the model calculates the probability of progression to sepsis using a logistic model. The daily log-odds of sepsis onset combine organism virulence, bacterial burden, infection duration, infection site/syndrome, host vulnerability, region, and medical-care status. These terms are additive on the log-odds scale.

| Component | Parameter(s) | Current value / pattern | Interpretation |
|-----------|--------------|-------------------------|----------------|
| Per-bacterium baseline sepsis propensity | `<bacterium>_sepsis_baseline_log_odds`; fallback `sepsis_baseline_log_odds` | Fallback −12.0. Explicit organism values range from effectively non-septic sentinels such as *H. pylori* (−500.0) and MDR-TB (−38.0), through low-sepsis organisms, to high-risk invasive pathogens such as *P. aeruginosa* and *S. pyogenes* (both −5.7). | Captures organism-level invasive potential and virulence at a given bacterial burden. This is separate from the syndrome/site term below: it says how septic this organism tends to be, not where the infection is. |
| Bacterial burden | `<bacterium>_log_odds_sepsis_infection_level`; fallback `log_odds_sepsis_infection_level` | Default +0.93 per unit of infection level. Current organism-specific overrides include smaller level effects for *S. epidermidis* (+0.04) and *S. maltophilia* (+0.08). | Higher simulated bacterial burden increases daily sepsis-onset risk. |
| Duration of infection | `<bacterium>_log_odds_sepsis_infection_duration`; fallback `log_odds_sepsis_infection_duration` | Default +0.005 per day since acquisition. Current organism-specific overrides include *S. epidermidis* (+0.005) and *S. maltophilia* (+0.012). | Longer-standing infections gradually become more dangerous, especially if not brought under control. |
| Age and bacterium-age interaction | `sepsis_age_log_odds_baseline`, `sepsis_age_log_odds_neonatal`, `sepsis_age_log_odds_pediatric`, `sepsis_age_log_odds_young_adult`, `sepsis_age_log_odds_elderly`, plus `<bacterium>_<age_category>_sepsis_log_odds` | General age terms: baseline 0.0; neonatal (≤28 days) +1.10; pediatric (>28 days to 18 years) +0.18; young adult (>18 to 65 years) 0.0; elderly (>65 years) +0.69. Selected organisms add extra age-specific deltas, e.g. neonatal GBS, neonatal *E. coli*, pediatric pneumococcus/*H. influenzae*/*N. meningitidis*, elderly pneumococcus/*E. coli*/*Klebsiella*/*Pseudomonas*/*Acinetobacter*/VRE/*S. aureus*, and young-adult meningococcus/*S. aureus*. | Captures both general host vulnerability and pathogen-specific age patterns. The detailed organism-age terms are calibration terms rather than direct empirical case-fatality estimates. |
| Syndrome / infection site | `log_odds_syndrome_<id>_sepsis` | No-syndrome/unspecified (0): 0.0 by code path. UTI (1): −2.0; skin/soft tissue (2): −1.0; respiratory (3): 0.0; bloodstream (4): +1.5; intra-abdominal (5): +0.8; CNS (6): +1.2; gastrointestinal (7): −0.5; genital (8): −1.5; bone/joint (9): +0.5; other (10): defaults to 0.0 because no explicit sepsis value is configured. | Allows the same bacterium to carry different sepsis risk depending on the clinical syndrome: for example, *E. coli* bacteraemia is very different from uncomplicated lower UTI. |
| Region | `log_odds_sepsis_onset_region_<region>`; `Region::Home` hard-coded neutral | North America −0.5; Europe −0.6; Oceania −0.5; Asia −0.1; South America 0.0; Africa +0.1; Home 0.0. | Represents broad differences in early recognition, healthcare access, sanitation, treatment delay, and resource availability. These are calibrated regional modifiers, not direct measured sepsis rates. |
| Immunodeficiency | `log_odds_sepsis_onset_immunosuppressed` | +0.7 if the person currently has an immunodeficiency state; otherwise 0.0. | Immunocompromised people have higher daily risk of progressing from infection to sepsis. |
| Hospitalisation | `log_odds_sepsis_onset_hospitalized` | +0.5 if the person is currently hospitalised; otherwise 0.0. | Hospitalised people are a higher-risk case mix and often have more severe or device-associated infections. |
| Not currently under medical care | `log_odds_sepsis_onset_not_under_care` | +1.0 unless the person is taking any antibiotic, is hospitalised, or has a bacterium identified for a still-active infection; otherwise 0.0. This is a person-level supportive-care proxy and does not require the antibiotic or test to concern the particular bacterium whose sepsis risk is being evaluated. | Represents the nonspecific benefits of clinical attention, monitoring, oxygen, fluids, and other supportive care in addition to antimicrobial treatment. |
| Special *H. pylori* handling | Hard-coded onset override in `apply_rules()` plus `helicobacter_pylori_sepsis_baseline_log_odds` | If *H. pylori* is the only active infection, its daily sepsis-onset probability is forced to 0.0. Its baseline is also set to −500.0, so it has effectively no intrinsic sepsis propensity in normal model use. | Keeps *H. pylori* as a localized chronic disease/gastric-cancer organism rather than a classic invasive sepsis pathogen. Other active organisms in polymicrobial infection can still cause sepsis through their own calculations. |

The result is converted to a daily probability using the logistic function. These terms are calibrated model components designed to preserve plausible clinical ordering; the exact numeric values should not be read as direct empirical estimates of portable sepsis risk.


### 4.4 Natural clearance and microbiome dynamics

This part of the model contains two related but distinct processes:

- **Microbiome or carriage clearance**: `default_microbiome_clearance_probability_per_day` = 0.01 is the default daily chance of losing asymptomatic carriage from the microbiome reservoir, with bacteria-specific overrides for organisms that are known to persist much longer or clear more quickly. This baseline 1% per-day clearance corresponds to an average carriage duration of ~100 days before spontaneous clearance, though many pathogens are configured with longer duration times.
- **Duration penalty on carriage clearance**: `carriage_duration_log_odds_coefficient` = −0.01 per day, capped by `carriage_duration_max_log_odds_effect` = −2.0, applies to microbiome carriage. The rationale is that long-established colonization becomes harder to dislodge because organisms have had time to occupy a stable niche, form biofilms, and adapt to the host environment (Trampuz A et al., 2005), reflecting that a bacterium that has been carried for 200 days is substantially harder to clear spontaneously than one only acquired 1 day ago.
- **Drug-assisted carriage clearance**: `microbiome_clearance_probability_on_drug_treatment` = 0.80 is the probability that effective treatment also clears carriage once a drug-treated infection resolves. This means 80% of the time, when an infection is successfully treated with an antibiotic, the asymptomatic carriage compartment is also cleared; 20% of the time, the microbiome is not cleared and persists.
- **Microbiome acquisition on active drug therapy**: The model tracks a separate flag `microbiome_acquired_on_drug_today` that identifies when carriage is acquired while the individual is already on antibiotic treatment. This captures "opportunistic overgrowth" — bacteria that would normally be outcompeted but thrive in the antibiotic-disrupted microbiome (example: *Clostridioides difficile* overgrowth during broad-spectrum therapy). This is distinct from acquisition in the absence of treatment.

**Infection resolution** Infection level changes each day according to bacterial growth, host-driven suppression, and any active antibiotic effect. An infection resolves when the simulated bacterial level is driven down to a near-zero threshold in the rules engine, or when an immune-clearance event is triggered; this is not controlled by `default_microbiome_clearance_probability_per_day`.

This distinction matters for AMR because there can be a delay between infection acquisition and symptom-driven treatment. During that untreated interval, bacteria continue replicating, and resistant subclones can emerge or expand within the infecting population before antibiotics are started. Before treatment begins, the dominant process is therefore untreated growth and diversification rather than antibiotic selection.

---



## 5. Diagnostic Testing

Since the transition from empiric to targeted prescribing depends on laboratory turnaround — culture followed by AST, often taking days during which empiric therapy continues — the model simulates the decision to send a test, the delay in getting results, the possibility of laboratory errors, and the historical availability of testing technology. In modern laboratories, species identification from a blood-culture bottle that has flagged positive and some genotypic resistance calls can often be available within hours rather than days, but the current model collapses that heterogeneity into a single simplified turnaround parameter.

We do not attempt to reproduce the full heterogeneity of specimen quality, breakpoint revision, platform-specific AST performance, MIC levels, or local reporting conventions; instead we include the parts of the laboratory pathway most likely to alter prescribing and therefore policy-relevant resistance dynamics.


### 5.1 Historical introduction

Since diagnostic microbiology did not exist in 1930 the framework introduces testing capabilities at historically appropriate time points:

| Technology | Available from | ~ Calendar year | Clinical context |
|------------|---------------|-----------------|-----------------|
| **Bacterial culture** | Day 5,478 | ~1945 | Basic culture techniques became routine in the mid-20th century |
| **Antimicrobial susceptibility testing (AST)** | Day 9,131 | ~1955 | Standardised AST methods (e.g., disc diffusion) followed about a decade later (Bauer AW et al., 1966) |



Before these dates, all prescribing in the model is entirely empiric — clinicians have no laboratory information to guide drug choice. This represents the early antibiotic era, when penicillin was prescribed without knowing the susceptibility of the infecting organism.


### 5.2 The testing process

Once testing is available and ordered, the model simulates a laboratory workflow:

| Step | Parameter | Value | Interpretation |
|------|-----------|-------|-------------|
| **Bacterial identification delay** | `test_delay_days` | 3 days | Bacterial identification cannot be recorded until 3 days after infection acquisition. This is a simplified aggregate representation of specimen collection, culture, and identification. |
| **AST result delay** | `resistance_test_result_delay_days` | 2 days | Once AST is ordered following bacterial identification, its result panel remains pending for 2 days. Treatment selection cannot use that panel during the pending interval. |
| **Reporting error rate** | `test_r_error_probability` | 2% | AST results are wrong 2% of the time — the lab reports a resistant organism as susceptible or vice versa. This aims to reflect issues with breakpoint interpretation, contaminated samples, and technical failures. Error rates in disc-diffusion and gradient-strip AST methods are typically in the 1–5% range depending on organism and drug class (International Organization for Standardization, 2021; EUCAST, 2023) |



AST has three explicit states: not ordered, ordered but pending, and result ready. A ready fully susceptible panel contains zero resistance values, so panel readiness is stored separately and is not inferred from whether any reported resistance value is positive. Reporting error is sampled once when the complete panel becomes ready, and the reported panel remains a snapshot until the infection clears. Until bacterial identification and the subsequent AST delay have elapsed, treatment proceeds without patient-specific susceptibility information.


### 5.3 Testing criteria and rates

Since culture ordering rates vary by care setting and clinical urgency — from selective but still often clinically useful urine cultures in outpatient UTIs to higher, though still imperfect, blood-culture use in sepsis — the model captures these differences:

| Factor | Parameter | Value | Clinical meaning |
|--------|-----------|-------|-----------------|
| **Baseline culture rate** | `bacterial_testing_base_rate_per_day` | 15% per day | A symptomatic outpatient has a 15% daily chance of having a culture sent |
| **AST reflex rate** | `resistance_testing_base_rate_per_day` | 95% per day | Once a culture is positive, AST is almost always performed |
| **Sepsis** | `testing_sepsis_multiplier` | ×4.0 | Sepsis makes clinicians much more likely to send cultures, although recommended blood cultures are still not obtained in every real-world case; this changes testing probability, not laboratory turnaround time |
| **Immunosuppressed** | `testing_immunosuppressed_multiplier` | ×2.5 | Clinicians investigate more aggressively |
| **Hospitalised (culture)** | `bacterial_testing_hospital_multiplier` | ×8.0 | Hospital patients have far greater access to microbiology labs |
| **Hospitalised (AST)** | `resistance_testing_hospital_multiplier` | ×5.0 | Hospital settings are modelled as more likely to complete AST when resources are constrained, especially in earlier eras and lower-capacity settings; in modern high-income systems this difference may be small because cultured clinically relevant isolates often get AST regardless of referral source |



**Regional differences:** Laboratory capacity varies dramatically around the world. Many hospitals in sub-Saharan Africa lack the microbiological infrastructure that is routine in European hospitals (Jacobs J et al., 2019). The model captures this with regional testing multipliers:

| Region | Testing multiplier | Context | Citation / source |
|--------|-------------------|---------|-------------------|
| Europe | ×1.2 | Highest testing density | Jacobs J et al., 2019; World Health Organization GLASS dashboard, accessed 2026 |
| North America | ×1.1 | High infrastructure | Jacobs J et al., 2019; World Health Organization GLASS dashboard, accessed 2026 |
| Oceania | ×0.8 | Geographically dispersed | Jacobs J et al., 2019; World Health Organization GLASS dashboard, accessed 2026 |
| Asia | ×0.7 | Highly variable by country | Jacobs J et al., 2019; World Health Organization GLASS dashboard, accessed 2026 |
| South America | ×0.6 | Variable access | Jacobs J et al., 2019; World Health Organization GLASS dashboard, accessed 2026 |
| Africa | ×0.3 | Very limited lab infrastructure in many settings | Jacobs J et al., 2019; World Health Organization GLASS dashboard, accessed 2026 |



These regional differences have direct consequences for AMR: in settings where testing is rare, patients are more likely to continue on ineffective empiric therapy, creating selection pressure for resistance without the feedback loop of culture results to guide narrower prescribing.

As with the admission and travel modifiers above, these testing multipliers should be read as qualitative effective-capacity terms rather than literal claims about national culture rates. They combine laboratory availability, specimen transport, clinician ordering behaviour, turnaround reliability, and AST reporting infrastructure, which is the same bundle of constraints emphasized by WHO's GLASS laboratory-strengthening programme and reviews of district-level bacteriology capacity in resource-limited settings (Jacobs J et al., 2019; World Health Organization GLASS dashboard, accessed 2026).


### 5.4 Serious resistance marker drugs

The code lookup `serious_resistance_marker_drugs` defines one or more sentinel drug-resistance results for each organism. These markers are used when interpreting serious-resistance burden and in admission logic once a relevant resistant test result is known.

| Organism group or species | Marker drug |
|---------------------------|-------------|
| Enterobacterales and selected Gram-negative nosocomial pathogens: *E. coli*, *K. pneumoniae*, *Enterobacter*, *Citrobacter*, *Serratia*, *Morganella*, *Proteus*, *P. stuartii*, *A. baumannii*, *P. aeruginosa*, *B. cepacia* complex, *B. fragilis* | `meropenem` |
| *S. maltophilia* | `trim_sulf` |
| Staphylococci | `flucloxacillin` |
| Enterococci and *C. difficile* | `vancomycin` |
| *S. pneumoniae*, *S. agalactiae*, *T. pallidum* | `penicillin_g` |
| *S. pyogenes* | `erythromycin` |
| *H. influenzae* | `amoxicillin_clavulanate` |
| Macrolide-treated respiratory, atypical, STI, and cholera organisms: *M. catarrhalis*, *M. pneumoniae*, *L. pneumophila*, *B. pertussis*, *V. cholerae*, *C. trachomatis*, *M. genitalium* | `azithromycin` |
| Fluoroquinolone-treated enteric organisms: *C. jejuni*, *S. Typhi*, *S. Paratyphi A*, *Shigella* spp., *Y. enterocolitica* | `ciprofloxacin` |
| iNTS, *N. gonorrhoeae*, *N. meningitidis* | `ceftriaxone` |
| *L. monocytogenes* | `ampicillin` |
| *H. pylori* | `clarithromycin` |
| MDR *M. tuberculosis* | `rifampicin` |

---



## 6. Antibiotic Treatment

This section covers the entire antibiotic prescribing process as the model simulates it — from the decision to start an antibiotic, through drug selection and dosing, to stopping the course. Antibiotic use drives the selection pressure that causes resistance to emerge and spread.

The model aims to reproduce how antibiotics are prescribed in clinical practice — including imperfect decisions, regional variation in drug access, and the distinction between empiric therapy (before microbiology results are available) and targeted therapy (guided by culture and susceptibility results). Here especially, the intention is not to encode every clinical nuance of antimicrobial decision-making, but to represent the prescribing features most likely to change AMR trajectories under different policy environments.


### 6.1 Treatment initiation — deciding to start antibiotics

Each day, the model decides whether to start a new antibiotic course for each person, using a logistic model (see Section 1.2). The probability of starting antibiotics depends on the person's clinical state:

| Factor | Log-odds | Approximate effect | Clinical rationale |
|--------|----------|-------------------|-------------------|
| Baseline (no symptoms) | −5.5 | ~0.4% daily chance | Represents the fact that the person may not have presented for care, so represents background prescribing without a clear indication, including non-specific or precautionary use seen in ambulatory care (Fleming-Dutra KE et al., 2016) |
| Symptomatic infection | +6.0 | Jumps to ~62% | Once a patient has obvious symptoms (fever, pain, etc.), prescribing becomes likely |
| Sepsis | +6.0 | Near-certain | Sepsis is a medical emergency requiring immediate antibiotics |
| Immunodeficiency | +0.2 | weak positive effect in isolation | Immunodeficiency alone is not usually a stand-alone treatment indication; symptoms, sepsis, test confirmation, and the constrained prophylaxis pool drive most starts |
| No clinical indication | −1.05 | ~3× less likely | A protective factor: if there is no active bacterial infection prescribing is non-zero (due to symptoms that are not caused by bacteria infection |
| Lab-confirmed infection | +0.92 | ~2.5× more likely | Positive culture results prompt targeted therapy |
| Already on an antibiotic | +0.18 | ~1.2× more likely | People taking an antibiotic may be started on an additional agent (combination therapy) |


**Syndrome antibiotic-initiation odds ratios.** The parameters named `syndrome_<id>_initiation_multiplier` belong to this start-any-antibiotic step. They are odds ratios, not drug-choice scores. In the Rust code, the model looks across the person's active infectious syndromes, takes the largest configured syndrome multiplier, converts it to a log-odds increment with `ln(multiplier)`, and adds that to the treatment-initiation log-odds. Taking the maximum rather than multiplying all syndromes together prevents polymicrobial or multi-site infections from inflating the initiation probability simply because several syndrome records are present.

| Syndrome parameter | Live value | Log-odds contribution | Interpretation |
|--------------------|------------|----------------------|----------------|
| `syndrome_1_initiation_multiplier` (UTI) | 6.0 | +1.79 | Common symptomatic outpatient presentation |
| `syndrome_2_initiation_multiplier` (skin/soft tissue) | 6.0 | +1.79 | Cellulitis, wounds, and abscesses commonly prompt treatment |
| `syndrome_3_initiation_multiplier` (respiratory) | 10.0 | +2.30 | Pneumonia-like presentations often prompt antibiotics |
| `syndrome_4_initiation_multiplier` (bloodstream) | 16.0 | +2.77 | Bacteraemia creates a strong treatment imperative |
| `syndrome_5_initiation_multiplier` (intra-abdominal) | 10.0 | +2.30 | Severe abdominal infection usually requires prompt systemic treatment |
| `syndrome_6_initiation_multiplier` (CNS) | 14.0 | +2.64 | Meningitis or brain abscess should trigger urgent antibiotics |
| `syndrome_7_initiation_multiplier` (GI) | 8.0 | +2.08 | Severe diarrhoeal illness and dehydration drive presentation |
| `syndrome_8_initiation_multiplier` (genital/pelvic) | 12.0 | +2.48 | STI/PID syndromes often have guideline-driven treatment |
| `syndrome_9_initiation_multiplier` (bone/joint) | 4.0 | +1.39 | Important but often less immediately explosive |
| `syndrome_10_initiation_multiplier` (other) | 4.0 | +1.39 | Modest boost for recognised infectious presentations |

These syndrome-level initiation multipliers are separate from the drug-bacterium `initiation_multiplier` values described in Section 6.2. The former answers "how likely is any antibiotic to be started today?"; the latter answers "if a drug is being selected for an identified organism, how much does prescribing practice favour this drug-organism pairing?"


For a 65-year-old immunosuppressed inpatient with a symptomatic *E. coli* UTI and no lab results yet, the daily initiation log-odds would be roughly: -5.5 (baseline) + 6.0 (symptomatic) + 0.7 (hospitalised) + 0.2 (immunodeficiency) + ln(6.0) for the UTI syndrome antibiotic-initiation odds ratio = odds of about +3.2 before regional access modifiers, which is a high probability of starting antibiotics that day. If the same patient also has sepsis, a further +6.0 is added.

**Regional variation in antibiotic access:**

Not everyone who needs antibiotics can get them. The model captures the large global gradients in antibiotic access documented by consumption surveys (Klein EY et al., 2018):

| Region | Log-odds modifier | Effect on prescribing | Rationale |
|--------|------------------|----------------------|-----------|
| North America, Europe, Oceania | 0.0 | Reference | Good pharmaceutical access |
| Asia | −0.5 | ~38% reduction | Variable access across countries |
| South America | −0.8 | ~55% reduction | Limited access in some settings |
| Africa | −1.4 | ~75% reduction | Major access barriers in many countries |



These access barriers produce a well-recognised tension: in settings where antibiotics are hard to obtain, selection pressure for resistance is lower, but people die of treatable infections. The model captures both sides. The regional prescribing modifiers should therefore be read as **effective access-and-behaviour terms** combining healthcare access, affordability, dispensing practice, and care-seeking, rather than as pure pharmacy-supply measurements (Klein EY et al., 2018).


### 6.2 Drug selection — choosing which antibiotic to use

Once it is decided that a person is being started on an antibiotic it must choose *which* antibiotic. The choice depends on the information available at the time of prescribing and on the reason the start was triggered.

**Two modes of prescribing:**

1. **Empiric therapy** — the treater (or individual self-treating in some situations) has no lab results (yet) and must choose a drug on syndromic grounds. The model uses syndrome-specific scoring templates (see below) to approximate this guideline-anchored prescribing, and does not use hidden knowledge of the true causative bacterium or the bacterium-drug potency matrix at this stage. If there is no positive syndrome-specific empiric prescribing signal, the candidate drug is heavily penalised as having no usual empiric role for that syndrome, rather than being allowed to compete only because it is broad-spectrum or widely available. In the current implementation, a positive empiric signal means a syndrome-template score above the neutral baseline; a score of 1.0 is neutral/background and is not treated as positive. This is distinct from intrinsic bacterium-drug potency, which is used once the organism is identified.

2. **Targeted therapy** — lab results have identified the bacterium and its susceptibility profile. The model can now use intrinsic bacterium-drug potency and the simulated susceptibility profile to identify drugs with adequate activity. It strongly rewards narrow-spectrum choices at this stage (×5.0 bonus for narrow-spectrum drugs) and penalises unnecessary broad-spectrum use (×0.1 penalty), reflecting the principle of antibiotic de-escalation that sits at the core of antimicrobial stewardship guidance and is supported by hospital stewardship evidence from Europe and Asia (Barlam TF et al., 2016; Schuts EC et al., 2016; Lee CF et al., 2018).

**Drug scoring algorithm:**

For each candidate drug, the model calculates a score based on several factors. The final candidate scores are converted into weighted probabilities rather than choosing the single highest-scoring drug every time. The config parameter `drug_selection_temperature` is best read as a **prescribing-variability setting**. With the baseline value of **0.55**, the model usually selects one of the highest-scoring drugs, but it can still choose another clinically plausible option. Smaller values reduce this variability; larger values make prescribing more dispersed across available options.

The live implementation also records the **context** in which a course starts. This is important for interpreting drug-use outputs:

- A symptomatic infection without organism identification is recorded as `Empiric` and uses the active syndrome templates.
- A symptomatic, identified infection is recorded as `Targeted` and uses potency, susceptibility, and organism-specific prescribing preferences.
- An immunodeficient person without symptomatic infection enters a small `Prophylaxis` pool rather than the general empiric pool.
- A person with an active modelled bacterial infection that has not caused symptoms can still start treatment through the background syndrome-0 scoring path, but the course is recorded as `OtherActiveAsymptomaticModelledBacterialInfection`. The unobserved infection syndrome is not supplied to drug selection.
- A start with **no active modelled bacterial infection** and no prophylaxis indication is recorded as `OtherNoActiveModelledInfection` and uses syndrome 0 as a background prescribing template.

Syndrome 0 is therefore not an eleventh clinical infection syndrome. It is a background / off-model prescribing distribution for situations such as diagnostic uncertainty, viral-like illness, non-modelled infections, dental or procedural use, self-medication, and pharmacy supply without prescription. Its per-drug scores are relative weights, not probabilities. They shape drug choice only after the model has already sampled a start-any-antibiotic event, and they still pass through the usual era, regional availability, allergy, age contraindication, reserve-drug, and compartment/niche-drug restrictions. The current syndrome 0 template strongly favours common community oral agents and assigns trace scores to hospital-only, reserve, and novel agents so these do not compete meaningfully in no-active-infection starts.

**Historical empiric syndrome-era scores.** Empiric prescribing templates can also have time-varying syndrome-drug overrides of the form `syndrome_<id>_empiric_drug_<drug>_score_before_<YYYY>`. These are used only in the empiric branch, before organism identification. They allow the model to represent historical syndromic prescribing practice without giving the clinician hidden knowledge of the true organism. For example, genital/STI syndrome prescribing can favour sulfonamides, penicillin, tetracyclines, fluoroquinolones, cephalosporins, or azithromycin in the eras when those drugs were used empirically for urethritis/cervicitis/PID or suspected gonorrhoea. The same syntax can be used for other syndromes in future if their empiric guideline history needs to vary by era.

These syndrome-era scores are separate from drug-bacterium `initiation_multiplier_before_<YYYY>` values. Syndrome-era scores answer "what was commonly chosen empirically for this clinical syndrome before microbiology results?"; drug-bacterium multipliers answer "once a bacterium is identified, how strongly does prescribing favour this drug-organism pair?". Drug introduction dates and regional availability gates still apply to both pathways, so an era score cannot make a drug available before it exists in the simulation.

| Scoring factor | Empiric phase | Targeted phase | What it captures |
|---------------|---------------|----------------|-----------------|
| Syndrome-specific template score | Primary driver; may be static or era-specific; scores >1.0 count as positive empiric signals | Not used directly | How well this drug matches usual empiric prescribing for the infection site and calendar era |
| Spectrum width | Slight bonus (×0.85) for broad-spectrum | Strong penalty (×0.1) for broad-spectrum | Empiric phase favours broader coverage; targeted phase rewards spectrum minimisation |
| No usual empiric role for syndrome | Strong penalty (×0.001) when no positive syndrome-template signal is present | — | Prevents drugs from competing empirically on generic breadth or availability alone; this is not a claim about intrinsic activity against the hidden true organism |
| Insufficient intrinsic bacterium-drug activity | — | Near-zero score (×0.001) or exclusion when potency is below the targeted threshold | Once the organism is identified, potency and susceptibility guide whether the drug can work |
| Drug-bacterium initiation multiplier | — | Score multiplier, default 1.0 | Organism-specific prescribing preference for identified infections; distinct from potency and susceptibility |
| Narrow-spectrum bonus | — | ×5.0 | Reward de-escalation to targeted therapy |



**Restricted niche agents:** Some drugs are hard-blocked outside their clinically plausible niche, reflecting that such drugs are not used broadly despite their activity in vitro.

- **Retapamulin** is restricted to **skin/soft-tissue prescribing contexts only**. It is excluded from undifferentiated prophylaxis, no-syndrome empiric starts, sepsis, bloodstream infections, and all non-skin systemic syndromes. In targeted (organism-identified) therapy, it is only allowed when the pathogen is *Staphylococcus aureus* or *Streptococcus pyogenes* and the syndrome is skin/soft-tissue. This reflects retapamulin's approval and clinical use as a topical agent for impetigo and other superficial skin infections (Stevens DL et al., 2014).

- **Fusidic acid** is excluded from sepsis, bloodstream infection, and undifferentiated/no-syndrome starts, but in targeted therapy it is allowed for anti-staphylococcal **skin/soft-tissue infections** and also for **bone/joint infections** (reflecting its use in osteomyelitis when combined with other agents). This gives fusidic acid a somewhat broader niche than retapamulin while still preventing it from competing as a generic systemic therapy (Koning S et al., 2012).

- **Nitrofurantoin** is restricted to genuine lower-UTI (urinary tract infection) contexts. It is excluded from sepsis, bloodstream infection, systemic infections, and non-UTI syndromes. This reflects its renal concentration and low systemic bioavailability (Brunton LL et al., 2018).

- **Fosfomycin** is also kept within lower-UTI prescribing contexts in the current model, including situations where prior cultures or resistance history would make ESBL-active oral cover attractive (e.g., ESBL-producing *E. coli* UTI where oral options are limited). Both nitrofurantoin and fosfomycin remain excluded from sepsis, bloodstream infection, and undifferentiated/no-syndrome starts (Gupta K et al., 2011).

- **Furazolidone** is modeled as a **GI-local agent** rather than a urinary agent. It is only eligible in gastrointestinal syndromes and is excluded from sepsis, bloodstream infection, and all non-GI prescribing contexts. This keeps compartment-limited agents from competing as generic systemic therapy when their clinical role is anatomically narrow.

These restrictions substantially reduce empiric broad-spectrum therapy by preventing niche agents from out-competing more appropriate choices in syndromes where they are never used clinically.

**Immunodeficiency prophylaxis.** A separate constrained prophylaxis path is evaluated when a person has immunodeficiency but no symptomatic infection. It is not a general empiric-treatment pool and it does not use syndrome 0. The only candidate drugs are `trim_sulf` (score 1.2), `azithromycin` (3.5), `ciprofloxacin` (2.0), and `levofloxacin` (1.5). Reserve drugs are excluded, broad-spectrum drugs are penalised, lower-spectrum drugs receive a modest bonus, and fluoroquinolone prophylaxis is blocked for people younger than 18 years.

**Regional resistance surveillance:** If population-level resistance data shows that a drug class is failing frequently in the region, the model penalises empiric use of that drug — mimicking real-world guideline updates when local resistance rates exceed thresholds. For symptomatic empiric treatment, surveillance is restricted to bacteria in the model's clinical-presentation table that can plausibly cause at least one of the person's observed syndromes; it never uses the person's hidden infecting organism. With multiple observed syndromes, the plausible-organism sets are combined. The strongest applicable penalty is retained. Prophylaxis, asymptomatic active infection, and background/no-active-infection starts have no observed syndrome and therefore do not receive this syndrome-based resistance adjustment. Once a bacterium is identified but AST is still pending, surveillance instead uses the identified active organism or organisms. The same information boundary is used when regional resistance is considered as justification for empiric reserve-drug escalation.

| Local resistance rate | Empiric score penalty | Clinical parallel |
|----------------------|----------------------|------------------|
| >60% resistant | ×0.3 | Drug dropped from guidelines (e.g., ciprofloxacin for *E. coli* UTI in South-East Asia) |
| >45% resistant | ×0.5 | Drug used cautiously, alternatives preferred |
| >10% resistant | ×0.8 | Drug still used but with awareness of resistance risk |



The syndrome scoring tables below are therefore stylised prescribing-preference weights, not literal market-share estimates for each antibiotic. They are designed to preserve broad world-recognisable clinical tendencies such as narrower outpatient UTI therapy, broader empiric treatment for sepsis and intra-abdominal infection, and de-escalation after microbiology results, while allowing the realised prescribing mix to emerge from access constraints, testing availability, and resistance feedback.


#### Treatment cessation — stopping antibiotics

People stop antibiotics based on several factors.

These values are best interpreted as **daily probabilities of prematurely stopping treatment**, not as the inverse of total course length. They are dropout hazards calibrated so that most patients remain on therapy through a guideline-like treatment window, recognising that shorter courses may often be preferred (Llewelyn MJ et al., 2017).

| Scenario | Daily stop probability | Approximate implication | Real-world parallel |
|----------|----------------------|-------------------------|-------------------|
| Default course | 0.45% per day | About 94% of patients are still on treatment by day 14 | Standard course for many infections |
| No relevant active infection remains | 15% per day | Rapid discontinuation over the next few days once the presentation no longer seems to reflect an ongoing bacterial infection | Antibiotics stopped when the patient improves and ongoing bacterial infection is no longer thought likely, even if no bacterium was identified |
| Cholera / *E. coli* GI | 2.5% per day | Supports short-course therapy, with most patients still on treatment through about 3-5 days | Short courses per guidelines |
| *S. aureus* / *S. pneumoniae* | 1.5% per day | About 90% of patients are still on treatment by day 7 | Representative shorter courses for milder skin/soft-tissue or respiratory infection; more serious invasive infections can still require longer treatment |
| MDR-TB | 0.06% per day | About 90% of patients are still on treatment by 6 months before regional adherence modifiers | Prolonged anti-TB regimens |

A constant daily stop probability of 0.45% does not imply an average course of 14 days; rather, it means treatment is only rarely interrupted on any given day, so most courses extend to approximately two weeks.



#### Syndrome-specific empiric scoring templates

The exact empiric score templates are source-of-truth in `src/config.rs` under `empiric_syndrome_templates`. The current UTI template is shown below as an example:

| Drug | Current UTI empiric score |
|------|--------------------------:|
| meropenem | 90.0 |
| imipenem_c | 90.0 |
| ertapenem | 90.0 |
| meropenem_vaborbactam | 60.0 |
| amoxicillin_clavulanate | 30.0 |
| piperacillin_tazobactam | 30.0 |
| ceftazidime_avibactam | 30.0 |
| nitrofurantoin | 13.0 |
| fosfomycin | 10.0 |
| cephalexin | 8.0 |
| cefazolin | 7.0 |
| cefuroxime | 7.0 |
| cefixime | 7.0 |
| ciprofloxacin | 6.0 |
| levofloxacin | 6.0 |
| cefepime | 4.0 |
| ceftazidime | 4.0 |
| ceftriaxone | 3.0 |
| amoxicillin | 2.0 |
| ampicillin | 1.0 |
| trim_sulf | 1.0 |
| aztreonam_avibactam | 1.0 |
| colistin | 0.2 |
| vancomycin | 0.1 |
| linezolid | 0.1 |

These are scores before availability, age, syndrome restrictions, regional access, resistance surveillance, targeted-therapy modifiers, reserve-drug gates, toxicity avoidance, and the final weighted probabilistic draw are applied.

#### Species-specific and time-varying prescribing multipliers

The syndrome-level empiric tables above describe population-wide prescribing tendencies; they are agnostic to which organism is causing the infection because this is unknown to the treating clinician at the time. Once an organism is identified, targeted therapy uses a second parameter family named `drug_<drug>_for_bacteria_<bacterium>_initiation_multiplier`. These are **drug-choice score multipliers** for identified infections. They do not decide whether any antibiotic is started, and they are not measures of intrinsic drug activity.

Conceptually, bacterium-drug potency answers "can this drug work against this organism if there is no acquired resistance?", while the drug-bacterium initiation multiplier answers "given that this organism has been identified, how strongly would prescribing practice favour this drug?" This lets the model represent first-line guideline choices, fallback drugs, combination-regimen components, historical prescribing practice, and stewardship suppression separately from MIC-derived potency and simulated susceptibility results.

In the targeted-selection code, each candidate drug receives the largest applicable drug-bacterium initiation multiplier across the currently identified active infections. That maximum is multiplied into the candidate drug score after potency checks and organism-specific guideline heuristics. The maximum is used rather than multiplying all identified organisms together, so polymicrobial infections do not automatically create very large score inflation. The default for any drug-organism pair with no explicit entry is 1.0, meaning no modification.

Multipliers greater than 1 boost a drug's score for that organism; values less than 1 suppress it. Very small values such as 0.01 make a drug very unlikely to be chosen unless other logic has already hard-blocked it. A high initiation multiplier cannot by itself make an unavailable drug available, rescue a drug with insufficient intrinsic potency, override a known resistant susceptibility result, or bypass hard-coded niche/reserve-drug restrictions. It is a prescribing-preference term, not a biological efficacy term.

**Time-varying era overrides.** For many organisms the preferred drug changed substantially during the simulation period — not because of drug licensing (that is handled separately by `DRUG_INTRODUCTION_DATES`) but because of guideline shifts driven by accumulating clinical evidence or emerging resistance. The model encodes this using `_before_YYYY` suffix keys, such as `drug_ciprofloxacin_for_bacteria_neisseria_gonorrhoeae_initiation_multiplier_before_2007`. For a given drug-organism pair, the code uses the smallest cutoff year that is still later than the current simulation year; if no cutoff applies, the base multiplier is used. This allows a single config to describe a continuous temporal arc without adding time-step logic to the core prescribing loop.

The species-guideline shifts currently encoded are:

| Organism | Drug | Pre-era multiplier | Era | Base (post) multiplier | Historical rationale |
|---|---|---|---|---|---|
| *N. gonorrhoeae* | penicillin G | 14.0 | → 1987 | 2.0 | Penicillin sole first-line before FQs licenced |
| *N. gonorrhoeae* | doxycycline | 12.0 | → 1987 | 4.0 | Tetracyclines were a major alternative to penicillin G from the 1950s; doxycycline increasingly dominant from ~1967 |
| *N. gonorrhoeae* | TMP-SMX | 10.0 | → 1990 | 0.5 | Sulfonamides were the dominant GC treatment 1937–~1975 (near-universal resistance developed); TMP-SMX continued the selection pressure into the 1980s; proxy for the entire sulfonamide era |
| *N. gonorrhoeae* | ciprofloxacin | 0.5 | → 1987 | — | FQ not yet in gonorrhoea guidelines |
| *N. gonorrhoeae* | ciprofloxacin | 14.0 | 1987–2007 | 2.0 | Sole first-line in CDC and WHO guidance; resistance then ended use |
| *N. gonorrhoeae* | ofloxacin | 8.0 | → 2007 | 1.0 | Co-first-line FQ option in European/Asian guidelines 1990–2007; additional selection pressure alongside ciprofloxacin |
| *N. gonorrhoeae* | ceftriaxone | 2.0 | → 2007 | 12.0 | Adopted as first-line following FQ resistance |
| *S.* Typhi / *S.* Paratyphi A | chloramphenicol | 14.0 | → 1990 | 2.0 | Dominant first-line until FQ era |
| *S.* Typhi / *S.* Paratyphi A | ciprofloxacin | 0.5 → 14.0 | → 1990 / 1990–2010 | 2.0 | FQ first-line 1990–2010; declining after XDR emergence |
| *S.* Typhi / *S.* Paratyphi A | ceftriaxone / azithromycin | — | 2010+ | 10.0 / 8.0 | XDR typhoid era first-line |
| *S.* Typhi / *S.* Paratyphi A | ampicillin | 6.0 | → 2000 | 1.0 | Oral alternative to chloramphenicol until widespread resistance |
| *S.* Typhi / *S.* Paratyphi A | TMP-SMX | 7.0 | → 2000 | 1.0 | Widely used alternative until resistance spread |
| *Shigella* spp. | ampicillin + TMP-SMX | 7.0 | → 2000 | 1.0 | WHO first-line until multi-drug resistance |
| *Shigella* spp. | ciprofloxacin | 14.0 | 1990–2010 | 2.0 | Global first-line until FQ resistance |
| *Shigella* spp. | azithromycin | 2.0 → 10.0 | 1991–2010 / 2010+ | — | Preferred for FQ-resistant strains |
| *Campylobacter jejuni* | ciprofloxacin | 0.5 | → 1990 | — | Not established for *Campylobacter* before FQ era |
| *Campylobacter jejuni* | ciprofloxacin | 10.0 | 1990–2010 | 2.0 | Widely used for severe/traveller disease |
| *Campylobacter jejuni* | azithromycin | 3.0 | 1991–2010 | 5.0 | Preferred once FQ resistance prevalent |
| *E. faecalis* / *E. faecium* | vancomycin | 0.3 | → 1985 | 4.0 / 3.5 | Early vancomycin abandoned due to nephrotoxicity; reintroduced ~1985 for MRSA cover |
| *C. difficile* | metronidazole | 12.0 | 1977–2017 | 4.0 | IDSA dominant first-line; downgraded in 2017 IDSA/SHEA guidelines |
| *C. difficile* | vancomycin (oral) | 4.0 | 1977–2017 | 10.0 | Severe/refractory only before 2017; universal first-line after |
| *S. aureus* | ciprofloxacin | 10.0 | → 2000 | 2.0 | 1988–2000: heavily used empirically for MRSA before FQ resistance precluded guideline use; modern era restricted to MSSA indications |
| *M. genitalium* | doxycycline | 8.0 | → 1991 | 1.5 | Pre-PCR diagnosis era (before 1991): doxycycline was sole empiric first-line for all non-gonococcal urethritis; now used only for debulking |
| *T. pallidum* | erythromycin / azithromycin / clarithromycin | 3.0 / 3.5 / 1.5 | → 2010 | 0.75 / 0.5 / 0.25 | Historical oral alternatives for syphilis when penicillin was avoided or unavailable; modern use sharply reduced after widespread 23S-mediated macrolide resistance |
| *E. coli* | **nalidixic_acid** | 7.0 | **1963–1990** | 0.0 | First-line uncomplicated UTI 1963–1990; GyrA mutations seeded in *E. coli* reservoir for 24 years before ciprofloxacin was licensed |
| *Shigella* spp. | **nalidixic_acid** | 12.0 | **1963–1990** | 0.0 | First-line dysentery in LMIC 1963–1990; pre-existing GyrA profiles explain rapid FQ resistance emergence post-1990 |
| *Campylobacter jejuni* | **nalidixic_acid** | 8.0 | **1963–1990** | 0.0 | Gastroenteritis treatment; GyrA Thr86Ile selected before veterinary fluoroquinolone era |
| *S.* Typhi / Paratyphi A / iNTS | **nalidixic_acid** | 8.0 | **1963–1990** | 0.0 | First-line enteric fever in South/SE Asia and sub-Saharan Africa 1963–1990 |

*C. difficile* also has a `treatment_recognition_year` of 1977 — before the Bartlett et al. description of antibiotic-associated colitis, no antibiotic pressure is applied to this organism regardless of multiplier values.

For *M. genitalium*, the live calibration now treats macrolide resistance as primarily a 23S-rRNA mutation problem, with `erm_b` retained only as a low-probability MLS co-resistance fallback rather than a dominant pathway. This avoids forcing an implausibly strong mobile-element-style mechanism in an organism whose macrolide resistance is usually explained by chromosomal 23S substitutions.
The tetracycline side is handled similarly: `tet_m` is retained only as a small historical non-gonococcal-urethritis / doxycycline-era signal, not as a dominant modern resistance mechanism.

**Additional live prescribing heuristics.** Targeted therapy uses organism-specific first/second-line allow-lists for common pathogens such as *P. aeruginosa*, *S. aureus*, *S. maltophilia*, pneumococcus, *E. coli*, *K. pneumoniae*, Enterococcus spp., *A. baumannii*, AmpC-risk Enterobacterales, gonorrhoea, Shigella, MDR-TB, and others. Drugs outside those lists are heavily penalised once an organism is identified. High-potency drugs receive graded boosts, organism-drug initiation multipliers add disease-specific preferences, and regional resistance surveillance penalises choices with high local failure signals unless a susceptibility result supports them.

Reserve agents are gated separately. In targeted therapy, reserve drugs are strongly discouraged unless there is recent treatment failure or a severe hospitalised Gram-negative context. In empiric therapy, reserve use is stricter: recent failure or severe hospitalised Gram-negative disease must be present, and population resistance pressure can still block escalation. After this gate, the global reserve penalty is applied (`reserve_drug_score_penalty = 0.35`).

Drug choice also remembers recent drug toxicities. If toxicity stops a drug (toxicity as modelled is non-specific), that drug receives an approximately 1000-fold avoidance penalty for the configured avoidance window. When targeted therapy starts, ineffective existing drugs are stopped if they have been found to have no activity against any identified active infection. In non-severe empiric starts, a new empiric drug generally replaces existing empiric therapy rather than stacking indefinitely; people with sepsis can still retain combination therapy.

Starting hospital-managed therapy can itself trigger admission. Always-inpatient drugs and hospital-restricted reserve drugs force admission and block discharge while active. OPAT-eligible drugs such as selected cephalosporins, ertapenem, dalbavancin, and vancomycin trigger admission probabilistically (`opat_admission_probability = 0.70`) but do not block discharge once the patient is otherwise eligible.

### 6.3 Drug pharmacokinetics

The model uses a simplified pharmacokinetic representation in which each drug has a **half-life** and a **starting level** at administration. Since the mutant selection window — where sub-therapeutic concentrations select for resistance rather than clearing it — is a key driver of emergence (see Section 7.3), the shape of the drug-level decay matters for downstream resistance dynamics.

| Parameter | Default | What it represents |
|-----------|---------|-------------------|
| `drug_{name}_half_life_days` | Drug-specific | How quickly the drug is cleared from the body |
| `drug_{name}_initial_level` | 10.0 | Drug level immediately after dosing |
| `drug_{name}_double_dose_multiplier` | 2.0 | Level when a double dose is given |
| `drug_{name}_spectrum_breadth` | 3.0 | How broadly the drug disrupts the microbiome (higher = kills more bystander bacteria = more collateral damage) |



#### Selected drug half-lives

Half-lives vary enormously — from penicillin G (cleared within an hour, needing frequent dosing) to dalbavancin (which persists for two weeks, enabling single-dose therapy):

| Drug | Half-life (days) | Clinical note | Citation / source |
|------|-----------------|---------------|-------------------|
| colistin | 0.08 (~2 hours) | Very short — rapidly cleared by kidneys; poor CNS penetration (0.05) except in meningitis; reserved for pan-resistant Gram-negatives; toxicity includes nephrotoxicity and neurotoxicity (rare dosing-related; rare pre-existing) | Li J et al., 2006 |
| penicillin_g | 0.042 (~1 hour) | Very short — needs IV infusion or frequent dosing | Brunton LL et al., 2018 |
| ampicillin | 0.063 (~1.5 hours) | Short-acting penicillin | Brunton LL et al., 2018 |
| meropenem | 0.042 (~1 hour) | Short — given as IV infusion TDS | Brunton LL et al., 2018 |
| cefiderocol | 0.10 (~2.4 hours) | Short-acting novel siderophore cephalosporin | Wunderink RG et al., 2021 |
| ciprofloxacin | 0.17 (~4 hours) | Moderate — allows twice-daily oral dosing | Brunton LL et al., 2018 |
| nalidixic_acid | 0.08 (~2 hours) | Short — first-generation quinolone; rapid renal elimination; only renally/GI concentrated | Brunton LL et al., 2018 |
| linezolid | 0.21 (~5 hours) | Moderate | Brunton LL et al., 2018 |
| vancomycin | 0.25 (~6 hours) | Requires therapeutic drug monitoring | Rybak MJ et al., 2020 |
| sulfanilamide | 0.29 (~7 hours) | Historical agent | Brunton LL et al., 2018 |
| ceftriaxone | 0.33 (~8 hours) | Long enough for once-daily dosing | Brunton LL et al., 2018 |
| doxycycline | 0.75 (~18 hours) | Long — convenient once or twice-daily oral | Brunton LL et al., 2018 |
| azithromycin | 2.92 (~70 hours) | Very long tissue half-life — enables 3–5 day courses | Brunton LL et al., 2018 |
| dalbavancin | 14.0 (2 weeks) | Ultra-long — allows single-dose outpatient treatment | Dunne MW et al., 2016 |



#### Spectrum breadth — collateral damage to the microbiome

Since broad-spectrum agents exert collateral selection pressure on the commensal microbiome — creating ecological niches for resistant organisms — the model represents spectrum breadth in two related but distinct ways.

Specifically:

1. `spectrum_breadth` is a stewardship-facing drug property used when scoring treatment choices. In empiric therapy it favors broader agents when coverage is uncertain, while in targeted therapy it rewards de-escalation toward narrower agents once the pathogen is identified.
2. The longer ecological consequence is handled through each drug's `microbiome_disruption_log_odds`, which accumulates into a persistent `microbiome_disruption_level` reservoir. That reservoir decays over time rather than disappearing immediately when treatment stops, and it directly raises the log-odds of later microbiome acquisition events.

The model consequence is not limited to broader initial coverage. Broader therapy influences prescribing behavior up front, and microbiome disruption leaves a persistent ecological effect that can increase later carriage risk even after the course has finished.

Illustrative `spectrum_breadth` values:

| Drug | Breadth | Meaning |
|------|---------|---------|
| nitrofurantoin | 1.0 (Minimal) | Renally concentrated and rapidly metabolised; negligible gut microbiome disruption |
| penicillin_g | 2.0 (Narrow) | Minimal disruption to the microbiome |
| linezolid | 2.0 (Narrow) | Targets Gram-positives only |
| vancomycin | 2.5 (Narrow-medium) | Mainly Gram-positive spectrum |
| trim_sulf | 3.5 (Medium-broad) | Moderate disruption |
| azithromycin | 4.0 (Broad) | Significant microbiome disruption |
| ceftriaxone | 4.0 (Broad) | Major disruption; linked to *C. difficile* risk (Slimings C et al., 2021) |
| ciprofloxacin | 4.5 (Very broad) | Extensive gut microbiome disruption |
| meropenem | 5.0 (Very broad) | Maximum disruption — the broadest-spectrum agent |

Operationally, this means broad-spectrum therapy can affect the simulation in two downstream places: first by making a drug more attractive for empirical cover but less attractive for narrow targeted de-escalation, and second by increasing later colonization pressure through the microbiome-disruption reservoir that feeds carriage acquisition.



### 6.4 Drug penetration by syndrome

Since tissue penetration determines whether an antibiotic achieves adequate site concentrations, the model assigns penetration coefficients for each drug–syndrome pair. The pharmacokinetic distinctions most relevant to AMR involve:

- **CNS (meningitis):** The blood-brain barrier normally blocks most antibiotics, but bacterial meningitis causes substantial BBB inflammation that increases drug permeability. The penetration coefficients for CNS syndrome therefore reflect the *inflamed* BBB state rather than healthy-CNS values. Even so, drugs with very poor lipid solubility, large molecular weight, or active efflux transport (particularly aminoglycosides, polymyxins, and lipopeptides) remain inadequate at the site, while agents such as ceftriaxone, metronidazole, chloramphenicol, and linezolid achieve therapeutic CSF levels under these conditions.
- **Bone/joint:** Drugs must penetrate dense, poorly vascularised tissue. `rifamycins` and `fq` agents penetrate well; `ag_group1` and `ag_group2` do not.
- **Bloodstream:** By definition, any IV drug achieves full levels here (penetration = 1.0 for all drugs).

Penetration values range from 0.0 (no drug reaches the site) to 1.0 (full systemic concentration available). The table uses the model's internal drug-class keys; the full mapping is in Appendix A.3, but the keys used here are:

| Key | Meaning | Examples |
|-----|---------|----------|
| `fq` | Fluoroquinolones | ciprofloxacin, levofloxacin, moxifloxacin, ofloxacin |
| `sulf` | Sulfonamides | sulfanilamide, trim_sulf |
| `nitrofurans` | Nitrofurans | nitrofurantoin, furazolidone |
| `phosphonic_acids` | Phosphonic acids | fosfomycin |
| `mls` | Macrolides | erythromycin, azithromycin, clarithromycin |
| `lincosamides` | Lincosamides | clindamycin |
| `lipopeptides` | Lipopeptides | daptomycin |
| `oxa` | Oxazolidinones, not oxacillin | linezolid, tedizolid |
| `ag_group1` / `ag_group2` | Aminoglycosides | gentamicin/tobramycin; amikacin |
| `nitroimidazoles` | Nitroimidazoles | metronidazole |
| `carb_group1` / `carb_group2` | Carbapenems | ertapenem; meropenem/imipenem_c |
| `chl` | Chloramphenicol class | chloramphenicol |
| `poly` | Polymyxins | colistin |
| `macrocycles` | Macrocycles | fidaxomicin |
| `glyc` / `lipoglycopeptides` | Glycopeptides and lipoglycopeptides | vancomycin; teicoplanin/dalbavancin |
| `rifamycins` | Rifamycins | rifampicin |

| Syndrome | Best penetration | Poorest penetration |
|----------|-----------------|---------------------|
| UTI (1) | `fq`, `sulf`, `nitrofurans`, `phosphonic_acids` (up to 1.0) | `mls` (0.4), `lincosamides` (0.3), `lipopeptides` (0.1) |
| Skin (2) | `lipopeptides` (0.95), `fq` (0.9), `oxa` (0.9) | `nitrofurans` (0.2) |
| Respiratory (3) | `mls` (0.95), `fq` (0.95), `oxa` (0.9) | `lipopeptides` (0.0), `ag_group1`/`ag_group2` (0.4) |
| Bloodstream (4) | All 1.0 (reference compartment) | — |
| Intra-abdominal (5) | `nitroimidazoles` (0.9), `fq` (0.75), `carb_group1`/`carb_group2` (0.75) | `ag_group1`/`ag_group2` (0.3) |
| CNS (6) | `nitroimidazoles` (0.80), `oxa` (0.70), `chl` (0.70) | `ag_group1`/`ag_group2` (0.05), `poly` (0.05), `lipopeptides` (0.05) |
| GI (7) | `macrocycles` (1.0), `nitroimidazoles` (0.95), oral `glyc` (0.90) | IV `glyc`/`lipoglycopeptides` (0.35) |
| Genital (8) | `fq` (0.9), `nitroimidazoles` (0.8), `sulf` (0.8) | `ag_group1`/`ag_group2` (0.35) |
| Bone/joint (9) | `rifamycins` (0.80), `oxa` (0.75), `fq` (0.70) | `ag_group1`/`ag_group2` (0.25), `poly` (0.2) |



These penetration values directly affect treatment outcomes in the model: a drug with 0.05 penetration to the CNS will be nearly ineffective for meningitis even if the bacterium is fully susceptible.


### 6.5 Drug potency matrix

Since intrinsic susceptibility differs by organism, the model encodes a **potency matrix** — a 42×62 table (42 bacteria × 62 named drugs) where each cell represents the baseline activity of that drug against that bacterium when no acquired resistance is present. Resistance mechanisms are then applied on top of that baseline through the separate 39-class enhancement system described in Section 7.2.

Values range from 0.0 (no activity) to 1.0 (maximum activity). These potency values are based on published MIC data and clinical breakpoints. If an organism is intrinsically resistant to a drug (baseline potency $< 0.15$, the `minimal_potency_threshold_for_drug_selection` parameter), the model strictly prevents any *acquired* resistance mechanisms from being assigned to that organism-drug pair — for example, *Mycoplasma*, which lacks a cell wall, cannot acquire PBP mutations against penicillins.

Prescribing preference signals are represented through the `initiation_multiplier` parameter (e.g., fidaxomicin for *C. difficile* receives `initiation_multiplier = 1.05`).

A key modelling principle is that intrinsic resistance must be represented exclusively through low or zero potency.

| Drug class | Organisms | Basis for zeroing |
|---|---|---|
| Vancomycin | All Gram-negative bacteria | Glycopeptide molecule (~1450 Da) cannot penetrate the Gram-negative outer membrane; no acquired mechanism can confer susceptibility |
| Metronidazole | All aerobic and facultative organisms | Requires anaerobic nitroreductase activation; cytotoxic radical intermediates cannot form under aerobic conditions |
| Aztreonam | All Gram-positive organisms | Monobactam PBP3 target requires outer-membrane permeation absent in Gram-positive cell walls |
| Aminoglycosides | *C. difficile*, *B. fragilis* | Obligate anaerobes — AG uptake depends on oxygen-dependent active transport, which is completely abolished anaerobically |
| TMP-SMX | *B. fragilis* | Constitutively encoded insensitivity to sulfonamides (chromosomal folate pathway) |
| Nitrofurantoin | *S. maltophilia* | Intrinsic non-fermenter resistance; nitrofurantoin is never used for *Stenotrophomonas* infections |
| Penicillins, ceph 1–2G, carbapenems, macrolides, clindamycin, aztreonam | *S. maltophilia* | Chromosomally encoded L1 metallo-β-lactamase, L2 serine-β-lactamase, and SmeABC/SmeDEF efflux pumps render these drug classes intrinsically inactive; `potency_when_no_r` ≤ 0.05 across all affected classes (see Section 7.5) |

These zeroing values are encoded directly as `potency_when_no_r = 0.0` entries, covering every affected organism-drug pair explicitly.

Key examples:
- Meropenem vs *E. coli*: 0.95 (very high potency — carbapenem against susceptible Gram-negative)
- Vancomycin vs *S. aureus*: 0.95 (first-line MRSA therapy)
- Vancomycin vs *E. coli*: 0.0 (outer membrane blocks access — intrinsic, not acquired)
- Metronidazole vs *C. difficile*: 0.90 (obligate anaerobe — activated drug reaches target)
- Metronidazole vs *S. aureus*: 0.0 (aerobe — drug cannot be activated)
- Ceftriaxone vs *S. pneumoniae*: 0.95 (standard treatment for pneumococcal meningitis)
- Aztreonam vs *P. aeruginosa*: 0.80 (monobactam active against Gram-negatives including Pseudomonas)
- Aztreonam vs *S. aureus*: 0.0 (Gram-positive — outer-membrane route absent)



### 6.6 Drug availability by region and era

The model simulates antibiotics becoming available at their historical introduction dates. Before sulfanilamide was introduced in 1937, there were no antibiotics in the framework. Before penicillin G was introduced in 1942, the `pen` era had not yet begun. Before ciprofloxacin was introduced in 1987, the framework had no `fq` agents. This historical layering is essential for reproducing the sequential emergence of resistance over the 20th century.

**Regional availability:** Even after a drug is introduced globally, not all regions have equal access. Newer, more expensive drugs may be unavailable or rarely used in lower-income settings:

| Region | Access pattern |
|--------|---------------|
| North America | Full access to all drugs |
| Europe | Full access to all drugs |
| Asia | Most drugs available; limited access to tedizolid, ceftaroline (30%) |
| Oceania | Good access; limited novel agents (50%) |
| South America | Limited newer drugs (tedizolid 10%, linezolid 50%, carbapenems 60–70%) |
| Africa | Basic antibiotics available (80–100%); ceftriaxone 60%; vancomycin 30%; carbapenems 10–20%; most novel drugs 0–10% |



This has major implications for AMR: in Africa, where carbapenems are rarely available, carbapenem resistance may emerge more slowly — but when it does arrive (via travel or HGT), there are no last-resort drugs available to treat it.

These availability tiers should be interpreted as qualitative access strata rather than audited procurement shares. They summarize broad world patterns in which older essential antibiotics are much more widely available than newer reserve agents, and in which stewardship, financing, regulatory approval, supply-chain reliability, and laboratory support jointly determine whether a drug is realistically usable in practice (World Health Organization GLASS-AMU report, 2025; World Health Organization UHC fact sheet, 2025).



#### Drug introduction dates

The 62 antibiotics in the model span 88 years of pharmaceutical development. Antibiotic availability is gated by these introduction dates: before the date, no clinician can prescribe the drug; from the date forward, it becomes available (though subject to regional access constraints):

| Drug | Year | Day step | Drug | Year | Day step |
|------|------|----------|------|------|----------|
| sulfanilamide | 1937 | ~2,555 | ceftriaxone | 1984 | ~19,723 |
| penicillin_g | 1942 | ~4,383 | piperacillin_tazobactam | 1984 | ~19,723 |
| tetracycline | 1948 | ~6,575 | ceftazidime | 1985 | ~20,088 |
| chloramphenicol | 1949 | ~6,941 | imipenem_c | 1985 | ~20,088 |
| colistin | 1952 | ~8,036 | amoxicillin_clavulanate | 1985 | ~20,088 |
| erythromycin | 1952 | ~8,036 | aztreonam | 1986 | ~20,453 |
| nitrofurantoin | 1953 | ~8,401 | ciprofloxacin | 1987 | ~20,818 |
| furazolidone | 1955 | ~9,132 | teicoplanin | 1988 | ~21,183 |
| vancomycin | 1958 | ~10,227 | ampicillin_sulbactam | 1990 | ~21,913 |
| fosfomycin | 1959 | ~10,592 | clarithromycin | 1990 | ~21,913 |
| metronidazole | 1960 | ~10,957 | ofloxacin | 1990 | ~21,913 |
| ampicillin | 1961 | ~11,323 | azithromycin | 1991 | ~22,278 |
| fusidic_acid | 1962 | ~11,688 | cefepime | 1996 | ~24,106 |
| gentamicin | 1963 | ~12,053 | meropenem | 1996 | ~24,106 |
| rifampicin | 1966 | ~13,149 | levofloxacin | 1996 | ~24,106 |
| doxycycline | 1967 | ~13,514 | moxifloxacin | 1999 | ~25,201 |
| clindamycin | 1968 | ~13,879 | linezolid | 2000 | ~25,566 |
| trim_sulf | 1968 | ~13,879 | ertapenem | 2001 | ~25,931 |
| cephalexin | 1970 | ~14,610 | daptomycin | 2005 | ~27,394 |
| minocycline | 1971 | ~14,976 | ceftazidime_avibactam | 2006 | ~27,759 |
| amoxicillin | 1972 | ~15,341 | tigecycline | 2007 | ~28,124 |
| cefazolin | 1973 | ~15,706 | retapamulin | 2007 | ~28,124 |
| tobramycin | 1975 | ~16,441 | ceftaroline | 2010 | ~29,219 |
| amikacin | 1976 | ~16,806 | fidaxomicin | 2011 | ~29,584 |
| ticarcillin | 1977 | ~17,172 | tedizolid | 2014 | ~30,679 |
| cefuroxime | 1978 | ~17,537 | dalbavancin | 2014 | ~30,679 |
| piperacillin | 1981 | ~18,632 | ceftolozane_tazobactam | 2014 | ~30,679 |
| nalidixic_acid | 1963 | ~12,053 | meropenem_vaborbactam | 2018 | ~32,142 |
| quinu_dalfo | 1999 | ~25,201 | cefiderocol | 2019 | ~32,507 |
| flucloxacillin | 1970 | ~14,610 | (aztreonam_avibactam) | 2025 | ~34,675 |
| cefixime | 1989 | ~21,548 | — | — | — |

The canonical list also includes reference entries for less-modeled agents.

**Historical special case — colistin withdrawal and reintroduction.** Colistin was introduced in 1952 but experienced a clinical eclipse between approximately 1970 and 1995 due to severe nephrotoxicity and the availability of better-tolerated agents. The framework reflects this by reducing colistin availability to **5% of normal** during the withdrawal window (1970–1995 in calendar terms), modelling the situation where colistin was rarely prescribed outside select research and compassionate-use contexts. From 1995 onward, colistin availability returns to its regional baseline, reflecting its reintroduction as a last-resort agent for multi-drug-resistant Gram-negative bacteria (Li J et al., 2006).


### 6.7 Drug toxicity

Some drugs — particularly `ag_group1`/`ag_group2` agents (nephrotoxicity, ototoxicity) and `poly` (`colistin`, nephrotoxicity) — carry significant toxicity risks. The model simulates drug toxicity as a **reservoir** that accumulates with continued use and decays when the drug is stopped.

Toxicity can cause two outcomes:

**1. Drug discontinuation (sub-lethal toxicity):** When toxicity accumulates, the treating clinician may stop the drug. This is the more common outcome; for example, rising creatinine during gentamicin exposure may prompt a switch to a less nephrotoxic alternative. The model implements this as a **threshold check**: each day, the combined daily toxicity death risk (see below) is computed; if it exceeds a sub-lethal threshold, the drug with the highest toxicity reservoir is discontinued.

| Factor | Parameter | Value | Effect |
|--------|-----------|-------|--------|
| Sub-lethal threshold | `toxicity_discontinuation_threshold` | 0.00001 | When the daily toxicity death risk exceeds this level, the (most-toxic) drug is stopped |
| Recent toxicity avoidance | Avoidance penalty | ×0.001 (1000× penalty) | After stopping a drug for toxicity, it receives a strong prescribing penalty during the avoidance window |
| Avoidance window | `toxicity_discontinuation_avoidance_days` | 30 days | How long the prescriber avoids re-prescribing the toxicity-stopped drug |


**2. Drug-related death (lethal toxicity):** Rarely, severe drug toxicity can be fatal — for example, acute kidney injury from colistin leading to multiorgan failure.  The model uses a **multiplicative hazard** model: each drug has a per-unit daily hazard rate (typically in the 10⁻⁸ range), and the total risk is the sum of (drug level × drug-specific hazard) across all active drugs, multiplied by patient-specific vulnerability factors.

| Factor | Parameter | Value | Effect |
|--------|-----------|-------|--------|
| Drug-specific hazard | Per-drug hazard rate | ~10⁻⁸ (varies by drug) | Colistin and aminoglycosides carry the highest per-unit hazard |
| Infant vulnerability | `toxicity_age_multiplier_infant` | ×1.8 | Neonates more vulnerable to drug toxicity |
| Child vulnerability | `toxicity_age_multiplier_child` | ×1.2 | Moderate additional risk |
| Adult (reference) | `toxicity_age_multiplier_adult` | ×1.0 | Reference group |
| Elderly vulnerability | `toxicity_age_multiplier_elderly` | ×2.2 | Highest toxicity vulnerability — reduced renal clearance, polypharmacy |
| Immunosuppressed | `toxicity_immunosuppressed_multiplier` | ×2.5 | Immune compromise increases toxicity risk |
| Hospitalised | `toxicity_hospital_multiplier` | ×1.3 | Hospitalised patients are often sicker but also monitored |


### 6.8 Antibiotic infection prevention

Patients who are already receiving an effective antibiotic are partially protected against acquiring new infections (Bratzler DW et al., 2013). The model applies a 70% reduction in new infection risk for susceptible organisms when the patient is already on an active antibiotic (`antibiotic_infection_prevention_efficacy` = 0.7). This does *not* protect against resistant organisms — a crucial point, because it means patients on antibiotics are selectively more likely to acquire resistant infections relative to susceptible ones.


---


## 7. Resistance Dynamics

This section describes how the model represents the biology of resistance emergence and spread. The model tracks resistance at the level of individual **mechanisms**. This matters because the same phenotype (e.g., "carbapenem-resistant *K. pneumoniae*") can arise from very different mechanisms (KPC, NDM, OXA-48), each with different implications for treatment, spread, and even which novel drugs might still work.

**Mechanism-centred architecture.** For each person and bacterium, the model records which resistance mechanisms are present in the active infection and in any bacteria being carried / in the microbiome. Drug-level resistance is then calculated from those mechanisms. In the code these three records are `mechanism_any`, `mechanism_majority`, and `mechanism_microbiome`, while the derived summary resistance measures are `any_r` and `microbiome_r`. A single stored profile library (`MechanismCache` in the code) holds up to 1000 retained majority-mechanism profiles per region × care setting × bacterium for use in profile-based acquisition (Section 3.4). It is assembled from active infection-days rather than unique isolates or people. This library keeps separate hospital and community collections, retains profiles over time, preserves one resistant example at the hospital extinction boundary, derives cache prevalence directly from the stored profiles, and applies per-bacterium hospital enrichment when configured. It also maintains a `peak_mechanism_prevalence` table (Section 7.8) containing the highest annually sampled pooled-cache marginal prevalence for each bacterium-mechanism pair; stepped values derived from that table can support later exogenous reseeding.

The architecture is primarily mechanism-by-mechanism for **de novo** emergence under active drug pressure.


### 7.1 Resistance mechanisms

The model explicitly tracks **46** distinct resistance mechanisms. Each mechanism represents a specific biological pathway: an enzyme that destroys the drug, a mutation that changes the drug's target, a pump that ejects the drug from the cell, or a barrier that prevents the drug entering.

The table below lists every mechanism, the drugs it affects, and which bacterial groups can acquire it. It is intended as a reference table. The key point is that each mechanism has a defined scope: ESBL enzymes (rows 1–3) relate to `pen`, `c1_2g`, `c3g`, `c4g`, and related monobactam-active entries but not `carb_group1`/`carb_group2`, while KPC and NDM/VIM (rows 6–7) compromise the carbapenem classes as well, for example.


  | Mechanism | Variable name | Description | Explicit Drugs Affected | Bacterial Classes Affected |
  |-----------|--------------|-------------|-------------------------|----------------------------|
   | ESBL CTX-M | `esbl_ctx_m` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
   | ESBL TEM | `esbl_tem` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
   | ESBL SHV | `esbl_shv` | Extended-spectrum β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
   | AmpC CMY | `ampc_cmy` | Plasmid-mediated AmpC β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
   | AmpC DHA | `ampc_dha` | Plasmid-mediated AmpC β-lactamase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AmpC derepression | `mutation_ampc_derepression` | Chromosomal AmpC derepression (*ampC* regulatory mutations in SPACE organisms: *Enterobacter*, *Citrobacter*, *Serratia*, *Morganella*, *Providencia*) — chromosomal, non-transferable | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `aztreonam` | Enterobacterales |
   | KPC | `kpc` | *K. pneumoniae* carbapenemase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `ceftazidime_avibactam`, `meropenem_vaborbactam`, `aztreonam_avibactam`, `aztreonam`, `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
   | NDM/VIM | `ndm_vim` | Metallo-β-lactamases | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `ceftazidime_avibactam`, `meropenem_vaborbactam`, `aztreonam_avibactam`, `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
   | OXA-48 | `oxa_48` | Oxacillinase-type carbapenemase | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `ceftazidime_avibactam`, `aztreonam_avibactam`, `meropenem`, `imipenem_c`, `ertapenem`, `meropenem_vaborbactam` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | OXA-Acinetob. | `oxa_acinetobacter` | OXA-23/40/58 carbapenemases (A. baumannii) | `meropenem`, `imipenem_c`, `ertapenem`, `ceftazidime`, `cefepime`, `ceftazidime_avibactam` | Nonfermenters |
  | blaZ | `blaz` | Staphylococcal penicillinase | `penicillin_g`, `ampicillin`, `amoxicillin` | Staphylococci |
  | PBP2a/MecA | `pbp2a_meca` | PBP alteration (MRSA) | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `piperacillin_tazobactam`, `ampicillin_sulbactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `cefiderocol`, `ceftazidime_avibactam`, `meropenem_vaborbactam`, `aztreonam`, `aztreonam_avibactam`, `meropenem`, `imipenem_c`, `ertapenem` | Staphylococci, Helicobacter |
  | VanA | `vana` | High-level vancomycin resistance | `vancomycin`, `teicoplanin`, `dalbavancin` | Staphylococci, Streptococci, Helicobacter |
  | VanB | `vanb` | Variable-level vancomycin resistance | `vancomycin` | Staphylococci, Streptococci, Helicobacter |
  | GyrA (pri.) | `gyra_primary` | DNA gyrase mutation (step 1) | `ciprofloxacin`, `ofloxacin` | All |
  | GyrA + ParC | `gyra_parc` | Additional topoisomerase mutation | `ciprofloxacin`, `ofloxacin`, `levofloxacin`, `moxifloxacin` | All |
  | Qnr | `qnr` | Quinolone resistance protein | `ciprofloxacin`, `ofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | 16S rRMT | `16s_rrmt` | 16S rRNA methyltransferase | `gentamicin`, `tobramycin`, `amikacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | AAC/APH/ANT | `aac_aph` | Aminoglycoside-modifying enzymes | `gentamicin`, `tobramycin`, `amikacin`, `streptomycin`, `neomycin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Staphylococci, Streptococci |
  | ErmB | `ermb` | Erythromycin ribosome methylase | `erythromycin`, `azithromycin`, `clarithromycin`, `clindamycin`, `quinu_dalfo` | Staphylococci, Streptococci, Anaerobes, Fastidious, Helicobacter |
  | 23S rRNA | `23s_rrna` | 23S rRNA point mutation | `erythromycin`, `azithromycin`, `clarithromycin` | Helicobacter, Enteric Pathogens, Fastidious, Streptococci |
  | 23S rRNA (oxa) | `mutation_23s_rrna_oxazolidinone` | 23S rRNA domain V mutation conferring linezolid/tedizolid resistance — chromosomal, non-transferable | `linezolid`, `tedizolid` | Staphylococci, Streptococci |
  | Cfr | `cfr` | 23S rRNA methyltransferase | `linezolid`, `tedizolid`, `chloramphenicol`, `clindamycin`, `retapamulin` | Staphylococci, Streptococci, Anaerobes, Fastidious, Helicobacter |
  | CAT | `cat` | Chloramphenicol acetyltransferase | `chloramphenicol` | All |
  | MCR-1 | `mcr_1` | Mobilised colistin resistance | `colistin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | Polymyxin regulatory | `mutation_polymyxin_regulatory` | Chromosomal polymyxin resistance (*mgrB*, *pmrAB*, *phoPQ*, *lpx* mutations) — chromosomal, non-transferable | `colistin` | Enterobacterales, Nonfermenters |
  | AcrAB-TolC | `acrab_tolc` | Gram-negative efflux pump | `tetracycline`, `doxycycline`, `minocycline`, `tigecycline`, `chloramphenicol`, `ciprofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | MexXY-OprM | `mexxy_oprm` | Pseudomonas-specific efflux pump | `tetracycline`, `doxycycline`, `minocycline`, `gentamicin`, `tobramycin`, `amikacin`, `chloramphenicol`, `ciprofloxacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | Global eff. | `global_efflux` | Non-specific efflux upregulation | `tetracycline`, `doxycycline`, `minocycline`, `tigecycline`, `chloramphenicol`, `ciprofloxacin` | All |
  | TetA/B/C | `tet_abc` | Gram-negative tetracycline efflux | `tetracycline`, `doxycycline` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious |
  | TetM/TetO | `tetm` | Ribosomal protection | `tetracycline`, `doxycycline`, `minocycline` | All |
  | OmpK35/36 | `ompk35_36` | Outer membrane porin loss (Klebsiella) | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `ceftaroline`, `cefiderocol`, `aztreonam`, `aztreonam_avibactam`, `meropenem`, `imipenem_c`, `ertapenem`, `ciprofloxacin`, `levofloxacin`, `moxifloxacin`, `ofloxacin`, `gentamicin`, `tobramycin`, `amikacin` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | OprD | `oprd` | Outer membrane porin loss (Pseudomonas) | `meropenem`, `imipenem_c`, `ertapenem` | Enterobacterales, Nonfermenters, Enteric Pathogens, Fastidious, Anaerobes |
  | Global por. | `global_porin_loss` | Non-specific porin downregulation | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `ceftriaxone`, `ceftazidime`, `cefepime`, `ceftolozane_tazobactam`, `ceftaroline`, `cefiderocol`, `aztreonam`, `aztreonam_avibactam`, `meropenem`, `imipenem_c`, `ertapenem`, `ciprofloxacin`, `levofloxacin`, `moxifloxacin`, `ofloxacin`, `gentamicin`, `tobramycin`, `amikacin` | All |
  | Folate path | `folate_pathway` | Altered dihydrofolate reductase | `sulfanilamide`, `trim_sulf` | All |
  | Nitroreduct | `nitroreductase` | Nitroreductase loss | `metronidazole`, `nitrofurantoin`, `furazolidone` | Staphylococci, Streptococci, Enterobacterales, Enteric Pathogens, Anaerobes, Fastidious, Helicobacter |
  | FosA/FosB | `fosa` | Fosfomycin-modifying enzyme (FosA: Gram-negative; FosB: Gram-positive) | `fosfomycin` | Staphylococci, Streptococci, Enterobacterales, Nonfermenters, Enteric Pathogens |
  | MprF | `mprf` | Membrane charge modification | `daptomycin` | Staphylococci |
  | LiaFSR/Cls | `mutation_liafsr_cls` | Enterococcal daptomycin resistance via *liaFSR/cls* cell-envelope remodeling — chromosomal, non-transferable | `daptomycin` | Streptococci |
  | RpoB | `rpob` | RNA polymerase mutation | `fidaxomicin` | All |
  | FusB | `fusb` | Fusidic acid resistance determinant | `fusidic_a` | Staphylococci |
  | PBP mosaic | `mutation_pbp_mosaic` | Penicillin-binding protein mosaic mutations (PBP2x/2b/1a in pneumococcus, penA in gonococci, PBP3 in *H. influenzae*) — reduced β-lactam affinity | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `ceftazidime_avibactam`, `aztreonam`, `aztreonam_avibactam` | All |
  | mtrCDE efflux | `efflux_mtr_cde` | mtrCDE-type broad efflux pump (Neisseria, Haemophilus, Campylobacter CmeABC) | `erythromycin`, `azithromycin`, `clarithromycin`, `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `tetracycline`, `doxycycline`, `minocycline`, `chloramphenicol` | Fastidious, Enteric Pathogens |
   | Unknown | `as_yet_unknown` | Placeholder mechanism (dormant) | Evaluates `true` dynamically for all applied overrides | All (Calibration Placeholder) |



### 7.2 Mechanism–drug-class enhancement multipliers

Each mechanism a bacteria has **reduces** drug efficacy by a specific amount. The "enhancement multiplier" (0.0–1.0) represents by **how much** a drug's effectiveness is reduced:

- **0.0** = the mechanism has no effect on this drug (e.g., a tetracycline efflux pump does nothing against meropenem)
- **0.95** = the mechanism eliminates 95% of the drug's activity (e.g., NDM metallo-β-lactamase virtually destroys carbapenem efficacy)
- **1.0** = complete resistance - the mechanism being present means the drug has no effect (the same as if it had potency 0 (ie. maximal "intrinsic resistance"))

There are 46 mechanisms × 39 drug classes = 1,794 individual values. The table below shows the **global default** multiplier for the major mechanisms discussed most often in the text (used when a specific per-class value has not been directly specified):

These enhancement multipliers should be interpreted as qualitative within-model effect sizes rather than literal MIC shifts or breakpoint translations. Their role is to preserve the clinically familiar ordering in which carbapenemases, van genes, and key target-site alterations have very large effects, whereas efflux and permeability mechanisms are usually weaker on their own, while final residual drug activity still depends on baseline potency, site penetration, and combination with other mechanisms.

| Mechanism | Multiplier | Clinical interpretation |
|-----------|-----------|----------------------|
| NDM/VIM | 0.95 | Near-complete resistance — these metallo-β-lactamases destroy almost all β-lactams |
| VanA | 0.99 | Near-complete vancomycin resistance |
| KPC | 0.95 | Very high — KPC carbapenemases severely compromise carbapenems |
| PBP2a/MecA | 0.99 | Very high — defines MRSA; eliminates nearly all β-lactam activity |
| ESBL CTX-M | 0.80 | High — but β-lactamase inhibitor combinations retain partial activity |
| VanB | 0.99 | Very high vancomycin resistance |
| GyrA + ParC | 0.95 | High-level fluoroquinolone resistance (double mutation) |
| 16S rRMT | 0.95 | High-level aminoglycoside resistance |
| ESBL TEM | 0.60 | Moderate-high |
| OXA-48 | 0.60 | Moderate-high — but with variable carbapenem MICs |
| ErmB | 0.90 | MLS-B resistance (macrolides, lincosamides) |
| RpoB | 0.95 | Rifampicin resistance |
| ESBL SHV | 0.60 | Moderate-high |
| Cfr | 0.95 | Cross-resistance to oxazolidinones and phenicols |
| AmpC CMY/DHA | 0.70 | Moderate-high — overcomes β-lactamase inhibitors too |
| CAT | 0.90 | Chloramphenicol resistance |
| GyrA primary | 0.40 | First-step fluoroquinolone resistance (partial) |
| Folate pathway | 0.85 | Trimethoprim-sulfamethoxazole resistance |
| FusB | 0.70 | Fusidic acid resistance |
| FosA | 0.80 | Fosfomycin resistance |
| MCR-1 | 0.85 | Colistin resistance — critically important as colistin is the last resort |
| Nitroreductase | 0.70 | Nitrofurantoin resistance |
| OprD | 0.80 | Porin loss — carbapenem resistance (mainly in *Pseudomonas*) |
| MprF | 0.60 | Daptomycin resistance |
| OmpK35/36 | 0.80 | Porin loss — broad resistance in Enterobacterales |
| Qnr | 0.20 | Low-level quinolone resistance (facilitates further mutation) |
| Global porin loss | 0.20 | Broad, non-specific resistance via reduced permeability |
| MexXY-OprM | 0.30 | Efflux pump — aminoglycoside/FQ resistance in *Pseudomonas* |
| AcrAB-TolC | 0.30 | Gram-negative efflux — modest broad-spectrum resistance |
| Global efflux | 0.20 | Non-specific efflux — weakest single mechanism |
| AmpC derepression | 0.75 | Derepressed chromosomal AmpC — β-lactamase inhibitor combinations show reduced activity compared to plasmid AmpC |
| Polymyxin regulatory | 0.90 | Chromosomal colistin resistance (*mgrB*, *pmrAB*, *phoPQ*, *lpx*) — high-level, comparable to MCR-1 |
| LiaFSR/Cls | 0.75 | Enterococcal daptomycin resistance via *liaFSR/cls* cell-envelope remodeling |
| 23S rRNA (oxa) | 0.85 | Domain V 23S rRNA mutation conferring linezolid/tedizolid resistance in staphylococci and enterococci |
| MphA | 0.85 | Macrolide phosphotransferase (MphA) — confers high-level resistance (8–32× MIC increase) to azithromycin, erythromycin, and clarithromycin in *Shigella* and other Enterobacterales|



### 7.3 Resistance emergence

This subsection concerns **de novo resistance emergence during treatment**. In the current model, this pathway is only evaluated when a patient is actively exposed to antibiotics.

**Sub-therapeutic exposure and resistance emergence:**

Given the familiar mutant selection window framework (Drlica K et al., 2007), the model parameterises emergence probability as a function of drug concentration that peaks at intermediate exposure:

- **Very low drug levels:** Minimal selective pressure; resistant and susceptible subpopulations have little differential advantage.
- **Sub-therapeutic levels:** Susceptible bacteria are differentially suppressed while resistant mutants retain a survival advantage — the peak of the emergence curve.
- **Full therapeutic levels:** Both susceptible and resistant bacteria are strongly suppressed.

Within this framework, incomplete courses, poor adherence, and underdosing matter because they extend the time spent in the sub-therapeutic selection window.

**The emergence formula:**

For de novo emergence in an active infection under antibiotic exposure, the model calculates:

```
emergence_rate = mechanism_rate
               × run_pathway_infection_de_novo_multiplier
               × counterfactual_resistance_multiplier
               × (1 + bacteria_level_factor)
               × max_emergence_drug_factor
               × multi_drug_penalty_factor
```

| Factor | What it represents |
|--------|-------------------|------------------|
| `mechanism_rate` | How biologically likely this bacterium is to acquire this specific mechanism |
| `run_pathway_infection_de_novo_multiplier` / `counterfactual_resistance_multiplier` | Run-level and scenario-level scaling applied on top of the organism-specific baseline |
| `bacteria_level_factor` | Logarithmic scaling by bacterial load |
| `max_emergence_drug_factor` | Drug-exposure term based on the highest relevant site-level exposure window across active drugs |
| `multi_drug_penalty_factor` | Suppression when multiple relevant drugs are used together |



#### Organism-specific emergence calibration

De novo resistance emergence is parameterised directly at the **bacterium-mechanism** level. Each organism therefore has its own baseline emergence profile across the various mechanisms, with zeros retained for biologically implausible combinations (for example, *S. pyogenes* acquiring NDM-type carbapenemase remains disallowed).

For active infection, the daily hazard for a given mechanism is:

```
mechanism_emergence_rate = mechanism_rate
                          × run_pathway_infection_de_novo_multiplier
                          × counterfactual_resistance_multiplier
                          × (1 + bacteria_level_factor)
                          × max_emergence_drug_factor
                          × multi_drug_penalty_factor
```

The additional terms do distinct jobs:

| Term | Role in the model |
|------|---------------------------|
| `mechanism_rate` | Organism-specific baseline for that exact mechanism |
| `bacteria_level_factor` | Log-scaled increase with within-host bacterial burden, bounded by the configured organism maximum |
| `max_emergence_drug_factor` | Drug-exposure effect, highest at intermediate site concentrations and low at both minimal and fully suppressive exposure |
| `multi_drug_penalty_factor` | Suppression of emergence when two or more relevant drugs are active and the candidate mechanism covers only part of the regimen |
| `run_pathway_infection_de_novo_multiplier` / `counterfactual_resistance_multiplier` | Run-level or scenario-level scaling applied without changing the organism-specific baseline table |

**Minority-to-majority evolution.** A mechanism newly present in `mechanism_any` but not yet in `mechanism_majority` receives one daily promotion attempt with probability `majority_r_evolution_rate_per_day_when_drug_present` (default 0.18) whenever at least one drug with a positive current level is applicable to that bacterium-mechanism pair. Concurrent applicable drugs do not create additional attempts. A successful transition affects the dominant profile contributed to the circulating-strain cache and the mechanism's HGT donor strength; it does not change the already mechanism-derived `any_r` value.

Population realism is distributed across three parts of the model: organism-specific infection acquisition parameters determine how often each pathogen appears, organism-specific mechanism baselines determine which resistance pathways are plausible and how readily they arise, and within-host modifiers determine whether the ecological conditions for emergence are present on a given day.

**Microbiome pathway.** While a bacterium is carried, each absent applicable mechanism receives one daily emergence attempt whenever at least one active drug has a positive current level and is applicable to that bacterium-mechanism pair. Concurrent applicable drugs do not create additional attempts. The probability uses the same organism-mechanism baseline table and applies `run_pathway_microbiome_de_novo_multiplier` (multiplied by the compatibility multiplier `microbiome_de_novo_multiplier`, currently 1.0) and counterfactual scaling, but it does not use the infection-burden, concentration-window, or multidrug-penalty terms described above. Current drug pressure is therefore a binary trigger in this pathway.

These parameters should therefore be read as **effective emergence hazards** rather than literal mutation-rate measurements. They absorb biology, treatment ecology, and calibration targets jointly through explicit organism-mechanism parameterisation rather than through a separate incidence-band layer.

**Run-level pathway multipliers.** The calibrated run uses a small set of `run_pathway_*` multipliers to vary whole resistance pathways without rewriting organism-specific baselines. The current defaults are all 1.0. `run_pathway_infection_de_novo_multiplier` scales de novo emergence during active infection; `run_pathway_microbiome_de_novo_multiplier` scales de novo emergence in carriage; `run_pathway_hgt_multiplier` scales horizontal gene transfer; `run_pathway_reversion_rate_multiplier` scales mechanism loss in infection and carriage; `run_pathway_microbiome_acquisition_multiplier` gates resistance-profile inheritance when new carriage is acquired; `run_pathway_carrier_inheritance_multiplier` scales microbiome-to-infection inheritance; `run_pathway_community_dilution_multiplier` scales the community dilution term applied to acquisition; and `run_pathway_microbiome_disruption_multiplier` scales the antibiotic-driven microbiome disruption reservoir. Separately, policy branches can set `counterfactual_resistance_multiplier` to 0.0, which blocks profile-sampled resistance, de novo emergence, microbiome-acquisition seeding, and HGT transfer in the counterfactual branch. Acquisition-time floors and the MDR-TB rifampicin rule are separate pathways described in Section 3.4.




### 7.4 Resistance reversion and fitness costs

Since fitness costs mean resistant bacteria often replicate more slowly than susceptible competitors in the absence of antibiotic pressure (Andersson DI & Hughes D, 2010), resistance can gradually decline when drug use is reduced. The model assigns each mechanism a daily **reversion rate**, used as an effective probability of removing the mechanism from carriage or demoting it from the dominant active-infection strain when no selecting antibiotic is present. Higher rates represent faster effective turnover; lower rates represent greater persistence. This is a compressed population/within-host approximation rather than a literal genetic back-mutation rate. In the rules engine, per-mechanism reversion rates are scaled by `run_pathway_reversion_rate_multiplier` (default 1.0), with additional community-specific reversion multipliers for some bacteria.

Reversion operates in **both** compartments, but not in exactly the same way. In the active infection, fitness-cost loss removes a mechanism from `mechanism_majority`, so it no longer contributes to majority-strain surveillance or seeding of newly acquired infections; `mechanism_any` is retained for the currently infected individual. In the microbiome compartment, reversion removes the mechanism from `mechanism_microbiome`, after which `microbiome_r` is re-derived from the updated carriage flags. In each compartment, a mechanism can only revert on a given day if no antibiotic with selective pressure for that mechanism is currently present. Carriage mechanisms are evaluated independently, so a drug selecting one mechanism does not preserve unrelated mechanisms. Both compartments use the same cached bacterium-drug-mechanism applicability predicate and a positive-level drug trigger. The cached predicate starts from the molecular mechanism-to-drug map, applies any explicit overrides, and excludes pairs whose susceptible-organism potency is below the model's non-negligible threshold unless explicitly retained. A class-mapped but intrinsically inactive drug therefore does not preserve a mechanism that provides no additional modelled advantage.

Key patterns:
- **Most stable:** Single point mutations (e.g., *gyrA* fluoroquinolone resistance, reversion 0.0001/day) — the mutation barely affects the bacterium's fitness, so it persists for years even without ciprofloxacin pressure
- **Least stable:** Complex multi-gene cassettes (e.g., VanA/VanB vancomycin resistance, reversion 0.002/day; *rpoB* rifampicin resistance, 0.002/day) — these impose significant metabolic costs and are lost relatively quickly without glycopeptide or rifampicin exposure
- **Default** for non-mechanism-specific resistance: 0.0004/day

The full reversion rates by mechanism category:

### Enzymatic Inactivation
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| **KPC** (*bla*KPC) | `0.001` | Plasmid-mediated carbapenemase; moderate maintenance cost. |
| **NDM / VIM** | `0.0015` | Metallo-β-lactamases, frequently on large, high-burden mobile genetic elements. |
| **OXA-48** | `0.0005` | Class D carbapenemase; comparatively lower fitness burden. |
| **ESBL CTX-M / TEM / SHV** | `0.0006` | Standard extended-spectrum β-lactamases. |
| **AmpC DHA** | `0.0006` | Plasmid-mediated AmpC; typical cost profile. |
| **AmpC CMY** | `0.0001` | Often native gene upregulation; minimal fitness loss to maintain. |
| **AmpC derepression** | `0.0002` | Chromosomal *ampC* regulatory mutation (SPACE organisms); low-moderate maintenance cost. |
| **FosA** | `0.0005` | Plasmid-mediated fosfomycin resistance; moderate cost. |
| **CAT** | `0.0005` | Chloramphenicol acetyltransferase. |
| **16S rRMTase** | `0.0005` | Ribosomal RNA methyltransferases conferring high-level aminoglycoside resistance. |



### Target Site Alterations
| Mechanism | Reversion Rate (per day) | Clinical Notes |
| :--- | :--- | :--- |
| **PBP2a / *mecA*** | `0.0009` | High energetic cost associated with maintaining the staphylococcal cassette chromosome *mec* (SCC*mec*). |
| ***erm(B)*** | `0.002` | High reversion rate; target methylation for macrolide-lincosamide-streptogramin B (MLS-B) resistance. |
| **VanA / VanB** | `0.002` | Highly complex target reprograming (D-Ala-D-Ala to D-Ala-D-Lac); significant energetic drain in the absence of glycopeptide exposure. |
| **CFR** | `0.0005` | RNA methyltransferase (oxazolidinone/phenicol cross-resistance). |
| **LiaFSR/Cls** | `0.0015` | Enterococcal *liaFSR/cls* cell-envelope remodeling for daptomycin resistance; significant metabolic burden without daptomycin pressure. |
| **23S rRNA (oxazolidinone)** | `0.001` | Domain V target-site mutation conferring linezolid/tedizolid resistance; measurable fitness cost. |



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
| **Polymyxin regulatory** | `0.0015` | *mgrB/pmrAB/phoPQ/lpx* mutations; chromosomal cell-envelope modification imposes fitness cost similar to MCR-1 plasmid burden. |
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

#### 7.4.1 Community-setting accelerated reversion for hospital-adapted organisms

The baseline reversion rates above are model-wide effective values. For a subset of predominantly nosocomial organisms, the model uses faster effective turnover outside hospital to represent loss of hospital-adapted resistant populations when acute-care selection and transmission are absent. This parameter should not be read as asserting literal genetic reversion of each mechanism. Three considerations motivate the approximation:

1. **Fitness cost without compensatory selection.** The high-level resistance cassettes that define MDR nosocomial clones — carbapenemases (OXA-type, KPC, NDM), VanA/VanB glycopeptide resistance operons, and acquired efflux overexpression — impose a significant metabolic burden. In the hospital, that cost is offset by continuous antibiotic pressure and clonal amplification under selection.  In community carriers who are no longer receiving antibiotics, susceptible wild-type competitors may displace the resistant clone, although clearance timing varies substantially by organism, host, and re-exposure risk. The model therefore treats this as an effective reversion term rather than a literal decolonisation half-life (Andersson DI & Hughes D, 2010; San Millán A & MacLean RC, 2017; Arcilla MS et al., 2017).

2. **Absence of the hospital ecological niche.** MDR hospital-adapted clones (e.g., OXA-23/40/58 *A. baumannii* international clones, VRE *E. faecium* clonal complex 17, MDR *Stenotrophomonas* lineages) are ecologically distinct from their community congeners. They are disproportionately maintained by healthcare exposure, devices, high antibiotic pressure, and repeated admissions rather than by stable community transmission for most people. Community-dwelling individuals who are not recent hospital discharges are therefore less likely to provide the conditions that maintain these clones (Magill SS et al., 2018; Bassetti M et al., 2018; Brooke JS, 2012; Werner G et al., 2008).

3. **Consequence for resistance pool dynamics.** If the same effective turnover is applied in hospital and community settings, discharged patients can repeatedly seed resistant majority/carriage profiles into the community cache and reduce the intended hospital-community contrast. The multiplier is a modelling device for that setting transition; its numerical value is calibrated within this compressed representation rather than measured as a genetic reversion rate.

To address this, the model assigns a per-bacterium `community_mechanism_reversion_multiplier` (default for most bacteria: 1.0, no effect). It multiplies the reversion probability for **every mechanism** of that bacterium when the individual is not currently hospitalised and no active drug selects for that mechanism. On the active-infection side a successful event removes the mechanism from `mechanism_majority` but leaves minority persistence in `mechanism_any`; in carriage it removes the mechanism from `mechanism_microbiome`. The parameter therefore represents effective setting-specific loss from transmissible/dominant or carriage states, not necessarily loss of every resistant cell or literal reversal of the underlying determinant.

The following organisms are assigned a community multiplier of **200×**:

| Organism | Multiplier | Effective community carbapenemase reversion | Effective community VanA/B reversion | Rationale |
|---|---|---|---|---|
| *Acinetobacter baumannii* | 200 | ~0.2/day (50% in ~3 days) | — | Strongly healthcare-associated MDR ecology; multiplier produces rapid effective community turnover in the compressed model |
| *Stenotrophomonas maltophilia* | 200 | ~0.2/day | — | Predominantly healthcare/device-associated infection ecology; multiplier limits persistent community seeding by hospital-adapted profiles |
| *Enterococcus faecium* | 200 | — | ~0.4/day (50% in ~1.5 days) | VRE is strongly hospital/LTCF concentrated; multiplier preserves a marked setting contrast in the model |

These are very strong effective multipliers: in the absence of selecting drugs they can remove majority or carriage mechanism flags within days. They apply to all mechanisms for the listed bacterium, not only the illustrative carbapenemase or VanA/B mechanisms shown in the table. They should therefore be interpreted as calibration approximations for rapid loss of hospital-adapted resistant populations from community-transmissible states, not as direct estimates of genetic determinant loss.

The multiplier can also be set **below 1.0** to slow effective community turnover for organisms whose resistant lineages remain stable after the original drug pressure falls. *N. gonorrhoeae* fluoroquinolone resistance is the configured example: surveillance data from multiple countries show no meaningful rebound in ciprofloxacin susceptibility after fluoroquinolones were withdrawn from gonorrhoea treatment guidelines in 2007 (Unemo M & Shafer WM, 2014; Unemo M et al., 2021; World Health Organization STI fact sheet, 2025). The responsible *gyrA*/*parC* mutations and associated resistant lineages can persist with comparatively low effective cost. The current multiplier is 0.01 (100× slower than the unmodified mechanism rate). As with the 200× values, this bacterium-level multiplier applies to every eligible mechanism, although gyrA is the illustrative calibration target.

| Organism | Multiplier | Effective gyrA reversion | Rationale |
|---|---|---|---|
| *N. gonorrhoeae* | 0.01 | ~0.000001/day for primary gyrA (nominal half-life ~1,900 years when other multipliers are 1) | Strong effective persistence of gonococcal FQ-resistant lineages; surveillance shows no susceptibility rebound after FQ withdrawal |


*Note: The system reserves one remaining placeholder variable (`as_yet_unknown`, baseline rate `0.001`) designated for future empirical calibration.

### 7.5 Resistance floors

At finite simulated population sizes, uncommon infections and mechanisms can be represented by very few active infection-days and retained profiles. The current model therefore includes a small **universal acquisition reseeding guard**, configured under the resistance-floor parameter names, separate from the environmental and ratchet pathways described in Sections 7.7 and 7.8. This guard reduces the chance of losing a mechanism while it remains represented in the retained profile cache; it is not a hard lower bound on resistance prevalence.

| Parameter | Value | Function |
|-----------|-------|----------|
| `resistance_floor_feature_enabled` | 1.0 (on) | Master switch for the floor system |
| `resistance_floor_all_bacteria_enabled` | 1.0 (on) | Enable floors for all bacteria simultaneously |
| `resistance_floor_default_level` | 0.01 | Per-drug probability of entering the acquisition reseeding path |
| `bacteria_{name}_resistance_floor_enabled` | Per-organism | Per-species override when all-bacteria switch is 0.0 |

**How the acquisition reseeding guard works.** Every time a person acquires an infection, the model runs a four-step guard separately for each drug:

1. **Check whether the guard is live.** `calculate_resistance_floor(bacteria, drug, current_day)` checks that (a) floors are enabled for this organism, and (b) the drug class has already been introduced. If both conditions are met it returns `resistance_floor_default_level` (0.01).

2. **Random probability check — does this drug-specific attempt proceed?** With `max_resistance_level = 1.0`, the current implementation gives each eligible drug a 1% probability of entering the mechanism-assignment path. This is not a single 1% classification of the infection: one acquisition can receive zero, one, or several successful per-drug attempts.

3. **Select a mechanism — subject to the causal correctness guard.** The code iterates through applicable resistance mechanisms for that drug and organism. A candidate is eligible only if it is currently present somewhere in the retained global `MechanismCache`; otherwise it is skipped. The method name `mechanism_has_ever_emerged_globally` should therefore be read as a current-cache test, not as a permanent historical registry. Each eligible candidate is assigned with probability `mechanism_assignment_probability_on_any_r_gain = 0.8`, and one successful assignment is sufficient for that drug-specific attempt.

4. **Derive `any_r` from the mechanism.** The acquisition guard does not inject a drug-level resistance value directly. Instead, `propagate_mechanism_resistance()` converts the assigned mechanism into per-drug `any_r` using the same enhancement multipliers as every other resistance pathway.

**Design principle — the guard is deliberately weak.** The 1% value provides repeated opportunities to reseed mechanisms that are already represented in the retained cache. It is far too small to reproduce the observed prevalence of common resistance phenotypes by itself, and it does not guarantee persistence: after a mechanism disappears from every retained profile, the causal guard blocks this pathway from recreating it. Because attempts are made by drug, a mechanism affecting several drugs can also be encountered through more than one attempt. Larger persistence effects are handled through transmission of sampled profiles, explicit environmental reseeding probabilities, and the ratchet.

**Default configuration — all organisms receive reseeding attempts.** With `resistance_floor_all_bacteria_enabled = 1.0` and `resistance_floor_default_level = 0.01`, every organism receives a 1% attempt for each introduced drug, subject to mechanism applicability and the current-cache causal check.

**Alternative configuration — guard disabled.** Setting `resistance_floor_all_bacteria_enabled` to 0.0 removes the universal guard unless an organism-specific override enables it. Resistance is then generated and propagated through the remaining emergence, cache, exogenous, ratchet, and microbiome pathways.

**Details for specific example bacteria**

Stenotrophomonas maltophilia. Intrinsic non-susceptibility to carbapenems, unprotected penicillins, 1st/2nd-generation cephalosporins, macrolides, and most aminoglycosides is encoded directly as near-zero potency values (`potency_when_no_r` ≤ 0.05). Those drug classes never qualify for floor assignment (negligible potency → mechanism not applicable). Acquired resistance to TMP-SMX, fluoroquinolones, and tetracyclines emerges through the standard treatment-selection pathway and receives the same repeated 1% acquisition attempts as other organisms.

Helicobacter pylori. Resistance is driven primarily by chromosomal mutation during treatment courses. Once resistance has been established via de novo selection and while it remains in the retained cache, the universal guard contributes additional reseeding opportunities between relatively infrequent *H. pylori* infection events. The dominant mechanism generating above-guard resistance prevalence is treatment failure without second-line rescue therapy (`helicobacter_pylori_treatment_failure_no_second_line_probability = 0.80`). When first-line clarithromycin-based triple therapy fails, approximately 80% of patients in the model do not proceed to bismuth quadruple rescue therapy — reflecting real-world loss to follow-up, access barriers, poor tolerability of the quadruple regimen, and incomplete prescriber uptake. These patients are left chronically infected with an *H. pylori* strain that has been under macrolide selection pressure; each minority mechanism receives one 0.18/day promotion roll whenever any applicable drug pressure is present, and successful transitions enter the `mechanism_majority` profile recorded in the community mechanism cache. Over a 35-year simulation, this persistent-carrier pool builds up the circulating resistant-strain reservoir, generating realistic clarithromycin and metronidazole resistance prevalence through the standard mechanism-inheritance pathway at acquisition rather than through any floor override. Incidental macrolide and fluoroquinolone exposure from courses prescribed for other organisms contributes additional selection pressure on top of this foundation.

Enterococcus faecium. VRE clonal lineages (CC17) are globally disseminated hospital-adapted strains. While glycopeptide resistance remains represented in the retained cache, the 1% guard provides additional acquisition-side reseeding opportunities between sparse hospitalisation episodes; higher prevalence is generated by the standard selection and profile-transmission pathways.



### 7.6 Cross-resistance groups

Section 7.1 describes the primary biological mechanism-to-drug map: when a resistance mechanism is recorded in a bacterium, the model calculates which drugs that mechanism affects. After that mechanism-based calculation, the model applies a smaller, secondary phenotype-harmonisation step for selected bacteria-drug groups.

These cross-resistance groups are used where resistance to closely related drugs is assumed to move together at the level of the observed drug-resistance phenotype. For each configured group, if one drug has a non-zero `any_r`, the model raises the other drugs in that group to the same maximum `any_r` for that bacterium. This changes the drug-level phenotype only. It does not add, remove, or synchronise the underlying recorded resistance mechanisms.

The groups are therefore not a second mechanism map and should not be read as one-to-one copies of the "drugs affected" entries in Section 7.1. They are a pragmatic phenotype-coherence layer used after mechanism effects have been calculated. Because `any_r` is used downstream in drug activity calculations and phenotypic resistance testing, this layer can influence model behaviour; it is not merely an output-smoothing step.

Some mechanisms affect several related drugs, but not necessarily an entire therapeutic class. For example, ESBLs commonly affect extended-spectrum cephalosporins and aztreonam, while beta-lactamase inhibitor combinations, cefepime, carbapenems, cefiderocol, ceftazidime-avibactam, aztreonam-avibactam, and newer beta-lactam/beta-lactamase inhibitor agents may require separate handling depending on the organism and mechanism. For that reason, cross-resistance groups are deliberately kept narrower than broad therapeutic classes where clinically important exceptions exist.

The full source of truth is `CROSS_RESISTANCE_GROUPS` in `src/config.rs`, applied by `apply_cross_resistance()` in `src/rules/mod.rs`. The examples below are verified against the current configuration and are intended as a readable guide, not a complete duplicated table:

| Bacterium or group | Phenotype groups in current code | Cautious interpretation / edge cases |
|--------------------|-------------------------------------------|-------------------------------------|
| *E. coli* | Plain penicillins (`penicillin_g`, `ampicillin`, `amoxicillin`); early/third-generation cephalosporins (`cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`); fluoroquinolones; aminoglycosides | BL/BLI combinations are excluded from the penicillin and cephalosporin phenotype groups; cefepime, ceftazidime, carbapenems, cefiderocol, and novel BL/BLI agents are not included into the broad cephalosporin group. |
| *K. pneumoniae* | Penicillins; early/third-generation cephalosporins; carbapenems (`meropenem`, `imipenem_c`, `ertapenem`, `meropenem_vaborbactam`); fluoroquinolones | BL/BLI combinations are not grouped with the plain penicillins or ESBL cephalosporin group; newer agents require mechanism-specific handling. |
| *P. aeruginosa* | `piperacillin`, `piperacillin_tazobactam`, `ceftazidime`, and `cefepime`; carbapenems (`meropenem`, `meropenem_vaborbactam`, `imipenem_c`); fluoroquinolones; aminoglycosides | `ceftazidime_avibactam` is explicitly excluded from the AmpC/beta-lactamase group; cefiderocol, ceftolozane-tazobactam, and aztreonam-avibactam are not broad-grouped here. |
| AmpC-risk Enterobacterales | *Enterobacter*, *Citrobacter*, *E. cloacae*, *Morganella*, and *Serratia* use narrow extended-spectrum cephalosporin phenotype groups plus fluoroquinolone groups; some also have aminoglycoside groups | Penicillin+BLI+early-cephalosporin phenotype groups were removed for several AmpC-dominant organisms; AmpC and ESBL mechanisms handle those drug effects through the primary mechanism map. |
| *Shigella* spp. | `ampicillin`/`amoxicillin`; fluoroquinolones; `ceftriaxone`/`cefixime`/`ceftazidime`; tetracycline/doxycycline; macrolides (`azithromycin`, `erythromycin`, `clarithromycin`) | BL/BLI combinations are not included in the penicillin group to avoid levelling them up from ESBL-only resistance. |
| Selected Gram-positive cocci | Staphylococci and enterococci include macrolide-lincosamide groups where appropriate; enterococci include glycopeptide groups; streptococci include narrower macrolide and tetracycline/FQ groupings depending on species | Membership differs by organism. For example, pneumococcal macrolide grouping does not mirror every staphylococcal or enterococcal MLS pattern. |
| Fastidious respiratory/STI organisms | *H. influenzae*, *N. gonorrhoeae*, and *N. meningitidis* use narrower penicillin/FQ/macrolide/tetracycline or ceftriaxone-cefixime groups depending on organism | Inhibitor combinations are excluded from some penicillin groups where inhibitor-susceptible beta-lactamase biology or PBP-mediated resistance requires separate handling. |



---

### 7.7 Environmental and Exogenous Mechanism Floors

The model includes a dedicated family of **environmental / exogenous mechanism floors** that apply only to the non-human fraction of community acquisition. These are distinct from the universal 1% acquisition reseeding guard in Section 7.5.

**What they represent.** The environmental floors are used when resistance can be maintained outside the local human circulating pool. In the classic case, this means agricultural, food-chain, wastewater, or other exogenous reservoirs that keep resistance present even when direct human treatment of that bacterium is weak. In the current configuration they also cover a small number of explicitly modelled non-agricultural exogenous pathways, most notably the rifampicin `rpoB` floors discussed below.

**Where they act.** Each community active-infection acquisition is split by `community_resistance_dilution_factor` into a human-circulating fraction and an exogenous fraction. If the infection is drawn from the exogenous fraction, the model rolls each mechanism independently against the configured floor probability for that bacteria-mechanism pair. In the current implementation this environmental/ratchet floor loop is applied at active-infection acquisition; carriage acquisition uses its own profile-cache gate and does not directly run this exogenous floor loop.

Despite the parameter name, a configured value is therefore a Bernoulli mechanism-assignment probability on an exogenous active-infection acquisition, not a guaranteed lower bound on prevalence in infections, carriage, or the population. The resulting prevalence also depends on how often acquisition follows the exogenous path, the competing human-cache pathway, subsequent transmission and selection, and mechanism loss.

**How they are configured.** These floors use the keys:

```
bacteria_{slug}_mechanism_{mechanism_key}_environmental_floor
bacteria_{slug}_mechanism_{mechanism_key}_environmental_floor_before_{YYYY}
```

Era-specific `_before_{YYYY}` overrides allow the exogenous reservoir to change through time. The parser also recognises `_coselection_floor` and `_coselection_floor_before_{YYYY}` keys for potential cross-organism gut co-selection, but there are no active entries with those suffixes in the default configuration.

**Which organisms currently use them.** The live configuration now contains explicit environmental-floor entries for:

1. Agricultural or food-chain dominated organisms: *Escherichia coli*, *Campylobacter jejuni*, and invasive non-typhoidal *Salmonella* spp.
2. Globally circulating MDR-plasmid support in *Shigella* spp., where the exogenous fraction helps preserve historically established resistant clone backbones.
3. A broader `rpoB` block representing rifampicin-associated bystander or exogenous maintenance pathways for a wider set of organisms, including several Enterobacterales and selected respiratory or enteric pathogens.

So while the original motivation was agricultural resistance ecology, the live implementation is broader: it now represents any explicitly parameterised exogenous pathway that should continue to reseed resistance even when the local human profile cache alone would be insufficient.

---

### 7.8 Dynamic ratchet floor

**Motivation.** Even with long retention of previously observed profiles, a finite stored profile library can under-represent historically established resistance after current selection pressure falls. The ratchet preserves a coarse memory of high previously sampled cache prevalence and uses that memory to support later exogenous acquisition of mechanisms with sufficiently low configured reversion rates.

**Design principle.** The dynamic ratchet is *self-calibrating from the simulation's own sampled history*. It records the highest annually sampled prevalence of each mechanism in the retained profile cache, rounds that peak down to a 10-percentage-point step, and uses the result as a mechanism-assignment probability on subsequent exogenous acquisitions. It does not prevent active-infection, carriage, or population prevalence from falling below the recorded percentage. The biological motivation is that low-effective-cost mechanisms can remain available for reseeding through persistent lineages, co-selection, HGT, and reservoirs after use of the originally selecting drug declines. Integron-borne *sul1* and gonococcal *gyrA* mutations are motivating examples, but the implemented ratchet is a deliberately simplified population-level approximation of those processes.

**Implementation.** Once per simulated year, after that day's profile-cache update, the model pools retained profile entries across all regions and both community and hospital strata and computes the marginal cache prevalence of each (bacteria, mechanism) pair. Every retained entry contributes equally to this statistic; it is not weighted to reproduce regional population sizes. The `peak_mechanism_prevalence[bacteria][mechanism]` table is updated upward-only from these annual snapshots: once a sampled peak has been recorded, it is not allowed to fall.

In the exogenous acquisition path (`!from_human_reservoir`), the effective floor for each mechanism is:

$$\text{effective\_floor} = \max(\text{static\_floor}, \text{ratchet\_floor})$$

where:

$$\text{ratchet\_floor} = \begin{cases} 0 & \text{if reversion\_rate} > 0.003 \text{ /day} \\ \lfloor \text{peak\_prev} / 0.10 \rfloor \times 0.10 & \text{otherwise} \end{cases}$$

capped at 0.50. Here `ratchet_floor` is a probability in the exogenous assignment loop. For example, a recorded peak of 12% yields a 10% assignment probability and a peak of 23% yields a 20% probability. Annual sampling and 10-point steps reduce sensitivity to small fluctuations, but they do not eliminate sampling noise; because the peak is upward-only, a sampled threshold crossing has a lasting effect.

**Eligibility gate — reversion rate ≤ 0.003/day.** In the current code the ratchet applies to mechanisms whose configured base reversion rate is at most 0.003/day. This is a numerical eligibility rule rather than a guarantee that the phenotype is stable in every bacterium and setting; it includes both very-low-reversion mechanisms and some modest-cost mechanisms such as `rpoB`, `ermB`, `vanA`, and `vanB`.

| Mechanism | Reversion rate | Ratchet-eligible y/n |
|-----------|---------------|-----------------------------------|
| GyrA primary (*gyrA*) | 0.0001/day | ✓ |
| GyrA+ParC secondary | 0.0002/day | ✓ |
| Qnr | 0.0001/day | ✓ |
| AcrAB-TolC | 0.0005/day | ✓ |
| Global efflux | 0.0005/day | ✓ |
| CAT | 0.0005/day | ✓ |
| 16S rRMT | 0.0005/day | ✓ |
| TetM/TetO | 0.0005/day | ✓ |
| TetA/B/C efflux | ~0.0001/day | ✓ |
| Folate pathway (sul/dfrA) | 0.0001/day | ✓ |
| AAC/APH | ~0.0001/day | ✓ |
| MphA | ~0.0001/day | ✓ |
| Polymyxin regulatory | 0.0015/day | ✓ |
| ErmB | 0.002/day | ✓ |
| VanA/VanB | 0.002/day | ✓ |
| RpoB | 0.002/day | ✓ |

In conceptual terms, a static environmental floor is a configured exogenous reseeding probability, whereas the ratchet is a history-derived exogenous reseeding probability based on prior cache prevalence.

**Relationship to static floors.** Era-specific static environmental probabilities can provide exogenous seeding during configured historical periods, including periods when the finite profile cache might otherwise build a resistant reservoir slowly. The ratchet has no universal takeover date. At each eligible exogenous acquisition, the model uses the larger of the static and ratchet probabilities. For example, if a *Shigella* mechanism has a current static probability of 0.45 and a recorded cache peak that produces a ratchet probability of 0.40, the effective assignment probability remains 0.45. Neither value implies 45% or 40% total resistance prevalence.

**Diagnostic switch.** `run_pathway_ratchet_enabled` is 1.0 by default and preserves the ratchet behavior described above. Setting it to 0.0 disables only the dynamic ratchet contribution for component-ablation runs; static environmental floors and the separate universal acquisition guard remain active.

**Practical calibration framework for resistance shortfalls.** When the model under-shoots an observed resistance prevalence, the preferred response is *not* to keep inflating mechanism emergence rates above biologically plausible values. For stewardship-policy purposes, the model instead treats most persistent shortfalls as belonging to one of four interpretable categories:

1. **Ratchet reseeding:** use when a mechanism plausibly built up historically under human selection, meets the configured reversion-rate gate, and should retain an acquisition-side historical memory after selection falls.
2. **Environmental / exogenous floor:** use when resistance is continually re-seeded from outside the local human transmission pool, for example through agriculture, food-chain exposure, wastewater, travel-associated importation, or other exogenous reservoirs.
3. **Additional off-model drug exposure:** use when the real selecting drug pressure exists but is not explicitly represented in the treatment pathways already simulated. A canonical example is rifampicin exposure arising from drug-sensitive tuberculosis treatment or other non-modelled indications, which can impose bystander selection on carriage organisms even if those organisms are not themselves the treatment target in the simulation.
4. **Co-selection / linked-package enrichment:** use when the missing feature is not the amount of resistance generated, but the fact that mechanisms tend to arrive together on lineages, integrons, cassettes, or plasmids. A *S. aureus* acquisition-side package completion logic is the template for this class of solution.  (A *S. aureus* enrichment step sits specifically on the **acquisition** side of the model, after a profile has already been drawn from the human circulating reservoir and after `mecA` has already been observed in that sampled profile. This placement is deliberate. It separates two biological ideas that would otherwise be conflated: (1) a mechanism first emerging within a host under selection pressure, and (2) a newly infected person acquiring a strain that belongs to a wider successful MRSA lineage already circulating with a partially linked resistance package.)

**Which lever to prefer.** The intended order of thought is: first ask whether the mechanism needs a history-derived exogenous reseeding probability (ratchet); if not, ask whether it is repeatedly imported from outside the modelled human pool (static environmental probability); if not, ask whether the model is missing a real source of selection pressure (additional off-model drug use); and only then ask whether the main problem is missing correlation structure between mechanisms (co-selection or lineage/package completion). This keeps the model aligned with its purpose as a policy model rather than pushing it towards unnecessary (and costly) molecular detail.

**Implication for future calibration work.** As a working rule, mechanism emergence rates above 1.0 should be treated as a signal to review whether one or more of these four levers is the more appropriate explanation for the observed prevalence. High emergence can still be used where justified, but it should not be the default catch-all for persistence, reseeding, missing selection pressure, or linked resistance ecology.

---

## 8. Microbiome and Carriage

Since the commensal microbiome is the principal reservoir in which resistance is stored, selected by bystander antibiotic exposure, and exchanged between species (Werner G et al., 2008; van Schaik W, 2015; McInnes RS et al., 2020), the model tracks microbiome carriage as a distinct compartment from active infection. A patient treated with ciprofloxacin weeks earlier may still carry fluoroquinolone-resistant *E. coli* in the gut; if that strain subsequently causes a UTI, empiric therapy may fail.

As throughout the model, the microbiome layer is intentionally simplified. We represent the main ecological reservoirs and the policy-relevant consequences of bystander selection, endogenous infection, and within-host persistence, but not the full organism-by-organism spatial ecology that would be required for a dedicated colonisation model.


### 8.1 Carriage compartments

Each bacterium in the model has a designated ecological niche — where it naturally lives in (or on) the body:

| Compartment | Example bacteria | Clinical relevance |
|-------------|-----------------|-------------------|
| Gut | *E. coli*, *K. pneumoniae*, *Enterococcus spp.*, *Shigella*, *Salmonella*, *C. difficile* | Largest reservoir; disrupted by broad-spectrum antibiotics |
| Respiratory | *S. pneumoniae*, *H. influenzae*, *P. aeruginosa*, *A. baumannii*, *M. catarrhalis* | Carriage often precedes pneumonia |
| Respiratory (latent) | *M. tuberculosis* | Modelled as latent infection with stochastic reactivation (~0.1% daily hazard); this represents LTBI rather than conventional mucosal carriage |
| Skin/Soft tissue | *S. aureus*, *S. epidermidis* | Nasal/skin MRSA carriage drives surgical wound infections |
| Genitourinary | *N. gonorrhoeae*, *C. trachomatis*, *M. genitalium*, *T. pallidum*, *S. agalactiae* | Modelled as asymptomatic infection rather than commensal carriage — these are obligate pathogens present without causing symptoms, which enables onward transmission |



These compartment assignments are simplified ecological defaults rather than a full atlas of colonisation niches. They mainly provide the model with the right qualitative reservoirs for bystander selection, endogenous infection, and HGT opportunity.


### 8.2 Resistance in the microbiome

The microbiome serves as a hidden reservoir of resistance. For each organism carried asymptomatically, the model records which resistance mechanisms are present in the carriage compartment. It then calculates a drug-level microbiome resistance measure from those mechanisms using the same multiplicative susceptibility formula used for active infection (Section 7.2). In the code, the carriage mechanism record is `mechanism_microbiome`, the active-infection record is `mechanism_any`, and the derived carriage resistance measure is `microbiome_r`. This keeps carriage and infection resistance biologically aligned rather than tracking them as separate unrelated numerical scores. Clearing a carriage episode resets `mechanism_microbiome` and `microbiome_r` together without altering active-infection mechanisms or `any_r`.

Key dynamics:

| Process | Parameter | Value | Effect |
|---------|-----------|-------|---------------|
| Resistance seeding on acquisition | `microbiome_resistance_multiplier_on_acquisition` | 0.50 | When a person acquires a new carriage episode, there is a 50% probability that the colonising strain inherits the circulating resistance profile (sampled from the stored profile library, just as for infection acquisition — see Section 3.4). If the draw fails, the strain arrives susceptible. The 0.50 value reflects the **colonisation bottleneck**: when a resistant strain from the community pool colonises a new host, only a fraction of transmission events successfully establish a resistant lineage, because susceptible strains in the incoming inoculum can outcompete resistant ones when antibiotic pressure is absent, and because small founding populations are subject to stochastic loss. Carriage acquisition studies consistently show that post-travel or post-admission ESBL carriage rates reach only 20–50% of the prevalence suggested by source-population data, supporting a sub-unity transfer probability (Arcilla MS et al., 2017; Buelow E et al., 2017). The parameter is therefore best interpreted as a colonisation-efficiency discount applied to the profile-sampling step. |
| Established colonies harder to clear | `carriage_duration_log_odds_coefficient` | −0.01/day (caps at −2.0) | The longer a resistant strain has been carried, the harder it is to eradicate — mature colonies are ~7× harder to clear than newly acquired ones |
| Mechanism-level reversion | `run_pathway_reversion_rate_multiplier` | 1.0 | Per-mechanism reversion operates in the microbiome compartment using the same rates, pathway multiplier, and potency-filtered applicability predicate as in the infection compartment (Section 7.4). Each mechanism can only revert when no positive-level active drug is clinically applicable to that bacterium-mechanism pair; selection for another mechanism does not block reversion. |
| De-novo emergence under treatment | `run_pathway_microbiome_de_novo_multiplier` | 1.0 | When at least one positive-level active drug applies to an absent carriage mechanism, that mechanism receives one daily emergence attempt via the microbiome pathway (Section 7.3), scaled by this multiplier and the compatibility multiplier `microbiome_de_novo_multiplier`. Concurrent applicable drugs do not add attempts. Emergence writes directly to `mechanism_microbiome`; transfer into infection occurs through the separate bridge pathways below. |
| Carrier → infection bridge | `carrier_resistance_inheritance_probability` | 0.50 | When a new same-organism infection is established in a carrier, each mechanism in `mechanism_microbiome` is independently considered for transfer to `mechanism_any` (see Section 3.3) |
| Infection ↔ microbiome transfer | `microbiome_resistance_transfer_probability_per_day` | 0.0001 | When both compartments are present and contain different mechanisms, a daily random probability check can trigger copying in both directions so that `mechanism_any` and `mechanism_microbiome` end up containing the combined profile. The drug-level resistance measures are then recalculated. |
| HGT into the microbiome | (see Section 9) | — | When a horizontal gene transfer event fires, the mechanism is written only to recipient compartments where that bacterium is present: `mechanism_microbiome` for carriage and `mechanism_any` for active infection |



## 9. Horizontal Gene Transfer (HGT)

Horizontal gene transfer (HGT) — the interspecies sharing of resistance determinants, as seen when the same ESBL plasmids appear across *E. coli*, *Klebsiella*, and *Proteus*  — is a major driver of resistance spread and is modelled explicitly.

The HGT layer is necessarily schematic. We preserve the major ecological compatibilities and the main amplifiers of transfer risk, but we do not attempt plasmid-by-plasmid reconstruction, incompatibility typing, or ward-level contact-network modelling. At the scale of the present model, that additional detail would be difficult to support empirically and would add substantial runtime and calibration burden without clearly improving the policy comparisons of interest.


### 9.1 Transfer compatibility

Not all bacteria can exchange genes equally. Transfer compatibility is not represented as a single species-to-species lookup table. Instead, each bacterium group is assigned to a **plasmid pool**, and the baseline pairwise HGT hazard is generated from that pool structure before the Section 9.2 multipliers are applied.

The pool mapping is:

- **GramPositive pool**: Staphylococci and Streptococci — both gram-positive groups share this pool; cross-group (Staph↔Strep) HGT operates at the lower cross-group rate (10× below within-group), and `TargetSitePbp2aMecA`/`EnzymeBlaZ` are restricted to Staphylococci so mecA cannot transfer into Streptococci regardless of HGT rate
- **EntericGramNegative pool**: Enterobacterales, non-fermenters, and enteric pathogens
- **RespiratoryGramNegative pool**: fastidious respiratory/genitourinary organisms
- **Anaerobe pool**: anaerobes
- **No-transfer structural exclusion**: spirochetes, helicobacters, and mycobacteria are assigned to `None` and therefore have baseline HGT probability `0.0`. This reflects well-established biological constraints: *Treponema pallidum* and other spirochetes lack classical conjugation machinery and have not been shown to exchange resistance plasmids under clinical conditions; *Helicobacter* and *Campylobacter* acquire resistance primarily through chromosomal mutation rather than plasmid-borne transfer; and *Mycobacterium tuberculosis* has a cell wall that is largely impermeable to conjugative pili and lacks the surface-exposed mating-pair formation proteins required for classical plasmid transfer (Carattoli A, 2009; Borger AL et al., 2023). Excluding these lineages from the HGT pool does not mean they never acquire new resistance — their emergence rates (Section 7.3) handle de novo mutation — but they do not participate in the plasmid-sharing network that connects Gram-negatives and Gram-positives in the model

The baseline compatibility assignments are:

| Donor-recipient relationship | Baseline pairwise HGT probability |
|-----------------------------|-----------------------------------|
| Same plasmid pool, same bacteria group | `1e-5` |
| Same plasmid pool, different bacteria group | `1e-6` |
| Enteric Gram-negative <-> respiratory Gram-negative | `3e-8` |
| Enteric Gram-negative <-> anaerobe | `3e-7` |
| Anaerobe <-> anaerobe | `1e-5` |
| All other cross-pool combinations | `0.0` |

These values are effective within-model hazards and their purpose is to preserve the current ordering in the code: transfer is easiest within the same ecological/plasmid pool, much weaker across the small set of allowed cross-pool bridges, and structurally absent for excluded groups.


### 9.2 The HGT process

Each day, after the per-bacterium updates, the model evaluates potential gene transfer events among bacteria present in the same person, whether as active infection, asymptomatic carriage, or both. The model evaluates HGT separately for each resistance mechanism, allowing different plasmids or chromosomal determinants (for example KPC and *mcr-1*) to transfer independently rather than as a single all-or-nothing package. Donor eligibility is determined from the compartments in which the organism is actually present: `mechanism_any` supplies active-infection mechanisms and `mechanism_microbiome` supplies carriage mechanisms. Bacteria do not restrict donation only to the dominant active-infection strain; minority active-infection mechanisms and carriage mechanisms can also donate, but they receive the configured minority-donor multiplier because carriage has no separate majority-mechanism flag.

When an HGT event occurs, the transferred mechanism is recorded in each compartment where the recipient organism is present: `mechanism_any` for an active infection and `mechanism_microbiome` for asymptomatic carriage. If the recipient is present in both, both records are updated. A carriage-only transfer does not create phantom active-infection mechanism state. All HGT rates are scaled by the run-level `run_pathway_hgt_multiplier` (default 1.0) and by `counterfactual_resistance_multiplier`.

| Step | Parameter | Value | Clinical parallel |
|------|-----------|-------|-------------------|
| Run-level HGT scaling | `run_pathway_hgt_multiplier` | 1.0 | Calibration knob — scales all HGT rates up or down uniformly |
| Base pairwise transfer rate | `hgt_prob_{donor}_to_{recipient}` matrix (defaults from Section 9.1) | `0.0` to `1e-5` | Baseline compatibility hazard determined by plasmid-pool ecology and donor-recipient pairing before context multipliers are applied |
| Amplification during antibiotic therapy | hgt_antibiotic_pressure_multiplier | 1.50 (×1.5) | Antibiotic use leads to activation of mobile genetic elements and increased conjugation rates by 50% (Beaber JW et al., 2004) — one of the reasons antibiotic use drives resistance even beyond the target pathogen |
| Hospitalization boost | hgt_hospital_multiplier | 4.0 (×4.0) | Captures increased transmission risks in clinical environments where close physical proximity and shared infrastructure elevate exchange. |
| Co-infection baseline | hgt_coinfection_multiplier | 1.25 (×1.25) | Active multi-pathogen infections slightly increase the probability of genetic collision. |
| Microbiome-only penalty | hgt_microbiome_only_penalty | 0.65 (×0.65) | Asymptomatic carriage interactions are less frequent than active infection environments. |
| Gut compartment boost | hgt_gut_compartment_multiplier | 2.0 (×2.0) | The gut has higher bacterial density and provides more conjugation opportunities compared to skin or respiratory tracts. |
| Minority donor penalty | hgt_minority_donor_multiplier | 0.20 (×0.20) | If an active-infection donor carries the mechanism outside its dominant strain, or if the mechanism is supplied only by carriage where no majority flag exists, its probability of successful conjugation is penalized by 80%. |

`microbiome_resistance_transfer_probability_per_day` is a separate parameter used for within-host infection↔microbiome mechanism exchange (Section 8.2), not for inter-species HGT.

The absolute HGT probabilities are intentionally low and should be interpreted as effective daily hazards at the model scale. Their main purpose is to preserve plausible relative ordering between low-contact community settings, antibiotic-stressed microbiomes, and high-contact hospital environments.

## 10. Mortality

The model tracks mortality from three sources: background (non-infection) causes, **bacterial infection induced sepsis**, and **non-sepsis (bacterial) infection death** (direct tissue damage, toxin production, or chronic complications of infection that do not involve the sepsis cascade). This dual-pathway architecture reflects the clinical reality that different pathogens kill through fundamentally different mechanisms (Rudd KE et al., 2020).

### 10.1 Background mortality

Every individual faces a baseline daily death risk shaped by age, sex, region, immune status, and the simulated calendar year. The probability is computed via a logistic model whose total log-odds sum the following components:

| Factor | Parameter | Default value | Effect |
|--------|-----------|--------------|-------|
| Baseline intercept | `background_mortality_baseline_log_odds` | -14.3 | Global anchor for the daily risk |
| Historical improvement | `mortality_baseline_1930_multiplier` / `mortality_baseline_2035_multiplier` / `mortality_improvement_half_life_years` | ×3 / ×1 / 35 yrs | Normalized exponential decline from a 3× higher 1930 rate to the configured 2035 reference rate exactly; half-life controls how front-loaded that improvement is |
| Linear age effect | `log_odds_mortality_per_year_of_age` | 0.055 | Each year of age adds a constant increment to log-odds (≈ ×1.06/year on the odds scale) |
| Elderly frailty acceleration | `log_odds_mortality_per_year_of_age_squared` | 0.008 | Quadratic term applied **only above age 80** — steepens mortality in the very elderly without making age-90 mortality implausibly extreme |
| Region | `log_odds_mortality_region_{name}` | N. America 0; S. America +0.26; Africa +0.69; Asia +0.18; Europe −0.11; Oceania 0 | Reflects broad differences in background mortality environment, healthcare access, and non-communicable disease burden |
| Sex | `log_odds_mortality_sex_male` / `_female` | +0.095 / −0.105 | Male ≈ ×1.1, female ≈ ×0.9 all-cause mortality differential |
| Immunosuppression | `log_odds_mortality_immunosuppressed` | +0.916 | ≈ ×2.5 higher risk when `immunodeficiency_type` is set |
| Hospital status | `log_odds_mortality_hospitalized` | +0.262 | ≈ ×1.3 higher risk while in hospital (captures inpatient case-mix and residual non-infectious acuity rather than HCAI, which is modelled separately) |

All parameters operate on a log-odds scale and sum additively before the logistic transform, so their effects multiply on the probability scale.

They should be read as effective demographic mortality-shape terms rather than direct life-table fits for any single country or year. Their role is to preserve the globally familiar pattern of sharply rising all-cause mortality with age and frailty while allowing the simulation's infection-specific pathways to add the AMR-relevant excess risk on top.

In the current implementation, the sex term is a lifelong multiplicative shift rather than an age-specific late-life modifier. That is a simplification: real male-female mortality gaps vary by age, cause, and setting.

Background mortality is treated as a competing risk alongside infection-specific death pathways rather than being added on top of them. Each day, the model checks for death in the following order: sepsis, drug toxicity, non-sepsis infection death, then background mortality. This means acute infectious deaths displace some deaths that would otherwise have been labelled as background mortality, ensuring each person receives at most one cause of death per time step.


### 10.2 Sepsis mortality

Sepsis is the primary death pathway for classic invasive bacterial pathogens. When an individual's infection progresses to sepsis (see Section 4.3), the model applies an escalated daily death risk using a logistic model. The probability of dying from sepsis each day depends on age, immune status, bacterial burden, and access to hospital care. Without effective antibiotics, sepsis is rapidly fatal, and resistant organisms that are untreatable with empiric therapy represent the principal scenario of concern (Murray CJL et al., 2022).

The mortality model also applies the existing `sepsis_death_log_odds_not_under_care` penalty (+1.4) when none of the three person-level care signals is present: current antibiotic use, hospitalisation, or bacterial identification for a still-active infection. This term represents nonspecific medical and supportive care rather than organism-specific antimicrobial effectiveness. Effective antibiotics additionally alter bacterial activity, burden, and clearance through the separate treatment pathways described in Section 6.

Since sepsis mortality varies enormously by organism — from near-zero for non-invasive STI pathogens to >30% for *S. aureus* bacteraemia (Tong SYC et al., 2015) — the model assigns per-bacterium sepsis baseline log-odds. The current configured values are:

| Bacterium | Sepsis baseline log odds | Clinical rationale |
|-----------|----------------|-------------------|
| *Acinetobacter baumannii* | -6.6 | Healthcare-associated invasive pathogen, especially ventilator-associated pneumonia and ICU bacteraemia; often severe and drug-resistant |
| *Citrobacter* spp. | -8.6 | Opportunistic Enterobacterales; can cause invasive infection but less commonly than *E. coli* or *Klebsiella* |
| *Enterobacter* spp. | -7.1 | Opportunistic hospital-associated Enterobacterales with meaningful bloodstream-infection potential |
| *Enterococcus faecalis* | -6.9 | Endocarditis and line-related bacteraemia |
| *Enterococcus faecium* | -6.1 | Hospital-acquired bloodstream infections, especially VRE |
| *Escherichia coli* | -8.5 | Most common Gram-negative bloodstream isolate; UTI-source sepsis usually less severe than highly invasive ICU pathogens (Poolman JT et al., 2016) |
| *Klebsiella pneumoniae* | -6.4 | Gram-negative sepsis; carbapenem-resistant strains carry high mortality (Xu L et al., 2017) |
| *Morganella* spp. | -7.1 | Opportunistic Enterobacterales associated with urinary, wound, and healthcare-associated invasive infection |
| *Proteus* spp. | -7.1 | UTI-associated Enterobacterales with potential for urosepsis, especially in older or catheterised patients |
| *Serratia* spp. | -7.3 | Opportunistic healthcare-associated Enterobacterales with bloodstream-infection potential |
| *Providencia stuartii* | -12.0 (fallback) | No organism-specific override configured; handled using the global fallback sepsis baseline |
| *Pseudomonas aeruginosa* | -5.8 | High mortality in ICU infections; often in immunocompromised hosts (Bassetti M et al., 2018) |
| *Stenotrophomonas maltophilia* | -7.3 | Opportunistic non-fermenter, mainly in severely ill or immunocompromised patients |
| *Staphylococcus aureus* | -7.3 | Aggressive bloodstream pathogen; 20-30% mortality in bacteraemia (Tong SYC et al., 2015) |
| *Staphylococcus epidermidis* | -7.3 | Device-associated and line-related infection; lower virulence than *S. aureus* but clinically important in hospitalised patients |
| *Streptococcus pneumoniae* | -9.8 | Invasive pneumococcal disease can cause sepsis, but many infections are respiratory and non-bacteraemic |
| *Salmonella enterica* serovar Typhi | -7.3 | Enteric fever with potential for systemic invasive disease |
| *Salmonella enterica* serovar Paratyphi A | -8.3 | Enteric fever with occasional septic complications |
| Invasive non-typhoidal *Salmonella* spp. | -8.5 | Invasive non-typhoidal salmonellosis; high mortality in sub-Saharan Africa (Stanaway JD et al., 2019) |
| *Shigella* spp. | -11.3 | Primarily dysentery/dehydration mortality; sepsis is not the dominant pathway |
| *Neisseria gonorrhoeae* | -23.3 | Disseminated gonococcal infection is rare |
| *Streptococcus pyogenes* | -5.8 | Invasive GAS disease including necrotising fasciitis and toxic shock; STSS can be highly lethal (Carapetis JR et al., 2005) |
| *Streptococcus agalactiae* | -6.3 | Neonatal and pregnancy-associated sepsis (Seale AC et al., 2013) |
| *Haemophilus influenzae* | -8.5 | Invasive respiratory pathogen, especially in young children, older adults, and unvaccinated populations |
| *Chlamydia trachomatis* | -18.3 | STI; essentially never causes classic bacterial sepsis |
| *Mycoplasma genitalium* | -12.0 (fallback) | No organism-specific override configured; handled using the global fallback sepsis baseline |
| *Vibrio cholerae* | -8.3 | Mortality is usually dehydration-mediated, but severe systemic illness is possible in vulnerable hosts |
| *Neisseria meningitidis* | -7.2 | Meningococcal disease; rapid sepsis progression with purpura fulminans and DIC (Tunkel AR et al., 2004; van de Beek D et al., 2012) |
| *Listeria monocytogenes* | -7.3 | Invasive disease in neonates, pregnancy, older adults, and immunocompromised patients |
| *Clostridioides difficile* | -10.3 | Deaths are often toxin-mediated colitis rather than classic bloodstream sepsis |
| *Bacteroides fragilis* | -12.0 (fallback) | No organism-specific override configured; handled using the global fallback sepsis baseline |
| *Campylobacter jejuni* | -19.3 | Usually enteritis; bacteraemia/sepsis is rare |
| *Enterobacter cloacae* | -7.0 | Opportunistic hospital-associated Enterobacterales with bloodstream-infection potential |
| *Yersinia enterocolitica* | -8.8 | Rare sepsis, mainly in iron-overload or immunosuppressed patients |
| *Moraxella catarrhalis* | -10.4 | Usually respiratory mucosal infection; invasive sepsis is uncommon |
| *Treponema pallidum* | -10.3 | Syphilis mortality is typically chronic, congenital, or cardiovascular/neurologic rather than acute sepsis |
| *Bordetella pertussis* | -11.1 | Deaths are usually respiratory, especially in infants, rather than classic sepsis |
| *Helicobacter pylori* | -500.0 | Gastric pathogen; deaths are chronic/cancer-related, and the implementation also forces zero sepsis risk when this is the sole active infection |
| MDR *Mycobacterium tuberculosis* | -38.0 | Chronic mycobacterial disease rather than acute bacterial sepsis |
| *Mycoplasma pneumoniae* | -12.0 (fallback) | No organism-specific override configured; handled using the global fallback sepsis baseline |
| *Legionella pneumophila* | -12.0 (fallback) | No organism-specific override configured; handled using the global fallback sepsis baseline |
| *Burkholderia cepacia* complex | -12.0 (fallback) | No organism-specific override configured; handled using the global fallback sepsis baseline |

Rows marked as fallback use the global `sepsis_baseline_log_odds = -12.0` because no organism-specific `{bacterium}_sepsis_baseline_log_odds` override is configured in `src/config.rs`.



These per-bacterium sepsis baselines are qualitative severity orderings anchored to widely observed differences between invasive and non-invasive pathogens, not claims of portable case-fatality estimates across all settings. Real-world sepsis mortality depends heavily on time-to-treatment, ICU access, comorbidity structure, and health-system capacity, so the model uses these terms mainly to maintain defensible ranking and then lets care access, treatment effectiveness, and syndrome site shape realised mortality in each branch (Rudd KE et al., 2020; Murray CJL et al., 2022).


### 10.2.1 Per-organism sepsis case-fatality adjustment

In addition to the per-bacterium sepsis entry baseline (Section 10.2), the model supports an **additive per-organism log-odds adjustment to the daily death probability given sepsis** (parameter name: `{organism}_sepsis_death_log_odds_override`). This term is added on top of all other factors in the sepsis death calculation — age, region, bacterial burden, medical-care status, and immunosuppression. Where multiple bacteria are simultaneously septic, the largest override across all septic organisms takes effect.

Three organisms currently receive non-zero adjustments:

| Bacterium | CFR adjustment | Relative CFR | Clinical rationale |
|-----------|---------------|--------------|-------------------|
| *N. meningitidis* | +0.69 | ≈×2 | Purpura fulminans and DIC; meningococcal sepsis has among the highest 24-hour CFR of any bacterial pathogen (Tunkel AR et al., 2004; van de Beek D et al., 2012) |
| *S. aureus* | +0.41 | ≈×1.5 | Infective endocarditis and MRSA bacteraemia; 30-day mortality 20–30% even with appropriate therapy (Tong SYC et al., 2015) |
| *A. baumannii* | +0.69 | ≈×2 | XDR ventilator-associated pneumonia and bloodstream infection; attributable mortality >30% in carbapenem-resistant strains (Bassetti M et al., 2018) |

All other organisms default to 0.0 (no adjustment).


### 10.3 Non-sepsis infection death

Not all infection-related deaths involve sepsis. Many pathogens kill through tissue-specific mechanisms: *V. cholerae* through fatal dehydration (Ali M et al., 2015), *B. pertussis* through infantile respiratory failure (Yeung KHT et al., 2017), *H. pylori* through gastric adenocarcinoma (Plummer M et al., 2015), *T. pallidum* through tertiary and congenital syphilis (Korenromp EL et al., 2019), and *C. difficile* through toxic megacolon (Guh AY et al., 2020). These deaths would not be captured by the sepsis pathway alone.

The model evaluates a **daily non-sepsis infection death probability** for every active infection that is *not* already progressing through the sepsis pathway. The probability is computed via a logistic model:

$$P(\text{non-sepsis death}) = \frac{1}{1 + \exp(-\text{log-odds})}$$

where the total log-odds combines:

| Component | Default | What it captures |
|-----------|---------|-----------------|
| Base (`infection_non_sepsis_base_log_odds`) | −9.0 | Global intercept — very low daily risk |
| Per-bacterium adjustment | 0.0 (default) | How lethal this organism is via non-sepsis mechanisms |
| Per-syndrome adjustment | 0.0 (default) | How dangerous this body site is |
| Bacterial level × coefficient | level × 0.0 | Higher burden → higher risk |
| Age adjustment | varies | Infants and elderly at higher risk |
| Hospital adjustment | 0.0 | Modified risk in hospital |
| Immunosuppression | 0.0 | Additional risk for immunocompromised |



The per-bacterium adjustments are the primary calibration lever. **Negative values** reduce non-sepsis death (used for organisms whose deaths are over-represented at the base rate), while **positive values** increase it (used for organisms whose real-world deaths come from non-sepsis mechanisms that would otherwise be invisible to the model):

| Bacterium | Adjustment | Rationale |
|-----------|-----------|-----------|
| *C. trachomatis* | −5.0 | STI with near-zero real-world mortality; base rate produced 128× over-death |
| *M. genitalium* | −4.5 | STI with essentially no deaths; base rate produced 66× over-death |
| *N. gonorrhoeae* | −2.5 | Gonorrhoea is rarely fatal; base rate produced 11.6× over-death |
| *M. pneumoniae* | −0.7 | Low-mortality respiratory pathogen |
| *C. jejuni* | −0.5 | Self-limiting gastroenteritis in most cases |
| *S. epidermidis* | −8.0 | Predominantly indolent, device-associated infection; acute sepsis is uncommon. Untreated prosthetic valve endocarditis carries ~20–30% mortality, but this is a subacute process partially captured by non-sepsis death parameters rather than the acute sepsis baseline |
| *S. maltophilia* | −4.0 | Some mortality via pneumonia progression, but limited |
| *B. pertussis* | +4.0 | Deaths from respiratory failure in infants, not sepsis (Yeung KHT et al., 2017) |
| *T. pallidum* | +3.5 | Tertiary/congenital syphilis deaths (Korenromp EL et al., 2019) |
| *V. cholerae* | +2.5 | Death from dehydration, not bacteraemia (Ali M et al., 2015) |
| *C. difficile* | +2.0 | Colitis and toxic megacolon deaths (Guh AY et al., 2020) |
| *S. pyogenes* | +3.0 | STSS and superantigen (SPE-A/C/SMEZ)-mediated rapid death independent of bacterial burden, plus rheumatic heart disease and post-streptococcal complications (Carapetis JR et al., 2005; Watkins DA et al., 2017) |
| *B. fragilis* | +1.5 | Intra-abdominal abscess mortality |
| *H. pylori* | +1.7 | Gastric cancer deaths; essentially zero sepsis risk (Plummer M et al., 2015) |
| *Shigella* spp. | +1.0 | Dysentery deaths in children; sepsis pathway contributes minimally (Troeger C et al., 2018) |



This dual-pathway design means that the model can reproduce both the typical sepsis mortality pattern (where broad-spectrum antibiotics and ICU care determine survival) and the non-sepsis mortality pattern (where the primary driver may be dehydration, organ-specific damage, or chronic sequelae).

The non-sepsis adjustments are therefore best viewed as compensating structural terms for important death pathways that a pure sepsis model would miss, rather than as direct organism-specific fatality estimates. That is especially important for globally important syndromes such as cholera, pertussis, diarrhoeal disease, and chronic sequelae-associated infections, where the pathway to death is real but not well represented by bloodstream invasion alone.


### 10.4 Infection mortality — syndrome multipliers

Both death pathways are modulated by the anatomical site of infection. The syndrome multipliers reflect how dangerous each body site is:

| Syndrome | Multiplier | Rationale |
|----------|-----------|-----------|
| Genital | 0.05 | Rarely fatal (localised mucosal infections) |
| Skin / Ear | 0.1 | Low systemic risk unless secondary bacteraemia |
| UTI | 0.5 | Usually self-limiting, but untreated or inadequately treated UTI — particularly in the elderly, pregnant women, or those with structural abnormalities — can ascend to urosepsis and carries meaningful mortality in that context |
| Bone/Joint | 0.8 | Serious but often slow-progressing; mortality arises from surgical complications, but also from haematogenous spread — bacteraemia, sepsis, and seeding of distant sites (e.g. vertebral osteomyelitis seeding cardiac valves) |
| Intra-abdominal | 1.5 | Peritonitis carries high mortality even with surgery |
| Respiratory | 1.5 | Pneumonia — leading infectious cause of death globally (GBD 2019 Lower Respiratory Infections Collaborators, 2022) |
| CNS | 3.0 | Meningitis/brain abscess — poor penetration of many antibiotics (Tunkel AR et al., 2004) |
| Bloodstream | 4.0 | Bacteraemia/sepsis — the most immediately life-threatening |



These syndrome multipliers encode the broad global ordering in which bloodstream and CNS infections are most lethal, respiratory and intra-abdominal infections are high-risk, and genital or superficial infections are usually much less fatal unless they progress, which is the main pattern needed for policy comparisons in the model.

---

## 11. Potential future model uses

Although the initial application focuses exclusively on reconstruction of global emergence of resistance, the model architecture is designed to support **policy comparison** whereby the simulation can branch off between different policy options in order to predict the effects of the alternative policies.  Potential future applications include comparing antibiotic stewardship packages (e.g., narrower empiric prescribing, expanded susceptibility testing, shorter course durations), evaluating the trade-off between restricting reserve drugs and preserving last-resort efficacy, and quantifying the projected impact of improved point-of-care diagnostics on resistance trajectories and mortality over multi-decade horizons.

---



## 12. Limitations

The central design judgement has been to retain the features most likely to matter for stewardship, diagnostics, access, and mortality questions, while omitting layers of nuance that would make a model of this scope difficult to calibrate, computationally too burdensome, or unnecessarily difficult to interpret. The main limitations are therefore not incidental omissions but deliberate trade-offs made in order to keep the model usable for the policy questions it is intended to address:

Several of the appendices that follow list exact configuration values and enum definitions. Those tables are included for transparency and reproducibility, but they should still be read in the context established above: many entries are implementation defaults, calibration targets, or structural coding choices rather than direct empirical measurements. Where this document presents an exact value, that should not automatically be interpreted as implying an equivalent degree of empirical certainty.

### 12.1 Calibration target provenance and interpretation

The calibration summary generated by the code reports a small number of headline targets plus larger per-organism and per-drug target families. These values are used to set the order of magnitude and relative structure of the model. Some are close to directly observed surveillance quantities, but many are derived calibration anchors that combine multiple evidence streams and should not be read as single observed global statistics.

| Target block | Target value or target family | Main source(s) | Direct estimate vs derived calibration anchor | Important caveat |
|--------------|-------------------------------|----------------|----------------------------------------------|------------------|
| Infection deaths | 6.4 million model-scope bacterial infection deaths per year | Ikuta KS et al., 2022; Murray CJL et al., 2022; GBD 2021 AMR Collaborators, 2024 | Derived calibration anchor | Like-for-like target for the Figure 1 simulation numerator: sepsis plus non-sepsis bacterial infection deaths, excluding H. pylori and MDR-TB. It is the rounded sum of per-organism mortality targets after those exclusions, so broader bacterial mortality estimates can be higher. |
| People on antibiotics on an average day | 100 million people | WHO Collaborating Centre for Drug Statistics Methodology, accessed 2026; Browne AJ et al., 2021; Klein EY et al., 2018; ECDC ESAC-Net, 2024; World Health Organization GLASS-AMU report, 2025 | Derived person-day calibration anchor | Derived from DDD-based consumption evidence and class-use patterns, but DDDs are not identical to unique people treated because dose, duration, combination therapy, stock buffering, wastage, and non-human channels differ by setting. |
| Annual bacterial infection incidence | 15% of world population per year | GBD/GRAM burden estimates, WHO STI estimates, WHO FERG/Havelaar et al., 2015, WHO Global Tuberculosis Report 2025 | Derived calibration anchor with high uncertainty | No single global surveillance system directly counts all bacterial infections across syndromes; this target is used to keep the all-organism acquisition scale in a plausible range. |
| Incident bacterial sepsis | 30 million cases per year | Rudd KE et al., 2020 | Derived bacterial-subset assumption | Rudd et al. estimate all-cause sepsis; the model simulates bacterial infections only, so this is a bacterial subset anchor rather than a direct observed global count. |
| Per-bacteria infection incidence | Per-organism annual infection targets in the calibration summary | Ikuta KS et al., 2022; Murray CJL et al., 2022; Rowley J et al., 2019; WHO STI fact sheet, 2025; Havelaar AH et al., 2015; WHO Global Tuberculosis Report 2025 | Mix of direct disease estimates and derived calibration anchors | Common pathogens and named diseases are better anchored; smaller healthcare-associated organisms use placeholder estimates with high uncertainty. |
| Per-bacteria mortality | Per-organism death targets in the calibration summary | Ikuta KS et al., 2022; Murray CJL et al., 2022; WHO Global Tuberculosis Report 2025; Plummer M et al., 2015; Troeger C et al., 2018; Guh AY et al., 2020 | Mix of direct estimates, modelled burden estimates, and derived anchors | Pathogen attribution differs between sources, and polymicrobial attribution can make summed organism deaths exceed headline infection deaths. |
| Carriage prevalence | Per-organism carriage target percentages | Human Microbiome Project Consortium, 2012; Wertheim HFL et al., 2005; Bogaert D et al., 2004; Gorrie CL et al., 2017; Verani JR et al., 2010; Forslund K et al., 2013 | Mix of colonisation estimates and model-design anchors | Some entries are model-design zeros or structural reservoirs rather than direct carriage prevalence estimates; microbiome resistance prevalence is a simulated state, not a direct target for every species. |
| Drug-class shares | Drug-use and class-share targets | Browne AJ et al., 2021; Klein EY et al., 2018; van Boeckel TP et al., 2014; ECDC ESAC-Net, 2024; World Health Organization GLASS-AMU report, 2025 | Derived consumption anchor | Sales and DDD/1000/day data inform relative exposure but do not map one-to-one to active simulated users or syndrome-specific prescribing decisions. |
| Bacteria/drug resistance targets | Surveillance isolate resistance targets by organism, drug, and care setting | World Health Organization Global antibiotic resistance surveillance report, 2025; World Health Organization GLASS dashboard, accessed 2026; ECDC EARS-Net, 2024; CDC AR Threats, 2019; CDC COVID-19 Special Report, 2022; Murray CJL et al., 2022; WHO Global Tuberculosis Report 2025; Unemo M et al., 2021; WHO STI fact sheet, 2025 | Mostly surveillance-isolate targets, with some derived anchors | Routine clinical isolates overrepresent tested, severe, and healthcare-associated infections compared with all simulated infections; breakpoint methods, sampling frames, and laboratory coverage vary across time and region. |

1. **Abstract drug levels**: Antibiotic concentrations are modelled as dimensionless units rather than true pharmacokinetic concentrations (mg/L). This allows the model to capture the *relative* dynamics of drug accumulation and clearance, but it means model values cannot be compared directly with MIC breakpoints, therapeutic drug monitoring results, or compartment-specific pharmacokinetic measurements from clinical microbiology or pharmacology practice. In particular, the model does not implement pharmacokinetic/pharmacodynamic (PK/PD) target-attainment analysis — it does not compute AUC/MIC or T>MIC indices, nor does it model the Cmax and distribution volume differences between patient subgroups (e.g., critically ill patients with altered volumes of distribution, or renal impairment affecting aminoglycoside and vancomycin clearance). Full mechanistic PK/PD frameworks can generate organism-specific probability-of-target-attainment curves and inform optimal dosing regimens (Nielsen EI & Friberg LE, 2013), which is beyond the scope of this policy-comparison model. The practical consequence is that the model's drug-level dynamics can reproduce the broad direction of resistance selection associated with sub-therapeutic exposure, but cannot support dosing-optimisation analyses or precisely model regimens where PK/PD target attainment drives clinical outcome.

2. **No explicit strain competition**: Within the microbiome, resistant and susceptible strains do not explicitly compete for ecological resources. The model therefore cannot represent scenarios in which clonal replacement, compensatory evolution, or near-cost-free resistance leads to durable dominance of resistant strains in the absence of ongoing antibiotic selection. That said, the model does capture several distinct mechanisms by which antibiotic use promotes resistance in the microbiome: (i) a *microbiome disruption reservoir* that accumulates while drugs are active and decays with a configurable half-life (`antibiotic_disruption_decay_half_life_days`), raising future colonisation risk; (ii) *de novo resistance emergence* in the microbiome, where each absent applicable mechanism receives one daily attempt under binary positive-level drug pressure, scaled by `run_pathway_microbiome_de_novo_multiplier` and the compatibility multiplier `microbiome_de_novo_multiplier`; (iii) *selective maintenance* of existing resistance — mechanisms only revert when no selecting drug is active, so ongoing treatment blocks loss of resistance; (iv) daily bidirectional *infection–microbiome resistance spillover* governed by `microbiome_resistance_transfer_probability_per_day`; and (v) *horizontal gene transfer amplified by antibiotic pressure* through `hgt_antibiotic_pressure_multiplier`. Together these five pathways mean that antibiotic exposure promotes and sustains microbiome resistance through multiple complementary routes, even though the model does not track explicit clonal competition between resistant and susceptible lineages.

3. **No within-host spatial structure**: Infections are treated as homogeneous within a body compartment. Biofilm formation, abscess walling-off, source control, and planktonic-versus-sessile distinctions are not modelled. The model therefore cannot reproduce the full treatment implications of deep-seated infection architecture, even though such structure is often decisive in real clinical microbiology and infectious diseases practice.

4. **Static vaccine model**: Vaccinated individuals have a fixed proportional reduction in infection risk. Vaccine effects do not depend on background prevalence (no herd immunity dynamics), and vaccine-driven serotype or lineage replacement is not captured. The vaccine layer should therefore be interpreted as a simplified background modifier on acquisition risk rather than a full transmission model of vaccine ecology. The current implementation improves on the earlier dormant design by assigning vaccination to birth cohorts, but it still does not model herd effects, serotype replacement, waning, boosters, or catch-up campaigns.

5. **Broad regional groupings**: The model uses continental-level regions (e.g., "Europe", "Africa") rather than country-level or hospital-level variation. Antibiotic consumption patterns, testing capacity, pathogen mix, and resistance rates can vary dramatically between countries and institutions within the same region. The regional layer should therefore be read as a coarse structuring device for global comparisons, not as a substitute for country-specific or centre-specific epidemiology.

6. **No person-to-person transmission network**: Community infection rates are driven by organism-specific log-odds parameters calibrated to match observed incidence, not by direct contacts between simulated individuals. There is no explicit transmission network, no basic reproduction number (R₀), and no herd-immunity dynamic. Hospital acquisition is the one partial exception: nosocomial infection rates scale with the current hospital census, creating an implicit density-dependence within the inpatient population. The absence of a transmission model means the simulation cannot reproduce epidemic waves, outbreak amplification, or the impact of interventions — such as isolation, contact tracing, or infection-control procedures — that primarily work through blocking transmission chains. It also means community resistance prevalence is driven by selection, reversion, HGT, and calibrated acquisition rates rather than by strain spread from person to person. This is a deliberate trade-off: adding a full population-transmission layer for 42 organisms would require extensive additional parameterisation and would substantially increase runtime, while the primary policy questions addressed here (prescribing, stewardship, diagnostics, and access) are primarily mediated through selection pressure rather than transmission dynamics.

7. **Constant infection acquisition rates for most organisms**: With the exception of a global sanitation-improvement adjustment that raises all community-acquisition log-odds by approximately +1 in 1930 and declines to zero by 1950, bacterial infection acquisition rates are treated as constant over the simulation period. This is defensible for most of the 42 modelled organisms, where secular incidence trends have been small relative to the magnitude of resistance-selection effects. Two organisms are worth noting as caveats. First, *Helicobacter pylori*: birth-cohort prevalence in high-income countries was substantially higher in populations that grew up before 1960 (approximately 70–80%) than in those born after 1970 (approximately 30%), a consequence of improved sanitation and reduced intrafamilial transmission. Since clarithromycin-based triple therapy was introduced in the 1990s and most patients treated in that era were adults who had acquired *H. pylori* in childhood, the model's constant acquisition rate may underestimate the size of the treated pool that generated early clarithromycin and metronidazole resistance. Second, *Campylobacter jejuni*: community incidence in high-income countries rose approximately 4–6-fold between 1960 and 1990 as a consequence of intensification of industrial poultry farming (Blaser MJ, 1997). The model therefore applies nalidixic acid selection pressure to a *Campylobacter* population that is implicitly larger than it would have been in reality during the 1963–1975 period, which may slightly overestimate early quinolone-driven GyrA selection for this organism. Neither effect is expected to alter the model's resistance trajectories in a clinically meaningful way given the logarithmic relationship between drug-pressure volume and equilibrium resistance prevalence, but both are limitations to bear in mind when interpreting *H. pylori* clarithromycin resistance trends and the early phase of *Campylobacter* fluoroquinolone resistance.

8. **Simplified historical gonorrhoea incidence decline**: Gonorrhoea incidence in high-income countries was approximately 3-fold higher in the 1960s–1970s than in the present day, before sexual health programmes, contact tracing, and partner notification services substantially reduced transmission. The model now applies coarse historical acquisition multipliers for *N. gonorrhoeae* so that the pre-1980 and transitional pre-2000 eras generate a larger organism pool than the modern era, which better aligns the volume of penicillin, tetracycline, sulfonamide, and early fluoroquinolone treatment with the historical period that seeded the major resistance layers. This is still a compressed approximation rather than a full epidemiological reconstruction: the era multipliers are global and intentionally simple, so resistance trajectories for *N. gonorrhoeae* should still be interpreted as reflecting the broad timing and scale of historical selection rather than an exact quantitative reconstruction of country-specific gonorrhoea incidence.

9. **Exogenous resistance reservoirs are approximated rather than explicitly modelled**: The environmental/exogenous floor system (§7.7) applies background mechanism probabilities to the non-human fraction of community acquisition instead of explicitly simulating livestock populations, wastewater ecology, food-chain contamination, or other external reservoirs. In the current configuration this includes both classic agricultural pathways and a few explicitly parameterised non-agricultural exogenous pathways such as the rifampicin `rpoB` block. The floor values are therefore calibrated surrogates rather than direct measurements of source-reservoir prevalence. They can encode historically realistic persistence, but they do not respond endogenously to interventions targeted at those outside-human reservoirs unless the corresponding parameters are changed.

10. **No explicit background rifampicin exposure pathway for non-target organisms**: Outside directly modelled treatment-linked uses such as staphylococcal combination therapy and selected ICU combination regimens, the simulation does not currently place people on rifampicin simply because they are receiving TB treatment or another non-modelled indication. As a result, bystander selection of `rpoB` resistance in coincident gut or enteric carriage is only represented indirectly through the existing floor and persistence machinery, not through an explicit infection-independent rifampicin-use process. This can lead the model to understate rifampicin resistance prevalence for some non-target organisms, but that limitation is currently accepted because it is unlikely to materially alter the main stewardship-policy comparisons the model is being used for.

11. **Residual high-emergence-rate mechanisms are now treated as reviewed coarse-graining rather than open calibration defects**: The remaining bacteria–mechanism pairs with nominal `emergence_rate > 1.0` have been re-audited using the four-lever framework. Where the audit exposed a clean local miss, the model was changed: for example, gonorrhoea now has a simplified historical acquisition-volume correction, syphilis now has explicit historical macrolide alternative-treatment pressure, and *M. genitalium* no longer treats `erm_b` or `tet_m` as dominant pathways. The remaining high values are therefore not currently being interpreted as “obvious next bugs to fix”, but as coarse surrogate parameters sitting on top of other accepted abstractions already in the model, such as simplified exogenous reservoirs, coarse historical incidence trajectories, and omitted off-model drug exposures. Unless a future use case reveals another equally local, mechanistically defensible missing pathway, these residual high-emergence settings should be read as accepted limitations of the current policy-scale representation rather than as active pending calibration work.

---



## Appendix A — Bacteria, Drugs, Mechanisms and Enums

This appendix lists every entity in the model. Use it as a lookup reference when you encounter a specific bacterium, drug, or mechanism code in the main text.

The appendix is implementation-facing. Names, groupings, and enum labels are the simulation's internal vocabulary for representing major clinical categories; they are not meant to imply that every organism, drug, or ecological niche is exhaustively or uniquely represented by a single real-world classification scheme. They are included so that readers can see exactly how clinically familiar categories were operationalised inside a policy-scale simulation.



### A.1 Bacteria (42 species)

| Index | Species | Group | Carriage compartment |
|-------|---------|-------|---------------------|
| 0 | Acinetobacter baumannii | NonFermenter | Respiratory |
| 1 | Citrobacter spp. | Enterobacterales | Gut |
| 2 | Enterobacter spp. | Enterobacterales | Gut |
| 3 | Enterococcus faecalis | Streptococci | Gut |
| 4 | Enterococcus faecium | Streptococci | Gut |
| 5 | Escherichia coli | Enterobacterales | Gut |
| 6 | Klebsiella pneumoniae | Enterobacterales | Gut |
| 7 | Morganella spp. | Enterobacterales | Gut |
| 8 | Proteus spp. | Enterobacterales | Gut |
| 9 | Serratia spp. | Enterobacterales | Gut |
| 10 | Providencia stuartii | Enterobacterales | Genitourinary |
| 11 | Pseudomonas aeruginosa | NonFermenter | Respiratory |
| 12 | Stenotrophomonas maltophilia | NonFermenter | Respiratory |
| 13 | Staphylococcus aureus | Staphylococci | Skin/Soft Tissue |
| 14 | Staphylococcus epidermidis | Staphylococci | Skin/Soft Tissue |
| 15 | Streptococcus pneumoniae | Streptococci | Respiratory |
| 16 | Salmonella enterica serovar Typhi | Enterobacterales | Gut |
| 17 | Salmonella enterica serovar Paratyphi A | Enterobacterales | Gut |
| 18 | Invasive non-typhoidal Salmonella spp. | Enterobacterales | Gut |
| 19 | Shigella spp. | Enterobacterales | Gut |
| 20 | Neisseria gonorrhoeae | Fastidious | Genitourinary |
| 21 | Streptococcus pyogenes | Streptococci | Respiratory |
| 22 | Streptococcus agalactiae | Streptococci | Genitourinary |
| 23 | Haemophilus influenzae | Fastidious | Respiratory |
| 24 | Chlamydia trachomatis | Fastidious | Genitourinary |
| 25 | Mycoplasma genitalium | Fastidious | Genitourinary |
| 26 | Vibrio cholerae | EntericPathogen | Gut |
| 27 | Neisseria meningitidis | Fastidious | Respiratory |
| 28 | Listeria monocytogenes | Streptococci | Gut |
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



### A.2 Antibiotics (62 drugs)

The class labels in this table mirror the model's current internal `DrugClass` mapping rather than broader textbook umbrella categories.

| Drug | Class |
|------|-------|
| sulfanilamide | Sulfonamides |
| penicillin_g | Penicillins |
| ampicillin | Penicillins |
| amoxicillin | Penicillins |
| piperacillin | Penicillins |
| ticarcillin | Penicillins |
| flucloxacillin | Penicillins |
| cephalexin | Cephalosporins 1–2G |
| cefazolin | Cephalosporins 1–2G |
| cefuroxime | Cephalosporins 1–2G |
| ceftriaxone | Cephalosporins 3G |
| cefixime | Cephalosporins 3G |
| ceftazidime | Cephalosporins 3G |
| cefepime | Cephalosporins 4G |
| ceftaroline | Anti-MRSA Cephalosporins (5G) |
| ceftolozane_tazobactam | Cephalosporins 3G/BLI |
| cefiderocol | Siderophore Cephalosporins |
| meropenem | Carbapenems Group 2 |
| imipenem_c | Carbapenems Group 2 |
| ertapenem | Carbapenems Group 1 |
| aztreonam | Monobactams |
| erythromycin | Macrolides |
| azithromycin | Macrolides |
| clarithromycin | Macrolides |
| clindamycin | Lincosamides |
| gentamicin | Aminoglycosides Group 1 |
| tobramycin | Aminoglycosides Group 1 |
| amikacin | Aminoglycosides Group 2 |
| ciprofloxacin | Fluoroquinolones |
| levofloxacin | Fluoroquinolones |
| moxifloxacin | Fluoroquinolones |
| ofloxacin | Fluoroquinolones |
| tetracycline | Tetracyclines |
| doxycycline | Tetracyclines |
| minocycline | Tetracyclines |
| tigecycline | Glycylcyclines |
| vancomycin | Glycopeptides |
| teicoplanin | Lipoglycopeptides |
| dalbavancin | Lipoglycopeptides |
| linezolid | Oxazolidinones |
| tedizolid | Oxazolidinones |
| daptomycin | Lipopeptides |
| quinu_dalfo | Streptogramins |
| trim_sulf | Sulfonamides |
| chloramphenicol | Chloramphenicol |
| nitrofurantoin | Nitrofurans |
| fosfomycin | Phosphonic Acids |
| retapamulin | Pleuromutilins |
| fusidic_a | Steroid Antibacterials |
| metronidazole | Nitroimidazoles |
| fidaxomicin | Macrocycles |
| furazolidone | Nitrofurans |
| rifampicin | Rifamycins |
| amoxicillin_clavulanate | BLI Combinations |
| piperacillin_tazobactam | BLI Anti-Pseudomonal |
| ampicillin_sulbactam | BLI Sulbactam |
| ticarcillin_clavulanate | BLI Combinations |
| ceftazidime_avibactam | Ceftazidime-Avibactam |
| meropenem_vaborbactam | Meropenem-Vaborbactam |
| aztreonam_avibactam | Aztreonam-Avibactam |
| colistin | Polymyxins |



### A.3 Drug Classes (39 internal classes, mirroring the live `DrugClass` enum)

| Code | Enum variant | Meaning | Canonical drugs |
|------|--------------|---------|-----------------|
| `pen` | `Penicillins` | Penicillins | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin` |
| `bli` | `BliCombinations` | Beta-lactam/beta-lactamase inhibitor combinations | `amoxicillin_clavulanate`, `ticarcillin_clavulanate` |
| `bli_anti_pseudomonal` | `BliAntiPseudomonal` | Anti-pseudomonal beta-lactam/BLI combinations | `piperacillin_tazobactam` |
| `bli_sulbactam` | `BliSulbactam` | Sulbactam-containing beta-lactam/BLI combinations | `ampicillin_sulbactam` |
| `c1_2g` | `Cephalosporins1_2` | First- and second-generation cephalosporins | `cephalexin`, `cefazolin`, `cefuroxime` |
| `c3g` | `Cephalosporins3` | Third-generation cephalosporins | `ceftriaxone`, `ceftazidime`, `cefixime` |
| `c3g_bli` | `Cephalosporins3Bli` | Third-generation cephalosporin/BLI combinations | `ceftolozane_tazobactam` |
| `c4g` | `Cephalosporins4` | Fourth-generation cephalosporins | `cefepime` |
| `anti_mrsa_ceph` | `AntiMrsaCephalosporins` | Anti-MRSA cephalosporins | `ceftaroline` |
| `siderophore_ceph` | `SiderophoreCephalosporins` | Siderophore cephalosporins | `cefiderocol` |
| `cft_avi` | `CeftazidimeAvibactam` | Ceftazidime-avibactam class | `ceftazidime_avibactam` |
| `mer_vab` | `MeropenemVaborbactam` | Meropenem-vaborbactam class | `meropenem_vaborbactam` |
| `azt_avi` | `AztreonamAvibactam` | Aztreonam-avibactam class | `aztreonam_avibactam` |
| `carb_group1` | `CarbapenemsGroup1` | Carbapenems lacking non-fermenter coverage | `ertapenem` |
| `carb_group2` | `CarbapenemsGroup2` | Broad carbapenems with non-fermenter coverage | `meropenem`, `imipenem_c` |
| `mono` | `Monobactams` | Monobactams | `aztreonam` |
| `fq` | `Fluoroquinolones` | Fluoroquinolones | `ciprofloxacin`, `levofloxacin`, `moxifloxacin`, `ofloxacin` |
| `ag_group1` | `AminoglycosidesGroup1` | Aminoglycosides group 1 | `gentamicin`, `tobramycin` |
| `ag_group2` | `AminoglycosidesGroup2` | Aminoglycosides group 2 | `amikacin` |
| `mls` | `Macrolides` | Macrolides | `erythromycin`, `azithromycin`, `clarithromycin` |
| `lincosamides` | `Lincosamides` | Lincosamides | `clindamycin` |
| `glyc` | `Glycopeptides` | Glycopeptides | `vancomycin` |
| `lipoglycopeptides` | `Lipoglycopeptides` | Lipoglycopeptides | `teicoplanin`, `dalbavancin` |
| `tet` | `Tetracyclines` | Tetracyclines | `tetracycline`, `doxycycline`, `minocycline` |
| `glycylcyclines` | `Glycylcyclines` | Glycylcyclines | `tigecycline` |
| `poly` | `Polymyxins` | Polymyxins | `colistin` |
| `oxa` | `Oxazolidinones` | Oxazolidinones | `linezolid`, `tedizolid` |
| `chl` | `Chloramphenicol` | Chloramphenicol class | `chloramphenicol` |
| `sulf` | `Sulfonamides` | Sulfonamides | `sulfanilamide`, `trim_sulf` |
| `lipopeptides` | `Lipopeptides` | Lipopeptides | `daptomycin` |
| `streptogramins` | `Streptogramins` | Streptogramins | `quinu_dalfo` |
| `nitrofurans` | `Nitrofurans` | Nitrofurans | `nitrofurantoin`, `furazolidone` |
| `phosphonic_acids` | `PhosphonicAcids` | Phosphonic acids | `fosfomycin` |
| `nitroimidazoles` | `Nitroimidazoles` | Nitroimidazoles | `metronidazole` |
| `rifamycins` | `Rifamycins` | Rifamycins | `rifampicin` |
| `macrocycles` | `Macrocycles` | Macrocycles | `fidaxomicin` |
| `steroid_antibacterials` | `SteroidAntibacterials` | Steroid antibacterials | `fusidic_a` |
| `pleuromutilins` | `Pleuromutilins` | Pleuromutilins | `retapamulin` |
| `other` | `Other` | Fallback catch-all class | none currently in `DRUG_SHORT_NAMES`; used only if a future drug lacks an explicit mapping |



### A.4 Resistance Mechanisms (46)

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
| `AtInfectionTB` | Acquired at MDR-TB infection event; rifampicin resistance (`rpoB`) is pre-seeded deterministically because MDR-TB is by definition rifampicin-resistant |
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

This appendix is a documentation snapshot derived from the Rust configuration. The Rust source remains authoritative, and values should be regenerated or checked against `src/config.rs` before publication because calibration edits can move faster than this reference table.

### B.1 Global Parameters

Global parameters that govern cross-cutting model behaviour. Grouped thematically; each row gives the parameter name and its default value.

See: [§6.1 Treatment initiation](#61-treatment-initiation-deciding-to-start-antibiotics), [§6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use), [§6.3 Drug pharmacokinetics](#63-drug-pharmacokinetics), [§6.7 Drug toxicity](#67-drug-toxicity), [§2.4 Hospitalisation](#24-hospitalisation), [§2.5 Travel](#25-travel), [§4.3 Sepsis](#43-sepsis), [§7.3 Resistance emergence](#73-resistance-emergence), [§7.4 Resistance reversion](#74-resistance-reversion-and-fitness-costs), [§8 Microbiome and Carriage](#8-microbiome-and-carriage), [§9 Horizontal Gene Transfer](#9-horizontal-gene-transfer-hgt), [§10 Mortality](#10-mortality).

#### Treatment Initiation (logistic model)

| Parameter | Value |
| --- | ---: |
| antibiotic_initiation_base_log_odds | -5.5 |
| antibiotic_initiation_log_odds_symptomatic_infection | 6 |
| antibiotic_initiation_log_odds_test_identified | 0.92 |
| antibiotic_initiation_log_odds_already_on_drug | 0.18 |
| antibiotic_initiation_log_odds_immunodeficiency | 0.2 |
| antibiotic_initiation_log_odds_hospitalized | 0.7 |
| antibiotic_initiation_log_odds_sepsis | 6 |
| antibiotic_initiation_log_odds_no_indication | -1.05 |

#### Drug Activity and Cessation

| Parameter | Value |
| --- | ---: |
| drug_activity_to_bacteria_level_multiplier | 0.75 |
| drug_activity_slow_clearance_probability | 0.25 |
| drug_activity_slow_clearance_multiplier | 0.2 |
| double_dose_probability_if_identified_infection | 0.25 |
| random_drug_cessation_probability | 0.0045 |
| random_drug_cessation_probability_if_no_active_infection | 0.15 |
| antibiotic_infection_prevention_efficacy | 0.7 |

#### Drug Selection

| Parameter | Value |
| --- | ---: |
| minimal_potency_threshold_for_drug_selection | 0.15 |
| drug_selection_temperature | 0.55 |
| reserve_drug_score_penalty | 0.35 |

Here, `drug_selection_temperature` is the prescribing-variability setting described in Section 6.2: smaller values concentrate choices on the highest-scoring drug, while larger values spread choices across clinically plausible alternatives.

#### Treatment Failure and Restart

| Parameter | Value |
| --- | ---: |
| treatment_failure_assessment_day | 4 |
| treatment_failure_threshold | 0.5 |
| drug_failure_memory_days | 14 |
| restart_window_days | 5 |
| restart_bacteria_level_threshold | 1.5 |
| restart_window_probability | 0.3 |

For organisms where rescue therapy following first-line failure is incomplete in practice, a per-bacteria `{bacteria}_treatment_failure_no_second_line_probability` parameter (default 0.0) controls the fraction of treatment failures that do not receive any alternative drug. When this probability fires, the current (failed) drugs are stopped but no second-line agent is started, leaving the patient with a persistent resistant infection. Because the community mechanism cache records `mechanism_majority` profiles from active infections, these chronic carriers continuously seed the pool of circulating resistant strains at acquisition, generating long-run resistance prevalence without a microbiome carriage reservoir and without floor overrides.

*H. pylori* is the only organism currently configured with a non-zero value (`helicobacter_pylori_treatment_failure_no_second_line_probability = 0.80`). This reflects the substantial real-world gap between the proportion of clarithromycin-triple failures that are confirmed and the proportion that actually complete bismuth quadruple rescue therapy, driven by side-effect burden (nausea, metallic taste, darkening stools), cost, and loss to follow-up.

#### Hospitalization

| Parameter | Value |
| --- | ---: |
| hospitalization_base_log_odds | -10.4 |
| hospitalization_log_odds_per_age_year | 0.02 |
| hospitalization_log_odds_sepsis | 13 |
| hospitalization_log_odds_symptomatic_infection | 9.5 |
| hospitalization_symptomatic_infection_level_threshold | 3 |
| hospital_recovery_rate_per_day | 0.28 |
| hospital_max_days | 30 |
| hospital_prevent_discharge_with_sepsis | 1 |

#### Resistance Emergence and Decay

| Parameter | Value |
| --- | ---: |
| max_resistance_level | 1 |
| resistance_emergence_bacteria_level_multiplier | 9 |
| any_r_emergence_level_on_first_emergence | 0.5 |
| multi_drug_penalty_threshold_num_drugs | 2 |
| resistance_development_inhibition_single_drug | 0.05 |
| resistance_development_inhibition_partial_cross | 0.3 |
| mechanism_assignment_probability_on_any_r_gain | 0.8 |
| community_profile_cache_retention | 0.999 |
| hospital_profile_cache_retention | 0.999 |
| hospital_resistance_concentration_factor | per-bacteria (1.0–2.25; see Section 3.4) |
| run_pathway_reversion_rate_multiplier | 1 |
| mechanism_reversion_rate_global_multiplier | 1 |
| majority_r_memory_retention_per_day | 0.93 |

#### Microbiome Dynamics

| Parameter | Value |
| --- | ---: |
| microbiome_resistance_transfer_probability_per_day | 1e-4 |
| antibiotic_disruption_decay_half_life_days | 30 |
| microbiome_resistance_multiplier_on_acquisition | 0.5 |
| infection_from_microbiome_dampening | 0.7 |
| carriage_duration_log_odds_coefficient | -0.01 |
| carriage_duration_max_log_odds_effect | -2 |
| antibiotic_clearance_log_odds_per_unit_activity | 0.5 |
| carrier_resistance_inheritance_probability | 0.5 |
| community_resistance_dilution_factor | per-bacteria (0.03–1.00; see Section 3.4) |
| microbiome_majority_decay_half_life_days | 60 |
| microbiome_minority_decay_half_life_days | 18 |
| microbiome_majority_promotion_rate_per_day | 0.02 |

#### De Novo and HGT Multipliers

| Parameter | Value |
| --- | ---: |
| run_pathway_infection_de_novo_multiplier | 1 |
| run_pathway_microbiome_de_novo_multiplier | 1 |
| run_pathway_hgt_multiplier | 1 |
| run_pathway_reversion_rate_multiplier | 1 |
| run_pathway_microbiome_acquisition_multiplier | 1 |
| run_pathway_carrier_inheritance_multiplier | 1 |
| run_pathway_community_dilution_multiplier | 1 |
| run_pathway_microbiome_disruption_multiplier | 1 |
| counterfactual_resistance_multiplier | 1 by default; counterfactual branches set 0 |
| microbiome_de_novo_multiplier | 1 compatibility multiplier applied within the microbiome de novo pathway |

#### Horizontal Gene Transfer Modifiers

| Parameter | Value |
| --- | ---: |
| hgt_hospital_multiplier | 4 |
| hgt_antibiotic_pressure_multiplier | 1.5 |
| hgt_coinfection_multiplier | 1.25 |
| hgt_microbiome_only_penalty | 0.65 |
| hgt_gut_compartment_multiplier | 2 |
| hgt_minority_donor_multiplier | 0.2 |

#### Travel

| Parameter | Value |
| --- | ---: |
| travel_probability_per_day | 0.0002 |

#### Bacteria Growth Age Multipliers

| Parameter | Value |
| --- | ---: |
| bacteria_growth_age_multiplier_infant | 1.3 |
| bacteria_growth_age_multiplier_child | 1 |
| bacteria_growth_age_multiplier_adult | 1 |
| bacteria_growth_age_multiplier_elderly | 1.2 |
| bacteria_growth_immunodeficiency_multiplier | 1.5 |

#### Sepsis Onset

| Parameter | Value |
| --- | ---: |
| sepsis_minimum_duration_days | 1 |
| log_odds_sepsis_onset_immunosuppressed | 0.7 |
| log_odds_sepsis_onset_hospitalized | 0.5 |
| log_odds_sepsis_onset_not_under_care | 1 |
| log_odds_sepsis_onset_region_north_america | -0.5 |
| log_odds_sepsis_onset_region_europe | -0.6 |
| log_odds_sepsis_onset_region_oceania | -0.5 |
| log_odds_sepsis_onset_region_asia | -0.1 |
| log_odds_sepsis_onset_region_south_america | 0 |
| log_odds_sepsis_onset_region_africa | 0.1 |

#### Sepsis Recovery

| Parameter | Value |
| --- | ---: |
| sepsis_recovery_base_log_odds_per_day | 0 |
| sepsis_recovery_log_odds_bacteria_level | -0.3 |
| sepsis_recovery_log_odds_in_hospital | 0.8 |
| sepsis_recovery_log_odds_age_infant | -0.5 |
| sepsis_recovery_log_odds_age_child | 0.4 |
| sepsis_recovery_log_odds_age_adult | 0 |
| sepsis_recovery_log_odds_age_elderly | -0.7 |
| sepsis_recovery_log_odds_immunosuppressed | -1 |

#### Sepsis Death

| Parameter | Value |
| --- | ---: |
| sepsis_death_base_log_odds | -5 |
| sepsis_death_log_odds_age_infant | 1.1 |
| sepsis_death_log_odds_age_child | -0.7 |
| sepsis_death_log_odds_age_adult | 0 |
| sepsis_death_log_odds_age_elderly | 0.9 |
| sepsis_death_log_odds_immunosuppressed | 1.5 |
| sepsis_death_log_odds_bacteria_level | 0.35 |
| sepsis_death_log_odds_duration | 0.04 |
| sepsis_death_log_odds_early_phase | 0.8 |
| sepsis_death_early_phase_days | 3 |
| sepsis_death_log_odds_not_under_care | 1.4 |

#### Non-Sepsis Infection Mortality

| Parameter | Value |
| --- | ---: |
| infection_non_sepsis_base_log_odds | -9 |
| infection_non_sepsis_log_odds_per_level | 0 |
| infection_non_sepsis_log_odds_age_infant | 0 |
| infection_non_sepsis_log_odds_age_child | 0 |
| infection_non_sepsis_log_odds_age_adult | 0 |
| infection_non_sepsis_log_odds_age_elderly | 0 |
| infection_non_sepsis_log_odds_immunosuppressed | 0 |
| infection_non_sepsis_log_odds_in_hospital | 0 |
| infection_non_sepsis_minimum_bacteria_level | 0.5 |

#### Background Mortality

| Parameter | Value |
| --- | ---: |
| background_mortality_baseline_log_odds | -14.3 |
| mortality_baseline_1930_multiplier | 3 |
| mortality_baseline_2035_multiplier | 1 |
| mortality_improvement_half_life_years | 35 |
| log_odds_mortality_per_year_of_age | 0.055 |
| log_odds_mortality_per_year_of_age_squared | 0.008 |
| log_odds_mortality_immunosuppressed | 0.916 |
| log_odds_mortality_hospitalized | 0.262 |

#### Drug Toxicity

| Parameter | Value |
| --- | ---: |
| default_toxicity_reservoir_half_life_days | 1.5 |
| toxicity_age_multiplier_infant | 1.8 |
| toxicity_age_multiplier_child | 1.2 |
| toxicity_age_multiplier_adult | 1 |
| toxicity_age_multiplier_elderly | 2.2 |
| toxicity_immunosuppressed_multiplier | 2.5 |
| toxicity_hospital_multiplier | 1.3 |
| toxicity_discontinuation_threshold | 1e-5 |
| toxicity_discontinuation_avoidance_days | 30 |

#### Regional Resistance Scoring

| Parameter | Value |
| --- | ---: |
| regional_resistance_threshold_very_high | 0.6 |
| regional_resistance_threshold_high | 0.45 |
| regional_resistance_threshold_moderate | 0.1 |
| regional_resistance_penalty_very_high | 0.3 |
| regional_resistance_penalty_high | 0.5 |
| regional_resistance_penalty_moderate | 0.8 |

#### Therapy Scoring

| Parameter | Value |
| --- | ---: |
| targeted_therapy_narrow_spectrum_bonus | 5 |
| targeted_therapy_broad_spectrum_penalty | 0.1 |
| targeted_therapy_ineffective_drug_penalty | 0.001 |
| effective_potency_threshold_for_targeted_therapy | 0.1 |
| empiric_therapy_broad_spectrum_bonus | 0.85 |
| empiric_therapy_ineffective_penalty | 0.001 |

#### MDR-TB Era Multipliers

| Parameter | Value |
| --- | ---: |
| mdr_tb_pre_antibiotic_era_multiplier | 0 |
| mdr_tb_early_antibiotic_era_multiplier | 0 |
| mdr_tb_modern_era_multiplier | 1 |

### B.2 Drug Properties

Pharmacokinetic and clinical properties for each of the 62 modelled antimicrobial agents. The introduction time step is measured in days from 1 January 1930.

See: [§6.3 Drug pharmacokinetics](#63-drug-pharmacokinetics), [§6.5 Drug potency matrix](#65-drug-potency-matrix), [§6.6 Drug availability](#66-drug-availability-by-region-and-era), [§6.7 Drug toxicity](#67-drug-toxicity), [§6.8 Antibiotic infection prevention](#68-antibiotic-infection-prevention).

| Drug | Class | Intro (days) | Init level | t½ (days) | 2× dose mult | Spectrum | Tox hazard | Tox t½ (days) | Microbiome disrupt |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| sulfanilamide | folate_antagonists | 2555 | 10 | 0.45 | 2 | 3 | 0 | 1.5 | 0.3 |
| penicillin_g | penicillins | 3555 | 10 | 0.04 | 2 | 2 | 5e-10 | 1.5 | 0.3 |
| ampicillin | penicillins | 11315 | 10 | 0.04 | 2 | 3 | 5e-10 | 1.5 | 0.3 |
| amoxicillin | penicillins | 13780 | 10 | 0.04 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| piperacillin | penicillins | 16065 | 10 | 0.04 | 2 | 3 | 5e-10 | 1.5 | 0.3 |
| ticarcillin | penicillins | 14600 | 10 | 0.046 | 2 | 3 | 0 | 1.5 | 0.3 |
| cephalexin | cephalosporins_1_2 | 14605 | 10 | 0.04 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| cefazolin | cephalosporins_1_2 | 15700 | 10 | 0.08 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| cefuroxime | cephalosporins_1_2 | 17525 | 10 | 0.05 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| ceftriaxone | cephalosporins_3_4 | 19715 | 10 | 0.33 | 2 | 4 | 5e-10 | 1.5 | 0.3 |
| ceftazidime | cephalosporins_3_4 | 20080 | 10 | 0.08 | 2 | 3 | 4e-10 | 1.5 | 0.3 |
| cefepime | cephalosporins_3_4 | 24195 | 10 | 0.08 | 2 | 4 | 1e-9 | 1.5 | 0.3 |
| ceftaroline | cephalosporins_3_4 | 29305 | 10 | 0.11 | 2 | 3 | 5e-10 | 1.5 | 0.3 |
| ceftolozane_tazobactam | cephalosporins_3_4 | 30295 | 10 | 0.125 | 2 | 3 | 0 | 1.5 | 0.3 |
| cefiderocol | unknown | 33510 | 10 | 0.1 | 2 | 3 | 0 | 1.5 | 0.3 |
| meropenem | carbapenems | 24195 | 10 | 0.04 | 2 | 5 | 6e-10 | 1.5 | 0.3 |
| imipenem_c | carbapenems | 20080 | 10 | 0.04 | 2 | 3 | 1e-9 | 1.5 | 0.3 |
| ertapenem | carbapenems | 25920 | 10 | 0.17 | 2 | 3 | 6e-10 | 1.5 | 0.3 |
| aztreonam | cephalosporins_3_4 | 20445 | 10 | 0.08 | 2 | 3 | 0 | 1.5 | 0.3 |
| erythromycin | macrolides | 8025 | 10 | 0.08 | 2 | 3 | 2e-9 | 1.5 | 0.3 |
| azithromycin | macrolides | 22260 | 10 | 2.8 | 2 | 4 | 1.5e-9 | 1.5 | 0.3 |
| clarithromycin | macrolides | 21895 | 10 | 0.25 | 2 | 3 | 1.5e-9 | 1.5 | 0.3 |
| clindamycin | macrolides | 13870 | 10 | 0.125 | 2 | 3 | 1e-9 | 1.5 | 0.3 |
| gentamicin | aminoglycosides | 12045 | 10 | 0.08 | 2 | 3 | 1.5e-8 | 1.5 | 0.3 |
| tobramycin | aminoglycosides | 16325 | 10 | 0.08 | 2 | 3 | 1.3e-8 | 1.5 | 0.3 |
| amikacin | aminoglycosides | 16690 | 10 | 0.08 | 2 | 3 | 1.7e-8 | 1.5 | 0.3 |
| ciprofloxacin | fluoroquinolones | 20805 | 10 | 0.17 | 2 | 4.5 | 3e-9 | 1.5 | 0.3 |
| levofloxacin | fluoroquinolones | 24195 | 10 | 0.33 | 2 | 3 | 3e-9 | 1.5 | 0.3 |
| moxifloxacin | fluoroquinolones | 25290 | 10 | 0.5 | 2 | 3 | 5e-9 | 1.5 | 0.3 |
| ofloxacin | fluoroquinolones | 21895 | 10 | 0.25 | 2 | 3 | 3e-9 | 1.5 | 0.3 |
| tetracycline | tetracyclines | 6575 | 10 | 0.33 | 2 | 3 | 1e-9 | 1.5 | 0.3 |
| doxycycline | tetracyclines | 13505 | 10 | 0.75 | 2 | 3 | 1e-9 | 1.5 | 0.3 |
| minocycline | tetracyclines | 14965 | 10 | 0.67 | 2 | 3 | 1.5e-9 | 1.5 | 0.3 |
| tigecycline | unknown | 28040 | 10 | 1.75 | 2 | 3 | 0 | 1.5 | 0.3 |
| vancomycin | glycopeptides | 10215 | 10 | 0.25 | 2 | 2.5 | 6e-9 | 1.5 | 0.3 |
| teicoplanin | lipoglycopeptides | 21170 | 10 | 3.5 | 2 | 3 | 0 | 1.5 | 0.3 |
| dalbavancin | lipoglycopeptides | 30660 | 10 | 10 | 2 | 3 | 0 | 1.5 | 0.3 |
| linezolid | oxazolidinones | 25550 | 10 | 0.21 | 2 | 2 | 8e-9 | 1.5 | 0.3 |
| tedizolid | oxazolidinones | 30660 | 10 | 0.5 | 2 | 3 | 4e-9 | 1.5 | 0.3 |
| daptomycin | unknown | 27375 | 10 | 0.33 | 2 | 3 | 0 | 1.5 | 0.3 |
| quinu_dalfo | unknown | 25290 | 10 | 0.5 | 2 | 3 | 0 | 1.5 | 0.3 |
| trim_sulf | folate_antagonists | 13870 | 10 | 0.5 | 2 | 3.5 | 2e-9 | 1.5 | 0.3 |
| chloramphenicol | unknown | 6935 | 10 | 0.125 | 2 | 3 | 1e-8 | 1.5 | 0.3 |
| nitrofurantoin | unknown | 8395 | 10 | 0.017 | 2 | 1 | 3e-9 | 1.5 | 0.3 |
| fosfomycin | unknown | 10590 | 10 | 0.15 | 2 | 3 | 0 | 1.5 | 0.3 |
| retapamulin | unknown | 28405 | 10 | 0.25 | 2 | 3 | 0 | 1.5 | 0.3 |
| fusidic_a | unknown | 11680 | 10 | 0.375 | 2 | 3 | 0 | 1.5 | 0.3 |
| metronidazole | nitroimidazoles | 10965 | 10 | 0.33 | 2 | 3 | 2e-9 | 1.5 | 0.3 |
| fidaxomicin | unknown | 29565 | 10 | 0.5 | 2 | 3 | 0 | 1.5 | 0.3 |
| furazolidone | unknown | 9125 | 10 | 0.25 | 2 | 3 | 0 | 1.5 | 0.3 |
| rifampicin | unknown | 13140 | 10 | 0.25 | 2 | 3 | 4e-9 | 1.5 | 0.3 |
| amoxicillin_clavulanate | penicillins | 16425 | 10 | 0.04 | 2 | 3 | 6e-10 | 1.5 | 0.3 |
| piperacillin_tazobactam | penicillins | 19715 | 10 | 0.04 | 2 | 3 | 6e-10 | 1.5 | 0.3 |
| ampicillin_sulbactam | penicillins | 18250 | 10 | 0.04 | 2 | 3 | 0 | 1.5 | 0.3 |
| ticarcillin_clavulanate | penicillins | 18250 | 10 | 0.046 | 2 | 3 | 0 | 1.5 | 0.3 |
| ceftazidime_avibactam | cephalosporins_3_4 | 27740 | 10 | 0.08 | 2 | 3 | 0 | 1.5 | 0.3 |
| meropenem_vaborbactam | carbapenems | 32045 | 10 | 0.04 | 2 | 3 | 0 | 1.5 | 0.3 |
| colistin | polymyxins | 8020 | 10 | 0.08 | 2 | 4 | 2.5e-8 | 1.5 | 0.3 |
| flucloxacillin | penicillins | 14600 | 10 | 0.04 | 2 | 1.6 | 1e-8 | 1.5 | 0.3 |
| aztreonam_avibactam | cephalosporins_3_4 | 34675 | 10 | 0.08 | 2 | 3 | 0 | 1.5 | 0.3 |
| cefixime | cephalosporins_3_4 | 21535 | 10 | 0.125 | 2 | 2.8 | 5e-10 | 1.5 | 0.3 |
| nalidixic_acid | fluoroquinolones | 12045 | 10 | 0.08 | 2 | 2 | 4e-9 | 1.5 | 0.3 |

### B.3 Bacteria Properties

Per-bacteria parameters governing acquisition, growth, symptom onset, and clinical outcomes for each of the 42 bacterial species. Resistance-dilution and hospital-concentration settings are discussed separately in Section 3.4 because they are calibrated as part of the acquisition-reservoir logic rather than the core clinical trajectory table below.

See: [§3.1 Community acquisition](#31-community-acquisition), [§4.2 Infection dynamics](#42-infection-dynamics), [§4.3 Sepsis](#43-sepsis), [§4.4 Natural clearance](#44-natural-clearance), [§8.1 Carriage compartments](#81-carriage-compartments).

| Bacteria | Acq log-odds | Init level | Delta level/day | Max level | Microb clr/day | Microb vs inf | Drug cess prob | Sx threshold | Sx delay (d) | Sepsis log-odds | Mech-less rev rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| acinetobacter_baumannii | -18.7 | 0.01 | 0.55 | 5 | 0.1 | 8 | 0.0075 | 0.5 | 1 | -7.1 | 4e-4 |
| citrobacter_spp. | -16.3 | 0.01 | 0.5 | 5 | 0.08 | 11.3 | 0.0045 | 0.5 | 1 | -8.6 | 4e-4 |
| enterobacter_spp. | -16.8 | 0.01 | 0.5 | 5 | 0.07 | 11.5 | 0.0045 | 0.5 | 1 | -7.1 | 4e-4 |
| enterococcus_faecalis | -17.1 | 0.01 | 0.48 | 5 | 0.008 | 11.5 | 0.0075 | 0.5 | 1 | -6.9 | 4e-4 |
| enterococcus_faecium | -17.6 | 0.01 | 0.48 | 5 | 0.06 | 13 | 0.0075 | 0.5 | 1 | -6.4 | 4e-4 |
| escherichia_coli | -11.6 | 0.01 | 0.5 | 5 | 0.005 | 6 | 0.025 | 0.5 | 1 | -8.9 | 4e-4 |
| klebsiella_pneumoniae | -14.7 | 0.01 | 0.52 | 5 | 0.03 | 9.2 | 0.0075 | 0.5 | 1 | -6.9 | 4e-4 |
| morganella_spp. | -17.2 | 0.01 | 0.48 | 5 | 0.1 | 10 | 0.0045 | 0.5 | 1 | -7.2 | 4e-4 |
| proteus_spp. | -16.1 | 0.01 | 0.5 | 5 | 0.08 | 8.5 | 0.0045 | 0.5 | 1 | -7.2 | 4e-4 |
| serratia_spp. | -17.3 | 0.01 | 0.48 | 5 | 0.1 | 10 | 0.0045 | 0.5 | 1 | -7.4 | 4e-4 |
| p_stuartii | -17.5 | 0.01 | 0.5 | 5 | 0.09 | 8.7 | 0.0045 | 0.75 | 1 | -12 | 4e-4 |
| pseudomonas_aeruginosa | -16 | 0.01 | 0.55 | 5 | 0.12 | 9 | 0.0075 | 0.8 | 1 | -5.9 | 4e-4 |
| stenotrophomonas_maltophilia | -19 | 0.01 | 0.45 | 5 | 0.06 | 7 | 0.0045 | 0.9 | 2.5 | -7.4 | 4e-4 |
| staphylococcus_aureus | -12.9 | 0.01 | 0.6 | 5 | 0.05 | 7.5 | 0.015 | 0.5 | 1 | -7.4 | 4e-4 |
| staphylococcus_epidermidis | -16 | 0.01 | 0.35 | 4 | 0.015 | 13.5 | 0.0045 | 1 | 3 | -7.4 | 4e-4 |
| streptococcus_pneumoniae | -12.6 | 0.01 | 0.6 | 5 | 0.05 | 7 | 0.015 | 0.5 | 1 | -9.9 | 4e-4 |
| salmonella_enterica_serovar_typhi | -17.3 | 0.01 | 0.45 | 5 | 0.003 | -7 | 0.0045 | 0.5 | 1 | -7.4 | 4e-4 |
| salmonella_enterica_serovar_paratyphi_a | -16.8 | 0.01 | 0.45 | 5 | 0.15 | -0.5 | 0.0045 | 0.5 | 1 | -8.4 | 4e-4 |
| invasive_non-typhoidal_salmonella_spp. | -17.8 | 0.01 | 0.5 | 5 | 0.12 | 3.2 | 0.0045 | 0.5 | 1 | -8.6 | 4e-4 |
| shigella_spp. | -12.6 | 0.01 | 0.55 | 5 | 0.15 | -0.5 | 0.0045 | 0.5 | 1 | -11.4 | 4e-4 |
| neisseria_gonorrhoeae | -13.5 | 0.01 | 0.55 | 5 | 0.2 | 3 | 0.0045 | 0.5 | 1 | -23.4 | 4e-4 |
| streptococcus_pyogenes | -14.4 | 0.01 | 0.7 | 5 | 0.08 | 8 | 0.015 | 0.5 | 1 | -5.9 | 4e-4 |
| streptococcus_agalactiae | -15.9 | 0.01 | 0.52 | 5 | 0.06 | 11 | 0.015 | 0.5 | 1 | -6.4 | 4e-4 |
| haemophilus_influenzae | -17.4 | 0.01 | 0.55 | 5 | 0.06 | 12.5 | 0.015 | 0.5 | 1 | -8.6 | 4e-4 |
| chlamydia_trachomatis | -12.8 | 0.01 | 0.25 | 5 | 0.2 | 4.5 | 0.007 | 0.8 | 1 | -18.4 | 4e-4 |
| mycoplasma_genitalium | -12.1 | 0.01 | 0.28 | 5 | 0.18 | 4.7 | 0.0045 | 0.9 | 5 | -12 | 4e-4 |
| vibrio_cholerae | -18.65 | 0.01 | 0.7 | 5 | 0.15 | 1 | 0.025 | 0.5 | 1 | -8.4 | 4e-4 |
| neisseria_meningitidis | -18.5 | 0.01 | 0.65 | 5 | 0.05 | 10.5 | 0.01 | 3 | 1 | -7.3 | 4e-4 |
| listeria_monocytogenes | -19 | 0.01 | 0.25 | 5 | 0.1 | 12.5 | 0.0045 | 0.5 | 1 | -7.4 | 4e-4 |
| clostridioides_difficile | -15.15 | 0.01 | 0.55 | 5 | 0.02 | 7 | 0.005 | 0.5 | 1 | -10.4 | 4e-4 |
| bacteroides_fragilis | -15.1 | 0.01 | 0.42 | 5 | 0.004 | 11 | 0.0045 | 1.2 | 2 | -12 | 4e-4 |
| campylobacter_jejuni | -13 | 0.01 | 0.52 | 5 | 0.12 | 2.5 | 0.015 | 0.5 | 1 | -19.4 | 4e-4 |
| enterobacter_cloacae | -16.8 | 0.01 | 0.5 | 5 | 0.04 | 13 | 0.0045 | 0.5 | 1 | -7.2 | 4e-4 |
| yersinia_enterocolitica | -16.6 | 0.01 | 0.45 | 5 | 0.25 | 5.5 | 0.0045 | 0.5 | 1 | -8.9 | 4e-4 |
| moraxella_catarrhalis | -14.6 | 0.01 | 0.55 | 5 | 0.05 | 11 | 0.0045 | 2 | 1 | -10.2 | 4e-4 |
| treponema_pallidum | -12.7 | 0.01 | 0.18 | 5 | 0.35 | 5.5 | 0.0045 | 0.6 | 1 | -10.4 | 4e-4 |
| bordetella_pertussis | -12.85 | 0.01 | 0.42 | 5 | 0.2 | 2.5 | 0.0075 | 0.5 | 1 | -10.4 | 4e-4 |
| helicobacter_pylori | -13.5 | 0.01 | 0.2 | 5 | 0.001 | 6.65 | 0.005 | 1.5 | 30 | -250 | 4e-4 |
| mdr_mycobacterium_tuberculosis | -16.9 | 0.01 | 0.15 | 5 | 0.0015 | -1 | 6e-4 | 2 | 1 | -38 | 4e-4 |
| mycoplasma_pneumoniae | -12 | 0.01 | 0.35 | 5 | 0.01 | 0.5 | 0.015 | 0.5 | 1 | -12 | 4e-4 |
| legionella_pneumophila | -15.5 | 0.01 | 0.55 | 5 | 0.01 | -2 | 0.0085 | 0.5 | 1 | -12 | 4e-4 |
| burkholderia_cepacia_complex | -20 | 0.01 | 0.45 | 5 | 0.01 | 0.5 | 0.0075 | 0.5 | 1 | -12 | 4e-4 |

### B.4 Drug–Bacteria Potency Matrix

Baseline potency (MIC-derived effectiveness when no resistance is present) and initiation multiplier (stewardship weighting for drug selection) for each drug–bacteria pair. 42 bacteria × 62 drugs = 2604 entries.

See: [§6.5 Drug potency matrix](#65-drug-potency-matrix), [§6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use).

| Bacteria | Drug | Potency (no R) | Init multiplier |
| --- | ---: | ---: | ---: |
| acinetobacter_baumannii | sulfanilamide | 0.1 | 1 |
| acinetobacter_baumannii | penicillin_g | 0.05 | 0.01 |
| acinetobacter_baumannii | ampicillin | 0.05 | 0.01 |
| acinetobacter_baumannii | amoxicillin | 0.05 | 0.01 |
| acinetobacter_baumannii | piperacillin | 0.6 | 1 |
| acinetobacter_baumannii | ticarcillin | 0.5 | 1 |
| acinetobacter_baumannii | cephalexin | 0.05 | 1 |
| acinetobacter_baumannii | cefazolin | 0.05 | 1 |
| acinetobacter_baumannii | cefuroxime | 0.1 | 1 |
| acinetobacter_baumannii | ceftriaxone | 0.1 | 1 |
| acinetobacter_baumannii | ceftazidime | 0.6 | 1 |
| acinetobacter_baumannii | cefepime | 0.7 | 1 |
| acinetobacter_baumannii | ceftaroline | 0.1 | 1 |
| acinetobacter_baumannii | ceftolozane_tazobactam | 0.1 | 1 |
| acinetobacter_baumannii | cefiderocol | 0.1 | 1 |
| acinetobacter_baumannii | meropenem | 0.85 | 0.005 |
| acinetobacter_baumannii | imipenem_c | 0.8 | 0.005 |
| acinetobacter_baumannii | ertapenem | 0.1 | 0.005 |
| acinetobacter_baumannii | aztreonam | 0.1 | 1 |
| acinetobacter_baumannii | gentamicin | 0.75 | 1 |
| acinetobacter_baumannii | tobramycin | 0.7 | 1 |
| acinetobacter_baumannii | amikacin | 0.8 | 1 |
| acinetobacter_baumannii | ciprofloxacin | 0.7 | 1 |
| acinetobacter_baumannii | levofloxacin | 0.7 | 1 |
| acinetobacter_baumannii | moxifloxacin | 0.6 | 1 |
| acinetobacter_baumannii | ofloxacin | 0.6 | 1 |
| acinetobacter_baumannii | tetracycline | 0.6 | 1 |
| acinetobacter_baumannii | doxycycline | 0.7 | 1 |
| acinetobacter_baumannii | minocycline | 0.8 | 1 |
| acinetobacter_baumannii | tigecycline | 0.1 | 1 |
| acinetobacter_baumannii | dalbavancin | 0 | 0.005 |
| acinetobacter_baumannii | linezolid | 0 | 0.005 |
| acinetobacter_baumannii | tedizolid | 0 | 0.005 |
| acinetobacter_baumannii | daptomycin | 0.1 | 1 |
| acinetobacter_baumannii | quinu_dalfo | 0 | 0.005 |
| acinetobacter_baumannii | trim_sulf | 0.6 | 1 |
| acinetobacter_baumannii | chloramphenicol | 0.7 | 1 |
| acinetobacter_baumannii | nitrofurantoin | 0.1 | 1 |
| acinetobacter_baumannii | fosfomycin | 0.1 | 1 |
| acinetobacter_baumannii | fidaxomicin | 0.1 | 1 |
| acinetobacter_baumannii | furazolidone | 0.1 | 1 |
| acinetobacter_baumannii | rifampicin | 0.6 | 1 |
| acinetobacter_baumannii | amoxicillin_clavulanate | 0.05 | 1 |
| acinetobacter_baumannii | piperacillin_tazobactam | 0.7 | 1 |
| acinetobacter_baumannii | ampicillin_sulbactam | 0.7 | 1 |
| acinetobacter_baumannii | ticarcillin_clavulanate | 0.6 | 1 |
| acinetobacter_baumannii | ceftazidime_avibactam | 0.7 | 0.005 |
| acinetobacter_baumannii | meropenem_vaborbactam | 0.8 | 0.005 |
| acinetobacter_baumannii | colistin | 0.9 | 0.005 |
| acinetobacter_baumannii | flucloxacillin | 0.01 | 1 |
| acinetobacter_baumannii | aztreonam_avibactam | 0.1 | 1 |
| acinetobacter_baumannii | cefixime | 0.1 | 1 |
| citrobacter_spp. | sulfanilamide | 0.5 | 1 |
| citrobacter_spp. | penicillin_g | 0.1 | 1 |
| citrobacter_spp. | ampicillin | 0.1 | 1 |
| citrobacter_spp. | amoxicillin | 0.1 | 1 |
| citrobacter_spp. | piperacillin | 0.8 | 1 |
| citrobacter_spp. | ticarcillin | 0.75 | 1 |
| citrobacter_spp. | cephalexin | 0.1 | 1 |
| citrobacter_spp. | cefazolin | 0.1 | 1 |
| citrobacter_spp. | cefuroxime | 0.8 | 1 |
| citrobacter_spp. | ceftriaxone | 0.85 | 1 |
| citrobacter_spp. | ceftazidime | 0.8 | 1 |
| citrobacter_spp. | cefepime | 0.9 | 1 |
| citrobacter_spp. | ceftaroline | 0.6 | 1 |
| citrobacter_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| citrobacter_spp. | cefiderocol | 0.8 | 1 |
| citrobacter_spp. | meropenem | 0.95 | 0.005 |
| citrobacter_spp. | imipenem_c | 0.95 | 0.005 |
| citrobacter_spp. | ertapenem | 0.9 | 0.005 |
| citrobacter_spp. | aztreonam | 0.85 | 1 |
| citrobacter_spp. | erythromycin | 0.1 | 1 |
| citrobacter_spp. | azithromycin | 0.1 | 1 |
| citrobacter_spp. | clarithromycin | 0.1 | 1 |
| citrobacter_spp. | clindamycin | 0.1 | 1 |
| citrobacter_spp. | gentamicin | 0.85 | 1 |
| citrobacter_spp. | tobramycin | 0.8 | 1 |
| citrobacter_spp. | amikacin | 0.9 | 1 |
| citrobacter_spp. | ciprofloxacin | 0.9 | 1 |
| citrobacter_spp. | levofloxacin | 0.85 | 1 |
| citrobacter_spp. | moxifloxacin | 0.7 | 1 |
| citrobacter_spp. | ofloxacin | 0.8 | 1 |
| citrobacter_spp. | tetracycline | 0.8 | 1 |
| citrobacter_spp. | doxycycline | 0.85 | 1 |
| citrobacter_spp. | minocycline | 0.85 | 1 |
| citrobacter_spp. | tigecycline | 0.1 | 1 |
| citrobacter_spp. | vancomycin | 0.1 | 1 |
| citrobacter_spp. | teicoplanin | 0.1 | 1 |
| citrobacter_spp. | dalbavancin | 0.1 | 0.005 |
| citrobacter_spp. | linezolid | 0.1 | 0.005 |
| citrobacter_spp. | tedizolid | 0.1 | 0.005 |
| citrobacter_spp. | daptomycin | 0.1 | 1 |
| citrobacter_spp. | quinu_dalfo | 0.1 | 0.005 |
| citrobacter_spp. | trim_sulf | 0.9 | 1 |
| citrobacter_spp. | chloramphenicol | 0.85 | 1 |
| citrobacter_spp. | nitrofurantoin | 0.8 | 1 |
| citrobacter_spp. | fosfomycin | 0.1 | 1 |
| citrobacter_spp. | retapamulin | 0.05 | 1 |
| citrobacter_spp. | fusidic_a | 0.05 | 1 |
| citrobacter_spp. | metronidazole | 0.05 | 1 |
| citrobacter_spp. | fidaxomicin | 0.1 | 1 |
| citrobacter_spp. | furazolidone | 0.1 | 1 |
| citrobacter_spp. | rifampicin | 0.7 | 1 |
| citrobacter_spp. | amoxicillin_clavulanate | 0.9 | 1 |
| citrobacter_spp. | piperacillin_tazobactam | 0.9 | 1 |
| citrobacter_spp. | ampicillin_sulbactam | 0.85 | 1 |
| citrobacter_spp. | ticarcillin_clavulanate | 0.8 | 1 |
| citrobacter_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| citrobacter_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| citrobacter_spp. | colistin | 0.7 | 0.005 |
| citrobacter_spp. | flucloxacillin | 0.01 | 1 |
| citrobacter_spp. | aztreonam_avibactam | 1 | 1 |
| citrobacter_spp. | cefixime | 0.8 | 1 |
| enterobacter_spp. | sulfanilamide | 0.5 | 1 |
| enterobacter_spp. | penicillin_g | 0.1 | 1 |
| enterobacter_spp. | ampicillin | 0.1 | 1 |
| enterobacter_spp. | amoxicillin | 0.1 | 1 |
| enterobacter_spp. | piperacillin | 0.75 | 1 |
| enterobacter_spp. | ticarcillin | 0.7 | 1 |
| enterobacter_spp. | cephalexin | 0.1 | 1 |
| enterobacter_spp. | cefazolin | 0.1 | 1 |
| enterobacter_spp. | cefuroxime | 0.6 | 1 |
| enterobacter_spp. | ceftriaxone | 0.5 | 1 |
| enterobacter_spp. | ceftazidime | 0.8 | 1 |
| enterobacter_spp. | cefepime | 0.85 | 1 |
| enterobacter_spp. | ceftaroline | 0.4 | 1 |
| enterobacter_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| enterobacter_spp. | cefiderocol | 0.8 | 1 |
| enterobacter_spp. | meropenem | 0.95 | 0.005 |
| enterobacter_spp. | imipenem_c | 0.95 | 0.005 |
| enterobacter_spp. | ertapenem | 0.9 | 0.005 |
| enterobacter_spp. | aztreonam | 0.8 | 1 |
| enterobacter_spp. | gentamicin | 0.85 | 1 |
| enterobacter_spp. | tobramycin | 0.8 | 1 |
| enterobacter_spp. | amikacin | 0.9 | 1 |
| enterobacter_spp. | ciprofloxacin | 0.9 | 1 |
| enterobacter_spp. | levofloxacin | 0.85 | 1 |
| enterobacter_spp. | moxifloxacin | 0.7 | 1 |
| enterobacter_spp. | ofloxacin | 0.8 | 1 |
| enterobacter_spp. | tetracycline | 0.8 | 1 |
| enterobacter_spp. | doxycycline | 0.85 | 1 |
| enterobacter_spp. | minocycline | 0.85 | 1 |
| enterobacter_spp. | tigecycline | 0.1 | 1 |
| enterobacter_spp. | dalbavancin | 0 | 0.005 |
| enterobacter_spp. | linezolid | 0 | 0.005 |
| enterobacter_spp. | tedizolid | 0 | 0.005 |
| enterobacter_spp. | daptomycin | 0.1 | 1 |
| enterobacter_spp. | quinu_dalfo | 0 | 0.005 |
| enterobacter_spp. | trim_sulf | 0.85 | 1 |
| enterobacter_spp. | chloramphenicol | 0.8 | 1 |
| enterobacter_spp. | nitrofurantoin | 0.7 | 1 |
| enterobacter_spp. | fosfomycin | 0.1 | 1 |
| enterobacter_spp. | fidaxomicin | 0.1 | 1 |
| enterobacter_spp. | furazolidone | 0.1 | 1 |
| enterobacter_spp. | rifampicin | 0.6 | 1 |
| enterobacter_spp. | amoxicillin_clavulanate | 0.7 | 1 |
| enterobacter_spp. | piperacillin_tazobactam | 0.85 | 1 |
| enterobacter_spp. | ampicillin_sulbactam | 0.7 | 1 |
| enterobacter_spp. | ticarcillin_clavulanate | 0.8 | 1 |
| enterobacter_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| enterobacter_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| enterobacter_spp. | colistin | 0.7 | 0.005 |
| enterobacter_spp. | flucloxacillin | 0.01 | 1 |
| enterobacter_spp. | aztreonam_avibactam | 1 | 1 |
| enterobacter_spp. | cefixime | 0.8 | 1 |
| enterococcus_faecalis | sulfanilamide | 0.1 | 1 |
| enterococcus_faecalis | penicillin_g | 0.8 | 1 |
| enterococcus_faecalis | ampicillin | 0.9 | 1 |
| enterococcus_faecalis | amoxicillin | 0.9 | 1 |
| enterococcus_faecalis | piperacillin | 0.75 | 1 |
| enterococcus_faecalis | ticarcillin | 0.7 | 1 |
| enterococcus_faecalis | cephalexin | 0.1 | 1 |
| enterococcus_faecalis | cefazolin | 0.1 | 1 |
| enterococcus_faecalis | cefuroxime | 0.1 | 1 |
| enterococcus_faecalis | ceftriaxone | 0.1 | 1 |
| enterococcus_faecalis | ceftazidime | 0.1 | 1 |
| enterococcus_faecalis | cefepime | 0.1 | 1 |
| enterococcus_faecalis | ceftaroline | 0.1 | 1 |
| enterococcus_faecalis | ceftolozane_tazobactam | 0.05 | 1 |
| enterococcus_faecalis | cefiderocol | 0.05 | 1 |
| enterococcus_faecalis | meropenem | 0.7 | 0.005 |
| enterococcus_faecalis | imipenem_c | 0.7 | 0.005 |
| enterococcus_faecalis | ertapenem | 0.1 | 0.005 |
| enterococcus_faecalis | erythromycin | 0.7 | 0.01 |
| enterococcus_faecalis | azithromycin | 0.7 | 0.01 |
| enterococcus_faecalis | clarithromycin | 0.7 | 1 |
| enterococcus_faecalis | clindamycin | 0.7 | 1 |
| enterococcus_faecalis | gentamicin | 0.1 | 1 |
| enterococcus_faecalis | tobramycin | 0.1 | 1 |
| enterococcus_faecalis | amikacin | 0.1 | 1 |
| enterococcus_faecalis | ciprofloxacin | 0.7 | 1 |
| enterococcus_faecalis | levofloxacin | 0.7 | 1 |
| enterococcus_faecalis | moxifloxacin | 0.7 | 1 |
| enterococcus_faecalis | ofloxacin | 0.7 | 1 |
| enterococcus_faecalis | tetracycline | 0.8 | 1 |
| enterococcus_faecalis | doxycycline | 0.8 | 1 |
| enterococcus_faecalis | minocycline | 0.85 | 1 |
| enterococcus_faecalis | tigecycline | 0.1 | 1 |
| enterococcus_faecalis | vancomycin | 0.95 | 4 |
| enterococcus_faecalis | teicoplanin | 0.9 | 1 |
| enterococcus_faecalis | dalbavancin | 0.9 | 0.005 |
| enterococcus_faecalis | linezolid | 0.9 | 0.25 |
| enterococcus_faecalis | tedizolid | 0.9 | 0.005 |
| enterococcus_faecalis | daptomycin | 0.1 | 1 |
| enterococcus_faecalis | quinu_dalfo | 0.1 | 0.005 |
| enterococcus_faecalis | trim_sulf | 0.1 | 1 |
| enterococcus_faecalis | chloramphenicol | 0.8 | 1 |
| enterococcus_faecalis | nitrofurantoin | 0.9 | 1 |
| enterococcus_faecalis | fosfomycin | 0.1 | 1 |
| enterococcus_faecalis | retapamulin | 0.1 | 1 |
| enterococcus_faecalis | fusidic_a | 0.1 | 1 |
| enterococcus_faecalis | metronidazole | 0.1 | 1 |
| enterococcus_faecalis | fidaxomicin | 0.1 | 1 |
| enterococcus_faecalis | furazolidone | 0.1 | 1 |
| enterococcus_faecalis | rifampicin | 0.1 | 1 |
| enterococcus_faecalis | amoxicillin_clavulanate | 0.9 | 1 |
| enterococcus_faecalis | piperacillin_tazobactam | 0.75 | 1 |
| enterococcus_faecalis | ampicillin_sulbactam | 0.9 | 1 |
| enterococcus_faecalis | ticarcillin_clavulanate | 0.7 | 1 |
| enterococcus_faecalis | ceftazidime_avibactam | 0.1 | 0.005 |
| enterococcus_faecalis | meropenem_vaborbactam | 0.75 | 0.005 |
| enterococcus_faecalis | colistin | 0 | 0.005 |
| enterococcus_faecalis | flucloxacillin | 0.05 | 1 |
| enterococcus_faecalis | aztreonam_avibactam | 0.01 | 1 |
| enterococcus_faecalis | cefixime | 0.05 | 1 |
| enterococcus_faecium | sulfanilamide | 0.1 | 1 |
| enterococcus_faecium | penicillin_g | 0.1 | 1 |
| enterococcus_faecium | ampicillin | 0.3 | 1 |
| enterococcus_faecium | amoxicillin | 0.3 | 1 |
| enterococcus_faecium | piperacillin | 0.1 | 1 |
| enterococcus_faecium | ticarcillin | 0.1 | 1 |
| enterococcus_faecium | cephalexin | 0.1 | 1 |
| enterococcus_faecium | cefazolin | 0.1 | 1 |
| enterococcus_faecium | cefuroxime | 0.1 | 1 |
| enterococcus_faecium | ceftriaxone | 0.1 | 1 |
| enterococcus_faecium | ceftazidime | 0.1 | 1 |
| enterococcus_faecium | cefepime | 0.1 | 1 |
| enterococcus_faecium | ceftaroline | 0.1 | 1 |
| enterococcus_faecium | ceftolozane_tazobactam | 0.05 | 1 |
| enterococcus_faecium | cefiderocol | 0.05 | 1 |
| enterococcus_faecium | meropenem | 0.1 | 0.005 |
| enterococcus_faecium | imipenem_c | 0.1 | 0.005 |
| enterococcus_faecium | ertapenem | 0.1 | 0.005 |
| enterococcus_faecium | erythromycin | 0.7 | 0.01 |
| enterococcus_faecium | azithromycin | 0.7 | 0.01 |
| enterococcus_faecium | clarithromycin | 0.7 | 1 |
| enterococcus_faecium | clindamycin | 0.7 | 1 |
| enterococcus_faecium | gentamicin | 0.75 | 1 |
| enterococcus_faecium | tobramycin | 0.1 | 1 |
| enterococcus_faecium | amikacin | 0.1 | 1 |
| enterococcus_faecium | ciprofloxacin | 0.7 | 1 |
| enterococcus_faecium | levofloxacin | 0.7 | 1 |
| enterococcus_faecium | moxifloxacin | 0.7 | 1 |
| enterococcus_faecium | ofloxacin | 0.7 | 1 |
| enterococcus_faecium | tetracycline | 0.8 | 1 |
| enterococcus_faecium | doxycycline | 0.8 | 1 |
| enterococcus_faecium | minocycline | 0.85 | 1 |
| enterococcus_faecium | tigecycline | 0.1 | 1 |
| enterococcus_faecium | vancomycin | 0.9 | 3.5 |
| enterococcus_faecium | teicoplanin | 0.85 | 1 |
| enterococcus_faecium | dalbavancin | 0.85 | 0.005 |
| enterococcus_faecium | linezolid | 0.9 | 0.3 |
| enterococcus_faecium | tedizolid | 0.9 | 0.005 |
| enterococcus_faecium | daptomycin | 0.1 | 1 |
| enterococcus_faecium | quinu_dalfo | 0.7 | 0.005 |
| enterococcus_faecium | trim_sulf | 0.6 | 1 |
| enterococcus_faecium | chloramphenicol | 0.7 | 1 |
| enterococcus_faecium | nitrofurantoin | 0.7 | 1 |
| enterococcus_faecium | fosfomycin | 0.1 | 1 |
| enterococcus_faecium | retapamulin | 0.1 | 1 |
| enterococcus_faecium | fusidic_a | 0.1 | 1 |
| enterococcus_faecium | metronidazole | 0.1 | 1 |
| enterococcus_faecium | fidaxomicin | 0.1 | 1 |
| enterococcus_faecium | furazolidone | 0.1 | 1 |
| enterococcus_faecium | rifampicin | 0.1 | 1 |
| enterococcus_faecium | amoxicillin_clavulanate | 0.8 | 1 |
| enterococcus_faecium | piperacillin_tazobactam | 0.1 | 1 |
| enterococcus_faecium | ampicillin_sulbactam | 0.85 | 1 |
| enterococcus_faecium | ticarcillin_clavulanate | 0.1 | 1 |
| enterococcus_faecium | ceftazidime_avibactam | 0.1 | 0.005 |
| enterococcus_faecium | meropenem_vaborbactam | 0.1 | 0.005 |
| enterococcus_faecium | colistin | 0 | 0.005 |
| enterococcus_faecium | flucloxacillin | 0.05 | 1 |
| enterococcus_faecium | aztreonam_avibactam | 0.01 | 1 |
| enterococcus_faecium | cefixime | 0.05 | 1 |
| escherichia_coli | sulfanilamide | 0.5 | 1 |
| escherichia_coli | penicillin_g | 0.1 | 1 |
| escherichia_coli | ampicillin | 0.8 | 1 |
| escherichia_coli | amoxicillin | 0.8 | 1 |
| escherichia_coli | piperacillin | 0.85 | 1 |
| escherichia_coli | ticarcillin | 0.8 | 1 |
| escherichia_coli | cephalexin | 0.7 | 1 |
| escherichia_coli | cefazolin | 0.75 | 1 |
| escherichia_coli | cefuroxime | 0.8 | 1 |
| escherichia_coli | ceftriaxone | 0.9 | 1 |
| escherichia_coli | ceftazidime | 0.9 | 1 |
| escherichia_coli | cefepime | 0.9 | 1 |
| escherichia_coli | ceftaroline | 0.7 | 1 |
| escherichia_coli | ceftolozane_tazobactam | 0.8 | 1 |
| escherichia_coli | cefiderocol | 0.8 | 1 |
| escherichia_coli | meropenem | 0.95 | 0.05 |
| escherichia_coli | imipenem_c | 0.95 | 0.005 |
| escherichia_coli | ertapenem | 0.95 | 0.05 |
| escherichia_coli | aztreonam | 0.9 | 1 |
| escherichia_coli | gentamicin | 0.9 | 1 |
| escherichia_coli | tobramycin | 0.85 | 1 |
| escherichia_coli | amikacin | 0.9 | 1 |
| escherichia_coli | ciprofloxacin | 0.95 | 1 |
| escherichia_coli | levofloxacin | 0.9 | 1 |
| escherichia_coli | moxifloxacin | 0.8 | 1 |
| escherichia_coli | ofloxacin | 0.9 | 1 |
| escherichia_coli | tetracycline | 0.8 | 1 |
| escherichia_coli | doxycycline | 0.8 | 1 |
| escherichia_coli | minocycline | 0.85 | 1 |
| escherichia_coli | tigecycline | 0.1 | 1 |
| escherichia_coli | dalbavancin | 0 | 0.005 |
| escherichia_coli | linezolid | 0 | 0.005 |
| escherichia_coli | tedizolid | 0 | 0.005 |
| escherichia_coli | daptomycin | 0.1 | 1 |
| escherichia_coli | quinu_dalfo | 0 | 0.005 |
| escherichia_coli | trim_sulf | 0.9 | 1 |
| escherichia_coli | chloramphenicol | 0.85 | 1 |
| escherichia_coli | nitrofurantoin | 0.95 | 1 |
| escherichia_coli | fosfomycin | 0.1 | 1 |
| escherichia_coli | fidaxomicin | 0.1 | 1 |
| escherichia_coli | furazolidone | 0.1 | 1 |
| escherichia_coli | rifampicin | 0.7 | 1 |
| escherichia_coli | amoxicillin_clavulanate | 0.9 | 1 |
| escherichia_coli | piperacillin_tazobactam | 0.97 | 1 |
| escherichia_coli | ampicillin_sulbactam | 0.9 | 1 |
| escherichia_coli | ticarcillin_clavulanate | 0.9 | 1 |
| escherichia_coli | ceftazidime_avibactam | 0.95 | 0.005 |
| escherichia_coli | meropenem_vaborbactam | 0.95 | 0.005 |
| escherichia_coli | colistin | 0.7 | 0.005 |
| escherichia_coli | flucloxacillin | 0.01 | 1 |
| escherichia_coli | aztreonam_avibactam | 1 | 1 |
| escherichia_coli | cefixime | 0.8 | 1 |
| klebsiella_pneumoniae | sulfanilamide | 0.5 | 1 |
| klebsiella_pneumoniae | penicillin_g | 0.1 | 1 |
| klebsiella_pneumoniae | ampicillin | 0.1 | 1 |
| klebsiella_pneumoniae | amoxicillin | 0.1 | 1 |
| klebsiella_pneumoniae | piperacillin | 0.8 | 1 |
| klebsiella_pneumoniae | ticarcillin | 0.75 | 1 |
| klebsiella_pneumoniae | cephalexin | 0.5 | 1 |
| klebsiella_pneumoniae | cefazolin | 0.5 | 1 |
| klebsiella_pneumoniae | cefuroxime | 0.7 | 1 |
| klebsiella_pneumoniae | ceftriaxone | 0.9 | 1 |
| klebsiella_pneumoniae | ceftazidime | 0.85 | 1 |
| klebsiella_pneumoniae | cefepime | 0.92 | 1 |
| klebsiella_pneumoniae | ceftaroline | 0.5 | 1 |
| klebsiella_pneumoniae | ceftolozane_tazobactam | 0.8 | 1 |
| klebsiella_pneumoniae | cefiderocol | 0.8 | 1 |
| klebsiella_pneumoniae | meropenem | 0.94 | 0.05 |
| klebsiella_pneumoniae | imipenem_c | 0.95 | 0.05 |
| klebsiella_pneumoniae | ertapenem | 0.94 | 0.05 |
| klebsiella_pneumoniae | aztreonam | 0.85 | 1 |
| klebsiella_pneumoniae | gentamicin | 0.9 | 1 |
| klebsiella_pneumoniae | tobramycin | 0.85 | 1 |
| klebsiella_pneumoniae | amikacin | 0.9 | 1 |
| klebsiella_pneumoniae | ciprofloxacin | 0.9 | 1 |
| klebsiella_pneumoniae | levofloxacin | 0.85 | 1 |
| klebsiella_pneumoniae | moxifloxacin | 0.7 | 1 |
| klebsiella_pneumoniae | ofloxacin | 0.8 | 1 |
| klebsiella_pneumoniae | tetracycline | 0.8 | 1 |
| klebsiella_pneumoniae | doxycycline | 0.8 | 1 |
| klebsiella_pneumoniae | minocycline | 0.85 | 1 |
| klebsiella_pneumoniae | tigecycline | 0.1 | 1 |
| klebsiella_pneumoniae | dalbavancin | 0 | 0.005 |
| klebsiella_pneumoniae | linezolid | 0 | 0.005 |
| klebsiella_pneumoniae | tedizolid | 0 | 0.005 |
| klebsiella_pneumoniae | daptomycin | 0.1 | 1 |
| klebsiella_pneumoniae | quinu_dalfo | 0 | 0.005 |
| klebsiella_pneumoniae | trim_sulf | 0.9 | 1 |
| klebsiella_pneumoniae | chloramphenicol | 0.85 | 1 |
| klebsiella_pneumoniae | nitrofurantoin | 0.8 | 1 |
| klebsiella_pneumoniae | fosfomycin | 0.1 | 1 |
| klebsiella_pneumoniae | fidaxomicin | 0.1 | 1 |
| klebsiella_pneumoniae | furazolidone | 0.1 | 1 |
| klebsiella_pneumoniae | rifampicin | 0.6 | 1 |
| klebsiella_pneumoniae | amoxicillin_clavulanate | 0.85 | 1 |
| klebsiella_pneumoniae | piperacillin_tazobactam | 0.92 | 1 |
| klebsiella_pneumoniae | ampicillin_sulbactam | 0.75 | 1 |
| klebsiella_pneumoniae | ticarcillin_clavulanate | 0.75 | 1 |
| klebsiella_pneumoniae | ceftazidime_avibactam | 0.95 | 0.005 |
| klebsiella_pneumoniae | meropenem_vaborbactam | 0.95 | 0.005 |
| klebsiella_pneumoniae | colistin | 0.7 | 0.005 |
| klebsiella_pneumoniae | flucloxacillin | 0.01 | 1 |
| klebsiella_pneumoniae | aztreonam_avibactam | 1 | 1 |
| klebsiella_pneumoniae | cefixime | 0.8 | 1 |
| morganella_spp. | sulfanilamide | 0.5 | 1 |
| morganella_spp. | penicillin_g | 0.1 | 1 |
| morganella_spp. | ampicillin | 0.5 | 1 |
| morganella_spp. | amoxicillin | 0.5 | 1 |
| morganella_spp. | piperacillin | 0.75 | 1 |
| morganella_spp. | ticarcillin | 0.7 | 1 |
| morganella_spp. | cephalexin | 0.5 | 1 |
| morganella_spp. | cefazolin | 0.5 | 1 |
| morganella_spp. | cefuroxime | 0.6 | 1 |
| morganella_spp. | ceftriaxone | 0.8 | 1 |
| morganella_spp. | ceftazidime | 0.8 | 1 |
| morganella_spp. | cefepime | 0.85 | 1 |
| morganella_spp. | ceftaroline | 0.4 | 1 |
| morganella_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| morganella_spp. | cefiderocol | 0.8 | 1 |
| morganella_spp. | meropenem | 0.95 | 0.005 |
| morganella_spp. | imipenem_c | 0.95 | 0.005 |
| morganella_spp. | ertapenem | 0.9 | 0.005 |
| morganella_spp. | aztreonam | 0.8 | 1 |
| morganella_spp. | erythromycin | 0.1 | 1 |
| morganella_spp. | azithromycin | 0.1 | 1 |
| morganella_spp. | clarithromycin | 0.1 | 1 |
| morganella_spp. | clindamycin | 0.1 | 1 |
| morganella_spp. | gentamicin | 0.85 | 1 |
| morganella_spp. | tobramycin | 0.8 | 1 |
| morganella_spp. | amikacin | 0.9 | 1 |
| morganella_spp. | ciprofloxacin | 0.9 | 1 |
| morganella_spp. | levofloxacin | 0.85 | 1 |
| morganella_spp. | moxifloxacin | 0.7 | 1 |
| morganella_spp. | ofloxacin | 0.8 | 1 |
| morganella_spp. | tetracycline | 0.8 | 1 |
| morganella_spp. | doxycycline | 0.85 | 1 |
| morganella_spp. | minocycline | 0.85 | 1 |
| morganella_spp. | tigecycline | 0.1 | 1 |
| morganella_spp. | vancomycin | 0.1 | 1 |
| morganella_spp. | teicoplanin | 0.1 | 1 |
| morganella_spp. | dalbavancin | 0.1 | 0.005 |
| morganella_spp. | linezolid | 0.1 | 0.005 |
| morganella_spp. | tedizolid | 0.1 | 0.005 |
| morganella_spp. | daptomycin | 0.1 | 1 |
| morganella_spp. | quinu_dalfo | 0.1 | 0.005 |
| morganella_spp. | trim_sulf | 0.8 | 1 |
| morganella_spp. | chloramphenicol | 0.85 | 1 |
| morganella_spp. | nitrofurantoin | 0.7 | 1 |
| morganella_spp. | fosfomycin | 0.1 | 1 |
| morganella_spp. | retapamulin | 0.05 | 1 |
| morganella_spp. | fusidic_a | 0.05 | 1 |
| morganella_spp. | metronidazole | 0.05 | 1 |
| morganella_spp. | fidaxomicin | 0.1 | 1 |
| morganella_spp. | furazolidone | 0.1 | 1 |
| morganella_spp. | rifampicin | 0.6 | 1 |
| morganella_spp. | amoxicillin_clavulanate | 0.7 | 1 |
| morganella_spp. | piperacillin_tazobactam | 0.85 | 1 |
| morganella_spp. | ampicillin_sulbactam | 0.7 | 1 |
| morganella_spp. | ticarcillin_clavulanate | 0.8 | 1 |
| morganella_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| morganella_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| morganella_spp. | colistin | 0.7 | 0.005 |
| morganella_spp. | flucloxacillin | 0.01 | 1 |
| morganella_spp. | aztreonam_avibactam | 1 | 1 |
| morganella_spp. | cefixime | 0.8 | 1 |
| proteus_spp. | sulfanilamide | 0.5 | 1 |
| proteus_spp. | penicillin_g | 0.1 | 1 |
| proteus_spp. | ampicillin | 0.8 | 1 |
| proteus_spp. | amoxicillin | 0.8 | 1 |
| proteus_spp. | piperacillin | 0.85 | 1 |
| proteus_spp. | ticarcillin | 0.8 | 1 |
| proteus_spp. | cephalexin | 0.7 | 1 |
| proteus_spp. | cefazolin | 0.75 | 1 |
| proteus_spp. | cefuroxime | 0.8 | 1 |
| proteus_spp. | ceftriaxone | 0.95 | 1 |
| proteus_spp. | ceftazidime | 0.9 | 1 |
| proteus_spp. | cefepime | 0.9 | 1 |
| proteus_spp. | ceftaroline | 0.7 | 1 |
| proteus_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| proteus_spp. | cefiderocol | 0.8 | 1 |
| proteus_spp. | meropenem | 0.95 | 0.005 |
| proteus_spp. | imipenem_c | 0.95 | 0.005 |
| proteus_spp. | ertapenem | 0.95 | 0.005 |
| proteus_spp. | aztreonam | 0.9 | 1 |
| proteus_spp. | gentamicin | 0.8 | 1 |
| proteus_spp. | tobramycin | 0.75 | 1 |
| proteus_spp. | amikacin | 0.85 | 1 |
| proteus_spp. | ciprofloxacin | 0.9 | 1 |
| proteus_spp. | levofloxacin | 0.85 | 1 |
| proteus_spp. | moxifloxacin | 0.7 | 1 |
| proteus_spp. | ofloxacin | 0.8 | 1 |
| proteus_spp. | tetracycline | 0.8 | 1 |
| proteus_spp. | doxycycline | 0.8 | 1 |
| proteus_spp. | minocycline | 0.85 | 1 |
| proteus_spp. | tigecycline | 0.1 | 1 |
| proteus_spp. | dalbavancin | 0 | 0.005 |
| proteus_spp. | linezolid | 0 | 0.005 |
| proteus_spp. | tedizolid | 0 | 0.005 |
| proteus_spp. | daptomycin | 0.1 | 1 |
| proteus_spp. | quinu_dalfo | 0 | 0.005 |
| proteus_spp. | trim_sulf | 0.9 | 1 |
| proteus_spp. | chloramphenicol | 0.85 | 1 |
| proteus_spp. | nitrofurantoin | 0.8 | 1 |
| proteus_spp. | fosfomycin | 0.1 | 1 |
| proteus_spp. | fidaxomicin | 0.1 | 1 |
| proteus_spp. | furazolidone | 0.1 | 1 |
| proteus_spp. | rifampicin | 0.7 | 1 |
| proteus_spp. | amoxicillin_clavulanate | 0.9 | 1 |
| proteus_spp. | piperacillin_tazobactam | 0.95 | 1 |
| proteus_spp. | ampicillin_sulbactam | 0.9 | 1 |
| proteus_spp. | ticarcillin_clavulanate | 0.9 | 1 |
| proteus_spp. | ceftazidime_avibactam | 0.95 | 0.005 |
| proteus_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| proteus_spp. | colistin | 0.7 | 0.005 |
| proteus_spp. | flucloxacillin | 0.01 | 1 |
| proteus_spp. | aztreonam_avibactam | 1 | 1 |
| proteus_spp. | cefixime | 0.8 | 1 |
| serratia_spp. | sulfanilamide | 0.5 | 1 |
| serratia_spp. | penicillin_g | 0.1 | 1 |
| serratia_spp. | ampicillin | 0.1 | 1 |
| serratia_spp. | amoxicillin | 0.1 | 1 |
| serratia_spp. | piperacillin | 0.75 | 1 |
| serratia_spp. | ticarcillin | 0.7 | 1 |
| serratia_spp. | cephalexin | 0.1 | 1 |
| serratia_spp. | cefazolin | 0.1 | 1 |
| serratia_spp. | cefuroxime | 0.6 | 1 |
| serratia_spp. | ceftriaxone | 0.8 | 1 |
| serratia_spp. | ceftazidime | 0.85 | 1 |
| serratia_spp. | cefepime | 0.85 | 1 |
| serratia_spp. | ceftaroline | 0.5 | 1 |
| serratia_spp. | ceftolozane_tazobactam | 0.8 | 1 |
| serratia_spp. | cefiderocol | 0.8 | 1 |
| serratia_spp. | meropenem | 0.95 | 0.005 |
| serratia_spp. | imipenem_c | 0.95 | 0.005 |
| serratia_spp. | ertapenem | 0.9 | 0.005 |
| serratia_spp. | aztreonam | 0.85 | 1 |
| serratia_spp. | erythromycin | 0.1 | 1 |
| serratia_spp. | azithromycin | 0.1 | 1 |
| serratia_spp. | clarithromycin | 0.1 | 1 |
| serratia_spp. | clindamycin | 0.1 | 1 |
| serratia_spp. | gentamicin | 0.85 | 1 |
| serratia_spp. | tobramycin | 0.8 | 1 |
| serratia_spp. | amikacin | 0.9 | 1 |
| serratia_spp. | ciprofloxacin | 0.85 | 1 |
| serratia_spp. | levofloxacin | 0.8 | 1 |
| serratia_spp. | moxifloxacin | 0.7 | 1 |
| serratia_spp. | ofloxacin | 0.75 | 1 |
| serratia_spp. | tetracycline | 0.8 | 1 |
| serratia_spp. | doxycycline | 0.8 | 1 |
| serratia_spp. | minocycline | 0.85 | 1 |
| serratia_spp. | tigecycline | 0.1 | 1 |
| serratia_spp. | vancomycin | 0.1 | 1 |
| serratia_spp. | teicoplanin | 0.1 | 1 |
| serratia_spp. | dalbavancin | 0.1 | 0.005 |
| serratia_spp. | linezolid | 0.1 | 0.005 |
| serratia_spp. | tedizolid | 0.1 | 0.005 |
| serratia_spp. | daptomycin | 0.1 | 1 |
| serratia_spp. | quinu_dalfo | 0.1 | 0.005 |
| serratia_spp. | trim_sulf | 0.85 | 1 |
| serratia_spp. | chloramphenicol | 0.8 | 1 |
| serratia_spp. | nitrofurantoin | 0.7 | 1 |
| serratia_spp. | fosfomycin | 0.1 | 1 |
| serratia_spp. | retapamulin | 0.05 | 1 |
| serratia_spp. | fusidic_a | 0.05 | 1 |
| serratia_spp. | metronidazole | 0.05 | 1 |
| serratia_spp. | fidaxomicin | 0.1 | 1 |
| serratia_spp. | furazolidone | 0.1 | 1 |
| serratia_spp. | rifampicin | 0.6 | 1 |
| serratia_spp. | amoxicillin_clavulanate | 0.7 | 1 |
| serratia_spp. | piperacillin_tazobactam | 0.85 | 1 |
| serratia_spp. | ampicillin_sulbactam | 0.7 | 1 |
| serratia_spp. | ticarcillin_clavulanate | 0.75 | 1 |
| serratia_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| serratia_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| serratia_spp. | colistin | 0.7 | 0.005 |
| serratia_spp. | flucloxacillin | 0.01 | 1 |
| serratia_spp. | aztreonam_avibactam | 1 | 1 |
| serratia_spp. | cefixime | 0.8 | 1 |
| p_stuartii | sulfanilamide | 0.1 | 1 |
| p_stuartii | penicillin_g | 0.05 | 1 |
| p_stuartii | ampicillin | 0.05 | 1 |
| p_stuartii | amoxicillin | 0.05 | 1 |
| p_stuartii | piperacillin | 0.35 | 1 |
| p_stuartii | ticarcillin | 0.3 | 1 |
| p_stuartii | cephalexin | 0.1 | 1 |
| p_stuartii | cefazolin | 0.1 | 1 |
| p_stuartii | cefuroxime | 0.2 | 1 |
| p_stuartii | ceftriaxone | 0.45 | 1 |
| p_stuartii | ceftazidime | 0.75 | 1 |
| p_stuartii | cefepime | 0.85 | 1 |
| p_stuartii | ceftaroline | 0.2 | 1 |
| p_stuartii | ceftolozane_tazobactam | 0.8 | 1 |
| p_stuartii | cefiderocol | 0.8 | 1 |
| p_stuartii | meropenem | 0.9 | 0.005 |
| p_stuartii | imipenem_c | 0.9 | 0.005 |
| p_stuartii | ertapenem | 0.85 | 0.005 |
| p_stuartii | aztreonam | 0.65 | 1 |
| p_stuartii | erythromycin | 0.05 | 1 |
| p_stuartii | azithromycin | 0.05 | 1 |
| p_stuartii | clarithromycin | 0.05 | 1 |
| p_stuartii | clindamycin | 0.05 | 1 |
| p_stuartii | gentamicin | 0.45 | 1 |
| p_stuartii | tobramycin | 0.5 | 1 |
| p_stuartii | amikacin | 0.75 | 1 |
| p_stuartii | ciprofloxacin | 0.55 | 1 |
| p_stuartii | levofloxacin | 0.6 | 1 |
| p_stuartii | moxifloxacin | 0.6 | 1 |
| p_stuartii | ofloxacin | 0.55 | 1 |
| p_stuartii | tetracycline | 0.2 | 1 |
| p_stuartii | doxycycline | 0.3 | 1 |
| p_stuartii | minocycline | 0.35 | 1 |
| p_stuartii | tigecycline | 0.1 | 1 |
| p_stuartii | vancomycin | 0.05 | 1 |
| p_stuartii | teicoplanin | 0.05 | 1 |
| p_stuartii | dalbavancin | 0.05 | 0.005 |
| p_stuartii | linezolid | 0.05 | 0.005 |
| p_stuartii | tedizolid | 0.05 | 0.005 |
| p_stuartii | daptomycin | 0.1 | 1 |
| p_stuartii | quinu_dalfo | 0.05 | 0.005 |
| p_stuartii | trim_sulf | 0.3 | 1 |
| p_stuartii | chloramphenicol | 0.25 | 1 |
| p_stuartii | nitrofurantoin | 0.05 | 1 |
| p_stuartii | fosfomycin | 0.1 | 1 |
| p_stuartii | retapamulin | 0.05 | 1 |
| p_stuartii | fusidic_a | 0.05 | 1 |
| p_stuartii | metronidazole | 0.05 | 1 |
| p_stuartii | fidaxomicin | 0.1 | 1 |
| p_stuartii | furazolidone | 0.05 | 1 |
| p_stuartii | rifampicin | 0.1 | 1 |
| p_stuartii | amoxicillin_clavulanate | 0.2 | 1 |
| p_stuartii | piperacillin_tazobactam | 0.75 | 1 |
| p_stuartii | ampicillin_sulbactam | 0.2 | 1 |
| p_stuartii | ticarcillin_clavulanate | 0.45 | 1 |
| p_stuartii | ceftazidime_avibactam | 0.9 | 0.005 |
| p_stuartii | meropenem_vaborbactam | 0.95 | 0.005 |
| p_stuartii | colistin | 0.05 | 0.005 |
| p_stuartii | flucloxacillin | 0.01 | 1 |
| p_stuartii | aztreonam_avibactam | 1 | 1 |
| p_stuartii | cefixime | 0.8 | 1 |
| pseudomonas_aeruginosa | sulfanilamide | 0.1 | 1 |
| pseudomonas_aeruginosa | penicillin_g | 0.05 | 0.01 |
| pseudomonas_aeruginosa | ampicillin | 0.05 | 0.01 |
| pseudomonas_aeruginosa | amoxicillin | 0.05 | 0.01 |
| pseudomonas_aeruginosa | piperacillin | 0.8 | 1 |
| pseudomonas_aeruginosa | ticarcillin | 0.7 | 1 |
| pseudomonas_aeruginosa | cephalexin | 0.05 | 1 |
| pseudomonas_aeruginosa | cefazolin | 0.05 | 1 |
| pseudomonas_aeruginosa | cefuroxime | 0.1 | 1 |
| pseudomonas_aeruginosa | ceftriaxone | 0.1 | 1 |
| pseudomonas_aeruginosa | ceftazidime | 0.85 | 4.5 |
| pseudomonas_aeruginosa | cefepime | 0.9 | 4.5 |
| pseudomonas_aeruginosa | ceftaroline | 0.1 | 1 |
| pseudomonas_aeruginosa | ceftolozane_tazobactam | 0.1 | 1 |
| pseudomonas_aeruginosa | cefiderocol | 0.1 | 1 |
| pseudomonas_aeruginosa | meropenem | 0.9 | 0.05 |
| pseudomonas_aeruginosa | imipenem_c | 0.85 | 0.05 |
| pseudomonas_aeruginosa | ertapenem | 0.1 | 0.005 |
| pseudomonas_aeruginosa | aztreonam | 0.8 | 1 |
| pseudomonas_aeruginosa | gentamicin | 0.85 | 1 |
| pseudomonas_aeruginosa | tobramycin | 0.9 | 1 |
| pseudomonas_aeruginosa | amikacin | 0.9 | 1 |
| pseudomonas_aeruginosa | ciprofloxacin | 0.9 | 1 |
| pseudomonas_aeruginosa | levofloxacin | 0.8 | 1 |
| pseudomonas_aeruginosa | moxifloxacin | 0.5 | 1 |
| pseudomonas_aeruginosa | ofloxacin | 0.7 | 1 |
| pseudomonas_aeruginosa | tetracycline | 0.1 | 1 |
| pseudomonas_aeruginosa | doxycycline | 0.1 | 1 |
| pseudomonas_aeruginosa | minocycline | 0.1 | 1 |
| pseudomonas_aeruginosa | tigecycline | 0.1 | 1 |
| pseudomonas_aeruginosa | dalbavancin | 0 | 0.005 |
| pseudomonas_aeruginosa | linezolid | 0 | 0.005 |
| pseudomonas_aeruginosa | tedizolid | 0 | 0.005 |
| pseudomonas_aeruginosa | daptomycin | 0.1 | 1 |
| pseudomonas_aeruginosa | quinu_dalfo | 0 | 0.005 |
| pseudomonas_aeruginosa | trim_sulf | 0.1 | 1 |
| pseudomonas_aeruginosa | chloramphenicol | 0.1 | 1 |
| pseudomonas_aeruginosa | nitrofurantoin | 0.05 | 1 |
| pseudomonas_aeruginosa | fosfomycin | 0.1 | 1 |
| pseudomonas_aeruginosa | fidaxomicin | 0.1 | 1 |
| pseudomonas_aeruginosa | furazolidone | 0.05 | 1 |
| pseudomonas_aeruginosa | rifampicin | 0.1 | 1 |
| pseudomonas_aeruginosa | amoxicillin_clavulanate | 0.05 | 1 |
| pseudomonas_aeruginosa | piperacillin_tazobactam | 0.9 | 5 |
| pseudomonas_aeruginosa | ampicillin_sulbactam | 0.05 | 1 |
| pseudomonas_aeruginosa | ticarcillin_clavulanate | 0.8 | 1 |
| pseudomonas_aeruginosa | ceftazidime_avibactam | 0.95 | 0.005 |
| pseudomonas_aeruginosa | meropenem_vaborbactam | 0.9 | 0.005 |
| pseudomonas_aeruginosa | colistin | 0.85 | 0.05 |
| pseudomonas_aeruginosa | flucloxacillin | 0.01 | 1 |
| pseudomonas_aeruginosa | aztreonam_avibactam | 0.9 | 1 |
| pseudomonas_aeruginosa | cefixime | 0.1 | 1 |
| stenotrophomonas_maltophilia | sulfanilamide | 0.6 | 1 |
| stenotrophomonas_maltophilia | penicillin_g | 0.05 | 1 |
| stenotrophomonas_maltophilia | ampicillin | 0.05 | 1 |
| stenotrophomonas_maltophilia | amoxicillin | 0.05 | 1 |
| stenotrophomonas_maltophilia | piperacillin | 0.2 | 1 |
| stenotrophomonas_maltophilia | ticarcillin | 0.25 | 1 |
| stenotrophomonas_maltophilia | cephalexin | 0.05 | 1 |
| stenotrophomonas_maltophilia | cefazolin | 0.05 | 1 |
| stenotrophomonas_maltophilia | cefuroxime | 0.05 | 1 |
| stenotrophomonas_maltophilia | ceftriaxone | 0.05 | 1 |
| stenotrophomonas_maltophilia | ceftazidime | 0.35 | 0.1 |
| stenotrophomonas_maltophilia | cefepime | 0.15 | 1 |
| stenotrophomonas_maltophilia | ceftaroline | 0.05 | 1 |
| stenotrophomonas_maltophilia | ceftolozane_tazobactam | 0.1 | 1 |
| stenotrophomonas_maltophilia | cefiderocol | 0.1 | 1 |
| stenotrophomonas_maltophilia | meropenem | 0.05 | 0.01 |
| stenotrophomonas_maltophilia | imipenem_c | 0.05 | 0.01 |
| stenotrophomonas_maltophilia | ertapenem | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | aztreonam | 0.1 | 1 |
| stenotrophomonas_maltophilia | erythromycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | azithromycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | clarithromycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | clindamycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | gentamicin | 0.1 | 1 |
| stenotrophomonas_maltophilia | tobramycin | 0.1 | 1 |
| stenotrophomonas_maltophilia | amikacin | 0.1 | 1 |
| stenotrophomonas_maltophilia | ciprofloxacin | 0.4 | 1 |
| stenotrophomonas_maltophilia | levofloxacin | 0.75 | 3.5 |
| stenotrophomonas_maltophilia | moxifloxacin | 0.8 | 1 |
| stenotrophomonas_maltophilia | ofloxacin | 0.6 | 1 |
| stenotrophomonas_maltophilia | tetracycline | 0.35 | 1 |
| stenotrophomonas_maltophilia | doxycycline | 0.6 | 4.5 |
| stenotrophomonas_maltophilia | minocycline | 0.85 | 6 |
| stenotrophomonas_maltophilia | tigecycline | 0.1 | 1 |
| stenotrophomonas_maltophilia | vancomycin | 0.05 | 1 |
| stenotrophomonas_maltophilia | teicoplanin | 0.05 | 1 |
| stenotrophomonas_maltophilia | dalbavancin | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | linezolid | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | tedizolid | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | daptomycin | 0.1 | 1 |
| stenotrophomonas_maltophilia | quinu_dalfo | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | trim_sulf | 0.95 | 7 |
| stenotrophomonas_maltophilia | chloramphenicol | 0.4 | 1 |
| stenotrophomonas_maltophilia | nitrofurantoin | 0.05 | 1 |
| stenotrophomonas_maltophilia | fosfomycin | 0.1 | 1 |
| stenotrophomonas_maltophilia | retapamulin | 0.05 | 1 |
| stenotrophomonas_maltophilia | fusidic_a | 0.05 | 1 |
| stenotrophomonas_maltophilia | metronidazole | 0.05 | 1 |
| stenotrophomonas_maltophilia | fidaxomicin | 0.1 | 1 |
| stenotrophomonas_maltophilia | furazolidone | 0.05 | 1 |
| stenotrophomonas_maltophilia | rifampicin | 0.2 | 1 |
| stenotrophomonas_maltophilia | amoxicillin_clavulanate | 0.05 | 1 |
| stenotrophomonas_maltophilia | piperacillin_tazobactam | 0.3 | 0.05 |
| stenotrophomonas_maltophilia | ampicillin_sulbactam | 0.05 | 1 |
| stenotrophomonas_maltophilia | ticarcillin_clavulanate | 0.7 | 1 |
| stenotrophomonas_maltophilia | ceftazidime_avibactam | 0.4 | 0.005 |
| stenotrophomonas_maltophilia | meropenem_vaborbactam | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | colistin | 0.05 | 0.005 |
| stenotrophomonas_maltophilia | flucloxacillin | 0.01 | 1 |
| stenotrophomonas_maltophilia | aztreonam_avibactam | 0.75 | 1 |
| stenotrophomonas_maltophilia | cefixime | 0.1 | 1 |
| staphylococcus_aureus | sulfanilamide | 0.1 | 1 |
| staphylococcus_aureus | penicillin_g | 0.95 | 1 |
| staphylococcus_aureus | ampicillin | 0.1 | 1 |
| staphylococcus_aureus | amoxicillin | 0.1 | 1 |
| staphylococcus_aureus | piperacillin | 0.7 | 1 |
| staphylococcus_aureus | ticarcillin | 0.6 | 1 |
| staphylococcus_aureus | cephalexin | 0.8 | 1 |
| staphylococcus_aureus | cefazolin | 0.85 | 1 |
| staphylococcus_aureus | cefuroxime | 0.7 | 1 |
| staphylococcus_aureus | ceftriaxone | 0.7 | 1 |
| staphylococcus_aureus | ceftazidime | 0.1 | 1 |
| staphylococcus_aureus | cefepime | 0.6 | 1 |
| staphylococcus_aureus | ceftaroline | 0.95 | 1 |
| staphylococcus_aureus | ceftolozane_tazobactam | 0.75 | 1 |
| staphylococcus_aureus | cefiderocol | 0.75 | 1 |
| staphylococcus_aureus | meropenem | 0.7 | 0.005 |
| staphylococcus_aureus | imipenem_c | 0.7 | 0.005 |
| staphylococcus_aureus | ertapenem | 0.7 | 0.005 |
| staphylococcus_aureus | erythromycin | 0.8 | 1 |
| staphylococcus_aureus | azithromycin | 0.8 | 1 |
| staphylococcus_aureus | clarithromycin | 0.8 | 1 |
| staphylococcus_aureus | clindamycin | 0.8 | 1 |
| staphylococcus_aureus | gentamicin | 0.7 | 1 |
| staphylococcus_aureus | tobramycin | 0.7 | 1 |
| staphylococcus_aureus | amikacin | 0.7 | 1 |
| staphylococcus_aureus | ciprofloxacin | 0.7 | 1 |
| staphylococcus_aureus | levofloxacin | 0.7 | 1 |
| staphylococcus_aureus | moxifloxacin | 0.8 | 1 |
| staphylococcus_aureus | ofloxacin | 0.7 | 1 |
| staphylococcus_aureus | tetracycline | 0.8 | 1 |
| staphylococcus_aureus | doxycycline | 0.85 | 1 |
| staphylococcus_aureus | minocycline | 0.85 | 1 |
| staphylococcus_aureus | tigecycline | 0.1 | 1 |
| staphylococcus_aureus | vancomycin | 0.95 | 5 |
| staphylococcus_aureus | teicoplanin | 0.9 | 4 |
| staphylococcus_aureus | dalbavancin | 0.9 | 0.005 |
| staphylococcus_aureus | linezolid | 0.9 | 4 |
| staphylococcus_aureus | tedizolid | 0.9 | 0.005 |
| staphylococcus_aureus | daptomycin | 0.1 | 1 |
| staphylococcus_aureus | quinu_dalfo | 0.85 | 0.005 |
| staphylococcus_aureus | trim_sulf | 0.7 | 1 |
| staphylococcus_aureus | chloramphenicol | 0.8 | 1 |
| staphylococcus_aureus | nitrofurantoin | 0.1 | 1 |
| staphylococcus_aureus | fosfomycin | 0.1 | 1 |
| staphylococcus_aureus | retapamulin | 0.9 | 1 |
| staphylococcus_aureus | fusidic_a | 0.85 | 1 |
| staphylococcus_aureus | metronidazole | 0.1 | 1 |
| staphylococcus_aureus | fidaxomicin | 0.1 | 1 |
| staphylococcus_aureus | furazolidone | 0.1 | 1 |
| staphylococcus_aureus | rifampicin | 0.8 | 1 |
| staphylococcus_aureus | amoxicillin_clavulanate | 0.85 | 1 |
| staphylococcus_aureus | piperacillin_tazobactam | 0.7 | 1 |
| staphylococcus_aureus | ampicillin_sulbactam | 0.8 | 1 |
| staphylococcus_aureus | ticarcillin_clavulanate | 0.6 | 1 |
| staphylococcus_aureus | ceftazidime_avibactam | 0.1 | 0.005 |
| staphylococcus_aureus | meropenem_vaborbactam | 0.7 | 0.005 |
| staphylococcus_aureus | colistin | 0 | 0.005 |
| staphylococcus_aureus | flucloxacillin | 0.95 | 1 |
| staphylococcus_aureus | aztreonam_avibactam | 0.01 | 1 |
| staphylococcus_aureus | cefixime | 0.75 | 1 |
| staphylococcus_epidermidis | sulfanilamide | 0.1 | 1 |
| staphylococcus_epidermidis | penicillin_g | 0.15 | 1 |
| staphylococcus_epidermidis | ampicillin | 0.15 | 1 |
| staphylococcus_epidermidis | amoxicillin | 0.15 | 1 |
| staphylococcus_epidermidis | piperacillin | 0.2 | 1 |
| staphylococcus_epidermidis | ticarcillin | 0.2 | 1 |
| staphylococcus_epidermidis | cephalexin | 0.2 | 1 |
| staphylococcus_epidermidis | cefazolin | 0.2 | 1 |
| staphylococcus_epidermidis | cefuroxime | 0.2 | 1 |
| staphylococcus_epidermidis | ceftriaxone | 0.25 | 1 |
| staphylococcus_epidermidis | ceftazidime | 0.1 | 1 |
| staphylococcus_epidermidis | cefepime | 0.15 | 1 |
| staphylococcus_epidermidis | ceftaroline | 0.75 | 1 |
| staphylococcus_epidermidis | ceftolozane_tazobactam | 0.75 | 1 |
| staphylococcus_epidermidis | cefiderocol | 0.75 | 1 |
| staphylococcus_epidermidis | meropenem | 0.4 | 0.005 |
| staphylococcus_epidermidis | imipenem_c | 0.5 | 0.005 |
| staphylococcus_epidermidis | ertapenem | 0.4 | 0.005 |
| staphylococcus_epidermidis | aztreonam | 0.05 | 1 |
| staphylococcus_epidermidis | erythromycin | 0.45 | 1 |
| staphylococcus_epidermidis | azithromycin | 0.5 | 1 |
| staphylococcus_epidermidis | clarithromycin | 0.5 | 1 |
| staphylococcus_epidermidis | clindamycin | 0.6 | 1 |
| staphylococcus_epidermidis | gentamicin | 0.6 | 1 |
| staphylococcus_epidermidis | tobramycin | 0.65 | 1 |
| staphylococcus_epidermidis | amikacin | 0.7 | 1 |
| staphylococcus_epidermidis | ciprofloxacin | 0.5 | 1 |
| staphylococcus_epidermidis | levofloxacin | 0.55 | 1 |
| staphylococcus_epidermidis | moxifloxacin | 0.6 | 1 |
| staphylococcus_epidermidis | ofloxacin | 0.5 | 1 |
| staphylococcus_epidermidis | tetracycline | 0.5 | 1 |
| staphylococcus_epidermidis | doxycycline | 0.75 | 1 |
| staphylococcus_epidermidis | minocycline | 0.8 | 1 |
| staphylococcus_epidermidis | tigecycline | 0.1 | 1 |
| staphylococcus_epidermidis | vancomycin | 0.95 | 6 |
| staphylococcus_epidermidis | teicoplanin | 0.95 | 5 |
| staphylococcus_epidermidis | dalbavancin | 0.95 | 0.005 |
| staphylococcus_epidermidis | linezolid | 0.95 | 5 |
| staphylococcus_epidermidis | tedizolid | 0.95 | 0.005 |
| staphylococcus_epidermidis | daptomycin | 0.1 | 1 |
| staphylococcus_epidermidis | quinu_dalfo | 0.9 | 4 |
| staphylococcus_epidermidis | trim_sulf | 0.75 | 2.5 |
| staphylococcus_epidermidis | chloramphenicol | 0.6 | 1 |
| staphylococcus_epidermidis | nitrofurantoin | 0.2 | 1 |
| staphylococcus_epidermidis | fosfomycin | 0.1 | 1 |
| staphylococcus_epidermidis | retapamulin | 0.8 | 1 |
| staphylococcus_epidermidis | fusidic_a | 0.85 | 1 |
| staphylococcus_epidermidis | metronidazole | 0.05 | 1 |
| staphylococcus_epidermidis | fidaxomicin | 0.1 | 1 |
| staphylococcus_epidermidis | furazolidone | 0.05 | 1 |
| staphylococcus_epidermidis | rifampicin | 0.9 | 1 |
| staphylococcus_epidermidis | amoxicillin_clavulanate | 0.2 | 1 |
| staphylococcus_epidermidis | piperacillin_tazobactam | 0.4 | 1 |
| staphylococcus_epidermidis | ampicillin_sulbactam | 0.2 | 1 |
| staphylococcus_epidermidis | ticarcillin_clavulanate | 0.25 | 1 |
| staphylococcus_epidermidis | ceftazidime_avibactam | 0.1 | 0.005 |
| staphylococcus_epidermidis | meropenem_vaborbactam | 0.4 | 0.005 |
| staphylococcus_epidermidis | colistin | 0.05 | 0.005 |
| staphylococcus_epidermidis | flucloxacillin | 0.85 | 1 |
| staphylococcus_epidermidis | aztreonam_avibactam | 0.01 | 1 |
| staphylococcus_epidermidis | cefixime | 0.75 | 1 |
| streptococcus_pneumoniae | sulfanilamide | 0.1 | 1 |
| streptococcus_pneumoniae | penicillin_g | 0.95 | 1 |
| streptococcus_pneumoniae | ampicillin | 0.95 | 1 |
| streptococcus_pneumoniae | amoxicillin | 0.95 | 1 |
| streptococcus_pneumoniae | piperacillin | 0.9 | 1 |
| streptococcus_pneumoniae | ticarcillin | 0.9 | 1 |
| streptococcus_pneumoniae | cephalexin | 0.85 | 1 |
| streptococcus_pneumoniae | cefazolin | 0.9 | 1 |
| streptococcus_pneumoniae | cefuroxime | 0.9 | 1 |
| streptococcus_pneumoniae | ceftriaxone | 0.95 | 1 |
| streptococcus_pneumoniae | ceftazidime | 0.7 | 1 |
| streptococcus_pneumoniae | cefepime | 0.8 | 1 |
| streptococcus_pneumoniae | ceftaroline | 0.95 | 1 |
| streptococcus_pneumoniae | ceftolozane_tazobactam | 0.75 | 1 |
| streptococcus_pneumoniae | cefiderocol | 0.75 | 1 |
| streptococcus_pneumoniae | meropenem | 0.95 | 0.005 |
| streptococcus_pneumoniae | imipenem_c | 0.95 | 0.005 |
| streptococcus_pneumoniae | ertapenem | 0.95 | 0.005 |
| streptococcus_pneumoniae | erythromycin | 0.8 | 1 |
| streptococcus_pneumoniae | azithromycin | 0.85 | 1 |
| streptococcus_pneumoniae | clarithromycin | 0.85 | 1 |
| streptococcus_pneumoniae | clindamycin | 0.8 | 1 |
| streptococcus_pneumoniae | gentamicin | 0.1 | 1 |
| streptococcus_pneumoniae | tobramycin | 0.1 | 1 |
| streptococcus_pneumoniae | amikacin | 0.1 | 1 |
| streptococcus_pneumoniae | ciprofloxacin | 0.9 | 1 |
| streptococcus_pneumoniae | levofloxacin | 0.95 | 1 |
| streptococcus_pneumoniae | moxifloxacin | 0.95 | 1 |
| streptococcus_pneumoniae | ofloxacin | 0.9 | 1 |
| streptococcus_pneumoniae | tetracycline | 0.8 | 1 |
| streptococcus_pneumoniae | doxycycline | 0.85 | 1 |
| streptococcus_pneumoniae | minocycline | 0.85 | 1 |
| streptococcus_pneumoniae | tigecycline | 0.1 | 1 |
| streptococcus_pneumoniae | vancomycin | 0.95 | 1 |
| streptococcus_pneumoniae | teicoplanin | 0.9 | 1 |
| streptococcus_pneumoniae | dalbavancin | 0.9 | 0.005 |
| streptococcus_pneumoniae | linezolid | 0.9 | 0.005 |
| streptococcus_pneumoniae | tedizolid | 0.9 | 0.005 |
| streptococcus_pneumoniae | daptomycin | 0.1 | 1 |
| streptococcus_pneumoniae | quinu_dalfo | 0.85 | 0.005 |
| streptococcus_pneumoniae | trim_sulf | 0.7 | 1 |
| streptococcus_pneumoniae | chloramphenicol | 0.8 | 1 |
| streptococcus_pneumoniae | nitrofurantoin | 0.1 | 1 |
| streptococcus_pneumoniae | fosfomycin | 0.1 | 1 |
| streptococcus_pneumoniae | retapamulin | 0.1 | 1 |
| streptococcus_pneumoniae | fusidic_a | 0.1 | 1 |
| streptococcus_pneumoniae | metronidazole | 0.1 | 1 |
| streptococcus_pneumoniae | fidaxomicin | 0.1 | 1 |
| streptococcus_pneumoniae | furazolidone | 0.1 | 1 |
| streptococcus_pneumoniae | rifampicin | 0.8 | 1 |
| streptococcus_pneumoniae | amoxicillin_clavulanate | 0.95 | 1 |
| streptococcus_pneumoniae | piperacillin_tazobactam | 0.9 | 1 |
| streptococcus_pneumoniae | ampicillin_sulbactam | 0.95 | 1 |
| streptococcus_pneumoniae | ticarcillin_clavulanate | 0.9 | 1 |
| streptococcus_pneumoniae | ceftazidime_avibactam | 0.95 | 0.005 |
| streptococcus_pneumoniae | meropenem_vaborbactam | 0.95 | 0.005 |
| streptococcus_pneumoniae | colistin | 0 | 0.005 |
| streptococcus_pneumoniae | flucloxacillin | 0.8 | 1 |
| streptococcus_pneumoniae | aztreonam_avibactam | 0.01 | 1 |
| streptococcus_pneumoniae | cefixime | 0.75 | 1 |
| salmonella_enterica_serovar_typhi | sulfanilamide | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | penicillin_g | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | ampicillin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | amoxicillin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | piperacillin | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | ticarcillin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | cephalexin | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | cefazolin | 0.75 | 1 |
| salmonella_enterica_serovar_typhi | cefuroxime | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | ceftriaxone | 0.95 | 1 |
| salmonella_enterica_serovar_typhi | ceftazidime | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | cefepime | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | ceftaroline | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | ceftolozane_tazobactam | 0.75 | 1 |
| salmonella_enterica_serovar_typhi | cefiderocol | 0.75 | 1 |
| salmonella_enterica_serovar_typhi | meropenem | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | imipenem_c | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | ertapenem | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | aztreonam | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | erythromycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | azithromycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | clarithromycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | gentamicin | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | tobramycin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | amikacin | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | ciprofloxacin | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | levofloxacin | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | moxifloxacin | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | ofloxacin | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | tetracycline | 0.8 | 1 |
| salmonella_enterica_serovar_typhi | doxycycline | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | minocycline | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | tigecycline | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | dalbavancin | 0 | 0.005 |
| salmonella_enterica_serovar_typhi | linezolid | 0 | 0.005 |
| salmonella_enterica_serovar_typhi | tedizolid | 0 | 0.005 |
| salmonella_enterica_serovar_typhi | daptomycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | quinu_dalfo | 0 | 0.005 |
| salmonella_enterica_serovar_typhi | trim_sulf | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | chloramphenicol | 0.85 | 1 |
| salmonella_enterica_serovar_typhi | nitrofurantoin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | fosfomycin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | fidaxomicin | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | furazolidone | 0.1 | 1 |
| salmonella_enterica_serovar_typhi | rifampicin | 0.7 | 1 |
| salmonella_enterica_serovar_typhi | amoxicillin_clavulanate | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | piperacillin_tazobactam | 0.95 | 1 |
| salmonella_enterica_serovar_typhi | ampicillin_sulbactam | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | ticarcillin_clavulanate | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | ceftazidime_avibactam | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | meropenem_vaborbactam | 0.95 | 0.005 |
| salmonella_enterica_serovar_typhi | colistin | 0.7 | 0.005 |
| salmonella_enterica_serovar_typhi | flucloxacillin | 0.01 | 1 |
| salmonella_enterica_serovar_typhi | aztreonam_avibactam | 0.9 | 1 |
| salmonella_enterica_serovar_typhi | cefixime | 0.75 | 1 |
| salmonella_enterica_serovar_paratyphi_a | sulfanilamide | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | penicillin_g | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ampicillin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | amoxicillin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | piperacillin | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ticarcillin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | cephalexin | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | cefazolin | 0.75 | 1 |
| salmonella_enterica_serovar_paratyphi_a | cefuroxime | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ceftriaxone | 0.95 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ceftazidime | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | cefepime | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ceftaroline | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ceftolozane_tazobactam | 0.75 | 1 |
| salmonella_enterica_serovar_paratyphi_a | cefiderocol | 0.75 | 1 |
| salmonella_enterica_serovar_paratyphi_a | meropenem | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | imipenem_c | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | ertapenem | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | aztreonam | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | erythromycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | azithromycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | clarithromycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | gentamicin | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | tobramycin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | amikacin | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ciprofloxacin | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | levofloxacin | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | moxifloxacin | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ofloxacin | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | tetracycline | 0.8 | 1 |
| salmonella_enterica_serovar_paratyphi_a | doxycycline | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | minocycline | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | tigecycline | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | dalbavancin | 0 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | linezolid | 0 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | tedizolid | 0 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | daptomycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | quinu_dalfo | 0 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | trim_sulf | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | chloramphenicol | 0.85 | 1 |
| salmonella_enterica_serovar_paratyphi_a | nitrofurantoin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | fosfomycin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | fidaxomicin | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | furazolidone | 0.1 | 1 |
| salmonella_enterica_serovar_paratyphi_a | rifampicin | 0.7 | 1 |
| salmonella_enterica_serovar_paratyphi_a | amoxicillin_clavulanate | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | piperacillin_tazobactam | 0.95 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ampicillin_sulbactam | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ticarcillin_clavulanate | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | ceftazidime_avibactam | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | meropenem_vaborbactam | 0.95 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | colistin | 0.7 | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | flucloxacillin | 0.01 | 1 |
| salmonella_enterica_serovar_paratyphi_a | aztreonam_avibactam | 0.9 | 1 |
| salmonella_enterica_serovar_paratyphi_a | cefixime | 0.75 | 1 |
| invasive_non-typhoidal_salmonella_spp. | sulfanilamide | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | penicillin_g | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ampicillin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | amoxicillin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | piperacillin | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ticarcillin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | cephalexin | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | cefazolin | 0.75 | 1 |
| invasive_non-typhoidal_salmonella_spp. | cefuroxime | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ceftriaxone | 0.95 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ceftazidime | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | cefepime | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ceftaroline | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ceftolozane_tazobactam | 0.75 | 1 |
| invasive_non-typhoidal_salmonella_spp. | cefiderocol | 0.75 | 1 |
| invasive_non-typhoidal_salmonella_spp. | meropenem | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | imipenem_c | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | ertapenem | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | aztreonam | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | erythromycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | azithromycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | clarithromycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | clindamycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | gentamicin | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | tobramycin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | amikacin | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ciprofloxacin | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | levofloxacin | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | moxifloxacin | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ofloxacin | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | tetracycline | 0.8 | 1 |
| invasive_non-typhoidal_salmonella_spp. | doxycycline | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | minocycline | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | tigecycline | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | vancomycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | teicoplanin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | dalbavancin | 0.1 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | linezolid | 0.1 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | tedizolid | 0.1 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | daptomycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | quinu_dalfo | 0.1 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | trim_sulf | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | chloramphenicol | 0.85 | 1 |
| invasive_non-typhoidal_salmonella_spp. | nitrofurantoin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | fosfomycin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | retapamulin | 0.05 | 1 |
| invasive_non-typhoidal_salmonella_spp. | fusidic_a | 0.05 | 1 |
| invasive_non-typhoidal_salmonella_spp. | fidaxomicin | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | furazolidone | 0.1 | 1 |
| invasive_non-typhoidal_salmonella_spp. | rifampicin | 0.7 | 1 |
| invasive_non-typhoidal_salmonella_spp. | amoxicillin_clavulanate | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | piperacillin_tazobactam | 0.95 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ampicillin_sulbactam | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ticarcillin_clavulanate | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | ceftazidime_avibactam | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | meropenem_vaborbactam | 0.95 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | colistin | 0.7 | 0.005 |
| invasive_non-typhoidal_salmonella_spp. | flucloxacillin | 0.01 | 1 |
| invasive_non-typhoidal_salmonella_spp. | aztreonam_avibactam | 0.9 | 1 |
| invasive_non-typhoidal_salmonella_spp. | cefixime | 0.75 | 1 |
| shigella_spp. | sulfanilamide | 0.5 | 1 |
| shigella_spp. | penicillin_g | 0.1 | 1 |
| shigella_spp. | ampicillin | 0.7 | 1 |
| shigella_spp. | amoxicillin | 0.7 | 1 |
| shigella_spp. | piperacillin | 0.75 | 1 |
| shigella_spp. | ticarcillin | 0.7 | 1 |
| shigella_spp. | cephalexin | 0.6 | 1 |
| shigella_spp. | cefazolin | 0.65 | 1 |
| shigella_spp. | cefuroxime | 0.7 | 1 |
| shigella_spp. | ceftriaxone | 0.9 | 1 |
| shigella_spp. | ceftazidime | 0.85 | 1 |
| shigella_spp. | cefepime | 0.85 | 1 |
| shigella_spp. | ceftaroline | 0.6 | 1 |
| shigella_spp. | ceftolozane_tazobactam | 0.75 | 1 |
| shigella_spp. | cefiderocol | 0.75 | 1 |
| shigella_spp. | meropenem | 0.9 | 0.005 |
| shigella_spp. | imipenem_c | 0.9 | 0.005 |
| shigella_spp. | ertapenem | 0.9 | 0.005 |
| shigella_spp. | aztreonam | 0.8 | 1 |
| shigella_spp. | erythromycin | 0.7 | 1 |
| shigella_spp. | azithromycin | 0.85 | 1 |
| shigella_spp. | clarithromycin | 0.75 | 1 |
| shigella_spp. | gentamicin | 0.8 | 1 |
| shigella_spp. | tobramycin | 0.75 | 1 |
| shigella_spp. | amikacin | 0.85 | 1 |
| shigella_spp. | ciprofloxacin | 0.95 | 1 |
| shigella_spp. | levofloxacin | 0.9 | 1 |
| shigella_spp. | moxifloxacin | 0.8 | 1 |
| shigella_spp. | ofloxacin | 0.9 | 1 |
| shigella_spp. | tetracycline | 0.8 | 1 |
| shigella_spp. | doxycycline | 0.85 | 1 |
| shigella_spp. | minocycline | 0.85 | 1 |
| shigella_spp. | tigecycline | 0.7 | 1 |
| shigella_spp. | dalbavancin | 0 | 0.005 |
| shigella_spp. | linezolid | 0 | 0.005 |
| shigella_spp. | tedizolid | 0 | 0.005 |
| shigella_spp. | daptomycin | 0.1 | 1 |
| shigella_spp. | quinu_dalfo | 0 | 0.005 |
| shigella_spp. | trim_sulf | 0.9 | 1 |
| shigella_spp. | chloramphenicol | 0.85 | 1 |
| shigella_spp. | nitrofurantoin | 0.1 | 1 |
| shigella_spp. | fosfomycin | 0.1 | 1 |
| shigella_spp. | fidaxomicin | 0.1 | 1 |
| shigella_spp. | furazolidone | 0.1 | 1 |
| shigella_spp. | rifampicin | 0.7 | 1 |
| shigella_spp. | amoxicillin_clavulanate | 0.8 | 1 |
| shigella_spp. | piperacillin_tazobactam | 0.85 | 1 |
| shigella_spp. | ampicillin_sulbactam | 0.8 | 1 |
| shigella_spp. | ticarcillin_clavulanate | 0.85 | 1 |
| shigella_spp. | ceftazidime_avibactam | 0.9 | 0.005 |
| shigella_spp. | meropenem_vaborbactam | 0.9 | 0.005 |
| shigella_spp. | colistin | 0.7 | 0.005 |
| shigella_spp. | flucloxacillin | 0.01 | 1 |
| shigella_spp. | aztreonam_avibactam | 0.9 | 1 |
| shigella_spp. | cefixime | 0.75 | 1 |
| neisseria_gonorrhoeae | sulfanilamide | 0.1 | 1 |
| neisseria_gonorrhoeae | penicillin_g | 0.9 | 4 |
| neisseria_gonorrhoeae | ampicillin | 0.85 | 1 |
| neisseria_gonorrhoeae | amoxicillin | 0.85 | 2.5 |
| neisseria_gonorrhoeae | piperacillin | 0.8 | 1 |
| neisseria_gonorrhoeae | ticarcillin | 0.8 | 1 |
| neisseria_gonorrhoeae | cephalexin | 0.7 | 1 |
| neisseria_gonorrhoeae | cefazolin | 0.75 | 1 |
| neisseria_gonorrhoeae | cefuroxime | 0.85 | 1 |
| neisseria_gonorrhoeae | ceftriaxone | 0.95 | 12 |
| neisseria_gonorrhoeae | ceftazidime | 0.9 | 1 |
| neisseria_gonorrhoeae | cefepime | 0.9 | 1 |
| neisseria_gonorrhoeae | ceftaroline | 0.8 | 1 |
| neisseria_gonorrhoeae | ceftolozane_tazobactam | 0.8 | 1 |
| neisseria_gonorrhoeae | cefiderocol | 0.8 | 1 |
| neisseria_gonorrhoeae | meropenem | 0.9 | 0.005 |
| neisseria_gonorrhoeae | imipenem_c | 0.9 | 0.005 |
| neisseria_gonorrhoeae | ertapenem | 0.9 | 0.005 |
| neisseria_gonorrhoeae | aztreonam | 0.9 | 1 |
| neisseria_gonorrhoeae | erythromycin | 0.7 | 1 |
| neisseria_gonorrhoeae | azithromycin | 0.7 | 5 |
| neisseria_gonorrhoeae | clarithromycin | 0.7 | 1 |
| neisseria_gonorrhoeae | gentamicin | 0.7 | 2 |
| neisseria_gonorrhoeae | tobramycin | 0.7 | 1 |
| neisseria_gonorrhoeae | amikacin | 0.7 | 1 |
| neisseria_gonorrhoeae | ciprofloxacin | 0.9 | 5 |
| neisseria_gonorrhoeae | levofloxacin | 0.85 | 1 |
| neisseria_gonorrhoeae | moxifloxacin | 0.8 | 1 |
| neisseria_gonorrhoeae | ofloxacin | 0.85 | 1 |
| neisseria_gonorrhoeae | tetracycline | 0.8 | 1 |
| neisseria_gonorrhoeae | doxycycline | 0.9 | 4 |
| neisseria_gonorrhoeae | minocycline | 0.85 | 1 |
| neisseria_gonorrhoeae | tigecycline | 0.1 | 1 |
| neisseria_gonorrhoeae | dalbavancin | 0 | 0.005 |
| neisseria_gonorrhoeae | linezolid | 0 | 0.005 |
| neisseria_gonorrhoeae | tedizolid | 0 | 0.005 |
| neisseria_gonorrhoeae | daptomycin | 0.1 | 1 |
| neisseria_gonorrhoeae | quinu_dalfo | 0 | 0.005 |
| neisseria_gonorrhoeae | trim_sulf | 0.7 | 1 |
| neisseria_gonorrhoeae | chloramphenicol | 0.8 | 1 |
| neisseria_gonorrhoeae | nitrofurantoin | 0.1 | 1 |
| neisseria_gonorrhoeae | fosfomycin | 0.1 | 1 |
| neisseria_gonorrhoeae | fidaxomicin | 0.1 | 1 |
| neisseria_gonorrhoeae | furazolidone | 0.1 | 1 |
| neisseria_gonorrhoeae | rifampicin | 0.7 | 1 |
| neisseria_gonorrhoeae | amoxicillin_clavulanate | 0.85 | 1 |
| neisseria_gonorrhoeae | piperacillin_tazobactam | 0.85 | 1 |
| neisseria_gonorrhoeae | ampicillin_sulbactam | 0.8 | 1 |
| neisseria_gonorrhoeae | ticarcillin_clavulanate | 0.8 | 1 |
| neisseria_gonorrhoeae | ceftazidime_avibactam | 0.9 | 0.005 |
| neisseria_gonorrhoeae | meropenem_vaborbactam | 0.9 | 0.005 |
| neisseria_gonorrhoeae | colistin | 0.05 | 0.005 |
| neisseria_gonorrhoeae | flucloxacillin | 0.01 | 1 |
| neisseria_gonorrhoeae | aztreonam_avibactam | 0.8 | 1 |
| neisseria_gonorrhoeae | cefixime | 0.55 | 6 |
| streptococcus_pyogenes | sulfanilamide | 0.1 | 1 |
| streptococcus_pyogenes | penicillin_g | 1 | 1 |
| streptococcus_pyogenes | ampicillin | 0.95 | 1 |
| streptococcus_pyogenes | amoxicillin | 0.95 | 1 |
| streptococcus_pyogenes | piperacillin | 0.9 | 1 |
| streptococcus_pyogenes | ticarcillin | 0.9 | 1 |
| streptococcus_pyogenes | cephalexin | 0.9 | 1 |
| streptococcus_pyogenes | cefazolin | 0.9 | 1 |
| streptococcus_pyogenes | cefuroxime | 0.95 | 1 |
| streptococcus_pyogenes | ceftriaxone | 0.95 | 1 |
| streptococcus_pyogenes | ceftazidime | 0.7 | 1 |
| streptococcus_pyogenes | cefepime | 0.8 | 1 |
| streptococcus_pyogenes | ceftaroline | 0.95 | 1 |
| streptococcus_pyogenes | ceftolozane_tazobactam | 0.75 | 1 |
| streptococcus_pyogenes | cefiderocol | 0.75 | 1 |
| streptococcus_pyogenes | meropenem | 0.95 | 0.005 |
| streptococcus_pyogenes | imipenem_c | 0.95 | 0.005 |
| streptococcus_pyogenes | ertapenem | 0.95 | 0.005 |
| streptococcus_pyogenes | erythromycin | 0.9 | 1 |
| streptococcus_pyogenes | azithromycin | 0.9 | 1 |
| streptococcus_pyogenes | clarithromycin | 0.9 | 1 |
| streptococcus_pyogenes | clindamycin | 0.85 | 1 |
| streptococcus_pyogenes | gentamicin | 0.1 | 1 |
| streptococcus_pyogenes | tobramycin | 0.1 | 1 |
| streptococcus_pyogenes | amikacin | 0.1 | 1 |
| streptococcus_pyogenes | ciprofloxacin | 0.8 | 1 |
| streptococcus_pyogenes | levofloxacin | 0.9 | 1 |
| streptococcus_pyogenes | moxifloxacin | 0.9 | 1 |
| streptococcus_pyogenes | ofloxacin | 0.85 | 1 |
| streptococcus_pyogenes | tetracycline | 0.8 | 1 |
| streptococcus_pyogenes | doxycycline | 0.85 | 1 |
| streptococcus_pyogenes | minocycline | 0.85 | 1 |
| streptococcus_pyogenes | tigecycline | 0.1 | 1 |
| streptococcus_pyogenes | vancomycin | 0.95 | 1 |
| streptococcus_pyogenes | teicoplanin | 0.9 | 1 |
| streptococcus_pyogenes | dalbavancin | 0.9 | 0.005 |
| streptococcus_pyogenes | linezolid | 0.9 | 0.005 |
| streptococcus_pyogenes | tedizolid | 0.9 | 0.005 |
| streptococcus_pyogenes | daptomycin | 0.1 | 1 |
| streptococcus_pyogenes | quinu_dalfo | 0.85 | 0.005 |
| streptococcus_pyogenes | trim_sulf | 0.7 | 1 |
| streptococcus_pyogenes | chloramphenicol | 0.8 | 1 |
| streptococcus_pyogenes | nitrofurantoin | 0.1 | 1 |
| streptococcus_pyogenes | fosfomycin | 0.1 | 1 |
| streptococcus_pyogenes | retapamulin | 0.1 | 1 |
| streptococcus_pyogenes | fusidic_a | 0.1 | 1 |
| streptococcus_pyogenes | metronidazole | 0.1 | 1 |
| streptococcus_pyogenes | fidaxomicin | 0.1 | 1 |
| streptococcus_pyogenes | furazolidone | 0.1 | 1 |
| streptococcus_pyogenes | rifampicin | 0.8 | 1 |
| streptococcus_pyogenes | amoxicillin_clavulanate | 0.95 | 1 |
| streptococcus_pyogenes | piperacillin_tazobactam | 0.9 | 1 |
| streptococcus_pyogenes | ampicillin_sulbactam | 0.95 | 1 |
| streptococcus_pyogenes | ticarcillin_clavulanate | 0.9 | 1 |
| streptococcus_pyogenes | ceftazidime_avibactam | 0.95 | 0.005 |
| streptococcus_pyogenes | meropenem_vaborbactam | 0.95 | 0.005 |
| streptococcus_pyogenes | colistin | 0 | 0.005 |
| streptococcus_pyogenes | flucloxacillin | 0.8 | 1 |
| streptococcus_pyogenes | aztreonam_avibactam | 0.01 | 1 |
| streptococcus_pyogenes | cefixime | 0.75 | 1 |
| streptococcus_agalactiae | sulfanilamide | 0.1 | 1 |
| streptococcus_agalactiae | penicillin_g | 0.95 | 1 |
| streptococcus_agalactiae | ampicillin | 0.95 | 1 |
| streptococcus_agalactiae | amoxicillin | 0.95 | 1 |
| streptococcus_agalactiae | piperacillin | 0.9 | 1 |
| streptococcus_agalactiae | ticarcillin | 0.9 | 1 |
| streptococcus_agalactiae | cephalexin | 0.9 | 1 |
| streptococcus_agalactiae | cefazolin | 0.9 | 1 |
| streptococcus_agalactiae | cefuroxime | 0.95 | 1 |
| streptococcus_agalactiae | ceftriaxone | 0.95 | 1 |
| streptococcus_agalactiae | ceftazidime | 0.7 | 1 |
| streptococcus_agalactiae | cefepime | 0.8 | 1 |
| streptococcus_agalactiae | ceftaroline | 0.95 | 1 |
| streptococcus_agalactiae | ceftolozane_tazobactam | 0.75 | 1 |
| streptococcus_agalactiae | cefiderocol | 0.75 | 1 |
| streptococcus_agalactiae | meropenem | 0.95 | 0.005 |
| streptococcus_agalactiae | imipenem_c | 0.95 | 0.005 |
| streptococcus_agalactiae | ertapenem | 0.95 | 0.005 |
| streptococcus_agalactiae | erythromycin | 0.8 | 1 |
| streptococcus_agalactiae | azithromycin | 0.85 | 1 |
| streptococcus_agalactiae | clarithromycin | 0.85 | 1 |
| streptococcus_agalactiae | clindamycin | 0.8 | 1 |
| streptococcus_agalactiae | gentamicin | 0.1 | 1 |
| streptococcus_agalactiae | tobramycin | 0.1 | 1 |
| streptococcus_agalactiae | amikacin | 0.1 | 1 |
| streptococcus_agalactiae | ciprofloxacin | 0.8 | 1 |
| streptococcus_agalactiae | levofloxacin | 0.9 | 1 |
| streptococcus_agalactiae | moxifloxacin | 0.9 | 1 |
| streptococcus_agalactiae | ofloxacin | 0.85 | 1 |
| streptococcus_agalactiae | tetracycline | 0.8 | 1 |
| streptococcus_agalactiae | doxycycline | 0.85 | 1 |
| streptococcus_agalactiae | minocycline | 0.85 | 1 |
| streptococcus_agalactiae | tigecycline | 0.1 | 1 |
| streptococcus_agalactiae | vancomycin | 0.95 | 1 |
| streptococcus_agalactiae | teicoplanin | 0.9 | 1 |
| streptococcus_agalactiae | dalbavancin | 0.9 | 0.005 |
| streptococcus_agalactiae | linezolid | 0.9 | 0.005 |
| streptococcus_agalactiae | tedizolid | 0.9 | 0.005 |
| streptococcus_agalactiae | daptomycin | 0.1 | 1 |
| streptococcus_agalactiae | quinu_dalfo | 0.85 | 0.005 |
| streptococcus_agalactiae | trim_sulf | 0.7 | 1 |
| streptococcus_agalactiae | chloramphenicol | 0.8 | 1 |
| streptococcus_agalactiae | nitrofurantoin | 0.1 | 1 |
| streptococcus_agalactiae | fosfomycin | 0.1 | 1 |
| streptococcus_agalactiae | retapamulin | 0.1 | 1 |
| streptococcus_agalactiae | fusidic_a | 0.1 | 1 |
| streptococcus_agalactiae | metronidazole | 0.1 | 1 |
| streptococcus_agalactiae | fidaxomicin | 0.1 | 1 |
| streptococcus_agalactiae | furazolidone | 0.1 | 1 |
| streptococcus_agalactiae | rifampicin | 0.8 | 1 |
| streptococcus_agalactiae | amoxicillin_clavulanate | 0.95 | 1 |
| streptococcus_agalactiae | piperacillin_tazobactam | 0.9 | 1 |
| streptococcus_agalactiae | ampicillin_sulbactam | 0.95 | 1 |
| streptococcus_agalactiae | ticarcillin_clavulanate | 0.9 | 1 |
| streptococcus_agalactiae | ceftazidime_avibactam | 0.95 | 0.005 |
| streptococcus_agalactiae | meropenem_vaborbactam | 0.95 | 0.005 |
| streptococcus_agalactiae | colistin | 0 | 0.005 |
| streptococcus_agalactiae | flucloxacillin | 0.8 | 1 |
| streptococcus_agalactiae | aztreonam_avibactam | 0.01 | 1 |
| streptococcus_agalactiae | cefixime | 0.75 | 1 |
| haemophilus_influenzae | sulfanilamide | 0.1 | 1 |
| haemophilus_influenzae | penicillin_g | 0.7 | 1 |
| haemophilus_influenzae | ampicillin | 0.8 | 1 |
| haemophilus_influenzae | amoxicillin | 0.9 | 1 |
| haemophilus_influenzae | piperacillin | 0.85 | 1 |
| haemophilus_influenzae | ticarcillin | 0.8 | 1 |
| haemophilus_influenzae | cephalexin | 0.7 | 1 |
| haemophilus_influenzae | cefazolin | 0.75 | 1 |
| haemophilus_influenzae | cefuroxime | 0.85 | 1 |
| haemophilus_influenzae | ceftriaxone | 0.95 | 1 |
| haemophilus_influenzae | ceftazidime | 0.9 | 1 |
| haemophilus_influenzae | cefepime | 0.9 | 1 |
| haemophilus_influenzae | ceftaroline | 0.8 | 1 |
| haemophilus_influenzae | ceftolozane_tazobactam | 0.8 | 1 |
| haemophilus_influenzae | cefiderocol | 0.8 | 1 |
| haemophilus_influenzae | meropenem | 0.95 | 0.005 |
| haemophilus_influenzae | imipenem_c | 0.95 | 0.005 |
| haemophilus_influenzae | ertapenem | 0.95 | 0.005 |
| haemophilus_influenzae | aztreonam | 0.9 | 1 |
| haemophilus_influenzae | erythromycin | 0.7 | 4 |
| haemophilus_influenzae | azithromycin | 0.9 | 4 |
| haemophilus_influenzae | clarithromycin | 0.85 | 4 |
| haemophilus_influenzae | gentamicin | 0.7 | 1 |
| haemophilus_influenzae | tobramycin | 0.7 | 1 |
| haemophilus_influenzae | amikacin | 0.7 | 1 |
| haemophilus_influenzae | ciprofloxacin | 0.9 | 1 |
| haemophilus_influenzae | levofloxacin | 0.85 | 1 |
| haemophilus_influenzae | moxifloxacin | 0.8 | 1 |
| haemophilus_influenzae | ofloxacin | 0.85 | 1 |
| haemophilus_influenzae | tetracycline | 0.85 | 1 |
| haemophilus_influenzae | doxycycline | 0.85 | 1 |
| haemophilus_influenzae | minocycline | 0.85 | 1 |
| haemophilus_influenzae | tigecycline | 0.1 | 1 |
| haemophilus_influenzae | dalbavancin | 0 | 0.005 |
| haemophilus_influenzae | linezolid | 0 | 0.005 |
| haemophilus_influenzae | tedizolid | 0 | 0.005 |
| haemophilus_influenzae | daptomycin | 0.1 | 1 |
| haemophilus_influenzae | quinu_dalfo | 0 | 0.005 |
| haemophilus_influenzae | trim_sulf | 0.85 | 1 |
| haemophilus_influenzae | chloramphenicol | 0.8 | 1 |
| haemophilus_influenzae | nitrofurantoin | 0.1 | 1 |
| haemophilus_influenzae | fosfomycin | 0.1 | 1 |
| haemophilus_influenzae | fidaxomicin | 0.1 | 1 |
| haemophilus_influenzae | furazolidone | 0.1 | 1 |
| haemophilus_influenzae | rifampicin | 0.7 | 1 |
| haemophilus_influenzae | amoxicillin_clavulanate | 0.9 | 1 |
| haemophilus_influenzae | piperacillin_tazobactam | 0.85 | 1 |
| haemophilus_influenzae | ampicillin_sulbactam | 0.9 | 1 |
| haemophilus_influenzae | ticarcillin_clavulanate | 0.8 | 1 |
| haemophilus_influenzae | ceftazidime_avibactam | 0.95 | 0.005 |
| haemophilus_influenzae | meropenem_vaborbactam | 0.95 | 0.005 |
| haemophilus_influenzae | colistin | 0.05 | 0.005 |
| haemophilus_influenzae | flucloxacillin | 0.01 | 1 |
| haemophilus_influenzae | aztreonam_avibactam | 0.8 | 1 |
| haemophilus_influenzae | cefixime | 0.8 | 1 |
| chlamydia_trachomatis | sulfanilamide | 0.1 | 1 |
| chlamydia_trachomatis | penicillin_g | 0.1 | 1 |
| chlamydia_trachomatis | ampicillin | 0.1 | 1 |
| chlamydia_trachomatis | amoxicillin | 0.1 | 1 |
| chlamydia_trachomatis | piperacillin | 0.1 | 1 |
| chlamydia_trachomatis | ticarcillin | 0.1 | 1 |
| chlamydia_trachomatis | cephalexin | 0.1 | 1 |
| chlamydia_trachomatis | cefazolin | 0.1 | 1 |
| chlamydia_trachomatis | cefuroxime | 0.1 | 1 |
| chlamydia_trachomatis | ceftriaxone | 0.1 | 1 |
| chlamydia_trachomatis | ceftazidime | 0.1 | 1 |
| chlamydia_trachomatis | cefepime | 0.1 | 1 |
| chlamydia_trachomatis | ceftaroline | 0.1 | 1 |
| chlamydia_trachomatis | ceftolozane_tazobactam | 0.01 | 1 |
| chlamydia_trachomatis | cefiderocol | 0.01 | 1 |
| chlamydia_trachomatis | meropenem | 0.1 | 0.005 |
| chlamydia_trachomatis | imipenem_c | 0.1 | 0.005 |
| chlamydia_trachomatis | ertapenem | 0.1 | 0.005 |
| chlamydia_trachomatis | aztreonam | 0.1 | 1 |
| chlamydia_trachomatis | erythromycin | 0.8 | 1 |
| chlamydia_trachomatis | azithromycin | 0.95 | 4 |
| chlamydia_trachomatis | clarithromycin | 0.9 | 1 |
| chlamydia_trachomatis | clindamycin | 0.7 | 1 |
| chlamydia_trachomatis | gentamicin | 0.1 | 1 |
| chlamydia_trachomatis | tobramycin | 0.1 | 1 |
| chlamydia_trachomatis | amikacin | 0.1 | 1 |
| chlamydia_trachomatis | ciprofloxacin | 0.8 | 1 |
| chlamydia_trachomatis | levofloxacin | 0.85 | 1 |
| chlamydia_trachomatis | moxifloxacin | 0.85 | 1 |
| chlamydia_trachomatis | ofloxacin | 0.8 | 1 |
| chlamydia_trachomatis | tetracycline | 0.95 | 4.5 |
| chlamydia_trachomatis | doxycycline | 0.95 | 5 |
| chlamydia_trachomatis | minocycline | 0.9 | 1 |
| chlamydia_trachomatis | tigecycline | 0.85 | 1 |
| chlamydia_trachomatis | vancomycin | 0.1 | 1 |
| chlamydia_trachomatis | teicoplanin | 0.1 | 1 |
| chlamydia_trachomatis | dalbavancin | 0.1 | 0.005 |
| chlamydia_trachomatis | linezolid | 0.1 | 0.005 |
| chlamydia_trachomatis | tedizolid | 0.1 | 0.005 |
| chlamydia_trachomatis | daptomycin | 0.1 | 1 |
| chlamydia_trachomatis | quinu_dalfo | 0.1 | 0.005 |
| chlamydia_trachomatis | trim_sulf | 0.1 | 1 |
| chlamydia_trachomatis | chloramphenicol | 0.8 | 1 |
| chlamydia_trachomatis | nitrofurantoin | 0.1 | 1 |
| chlamydia_trachomatis | fosfomycin | 0.1 | 1 |
| chlamydia_trachomatis | retapamulin | 0.1 | 1 |
| chlamydia_trachomatis | fusidic_a | 0.1 | 1 |
| chlamydia_trachomatis | metronidazole | 0.1 | 1 |
| chlamydia_trachomatis | fidaxomicin | 0.1 | 1 |
| chlamydia_trachomatis | furazolidone | 0.1 | 1 |
| chlamydia_trachomatis | rifampicin | 0.1 | 1 |
| chlamydia_trachomatis | amoxicillin_clavulanate | 0.1 | 1 |
| chlamydia_trachomatis | piperacillin_tazobactam | 0.1 | 1 |
| chlamydia_trachomatis | ampicillin_sulbactam | 0.1 | 1 |
| chlamydia_trachomatis | ticarcillin_clavulanate | 0.1 | 1 |
| chlamydia_trachomatis | ceftazidime_avibactam | 0.1 | 0.005 |
| chlamydia_trachomatis | meropenem_vaborbactam | 0.1 | 0.005 |
| chlamydia_trachomatis | colistin | 0.1 | 0.005 |
| chlamydia_trachomatis | flucloxacillin | 0.01 | 1 |
| chlamydia_trachomatis | aztreonam_avibactam | 0.01 | 1 |
| chlamydia_trachomatis | cefixime | 0.01 | 1 |
| mycoplasma_genitalium | sulfanilamide | 0.05 | 1 |
| mycoplasma_genitalium | penicillin_g | 0.05 | 1 |
| mycoplasma_genitalium | ampicillin | 0.05 | 1 |
| mycoplasma_genitalium | amoxicillin | 0.05 | 1 |
| mycoplasma_genitalium | piperacillin | 0.05 | 1 |
| mycoplasma_genitalium | ticarcillin | 0.05 | 1 |
| mycoplasma_genitalium | cephalexin | 0.05 | 1 |
| mycoplasma_genitalium | cefazolin | 0.05 | 1 |
| mycoplasma_genitalium | cefuroxime | 0.05 | 1 |
| mycoplasma_genitalium | ceftriaxone | 0.05 | 1 |
| mycoplasma_genitalium | ceftazidime | 0.05 | 1 |
| mycoplasma_genitalium | cefepime | 0.05 | 1 |
| mycoplasma_genitalium | ceftaroline | 0.05 | 1 |
| mycoplasma_genitalium | ceftolozane_tazobactam | 0.01 | 1 |
| mycoplasma_genitalium | cefiderocol | 0.01 | 1 |
| mycoplasma_genitalium | meropenem | 0.05 | 0.005 |
| mycoplasma_genitalium | imipenem_c | 0.05 | 0.005 |
| mycoplasma_genitalium | ertapenem | 0.05 | 0.005 |
| mycoplasma_genitalium | aztreonam | 0.05 | 1 |
| mycoplasma_genitalium | erythromycin | 0.8 | 1 |
| mycoplasma_genitalium | azithromycin | 0.9 | 8 |
| mycoplasma_genitalium | clarithromycin | 0.9 | 1 |
| mycoplasma_genitalium | clindamycin | 0.2 | 1 |
| mycoplasma_genitalium | gentamicin | 0.05 | 1 |
| mycoplasma_genitalium | tobramycin | 0.05 | 1 |
| mycoplasma_genitalium | amikacin | 0.05 | 1 |
| mycoplasma_genitalium | ciprofloxacin | 0.3 | 1 |
| mycoplasma_genitalium | levofloxacin | 0.5 | 2.5 |
| mycoplasma_genitalium | moxifloxacin | 0.85 | 4 |
| mycoplasma_genitalium | ofloxacin | 0.45 | 1 |
| mycoplasma_genitalium | tetracycline | 0.4 | 1 |
| mycoplasma_genitalium | doxycycline | 0.6 | 1.5 |
| mycoplasma_genitalium | minocycline | 0.7 | 1 |
| mycoplasma_genitalium | tigecycline | 0.85 | 1 |
| mycoplasma_genitalium | vancomycin | 0.05 | 1 |
| mycoplasma_genitalium | teicoplanin | 0.05 | 1 |
| mycoplasma_genitalium | dalbavancin | 0.05 | 0.005 |
| mycoplasma_genitalium | linezolid | 0.05 | 0.005 |
| mycoplasma_genitalium | tedizolid | 0.05 | 0.005 |
| mycoplasma_genitalium | daptomycin | 0.1 | 1 |
| mycoplasma_genitalium | quinu_dalfo | 0.05 | 0.005 |
| mycoplasma_genitalium | trim_sulf | 0.05 | 1 |
| mycoplasma_genitalium | chloramphenicol | 0.2 | 1 |
| mycoplasma_genitalium | nitrofurantoin | 0.05 | 1 |
| mycoplasma_genitalium | fosfomycin | 0.1 | 1 |
| mycoplasma_genitalium | retapamulin | 0.05 | 1 |
| mycoplasma_genitalium | fusidic_a | 0.05 | 1 |
| mycoplasma_genitalium | metronidazole | 0.05 | 1 |
| mycoplasma_genitalium | fidaxomicin | 0.1 | 1 |
| mycoplasma_genitalium | furazolidone | 0.05 | 1 |
| mycoplasma_genitalium | rifampicin | 0.1 | 1 |
| mycoplasma_genitalium | amoxicillin_clavulanate | 0.05 | 1 |
| mycoplasma_genitalium | piperacillin_tazobactam | 0.05 | 1 |
| mycoplasma_genitalium | ampicillin_sulbactam | 0.05 | 1 |
| mycoplasma_genitalium | ticarcillin_clavulanate | 0.05 | 1 |
| mycoplasma_genitalium | ceftazidime_avibactam | 0.05 | 0.005 |
| mycoplasma_genitalium | meropenem_vaborbactam | 0.05 | 0.005 |
| mycoplasma_genitalium | colistin | 0.05 | 0.005 |
| mycoplasma_genitalium | flucloxacillin | 0.01 | 1 |
| mycoplasma_genitalium | aztreonam_avibactam | 0.01 | 1 |
| mycoplasma_genitalium | cefixime | 0.01 | 1 |
| vibrio_cholerae | sulfanilamide | 0.5 | 1 |
| vibrio_cholerae | penicillin_g | 0.7 | 1 |
| vibrio_cholerae | ampicillin | 0.8 | 1 |
| vibrio_cholerae | amoxicillin | 0.8 | 1 |
| vibrio_cholerae | piperacillin | 0.85 | 1 |
| vibrio_cholerae | ticarcillin | 0.8 | 1 |
| vibrio_cholerae | cephalexin | 0.7 | 1 |
| vibrio_cholerae | cefazolin | 0.75 | 1 |
| vibrio_cholerae | cefuroxime | 0.8 | 1 |
| vibrio_cholerae | ceftriaxone | 0.9 | 1 |
| vibrio_cholerae | ceftazidime | 0.85 | 1 |
| vibrio_cholerae | cefepime | 0.85 | 1 |
| vibrio_cholerae | ceftaroline | 0.7 | 1 |
| vibrio_cholerae | ceftolozane_tazobactam | 0.75 | 1 |
| vibrio_cholerae | cefiderocol | 0.75 | 1 |
| vibrio_cholerae | meropenem | 0.9 | 0.005 |
| vibrio_cholerae | imipenem_c | 0.9 | 0.005 |
| vibrio_cholerae | ertapenem | 0.9 | 0.005 |
| vibrio_cholerae | aztreonam | 0.8 | 1 |
| vibrio_cholerae | erythromycin | 0.7 | 1 |
| vibrio_cholerae | azithromycin | 0.8 | 1 |
| vibrio_cholerae | clarithromycin | 0.75 | 1 |
| vibrio_cholerae | clindamycin | 0.1 | 1 |
| vibrio_cholerae | gentamicin | 0.85 | 1 |
| vibrio_cholerae | tobramycin | 0.8 | 1 |
| vibrio_cholerae | amikacin | 0.85 | 1 |
| vibrio_cholerae | ciprofloxacin | 0.9 | 1 |
| vibrio_cholerae | levofloxacin | 0.85 | 1 |
| vibrio_cholerae | moxifloxacin | 0.75 | 1 |
| vibrio_cholerae | ofloxacin | 0.85 | 1 |
| vibrio_cholerae | tetracycline | 0.95 | 1 |
| vibrio_cholerae | doxycycline | 0.95 | 1 |
| vibrio_cholerae | minocycline | 0.9 | 1 |
| vibrio_cholerae | tigecycline | 0.7 | 1 |
| vibrio_cholerae | vancomycin | 0.1 | 1 |
| vibrio_cholerae | teicoplanin | 0.1 | 1 |
| vibrio_cholerae | dalbavancin | 0.1 | 0.005 |
| vibrio_cholerae | linezolid | 0.1 | 0.005 |
| vibrio_cholerae | tedizolid | 0.1 | 0.005 |
| vibrio_cholerae | daptomycin | 0.1 | 1 |
| vibrio_cholerae | quinu_dalfo | 0.1 | 0.005 |
| vibrio_cholerae | trim_sulf | 0.8 | 1 |
| vibrio_cholerae | chloramphenicol | 0.8 | 1 |
| vibrio_cholerae | nitrofurantoin | 0.1 | 1 |
| vibrio_cholerae | fosfomycin | 0.1 | 1 |
| vibrio_cholerae | retapamulin | 0.05 | 1 |
| vibrio_cholerae | fusidic_a | 0.05 | 1 |
| vibrio_cholerae | metronidazole | 0.05 | 1 |
| vibrio_cholerae | fidaxomicin | 0.1 | 1 |
| vibrio_cholerae | furazolidone | 0.1 | 1 |
| vibrio_cholerae | rifampicin | 0.7 | 1 |
| vibrio_cholerae | amoxicillin_clavulanate | 0.85 | 1 |
| vibrio_cholerae | piperacillin_tazobactam | 0.9 | 1 |
| vibrio_cholerae | ampicillin_sulbactam | 0.85 | 1 |
| vibrio_cholerae | ticarcillin_clavulanate | 0.85 | 1 |
| vibrio_cholerae | ceftazidime_avibactam | 0.9 | 0.005 |
| vibrio_cholerae | meropenem_vaborbactam | 0.9 | 0.005 |
| vibrio_cholerae | colistin | 0.7 | 0.005 |
| vibrio_cholerae | flucloxacillin | 0.01 | 1 |
| vibrio_cholerae | aztreonam_avibactam | 0.9 | 1 |
| vibrio_cholerae | cefixime | 0.75 | 1 |
| neisseria_meningitidis | sulfanilamide | 0.1 | 1 |
| neisseria_meningitidis | penicillin_g | 0.95 | 25 |
| neisseria_meningitidis | ampicillin | 0.9 | 22 |
| neisseria_meningitidis | amoxicillin | 0.9 | 1 |
| neisseria_meningitidis | piperacillin | 0.85 | 1 |
| neisseria_meningitidis | ticarcillin | 0.8 | 1 |
| neisseria_meningitidis | cephalexin | 0.8 | 1 |
| neisseria_meningitidis | cefazolin | 0.85 | 1 |
| neisseria_meningitidis | cefuroxime | 0.9 | 1 |
| neisseria_meningitidis | ceftriaxone | 0.95 | 30 |
| neisseria_meningitidis | ceftazidime | 0.9 | 1 |
| neisseria_meningitidis | cefepime | 0.9 | 1 |
| neisseria_meningitidis | ceftaroline | 0.8 | 1 |
| neisseria_meningitidis | ceftolozane_tazobactam | 0.8 | 1 |
| neisseria_meningitidis | cefiderocol | 0.8 | 1 |
| neisseria_meningitidis | meropenem | 0.95 | 0.005 |
| neisseria_meningitidis | imipenem_c | 0.95 | 0.005 |
| neisseria_meningitidis | ertapenem | 0.95 | 0.005 |
| neisseria_meningitidis | aztreonam | 0.9 | 1 |
| neisseria_meningitidis | erythromycin | 0.7 | 1 |
| neisseria_meningitidis | azithromycin | 0.8 | 1 |
| neisseria_meningitidis | clarithromycin | 0.75 | 1 |
| neisseria_meningitidis | clindamycin | 0.1 | 1 |
| neisseria_meningitidis | gentamicin | 0.7 | 1 |
| neisseria_meningitidis | tobramycin | 0.7 | 1 |
| neisseria_meningitidis | amikacin | 0.7 | 1 |
| neisseria_meningitidis | ciprofloxacin | 0.9 | 15 |
| neisseria_meningitidis | levofloxacin | 0.85 | 1 |
| neisseria_meningitidis | moxifloxacin | 0.8 | 1 |
| neisseria_meningitidis | ofloxacin | 0.85 | 1 |
| neisseria_meningitidis | tetracycline | 0.8 | 1 |
| neisseria_meningitidis | doxycycline | 0.8 | 1 |
| neisseria_meningitidis | minocycline | 0.85 | 1 |
| neisseria_meningitidis | tigecycline | 0.1 | 1 |
| neisseria_meningitidis | vancomycin | 0.1 | 1 |
| neisseria_meningitidis | teicoplanin | 0.1 | 1 |
| neisseria_meningitidis | dalbavancin | 0.1 | 0.005 |
| neisseria_meningitidis | linezolid | 0.1 | 0.005 |
| neisseria_meningitidis | tedizolid | 0.1 | 0.005 |
| neisseria_meningitidis | daptomycin | 0.1 | 1 |
| neisseria_meningitidis | quinu_dalfo | 0.1 | 0.005 |
| neisseria_meningitidis | trim_sulf | 0.7 | 1 |
| neisseria_meningitidis | chloramphenicol | 0.85 | 18 |
| neisseria_meningitidis | nitrofurantoin | 0.1 | 1 |
| neisseria_meningitidis | fosfomycin | 0.1 | 1 |
| neisseria_meningitidis | retapamulin | 0.05 | 1 |
| neisseria_meningitidis | fusidic_a | 0.05 | 1 |
| neisseria_meningitidis | metronidazole | 0.05 | 1 |
| neisseria_meningitidis | fidaxomicin | 0.1 | 1 |
| neisseria_meningitidis | furazolidone | 0.1 | 1 |
| neisseria_meningitidis | rifampicin | 0.85 | 12 |
| neisseria_meningitidis | amoxicillin_clavulanate | 0.9 | 1 |
| neisseria_meningitidis | piperacillin_tazobactam | 0.85 | 1 |
| neisseria_meningitidis | ampicillin_sulbactam | 0.9 | 1 |
| neisseria_meningitidis | ticarcillin_clavulanate | 0.85 | 1 |
| neisseria_meningitidis | ceftazidime_avibactam | 0.95 | 0.005 |
| neisseria_meningitidis | meropenem_vaborbactam | 0.95 | 0.005 |
| neisseria_meningitidis | colistin | 0.05 | 0.005 |
| neisseria_meningitidis | flucloxacillin | 0.01 | 1 |
| neisseria_meningitidis | aztreonam_avibactam | 0.8 | 1 |
| neisseria_meningitidis | cefixime | 0.8 | 1 |
| listeria_monocytogenes | sulfanilamide | 0.1 | 1 |
| listeria_monocytogenes | penicillin_g | 0.7 | 1 |
| listeria_monocytogenes | ampicillin | 0.95 | 1 |
| listeria_monocytogenes | amoxicillin | 0.95 | 1 |
| listeria_monocytogenes | piperacillin | 0.7 | 1 |
| listeria_monocytogenes | ticarcillin | 0.6 | 1 |
| listeria_monocytogenes | cephalexin | 0.1 | 1 |
| listeria_monocytogenes | cefazolin | 0.1 | 1 |
| listeria_monocytogenes | cefuroxime | 0.1 | 1 |
| listeria_monocytogenes | ceftriaxone | 0.1 | 1 |
| listeria_monocytogenes | ceftazidime | 0.1 | 1 |
| listeria_monocytogenes | cefepime | 0.1 | 1 |
| listeria_monocytogenes | ceftaroline | 0.1 | 1 |
| listeria_monocytogenes | ceftolozane_tazobactam | 0.05 | 1 |
| listeria_monocytogenes | cefiderocol | 0.05 | 1 |
| listeria_monocytogenes | meropenem | 0.7 | 0.005 |
| listeria_monocytogenes | imipenem_c | 0.7 | 0.005 |
| listeria_monocytogenes | ertapenem | 0.7 | 0.005 |
| listeria_monocytogenes | aztreonam | 0.1 | 1 |
| listeria_monocytogenes | erythromycin | 0.8 | 1 |
| listeria_monocytogenes | azithromycin | 0.85 | 1 |
| listeria_monocytogenes | clarithromycin | 0.8 | 1 |
| listeria_monocytogenes | clindamycin | 0.1 | 1 |
| listeria_monocytogenes | gentamicin | 0.1 | 1 |
| listeria_monocytogenes | tobramycin | 0.1 | 1 |
| listeria_monocytogenes | amikacin | 0.1 | 1 |
| listeria_monocytogenes | ciprofloxacin | 0.8 | 1 |
| listeria_monocytogenes | levofloxacin | 0.85 | 1 |
| listeria_monocytogenes | moxifloxacin | 0.8 | 1 |
| listeria_monocytogenes | ofloxacin | 0.8 | 1 |
| listeria_monocytogenes | tetracycline | 0.8 | 1 |
| listeria_monocytogenes | doxycycline | 0.85 | 1 |
| listeria_monocytogenes | minocycline | 0.85 | 1 |
| listeria_monocytogenes | tigecycline | 0.1 | 1 |
| listeria_monocytogenes | vancomycin | 0.1 | 1 |
| listeria_monocytogenes | teicoplanin | 0.1 | 1 |
| listeria_monocytogenes | dalbavancin | 0.1 | 0.005 |
| listeria_monocytogenes | linezolid | 0.1 | 0.005 |
| listeria_monocytogenes | tedizolid | 0.1 | 0.005 |
| listeria_monocytogenes | daptomycin | 0.1 | 1 |
| listeria_monocytogenes | quinu_dalfo | 0.1 | 0.005 |
| listeria_monocytogenes | trim_sulf | 0.9 | 1 |
| listeria_monocytogenes | chloramphenicol | 0.85 | 1 |
| listeria_monocytogenes | nitrofurantoin | 0.1 | 1 |
| listeria_monocytogenes | fosfomycin | 0.1 | 1 |
| listeria_monocytogenes | retapamulin | 0.1 | 1 |
| listeria_monocytogenes | fusidic_a | 0.1 | 1 |
| listeria_monocytogenes | metronidazole | 0.1 | 1 |
| listeria_monocytogenes | fidaxomicin | 0.1 | 1 |
| listeria_monocytogenes | furazolidone | 0.1 | 1 |
| listeria_monocytogenes | rifampicin | 0.8 | 1 |
| listeria_monocytogenes | amoxicillin_clavulanate | 0.7 | 1 |
| listeria_monocytogenes | piperacillin_tazobactam | 0.95 | 1 |
| listeria_monocytogenes | ampicillin_sulbactam | 0.6 | 1 |
| listeria_monocytogenes | ticarcillin_clavulanate | 0.1 | 1 |
| listeria_monocytogenes | ceftazidime_avibactam | 0.7 | 0.005 |
| listeria_monocytogenes | meropenem_vaborbactam | 0.05 | 0.005 |
| listeria_monocytogenes | colistin | 0.1 | 0.005 |
| listeria_monocytogenes | flucloxacillin | 0.05 | 1 |
| listeria_monocytogenes | aztreonam_avibactam | 0.01 | 1 |
| listeria_monocytogenes | cefixime | 0.05 | 1 |
| clostridioides_difficile | sulfanilamide | 0.1 | 1 |
| clostridioides_difficile | penicillin_g | 0.1 | 1 |
| clostridioides_difficile | ampicillin | 0.1 | 1 |
| clostridioides_difficile | amoxicillin | 0.1 | 1 |
| clostridioides_difficile | piperacillin | 0.1 | 1 |
| clostridioides_difficile | ticarcillin | 0.1 | 1 |
| clostridioides_difficile | cephalexin | 0.1 | 1 |
| clostridioides_difficile | cefazolin | 0.1 | 1 |
| clostridioides_difficile | cefuroxime | 0.1 | 1 |
| clostridioides_difficile | ceftriaxone | 0.1 | 1 |
| clostridioides_difficile | ceftazidime | 0.1 | 1 |
| clostridioides_difficile | cefepime | 0.1 | 1 |
| clostridioides_difficile | ceftaroline | 0.1 | 1 |
| clostridioides_difficile | ceftolozane_tazobactam | 0.05 | 1 |
| clostridioides_difficile | cefiderocol | 0.05 | 1 |
| clostridioides_difficile | meropenem | 0.1 | 0.005 |
| clostridioides_difficile | imipenem_c | 0.1 | 0.005 |
| clostridioides_difficile | ertapenem | 0.1 | 0.005 |
| clostridioides_difficile | aztreonam | 0.1 | 1 |
| clostridioides_difficile | erythromycin | 0.7 | 1 |
| clostridioides_difficile | azithromycin | 0.75 | 1 |
| clostridioides_difficile | clarithromycin | 0.7 | 1 |
| clostridioides_difficile | clindamycin | 0.1 | 1 |
| clostridioides_difficile | gentamicin | 0.1 | 1 |
| clostridioides_difficile | tobramycin | 0.1 | 1 |
| clostridioides_difficile | amikacin | 0.1 | 1 |
| clostridioides_difficile | ciprofloxacin | 0.1 | 1 |
| clostridioides_difficile | levofloxacin | 0.1 | 1 |
| clostridioides_difficile | moxifloxacin | 0.1 | 1 |
| clostridioides_difficile | ofloxacin | 0.1 | 1 |
| clostridioides_difficile | tetracycline | 0.7 | 1 |
| clostridioides_difficile | doxycycline | 0.7 | 1 |
| clostridioides_difficile | minocycline | 0.7 | 1 |
| clostridioides_difficile | tigecycline | 0.1 | 1 |
| clostridioides_difficile | vancomycin | 0.95 | 5 |
| clostridioides_difficile | teicoplanin | 0.9 | 1 |
| clostridioides_difficile | dalbavancin | 0.9 | 0.005 |
| clostridioides_difficile | linezolid | 0.85 | 0.005 |
| clostridioides_difficile | tedizolid | 0.85 | 0.005 |
| clostridioides_difficile | daptomycin | 0.1 | 1 |
| clostridioides_difficile | quinu_dalfo | 0.1 | 0.005 |
| clostridioides_difficile | trim_sulf | 0.1 | 1 |
| clostridioides_difficile | chloramphenicol | 0.1 | 1 |
| clostridioides_difficile | nitrofurantoin | 0.1 | 1 |
| clostridioides_difficile | fosfomycin | 0.1 | 1 |
| clostridioides_difficile | retapamulin | 0.1 | 1 |
| clostridioides_difficile | fusidic_a | 0.1 | 1 |
| clostridioides_difficile | metronidazole | 0.9 | 6 |
| clostridioides_difficile | fidaxomicin | 0.1 | 1 |
| clostridioides_difficile | furazolidone | 0.1 | 1 |
| clostridioides_difficile | rifampicin | 0.1 | 1 |
| clostridioides_difficile | amoxicillin_clavulanate | 0.1 | 1 |
| clostridioides_difficile | piperacillin_tazobactam | 0.1 | 1 |
| clostridioides_difficile | ampicillin_sulbactam | 0.1 | 1 |
| clostridioides_difficile | ticarcillin_clavulanate | 0.1 | 1 |
| clostridioides_difficile | ceftazidime_avibactam | 0.1 | 0.005 |
| clostridioides_difficile | meropenem_vaborbactam | 0.1 | 0.005 |
| clostridioides_difficile | colistin | 0.05 | 0.005 |
| clostridioides_difficile | flucloxacillin | 0.01 | 1 |
| clostridioides_difficile | aztreonam_avibactam | 0.01 | 1 |
| clostridioides_difficile | cefixime | 0.05 | 1 |
| bacteroides_fragilis | sulfanilamide | 0.05 | 1 |
| bacteroides_fragilis | penicillin_g | 0.1 | 1 |
| bacteroides_fragilis | ampicillin | 0.2 | 1 |
| bacteroides_fragilis | amoxicillin | 0.25 | 1 |
| bacteroides_fragilis | piperacillin | 0.5 | 1 |
| bacteroides_fragilis | ticarcillin | 0.4 | 1 |
| bacteroides_fragilis | cephalexin | 0.05 | 1 |
| bacteroides_fragilis | cefazolin | 0.05 | 1 |
| bacteroides_fragilis | cefuroxime | 0.2 | 1 |
| bacteroides_fragilis | ceftriaxone | 0.2 | 1 |
| bacteroides_fragilis | ceftazidime | 0.25 | 1 |
| bacteroides_fragilis | cefepime | 0.25 | 1 |
| bacteroides_fragilis | ceftaroline | 0.2 | 1 |
| bacteroides_fragilis | ceftolozane_tazobactam | 0.45 | 1 |
| bacteroides_fragilis | cefiderocol | 0.45 | 1 |
| bacteroides_fragilis | meropenem | 0.95 | 0.005 |
| bacteroides_fragilis | imipenem_c | 0.95 | 0.005 |
| bacteroides_fragilis | ertapenem | 0.95 | 0.005 |
| bacteroides_fragilis | aztreonam | 0.05 | 1 |
| bacteroides_fragilis | erythromycin | 0.05 | 1 |
| bacteroides_fragilis | azithromycin | 0.05 | 1 |
| bacteroides_fragilis | clarithromycin | 0.05 | 1 |
| bacteroides_fragilis | clindamycin | 0.6 | 1 |
| bacteroides_fragilis | gentamicin | 0.05 | 1 |
| bacteroides_fragilis | tobramycin | 0.05 | 1 |
| bacteroides_fragilis | amikacin | 0.05 | 1 |
| bacteroides_fragilis | ciprofloxacin | 0.25 | 1 |
| bacteroides_fragilis | levofloxacin | 0.35 | 1 |
| bacteroides_fragilis | moxifloxacin | 0.5 | 1 |
| bacteroides_fragilis | ofloxacin | 0.25 | 1 |
| bacteroides_fragilis | tetracycline | 0.3 | 1 |
| bacteroides_fragilis | doxycycline | 0.5 | 1 |
| bacteroides_fragilis | minocycline | 0.5 | 1 |
| bacteroides_fragilis | tigecycline | 0.1 | 1 |
| bacteroides_fragilis | vancomycin | 0.05 | 1 |
| bacteroides_fragilis | teicoplanin | 0.05 | 1 |
| bacteroides_fragilis | dalbavancin | 0.05 | 0.005 |
| bacteroides_fragilis | linezolid | 0.05 | 0.005 |
| bacteroides_fragilis | tedizolid | 0.05 | 0.005 |
| bacteroides_fragilis | daptomycin | 0.1 | 1 |
| bacteroides_fragilis | quinu_dalfo | 0.05 | 0.005 |
| bacteroides_fragilis | trim_sulf | 0.3 | 1 |
| bacteroides_fragilis | chloramphenicol | 0.7 | 1 |
| bacteroides_fragilis | nitrofurantoin | 0.05 | 1 |
| bacteroides_fragilis | fosfomycin | 0.1 | 1 |
| bacteroides_fragilis | retapamulin | 0.05 | 1 |
| bacteroides_fragilis | fusidic_a | 0.05 | 1 |
| bacteroides_fragilis | metronidazole | 0.95 | 1 |
| bacteroides_fragilis | fidaxomicin | 0.1 | 1 |
| bacteroides_fragilis | furazolidone | 0.05 | 1 |
| bacteroides_fragilis | rifampicin | 0.2 | 1 |
| bacteroides_fragilis | amoxicillin_clavulanate | 0.75 | 1 |
| bacteroides_fragilis | piperacillin_tazobactam | 0.85 | 1 |
| bacteroides_fragilis | ampicillin_sulbactam | 0.75 | 1 |
| bacteroides_fragilis | ticarcillin_clavulanate | 0.8 | 1 |
| bacteroides_fragilis | ceftazidime_avibactam | 0.5 | 0.005 |
| bacteroides_fragilis | meropenem_vaborbactam | 0.95 | 0.005 |
| bacteroides_fragilis | colistin | 0.05 | 0.005 |
| bacteroides_fragilis | flucloxacillin | 0.01 | 1 |
| bacteroides_fragilis | aztreonam_avibactam | 0.01 | 1 |
| bacteroides_fragilis | cefixime | 0.45 | 1 |
| campylobacter_jejuni | sulfanilamide | 0.1 | 1 |
| campylobacter_jejuni | penicillin_g | 0.1 | 1 |
| campylobacter_jejuni | ampicillin | 0.1 | 1 |
| campylobacter_jejuni | amoxicillin | 0.1 | 1 |
| campylobacter_jejuni | piperacillin | 0.1 | 1 |
| campylobacter_jejuni | ticarcillin | 0.1 | 1 |
| campylobacter_jejuni | cephalexin | 0.1 | 1 |
| campylobacter_jejuni | cefazolin | 0.1 | 1 |
| campylobacter_jejuni | cefuroxime | 0.1 | 1 |
| campylobacter_jejuni | ceftriaxone | 0.1 | 1 |
| campylobacter_jejuni | ceftazidime | 0.1 | 1 |
| campylobacter_jejuni | cefepime | 0.1 | 1 |
| campylobacter_jejuni | ceftaroline | 0.1 | 1 |
| campylobacter_jejuni | ceftolozane_tazobactam | 0.75 | 1 |
| campylobacter_jejuni | cefiderocol | 0.75 | 1 |
| campylobacter_jejuni | meropenem | 0.1 | 0.005 |
| campylobacter_jejuni | imipenem_c | 0.1 | 0.005 |
| campylobacter_jejuni | ertapenem | 0.1 | 0.005 |
| campylobacter_jejuni | aztreonam | 0.1 | 1 |
| campylobacter_jejuni | erythromycin | 0.85 | 5 |
| campylobacter_jejuni | azithromycin | 0.9 | 4.5 |
| campylobacter_jejuni | clarithromycin | 0.85 | 1 |
| campylobacter_jejuni | clindamycin | 0.7 | 1 |
| campylobacter_jejuni | gentamicin | 0.7 | 1 |
| campylobacter_jejuni | tobramycin | 0.7 | 1 |
| campylobacter_jejuni | amikacin | 0.7 | 1 |
| campylobacter_jejuni | ciprofloxacin | 0.8 | 1 |
| campylobacter_jejuni | levofloxacin | 0.75 | 1 |
| campylobacter_jejuni | moxifloxacin | 0.7 | 1 |
| campylobacter_jejuni | ofloxacin | 0.75 | 1 |
| campylobacter_jejuni | tetracycline | 0.75 | 1 |
| campylobacter_jejuni | doxycycline | 0.8 | 1 |
| campylobacter_jejuni | minocycline | 0.8 | 1 |
| campylobacter_jejuni | tigecycline | 0.7 | 1 |
| campylobacter_jejuni | vancomycin | 0.1 | 1 |
| campylobacter_jejuni | teicoplanin | 0.1 | 1 |
| campylobacter_jejuni | dalbavancin | 0.1 | 0.005 |
| campylobacter_jejuni | linezolid | 0.1 | 0.005 |
| campylobacter_jejuni | tedizolid | 0.1 | 0.005 |
| campylobacter_jejuni | daptomycin | 0.1 | 1 |
| campylobacter_jejuni | quinu_dalfo | 0.1 | 0.005 |
| campylobacter_jejuni | trim_sulf | 0.1 | 1 |
| campylobacter_jejuni | chloramphenicol | 0.7 | 1 |
| campylobacter_jejuni | nitrofurantoin | 0.1 | 1 |
| campylobacter_jejuni | fosfomycin | 0.1 | 1 |
| campylobacter_jejuni | retapamulin | 0.05 | 1 |
| campylobacter_jejuni | fusidic_a | 0.05 | 1 |
| campylobacter_jejuni | metronidazole | 0.05 | 1 |
| campylobacter_jejuni | fidaxomicin | 0.1 | 1 |
| campylobacter_jejuni | furazolidone | 0.1 | 1 |
| campylobacter_jejuni | rifampicin | 0.1 | 1 |
| campylobacter_jejuni | amoxicillin_clavulanate | 0.1 | 1 |
| campylobacter_jejuni | piperacillin_tazobactam | 0.1 | 1 |
| campylobacter_jejuni | ampicillin_sulbactam | 0.1 | 1 |
| campylobacter_jejuni | ticarcillin_clavulanate | 0.1 | 1 |
| campylobacter_jejuni | ceftazidime_avibactam | 0.1 | 0.005 |
| campylobacter_jejuni | meropenem_vaborbactam | 0.1 | 0.005 |
| campylobacter_jejuni | colistin | 0.05 | 0.005 |
| campylobacter_jejuni | flucloxacillin | 0.01 | 1 |
| campylobacter_jejuni | aztreonam_avibactam | 0.9 | 1 |
| campylobacter_jejuni | cefixime | 0.75 | 1 |
| enterobacter_cloacae | sulfanilamide | 0.5 | 1 |
| enterobacter_cloacae | penicillin_g | 0.1 | 1 |
| enterobacter_cloacae | ampicillin | 0.5 | 1 |
| enterobacter_cloacae | amoxicillin | 0.5 | 1 |
| enterobacter_cloacae | piperacillin | 0.75 | 1 |
| enterobacter_cloacae | ticarcillin | 0.7 | 1 |
| enterobacter_cloacae | cephalexin | 0.5 | 1 |
| enterobacter_cloacae | cefazolin | 0.5 | 1 |
| enterobacter_cloacae | cefuroxime | 0.6 | 1 |
| enterobacter_cloacae | ceftriaxone | 0.4 | 1 |
| enterobacter_cloacae | ceftazidime | 0.8 | 1 |
| enterobacter_cloacae | cefepime | 0.85 | 1 |
| enterobacter_cloacae | ceftaroline | 0.4 | 1 |
| enterobacter_cloacae | ceftolozane_tazobactam | 0.8 | 1 |
| enterobacter_cloacae | cefiderocol | 0.8 | 1 |
| enterobacter_cloacae | meropenem | 0.95 | 0.005 |
| enterobacter_cloacae | imipenem_c | 0.95 | 0.005 |
| enterobacter_cloacae | ertapenem | 0.9 | 0.005 |
| enterobacter_cloacae | aztreonam | 0.8 | 1 |
| enterobacter_cloacae | erythromycin | 0.1 | 1 |
| enterobacter_cloacae | azithromycin | 0.1 | 1 |
| enterobacter_cloacae | clarithromycin | 0.1 | 1 |
| enterobacter_cloacae | clindamycin | 0.1 | 1 |
| enterobacter_cloacae | gentamicin | 0.85 | 1 |
| enterobacter_cloacae | tobramycin | 0.8 | 1 |
| enterobacter_cloacae | amikacin | 0.9 | 1 |
| enterobacter_cloacae | ciprofloxacin | 0.9 | 1 |
| enterobacter_cloacae | levofloxacin | 0.85 | 1 |
| enterobacter_cloacae | moxifloxacin | 0.7 | 1 |
| enterobacter_cloacae | ofloxacin | 0.8 | 1 |
| enterobacter_cloacae | tetracycline | 0.8 | 1 |
| enterobacter_cloacae | doxycycline | 0.85 | 1 |
| enterobacter_cloacae | minocycline | 0.85 | 1 |
| enterobacter_cloacae | tigecycline | 0.1 | 1 |
| enterobacter_cloacae | vancomycin | 0.1 | 1 |
| enterobacter_cloacae | teicoplanin | 0.1 | 1 |
| enterobacter_cloacae | dalbavancin | 0.1 | 0.005 |
| enterobacter_cloacae | linezolid | 0.1 | 0.005 |
| enterobacter_cloacae | tedizolid | 0.1 | 0.005 |
| enterobacter_cloacae | daptomycin | 0.1 | 1 |
| enterobacter_cloacae | quinu_dalfo | 0.1 | 0.005 |
| enterobacter_cloacae | trim_sulf | 0.85 | 1 |
| enterobacter_cloacae | chloramphenicol | 0.8 | 1 |
| enterobacter_cloacae | nitrofurantoin | 0.7 | 1 |
| enterobacter_cloacae | fosfomycin | 0.1 | 1 |
| enterobacter_cloacae | retapamulin | 0.05 | 1 |
| enterobacter_cloacae | fusidic_a | 0.05 | 1 |
| enterobacter_cloacae | metronidazole | 0.05 | 1 |
| enterobacter_cloacae | fidaxomicin | 0.1 | 1 |
| enterobacter_cloacae | furazolidone | 0.1 | 1 |
| enterobacter_cloacae | rifampicin | 0.6 | 1 |
| enterobacter_cloacae | amoxicillin_clavulanate | 0.7 | 1 |
| enterobacter_cloacae | piperacillin_tazobactam | 0.85 | 1 |
| enterobacter_cloacae | ampicillin_sulbactam | 0.7 | 1 |
| enterobacter_cloacae | ticarcillin_clavulanate | 0.8 | 1 |
| enterobacter_cloacae | ceftazidime_avibactam | 0.9 | 0.005 |
| enterobacter_cloacae | meropenem_vaborbactam | 0.95 | 0.005 |
| enterobacter_cloacae | colistin | 0.7 | 0.005 |
| enterobacter_cloacae | flucloxacillin | 0.01 | 1 |
| enterobacter_cloacae | aztreonam_avibactam | 1 | 1 |
| enterobacter_cloacae | cefixime | 0.8 | 1 |
| yersinia_enterocolitica | sulfanilamide | 0.5 | 1 |
| yersinia_enterocolitica | penicillin_g | 0.1 | 1 |
| yersinia_enterocolitica | ampicillin | 0.7 | 1 |
| yersinia_enterocolitica | amoxicillin | 0.7 | 1 |
| yersinia_enterocolitica | piperacillin | 0.75 | 1 |
| yersinia_enterocolitica | ticarcillin | 0.7 | 1 |
| yersinia_enterocolitica | cephalexin | 0.6 | 1 |
| yersinia_enterocolitica | cefazolin | 0.65 | 1 |
| yersinia_enterocolitica | cefuroxime | 0.7 | 1 |
| yersinia_enterocolitica | ceftriaxone | 0.9 | 1 |
| yersinia_enterocolitica | ceftazidime | 0.85 | 1 |
| yersinia_enterocolitica | cefepime | 0.85 | 1 |
| yersinia_enterocolitica | ceftaroline | 0.6 | 1 |
| yersinia_enterocolitica | ceftolozane_tazobactam | 0.75 | 1 |
| yersinia_enterocolitica | cefiderocol | 0.75 | 1 |
| yersinia_enterocolitica | meropenem | 0.95 | 0.005 |
| yersinia_enterocolitica | imipenem_c | 0.95 | 0.005 |
| yersinia_enterocolitica | ertapenem | 0.95 | 0.005 |
| yersinia_enterocolitica | aztreonam | 0.85 | 1 |
| yersinia_enterocolitica | erythromycin | 0.1 | 1 |
| yersinia_enterocolitica | azithromycin | 0.1 | 1 |
| yersinia_enterocolitica | clarithromycin | 0.1 | 1 |
| yersinia_enterocolitica | clindamycin | 0.1 | 1 |
| yersinia_enterocolitica | gentamicin | 0.85 | 1 |
| yersinia_enterocolitica | tobramycin | 0.8 | 1 |
| yersinia_enterocolitica | amikacin | 0.9 | 1 |
| yersinia_enterocolitica | ciprofloxacin | 0.9 | 1 |
| yersinia_enterocolitica | levofloxacin | 0.85 | 1 |
| yersinia_enterocolitica | moxifloxacin | 0.7 | 1 |
| yersinia_enterocolitica | ofloxacin | 0.8 | 1 |
| yersinia_enterocolitica | tetracycline | 0.8 | 1 |
| yersinia_enterocolitica | doxycycline | 0.85 | 1 |
| yersinia_enterocolitica | minocycline | 0.85 | 1 |
| yersinia_enterocolitica | tigecycline | 0.7 | 1 |
| yersinia_enterocolitica | vancomycin | 0.1 | 1 |
| yersinia_enterocolitica | teicoplanin | 0.1 | 1 |
| yersinia_enterocolitica | dalbavancin | 0.1 | 0.005 |
| yersinia_enterocolitica | linezolid | 0.1 | 0.005 |
| yersinia_enterocolitica | tedizolid | 0.1 | 0.005 |
| yersinia_enterocolitica | daptomycin | 0.1 | 1 |
| yersinia_enterocolitica | quinu_dalfo | 0.1 | 0.005 |
| yersinia_enterocolitica | trim_sulf | 0.95 | 1 |
| yersinia_enterocolitica | chloramphenicol | 0.85 | 1 |
| yersinia_enterocolitica | nitrofurantoin | 0.1 | 1 |
| yersinia_enterocolitica | fosfomycin | 0.1 | 1 |
| yersinia_enterocolitica | retapamulin | 0.05 | 1 |
| yersinia_enterocolitica | fusidic_a | 0.05 | 1 |
| yersinia_enterocolitica | metronidazole | 0.05 | 1 |
| yersinia_enterocolitica | fidaxomicin | 0.1 | 1 |
| yersinia_enterocolitica | furazolidone | 0.1 | 1 |
| yersinia_enterocolitica | rifampicin | 0.7 | 1 |
| yersinia_enterocolitica | amoxicillin_clavulanate | 0.85 | 1 |
| yersinia_enterocolitica | piperacillin_tazobactam | 0.85 | 1 |
| yersinia_enterocolitica | ampicillin_sulbactam | 0.8 | 1 |
| yersinia_enterocolitica | ticarcillin_clavulanate | 0.8 | 1 |
| yersinia_enterocolitica | ceftazidime_avibactam | 0.95 | 0.005 |
| yersinia_enterocolitica | meropenem_vaborbactam | 0.95 | 0.005 |
| yersinia_enterocolitica | colistin | 0.7 | 0.005 |
| yersinia_enterocolitica | flucloxacillin | 0.01 | 1 |
| yersinia_enterocolitica | aztreonam_avibactam | 0.9 | 1 |
| yersinia_enterocolitica | cefixime | 0.75 | 1 |
| moraxella_catarrhalis | sulfanilamide | 0.1 | 1 |
| moraxella_catarrhalis | penicillin_g | 0.9 | 1 |
| moraxella_catarrhalis | ampicillin | 0.9 | 1 |
| moraxella_catarrhalis | amoxicillin | 0.9 | 1 |
| moraxella_catarrhalis | piperacillin | 0.8 | 1 |
| moraxella_catarrhalis | ticarcillin | 0.8 | 1 |
| moraxella_catarrhalis | cephalexin | 0.8 | 1 |
| moraxella_catarrhalis | cefazolin | 0.85 | 1 |
| moraxella_catarrhalis | cefuroxime | 0.9 | 1 |
| moraxella_catarrhalis | ceftriaxone | 0.95 | 1 |
| moraxella_catarrhalis | ceftazidime | 0.9 | 1 |
| moraxella_catarrhalis | cefepime | 0.9 | 1 |
| moraxella_catarrhalis | ceftaroline | 0.8 | 1 |
| moraxella_catarrhalis | ceftolozane_tazobactam | 0.8 | 1 |
| moraxella_catarrhalis | cefiderocol | 0.8 | 1 |
| moraxella_catarrhalis | meropenem | 0.95 | 0.005 |
| moraxella_catarrhalis | imipenem_c | 0.95 | 0.005 |
| moraxella_catarrhalis | ertapenem | 0.95 | 0.005 |
| moraxella_catarrhalis | aztreonam | 0.9 | 1 |
| moraxella_catarrhalis | erythromycin | 0.8 | 1 |
| moraxella_catarrhalis | azithromycin | 0.85 | 1 |
| moraxella_catarrhalis | clarithromycin | 0.8 | 1 |
| moraxella_catarrhalis | clindamycin | 0.1 | 1 |
| moraxella_catarrhalis | gentamicin | 0.1 | 1 |
| moraxella_catarrhalis | tobramycin | 0.1 | 1 |
| moraxella_catarrhalis | amikacin | 0.1 | 1 |
| moraxella_catarrhalis | ciprofloxacin | 0.9 | 1 |
| moraxella_catarrhalis | levofloxacin | 0.85 | 1 |
| moraxella_catarrhalis | moxifloxacin | 0.8 | 1 |
| moraxella_catarrhalis | ofloxacin | 0.85 | 1 |
| moraxella_catarrhalis | tetracycline | 0.8 | 1 |
| moraxella_catarrhalis | doxycycline | 0.8 | 1 |
| moraxella_catarrhalis | minocycline | 0.85 | 1 |
| moraxella_catarrhalis | tigecycline | 0.1 | 1 |
| moraxella_catarrhalis | vancomycin | 0.1 | 1 |
| moraxella_catarrhalis | teicoplanin | 0.1 | 1 |
| moraxella_catarrhalis | dalbavancin | 0.1 | 0.005 |
| moraxella_catarrhalis | linezolid | 0.1 | 0.005 |
| moraxella_catarrhalis | tedizolid | 0.1 | 0.005 |
| moraxella_catarrhalis | daptomycin | 0.1 | 1 |
| moraxella_catarrhalis | quinu_dalfo | 0.1 | 0.005 |
| moraxella_catarrhalis | trim_sulf | 0.95 | 1 |
| moraxella_catarrhalis | chloramphenicol | 0.85 | 1 |
| moraxella_catarrhalis | nitrofurantoin | 0.1 | 1 |
| moraxella_catarrhalis | fosfomycin | 0.1 | 1 |
| moraxella_catarrhalis | retapamulin | 0.05 | 1 |
| moraxella_catarrhalis | fusidic_a | 0.05 | 1 |
| moraxella_catarrhalis | metronidazole | 0.05 | 1 |
| moraxella_catarrhalis | fidaxomicin | 0.1 | 1 |
| moraxella_catarrhalis | furazolidone | 0.1 | 1 |
| moraxella_catarrhalis | rifampicin | 0.7 | 1 |
| moraxella_catarrhalis | amoxicillin_clavulanate | 0.95 | 1 |
| moraxella_catarrhalis | piperacillin_tazobactam | 0.85 | 1 |
| moraxella_catarrhalis | ampicillin_sulbactam | 0.95 | 1 |
| moraxella_catarrhalis | ticarcillin_clavulanate | 0.85 | 1 |
| moraxella_catarrhalis | ceftazidime_avibactam | 0.95 | 0.005 |
| moraxella_catarrhalis | meropenem_vaborbactam | 0.95 | 0.005 |
| moraxella_catarrhalis | colistin | 0.05 | 0.005 |
| moraxella_catarrhalis | flucloxacillin | 0.01 | 1 |
| moraxella_catarrhalis | aztreonam_avibactam | 0.8 | 1 |
| moraxella_catarrhalis | cefixime | 0.8 | 1 |
| treponema_pallidum | sulfanilamide | 0.1 | 1 |
| treponema_pallidum | penicillin_g | 1 | 1 |
| treponema_pallidum | ampicillin | 0.95 | 1 |
| treponema_pallidum | amoxicillin | 0.95 | 1 |
| treponema_pallidum | piperacillin | 0.9 | 1 |
| treponema_pallidum | ticarcillin | 0.9 | 1 |
| treponema_pallidum | cephalexin | 0.9 | 1 |
| treponema_pallidum | cefazolin | 0.9 | 1 |
| treponema_pallidum | cefuroxime | 0.95 | 1 |
| treponema_pallidum | ceftriaxone | 0.95 | 1 |
| treponema_pallidum | ceftazidime | 0.9 | 1 |
| treponema_pallidum | cefepime | 0.9 | 1 |
| treponema_pallidum | ceftaroline | 0.9 | 1 |
| treponema_pallidum | ceftolozane_tazobactam | 0.1 | 1 |
| treponema_pallidum | cefiderocol | 0.1 | 1 |
| treponema_pallidum | meropenem | 0.95 | 0.005 |
| treponema_pallidum | imipenem_c | 0.95 | 0.005 |
| treponema_pallidum | ertapenem | 0.95 | 0.005 |
| treponema_pallidum | aztreonam | 0.9 | 1 |
| treponema_pallidum | erythromycin | 0.8 | 1 |
| treponema_pallidum | azithromycin | 0.85 | 1 |
| treponema_pallidum | clarithromycin | 0.8 | 1 |
| treponema_pallidum | clindamycin | 0.1 | 1 |
| treponema_pallidum | gentamicin | 0.1 | 1 |
| treponema_pallidum | tobramycin | 0.1 | 1 |
| treponema_pallidum | amikacin | 0.1 | 1 |
| treponema_pallidum | ciprofloxacin | 0.7 | 1 |
| treponema_pallidum | levofloxacin | 0.75 | 1 |
| treponema_pallidum | moxifloxacin | 0.75 | 1 |
| treponema_pallidum | ofloxacin | 0.7 | 1 |
| treponema_pallidum | tetracycline | 0.8 | 1 |
| treponema_pallidum | doxycycline | 0.8 | 1 |
| treponema_pallidum | minocycline | 0.85 | 1 |
| treponema_pallidum | tigecycline | 0.1 | 1 |
| treponema_pallidum | vancomycin | 0.1 | 1 |
| treponema_pallidum | teicoplanin | 0.1 | 1 |
| treponema_pallidum | dalbavancin | 0.1 | 0.005 |
| treponema_pallidum | linezolid | 0.1 | 0.005 |
| treponema_pallidum | tedizolid | 0.1 | 0.005 |
| treponema_pallidum | daptomycin | 0.1 | 1 |
| treponema_pallidum | quinu_dalfo | 0.1 | 0.005 |
| treponema_pallidum | trim_sulf | 0.1 | 1 |
| treponema_pallidum | chloramphenicol | 0.8 | 1 |
| treponema_pallidum | nitrofurantoin | 0.1 | 1 |
| treponema_pallidum | fosfomycin | 0.1 | 1 |
| treponema_pallidum | retapamulin | 0.05 | 1 |
| treponema_pallidum | fusidic_a | 0.05 | 1 |
| treponema_pallidum | metronidazole | 0.05 | 1 |
| treponema_pallidum | fidaxomicin | 0.1 | 1 |
| treponema_pallidum | furazolidone | 0.1 | 1 |
| treponema_pallidum | rifampicin | 0.1 | 1 |
| treponema_pallidum | amoxicillin_clavulanate | 0.95 | 1 |
| treponema_pallidum | piperacillin_tazobactam | 0.9 | 1 |
| treponema_pallidum | ampicillin_sulbactam | 0.95 | 1 |
| treponema_pallidum | ticarcillin_clavulanate | 0.9 | 1 |
| treponema_pallidum | ceftazidime_avibactam | 0.95 | 0.005 |
| treponema_pallidum | meropenem_vaborbactam | 0.95 | 0.005 |
| treponema_pallidum | colistin | 0.05 | 0.005 |
| treponema_pallidum | flucloxacillin | 0.01 | 1 |
| treponema_pallidum | aztreonam_avibactam | 0.9 | 1 |
| treponema_pallidum | cefixime | 0.1 | 1 |
| bordetella_pertussis | sulfanilamide | 0.1 | 1 |
| bordetella_pertussis | penicillin_g | 0.1 | 1 |
| bordetella_pertussis | ampicillin | 0.1 | 1 |
| bordetella_pertussis | amoxicillin | 0.1 | 1 |
| bordetella_pertussis | piperacillin | 0.1 | 1 |
| bordetella_pertussis | ticarcillin | 0.1 | 1 |
| bordetella_pertussis | cephalexin | 0.1 | 1 |
| bordetella_pertussis | cefazolin | 0.1 | 1 |
| bordetella_pertussis | cefuroxime | 0.1 | 1 |
| bordetella_pertussis | ceftriaxone | 0.1 | 1 |
| bordetella_pertussis | ceftazidime | 0.1 | 1 |
| bordetella_pertussis | cefepime | 0.1 | 1 |
| bordetella_pertussis | ceftaroline | 0.1 | 1 |
| bordetella_pertussis | ceftolozane_tazobactam | 0.8 | 1 |
| bordetella_pertussis | cefiderocol | 0.8 | 1 |
| bordetella_pertussis | meropenem | 0.1 | 0.005 |
| bordetella_pertussis | imipenem_c | 0.1 | 0.005 |
| bordetella_pertussis | ertapenem | 0.1 | 0.005 |
| bordetella_pertussis | aztreonam | 0.1 | 1 |
| bordetella_pertussis | erythromycin | 0.9 | 1 |
| bordetella_pertussis | azithromycin | 0.95 | 1 |
| bordetella_pertussis | clarithromycin | 0.9 | 1 |
| bordetella_pertussis | clindamycin | 0.1 | 1 |
| bordetella_pertussis | gentamicin | 0.7 | 1 |
| bordetella_pertussis | tobramycin | 0.7 | 1 |
| bordetella_pertussis | amikacin | 0.7 | 1 |
| bordetella_pertussis | ciprofloxacin | 0.7 | 1 |
| bordetella_pertussis | levofloxacin | 0.75 | 1 |
| bordetella_pertussis | moxifloxacin | 0.75 | 1 |
| bordetella_pertussis | ofloxacin | 0.7 | 1 |
| bordetella_pertussis | tetracycline | 0.7 | 1 |
| bordetella_pertussis | doxycycline | 0.75 | 1 |
| bordetella_pertussis | minocycline | 0.75 | 1 |
| bordetella_pertussis | tigecycline | 0.1 | 1 |
| bordetella_pertussis | vancomycin | 0.1 | 1 |
| bordetella_pertussis | teicoplanin | 0.1 | 1 |
| bordetella_pertussis | dalbavancin | 0.1 | 0.005 |
| bordetella_pertussis | linezolid | 0.1 | 0.005 |
| bordetella_pertussis | tedizolid | 0.1 | 0.005 |
| bordetella_pertussis | daptomycin | 0.1 | 1 |
| bordetella_pertussis | quinu_dalfo | 0.1 | 0.005 |
| bordetella_pertussis | trim_sulf | 0.7 | 1 |
| bordetella_pertussis | chloramphenicol | 0.8 | 1 |
| bordetella_pertussis | nitrofurantoin | 0.1 | 1 |
| bordetella_pertussis | fosfomycin | 0.1 | 1 |
| bordetella_pertussis | retapamulin | 0.05 | 1 |
| bordetella_pertussis | fusidic_a | 0.05 | 1 |
| bordetella_pertussis | metronidazole | 0.05 | 1 |
| bordetella_pertussis | fidaxomicin | 0.1 | 1 |
| bordetella_pertussis | furazolidone | 0.1 | 1 |
| bordetella_pertussis | rifampicin | 0.7 | 1 |
| bordetella_pertussis | amoxicillin_clavulanate | 0.1 | 1 |
| bordetella_pertussis | piperacillin_tazobactam | 0.1 | 1 |
| bordetella_pertussis | ampicillin_sulbactam | 0.1 | 1 |
| bordetella_pertussis | ticarcillin_clavulanate | 0.1 | 1 |
| bordetella_pertussis | ceftazidime_avibactam | 0.1 | 0.005 |
| bordetella_pertussis | meropenem_vaborbactam | 0.1 | 0.005 |
| bordetella_pertussis | colistin | 0.05 | 0.005 |
| bordetella_pertussis | flucloxacillin | 0.01 | 1 |
| bordetella_pertussis | aztreonam_avibactam | 0.8 | 1 |
| bordetella_pertussis | cefixime | 0.8 | 1 |
| helicobacter_pylori | sulfanilamide | 0.1 | 1 |
| helicobacter_pylori | penicillin_g | 0.1 | 1 |
| helicobacter_pylori | ampicillin | 0.7 | 1 |
| helicobacter_pylori | amoxicillin | 0.85 | 12 |
| helicobacter_pylori | piperacillin | 0.1 | 1 |
| helicobacter_pylori | ticarcillin | 0.1 | 1 |
| helicobacter_pylori | cephalexin | 0.1 | 1 |
| helicobacter_pylori | cefazolin | 0.1 | 1 |
| helicobacter_pylori | cefuroxime | 0.1 | 1 |
| helicobacter_pylori | ceftriaxone | 0.1 | 1 |
| helicobacter_pylori | ceftazidime | 0.1 | 1 |
| helicobacter_pylori | cefepime | 0.1 | 1 |
| helicobacter_pylori | ceftaroline | 0.1 | 1 |
| helicobacter_pylori | ceftolozane_tazobactam | 0.05 | 1 |
| helicobacter_pylori | cefiderocol | 0.05 | 1 |
| helicobacter_pylori | meropenem | 0.1 | 0.005 |
| helicobacter_pylori | imipenem_c | 0.1 | 0.005 |
| helicobacter_pylori | ertapenem | 0.1 | 0.005 |
| helicobacter_pylori | aztreonam | 0.1 | 1 |
| helicobacter_pylori | erythromycin | 0.8 | 1 |
| helicobacter_pylori | azithromycin | 0.85 | 1 |
| helicobacter_pylori | clarithromycin | 0.8 | 15 |
| helicobacter_pylori | clindamycin | 0.1 | 1 |
| helicobacter_pylori | gentamicin | 0.1 | 1 |
| helicobacter_pylori | tobramycin | 0.1 | 1 |
| helicobacter_pylori | amikacin | 0.1 | 1 |
| helicobacter_pylori | ciprofloxacin | 0.7 | 1 |
| helicobacter_pylori | levofloxacin | 0.7 | 5 |
| helicobacter_pylori | moxifloxacin | 0.75 | 1 |
| helicobacter_pylori | ofloxacin | 0.7 | 1 |
| helicobacter_pylori | tetracycline | 0.8 | 6 |
| helicobacter_pylori | doxycycline | 0.8 | 1 |
| helicobacter_pylori | minocycline | 0.85 | 1 |
| helicobacter_pylori | tigecycline | 0.1 | 1 |
| helicobacter_pylori | vancomycin | 0.1 | 1 |
| helicobacter_pylori | teicoplanin | 0.1 | 1 |
| helicobacter_pylori | dalbavancin | 0.1 | 0.005 |
| helicobacter_pylori | linezolid | 0.1 | 0.005 |
| helicobacter_pylori | tedizolid | 0.1 | 0.005 |
| helicobacter_pylori | daptomycin | 0.1 | 1 |
| helicobacter_pylori | quinu_dalfo | 0.1 | 0.005 |
| helicobacter_pylori | trim_sulf | 0.1 | 1 |
| helicobacter_pylori | chloramphenicol | 0.7 | 1 |
| helicobacter_pylori | nitrofurantoin | 0.1 | 1 |
| helicobacter_pylori | fosfomycin | 0.1 | 1 |
| helicobacter_pylori | retapamulin | 0.05 | 1 |
| helicobacter_pylori | fusidic_a | 0.05 | 1 |
| helicobacter_pylori | metronidazole | 0.8 | 8 |
| helicobacter_pylori | fidaxomicin | 0.1 | 1 |
| helicobacter_pylori | furazolidone | 0.1 | 1 |
| helicobacter_pylori | rifampicin | 0.1 | 1 |
| helicobacter_pylori | amoxicillin_clavulanate | 0.85 | 1 |
| helicobacter_pylori | piperacillin_tazobactam | 0.1 | 1 |
| helicobacter_pylori | ampicillin_sulbactam | 0.7 | 1 |
| helicobacter_pylori | ticarcillin_clavulanate | 0.1 | 1 |
| helicobacter_pylori | ceftazidime_avibactam | 0.1 | 0.005 |
| helicobacter_pylori | meropenem_vaborbactam | 0.1 | 0.005 |
| helicobacter_pylori | colistin | 0.05 | 0.005 |
| helicobacter_pylori | flucloxacillin | 0.01 | 1 |
| helicobacter_pylori | aztreonam_avibactam | 0.01 | 1 |
| helicobacter_pylori | cefixime | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | sulfanilamide | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | penicillin_g | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ampicillin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | amoxicillin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | piperacillin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ticarcillin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | cephalexin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | cefazolin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | cefuroxime | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ceftriaxone | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ceftazidime | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | cefepime | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ceftaroline | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ceftolozane_tazobactam | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | cefiderocol | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | meropenem | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | imipenem_c | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | ertapenem | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | aztreonam | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | erythromycin | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | azithromycin | 0.25 | 1 |
| mdr_mycobacterium_tuberculosis | clarithromycin | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | clindamycin | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | gentamicin | 0.25 | 1 |
| mdr_mycobacterium_tuberculosis | tobramycin | 0.25 | 1 |
| mdr_mycobacterium_tuberculosis | amikacin | 0.3 | 1 |
| mdr_mycobacterium_tuberculosis | ciprofloxacin | 0.4 | 1 |
| mdr_mycobacterium_tuberculosis | levofloxacin | 0.45 | 1 |
| mdr_mycobacterium_tuberculosis | moxifloxacin | 0.45 | 1 |
| mdr_mycobacterium_tuberculosis | ofloxacin | 0.4 | 1 |
| mdr_mycobacterium_tuberculosis | tetracycline | 0.3 | 1 |
| mdr_mycobacterium_tuberculosis | doxycycline | 0.35 | 1 |
| mdr_mycobacterium_tuberculosis | minocycline | 0.35 | 1 |
| mdr_mycobacterium_tuberculosis | tigecycline | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | vancomycin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | teicoplanin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | dalbavancin | 0.1 | 0.005 |
| mdr_mycobacterium_tuberculosis | linezolid | 0.1 | 0.005 |
| mdr_mycobacterium_tuberculosis | tedizolid | 0.1 | 0.005 |
| mdr_mycobacterium_tuberculosis | daptomycin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | quinu_dalfo | 0.1 | 0.005 |
| mdr_mycobacterium_tuberculosis | trim_sulf | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | chloramphenicol | 0.2 | 1 |
| mdr_mycobacterium_tuberculosis | nitrofurantoin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | fosfomycin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | retapamulin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | fusidic_a | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | metronidazole | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | fidaxomicin | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | furazolidone | 0.1 | 1 |
| mdr_mycobacterium_tuberculosis | rifampicin | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | amoxicillin_clavulanate | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | piperacillin_tazobactam | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ampicillin_sulbactam | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ticarcillin_clavulanate | 0.05 | 1 |
| mdr_mycobacterium_tuberculosis | ceftazidime_avibactam | 0.05 | 0.005 |
| mdr_mycobacterium_tuberculosis | meropenem_vaborbactam | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | colistin | 0.2 | 0.005 |
| mdr_mycobacterium_tuberculosis | flucloxacillin | 0.01 | 1 |
| mdr_mycobacterium_tuberculosis | aztreonam_avibactam | 0.9 | 1 |
| mdr_mycobacterium_tuberculosis | cefixime | 0.1 | 1 |
| mycoplasma_pneumoniae | sulfanilamide | 0.05 | 1 |
| mycoplasma_pneumoniae | penicillin_g | 0.05 | 0.001 |
| mycoplasma_pneumoniae | ampicillin | 0.05 | 0.001 |
| mycoplasma_pneumoniae | amoxicillin | 0.05 | 0.001 |
| mycoplasma_pneumoniae | piperacillin | 0.05 | 1 |
| mycoplasma_pneumoniae | ticarcillin | 0.05 | 1 |
| mycoplasma_pneumoniae | cephalexin | 0.05 | 0.001 |
| mycoplasma_pneumoniae | cefazolin | 0.05 | 0.001 |
| mycoplasma_pneumoniae | cefuroxime | 0.05 | 1 |
| mycoplasma_pneumoniae | ceftriaxone | 0.05 | 0.001 |
| mycoplasma_pneumoniae | ceftazidime | 0.05 | 1 |
| mycoplasma_pneumoniae | cefepime | 0.05 | 1 |
| mycoplasma_pneumoniae | ceftaroline | 0.05 | 1 |
| mycoplasma_pneumoniae | ceftolozane_tazobactam | 0.01 | 1 |
| mycoplasma_pneumoniae | cefiderocol | 0.01 | 1 |
| mycoplasma_pneumoniae | meropenem | 0.05 | 0.001 |
| mycoplasma_pneumoniae | imipenem_c | 0.05 | 0.005 |
| mycoplasma_pneumoniae | ertapenem | 0.05 | 0.001 |
| mycoplasma_pneumoniae | aztreonam | 0.05 | 1 |
| mycoplasma_pneumoniae | erythromycin | 0.8 | 1 |
| mycoplasma_pneumoniae | azithromycin | 0.85 | 1 |
| mycoplasma_pneumoniae | clarithromycin | 0.8 | 1 |
| mycoplasma_pneumoniae | clindamycin | 0.05 | 1 |
| mycoplasma_pneumoniae | gentamicin | 0.05 | 1 |
| mycoplasma_pneumoniae | tobramycin | 0.05 | 1 |
| mycoplasma_pneumoniae | amikacin | 0.05 | 1 |
| mycoplasma_pneumoniae | ciprofloxacin | 0.7 | 1 |
| mycoplasma_pneumoniae | levofloxacin | 0.75 | 1 |
| mycoplasma_pneumoniae | moxifloxacin | 0.8 | 1 |
| mycoplasma_pneumoniae | ofloxacin | 0.6 | 1 |
| mycoplasma_pneumoniae | tetracycline | 0.7 | 1 |
| mycoplasma_pneumoniae | doxycycline | 0.75 | 1 |
| mycoplasma_pneumoniae | minocycline | 0.8 | 1 |
| mycoplasma_pneumoniae | tigecycline | 0.85 | 1 |
| mycoplasma_pneumoniae | vancomycin | 0.05 | 1 |
| mycoplasma_pneumoniae | teicoplanin | 0.05 | 1 |
| mycoplasma_pneumoniae | dalbavancin | 0.05 | 0.005 |
| mycoplasma_pneumoniae | linezolid | 0.05 | 0.005 |
| mycoplasma_pneumoniae | tedizolid | 0.05 | 0.005 |
| mycoplasma_pneumoniae | daptomycin | 0.1 | 1 |
| mycoplasma_pneumoniae | quinu_dalfo | 0.05 | 0.005 |
| mycoplasma_pneumoniae | trim_sulf | 0.05 | 1 |
| mycoplasma_pneumoniae | chloramphenicol | 0.05 | 1 |
| mycoplasma_pneumoniae | nitrofurantoin | 0.05 | 1 |
| mycoplasma_pneumoniae | fosfomycin | 0.1 | 1 |
| mycoplasma_pneumoniae | retapamulin | 0.05 | 1 |
| mycoplasma_pneumoniae | fusidic_a | 0.05 | 1 |
| mycoplasma_pneumoniae | metronidazole | 0.05 | 1 |
| mycoplasma_pneumoniae | fidaxomicin | 0.1 | 1 |
| mycoplasma_pneumoniae | furazolidone | 0.05 | 1 |
| mycoplasma_pneumoniae | rifampicin | 0.05 | 1 |
| mycoplasma_pneumoniae | amoxicillin_clavulanate | 0.05 | 1 |
| mycoplasma_pneumoniae | piperacillin_tazobactam | 0.05 | 1 |
| mycoplasma_pneumoniae | ampicillin_sulbactam | 0.05 | 1 |
| mycoplasma_pneumoniae | ticarcillin_clavulanate | 0.05 | 1 |
| mycoplasma_pneumoniae | ceftazidime_avibactam | 0.05 | 0.005 |
| mycoplasma_pneumoniae | meropenem_vaborbactam | 0.05 | 0.005 |
| mycoplasma_pneumoniae | colistin | 0.05 | 0.005 |
| mycoplasma_pneumoniae | flucloxacillin | 0.01 | 1 |
| mycoplasma_pneumoniae | aztreonam_avibactam | 0.01 | 1 |
| mycoplasma_pneumoniae | cefixime | 0.01 | 1 |
| legionella_pneumophila | sulfanilamide | 0.05 | 1 |
| legionella_pneumophila | penicillin_g | 0.05 | 0.001 |
| legionella_pneumophila | ampicillin | 0.05 | 0.001 |
| legionella_pneumophila | amoxicillin | 0.05 | 0.001 |
| legionella_pneumophila | piperacillin | 0.05 | 1 |
| legionella_pneumophila | ticarcillin | 0.05 | 1 |
| legionella_pneumophila | cephalexin | 0.05 | 0.001 |
| legionella_pneumophila | cefazolin | 0.05 | 0.001 |
| legionella_pneumophila | cefuroxime | 0.05 | 1 |
| legionella_pneumophila | ceftriaxone | 0.05 | 0.001 |
| legionella_pneumophila | ceftazidime | 0.05 | 1 |
| legionella_pneumophila | cefepime | 0.05 | 1 |
| legionella_pneumophila | ceftaroline | 0.05 | 1 |
| legionella_pneumophila | ceftolozane_tazobactam | 0.8 | 1 |
| legionella_pneumophila | cefiderocol | 0.8 | 1 |
| legionella_pneumophila | meropenem | 0.05 | 0.001 |
| legionella_pneumophila | imipenem_c | 0.05 | 0.005 |
| legionella_pneumophila | ertapenem | 0.05 | 0.001 |
| legionella_pneumophila | aztreonam | 0.8 | 1 |
| legionella_pneumophila | erythromycin | 0.8 | 1 |
| legionella_pneumophila | azithromycin | 0.9 | 1 |
| legionella_pneumophila | clarithromycin | 0.8 | 1 |
| legionella_pneumophila | clindamycin | 0.05 | 1 |
| legionella_pneumophila | gentamicin | 0.05 | 1 |
| legionella_pneumophila | tobramycin | 0.05 | 1 |
| legionella_pneumophila | amikacin | 0.05 | 1 |
| legionella_pneumophila | ciprofloxacin | 0.9 | 1 |
| legionella_pneumophila | levofloxacin | 0.95 | 1 |
| legionella_pneumophila | moxifloxacin | 0.9 | 1 |
| legionella_pneumophila | ofloxacin | 0.7 | 1 |
| legionella_pneumophila | tetracycline | 0.8 | 1 |
| legionella_pneumophila | doxycycline | 0.85 | 1 |
| legionella_pneumophila | minocycline | 0.9 | 1 |
| legionella_pneumophila | tigecycline | 0.1 | 1 |
| legionella_pneumophila | vancomycin | 0.05 | 1 |
| legionella_pneumophila | teicoplanin | 0.05 | 1 |
| legionella_pneumophila | dalbavancin | 0.05 | 0.005 |
| legionella_pneumophila | linezolid | 0.05 | 0.005 |
| legionella_pneumophila | tedizolid | 0.05 | 0.005 |
| legionella_pneumophila | daptomycin | 0.1 | 1 |
| legionella_pneumophila | quinu_dalfo | 0.05 | 0.005 |
| legionella_pneumophila | trim_sulf | 0.05 | 1 |
| legionella_pneumophila | chloramphenicol | 0.05 | 1 |
| legionella_pneumophila | nitrofurantoin | 0.05 | 1 |
| legionella_pneumophila | fosfomycin | 0.1 | 1 |
| legionella_pneumophila | retapamulin | 0.05 | 1 |
| legionella_pneumophila | fusidic_a | 0.05 | 1 |
| legionella_pneumophila | metronidazole | 0.05 | 1 |
| legionella_pneumophila | fidaxomicin | 0.1 | 1 |
| legionella_pneumophila | furazolidone | 0.05 | 1 |
| legionella_pneumophila | rifampicin | 0.05 | 1 |
| legionella_pneumophila | amoxicillin_clavulanate | 0.05 | 1 |
| legionella_pneumophila | piperacillin_tazobactam | 0.05 | 1 |
| legionella_pneumophila | ampicillin_sulbactam | 0.05 | 1 |
| legionella_pneumophila | ticarcillin_clavulanate | 0.05 | 1 |
| legionella_pneumophila | ceftazidime_avibactam | 0.05 | 0.005 |
| legionella_pneumophila | meropenem_vaborbactam | 0.05 | 0.005 |
| legionella_pneumophila | colistin | 0.05 | 0.005 |
| legionella_pneumophila | flucloxacillin | 0.01 | 1 |
| legionella_pneumophila | aztreonam_avibactam | 0.01 | 1 |
| legionella_pneumophila | cefixime | 0.8 | 1 |
| burkholderia_cepacia_complex | sulfanilamide | 0.1 | 1 |
| burkholderia_cepacia_complex | penicillin_g | 0.05 | 1 |
| burkholderia_cepacia_complex | ampicillin | 0.05 | 1 |
| burkholderia_cepacia_complex | amoxicillin | 0.05 | 1 |
| burkholderia_cepacia_complex | piperacillin | 0.6 | 1 |
| burkholderia_cepacia_complex | ticarcillin | 0.5 | 1 |
| burkholderia_cepacia_complex | cephalexin | 0.05 | 1 |
| burkholderia_cepacia_complex | cefazolin | 0.05 | 1 |
| burkholderia_cepacia_complex | cefuroxime | 0.1 | 1 |
| burkholderia_cepacia_complex | ceftriaxone | 0.1 | 1 |
| burkholderia_cepacia_complex | ceftazidime | 0.7 | 1 |
| burkholderia_cepacia_complex | cefepime | 0.75 | 1 |
| burkholderia_cepacia_complex | ceftaroline | 0.1 | 1 |
| burkholderia_cepacia_complex | ceftolozane_tazobactam | 0.1 | 1 |
| burkholderia_cepacia_complex | cefiderocol | 0.1 | 1 |
| burkholderia_cepacia_complex | meropenem | 0.8 | 0.005 |
| burkholderia_cepacia_complex | imipenem_c | 0.8 | 0.005 |
| burkholderia_cepacia_complex | ertapenem | 0.1 | 0.005 |
| burkholderia_cepacia_complex | aztreonam | 0.1 | 1 |
| burkholderia_cepacia_complex | gentamicin | 0.7 | 1 |
| burkholderia_cepacia_complex | tobramycin | 0.65 | 1 |
| burkholderia_cepacia_complex | amikacin | 0.75 | 1 |
| burkholderia_cepacia_complex | ciprofloxacin | 0.6 | 1 |
| burkholderia_cepacia_complex | levofloxacin | 0.65 | 1 |
| burkholderia_cepacia_complex | moxifloxacin | 0.6 | 1 |
| burkholderia_cepacia_complex | ofloxacin | 0.6 | 1 |
| burkholderia_cepacia_complex | tetracycline | 0.6 | 1 |
| burkholderia_cepacia_complex | doxycycline | 0.65 | 1 |
| burkholderia_cepacia_complex | minocycline | 0.7 | 1 |
| burkholderia_cepacia_complex | tigecycline | 0.1 | 1 |
| burkholderia_cepacia_complex | dalbavancin | 0 | 0.005 |
| burkholderia_cepacia_complex | linezolid | 0 | 0.005 |
| burkholderia_cepacia_complex | tedizolid | 0 | 0.005 |
| burkholderia_cepacia_complex | daptomycin | 0.1 | 1 |
| burkholderia_cepacia_complex | quinu_dalfo | 0 | 0.005 |
| burkholderia_cepacia_complex | trim_sulf | 0.6 | 1 |
| burkholderia_cepacia_complex | chloramphenicol | 0.7 | 1 |
| burkholderia_cepacia_complex | nitrofurantoin | 0.1 | 1 |
| burkholderia_cepacia_complex | fosfomycin | 0.1 | 1 |
| burkholderia_cepacia_complex | fidaxomicin | 0.1 | 1 |
| burkholderia_cepacia_complex | furazolidone | 0.1 | 1 |
| burkholderia_cepacia_complex | rifampicin | 0.5 | 1 |
| burkholderia_cepacia_complex | amoxicillin_clavulanate | 0.05 | 1 |
| burkholderia_cepacia_complex | piperacillin_tazobactam | 0.65 | 1 |
| burkholderia_cepacia_complex | ampicillin_sulbactam | 0.65 | 1 |
| burkholderia_cepacia_complex | ticarcillin_clavulanate | 0.6 | 1 |
| burkholderia_cepacia_complex | ceftazidime_avibactam | 0.65 | 0.005 |
| burkholderia_cepacia_complex | meropenem_vaborbactam | 0.75 | 0.005 |
| burkholderia_cepacia_complex | colistin | 0.8 | 0.005 |
| burkholderia_cepacia_complex | flucloxacillin | 0.01 | 1 |
| burkholderia_cepacia_complex | aztreonam_avibactam | 0.6 | 1 |
| burkholderia_cepacia_complex | cefixime | 0.1 | 1 |

### B.5 Regional Parameters

Region-level parameters that apply across bacteria, plus the per-region per-bacteria acquisition log-odds adjustments.

See: [§2.5 Travel](#25-travel), [§3.1 Community acquisition](#31-community-acquisition).

#### Region-Level Parameters

| Region | Travel mult | Cessation mult | Mortality log-odds | Sepsis log-odds | Sepsis mort mult | Testing mult | Abx init log-odds | Hosp log-odds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| north_america | 3 | 0.85 | 0 | 0.4 | 0.5 | 1.1 | 0 | 0.5 |
| south_america | 0.8 | 1.25 | 0.26 | -0.3 | 1.1 | 0.6 | -0.8 | -0.2 |
| africa | 0.3 | 1.4 | 0.69 | -0.7 | 1.5 | 0.3 | -1.4 | -0.5 |
| asia | 1.5 | 1.15 | 0.18 | 0 | 0.9 | 0.7 | -0.5 | 0 |
| europe | 3.5 | 0.8 | -0.105 | 0.5 | 0.4 | 1.2 | 0 | 0.6 |
| oceania | 2.5 | 0.85 | 0 | 0.3 | 0.5 | 0.8 | 0 | 0.4 |
| home | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 0 |

#### Region–Bacteria Acquisition Log-Odds

| Region | Bacteria | Acquisition log-odds |
| --- | ---: | ---: |
| south_america | acinetobacter_baumannii | 1.6 |
| south_america | citrobacter_spp. | -0.3 |
| south_america | enterobacter_spp. | -0.2 |
| south_america | enterococcus_faecalis | 1 |
| south_america | enterococcus_faecium | 0.9 |
| south_america | escherichia_coli | 0.6 |
| south_america | klebsiella_pneumoniae | 0.8 |
| south_america | morganella_spp. | -0.5 |
| south_america | proteus_spp. | 0.5 |
| south_america | serratia_spp. | 0.4 |
| south_america | pseudomonas_aeruginosa | 1 |
| south_america | staphylococcus_aureus | 0.8 |
| south_america | streptococcus_pneumoniae | 1.2 |
| south_america | salmonella_enterica_serovar_typhi | 1.6 |
| south_america | salmonella_enterica_serovar_paratyphi_a | 1.5 |
| south_america | invasive_non-typhoidal_salmonella_spp. | 1.7 |
| south_america | neisseria_gonorrhoeae | 1.2 |
| south_america | streptococcus_pyogenes | 1 |
| south_america | streptococcus_agalactiae | 0.7 |
| south_america | haemophilus_influenzae | 1 |
| south_america | chlamydia_trachomatis | 1.1 |
| south_america | vibrio_cholerae | 1.5 |
| south_america | neisseria_meningitidis | -0.2 |
| south_america | listeria_monocytogenes | 0.1 |
| south_america | campylobacter_jejuni | 1.5 |
| south_america | enterobacter_cloacae | -0.2 |
| south_america | moraxella_catarrhalis | 0.2 |
| south_america | treponema_pallidum | 0.1 |
| south_america | helicobacter_pylori | 1.9 |
| africa | acinetobacter_baumannii | 2.3 |
| africa | citrobacter_spp. | 0.1 |
| africa | enterobacter_spp. | 0.2 |
| africa | enterococcus_faecalis | 1.6 |
| africa | enterococcus_faecium | 1.4 |
| africa | escherichia_coli | 1.2 |
| africa | klebsiella_pneumoniae | 1.6 |
| africa | morganella_spp. | -0.3 |
| africa | proteus_spp. | 1.1 |
| africa | serratia_spp. | 0.8 |
| africa | pseudomonas_aeruginosa | 1.6 |
| africa | staphylococcus_aureus | 1.5 |
| africa | streptococcus_pneumoniae | 2.2 |
| africa | salmonella_enterica_serovar_typhi | 3.3 |
| africa | salmonella_enterica_serovar_paratyphi_a | 3.2 |
| africa | invasive_non-typhoidal_salmonella_spp. | 4.5 |
| africa | neisseria_gonorrhoeae | 2.1 |
| africa | streptococcus_pyogenes | 1.8 |
| africa | streptococcus_agalactiae | 1.4 |
| africa | haemophilus_influenzae | 2 |
| africa | chlamydia_trachomatis | 2.1 |
| africa | vibrio_cholerae | 3.5 |
| africa | neisseria_meningitidis | 1.5 |
| africa | listeria_monocytogenes | 0.4 |
| africa | clostridioides_difficile | 0.2 |
| africa | campylobacter_jejuni | 2.5 |
| africa | enterobacter_cloacae | 0.2 |
| africa | yersinia_enterocolitica | 0.2 |
| africa | moraxella_catarrhalis | 0.9 |
| africa | treponema_pallidum | 0.3 |
| africa | helicobacter_pylori | 3.2 |
| asia | acinetobacter_baumannii | 2 |
| asia | citrobacter_spp. | -0.2 |
| asia | enterococcus_faecalis | 1.3 |
| asia | enterococcus_faecium | 1.2 |
| asia | escherichia_coli | 1 |
| asia | klebsiella_pneumoniae | 1.3 |
| asia | morganella_spp. | -0.4 |
| asia | proteus_spp. | 0.8 |
| asia | serratia_spp. | 0.6 |
| asia | pseudomonas_aeruginosa | 1.3 |
| asia | staphylococcus_aureus | 1.2 |
| asia | streptococcus_pneumoniae | 1.8 |
| asia | salmonella_enterica_serovar_typhi | 3 |
| asia | salmonella_enterica_serovar_paratyphi_a | 2.9 |
| asia | invasive_non-typhoidal_salmonella_spp. | 2.3 |
| asia | neisseria_gonorrhoeae | 1.6 |
| asia | streptococcus_pyogenes | 1.4 |
| asia | streptococcus_agalactiae | 1.1 |
| asia | haemophilus_influenzae | 1.6 |
| asia | chlamydia_trachomatis | 1.7 |
| asia | vibrio_cholerae | 2.8 |
| asia | neisseria_meningitidis | -0.2 |
| asia | listeria_monocytogenes | 0.3 |
| asia | clostridioides_difficile | 0.1 |
| asia | campylobacter_jejuni | 1.9 |
| asia | moraxella_catarrhalis | 0.5 |
| asia | helicobacter_pylori | 2.8 |
| europe | acinetobacter_baumannii | -0.3 |
| europe | citrobacter_spp. | -0.5 |
| europe | enterobacter_spp. | -0.5 |
| europe | enterococcus_faecalis | 0.2 |
| europe | enterococcus_faecium | 0.2 |
| europe | escherichia_coli | -0.3 |
| europe | klebsiella_pneumoniae | -0.1 |
| europe | morganella_spp. | -0.7 |
| europe | proteus_spp. | -0.1 |
| europe | serratia_spp. | -0.1 |
| europe | pseudomonas_aeruginosa | 0.3 |
| europe | staphylococcus_aureus | -0.2 |
| europe | streptococcus_pneumoniae | -0.1 |
| europe | salmonella_enterica_serovar_typhi | -1.8 |
| europe | salmonella_enterica_serovar_paratyphi_a | -2.1 |
| europe | invasive_non-typhoidal_salmonella_spp. | -0.8 |
| europe | neisseria_gonorrhoeae | 0.4 |
| europe | streptococcus_pyogenes | -0.1 |
| europe | haemophilus_influenzae | -0.3 |
| europe | chlamydia_trachomatis | 0.2 |
| europe | vibrio_cholerae | -3 |
| europe | neisseria_meningitidis | -0.8 |
| europe | listeria_monocytogenes | -0.2 |
| europe | clostridioides_difficile | -0.1 |
| europe | campylobacter_jejuni | 1 |
| europe | enterobacter_cloacae | -0.5 |
| europe | moraxella_catarrhalis | -0.3 |
| europe | treponema_pallidum | -0.2 |
| europe | helicobacter_pylori | -0.2 |
| oceania | acinetobacter_baumannii | -0.1 |
| oceania | citrobacter_spp. | -0.4 |
| oceania | enterobacter_spp. | -0.3 |
| oceania | enterococcus_faecalis | 0.2 |
| oceania | enterococcus_faecium | 0.2 |
| oceania | escherichia_coli | 0.1 |
| oceania | morganella_spp. | -0.6 |
| oceania | pseudomonas_aeruginosa | 0.4 |
| oceania | streptococcus_pneumoniae | 0.3 |
| oceania | salmonella_enterica_serovar_typhi | -1.1 |
| oceania | salmonella_enterica_serovar_paratyphi_a | -1.3 |
| oceania | invasive_non-typhoidal_salmonella_spp. | -0.3 |
| oceania | neisseria_gonorrhoeae | 0.5 |
| oceania | streptococcus_pyogenes | 0.2 |
| oceania | streptococcus_agalactiae | 0.1 |
| oceania | haemophilus_influenzae | -0.2 |
| oceania | chlamydia_trachomatis | 0.6 |
| oceania | vibrio_cholerae | -2.2 |
| oceania | neisseria_meningitidis | -0.6 |
| oceania | listeria_monocytogenes | -0.3 |
| oceania | clostridioides_difficile | -0.2 |
| oceania | campylobacter_jejuni | 0.6 |
| oceania | enterobacter_cloacae | -0.3 |
| oceania | yersinia_enterocolitica | -0.1 |
| oceania | moraxella_catarrhalis | -0.2 |
| oceania | treponema_pallidum | -0.3 |
| oceania | helicobacter_pylori | 0.4 |

### B.6 Age-Dependent Parameters

Log-odds adjustments by age category for bacteria acquisition and regional effects. Age categories: infant, preschool, school, young_adult, middle_age, elderly.

See: [§2.2 Ageing and age categories](#22-ageing-and-age-categories), [§3.1 Community acquisition](#31-community-acquisition).

#### Default Age Log-Odds

| Age category | Default log-odds |
| --- | ---: |
| infant | 1.5 |
| preschool | 0.8 |
| school | 0.3 |
| young_adult | 0 |
| middle_age | 0.2 |
| elderly | 0.9 |

#### Bacteria–Age Log-Odds

| Bacteria | Age category | Log-odds |
| --- | ---: | ---: |
| acinetobacter_baumannii | infant | 0.5 |
| acinetobacter_baumannii | preschool | -0.5 |
| acinetobacter_baumannii | school | -0.8 |
| acinetobacter_baumannii | young_adult | 0.2 |
| acinetobacter_baumannii | middle_age | 0.8 |
| acinetobacter_baumannii | elderly | 1.5 |
| citrobacter_spp. | infant | 1.5 |
| citrobacter_spp. | preschool | 0.8 |
| citrobacter_spp. | school | 0.3 |
| citrobacter_spp. | middle_age | 0.2 |
| citrobacter_spp. | elderly | 0.9 |
| enterobacter_spp. | infant | 1.5 |
| enterobacter_spp. | preschool | 0.8 |
| enterobacter_spp. | school | 0.3 |
| enterobacter_spp. | middle_age | 0.2 |
| enterobacter_spp. | elderly | 0.9 |
| enterococcus_faecalis | infant | 1.5 |
| enterococcus_faecalis | preschool | 0.8 |
| enterococcus_faecalis | school | 0.3 |
| enterococcus_faecalis | middle_age | 0.2 |
| enterococcus_faecalis | elderly | 0.9 |
| enterococcus_faecium | infant | 1.5 |
| enterococcus_faecium | preschool | 0.8 |
| enterococcus_faecium | school | 0.3 |
| enterococcus_faecium | middle_age | 0.2 |
| enterococcus_faecium | elderly | 0.9 |
| escherichia_coli | infant | 0.6 |
| escherichia_coli | preschool | -0.1 |
| escherichia_coli | school | -0.2 |
| escherichia_coli | young_adult | 0.3 |
| escherichia_coli | middle_age | 0.2 |
| escherichia_coli | elderly | 0.9 |
| klebsiella_pneumoniae | infant | 1.5 |
| klebsiella_pneumoniae | preschool | 0.8 |
| klebsiella_pneumoniae | school | 0.3 |
| klebsiella_pneumoniae | middle_age | 0.2 |
| klebsiella_pneumoniae | elderly | 0.9 |
| morganella_spp. | infant | 1.5 |
| morganella_spp. | preschool | 0.8 |
| morganella_spp. | school | 0.3 |
| morganella_spp. | middle_age | 0.2 |
| morganella_spp. | elderly | 0.9 |
| proteus_spp. | infant | 1.5 |
| proteus_spp. | preschool | 0.8 |
| proteus_spp. | school | 0.3 |
| proteus_spp. | middle_age | 0.2 |
| proteus_spp. | elderly | 0.9 |
| serratia_spp. | infant | 1.5 |
| serratia_spp. | preschool | 0.8 |
| serratia_spp. | school | 0.3 |
| serratia_spp. | middle_age | 0.2 |
| serratia_spp. | elderly | 0.9 |
| p_stuartii | infant | 1.5 |
| p_stuartii | preschool | 0.8 |
| p_stuartii | school | 0.3 |
| p_stuartii | middle_age | 0.2 |
| p_stuartii | elderly | 0.9 |
| pseudomonas_aeruginosa | infant | 1.5 |
| pseudomonas_aeruginosa | preschool | 0.8 |
| pseudomonas_aeruginosa | school | 0.3 |
| pseudomonas_aeruginosa | middle_age | 0.2 |
| pseudomonas_aeruginosa | elderly | 0.9 |
| stenotrophomonas_maltophilia | infant | 1.5 |
| stenotrophomonas_maltophilia | preschool | 0.8 |
| stenotrophomonas_maltophilia | school | 0.3 |
| stenotrophomonas_maltophilia | middle_age | 0.2 |
| stenotrophomonas_maltophilia | elderly | 0.9 |
| staphylococcus_aureus | infant | 1.5 |
| staphylococcus_aureus | preschool | 0.8 |
| staphylococcus_aureus | school | 0.3 |
| staphylococcus_aureus | middle_age | 0.2 |
| staphylococcus_aureus | elderly | 0.9 |
| staphylococcus_epidermidis | infant | 1.5 |
| staphylococcus_epidermidis | preschool | 0.8 |
| staphylococcus_epidermidis | school | 0.3 |
| staphylococcus_epidermidis | middle_age | 0.2 |
| staphylococcus_epidermidis | elderly | 0.9 |
| streptococcus_pneumoniae | infant | 1.7 |
| streptococcus_pneumoniae | preschool | 0.9 |
| streptococcus_pneumoniae | school | 0.2 |
| streptococcus_pneumoniae | young_adult | -0.4 |
| streptococcus_pneumoniae | middle_age | -0.1 |
| streptococcus_pneumoniae | elderly | 1.2 |
| salmonella_enterica_serovar_typhi | infant | 1 |
| salmonella_enterica_serovar_typhi | preschool | 0.8 |
| salmonella_enterica_serovar_typhi | school | 0.5 |
| salmonella_enterica_serovar_typhi | middle_age | 0.2 |
| salmonella_enterica_serovar_typhi | elderly | 0.8 |
| salmonella_enterica_serovar_paratyphi_a | infant | 1.5 |
| salmonella_enterica_serovar_paratyphi_a | preschool | 0.8 |
| salmonella_enterica_serovar_paratyphi_a | school | 0.3 |
| salmonella_enterica_serovar_paratyphi_a | middle_age | 0.2 |
| salmonella_enterica_serovar_paratyphi_a | elderly | 0.9 |
| invasive_non-typhoidal_salmonella_spp. | infant | 1.5 |
| invasive_non-typhoidal_salmonella_spp. | preschool | 0.8 |
| invasive_non-typhoidal_salmonella_spp. | school | 0.3 |
| invasive_non-typhoidal_salmonella_spp. | middle_age | 0.2 |
| invasive_non-typhoidal_salmonella_spp. | elderly | 0.9 |
| shigella_spp. | infant | 1.5 |
| shigella_spp. | preschool | 2 |
| shigella_spp. | school | 1.2 |
| shigella_spp. | young_adult | 0.3 |
| shigella_spp. | elderly | 0.5 |
| neisseria_gonorrhoeae | infant | -2.5 |
| neisseria_gonorrhoeae | preschool | -3.5 |
| neisseria_gonorrhoeae | school | -0.8 |
| neisseria_gonorrhoeae | young_adult | 2 |
| neisseria_gonorrhoeae | middle_age | 0.9 |
| neisseria_gonorrhoeae | elderly | -0.8 |
| streptococcus_pyogenes | infant | 1.5 |
| streptococcus_pyogenes | preschool | 0.8 |
| streptococcus_pyogenes | school | 0.3 |
| streptococcus_pyogenes | middle_age | 0.2 |
| streptococcus_pyogenes | elderly | 0.9 |
| streptococcus_agalactiae | infant | 1.5 |
| streptococcus_agalactiae | preschool | 0.8 |
| streptococcus_agalactiae | school | 0.3 |
| streptococcus_agalactiae | middle_age | 0.2 |
| streptococcus_agalactiae | elderly | 0.9 |
| haemophilus_influenzae | infant | 2.5 |
| haemophilus_influenzae | preschool | 1.5 |
| haemophilus_influenzae | school | 0.8 |
| haemophilus_influenzae | young_adult | -0.5 |
| haemophilus_influenzae | middle_age | -0.2 |
| haemophilus_influenzae | elderly | 1 |
| chlamydia_trachomatis | infant | -3 |
| chlamydia_trachomatis | preschool | -3.5 |
| chlamydia_trachomatis | school | -1 |
| chlamydia_trachomatis | young_adult | 1.7 |
| chlamydia_trachomatis | middle_age | 0.8 |
| chlamydia_trachomatis | elderly | -1 |
| mycoplasma_genitalium | infant | 1.5 |
| mycoplasma_genitalium | preschool | 0.8 |
| mycoplasma_genitalium | school | 0.3 |
| mycoplasma_genitalium | middle_age | 0.2 |
| mycoplasma_genitalium | elderly | 0.9 |
| vibrio_cholerae | infant | 1.5 |
| vibrio_cholerae | preschool | 0.8 |
| vibrio_cholerae | school | 0.3 |
| vibrio_cholerae | middle_age | 0.2 |
| vibrio_cholerae | elderly | 0.9 |
| neisseria_meningitidis | infant | 1.8 |
| neisseria_meningitidis | preschool | 0.4 |
| neisseria_meningitidis | school | 1 |
| neisseria_meningitidis | young_adult | 1.3 |
| neisseria_meningitidis | middle_age | -0.2 |
| neisseria_meningitidis | elderly | 0.2 |
| listeria_monocytogenes | infant | 1.8 |
| listeria_monocytogenes | preschool | -0.5 |
| listeria_monocytogenes | school | -1 |
| listeria_monocytogenes | young_adult | 0.5 |
| listeria_monocytogenes | elderly | 1.5 |
| clostridioides_difficile | infant | -1 |
| clostridioides_difficile | preschool | -1.5 |
| clostridioides_difficile | school | -2 |
| clostridioides_difficile | young_adult | -0.5 |
| clostridioides_difficile | middle_age | 0.5 |
| clostridioides_difficile | elderly | 2 |
| bacteroides_fragilis | infant | 1.5 |
| bacteroides_fragilis | preschool | 0.8 |
| bacteroides_fragilis | school | 0.3 |
| bacteroides_fragilis | middle_age | 0.2 |
| bacteroides_fragilis | elderly | 0.9 |
| campylobacter_jejuni | infant | 1.6 |
| campylobacter_jejuni | preschool | 1.2 |
| campylobacter_jejuni | school | 0.4 |
| campylobacter_jejuni | young_adult | 0.2 |
| campylobacter_jejuni | elderly | 0.5 |
| enterobacter_cloacae | infant | 1.5 |
| enterobacter_cloacae | preschool | 0.8 |
| enterobacter_cloacae | school | 0.3 |
| enterobacter_cloacae | middle_age | 0.2 |
| enterobacter_cloacae | elderly | 0.9 |
| yersinia_enterocolitica | infant | 1.5 |
| yersinia_enterocolitica | preschool | 0.8 |
| yersinia_enterocolitica | school | 0.3 |
| yersinia_enterocolitica | middle_age | 0.2 |
| yersinia_enterocolitica | elderly | 0.9 |
| moraxella_catarrhalis | infant | 1.8 |
| moraxella_catarrhalis | preschool | 1 |
| moraxella_catarrhalis | school | 0.3 |
| moraxella_catarrhalis | young_adult | -0.8 |
| moraxella_catarrhalis | middle_age | -0.3 |
| moraxella_catarrhalis | elderly | 1.2 |
| treponema_pallidum | infant | -2.2 |
| treponema_pallidum | preschool | -4.3 |
| treponema_pallidum | school | -2.4 |
| treponema_pallidum | young_adult | 0.4 |
| treponema_pallidum | middle_age | -0.6 |
| treponema_pallidum | elderly | -2.2 |
| bordetella_pertussis | infant | 1.5 |
| bordetella_pertussis | preschool | 0.8 |
| bordetella_pertussis | school | 0.3 |
| bordetella_pertussis | middle_age | 0.2 |
| bordetella_pertussis | elderly | 0.9 |
| helicobacter_pylori | infant | -2.5 |
| helicobacter_pylori | preschool | -1.5 |
| helicobacter_pylori | school | -0.5 |
| helicobacter_pylori | young_adult | 0.8 |
| helicobacter_pylori | middle_age | 1.4 |
| helicobacter_pylori | elderly | 1.8 |
| mdr_mycobacterium_tuberculosis | infant | 1.5 |
| mdr_mycobacterium_tuberculosis | preschool | 0.8 |
| mdr_mycobacterium_tuberculosis | school | 0.3 |
| mdr_mycobacterium_tuberculosis | middle_age | 0.2 |
| mdr_mycobacterium_tuberculosis | elderly | 0.9 |
| mycoplasma_pneumoniae | infant | 1.5 |
| mycoplasma_pneumoniae | preschool | 0.8 |
| mycoplasma_pneumoniae | school | 0.3 |
| mycoplasma_pneumoniae | middle_age | 0.2 |
| mycoplasma_pneumoniae | elderly | 0.9 |
| legionella_pneumophila | infant | 1.5 |
| legionella_pneumophila | preschool | 0.8 |
| legionella_pneumophila | school | 0.3 |
| legionella_pneumophila | middle_age | 0.2 |
| legionella_pneumophila | elderly | 0.9 |
| burkholderia_cepacia_complex | infant | 1.5 |
| burkholderia_cepacia_complex | preschool | 0.8 |
| burkholderia_cepacia_complex | school | 0.3 |
| burkholderia_cepacia_complex | middle_age | 0.2 |
| burkholderia_cepacia_complex | elderly | 0.9 |

#### Region–Age Log-Odds

| Region | Age category | Log-odds |
| --- | ---: | ---: |
| south_america | infant | 1.2 |
| south_america | preschool | 0.7 |
| south_america | school | 0.3 |
| south_america | young_adult | 0.2 |
| south_america | middle_age | 0.3 |
| south_america | elderly | 0.9 |
| africa | infant | 2 |
| africa | preschool | 1.2 |
| africa | school | 0.6 |
| africa | young_adult | 0.3 |
| africa | middle_age | 0.4 |
| africa | elderly | 1.3 |
| asia | infant | 1 |
| asia | preschool | 0.5 |
| asia | school | 0.2 |
| asia | young_adult | 0.1 |
| asia | middle_age | 0.2 |
| asia | elderly | 0.8 |
| europe | infant | -0.2 |
| europe | preschool | -0.1 |
| europe | elderly | 0.2 |
| oceania | infant | 0.1 |
| oceania | elderly | 0.3 |

### B.7 Syndrome Parameters

Infection-site (syndrome) specific parameters. Syndromes are: 1 = UTI, 2 = skin/soft tissue, 3 = respiratory, 4 = bloodstream, 5 = intra-abdominal, 6 = CNS/meningitis, 7 = gastrointestinal, 8 = genital/STI, 9 = bone/joint, 10 = other.

See: [§4.1 Syndrome assignment](#41-syndrome-assignment), [§6.2 Drug selection](#62-drug-selection-choosing-which-antibiotic-to-use), [§6.4 Drug penetration by syndrome](#64-drug-penetration-by-syndrome).

#### Syndrome Empiric Drug Scores

The full empiric syndrome-score templates are source-of-truth in `src/config.rs` under `empiric_syndrome_templates`. They are intentionally not duplicated here in full, because these weights are calibration inputs and stale copied tables are misleading. The current UTI template is included as a compact reference for the commonly inspected urinary-treatment outputs.

Empiric syndrome scores may also have era-specific overrides using keys of the form `syndrome_<id>_empiric_drug_<drug>_score_before_<YYYY>`. The parser applies the earliest cutoff year that is still greater than the simulated calendar year; otherwise it uses the base score in the static syndrome template. These overrides are currently used for syndrome 8 (genital/STI) to represent historical empiric treatment of genital discharge, PID, and suspected gonorrhoea before organism confirmation or AST. The current configured examples include sulfanilamide before 1945/1965, penicillin and aminopenicillins before 1987, tetracycline and doxycycline before 1987, TMP-SMX before 1990, ciprofloxacin and ofloxacin before 2007, ceftriaxone and cefixime before 2007, and azithromycin before 2020. Drug introduction and regional availability gates still apply.

| Syndrome | Drug | Current empiric score |
| --- | --- | ---: |
| uti | meropenem | 90.0 |
| uti | imipenem_c | 90.0 |
| uti | ertapenem | 90.0 |
| uti | meropenem_vaborbactam | 60.0 |
| uti | amoxicillin_clavulanate | 30.0 |
| uti | piperacillin_tazobactam | 30.0 |
| uti | ceftazidime_avibactam | 30.0 |
| uti | nitrofurantoin | 13.0 |
| uti | fosfomycin | 10.0 |
| uti | cephalexin | 8.0 |
| uti | cefazolin | 7.0 |
| uti | cefuroxime | 7.0 |
| uti | cefixime | 7.0 |
| uti | ciprofloxacin | 6.0 |
| uti | levofloxacin | 6.0 |
| uti | cefepime | 4.0 |
| uti | ceftazidime | 4.0 |
| uti | ceftriaxone | 3.0 |
| uti | amoxicillin | 2.0 |
| uti | ampicillin | 1.0 |
| uti | trim_sulf | 1.0 |
| uti | aztreonam_avibactam | 1.0 |
| uti | colistin | 0.2 |
| uti | vancomycin | 0.1 |
| uti | linezolid | 0.1 |
#### Syndrome Drug Penetration

| Syndrome | Drug | Penetration factor |
| --- | ---: | ---: |
| uti | penicillin_g | 0.8 |
| uti | ampicillin | 0.8 |
| uti | amoxicillin | 0.8 |
| uti | piperacillin | 0.8 |
| uti | ticarcillin | 0.8 |
| uti | cephalexin | 0.85 |
| uti | cefazolin | 0.85 |
| uti | cefuroxime | 0.85 |
| uti | ceftriaxone | 0.85 |
| uti | ceftazidime | 0.85 |
| uti | cefepime | 0.85 |
| uti | ceftaroline | 0.85 |
| uti | ceftolozane_tazobactam | 0.85 |
| uti | cefiderocol | 0.85 |
| uti | meropenem | 0.85 |
| uti | imipenem_c | 0.85 |
| uti | ertapenem | 0.85 |
| uti | aztreonam | 0.8 |
| uti | erythromycin | 0.4 |
| uti | azithromycin | 0.4 |
| uti | clarithromycin | 0.4 |
| uti | clindamycin | 0.3 |
| uti | gentamicin | 0.75 |
| uti | tobramycin | 0.75 |
| uti | amikacin | 0.75 |
| uti | tetracycline | 0.5 |
| uti | doxycycline | 0.5 |
| uti | minocycline | 0.5 |
| uti | tigecycline | 0.5 |
| uti | vancomycin | 0.6 |
| uti | teicoplanin | 0.6 |
| uti | dalbavancin | 0.6 |
| uti | linezolid | 0.7 |
| uti | tedizolid | 0.7 |
| uti | daptomycin | 0.1 |
| uti | chloramphenicol | 0.4 |
| uti | metronidazole | 0.5 |
| uti | fidaxomicin | 0 |
| uti | rifampicin | 0.4 |
| uti | amoxicillin_clavulanate | 0.8 |
| uti | piperacillin_tazobactam | 0.8 |
| uti | ampicillin_sulbactam | 0.8 |
| uti | ticarcillin_clavulanate | 0.8 |
| uti | ceftazidime_avibactam | 0.8 |
| uti | meropenem_vaborbactam | 0.8 |
| uti | colistin | 0.7 |
| uti | flucloxacillin | 0.8 |
| uti | aztreonam_avibactam | 0.8 |
| uti | cefixime | 0.9 |
| skin_soft_tissue | penicillin_g | 0.85 |
| skin_soft_tissue | ampicillin | 0.85 |
| skin_soft_tissue | amoxicillin | 0.85 |
| skin_soft_tissue | piperacillin | 0.85 |
| skin_soft_tissue | ticarcillin | 0.85 |
| skin_soft_tissue | cephalexin | 0.8 |
| skin_soft_tissue | cefazolin | 0.85 |
| skin_soft_tissue | cefuroxime | 0.85 |
| skin_soft_tissue | ceftriaxone | 0.85 |
| skin_soft_tissue | ceftazidime | 0.85 |
| skin_soft_tissue | cefepime | 0.85 |
| skin_soft_tissue | ceftaroline | 0.85 |
| skin_soft_tissue | ceftolozane_tazobactam | 0.85 |
| skin_soft_tissue | cefiderocol | 0.85 |
| skin_soft_tissue | meropenem | 0.85 |
| skin_soft_tissue | imipenem_c | 0.85 |
| skin_soft_tissue | ertapenem | 0.85 |
| skin_soft_tissue | aztreonam | 0.75 |
| skin_soft_tissue | erythromycin | 0.8 |
| skin_soft_tissue | azithromycin | 0.8 |
| skin_soft_tissue | clarithromycin | 0.8 |
| skin_soft_tissue | clindamycin | 0.85 |
| skin_soft_tissue | gentamicin | 0.6 |
| skin_soft_tissue | tobramycin | 0.6 |
| skin_soft_tissue | amikacin | 0.6 |
| skin_soft_tissue | ciprofloxacin | 0.9 |
| skin_soft_tissue | levofloxacin | 0.9 |
| skin_soft_tissue | moxifloxacin | 0.9 |
| skin_soft_tissue | ofloxacin | 0.9 |
| skin_soft_tissue | tetracycline | 0.8 |
| skin_soft_tissue | doxycycline | 0.8 |
| skin_soft_tissue | minocycline | 0.8 |
| skin_soft_tissue | tigecycline | 0.8 |
| skin_soft_tissue | vancomycin | 0.75 |
| skin_soft_tissue | teicoplanin | 0.75 |
| skin_soft_tissue | dalbavancin | 0.75 |
| skin_soft_tissue | linezolid | 0.9 |
| skin_soft_tissue | tedizolid | 0.9 |
| skin_soft_tissue | daptomycin | 0.95 |
| skin_soft_tissue | trim_sulf | 0.8 |
| skin_soft_tissue | chloramphenicol | 0.7 |
| skin_soft_tissue | nitrofurantoin | 0.2 |
| skin_soft_tissue | fosfomycin | 0.5 |
| skin_soft_tissue | metronidazole | 0.75 |
| skin_soft_tissue | fidaxomicin | 0 |
| skin_soft_tissue | rifampicin | 0.8 |
| skin_soft_tissue | amoxicillin_clavulanate | 0.85 |
| skin_soft_tissue | piperacillin_tazobactam | 0.85 |
| skin_soft_tissue | ampicillin_sulbactam | 0.85 |
| skin_soft_tissue | ticarcillin_clavulanate | 0.85 |
| skin_soft_tissue | ceftazidime_avibactam | 0.85 |
| skin_soft_tissue | meropenem_vaborbactam | 0.85 |
| skin_soft_tissue | colistin | 0.5 |
| skin_soft_tissue | flucloxacillin | 0.85 |
| skin_soft_tissue | aztreonam_avibactam | 0.75 |
| skin_soft_tissue | cefixime | 0.75 |
| respiratory | penicillin_g | 0.65 |
| respiratory | ampicillin | 0.65 |
| respiratory | amoxicillin | 0.65 |
| respiratory | piperacillin | 0.65 |
| respiratory | ticarcillin | 0.65 |
| respiratory | cephalexin | 0.55 |
| respiratory | cefazolin | 0.7 |
| respiratory | cefuroxime | 0.7 |
| respiratory | ceftriaxone | 0.7 |
| respiratory | ceftazidime | 0.7 |
| respiratory | cefepime | 0.7 |
| respiratory | ceftaroline | 0.7 |
| respiratory | ceftolozane_tazobactam | 0.7 |
| respiratory | cefiderocol | 0.7 |
| respiratory | meropenem | 0.75 |
| respiratory | imipenem_c | 0.75 |
| respiratory | ertapenem | 0.75 |
| respiratory | aztreonam | 0.6 |
| respiratory | erythromycin | 0.95 |
| respiratory | azithromycin | 0.95 |
| respiratory | clarithromycin | 0.95 |
| respiratory | clindamycin | 0.75 |
| respiratory | gentamicin | 0.4 |
| respiratory | tobramycin | 0.4 |
| respiratory | amikacin | 0.4 |
| respiratory | ciprofloxacin | 0.95 |
| respiratory | levofloxacin | 0.95 |
| respiratory | moxifloxacin | 0.95 |
| respiratory | ofloxacin | 0.95 |
| respiratory | tetracycline | 0.7 |
| respiratory | doxycycline | 0.7 |
| respiratory | minocycline | 0.7 |
| respiratory | tigecycline | 0.7 |
| respiratory | vancomycin | 0.5 |
| respiratory | teicoplanin | 0.5 |
| respiratory | dalbavancin | 0.5 |
| respiratory | linezolid | 0.9 |
| respiratory | tedizolid | 0.9 |
| respiratory | daptomycin | 0 |
| respiratory | trim_sulf | 0.8 |
| respiratory | chloramphenicol | 0.7 |
| respiratory | nitrofurantoin | 0.15 |
| respiratory | fosfomycin | 0.4 |
| respiratory | metronidazole | 0.6 |
| respiratory | fidaxomicin | 0 |
| respiratory | rifampicin | 0.85 |
| respiratory | amoxicillin_clavulanate | 0.65 |
| respiratory | piperacillin_tazobactam | 0.65 |
| respiratory | ampicillin_sulbactam | 0.65 |
| respiratory | ticarcillin_clavulanate | 0.65 |
| respiratory | ceftazidime_avibactam | 0.65 |
| respiratory | meropenem_vaborbactam | 0.65 |
| respiratory | colistin | 0.3 |
| respiratory | flucloxacillin | 0.65 |
| respiratory | aztreonam_avibactam | 0.6 |
| respiratory | cefixime | 0.6 |
| intra_abdominal | penicillin_g | 0.6 |
| intra_abdominal | ampicillin | 0.6 |
| intra_abdominal | amoxicillin | 0.6 |
| intra_abdominal | piperacillin | 0.6 |
| intra_abdominal | ticarcillin | 0.6 |
| intra_abdominal | cephalexin | 0.45 |
| intra_abdominal | cefazolin | 0.65 |
| intra_abdominal | cefuroxime | 0.65 |
| intra_abdominal | ceftriaxone | 0.65 |
| intra_abdominal | ceftazidime | 0.65 |
| intra_abdominal | cefepime | 0.65 |
| intra_abdominal | ceftaroline | 0.65 |
| intra_abdominal | ceftolozane_tazobactam | 0.65 |
| intra_abdominal | cefiderocol | 0.65 |
| intra_abdominal | meropenem | 0.75 |
| intra_abdominal | imipenem_c | 0.75 |
| intra_abdominal | ertapenem | 0.75 |
| intra_abdominal | aztreonam | 0.55 |
| intra_abdominal | erythromycin | 0.5 |
| intra_abdominal | azithromycin | 0.5 |
| intra_abdominal | clarithromycin | 0.5 |
| intra_abdominal | clindamycin | 0.65 |
| intra_abdominal | gentamicin | 0.3 |
| intra_abdominal | tobramycin | 0.3 |
| intra_abdominal | amikacin | 0.3 |
| intra_abdominal | ciprofloxacin | 0.75 |
| intra_abdominal | levofloxacin | 0.75 |
| intra_abdominal | moxifloxacin | 0.75 |
| intra_abdominal | ofloxacin | 0.75 |
| intra_abdominal | tetracycline | 0.55 |
| intra_abdominal | doxycycline | 0.55 |
| intra_abdominal | minocycline | 0.55 |
| intra_abdominal | tigecycline | 0.55 |
| intra_abdominal | vancomycin | 0.45 |
| intra_abdominal | teicoplanin | 0.45 |
| intra_abdominal | dalbavancin | 0.45 |
| intra_abdominal | linezolid | 0.7 |
| intra_abdominal | tedizolid | 0.7 |
| intra_abdominal | daptomycin | 0.6 |
| intra_abdominal | trim_sulf | 0.6 |
| intra_abdominal | chloramphenicol | 0.6 |
| intra_abdominal | nitrofurantoin | 0.15 |
| intra_abdominal | fosfomycin | 0.5 |
| intra_abdominal | metronidazole | 0.9 |
| intra_abdominal | fidaxomicin | 0.05 |
| intra_abdominal | rifampicin | 0.65 |
| intra_abdominal | amoxicillin_clavulanate | 0.65 |
| intra_abdominal | piperacillin_tazobactam | 0.65 |
| intra_abdominal | ampicillin_sulbactam | 0.65 |
| intra_abdominal | ticarcillin_clavulanate | 0.65 |
| intra_abdominal | ceftazidime_avibactam | 0.65 |
| intra_abdominal | meropenem_vaborbactam | 0.65 |
| intra_abdominal | colistin | 0.35 |
| intra_abdominal | flucloxacillin | 0.6 |
| intra_abdominal | aztreonam_avibactam | 0.55 |
| intra_abdominal | cefixime | 0.55 |
| cns_meningitis | penicillin_g | 0.15 |
| cns_meningitis | ampicillin | 0.15 |
| cns_meningitis | amoxicillin | 0.15 |
| cns_meningitis | piperacillin | 0.15 |
| cns_meningitis | ticarcillin | 0.15 |
| cns_meningitis | cephalexin | 0.05 |
| cns_meningitis | cefazolin | 0.2 |
| cns_meningitis | cefuroxime | 0.2 |
| cns_meningitis | ceftriaxone | 0.35 |
| cns_meningitis | ceftazidime | 0.2 |
| cns_meningitis | cefepime | 0.2 |
| cns_meningitis | ceftaroline | 0.2 |
| cns_meningitis | ceftolozane_tazobactam | 0.2 |
| cns_meningitis | cefiderocol | 0.2 |
| cns_meningitis | meropenem | 0.35 |
| cns_meningitis | imipenem_c | 0.25 |
| cns_meningitis | ertapenem | 0.25 |
| cns_meningitis | aztreonam | 0.1 |
| cns_meningitis | erythromycin | 0.15 |
| cns_meningitis | azithromycin | 0.15 |
| cns_meningitis | clarithromycin | 0.15 |
| cns_meningitis | clindamycin | 0.15 |
| cns_meningitis | gentamicin | 0.05 |
| cns_meningitis | tobramycin | 0.05 |
| cns_meningitis | amikacin | 0.05 |
| cns_meningitis | ciprofloxacin | 0.5 |
| cns_meningitis | levofloxacin | 0.5 |
| cns_meningitis | moxifloxacin | 0.6 |
| cns_meningitis | ofloxacin | 0.5 |
| cns_meningitis | tetracycline | 0.25 |
| cns_meningitis | doxycycline | 0.25 |
| cns_meningitis | minocycline | 0.4 |
| cns_meningitis | tigecycline | 0.25 |
| cns_meningitis | vancomycin | 0.15 |
| cns_meningitis | teicoplanin | 0.15 |
| cns_meningitis | dalbavancin | 0.15 |
| cns_meningitis | linezolid | 0.7 |
| cns_meningitis | tedizolid | 0.7 |
| cns_meningitis | daptomycin | 0.05 |
| cns_meningitis | trim_sulf | 0.5 |
| cns_meningitis | chloramphenicol | 0.7 |
| cns_meningitis | nitrofurantoin | 0.05 |
| cns_meningitis | fosfomycin | 0.3 |
| cns_meningitis | metronidazole | 0.8 |
| cns_meningitis | fidaxomicin | 0 |
| cns_meningitis | rifampicin | 0.5 |
| cns_meningitis | amoxicillin_clavulanate | 0.15 |
| cns_meningitis | piperacillin_tazobactam | 0.15 |
| cns_meningitis | ampicillin_sulbactam | 0.15 |
| cns_meningitis | ticarcillin_clavulanate | 0.15 |
| cns_meningitis | ceftazidime_avibactam | 0.15 |
| cns_meningitis | meropenem_vaborbactam | 0.15 |
| cns_meningitis | colistin | 0.05 |
| cns_meningitis | flucloxacillin | 0.15 |
| cns_meningitis | aztreonam_avibactam | 0.1 |
| cns_meningitis | cefixime | 0.1 |
| gastrointestinal | penicillin_g | 0.55 |
| gastrointestinal | ampicillin | 0.55 |
| gastrointestinal | amoxicillin | 0.55 |
| gastrointestinal | piperacillin | 0.55 |
| gastrointestinal | ticarcillin | 0.55 |
| gastrointestinal | cephalexin | 0.5 |
| gastrointestinal | cefazolin | 0.6 |
| gastrointestinal | cefuroxime | 0.6 |
| gastrointestinal | ceftriaxone | 0.6 |
| gastrointestinal | ceftazidime | 0.6 |
| gastrointestinal | cefepime | 0.6 |
| gastrointestinal | ceftaroline | 0.6 |
| gastrointestinal | ceftolozane_tazobactam | 0.6 |
| gastrointestinal | cefiderocol | 0.6 |
| gastrointestinal | meropenem | 0.65 |
| gastrointestinal | imipenem_c | 0.65 |
| gastrointestinal | ertapenem | 0.65 |
| gastrointestinal | aztreonam | 0.5 |
| gastrointestinal | erythromycin | 0.7 |
| gastrointestinal | azithromycin | 0.7 |
| gastrointestinal | clarithromycin | 0.7 |
| gastrointestinal | clindamycin | 0.65 |
| gastrointestinal | gentamicin | 0.4 |
| gastrointestinal | tobramycin | 0.4 |
| gastrointestinal | amikacin | 0.4 |
| gastrointestinal | ciprofloxacin | 0.85 |
| gastrointestinal | levofloxacin | 0.85 |
| gastrointestinal | moxifloxacin | 0.85 |
| gastrointestinal | ofloxacin | 0.85 |
| gastrointestinal | tetracycline | 0.6 |
| gastrointestinal | doxycycline | 0.6 |
| gastrointestinal | minocycline | 0.6 |
| gastrointestinal | tigecycline | 0.6 |
| gastrointestinal | vancomycin | 0.9 |
| gastrointestinal | teicoplanin | 0.35 |
| gastrointestinal | dalbavancin | 0.35 |
| gastrointestinal | linezolid | 0.75 |
| gastrointestinal | tedizolid | 0.75 |
| gastrointestinal | daptomycin | 0.3 |
| gastrointestinal | trim_sulf | 0.7 |
| gastrointestinal | chloramphenicol | 0.65 |
| gastrointestinal | nitrofurantoin | 0.25 |
| gastrointestinal | fosfomycin | 0.4 |
| gastrointestinal | metronidazole | 0.95 |
| gastrointestinal | furazolidone | 0.9 |
| gastrointestinal | rifampicin | 0.6 |
| gastrointestinal | amoxicillin_clavulanate | 0.55 |
| gastrointestinal | piperacillin_tazobactam | 0.55 |
| gastrointestinal | ampicillin_sulbactam | 0.55 |
| gastrointestinal | ticarcillin_clavulanate | 0.55 |
| gastrointestinal | ceftazidime_avibactam | 0.55 |
| gastrointestinal | meropenem_vaborbactam | 0.55 |
| gastrointestinal | colistin | 0.4 |
| gastrointestinal | flucloxacillin | 0.55 |
| gastrointestinal | aztreonam_avibactam | 0.5 |
| gastrointestinal | cefixime | 0.55 |
| genital_sti | penicillin_g | 0.55 |
| genital_sti | ampicillin | 0.55 |
| genital_sti | amoxicillin | 0.55 |
| genital_sti | piperacillin | 0.55 |
| genital_sti | ticarcillin | 0.55 |
| genital_sti | cephalexin | 0.45 |
| genital_sti | cefazolin | 0.55 |
| genital_sti | cefuroxime | 0.55 |
| genital_sti | ceftriaxone | 0.55 |
| genital_sti | ceftazidime | 0.55 |
| genital_sti | cefepime | 0.55 |
| genital_sti | ceftaroline | 0.55 |
| genital_sti | ceftolozane_tazobactam | 0.55 |
| genital_sti | cefiderocol | 0.55 |
| genital_sti | meropenem | 0.6 |
| genital_sti | imipenem_c | 0.6 |
| genital_sti | ertapenem | 0.6 |
| genital_sti | aztreonam | 0.45 |
| genital_sti | erythromycin | 0.75 |
| genital_sti | azithromycin | 0.75 |
| genital_sti | clarithromycin | 0.75 |
| genital_sti | clindamycin | 0.6 |
| genital_sti | gentamicin | 0.35 |
| genital_sti | tobramycin | 0.35 |
| genital_sti | amikacin | 0.35 |
| genital_sti | ciprofloxacin | 0.9 |
| genital_sti | levofloxacin | 0.9 |
| genital_sti | moxifloxacin | 0.9 |
| genital_sti | ofloxacin | 0.9 |
| genital_sti | tetracycline | 0.75 |
| genital_sti | doxycycline | 0.75 |
| genital_sti | minocycline | 0.75 |
| genital_sti | tigecycline | 0.75 |
| genital_sti | vancomycin | 0.4 |
| genital_sti | teicoplanin | 0.4 |
| genital_sti | dalbavancin | 0.4 |
| genital_sti | linezolid | 0.7 |
| genital_sti | tedizolid | 0.7 |
| genital_sti | daptomycin | 0.4 |
| genital_sti | trim_sulf | 0.8 |
| genital_sti | chloramphenicol | 0.55 |
| genital_sti | nitrofurantoin | 0.3 |
| genital_sti | fosfomycin | 0.5 |
| genital_sti | metronidazole | 0.8 |
| genital_sti | fidaxomicin | 0 |
| genital_sti | rifampicin | 0.6 |
| genital_sti | amoxicillin_clavulanate | 0.55 |
| genital_sti | piperacillin_tazobactam | 0.55 |
| genital_sti | ampicillin_sulbactam | 0.55 |
| genital_sti | ticarcillin_clavulanate | 0.55 |
| genital_sti | ceftazidime_avibactam | 0.55 |
| genital_sti | meropenem_vaborbactam | 0.55 |
| genital_sti | colistin | 0.3 |
| genital_sti | flucloxacillin | 0.55 |
| genital_sti | aztreonam_avibactam | 0.45 |
| genital_sti | cefixime | 0.5 |
| bone_joint | penicillin_g | 0.4 |
| bone_joint | ampicillin | 0.4 |
| bone_joint | amoxicillin | 0.4 |
| bone_joint | piperacillin | 0.4 |
| bone_joint | ticarcillin | 0.4 |
| bone_joint | cephalexin | 0.3 |
| bone_joint | cefazolin | 0.45 |
| bone_joint | cefuroxime | 0.45 |
| bone_joint | ceftriaxone | 0.45 |
| bone_joint | ceftazidime | 0.45 |
| bone_joint | cefepime | 0.45 |
| bone_joint | ceftaroline | 0.45 |
| bone_joint | ceftolozane_tazobactam | 0.45 |
| bone_joint | cefiderocol | 0.45 |
| bone_joint | meropenem | 0.5 |
| bone_joint | imipenem_c | 0.5 |
| bone_joint | ertapenem | 0.5 |
| bone_joint | aztreonam | 0.35 |
| bone_joint | erythromycin | 0.4 |
| bone_joint | azithromycin | 0.4 |
| bone_joint | clarithromycin | 0.4 |
| bone_joint | clindamycin | 0.6 |
| bone_joint | gentamicin | 0.25 |
| bone_joint | tobramycin | 0.25 |
| bone_joint | amikacin | 0.25 |
| bone_joint | ciprofloxacin | 0.7 |
| bone_joint | levofloxacin | 0.7 |
| bone_joint | moxifloxacin | 0.7 |
| bone_joint | ofloxacin | 0.7 |
| bone_joint | tetracycline | 0.5 |
| bone_joint | doxycycline | 0.5 |
| bone_joint | minocycline | 0.5 |
| bone_joint | tigecycline | 0.5 |
| bone_joint | vancomycin | 0.35 |
| bone_joint | teicoplanin | 0.35 |
| bone_joint | dalbavancin | 0.35 |
| bone_joint | linezolid | 0.75 |
| bone_joint | tedizolid | 0.75 |
| bone_joint | daptomycin | 0.5 |
| bone_joint | trim_sulf | 0.55 |
| bone_joint | chloramphenicol | 0.5 |
| bone_joint | nitrofurantoin | 0.1 |
| bone_joint | fosfomycin | 0.6 |
| bone_joint | metronidazole | 0.55 |
| bone_joint | fidaxomicin | 0 |
| bone_joint | rifampicin | 0.8 |
| bone_joint | amoxicillin_clavulanate | 0.4 |
| bone_joint | piperacillin_tazobactam | 0.4 |
| bone_joint | ampicillin_sulbactam | 0.4 |
| bone_joint | ticarcillin_clavulanate | 0.4 |
| bone_joint | ceftazidime_avibactam | 0.4 |
| bone_joint | meropenem_vaborbactam | 0.4 |
| bone_joint | colistin | 0.2 |
| bone_joint | flucloxacillin | 0.4 |
| bone_joint | aztreonam_avibactam | 0.35 |
| bone_joint | cefixime | 0.4 |
| other | gentamicin | 0.7 |
| other | tobramycin | 0.7 |
| other | amikacin | 0.7 |
| other | nitrofurantoin | 0.3 |

### B.8 Clearance Parameters

Infection clearance model parameters. The clearance hazard is a logistic function of base log-odds, per-bacteria adjustments, age effects, immunodeficiency, bacteria level, and treatment duration.

See: [§4.4 Natural clearance](#44-natural-clearance).

| Parameter | Value |
| --- | ---: |
| base_clearance_log_odds | -4.2 |
| immunodeficient_log_odds_adjustment | -0.69 |

#### Per-Bacteria Clearance Adjustments

| Bacteria | Log-odds adjustment |
| --- | ---: |

### B.9 Immunodeficiency, Sex, and Vaccination Parameters

See: [§2.3 Immunodeficiency](#23-immunodeficiency), [§10 Mortality](#10-mortality).

#### Immunodeficiency

| Parameter | Value |
| --- | ---: |
| startup_seed_fraction | 0.05 |
| temporary_onset_rate_per_day | 5e-5 |
| temporary_recovery_rate_per_day | 0.01 |
| chronic_onset_rate_per_day | 6e-5 |
| chronic_recovery_rate_per_day | 0.0012 |
| chronic_probability_age_0_1 | 0.3 |
| chronic_probability_age_1_18 | 0.2 |
| chronic_probability_age_18_65 | 0.4 |
| chronic_probability_age_65_plus | 0.6 |

#### Sex

| Sex | Mortality log-odds |
| --- | ---: |
| male | 0.095 |
| female | -0.105 |

#### Vaccination

Vaccination parameters are split into three parts:

- a vaccine-specific historical availability year,
- a target birth-cohort coverage reached over a configurable rollout period,
- and a bacterium-specific acquisition-effect term `log_odds_vaccinated` (default `-2.0`).

`vaccination_status` is stored per bacterium rather than per vaccine brand, and once acquired it is permanent within the current model. The active runtime mapping is pneumococcal → *S. pneumoniae*, meningococcal → *N. meningitidis*, Hib → *H. influenzae*, and pertussis → *B. pertussis*. Vaccination is assigned once at birth / first day alive, not as a repeated daily age-band hazard.

Under the default parameter map below, vaccination is active. Coverage ramps linearly from 0 at the availability year to the target birth-cohort coverage over `rollout_years`.

| Vaccine | Availability year | Target birth-cohort coverage | Rollout years |
| --- | ---: | ---: | ---: |
| pneumococcal | 2000 | 0.75 | 20 |
| meningococcal | 1981 | 0.55 | 20 |
| hib | 1985 | 0.85 | 15 |
| pertussis | 1948 | 0.82 | 20 |

### B.10 Resistance Mechanisms

Parameters for the 46 resistance mechanisms modelled. Each mechanism has a per-day reversion rate, per-drug-class enhancement multipliers, and per-bacteria emergence rates.

See: [§7.1 Resistance mechanisms](#71-resistance-mechanisms), [§7.2 Mechanism–drug-class enhancement](#72-mechanismdrug-class-enhancement-multipliers), [§7.3 Resistance emergence](#73-resistance-emergence), [§7.4 Resistance reversion](#74-resistance-reversion-and-fitness-costs).

#### Mechanism Reversion Rates

| Mechanism | Reversion rate/day |
| --- | ---: |
| enzyme_esbl_ctx_m | 6e-4 |
| enzyme_esbl_tem | 6e-4 |
| enzyme_esbl_shv | 6e-4 |
| enzyme_kpc | 0.001 |
| enzyme_ndm_vim | 0.0015 |
| enzyme_oxa_48 | 5e-4 |
| enzyme_ampc_cmy | 1e-4 |
| enzyme_ampc_dha | 6e-4 |
| target_site_pbp2a_meca | 9e-4 |
| target_site_van_a | 0.002 |
| target_site_van_b | 0.002 |
| mutation_gyra_primary | 1e-4 |
| mutation_gyra_parc_secondary | 2e-4 |
| protection_qnr | 1e-4 |
| enzyme_16s_rrmt | 5e-4 |
| target_site_erm_b | 0.002 |
| target_site_cfr | 5e-4 |
| enzyme_cat | 5e-4 |
| efflux_acrab_tolc | 5e-4 |
| efflux_mexxy_oprm | 5e-4 |
| porin_loss_ompk35_36 | 5e-4 |
| porin_loss_oprd | 5e-4 |
| modification_mcr_1 | 0.0015 |
| global_efflux_pump | 5e-4 |
| global_porin_loss | 5e-4 |
| mutation_folate_pathway | 1e-4 |
| mutation_nitroreductase | 3e-4 |
| enzyme_fos_a | 5e-4 |
| mutation_mpr_f | 0.001 |
| mutation_rpo_b | 0.002 |
| protection_fus_b | 5e-4 |
| protection_tet_m | 5e-4 |
| enzyme_aac_aph | 1e-4 |
| enzyme_bla_z | 1e-4 |
| enzyme_oxa_acinetobacter | 1e-4 |
| mutation_23s_rrna | 1e-4 |
| mutation_ampc_derepression | 2e-4 |
| mutation_polymyxin_regulatory | 0.0015 |
| mutation_liafsr_cls | 0.0015 |
| mutation_23s_rrna_oxazolidinone | 0.001 |
| efflux_tet_abc | 1e-4 |
| mutation_pbp_mosaic | 0.001 |
| efflux_mtr_cde | 0.001 |
| as_yet_unknown | 0.001 |

#### Mechanism Enhancement Multipliers by Drug Class

How much resistance each mechanism confers against each drug class. Only non-zero entries shown.

| Mechanism | Drug class | Enhancement multiplier |
| --- | ---: | ---: |
| enzyme_esbl_ctx_m | pen | 0.9 |
| enzyme_esbl_ctx_m | bli | 0.25 |
| enzyme_esbl_ctx_m | bli_anti_pseudomonal | 0.25 |
| enzyme_esbl_ctx_m | bli_sulbactam | 0.25 |
| enzyme_esbl_ctx_m | c1_2g | 0.9 |
| enzyme_esbl_ctx_m | c3g | 0.85 |
| enzyme_esbl_ctx_m | c3g_bli | 0.85 |
| enzyme_esbl_ctx_m | c4g | 0.35 |
| enzyme_esbl_ctx_m | anti_mrsa_ceph | 0.35 |
| enzyme_esbl_ctx_m | siderophore_ceph | 0.35 |
| enzyme_esbl_ctx_m | cft_avi | 0.1 |
| enzyme_esbl_ctx_m | mer_vab | 0.1 |
| enzyme_esbl_ctx_m | azt_avi | 0.1 |
| enzyme_esbl_ctx_m | mono | 0.8 |
| enzyme_esbl_ctx_m | fq | 0.8 |
| enzyme_esbl_ctx_m | ag_group1 | 0.8 |
| enzyme_esbl_ctx_m | ag_group2 | 0.8 |
| enzyme_esbl_ctx_m | mls | 0.8 |
| enzyme_esbl_ctx_m | lincosamides | 0.8 |
| enzyme_esbl_ctx_m | glyc | 0.8 |
| enzyme_esbl_ctx_m | lipoglycopeptides | 0.8 |
| enzyme_esbl_ctx_m | tet | 0.8 |
| enzyme_esbl_ctx_m | glycylcyclines | 0.8 |
| enzyme_esbl_ctx_m | poly | 0.8 |
| enzyme_esbl_ctx_m | oxa | 0.8 |
| enzyme_esbl_ctx_m | chl | 0.8 |
| enzyme_esbl_ctx_m | sulf | 0.8 |
| enzyme_esbl_ctx_m | lipopeptides | 0.8 |
| enzyme_esbl_ctx_m | streptogramins | 0.8 |
| enzyme_esbl_ctx_m | nitrofurans | 0.8 |
| enzyme_esbl_ctx_m | phosphonic_acids | 0.8 |
| enzyme_esbl_ctx_m | nitroimidazoles | 0.8 |
| enzyme_esbl_ctx_m | rifamycins | 0.8 |
| enzyme_esbl_ctx_m | macrocycles | 0.8 |
| enzyme_esbl_ctx_m | steroid_antibacterials | 0.8 |
| enzyme_esbl_ctx_m | pleuromutilins | 0.8 |
| enzyme_esbl_ctx_m | other | 0.8 |
| enzyme_esbl_tem | pen | 0.85 |
| enzyme_esbl_tem | bli | 0.2 |
| enzyme_esbl_tem | bli_anti_pseudomonal | 0.2 |
| enzyme_esbl_tem | bli_sulbactam | 0.2 |
| enzyme_esbl_tem | c1_2g | 0.85 |
| enzyme_esbl_tem | c3g | 0.65 |
| enzyme_esbl_tem | c3g_bli | 0.65 |
| enzyme_esbl_tem | c4g | 0.25 |
| enzyme_esbl_tem | anti_mrsa_ceph | 0.25 |
| enzyme_esbl_tem | siderophore_ceph | 0.25 |
| enzyme_esbl_tem | cft_avi | 0.1 |
| enzyme_esbl_tem | mer_vab | 0.1 |
| enzyme_esbl_tem | azt_avi | 0.1 |
| enzyme_esbl_tem | mono | 0.6 |
| enzyme_esbl_tem | fq | 0.6 |
| enzyme_esbl_tem | ag_group1 | 0.6 |
| enzyme_esbl_tem | ag_group2 | 0.6 |
| enzyme_esbl_tem | mls | 0.6 |
| enzyme_esbl_tem | lincosamides | 0.6 |
| enzyme_esbl_tem | glyc | 0.6 |
| enzyme_esbl_tem | lipoglycopeptides | 0.6 |
| enzyme_esbl_tem | tet | 0.6 |
| enzyme_esbl_tem | glycylcyclines | 0.6 |
| enzyme_esbl_tem | poly | 0.6 |
| enzyme_esbl_tem | oxa | 0.6 |
| enzyme_esbl_tem | chl | 0.6 |
| enzyme_esbl_tem | sulf | 0.6 |
| enzyme_esbl_tem | lipopeptides | 0.6 |
| enzyme_esbl_tem | streptogramins | 0.6 |
| enzyme_esbl_tem | nitrofurans | 0.6 |
| enzyme_esbl_tem | phosphonic_acids | 0.6 |
| enzyme_esbl_tem | nitroimidazoles | 0.6 |
| enzyme_esbl_tem | rifamycins | 0.6 |
| enzyme_esbl_tem | macrocycles | 0.6 |
| enzyme_esbl_tem | steroid_antibacterials | 0.6 |
| enzyme_esbl_tem | pleuromutilins | 0.6 |
| enzyme_esbl_tem | other | 0.6 |
| enzyme_esbl_shv | pen | 0.8 |
| enzyme_esbl_shv | bli | 0.2 |
| enzyme_esbl_shv | bli_anti_pseudomonal | 0.2 |
| enzyme_esbl_shv | bli_sulbactam | 0.2 |
| enzyme_esbl_shv | c1_2g | 0.85 |
| enzyme_esbl_shv | c3g | 0.65 |
| enzyme_esbl_shv | c3g_bli | 0.65 |
| enzyme_esbl_shv | c4g | 0.3 |
| enzyme_esbl_shv | anti_mrsa_ceph | 0.3 |
| enzyme_esbl_shv | siderophore_ceph | 0.3 |
| enzyme_esbl_shv | cft_avi | 0.1 |
| enzyme_esbl_shv | mer_vab | 0.1 |
| enzyme_esbl_shv | azt_avi | 0.1 |
| enzyme_esbl_shv | mono | 0.55 |
| enzyme_esbl_shv | fq | 0.6 |
| enzyme_esbl_shv | ag_group1 | 0.6 |
| enzyme_esbl_shv | ag_group2 | 0.6 |
| enzyme_esbl_shv | mls | 0.6 |
| enzyme_esbl_shv | lincosamides | 0.6 |
| enzyme_esbl_shv | glyc | 0.6 |
| enzyme_esbl_shv | lipoglycopeptides | 0.6 |
| enzyme_esbl_shv | tet | 0.6 |
| enzyme_esbl_shv | glycylcyclines | 0.6 |
| enzyme_esbl_shv | poly | 0.6 |
| enzyme_esbl_shv | oxa | 0.6 |
| enzyme_esbl_shv | chl | 0.6 |
| enzyme_esbl_shv | sulf | 0.6 |
| enzyme_esbl_shv | lipopeptides | 0.6 |
| enzyme_esbl_shv | streptogramins | 0.6 |
| enzyme_esbl_shv | nitrofurans | 0.6 |
| enzyme_esbl_shv | phosphonic_acids | 0.6 |
| enzyme_esbl_shv | nitroimidazoles | 0.6 |
| enzyme_esbl_shv | rifamycins | 0.6 |
| enzyme_esbl_shv | macrocycles | 0.6 |
| enzyme_esbl_shv | steroid_antibacterials | 0.6 |
| enzyme_esbl_shv | pleuromutilins | 0.6 |
| enzyme_esbl_shv | other | 0.6 |
| enzyme_kpc | pen | 0.95 |
| enzyme_kpc | bli | 0.85 |
| enzyme_kpc | bli_anti_pseudomonal | 0.85 |
| enzyme_kpc | bli_sulbactam | 0.85 |
| enzyme_kpc | c1_2g | 0.95 |
| enzyme_kpc | c3g | 0.95 |
| enzyme_kpc | c3g_bli | 0.95 |
| enzyme_kpc | c4g | 0.85 |
| enzyme_kpc | anti_mrsa_ceph | 0.85 |
| enzyme_kpc | siderophore_ceph | 0.85 |
| enzyme_kpc | cft_avi | 0.3 |
| enzyme_kpc | mer_vab | 0.3 |
| enzyme_kpc | azt_avi | 0.3 |
| enzyme_kpc | carb_group1 | 0.9 |
| enzyme_kpc | carb_group2 | 0.9 |
| enzyme_kpc | mono | 0.9 |
| enzyme_kpc | fq | 0.95 |
| enzyme_kpc | ag_group1 | 0.95 |
| enzyme_kpc | ag_group2 | 0.95 |
| enzyme_kpc | mls | 0.95 |
| enzyme_kpc | lincosamides | 0.95 |
| enzyme_kpc | glyc | 0.95 |
| enzyme_kpc | lipoglycopeptides | 0.95 |
| enzyme_kpc | tet | 0.95 |
| enzyme_kpc | glycylcyclines | 0.95 |
| enzyme_kpc | poly | 0.95 |
| enzyme_kpc | oxa | 0.95 |
| enzyme_kpc | chl | 0.95 |
| enzyme_kpc | sulf | 0.95 |
| enzyme_kpc | lipopeptides | 0.95 |
| enzyme_kpc | streptogramins | 0.95 |
| enzyme_kpc | nitrofurans | 0.95 |
| enzyme_kpc | phosphonic_acids | 0.95 |
| enzyme_kpc | nitroimidazoles | 0.95 |
| enzyme_kpc | rifamycins | 0.95 |
| enzyme_kpc | macrocycles | 0.95 |
| enzyme_kpc | steroid_antibacterials | 0.95 |
| enzyme_kpc | pleuromutilins | 0.95 |
| enzyme_kpc | other | 0.95 |
| enzyme_ndm_vim | pen | 0.95 |
| enzyme_ndm_vim | bli | 0.95 |
| enzyme_ndm_vim | bli_anti_pseudomonal | 0.95 |
| enzyme_ndm_vim | bli_sulbactam | 0.95 |
| enzyme_ndm_vim | c1_2g | 0.95 |
| enzyme_ndm_vim | c3g | 0.95 |
| enzyme_ndm_vim | c3g_bli | 0.95 |
| enzyme_ndm_vim | c4g | 0.9 |
| enzyme_ndm_vim | anti_mrsa_ceph | 0.9 |
| enzyme_ndm_vim | siderophore_ceph | 0.9 |
| enzyme_ndm_vim | cft_avi | 0.5 |
| enzyme_ndm_vim | mer_vab | 0.5 |
| enzyme_ndm_vim | azt_avi | 0.5 |
| enzyme_ndm_vim | carb_group1 | 0.95 |
| enzyme_ndm_vim | carb_group2 | 0.95 |
| enzyme_ndm_vim | mono | 0.1 |
| enzyme_ndm_vim | fq | 0.95 |
| enzyme_ndm_vim | ag_group1 | 0.95 |
| enzyme_ndm_vim | ag_group2 | 0.95 |
| enzyme_ndm_vim | mls | 0.95 |
| enzyme_ndm_vim | lincosamides | 0.95 |
| enzyme_ndm_vim | glyc | 0.95 |
| enzyme_ndm_vim | lipoglycopeptides | 0.95 |
| enzyme_ndm_vim | tet | 0.95 |
| enzyme_ndm_vim | glycylcyclines | 0.95 |
| enzyme_ndm_vim | poly | 0.95 |
| enzyme_ndm_vim | oxa | 0.95 |
| enzyme_ndm_vim | chl | 0.95 |
| enzyme_ndm_vim | sulf | 0.95 |
| enzyme_ndm_vim | lipopeptides | 0.95 |
| enzyme_ndm_vim | streptogramins | 0.95 |
| enzyme_ndm_vim | nitrofurans | 0.95 |
| enzyme_ndm_vim | phosphonic_acids | 0.95 |
| enzyme_ndm_vim | nitroimidazoles | 0.95 |
| enzyme_ndm_vim | rifamycins | 0.95 |
| enzyme_ndm_vim | macrocycles | 0.95 |
| enzyme_ndm_vim | steroid_antibacterials | 0.95 |
| enzyme_ndm_vim | pleuromutilins | 0.95 |
| enzyme_ndm_vim | other | 0.95 |
| enzyme_oxa_48 | pen | 0.8 |
| enzyme_oxa_48 | bli | 0.5 |
| enzyme_oxa_48 | bli_anti_pseudomonal | 0.5 |
| enzyme_oxa_48 | bli_sulbactam | 0.5 |
| enzyme_oxa_48 | c1_2g | 0.4 |
| enzyme_oxa_48 | c3g | 0.15 |
| enzyme_oxa_48 | c3g_bli | 0.15 |
| enzyme_oxa_48 | c4g | 0.1 |
| enzyme_oxa_48 | anti_mrsa_ceph | 0.1 |
| enzyme_oxa_48 | siderophore_ceph | 0.1 |
| enzyme_oxa_48 | cft_avi | 0.15 |
| enzyme_oxa_48 | mer_vab | 0.15 |
| enzyme_oxa_48 | azt_avi | 0.15 |
| enzyme_oxa_48 | carb_group1 | 0.7 |
| enzyme_oxa_48 | carb_group2 | 0.7 |
| enzyme_oxa_48 | fq | 0.6 |
| enzyme_oxa_48 | ag_group1 | 0.6 |
| enzyme_oxa_48 | ag_group2 | 0.6 |
| enzyme_oxa_48 | mls | 0.6 |
| enzyme_oxa_48 | lincosamides | 0.6 |
| enzyme_oxa_48 | glyc | 0.6 |
| enzyme_oxa_48 | lipoglycopeptides | 0.6 |
| enzyme_oxa_48 | tet | 0.6 |
| enzyme_oxa_48 | glycylcyclines | 0.6 |
| enzyme_oxa_48 | poly | 0.6 |
| enzyme_oxa_48 | oxa | 0.6 |
| enzyme_oxa_48 | chl | 0.6 |
| enzyme_oxa_48 | sulf | 0.6 |
| enzyme_oxa_48 | lipopeptides | 0.6 |
| enzyme_oxa_48 | streptogramins | 0.6 |
| enzyme_oxa_48 | nitrofurans | 0.6 |
| enzyme_oxa_48 | phosphonic_acids | 0.6 |
| enzyme_oxa_48 | nitroimidazoles | 0.6 |
| enzyme_oxa_48 | rifamycins | 0.6 |
| enzyme_oxa_48 | macrocycles | 0.6 |
| enzyme_oxa_48 | steroid_antibacterials | 0.6 |
| enzyme_oxa_48 | pleuromutilins | 0.6 |
| enzyme_oxa_48 | other | 0.6 |
| enzyme_ampc_cmy | pen | 0.7 |
| enzyme_ampc_cmy | bli | 0.6 |
| enzyme_ampc_cmy | bli_anti_pseudomonal | 0.6 |
| enzyme_ampc_cmy | bli_sulbactam | 0.6 |
| enzyme_ampc_cmy | c1_2g | 0.8 |
| enzyme_ampc_cmy | c3g | 0.8 |
| enzyme_ampc_cmy | c3g_bli | 0.8 |
| enzyme_ampc_cmy | c4g | 0.15 |
| enzyme_ampc_cmy | anti_mrsa_ceph | 0.15 |
| enzyme_ampc_cmy | siderophore_ceph | 0.15 |
| enzyme_ampc_cmy | cft_avi | 0.1 |
| enzyme_ampc_cmy | mer_vab | 0.1 |
| enzyme_ampc_cmy | azt_avi | 0.1 |
| enzyme_ampc_cmy | mono | 0.1 |
| enzyme_ampc_cmy | fq | 0.7 |
| enzyme_ampc_cmy | ag_group1 | 0.7 |
| enzyme_ampc_cmy | ag_group2 | 0.7 |
| enzyme_ampc_cmy | mls | 0.7 |
| enzyme_ampc_cmy | lincosamides | 0.7 |
| enzyme_ampc_cmy | glyc | 0.7 |
| enzyme_ampc_cmy | lipoglycopeptides | 0.7 |
| enzyme_ampc_cmy | tet | 0.7 |
| enzyme_ampc_cmy | glycylcyclines | 0.7 |
| enzyme_ampc_cmy | poly | 0.7 |
| enzyme_ampc_cmy | oxa | 0.7 |
| enzyme_ampc_cmy | chl | 0.7 |
| enzyme_ampc_cmy | sulf | 0.7 |
| enzyme_ampc_cmy | lipopeptides | 0.7 |
| enzyme_ampc_cmy | streptogramins | 0.7 |
| enzyme_ampc_cmy | nitrofurans | 0.7 |
| enzyme_ampc_cmy | phosphonic_acids | 0.7 |
| enzyme_ampc_cmy | nitroimidazoles | 0.7 |
| enzyme_ampc_cmy | rifamycins | 0.7 |
| enzyme_ampc_cmy | macrocycles | 0.7 |
| enzyme_ampc_cmy | steroid_antibacterials | 0.7 |
| enzyme_ampc_cmy | pleuromutilins | 0.7 |
| enzyme_ampc_cmy | other | 0.7 |
| enzyme_ampc_dha | pen | 0.7 |
| enzyme_ampc_dha | bli | 0.55 |
| enzyme_ampc_dha | bli_anti_pseudomonal | 0.55 |
| enzyme_ampc_dha | bli_sulbactam | 0.55 |
| enzyme_ampc_dha | c1_2g | 0.75 |
| enzyme_ampc_dha | c3g | 0.75 |
| enzyme_ampc_dha | c3g_bli | 0.75 |
| enzyme_ampc_dha | c4g | 0.15 |
| enzyme_ampc_dha | anti_mrsa_ceph | 0.15 |
| enzyme_ampc_dha | siderophore_ceph | 0.15 |
| enzyme_ampc_dha | cft_avi | 0.1 |
| enzyme_ampc_dha | mer_vab | 0.1 |
| enzyme_ampc_dha | azt_avi | 0.1 |
| enzyme_ampc_dha | mono | 0.1 |
| enzyme_ampc_dha | fq | 0.7 |
| enzyme_ampc_dha | ag_group1 | 0.7 |
| enzyme_ampc_dha | ag_group2 | 0.7 |
| enzyme_ampc_dha | mls | 0.7 |
| enzyme_ampc_dha | lincosamides | 0.7 |
| enzyme_ampc_dha | glyc | 0.7 |
| enzyme_ampc_dha | lipoglycopeptides | 0.7 |
| enzyme_ampc_dha | tet | 0.7 |
| enzyme_ampc_dha | glycylcyclines | 0.7 |
| enzyme_ampc_dha | poly | 0.7 |
| enzyme_ampc_dha | oxa | 0.7 |
| enzyme_ampc_dha | chl | 0.7 |
| enzyme_ampc_dha | sulf | 0.7 |
| enzyme_ampc_dha | lipopeptides | 0.7 |
| enzyme_ampc_dha | streptogramins | 0.7 |
| enzyme_ampc_dha | nitrofurans | 0.7 |
| enzyme_ampc_dha | phosphonic_acids | 0.7 |
| enzyme_ampc_dha | nitroimidazoles | 0.7 |
| enzyme_ampc_dha | rifamycins | 0.7 |
| enzyme_ampc_dha | macrocycles | 0.7 |
| enzyme_ampc_dha | steroid_antibacterials | 0.7 |
| enzyme_ampc_dha | pleuromutilins | 0.7 |
| enzyme_ampc_dha | other | 0.7 |
| target_site_pbp2a_meca | pen | 0.99 |
| target_site_pbp2a_meca | bli | 0.99 |
| target_site_pbp2a_meca | bli_anti_pseudomonal | 0.99 |
| target_site_pbp2a_meca | bli_sulbactam | 0.99 |
| target_site_pbp2a_meca | c1_2g | 0.99 |
| target_site_pbp2a_meca | c3g | 0.99 |
| target_site_pbp2a_meca | c3g_bli | 0.99 |
| target_site_pbp2a_meca | c4g | 0.7 |
| target_site_pbp2a_meca | anti_mrsa_ceph | 0.7 |
| target_site_pbp2a_meca | siderophore_ceph | 0.7 |
| target_site_pbp2a_meca | cft_avi | 0.99 |
| target_site_pbp2a_meca | mer_vab | 0.99 |
| target_site_pbp2a_meca | azt_avi | 0.99 |
| target_site_pbp2a_meca | carb_group1 | 0.85 |
| target_site_pbp2a_meca | carb_group2 | 0.85 |
| target_site_pbp2a_meca | mono | 0.99 |
| target_site_pbp2a_meca | fq | 0.99 |
| target_site_pbp2a_meca | ag_group1 | 0.99 |
| target_site_pbp2a_meca | ag_group2 | 0.99 |
| target_site_pbp2a_meca | mls | 0.99 |
| target_site_pbp2a_meca | lincosamides | 0.99 |
| target_site_pbp2a_meca | glyc | 0.99 |
| target_site_pbp2a_meca | lipoglycopeptides | 0.99 |
| target_site_pbp2a_meca | tet | 0.99 |
| target_site_pbp2a_meca | glycylcyclines | 0.99 |
| target_site_pbp2a_meca | poly | 0.99 |
| target_site_pbp2a_meca | oxa | 0.99 |
| target_site_pbp2a_meca | chl | 0.99 |
| target_site_pbp2a_meca | sulf | 0.99 |
| target_site_pbp2a_meca | lipopeptides | 0.99 |
| target_site_pbp2a_meca | streptogramins | 0.99 |
| target_site_pbp2a_meca | nitrofurans | 0.99 |
| target_site_pbp2a_meca | phosphonic_acids | 0.99 |
| target_site_pbp2a_meca | nitroimidazoles | 0.99 |
| target_site_pbp2a_meca | rifamycins | 0.99 |
| target_site_pbp2a_meca | macrocycles | 0.99 |
| target_site_pbp2a_meca | steroid_antibacterials | 0.99 |
| target_site_pbp2a_meca | pleuromutilins | 0.99 |
| target_site_pbp2a_meca | other | 0.99 |
| target_site_van_a | pen | 0.99 |
| target_site_van_a | bli | 0.99 |
| target_site_van_a | bli_anti_pseudomonal | 0.99 |
| target_site_van_a | bli_sulbactam | 0.99 |
| target_site_van_a | c1_2g | 0.99 |
| target_site_van_a | c3g | 0.99 |
| target_site_van_a | c3g_bli | 0.99 |
| target_site_van_a | c4g | 0.99 |
| target_site_van_a | anti_mrsa_ceph | 0.99 |
| target_site_van_a | siderophore_ceph | 0.99 |
| target_site_van_a | cft_avi | 0.99 |
| target_site_van_a | mer_vab | 0.99 |
| target_site_van_a | azt_avi | 0.99 |
| target_site_van_a | carb_group1 | 0.99 |
| target_site_van_a | carb_group2 | 0.99 |
| target_site_van_a | mono | 0.99 |
| target_site_van_a | fq | 0.99 |
| target_site_van_a | ag_group1 | 0.99 |
| target_site_van_a | ag_group2 | 0.99 |
| target_site_van_a | mls | 0.99 |
| target_site_van_a | lincosamides | 0.99 |
| target_site_van_a | glyc | 0.99 |
| target_site_van_a | lipoglycopeptides | 0.99 |
| target_site_van_a | tet | 0.99 |
| target_site_van_a | glycylcyclines | 0.99 |
| target_site_van_a | poly | 0.99 |
| target_site_van_a | oxa | 0.99 |
| target_site_van_a | chl | 0.99 |
| target_site_van_a | sulf | 0.99 |
| target_site_van_a | lipopeptides | 0.99 |
| target_site_van_a | streptogramins | 0.99 |
| target_site_van_a | nitrofurans | 0.99 |
| target_site_van_a | phosphonic_acids | 0.99 |
| target_site_van_a | nitroimidazoles | 0.99 |
| target_site_van_a | rifamycins | 0.99 |
| target_site_van_a | macrocycles | 0.99 |
| target_site_van_a | steroid_antibacterials | 0.99 |
| target_site_van_a | pleuromutilins | 0.99 |
| target_site_van_a | other | 0.99 |
| target_site_van_b | pen | 0.99 |
| target_site_van_b | bli | 0.99 |
| target_site_van_b | bli_anti_pseudomonal | 0.99 |
| target_site_van_b | bli_sulbactam | 0.99 |
| target_site_van_b | c1_2g | 0.99 |
| target_site_van_b | c3g | 0.99 |
| target_site_van_b | c3g_bli | 0.99 |
| target_site_van_b | c4g | 0.99 |
| target_site_van_b | anti_mrsa_ceph | 0.99 |
| target_site_van_b | siderophore_ceph | 0.99 |
| target_site_van_b | cft_avi | 0.99 |
| target_site_van_b | mer_vab | 0.99 |
| target_site_van_b | azt_avi | 0.99 |
| target_site_van_b | carb_group1 | 0.99 |
| target_site_van_b | carb_group2 | 0.99 |
| target_site_van_b | mono | 0.99 |
| target_site_van_b | fq | 0.99 |
| target_site_van_b | ag_group1 | 0.99 |
| target_site_van_b | ag_group2 | 0.99 |
| target_site_van_b | mls | 0.99 |
| target_site_van_b | lincosamides | 0.99 |
| target_site_van_b | glyc | 0.7 |
| target_site_van_b | lipoglycopeptides | 0.7 |
| target_site_van_b | tet | 0.99 |
| target_site_van_b | glycylcyclines | 0.99 |
| target_site_van_b | poly | 0.99 |
| target_site_van_b | oxa | 0.99 |
| target_site_van_b | chl | 0.99 |
| target_site_van_b | sulf | 0.99 |
| target_site_van_b | lipopeptides | 0.99 |
| target_site_van_b | streptogramins | 0.99 |
| target_site_van_b | nitrofurans | 0.99 |
| target_site_van_b | phosphonic_acids | 0.99 |
| target_site_van_b | nitroimidazoles | 0.99 |
| target_site_van_b | rifamycins | 0.99 |
| target_site_van_b | macrocycles | 0.99 |
| target_site_van_b | steroid_antibacterials | 0.99 |
| target_site_van_b | pleuromutilins | 0.99 |
| target_site_van_b | other | 0.99 |
| mutation_gyra_primary | pen | 0.4 |
| mutation_gyra_primary | bli | 0.4 |
| mutation_gyra_primary | bli_anti_pseudomonal | 0.4 |
| mutation_gyra_primary | bli_sulbactam | 0.4 |
| mutation_gyra_primary | c1_2g | 0.4 |
| mutation_gyra_primary | c3g | 0.4 |
| mutation_gyra_primary | c3g_bli | 0.4 |
| mutation_gyra_primary | c4g | 0.4 |
| mutation_gyra_primary | anti_mrsa_ceph | 0.4 |
| mutation_gyra_primary | siderophore_ceph | 0.4 |
| mutation_gyra_primary | cft_avi | 0.4 |
| mutation_gyra_primary | mer_vab | 0.4 |
| mutation_gyra_primary | azt_avi | 0.4 |
| mutation_gyra_primary | carb_group1 | 0.4 |
| mutation_gyra_primary | carb_group2 | 0.4 |
| mutation_gyra_primary | mono | 0.4 |
| mutation_gyra_primary | fq | 0.4 |
| mutation_gyra_primary | ag_group1 | 0.4 |
| mutation_gyra_primary | ag_group2 | 0.4 |
| mutation_gyra_primary | mls | 0.4 |
| mutation_gyra_primary | lincosamides | 0.4 |
| mutation_gyra_primary | glyc | 0.4 |
| mutation_gyra_primary | lipoglycopeptides | 0.4 |
| mutation_gyra_primary | tet | 0.4 |
| mutation_gyra_primary | glycylcyclines | 0.4 |
| mutation_gyra_primary | poly | 0.4 |
| mutation_gyra_primary | oxa | 0.4 |
| mutation_gyra_primary | chl | 0.4 |
| mutation_gyra_primary | sulf | 0.4 |
| mutation_gyra_primary | lipopeptides | 0.4 |
| mutation_gyra_primary | streptogramins | 0.4 |
| mutation_gyra_primary | nitrofurans | 0.4 |
| mutation_gyra_primary | phosphonic_acids | 0.4 |
| mutation_gyra_primary | nitroimidazoles | 0.4 |
| mutation_gyra_primary | rifamycins | 0.4 |
| mutation_gyra_primary | macrocycles | 0.4 |
| mutation_gyra_primary | steroid_antibacterials | 0.4 |
| mutation_gyra_primary | pleuromutilins | 0.4 |
| mutation_gyra_primary | other | 0.4 |
| mutation_gyra_parc_secondary | pen | 0.95 |
| mutation_gyra_parc_secondary | bli | 0.95 |
| mutation_gyra_parc_secondary | bli_anti_pseudomonal | 0.95 |
| mutation_gyra_parc_secondary | bli_sulbactam | 0.95 |
| mutation_gyra_parc_secondary | c1_2g | 0.95 |
| mutation_gyra_parc_secondary | c3g | 0.95 |
| mutation_gyra_parc_secondary | c3g_bli | 0.95 |
| mutation_gyra_parc_secondary | c4g | 0.95 |
| mutation_gyra_parc_secondary | anti_mrsa_ceph | 0.95 |
| mutation_gyra_parc_secondary | siderophore_ceph | 0.95 |
| mutation_gyra_parc_secondary | cft_avi | 0.95 |
| mutation_gyra_parc_secondary | mer_vab | 0.95 |
| mutation_gyra_parc_secondary | azt_avi | 0.95 |
| mutation_gyra_parc_secondary | carb_group1 | 0.95 |
| mutation_gyra_parc_secondary | carb_group2 | 0.95 |
| mutation_gyra_parc_secondary | mono | 0.95 |
| mutation_gyra_parc_secondary | fq | 0.95 |
| mutation_gyra_parc_secondary | ag_group1 | 0.95 |
| mutation_gyra_parc_secondary | ag_group2 | 0.95 |
| mutation_gyra_parc_secondary | mls | 0.95 |
| mutation_gyra_parc_secondary | lincosamides | 0.95 |
| mutation_gyra_parc_secondary | glyc | 0.95 |
| mutation_gyra_parc_secondary | lipoglycopeptides | 0.95 |
| mutation_gyra_parc_secondary | tet | 0.95 |
| mutation_gyra_parc_secondary | glycylcyclines | 0.95 |
| mutation_gyra_parc_secondary | poly | 0.95 |
| mutation_gyra_parc_secondary | oxa | 0.95 |
| mutation_gyra_parc_secondary | chl | 0.95 |
| mutation_gyra_parc_secondary | sulf | 0.95 |
| mutation_gyra_parc_secondary | lipopeptides | 0.95 |
| mutation_gyra_parc_secondary | streptogramins | 0.95 |
| mutation_gyra_parc_secondary | nitrofurans | 0.95 |
| mutation_gyra_parc_secondary | phosphonic_acids | 0.95 |
| mutation_gyra_parc_secondary | nitroimidazoles | 0.95 |
| mutation_gyra_parc_secondary | rifamycins | 0.95 |
| mutation_gyra_parc_secondary | macrocycles | 0.95 |
| mutation_gyra_parc_secondary | steroid_antibacterials | 0.95 |
| mutation_gyra_parc_secondary | pleuromutilins | 0.95 |
| mutation_gyra_parc_secondary | other | 0.95 |
| protection_qnr | pen | 0.2 |
| protection_qnr | bli | 0.2 |
| protection_qnr | bli_anti_pseudomonal | 0.2 |
| protection_qnr | bli_sulbactam | 0.2 |
| protection_qnr | c1_2g | 0.2 |
| protection_qnr | c3g | 0.2 |
| protection_qnr | c3g_bli | 0.2 |
| protection_qnr | c4g | 0.2 |
| protection_qnr | anti_mrsa_ceph | 0.2 |
| protection_qnr | siderophore_ceph | 0.2 |
| protection_qnr | cft_avi | 0.2 |
| protection_qnr | mer_vab | 0.2 |
| protection_qnr | azt_avi | 0.2 |
| protection_qnr | carb_group1 | 0.2 |
| protection_qnr | carb_group2 | 0.2 |
| protection_qnr | mono | 0.2 |
| protection_qnr | fq | 0.2 |
| protection_qnr | ag_group1 | 0.2 |
| protection_qnr | ag_group2 | 0.2 |
| protection_qnr | mls | 0.2 |
| protection_qnr | lincosamides | 0.2 |
| protection_qnr | glyc | 0.2 |
| protection_qnr | lipoglycopeptides | 0.2 |
| protection_qnr | tet | 0.2 |
| protection_qnr | glycylcyclines | 0.2 |
| protection_qnr | poly | 0.2 |
| protection_qnr | oxa | 0.2 |
| protection_qnr | chl | 0.2 |
| protection_qnr | sulf | 0.2 |
| protection_qnr | lipopeptides | 0.2 |
| protection_qnr | streptogramins | 0.2 |
| protection_qnr | nitrofurans | 0.2 |
| protection_qnr | phosphonic_acids | 0.2 |
| protection_qnr | nitroimidazoles | 0.2 |
| protection_qnr | rifamycins | 0.2 |
| protection_qnr | macrocycles | 0.2 |
| protection_qnr | steroid_antibacterials | 0.2 |
| protection_qnr | pleuromutilins | 0.2 |
| protection_qnr | other | 0.2 |
| enzyme_16s_rrmt | pen | 0.95 |
| enzyme_16s_rrmt | bli | 0.95 |
| enzyme_16s_rrmt | bli_anti_pseudomonal | 0.95 |
| enzyme_16s_rrmt | bli_sulbactam | 0.95 |
| enzyme_16s_rrmt | c1_2g | 0.95 |
| enzyme_16s_rrmt | c3g | 0.95 |
| enzyme_16s_rrmt | c3g_bli | 0.95 |
| enzyme_16s_rrmt | c4g | 0.95 |
| enzyme_16s_rrmt | anti_mrsa_ceph | 0.95 |
| enzyme_16s_rrmt | siderophore_ceph | 0.95 |
| enzyme_16s_rrmt | cft_avi | 0.95 |
| enzyme_16s_rrmt | mer_vab | 0.95 |
| enzyme_16s_rrmt | azt_avi | 0.95 |
| enzyme_16s_rrmt | carb_group1 | 0.95 |
| enzyme_16s_rrmt | carb_group2 | 0.95 |
| enzyme_16s_rrmt | mono | 0.95 |
| enzyme_16s_rrmt | fq | 0.95 |
| enzyme_16s_rrmt | ag_group1 | 0.95 |
| enzyme_16s_rrmt | ag_group2 | 0.95 |
| enzyme_16s_rrmt | mls | 0.95 |
| enzyme_16s_rrmt | lincosamides | 0.95 |
| enzyme_16s_rrmt | glyc | 0.95 |
| enzyme_16s_rrmt | lipoglycopeptides | 0.95 |
| enzyme_16s_rrmt | tet | 0.95 |
| enzyme_16s_rrmt | glycylcyclines | 0.95 |
| enzyme_16s_rrmt | poly | 0.95 |
| enzyme_16s_rrmt | oxa | 0.95 |
| enzyme_16s_rrmt | chl | 0.95 |
| enzyme_16s_rrmt | sulf | 0.95 |
| enzyme_16s_rrmt | lipopeptides | 0.95 |
| enzyme_16s_rrmt | streptogramins | 0.95 |
| enzyme_16s_rrmt | nitrofurans | 0.95 |
| enzyme_16s_rrmt | phosphonic_acids | 0.95 |
| enzyme_16s_rrmt | nitroimidazoles | 0.95 |
| enzyme_16s_rrmt | rifamycins | 0.95 |
| enzyme_16s_rrmt | macrocycles | 0.95 |
| enzyme_16s_rrmt | steroid_antibacterials | 0.95 |
| enzyme_16s_rrmt | pleuromutilins | 0.95 |
| enzyme_16s_rrmt | other | 0.95 |
| target_site_erm_b | pen | 0.9 |
| target_site_erm_b | bli | 0.9 |
| target_site_erm_b | bli_anti_pseudomonal | 0.9 |
| target_site_erm_b | bli_sulbactam | 0.9 |
| target_site_erm_b | c1_2g | 0.9 |
| target_site_erm_b | c3g | 0.9 |
| target_site_erm_b | c3g_bli | 0.9 |
| target_site_erm_b | c4g | 0.9 |
| target_site_erm_b | anti_mrsa_ceph | 0.9 |
| target_site_erm_b | siderophore_ceph | 0.9 |
| target_site_erm_b | cft_avi | 0.9 |
| target_site_erm_b | mer_vab | 0.9 |
| target_site_erm_b | azt_avi | 0.9 |
| target_site_erm_b | carb_group1 | 0.9 |
| target_site_erm_b | carb_group2 | 0.9 |
| target_site_erm_b | mono | 0.9 |
| target_site_erm_b | fq | 0.9 |
| target_site_erm_b | ag_group1 | 0.9 |
| target_site_erm_b | ag_group2 | 0.9 |
| target_site_erm_b | mls | 0.9 |
| target_site_erm_b | lincosamides | 0.9 |
| target_site_erm_b | glyc | 0.9 |
| target_site_erm_b | lipoglycopeptides | 0.9 |
| target_site_erm_b | tet | 0.9 |
| target_site_erm_b | glycylcyclines | 0.9 |
| target_site_erm_b | poly | 0.9 |
| target_site_erm_b | oxa | 0.9 |
| target_site_erm_b | chl | 0.9 |
| target_site_erm_b | sulf | 0.9 |
| target_site_erm_b | lipopeptides | 0.9 |
| target_site_erm_b | streptogramins | 0.9 |
| target_site_erm_b | nitrofurans | 0.9 |
| target_site_erm_b | phosphonic_acids | 0.9 |
| target_site_erm_b | nitroimidazoles | 0.9 |
| target_site_erm_b | rifamycins | 0.9 |
| target_site_erm_b | macrocycles | 0.9 |
| target_site_erm_b | steroid_antibacterials | 0.9 |
| target_site_erm_b | pleuromutilins | 0.9 |
| target_site_erm_b | other | 0.9 |
| target_site_cfr | pen | 0.95 |
| target_site_cfr | bli | 0.95 |
| target_site_cfr | bli_anti_pseudomonal | 0.95 |
| target_site_cfr | bli_sulbactam | 0.95 |
| target_site_cfr | c1_2g | 0.95 |
| target_site_cfr | c3g | 0.95 |
| target_site_cfr | c3g_bli | 0.95 |
| target_site_cfr | c4g | 0.95 |
| target_site_cfr | anti_mrsa_ceph | 0.95 |
| target_site_cfr | siderophore_ceph | 0.95 |
| target_site_cfr | cft_avi | 0.95 |
| target_site_cfr | mer_vab | 0.95 |
| target_site_cfr | azt_avi | 0.95 |
| target_site_cfr | carb_group1 | 0.95 |
| target_site_cfr | carb_group2 | 0.95 |
| target_site_cfr | mono | 0.95 |
| target_site_cfr | fq | 0.95 |
| target_site_cfr | ag_group1 | 0.95 |
| target_site_cfr | ag_group2 | 0.95 |
| target_site_cfr | mls | 0.7 |
| target_site_cfr | lincosamides | 0.7 |
| target_site_cfr | glyc | 0.95 |
| target_site_cfr | lipoglycopeptides | 0.95 |
| target_site_cfr | tet | 0.95 |
| target_site_cfr | glycylcyclines | 0.95 |
| target_site_cfr | poly | 0.95 |
| target_site_cfr | oxa | 0.9 |
| target_site_cfr | chl | 0.7 |
| target_site_cfr | sulf | 0.95 |
| target_site_cfr | lipopeptides | 0.95 |
| target_site_cfr | streptogramins | 0.95 |
| target_site_cfr | nitrofurans | 0.95 |
| target_site_cfr | phosphonic_acids | 0.95 |
| target_site_cfr | nitroimidazoles | 0.95 |
| target_site_cfr | rifamycins | 0.95 |
| target_site_cfr | macrocycles | 0.95 |
| target_site_cfr | steroid_antibacterials | 0.95 |
| target_site_cfr | pleuromutilins | 0.95 |
| target_site_cfr | other | 0.95 |
| enzyme_cat | pen | 0.9 |
| enzyme_cat | bli | 0.9 |
| enzyme_cat | bli_anti_pseudomonal | 0.9 |
| enzyme_cat | bli_sulbactam | 0.9 |
| enzyme_cat | c1_2g | 0.9 |
| enzyme_cat | c3g | 0.9 |
| enzyme_cat | c3g_bli | 0.9 |
| enzyme_cat | c4g | 0.9 |
| enzyme_cat | anti_mrsa_ceph | 0.9 |
| enzyme_cat | siderophore_ceph | 0.9 |
| enzyme_cat | cft_avi | 0.9 |
| enzyme_cat | mer_vab | 0.9 |
| enzyme_cat | azt_avi | 0.9 |
| enzyme_cat | carb_group1 | 0.9 |
| enzyme_cat | carb_group2 | 0.9 |
| enzyme_cat | mono | 0.9 |
| enzyme_cat | fq | 0.9 |
| enzyme_cat | ag_group1 | 0.9 |
| enzyme_cat | ag_group2 | 0.9 |
| enzyme_cat | mls | 0.9 |
| enzyme_cat | lincosamides | 0.9 |
| enzyme_cat | glyc | 0.9 |
| enzyme_cat | lipoglycopeptides | 0.9 |
| enzyme_cat | tet | 0.9 |
| enzyme_cat | glycylcyclines | 0.9 |
| enzyme_cat | poly | 0.9 |
| enzyme_cat | oxa | 0.9 |
| enzyme_cat | chl | 0.9 |
| enzyme_cat | sulf | 0.9 |
| enzyme_cat | lipopeptides | 0.9 |
| enzyme_cat | streptogramins | 0.9 |
| enzyme_cat | nitrofurans | 0.9 |
| enzyme_cat | phosphonic_acids | 0.9 |
| enzyme_cat | nitroimidazoles | 0.9 |
| enzyme_cat | rifamycins | 0.9 |
| enzyme_cat | macrocycles | 0.9 |
| enzyme_cat | steroid_antibacterials | 0.9 |
| enzyme_cat | pleuromutilins | 0.9 |
| enzyme_cat | other | 0.9 |
| efflux_acrab_tolc | c1_2g | 0.3 |
| efflux_acrab_tolc | cft_avi | 0.3 |
| efflux_acrab_tolc | mer_vab | 0.3 |
| efflux_acrab_tolc | azt_avi | 0.3 |
| efflux_acrab_tolc | carb_group1 | 0.3 |
| efflux_acrab_tolc | carb_group2 | 0.3 |
| efflux_acrab_tolc | mono | 0.3 |
| efflux_acrab_tolc | fq | 0.25 |
| efflux_acrab_tolc | glyc | 0.3 |
| efflux_acrab_tolc | lipoglycopeptides | 0.3 |
| efflux_acrab_tolc | tet | 0.25 |
| efflux_acrab_tolc | glycylcyclines | 0.25 |
| efflux_acrab_tolc | poly | 0.3 |
| efflux_acrab_tolc | oxa | 0.3 |
| efflux_acrab_tolc | chl | 0.2 |
| efflux_acrab_tolc | sulf | 0.3 |
| efflux_acrab_tolc | lipopeptides | 0.3 |
| efflux_acrab_tolc | streptogramins | 0.3 |
| efflux_acrab_tolc | nitrofurans | 0.3 |
| efflux_acrab_tolc | phosphonic_acids | 0.3 |
| efflux_acrab_tolc | nitroimidazoles | 0.3 |
| efflux_acrab_tolc | rifamycins | 0.3 |
| efflux_acrab_tolc | macrocycles | 0.3 |
| efflux_acrab_tolc | steroid_antibacterials | 0.3 |
| efflux_acrab_tolc | pleuromutilins | 0.3 |
| efflux_acrab_tolc | other | 0.3 |
| efflux_mexxy_oprm | pen | 0.3 |
| efflux_mexxy_oprm | bli | 0.3 |
| efflux_mexxy_oprm | bli_anti_pseudomonal | 0.3 |
| efflux_mexxy_oprm | bli_sulbactam | 0.3 |
| efflux_mexxy_oprm | c1_2g | 0.3 |
| efflux_mexxy_oprm | c3g | 0.3 |
| efflux_mexxy_oprm | c3g_bli | 0.3 |
| efflux_mexxy_oprm | c4g | 0.2 |
| efflux_mexxy_oprm | anti_mrsa_ceph | 0.2 |
| efflux_mexxy_oprm | siderophore_ceph | 0.2 |
| efflux_mexxy_oprm | cft_avi | 0.3 |
| efflux_mexxy_oprm | mer_vab | 0.3 |
| efflux_mexxy_oprm | azt_avi | 0.3 |
| efflux_mexxy_oprm | carb_group1 | 0.05 |
| efflux_mexxy_oprm | carb_group2 | 0.05 |
| efflux_mexxy_oprm | mono | 0.3 |
| efflux_mexxy_oprm | fq | 0.2 |
| efflux_mexxy_oprm | ag_group1 | 0.3 |
| efflux_mexxy_oprm | ag_group2 | 0.3 |
| efflux_mexxy_oprm | mls | 0.3 |
| efflux_mexxy_oprm | lincosamides | 0.3 |
| efflux_mexxy_oprm | glyc | 0.3 |
| efflux_mexxy_oprm | lipoglycopeptides | 0.3 |
| efflux_mexxy_oprm | tet | 0.3 |
| efflux_mexxy_oprm | glycylcyclines | 0.3 |
| efflux_mexxy_oprm | poly | 0.3 |
| efflux_mexxy_oprm | oxa | 0.3 |
| efflux_mexxy_oprm | chl | 0.3 |
| efflux_mexxy_oprm | sulf | 0.3 |
| efflux_mexxy_oprm | lipopeptides | 0.3 |
| efflux_mexxy_oprm | streptogramins | 0.3 |
| efflux_mexxy_oprm | nitrofurans | 0.3 |
| efflux_mexxy_oprm | phosphonic_acids | 0.3 |
| efflux_mexxy_oprm | nitroimidazoles | 0.3 |
| efflux_mexxy_oprm | rifamycins | 0.3 |
| efflux_mexxy_oprm | macrocycles | 0.3 |
| efflux_mexxy_oprm | steroid_antibacterials | 0.3 |
| efflux_mexxy_oprm | pleuromutilins | 0.3 |
| efflux_mexxy_oprm | other | 0.3 |
| porin_loss_ompk35_36 | pen | 0.3 |
| porin_loss_ompk35_36 | bli | 0.4 |
| porin_loss_ompk35_36 | bli_anti_pseudomonal | 0.4 |
| porin_loss_ompk35_36 | bli_sulbactam | 0.4 |
| porin_loss_ompk35_36 | c1_2g | 0.8 |
| porin_loss_ompk35_36 | c3g | 0.4 |
| porin_loss_ompk35_36 | c3g_bli | 0.4 |
| porin_loss_ompk35_36 | c4g | 0.3 |
| porin_loss_ompk35_36 | anti_mrsa_ceph | 0.3 |
| porin_loss_ompk35_36 | siderophore_ceph | 0.3 |
| porin_loss_ompk35_36 | cft_avi | 0.25 |
| porin_loss_ompk35_36 | mer_vab | 0.25 |
| porin_loss_ompk35_36 | azt_avi | 0.25 |
| porin_loss_ompk35_36 | mono | 0.8 |
| porin_loss_ompk35_36 | mls | 0.8 |
| porin_loss_ompk35_36 | lincosamides | 0.8 |
| porin_loss_ompk35_36 | glyc | 0.8 |
| porin_loss_ompk35_36 | lipoglycopeptides | 0.8 |
| porin_loss_ompk35_36 | tet | 0.8 |
| porin_loss_ompk35_36 | glycylcyclines | 0.8 |
| porin_loss_ompk35_36 | poly | 0.8 |
| porin_loss_ompk35_36 | oxa | 0.8 |
| porin_loss_ompk35_36 | chl | 0.8 |
| porin_loss_ompk35_36 | sulf | 0.8 |
| porin_loss_ompk35_36 | lipopeptides | 0.8 |
| porin_loss_ompk35_36 | streptogramins | 0.8 |
| porin_loss_ompk35_36 | nitrofurans | 0.8 |
| porin_loss_ompk35_36 | phosphonic_acids | 0.8 |
| porin_loss_ompk35_36 | nitroimidazoles | 0.8 |
| porin_loss_ompk35_36 | rifamycins | 0.8 |
| porin_loss_ompk35_36 | macrocycles | 0.8 |
| porin_loss_ompk35_36 | steroid_antibacterials | 0.8 |
| porin_loss_ompk35_36 | pleuromutilins | 0.8 |
| porin_loss_ompk35_36 | other | 0.8 |
| porin_loss_oprd | pen | 0.8 |
| porin_loss_oprd | bli | 0.8 |
| porin_loss_oprd | bli_anti_pseudomonal | 0.8 |
| porin_loss_oprd | bli_sulbactam | 0.8 |
| porin_loss_oprd | c1_2g | 0.8 |
| porin_loss_oprd | c3g | 0.8 |
| porin_loss_oprd | c3g_bli | 0.8 |
| porin_loss_oprd | c4g | 0.8 |
| porin_loss_oprd | anti_mrsa_ceph | 0.8 |
| porin_loss_oprd | siderophore_ceph | 0.8 |
| porin_loss_oprd | cft_avi | 0.8 |
| porin_loss_oprd | mer_vab | 0.8 |
| porin_loss_oprd | azt_avi | 0.8 |
| porin_loss_oprd | carb_group1 | 0.8 |
| porin_loss_oprd | carb_group2 | 0.8 |
| porin_loss_oprd | mono | 0.8 |
| porin_loss_oprd | fq | 0.8 |
| porin_loss_oprd | ag_group1 | 0.8 |
| porin_loss_oprd | ag_group2 | 0.8 |
| porin_loss_oprd | mls | 0.8 |
| porin_loss_oprd | lincosamides | 0.8 |
| porin_loss_oprd | glyc | 0.8 |
| porin_loss_oprd | lipoglycopeptides | 0.8 |
| porin_loss_oprd | tet | 0.8 |
| porin_loss_oprd | glycylcyclines | 0.8 |
| porin_loss_oprd | poly | 0.8 |
| porin_loss_oprd | oxa | 0.8 |
| porin_loss_oprd | chl | 0.8 |
| porin_loss_oprd | sulf | 0.8 |
| porin_loss_oprd | lipopeptides | 0.8 |
| porin_loss_oprd | streptogramins | 0.8 |
| porin_loss_oprd | nitrofurans | 0.8 |
| porin_loss_oprd | phosphonic_acids | 0.8 |
| porin_loss_oprd | nitroimidazoles | 0.8 |
| porin_loss_oprd | rifamycins | 0.8 |
| porin_loss_oprd | macrocycles | 0.8 |
| porin_loss_oprd | steroid_antibacterials | 0.8 |
| porin_loss_oprd | pleuromutilins | 0.8 |
| porin_loss_oprd | other | 0.8 |
| modification_mcr_1 | pen | 0.85 |
| modification_mcr_1 | bli | 0.85 |
| modification_mcr_1 | bli_anti_pseudomonal | 0.85 |
| modification_mcr_1 | bli_sulbactam | 0.85 |
| modification_mcr_1 | c1_2g | 0.85 |
| modification_mcr_1 | c3g | 0.85 |
| modification_mcr_1 | c3g_bli | 0.85 |
| modification_mcr_1 | c4g | 0.85 |
| modification_mcr_1 | anti_mrsa_ceph | 0.85 |
| modification_mcr_1 | siderophore_ceph | 0.85 |
| modification_mcr_1 | cft_avi | 0.85 |
| modification_mcr_1 | mer_vab | 0.85 |
| modification_mcr_1 | azt_avi | 0.85 |
| modification_mcr_1 | carb_group1 | 0.85 |
| modification_mcr_1 | carb_group2 | 0.85 |
| modification_mcr_1 | mono | 0.85 |
| modification_mcr_1 | fq | 0.85 |
| modification_mcr_1 | ag_group1 | 0.85 |
| modification_mcr_1 | ag_group2 | 0.85 |
| modification_mcr_1 | mls | 0.85 |
| modification_mcr_1 | lincosamides | 0.85 |
| modification_mcr_1 | glyc | 0.85 |
| modification_mcr_1 | lipoglycopeptides | 0.85 |
| modification_mcr_1 | tet | 0.85 |
| modification_mcr_1 | glycylcyclines | 0.85 |
| modification_mcr_1 | poly | 0.85 |
| modification_mcr_1 | oxa | 0.85 |
| modification_mcr_1 | chl | 0.85 |
| modification_mcr_1 | sulf | 0.85 |
| modification_mcr_1 | lipopeptides | 0.85 |
| modification_mcr_1 | streptogramins | 0.85 |
| modification_mcr_1 | nitrofurans | 0.85 |
| modification_mcr_1 | phosphonic_acids | 0.85 |
| modification_mcr_1 | nitroimidazoles | 0.85 |
| modification_mcr_1 | rifamycins | 0.85 |
| modification_mcr_1 | macrocycles | 0.85 |
| modification_mcr_1 | steroid_antibacterials | 0.85 |
| modification_mcr_1 | pleuromutilins | 0.85 |
| modification_mcr_1 | other | 0.85 |
| global_efflux_pump | c1_2g | 0.2 |
| global_efflux_pump | cft_avi | 0.2 |
| global_efflux_pump | mer_vab | 0.2 |
| global_efflux_pump | azt_avi | 0.2 |
| global_efflux_pump | carb_group1 | 0.2 |
| global_efflux_pump | carb_group2 | 0.2 |
| global_efflux_pump | mono | 0.2 |
| global_efflux_pump | fq | 0.15 |
| global_efflux_pump | mls | 0.1 |
| global_efflux_pump | lincosamides | 0.1 |
| global_efflux_pump | glyc | 0.2 |
| global_efflux_pump | lipoglycopeptides | 0.2 |
| global_efflux_pump | tet | 0.15 |
| global_efflux_pump | glycylcyclines | 0.15 |
| global_efflux_pump | poly | 0.2 |
| global_efflux_pump | oxa | 0.2 |
| global_efflux_pump | chl | 0.15 |
| global_efflux_pump | sulf | 0.2 |
| global_efflux_pump | lipopeptides | 0.2 |
| global_efflux_pump | streptogramins | 0.2 |
| global_efflux_pump | nitrofurans | 0.2 |
| global_efflux_pump | phosphonic_acids | 0.2 |
| global_efflux_pump | nitroimidazoles | 0.2 |
| global_efflux_pump | rifamycins | 0.2 |
| global_efflux_pump | macrocycles | 0.2 |
| global_efflux_pump | steroid_antibacterials | 0.2 |
| global_efflux_pump | pleuromutilins | 0.2 |
| global_efflux_pump | other | 0.2 |
| global_porin_loss | c1_2g | 0.2 |
| global_porin_loss | cft_avi | 0.2 |
| global_porin_loss | mer_vab | 0.2 |
| global_porin_loss | azt_avi | 0.2 |
| global_porin_loss | mono | 0.2 |
| global_porin_loss | mls | 0.2 |
| global_porin_loss | lincosamides | 0.2 |
| global_porin_loss | glyc | 0.2 |
| global_porin_loss | lipoglycopeptides | 0.2 |
| global_porin_loss | tet | 0.2 |
| global_porin_loss | glycylcyclines | 0.2 |
| global_porin_loss | poly | 0.2 |
| global_porin_loss | oxa | 0.2 |
| global_porin_loss | chl | 0.2 |
| global_porin_loss | sulf | 0.2 |
| global_porin_loss | lipopeptides | 0.2 |
| global_porin_loss | streptogramins | 0.2 |
| global_porin_loss | nitrofurans | 0.2 |
| global_porin_loss | phosphonic_acids | 0.2 |
| global_porin_loss | nitroimidazoles | 0.2 |
| global_porin_loss | rifamycins | 0.2 |
| global_porin_loss | macrocycles | 0.2 |
| global_porin_loss | steroid_antibacterials | 0.2 |
| global_porin_loss | pleuromutilins | 0.2 |
| global_porin_loss | other | 0.2 |
| mutation_folate_pathway | pen | 0.85 |
| mutation_folate_pathway | bli | 0.85 |
| mutation_folate_pathway | bli_anti_pseudomonal | 0.85 |
| mutation_folate_pathway | bli_sulbactam | 0.85 |
| mutation_folate_pathway | c1_2g | 0.85 |
| mutation_folate_pathway | c3g | 0.85 |
| mutation_folate_pathway | c3g_bli | 0.85 |
| mutation_folate_pathway | c4g | 0.85 |
| mutation_folate_pathway | anti_mrsa_ceph | 0.85 |
| mutation_folate_pathway | siderophore_ceph | 0.85 |
| mutation_folate_pathway | cft_avi | 0.85 |
| mutation_folate_pathway | mer_vab | 0.85 |
| mutation_folate_pathway | azt_avi | 0.85 |
| mutation_folate_pathway | carb_group1 | 0.85 |
| mutation_folate_pathway | carb_group2 | 0.85 |
| mutation_folate_pathway | mono | 0.85 |
| mutation_folate_pathway | fq | 0.85 |
| mutation_folate_pathway | ag_group1 | 0.85 |
| mutation_folate_pathway | ag_group2 | 0.85 |
| mutation_folate_pathway | mls | 0.85 |
| mutation_folate_pathway | lincosamides | 0.85 |
| mutation_folate_pathway | glyc | 0.85 |
| mutation_folate_pathway | lipoglycopeptides | 0.85 |
| mutation_folate_pathway | tet | 0.85 |
| mutation_folate_pathway | glycylcyclines | 0.85 |
| mutation_folate_pathway | poly | 0.85 |
| mutation_folate_pathway | oxa | 0.85 |
| mutation_folate_pathway | chl | 0.85 |
| mutation_folate_pathway | sulf | 0.85 |
| mutation_folate_pathway | lipopeptides | 0.85 |
| mutation_folate_pathway | streptogramins | 0.85 |
| mutation_folate_pathway | nitrofurans | 0.85 |
| mutation_folate_pathway | phosphonic_acids | 0.85 |
| mutation_folate_pathway | nitroimidazoles | 0.85 |
| mutation_folate_pathway | rifamycins | 0.85 |
| mutation_folate_pathway | macrocycles | 0.85 |
| mutation_folate_pathway | steroid_antibacterials | 0.85 |
| mutation_folate_pathway | pleuromutilins | 0.85 |
| mutation_folate_pathway | other | 0.85 |
| mutation_nitroreductase | pen | 0.7 |
| mutation_nitroreductase | bli | 0.7 |
| mutation_nitroreductase | bli_anti_pseudomonal | 0.7 |
| mutation_nitroreductase | bli_sulbactam | 0.7 |
| mutation_nitroreductase | c1_2g | 0.7 |
| mutation_nitroreductase | c3g | 0.7 |
| mutation_nitroreductase | c3g_bli | 0.7 |
| mutation_nitroreductase | c4g | 0.7 |
| mutation_nitroreductase | anti_mrsa_ceph | 0.7 |
| mutation_nitroreductase | siderophore_ceph | 0.7 |
| mutation_nitroreductase | cft_avi | 0.7 |
| mutation_nitroreductase | mer_vab | 0.7 |
| mutation_nitroreductase | azt_avi | 0.7 |
| mutation_nitroreductase | carb_group1 | 0.7 |
| mutation_nitroreductase | carb_group2 | 0.7 |
| mutation_nitroreductase | mono | 0.7 |
| mutation_nitroreductase | fq | 0.7 |
| mutation_nitroreductase | ag_group1 | 0.7 |
| mutation_nitroreductase | ag_group2 | 0.7 |
| mutation_nitroreductase | mls | 0.7 |
| mutation_nitroreductase | lincosamides | 0.7 |
| mutation_nitroreductase | glyc | 0.7 |
| mutation_nitroreductase | lipoglycopeptides | 0.7 |
| mutation_nitroreductase | tet | 0.7 |
| mutation_nitroreductase | glycylcyclines | 0.7 |
| mutation_nitroreductase | poly | 0.7 |
| mutation_nitroreductase | oxa | 0.7 |
| mutation_nitroreductase | chl | 0.7 |
| mutation_nitroreductase | sulf | 0.7 |
| mutation_nitroreductase | lipopeptides | 0.7 |
| mutation_nitroreductase | streptogramins | 0.7 |
| mutation_nitroreductase | nitrofurans | 0.7 |
| mutation_nitroreductase | phosphonic_acids | 0.7 |
| mutation_nitroreductase | nitroimidazoles | 0.7 |
| mutation_nitroreductase | rifamycins | 0.7 |
| mutation_nitroreductase | macrocycles | 0.7 |
| mutation_nitroreductase | steroid_antibacterials | 0.7 |
| mutation_nitroreductase | pleuromutilins | 0.7 |
| mutation_nitroreductase | other | 0.7 |
| enzyme_fos_a | pen | 0.8 |
| enzyme_fos_a | bli | 0.8 |
| enzyme_fos_a | bli_anti_pseudomonal | 0.8 |
| enzyme_fos_a | bli_sulbactam | 0.8 |
| enzyme_fos_a | c1_2g | 0.8 |
| enzyme_fos_a | c3g | 0.8 |
| enzyme_fos_a | c3g_bli | 0.8 |
| enzyme_fos_a | c4g | 0.8 |
| enzyme_fos_a | anti_mrsa_ceph | 0.8 |
| enzyme_fos_a | siderophore_ceph | 0.8 |
| enzyme_fos_a | cft_avi | 0.8 |
| enzyme_fos_a | mer_vab | 0.8 |
| enzyme_fos_a | azt_avi | 0.8 |
| enzyme_fos_a | carb_group1 | 0.8 |
| enzyme_fos_a | carb_group2 | 0.8 |
| enzyme_fos_a | mono | 0.8 |
| enzyme_fos_a | fq | 0.8 |
| enzyme_fos_a | ag_group1 | 0.8 |
| enzyme_fos_a | ag_group2 | 0.8 |
| enzyme_fos_a | mls | 0.8 |
| enzyme_fos_a | lincosamides | 0.8 |
| enzyme_fos_a | glyc | 0.8 |
| enzyme_fos_a | lipoglycopeptides | 0.8 |
| enzyme_fos_a | tet | 0.8 |
| enzyme_fos_a | glycylcyclines | 0.8 |
| enzyme_fos_a | poly | 0.8 |
| enzyme_fos_a | oxa | 0.8 |
| enzyme_fos_a | chl | 0.8 |
| enzyme_fos_a | sulf | 0.8 |
| enzyme_fos_a | lipopeptides | 0.8 |
| enzyme_fos_a | streptogramins | 0.8 |
| enzyme_fos_a | nitrofurans | 0.8 |
| enzyme_fos_a | phosphonic_acids | 0.8 |
| enzyme_fos_a | nitroimidazoles | 0.8 |
| enzyme_fos_a | rifamycins | 0.8 |
| enzyme_fos_a | macrocycles | 0.8 |
| enzyme_fos_a | steroid_antibacterials | 0.8 |
| enzyme_fos_a | pleuromutilins | 0.8 |
| enzyme_fos_a | other | 0.8 |
| mutation_mpr_f | pen | 0.6 |
| mutation_mpr_f | bli | 0.6 |
| mutation_mpr_f | bli_anti_pseudomonal | 0.6 |
| mutation_mpr_f | bli_sulbactam | 0.6 |
| mutation_mpr_f | c1_2g | 0.6 |
| mutation_mpr_f | c3g | 0.6 |
| mutation_mpr_f | c3g_bli | 0.6 |
| mutation_mpr_f | c4g | 0.6 |
| mutation_mpr_f | anti_mrsa_ceph | 0.6 |
| mutation_mpr_f | siderophore_ceph | 0.6 |
| mutation_mpr_f | cft_avi | 0.6 |
| mutation_mpr_f | mer_vab | 0.6 |
| mutation_mpr_f | azt_avi | 0.6 |
| mutation_mpr_f | carb_group1 | 0.6 |
| mutation_mpr_f | carb_group2 | 0.6 |
| mutation_mpr_f | mono | 0.6 |
| mutation_mpr_f | fq | 0.6 |
| mutation_mpr_f | ag_group1 | 0.6 |
| mutation_mpr_f | ag_group2 | 0.6 |
| mutation_mpr_f | mls | 0.6 |
| mutation_mpr_f | lincosamides | 0.6 |
| mutation_mpr_f | glyc | 0.6 |
| mutation_mpr_f | lipoglycopeptides | 0.6 |
| mutation_mpr_f | tet | 0.6 |
| mutation_mpr_f | glycylcyclines | 0.6 |
| mutation_mpr_f | poly | 0.6 |
| mutation_mpr_f | oxa | 0.6 |
| mutation_mpr_f | chl | 0.6 |
| mutation_mpr_f | sulf | 0.6 |
| mutation_mpr_f | lipopeptides | 0.6 |
| mutation_mpr_f | streptogramins | 0.6 |
| mutation_mpr_f | nitrofurans | 0.6 |
| mutation_mpr_f | phosphonic_acids | 0.6 |
| mutation_mpr_f | nitroimidazoles | 0.6 |
| mutation_mpr_f | rifamycins | 0.6 |
| mutation_mpr_f | macrocycles | 0.6 |
| mutation_mpr_f | steroid_antibacterials | 0.6 |
| mutation_mpr_f | pleuromutilins | 0.6 |
| mutation_mpr_f | other | 0.6 |
| mutation_rpo_b | pen | 0.95 |
| mutation_rpo_b | bli | 0.95 |
| mutation_rpo_b | bli_anti_pseudomonal | 0.95 |
| mutation_rpo_b | bli_sulbactam | 0.95 |
| mutation_rpo_b | c1_2g | 0.95 |
| mutation_rpo_b | c3g | 0.95 |
| mutation_rpo_b | c3g_bli | 0.95 |
| mutation_rpo_b | c4g | 0.95 |
| mutation_rpo_b | anti_mrsa_ceph | 0.95 |
| mutation_rpo_b | siderophore_ceph | 0.95 |
| mutation_rpo_b | cft_avi | 0.95 |
| mutation_rpo_b | mer_vab | 0.95 |
| mutation_rpo_b | azt_avi | 0.95 |
| mutation_rpo_b | carb_group1 | 0.95 |
| mutation_rpo_b | carb_group2 | 0.95 |
| mutation_rpo_b | mono | 0.95 |
| mutation_rpo_b | fq | 0.95 |
| mutation_rpo_b | ag_group1 | 0.95 |
| mutation_rpo_b | ag_group2 | 0.95 |
| mutation_rpo_b | mls | 0.95 |
| mutation_rpo_b | lincosamides | 0.95 |
| mutation_rpo_b | glyc | 0.95 |
| mutation_rpo_b | lipoglycopeptides | 0.95 |
| mutation_rpo_b | tet | 0.95 |
| mutation_rpo_b | glycylcyclines | 0.95 |
| mutation_rpo_b | poly | 0.95 |
| mutation_rpo_b | oxa | 0.95 |
| mutation_rpo_b | chl | 0.95 |
| mutation_rpo_b | sulf | 0.95 |
| mutation_rpo_b | lipopeptides | 0.95 |
| mutation_rpo_b | streptogramins | 0.95 |
| mutation_rpo_b | nitrofurans | 0.95 |
| mutation_rpo_b | phosphonic_acids | 0.95 |
| mutation_rpo_b | nitroimidazoles | 0.95 |
| mutation_rpo_b | rifamycins | 0.95 |
| mutation_rpo_b | macrocycles | 0.95 |
| mutation_rpo_b | steroid_antibacterials | 0.95 |
| mutation_rpo_b | pleuromutilins | 0.95 |
| mutation_rpo_b | other | 0.95 |
| protection_fus_b | pen | 0.7 |
| protection_fus_b | bli | 0.7 |
| protection_fus_b | bli_anti_pseudomonal | 0.7 |
| protection_fus_b | bli_sulbactam | 0.7 |
| protection_fus_b | c1_2g | 0.7 |
| protection_fus_b | c3g | 0.7 |
| protection_fus_b | c3g_bli | 0.7 |
| protection_fus_b | c4g | 0.7 |
| protection_fus_b | anti_mrsa_ceph | 0.7 |
| protection_fus_b | siderophore_ceph | 0.7 |
| protection_fus_b | cft_avi | 0.7 |
| protection_fus_b | mer_vab | 0.7 |
| protection_fus_b | azt_avi | 0.7 |
| protection_fus_b | carb_group1 | 0.7 |
| protection_fus_b | carb_group2 | 0.7 |
| protection_fus_b | mono | 0.7 |
| protection_fus_b | fq | 0.7 |
| protection_fus_b | ag_group1 | 0.7 |
| protection_fus_b | ag_group2 | 0.7 |
| protection_fus_b | mls | 0.7 |
| protection_fus_b | lincosamides | 0.7 |
| protection_fus_b | glyc | 0.7 |
| protection_fus_b | lipoglycopeptides | 0.7 |
| protection_fus_b | tet | 0.7 |
| protection_fus_b | glycylcyclines | 0.7 |
| protection_fus_b | poly | 0.7 |
| protection_fus_b | oxa | 0.7 |
| protection_fus_b | chl | 0.7 |
| protection_fus_b | sulf | 0.7 |
| protection_fus_b | lipopeptides | 0.7 |
| protection_fus_b | streptogramins | 0.7 |
| protection_fus_b | nitrofurans | 0.7 |
| protection_fus_b | phosphonic_acids | 0.7 |
| protection_fus_b | nitroimidazoles | 0.7 |
| protection_fus_b | rifamycins | 0.7 |
| protection_fus_b | macrocycles | 0.7 |
| protection_fus_b | steroid_antibacterials | 0.7 |
| protection_fus_b | pleuromutilins | 0.7 |
| protection_fus_b | other | 0.7 |
| protection_tet_m | pen | 0.9 |
| protection_tet_m | bli | 0.9 |
| protection_tet_m | bli_anti_pseudomonal | 0.9 |
| protection_tet_m | bli_sulbactam | 0.9 |
| protection_tet_m | c1_2g | 0.9 |
| protection_tet_m | c3g | 0.9 |
| protection_tet_m | c3g_bli | 0.9 |
| protection_tet_m | c4g | 0.9 |
| protection_tet_m | anti_mrsa_ceph | 0.9 |
| protection_tet_m | siderophore_ceph | 0.9 |
| protection_tet_m | cft_avi | 0.9 |
| protection_tet_m | mer_vab | 0.9 |
| protection_tet_m | azt_avi | 0.9 |
| protection_tet_m | carb_group1 | 0.9 |
| protection_tet_m | carb_group2 | 0.9 |
| protection_tet_m | mono | 0.9 |
| protection_tet_m | fq | 0.9 |
| protection_tet_m | ag_group1 | 0.9 |
| protection_tet_m | ag_group2 | 0.9 |
| protection_tet_m | mls | 0.9 |
| protection_tet_m | lincosamides | 0.9 |
| protection_tet_m | glyc | 0.9 |
| protection_tet_m | lipoglycopeptides | 0.9 |
| protection_tet_m | tet | 0.9 |
| protection_tet_m | glycylcyclines | 0.9 |
| protection_tet_m | poly | 0.9 |
| protection_tet_m | oxa | 0.9 |
| protection_tet_m | chl | 0.9 |
| protection_tet_m | sulf | 0.9 |
| protection_tet_m | lipopeptides | 0.9 |
| protection_tet_m | streptogramins | 0.9 |
| protection_tet_m | nitrofurans | 0.9 |
| protection_tet_m | phosphonic_acids | 0.9 |
| protection_tet_m | nitroimidazoles | 0.9 |
| protection_tet_m | rifamycins | 0.9 |
| protection_tet_m | macrocycles | 0.9 |
| protection_tet_m | steroid_antibacterials | 0.9 |
| protection_tet_m | pleuromutilins | 0.9 |
| protection_tet_m | other | 0.9 |
| enzyme_aac_aph | pen | 0.85 |
| enzyme_aac_aph | bli | 0.85 |
| enzyme_aac_aph | bli_anti_pseudomonal | 0.85 |
| enzyme_aac_aph | bli_sulbactam | 0.85 |
| enzyme_aac_aph | c1_2g | 0.85 |
| enzyme_aac_aph | c3g | 0.85 |
| enzyme_aac_aph | c3g_bli | 0.85 |
| enzyme_aac_aph | c4g | 0.85 |
| enzyme_aac_aph | anti_mrsa_ceph | 0.85 |
| enzyme_aac_aph | siderophore_ceph | 0.85 |
| enzyme_aac_aph | cft_avi | 0.85 |
| enzyme_aac_aph | mer_vab | 0.85 |
| enzyme_aac_aph | azt_avi | 0.85 |
| enzyme_aac_aph | carb_group1 | 0.85 |
| enzyme_aac_aph | carb_group2 | 0.85 |
| enzyme_aac_aph | mono | 0.85 |
| enzyme_aac_aph | fq | 0.85 |
| enzyme_aac_aph | ag_group1 | 0.85 |
| enzyme_aac_aph | ag_group2 | 0.85 |
| enzyme_aac_aph | mls | 0.85 |
| enzyme_aac_aph | lincosamides | 0.85 |
| enzyme_aac_aph | glyc | 0.85 |
| enzyme_aac_aph | lipoglycopeptides | 0.85 |
| enzyme_aac_aph | tet | 0.85 |
| enzyme_aac_aph | glycylcyclines | 0.85 |
| enzyme_aac_aph | poly | 0.85 |
| enzyme_aac_aph | oxa | 0.85 |
| enzyme_aac_aph | chl | 0.85 |
| enzyme_aac_aph | sulf | 0.85 |
| enzyme_aac_aph | lipopeptides | 0.85 |
| enzyme_aac_aph | streptogramins | 0.85 |
| enzyme_aac_aph | nitrofurans | 0.85 |
| enzyme_aac_aph | phosphonic_acids | 0.85 |
| enzyme_aac_aph | nitroimidazoles | 0.85 |
| enzyme_aac_aph | rifamycins | 0.85 |
| enzyme_aac_aph | macrocycles | 0.85 |
| enzyme_aac_aph | steroid_antibacterials | 0.85 |
| enzyme_aac_aph | pleuromutilins | 0.85 |
| enzyme_aac_aph | other | 0.85 |
| enzyme_bla_z | pen | 0.9 |
| enzyme_bla_z | bli | 0.9 |
| enzyme_bla_z | bli_anti_pseudomonal | 0.9 |
| enzyme_bla_z | bli_sulbactam | 0.9 |
| enzyme_bla_z | c1_2g | 0.9 |
| enzyme_bla_z | c3g | 0.9 |
| enzyme_bla_z | c3g_bli | 0.9 |
| enzyme_bla_z | c4g | 0.9 |
| enzyme_bla_z | anti_mrsa_ceph | 0.9 |
| enzyme_bla_z | siderophore_ceph | 0.9 |
| enzyme_bla_z | cft_avi | 0.9 |
| enzyme_bla_z | mer_vab | 0.9 |
| enzyme_bla_z | azt_avi | 0.9 |
| enzyme_bla_z | carb_group1 | 0.9 |
| enzyme_bla_z | carb_group2 | 0.9 |
| enzyme_bla_z | mono | 0.9 |
| enzyme_bla_z | fq | 0.9 |
| enzyme_bla_z | ag_group1 | 0.9 |
| enzyme_bla_z | ag_group2 | 0.9 |
| enzyme_bla_z | mls | 0.9 |
| enzyme_bla_z | lincosamides | 0.9 |
| enzyme_bla_z | glyc | 0.9 |
| enzyme_bla_z | lipoglycopeptides | 0.9 |
| enzyme_bla_z | tet | 0.9 |
| enzyme_bla_z | glycylcyclines | 0.9 |
| enzyme_bla_z | poly | 0.9 |
| enzyme_bla_z | oxa | 0.9 |
| enzyme_bla_z | chl | 0.9 |
| enzyme_bla_z | sulf | 0.9 |
| enzyme_bla_z | lipopeptides | 0.9 |
| enzyme_bla_z | streptogramins | 0.9 |
| enzyme_bla_z | nitrofurans | 0.9 |
| enzyme_bla_z | phosphonic_acids | 0.9 |
| enzyme_bla_z | nitroimidazoles | 0.9 |
| enzyme_bla_z | rifamycins | 0.9 |
| enzyme_bla_z | macrocycles | 0.9 |
| enzyme_bla_z | steroid_antibacterials | 0.9 |
| enzyme_bla_z | pleuromutilins | 0.9 |
| enzyme_bla_z | other | 0.9 |
| enzyme_oxa_acinetobacter | pen | 0.8 |
| enzyme_oxa_acinetobacter | bli | 0.8 |
| enzyme_oxa_acinetobacter | bli_anti_pseudomonal | 0.8 |
| enzyme_oxa_acinetobacter | bli_sulbactam | 0.8 |
| enzyme_oxa_acinetobacter | c1_2g | 0.8 |
| enzyme_oxa_acinetobacter | c3g | 0.8 |
| enzyme_oxa_acinetobacter | c3g_bli | 0.8 |
| enzyme_oxa_acinetobacter | c4g | 0.8 |
| enzyme_oxa_acinetobacter | anti_mrsa_ceph | 0.8 |
| enzyme_oxa_acinetobacter | siderophore_ceph | 0.8 |
| enzyme_oxa_acinetobacter | cft_avi | 0.8 |
| enzyme_oxa_acinetobacter | mer_vab | 0.8 |
| enzyme_oxa_acinetobacter | azt_avi | 0.8 |
| enzyme_oxa_acinetobacter | carb_group1 | 0.8 |
| enzyme_oxa_acinetobacter | carb_group2 | 0.8 |
| enzyme_oxa_acinetobacter | mono | 0.8 |
| enzyme_oxa_acinetobacter | fq | 0.8 |
| enzyme_oxa_acinetobacter | ag_group1 | 0.8 |
| enzyme_oxa_acinetobacter | ag_group2 | 0.8 |
| enzyme_oxa_acinetobacter | mls | 0.8 |
| enzyme_oxa_acinetobacter | lincosamides | 0.8 |
| enzyme_oxa_acinetobacter | glyc | 0.8 |
| enzyme_oxa_acinetobacter | lipoglycopeptides | 0.8 |
| enzyme_oxa_acinetobacter | tet | 0.8 |
| enzyme_oxa_acinetobacter | glycylcyclines | 0.8 |
| enzyme_oxa_acinetobacter | poly | 0.8 |
| enzyme_oxa_acinetobacter | oxa | 0.8 |
| enzyme_oxa_acinetobacter | chl | 0.8 |
| enzyme_oxa_acinetobacter | sulf | 0.8 |
| enzyme_oxa_acinetobacter | lipopeptides | 0.8 |
| enzyme_oxa_acinetobacter | streptogramins | 0.8 |
| enzyme_oxa_acinetobacter | nitrofurans | 0.8 |
| enzyme_oxa_acinetobacter | phosphonic_acids | 0.8 |
| enzyme_oxa_acinetobacter | nitroimidazoles | 0.8 |
| enzyme_oxa_acinetobacter | rifamycins | 0.8 |
| enzyme_oxa_acinetobacter | macrocycles | 0.8 |
| enzyme_oxa_acinetobacter | steroid_antibacterials | 0.8 |
| enzyme_oxa_acinetobacter | pleuromutilins | 0.8 |
| enzyme_oxa_acinetobacter | other | 0.8 |
| mutation_23s_rrna | pen | 0.8 |
| mutation_23s_rrna | bli | 0.8 |
| mutation_23s_rrna | bli_anti_pseudomonal | 0.8 |
| mutation_23s_rrna | bli_sulbactam | 0.8 |
| mutation_23s_rrna | c1_2g | 0.8 |
| mutation_23s_rrna | c3g | 0.8 |
| mutation_23s_rrna | c3g_bli | 0.8 |
| mutation_23s_rrna | c4g | 0.8 |
| mutation_23s_rrna | anti_mrsa_ceph | 0.8 |
| mutation_23s_rrna | siderophore_ceph | 0.8 |
| mutation_23s_rrna | cft_avi | 0.8 |
| mutation_23s_rrna | mer_vab | 0.8 |
| mutation_23s_rrna | azt_avi | 0.8 |
| mutation_23s_rrna | carb_group1 | 0.8 |
| mutation_23s_rrna | carb_group2 | 0.8 |
| mutation_23s_rrna | mono | 0.8 |
| mutation_23s_rrna | fq | 0.8 |
| mutation_23s_rrna | ag_group1 | 0.8 |
| mutation_23s_rrna | ag_group2 | 0.8 |
| mutation_23s_rrna | mls | 0.8 |
| mutation_23s_rrna | lincosamides | 0.8 |
| mutation_23s_rrna | glyc | 0.8 |
| mutation_23s_rrna | lipoglycopeptides | 0.8 |
| mutation_23s_rrna | tet | 0.8 |
| mutation_23s_rrna | glycylcyclines | 0.8 |
| mutation_23s_rrna | poly | 0.8 |
| mutation_23s_rrna | oxa | 0.8 |
| mutation_23s_rrna | chl | 0.8 |
| mutation_23s_rrna | sulf | 0.8 |
| mutation_23s_rrna | lipopeptides | 0.8 |
| mutation_23s_rrna | streptogramins | 0.8 |
| mutation_23s_rrna | nitrofurans | 0.8 |
| mutation_23s_rrna | phosphonic_acids | 0.8 |
| mutation_23s_rrna | nitroimidazoles | 0.8 |
| mutation_23s_rrna | rifamycins | 0.8 |
| mutation_23s_rrna | macrocycles | 0.8 |
| mutation_23s_rrna | steroid_antibacterials | 0.8 |
| mutation_23s_rrna | pleuromutilins | 0.8 |
| mutation_23s_rrna | other | 0.8 |
| efflux_tet_abc | pen | 0.7 |
| efflux_tet_abc | bli | 0.7 |
| efflux_tet_abc | bli_anti_pseudomonal | 0.7 |
| efflux_tet_abc | bli_sulbactam | 0.7 |
| efflux_tet_abc | c1_2g | 0.7 |
| efflux_tet_abc | c3g | 0.7 |
| efflux_tet_abc | c3g_bli | 0.7 |
| efflux_tet_abc | c4g | 0.7 |
| efflux_tet_abc | anti_mrsa_ceph | 0.7 |
| efflux_tet_abc | siderophore_ceph | 0.7 |
| efflux_tet_abc | cft_avi | 0.7 |
| efflux_tet_abc | mer_vab | 0.7 |
| efflux_tet_abc | azt_avi | 0.7 |
| efflux_tet_abc | carb_group1 | 0.7 |
| efflux_tet_abc | carb_group2 | 0.7 |
| efflux_tet_abc | mono | 0.7 |
| efflux_tet_abc | fq | 0.7 |
| efflux_tet_abc | ag_group1 | 0.7 |
| efflux_tet_abc | ag_group2 | 0.7 |
| efflux_tet_abc | mls | 0.7 |
| efflux_tet_abc | lincosamides | 0.7 |
| efflux_tet_abc | glyc | 0.7 |
| efflux_tet_abc | lipoglycopeptides | 0.7 |
| efflux_tet_abc | tet | 0.7 |
| efflux_tet_abc | glycylcyclines | 0.7 |
| efflux_tet_abc | poly | 0.7 |
| efflux_tet_abc | oxa | 0.7 |
| efflux_tet_abc | chl | 0.7 |
| efflux_tet_abc | sulf | 0.7 |
| efflux_tet_abc | lipopeptides | 0.7 |
| efflux_tet_abc | streptogramins | 0.7 |
| efflux_tet_abc | nitrofurans | 0.7 |
| efflux_tet_abc | phosphonic_acids | 0.7 |
| efflux_tet_abc | nitroimidazoles | 0.7 |
| efflux_tet_abc | rifamycins | 0.7 |
| efflux_tet_abc | macrocycles | 0.7 |
| efflux_tet_abc | steroid_antibacterials | 0.7 |
| efflux_tet_abc | pleuromutilins | 0.7 |
| efflux_tet_abc | other | 0.7 |
| mutation_pbp_mosaic | pen | 0.8 |
| mutation_pbp_mosaic | bli | 0.7 |
| mutation_pbp_mosaic | bli_anti_pseudomonal | 0.7 |
| mutation_pbp_mosaic | bli_sulbactam | 0.7 |
| mutation_pbp_mosaic | c1_2g | 0.6 |
| mutation_pbp_mosaic | c3g | 0.3 |
| mutation_pbp_mosaic | c3g_bli | 0.3 |
| mutation_pbp_mosaic | c4g | 0.15 |
| mutation_pbp_mosaic | anti_mrsa_ceph | 0.15 |
| mutation_pbp_mosaic | siderophore_ceph | 0.15 |
| mutation_pbp_mosaic | cft_avi | 0.1 |
| mutation_pbp_mosaic | azt_avi | 0.5 |
| mutation_pbp_mosaic | mono | 0.5 |
| mutation_pbp_mosaic | lipopeptides | 0.5 |
| mutation_pbp_mosaic | streptogramins | 0.5 |
| mutation_pbp_mosaic | nitrofurans | 0.5 |
| mutation_pbp_mosaic | phosphonic_acids | 0.5 |
| mutation_pbp_mosaic | nitroimidazoles | 0.5 |
| mutation_pbp_mosaic | rifamycins | 0.5 |
| mutation_pbp_mosaic | macrocycles | 0.5 |
| mutation_pbp_mosaic | steroid_antibacterials | 0.5 |
| mutation_pbp_mosaic | pleuromutilins | 0.5 |
| efflux_mtr_cde | pen | 0.3 |
| efflux_mtr_cde | mls | 0.5 |
| efflux_mtr_cde | tet | 0.4 |
| efflux_mtr_cde | chl | 0.4 |
| efflux_mtr_cde | lipopeptides | 0.4 |
| efflux_mtr_cde | streptogramins | 0.4 |
| efflux_mtr_cde | nitrofurans | 0.4 |
| efflux_mtr_cde | phosphonic_acids | 0.4 |
| efflux_mtr_cde | nitroimidazoles | 0.4 |
| efflux_mtr_cde | rifamycins | 0.4 |
| efflux_mtr_cde | macrocycles | 0.4 |
| efflux_mtr_cde | steroid_antibacterials | 0.4 |
| efflux_mtr_cde | pleuromutilins | 0.4 |
| as_yet_unknown | pen | 0.5 |
| as_yet_unknown | bli | 0.5 |
| as_yet_unknown | bli_anti_pseudomonal | 0.5 |
| as_yet_unknown | bli_sulbactam | 0.5 |
| as_yet_unknown | c1_2g | 0.5 |
| as_yet_unknown | c3g | 0.5 |
| as_yet_unknown | c3g_bli | 0.5 |
| as_yet_unknown | c4g | 0.5 |
| as_yet_unknown | anti_mrsa_ceph | 0.5 |
| as_yet_unknown | siderophore_ceph | 0.5 |
| as_yet_unknown | cft_avi | 0.5 |
| as_yet_unknown | mer_vab | 0.5 |
| as_yet_unknown | azt_avi | 0.5 |
| as_yet_unknown | carb_group1 | 0.5 |
| as_yet_unknown | carb_group2 | 0.5 |
| as_yet_unknown | mono | 0.5 |
| as_yet_unknown | fq | 0.5 |
| as_yet_unknown | ag_group1 | 0.5 |
| as_yet_unknown | ag_group2 | 0.5 |
| as_yet_unknown | mls | 0.5 |
| as_yet_unknown | lincosamides | 0.5 |
| as_yet_unknown | glyc | 0.5 |
| as_yet_unknown | lipoglycopeptides | 0.5 |
| as_yet_unknown | tet | 0.5 |
| as_yet_unknown | glycylcyclines | 0.5 |
| as_yet_unknown | poly | 0.5 |
| as_yet_unknown | oxa | 0.5 |
| as_yet_unknown | chl | 0.5 |
| as_yet_unknown | sulf | 0.5 |
| as_yet_unknown | lipopeptides | 0.5 |
| as_yet_unknown | streptogramins | 0.5 |
| as_yet_unknown | nitrofurans | 0.5 |
| as_yet_unknown | phosphonic_acids | 0.5 |
| as_yet_unknown | nitroimidazoles | 0.5 |
| as_yet_unknown | rifamycins | 0.5 |
| as_yet_unknown | macrocycles | 0.5 |
| as_yet_unknown | steroid_antibacterials | 0.5 |
| as_yet_unknown | pleuromutilins | 0.5 |
| as_yet_unknown | other | 0.5 |

#### Bacteria–Mechanism Emergence Rates

De novo emergence rate per day for each bacteria–mechanism pair. Only non-zero entries shown.

| Bacteria | Mechanism | Emergence rate/day |
| --- | ---: | ---: |
| acinetobacter_baumannii | enzyme_esbl_ctx_m | 0.006 |
| acinetobacter_baumannii | enzyme_esbl_tem | 0.006 |
| acinetobacter_baumannii | enzyme_esbl_shv | 0.006 |
| acinetobacter_baumannii | enzyme_kpc | 0.006 |
| acinetobacter_baumannii | enzyme_ndm_vim | 0.006 |
| acinetobacter_baumannii | enzyme_oxa_48 | 0.006 |
| acinetobacter_baumannii | enzyme_ampc_cmy | 0.006 |
| acinetobacter_baumannii | enzyme_ampc_dha | 0.006 |
| acinetobacter_baumannii | mutation_gyra_primary | 0.09 |
| acinetobacter_baumannii | mutation_gyra_parc_secondary | 0.09 |
| acinetobacter_baumannii | protection_qnr | 0.09 |
| acinetobacter_baumannii | enzyme_16s_rrmt | 10 |
| acinetobacter_baumannii | enzyme_cat | 0.01 |
| acinetobacter_baumannii | modification_mcr_1 | 0.01 |
| acinetobacter_baumannii | global_efflux_pump | 0.025 |
| acinetobacter_baumannii | global_porin_loss | 1e-4 |
| acinetobacter_baumannii | mutation_folate_pathway | 0.05 |
| acinetobacter_baumannii | enzyme_fos_a | 0.0035 |
| acinetobacter_baumannii | mutation_rpo_b | 0.045 |
| acinetobacter_baumannii | protection_tet_m | 0.05 |
| acinetobacter_baumannii | enzyme_aac_aph | 10 |
| acinetobacter_baumannii | enzyme_oxa_acinetobacter | 0.006 |
| acinetobacter_baumannii | efflux_tet_abc | 0.0025 |
| acinetobacter_baumannii | mutation_pbp_mosaic | 0.006 |
| citrobacter_spp. | enzyme_esbl_ctx_m | 0.001 |
| citrobacter_spp. | enzyme_esbl_tem | 0.001 |
| citrobacter_spp. | enzyme_esbl_shv | 0.001 |
| citrobacter_spp. | enzyme_kpc | 0.001 |
| citrobacter_spp. | enzyme_ndm_vim | 0.001 |
| citrobacter_spp. | enzyme_oxa_48 | 0.001 |
| citrobacter_spp. | enzyme_ampc_cmy | 0.001 |
| citrobacter_spp. | enzyme_ampc_dha | 0.001 |
| citrobacter_spp. | mutation_gyra_primary | 0.04 |
| citrobacter_spp. | mutation_gyra_parc_secondary | 3e-4 |
| citrobacter_spp. | protection_qnr | 0.04 |
| citrobacter_spp. | enzyme_16s_rrmt | 0.5 |
| citrobacter_spp. | enzyme_cat | 0.03 |
| citrobacter_spp. | efflux_acrab_tolc | 0.015 |
| citrobacter_spp. | efflux_mexxy_oprm | 0.0015 |
| citrobacter_spp. | modification_mcr_1 | 0.05 |
| citrobacter_spp. | global_efflux_pump | 0.003 |
| citrobacter_spp. | global_porin_loss | 5e-4 |
| citrobacter_spp. | mutation_folate_pathway | 0.003 |
| citrobacter_spp. | mutation_nitroreductase | 0.015 |
| citrobacter_spp. | enzyme_fos_a | 0.005 |
| citrobacter_spp. | mutation_rpo_b | 0.02 |
| citrobacter_spp. | protection_tet_m | 0.004 |
| citrobacter_spp. | enzyme_aac_aph | 0.5 |
| citrobacter_spp. | efflux_tet_abc | 0.004 |
| citrobacter_spp. | mutation_pbp_mosaic | 0.0015 |
| enterobacter_spp. | enzyme_esbl_ctx_m | 0.001 |
| enterobacter_spp. | enzyme_esbl_tem | 0.001 |
| enterobacter_spp. | enzyme_esbl_shv | 0.001 |
| enterobacter_spp. | enzyme_kpc | 3e-4 |
| enterobacter_spp. | enzyme_ndm_vim | 3e-4 |
| enterobacter_spp. | enzyme_oxa_48 | 3e-4 |
| enterobacter_spp. | enzyme_ampc_cmy | 0.001 |
| enterobacter_spp. | enzyme_ampc_dha | 0.001 |
| enterobacter_spp. | mutation_gyra_primary | 5e-4 |
| enterobacter_spp. | mutation_gyra_parc_secondary | 5e-4 |
| enterobacter_spp. | protection_qnr | 5e-4 |
| enterobacter_spp. | enzyme_16s_rrmt | 0.3 |
| enterobacter_spp. | enzyme_cat | 0.004 |
| enterobacter_spp. | efflux_acrab_tolc | 5e-4 |
| enterobacter_spp. | modification_mcr_1 | 0.03 |
| enterobacter_spp. | global_efflux_pump | 5e-4 |
| enterobacter_spp. | global_porin_loss | 5e-4 |
| enterobacter_spp. | mutation_folate_pathway | 0.001 |
| enterobacter_spp. | mutation_nitroreductase | 0.005 |
| enterobacter_spp. | enzyme_fos_a | 0.02 |
| enterobacter_spp. | mutation_rpo_b | 0.03 |
| enterobacter_spp. | protection_tet_m | 0.01 |
| enterobacter_spp. | enzyme_aac_aph | 0.3 |
| enterobacter_spp. | efflux_tet_abc | 0.01 |
| enterobacter_spp. | mutation_pbp_mosaic | 0.001 |
| enterococcus_faecalis | target_site_pbp2a_meca | 2e-5 |
| enterococcus_faecalis | target_site_van_a | 0.02 |
| enterococcus_faecalis | target_site_van_b | 0.02 |
| enterococcus_faecalis | mutation_gyra_primary | 0.002 |
| enterococcus_faecalis | mutation_gyra_parc_secondary | 3e-4 |
| enterococcus_faecalis | target_site_erm_b | 0.002 |
| enterococcus_faecalis | target_site_cfr | 0.003 |
| enterococcus_faecalis | enzyme_cat | 2e-4 |
| enterococcus_faecalis | global_efflux_pump | 3e-4 |
| enterococcus_faecalis | global_porin_loss | 1e-4 |
| enterococcus_faecalis | mutation_folate_pathway | 0.02 |
| enterococcus_faecalis | mutation_nitroreductase | 0.02 |
| enterococcus_faecalis | mutation_mpr_f | 0.005 |
| enterococcus_faecalis | mutation_rpo_b | 0.005 |
| enterococcus_faecalis | protection_fus_b | 5e-4 |
| enterococcus_faecalis | protection_tet_m | 0.002 |
| enterococcus_faecalis | enzyme_aac_aph | 2e-5 |
| enterococcus_faecalis | mutation_23s_rrna | 3e-4 |
| enterococcus_faecalis | mutation_pbp_mosaic | 1.5e-5 |
| enterococcus_faecium | target_site_van_a | 0.008 |
| enterococcus_faecium | target_site_van_b | 0.008 |
| enterococcus_faecium | mutation_gyra_primary | 5e-4 |
| enterococcus_faecium | mutation_gyra_parc_secondary | 0.002 |
| enterococcus_faecium | enzyme_16s_rrmt | 0.003 |
| enterococcus_faecium | target_site_erm_b | 0.008 |
| enterococcus_faecium | target_site_cfr | 0.008 |
| enterococcus_faecium | enzyme_cat | 0.0015 |
| enterococcus_faecium | global_efflux_pump | 0.005 |
| enterococcus_faecium | mutation_folate_pathway | 0.004 |
| enterococcus_faecium | mutation_nitroreductase | 0.3 |
| enterococcus_faecium | enzyme_fos_a | 0.3 |
| enterococcus_faecium | mutation_mpr_f | 0.012 |
| enterococcus_faecium | mutation_rpo_b | 0.12 |
| enterococcus_faecium | protection_fus_b | 0.012 |
| enterococcus_faecium | protection_tet_m | 0.004 |
| enterococcus_faecium | enzyme_aac_aph | 0.025 |
| enterococcus_faecium | mutation_23s_rrna | 0.0012 |
| enterococcus_faecium | mutation_pbp_mosaic | 0.012 |
| enterococcus_faecium | efflux_mtr_cde | 0.0012 |
| escherichia_coli | enzyme_esbl_ctx_m | 1e-6 |
| escherichia_coli | enzyme_esbl_tem | 1e-6 |
| escherichia_coli | enzyme_esbl_shv | 1e-6 |
| escherichia_coli | enzyme_kpc | 1e-6 |
| escherichia_coli | enzyme_ndm_vim | 1e-6 |
| escherichia_coli | enzyme_oxa_48 | 1e-6 |
| escherichia_coli | enzyme_ampc_cmy | 1e-6 |
| escherichia_coli | enzyme_ampc_dha | 1e-6 |
| escherichia_coli | mutation_gyra_primary | 0.003 |
| escherichia_coli | mutation_gyra_parc_secondary | 3e-4 |
| escherichia_coli | protection_qnr | 0.003 |
| escherichia_coli | enzyme_16s_rrmt | 1e-6 |
| escherichia_coli | enzyme_cat | 1e-6 |
| escherichia_coli | efflux_acrab_tolc | 1e-4 |
| escherichia_coli | modification_mcr_1 | 1e-5 |
| escherichia_coli | global_efflux_pump | 1e-4 |
| escherichia_coli | global_porin_loss | 1e-7 |
| escherichia_coli | mutation_folate_pathway | 1e-4 |
| escherichia_coli | mutation_nitroreductase | 1e-6 |
| escherichia_coli | enzyme_fos_a | 1e-6 |
| escherichia_coli | mutation_rpo_b | 1e-5 |
| escherichia_coli | protection_tet_m | 1e-6 |
| escherichia_coli | enzyme_aac_aph | 0.001 |
| escherichia_coli | efflux_tet_abc | 1e-6 |
| klebsiella_pneumoniae | enzyme_esbl_ctx_m | 1e-5 |
| klebsiella_pneumoniae | enzyme_esbl_tem | 1e-5 |
| klebsiella_pneumoniae | enzyme_esbl_shv | 1e-5 |
| klebsiella_pneumoniae | enzyme_kpc | 1e-5 |
| klebsiella_pneumoniae | enzyme_ndm_vim | 1e-5 |
| klebsiella_pneumoniae | enzyme_oxa_48 | 1e-5 |
| klebsiella_pneumoniae | enzyme_ampc_cmy | 1e-5 |
| klebsiella_pneumoniae | enzyme_ampc_dha | 1e-5 |
| klebsiella_pneumoniae | mutation_gyra_primary | 0.001 |
| klebsiella_pneumoniae | mutation_gyra_parc_secondary | 0.001 |
| klebsiella_pneumoniae | protection_qnr | 0.001 |
| klebsiella_pneumoniae | enzyme_16s_rrmt | 0.015 |
| klebsiella_pneumoniae | enzyme_cat | 0.002 |
| klebsiella_pneumoniae | efflux_acrab_tolc | 0.001 |
| klebsiella_pneumoniae | porin_loss_ompk35_36 | 1e-5 |
| klebsiella_pneumoniae | modification_mcr_1 | 0.01 |
| klebsiella_pneumoniae | global_efflux_pump | 0.001 |
| klebsiella_pneumoniae | global_porin_loss | 3e-9 |
| klebsiella_pneumoniae | mutation_folate_pathway | 0.001 |
| klebsiella_pneumoniae | mutation_nitroreductase | 0.002 |
| klebsiella_pneumoniae | enzyme_fos_a | 0.003 |
| klebsiella_pneumoniae | mutation_rpo_b | 2e-4 |
| klebsiella_pneumoniae | protection_tet_m | 2e-4 |
| klebsiella_pneumoniae | enzyme_aac_aph | 0.04 |
| klebsiella_pneumoniae | efflux_tet_abc | 3e-4 |
| morganella_spp. | enzyme_esbl_ctx_m | 2e-4 |
| morganella_spp. | enzyme_esbl_tem | 2e-4 |
| morganella_spp. | enzyme_esbl_shv | 2e-4 |
| morganella_spp. | enzyme_kpc | 5e-5 |
| morganella_spp. | enzyme_ndm_vim | 5e-5 |
| morganella_spp. | enzyme_oxa_48 | 5e-5 |
| morganella_spp. | enzyme_ampc_cmy | 1e-4 |
| morganella_spp. | enzyme_ampc_dha | 1e-4 |
| morganella_spp. | mutation_gyra_primary | 0.03 |
| morganella_spp. | mutation_gyra_parc_secondary | 0.03 |
| morganella_spp. | protection_qnr | 0.03 |
| morganella_spp. | enzyme_16s_rrmt | 1 |
| morganella_spp. | enzyme_cat | 0.1 |
| morganella_spp. | efflux_acrab_tolc | 0.003 |
| morganella_spp. | efflux_mexxy_oprm | 2e-4 |
| morganella_spp. | modification_mcr_1 | 0.02 |
| morganella_spp. | global_efflux_pump | 0.003 |
| morganella_spp. | global_porin_loss | 1e-5 |
| morganella_spp. | mutation_folate_pathway | 0.004 |
| morganella_spp. | mutation_nitroreductase | 0.03 |
| morganella_spp. | enzyme_fos_a | 0.003 |
| morganella_spp. | mutation_rpo_b | 0.01 |
| morganella_spp. | protection_tet_m | 0.003 |
| morganella_spp. | enzyme_aac_aph | 1 |
| morganella_spp. | efflux_tet_abc | 0.003 |
| morganella_spp. | mutation_pbp_mosaic | 2e-5 |
| proteus_spp. | enzyme_esbl_ctx_m | 2e-5 |
| proteus_spp. | enzyme_esbl_tem | 2e-5 |
| proteus_spp. | enzyme_esbl_shv | 2e-5 |
| proteus_spp. | enzyme_kpc | 1.5e-5 |
| proteus_spp. | enzyme_ndm_vim | 1.5e-5 |
| proteus_spp. | enzyme_oxa_48 | 1.5e-5 |
| proteus_spp. | enzyme_ampc_cmy | 2e-5 |
| proteus_spp. | enzyme_ampc_dha | 2e-5 |
| proteus_spp. | mutation_gyra_primary | 0.04 |
| proteus_spp. | mutation_gyra_parc_secondary | 0.01 |
| proteus_spp. | protection_qnr | 0.01 |
| proteus_spp. | enzyme_16s_rrmt | 0.02 |
| proteus_spp. | enzyme_cat | 3.5e-4 |
| proteus_spp. | efflux_acrab_tolc | 0.005 |
| proteus_spp. | modification_mcr_1 | 5e-4 |
| proteus_spp. | global_efflux_pump | 0.05 |
| proteus_spp. | global_porin_loss | 1.5e-4 |
| proteus_spp. | mutation_folate_pathway | 0.002 |
| proteus_spp. | mutation_nitroreductase | 1.5e-5 |
| proteus_spp. | enzyme_fos_a | 0.0015 |
| proteus_spp. | mutation_rpo_b | 1e-4 |
| proteus_spp. | protection_tet_m | 0.0035 |
| proteus_spp. | enzyme_aac_aph | 0.02 |
| proteus_spp. | efflux_tet_abc | 1e-6 |
| serratia_spp. | enzyme_esbl_ctx_m | 0.0018 |
| serratia_spp. | enzyme_esbl_tem | 0.0018 |
| serratia_spp. | enzyme_esbl_shv | 0.0018 |
| serratia_spp. | enzyme_kpc | 4e-4 |
| serratia_spp. | enzyme_ndm_vim | 4e-4 |
| serratia_spp. | enzyme_oxa_48 | 4e-4 |
| serratia_spp. | enzyme_ampc_cmy | 0.0018 |
| serratia_spp. | enzyme_ampc_dha | 0.0018 |
| serratia_spp. | mutation_gyra_primary | 0.02 |
| serratia_spp. | mutation_gyra_parc_secondary | 0.01 |
| serratia_spp. | protection_qnr | 0.01 |
| serratia_spp. | enzyme_16s_rrmt | 0.9 |
| serratia_spp. | enzyme_cat | 5e-5 |
| serratia_spp. | efflux_acrab_tolc | 0.008 |
| serratia_spp. | modification_mcr_1 | 0.03 |
| serratia_spp. | global_efflux_pump | 0.0015 |
| serratia_spp. | global_porin_loss | 3e-5 |
| serratia_spp. | mutation_folate_pathway | 0.001 |
| serratia_spp. | mutation_nitroreductase | 0.002 |
| serratia_spp. | enzyme_fos_a | 2e-4 |
| serratia_spp. | mutation_rpo_b | 0.03 |
| serratia_spp. | protection_tet_m | 0.0015 |
| serratia_spp. | enzyme_aac_aph | 1 |
| serratia_spp. | efflux_tet_abc | 2e-4 |
| p_stuartii | enzyme_esbl_ctx_m | 1.9e-4 |
| p_stuartii | enzyme_esbl_tem | 1.9e-4 |
| p_stuartii | enzyme_esbl_shv | 3.7e-5 |
| p_stuartii | enzyme_kpc | 1.9e-5 |
| p_stuartii | enzyme_ndm_vim | 1.9e-5 |
| p_stuartii | enzyme_oxa_48 | 1.9e-5 |
| p_stuartii | enzyme_ampc_cmy | 3.8e-4 |
| p_stuartii | enzyme_ampc_dha | 1.9e-4 |
| p_stuartii | mutation_gyra_primary | 3.8e-4 |
| p_stuartii | mutation_gyra_parc_secondary | 1.9e-4 |
| p_stuartii | protection_qnr | 1.9e-4 |
| p_stuartii | enzyme_16s_rrmt | 3.8e-6 |
| p_stuartii | enzyme_cat | 1.9e-4 |
| p_stuartii | efflux_acrab_tolc | 1.9e-4 |
| p_stuartii | modification_mcr_1 | 1.9e-5 |
| p_stuartii | global_efflux_pump | 3.7e-5 |
| p_stuartii | global_porin_loss | 3.7e-5 |
| p_stuartii | mutation_folate_pathway | 3.8e-4 |
| p_stuartii | mutation_nitroreductase | 3.7e-5 |
| p_stuartii | enzyme_fos_a | 3.7e-5 |
| p_stuartii | mutation_rpo_b | 3.8e-6 |
| p_stuartii | protection_tet_m | 3.8e-4 |
| p_stuartii | enzyme_aac_aph | 5e-11 |
| p_stuartii | efflux_tet_abc | 5e-11 |
| pseudomonas_aeruginosa | enzyme_esbl_ctx_m | 3e-5 |
| pseudomonas_aeruginosa | enzyme_esbl_tem | 3e-5 |
| pseudomonas_aeruginosa | enzyme_esbl_shv | 3e-5 |
| pseudomonas_aeruginosa | enzyme_kpc | 3e-5 |
| pseudomonas_aeruginosa | enzyme_ndm_vim | 3e-5 |
| pseudomonas_aeruginosa | enzyme_oxa_48 | 3e-5 |
| pseudomonas_aeruginosa | enzyme_ampc_cmy | 3e-5 |
| pseudomonas_aeruginosa | enzyme_ampc_dha | 3e-5 |
| pseudomonas_aeruginosa | mutation_gyra_primary | 0.001 |
| pseudomonas_aeruginosa | mutation_gyra_parc_secondary | 0.001 |
| pseudomonas_aeruginosa | protection_qnr | 0.001 |
| pseudomonas_aeruginosa | enzyme_16s_rrmt | 4e-4 |
| pseudomonas_aeruginosa | target_site_erm_b | 3e-5 |
| pseudomonas_aeruginosa | target_site_cfr | 3e-5 |
| pseudomonas_aeruginosa | enzyme_cat | 4e-4 |
| pseudomonas_aeruginosa | efflux_mexxy_oprm | 5e-4 |
| pseudomonas_aeruginosa | porin_loss_oprd | 5e-4 |
| pseudomonas_aeruginosa | modification_mcr_1 | 0.003 |
| pseudomonas_aeruginosa | global_efflux_pump | 3e-4 |
| pseudomonas_aeruginosa | global_porin_loss | 3e-4 |
| pseudomonas_aeruginosa | mutation_folate_pathway | 0.003 |
| pseudomonas_aeruginosa | mutation_nitroreductase | 3e-5 |
| pseudomonas_aeruginosa | enzyme_fos_a | 0.01 |
| pseudomonas_aeruginosa | mutation_rpo_b | 0.001 |
| pseudomonas_aeruginosa | protection_tet_m | 0.004 |
| pseudomonas_aeruginosa | enzyme_aac_aph | 0.001 |
| pseudomonas_aeruginosa | efflux_tet_abc | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_esbl_ctx_m | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_esbl_tem | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_esbl_shv | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_kpc | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_ndm_vim | 0.1 |
| stenotrophomonas_maltophilia | enzyme_oxa_48 | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_ampc_cmy | 0.05 |
| stenotrophomonas_maltophilia | enzyme_ampc_dha | 2e-4 |
| stenotrophomonas_maltophilia | mutation_gyra_primary | 5e-8 |
| stenotrophomonas_maltophilia | mutation_gyra_parc_secondary | 5e-8 |
| stenotrophomonas_maltophilia | protection_qnr | 2e-5 |
| stenotrophomonas_maltophilia | enzyme_16s_rrmt | 0.01 |
| stenotrophomonas_maltophilia | target_site_erm_b | 0.1 |
| stenotrophomonas_maltophilia | enzyme_cat | 5e-7 |
| stenotrophomonas_maltophilia | modification_mcr_1 | 0.05 |
| stenotrophomonas_maltophilia | global_efflux_pump | 1e-7 |
| stenotrophomonas_maltophilia | global_porin_loss | 1e-13 |
| stenotrophomonas_maltophilia | mutation_folate_pathway | 0.005 |
| stenotrophomonas_maltophilia | mutation_nitroreductase | 10 |
| stenotrophomonas_maltophilia | enzyme_fos_a | 0.03 |
| stenotrophomonas_maltophilia | mutation_rpo_b | 2e-4 |
| stenotrophomonas_maltophilia | protection_tet_m | 5e-7 |
| stenotrophomonas_maltophilia | enzyme_aac_aph | 0.05 |
| stenotrophomonas_maltophilia | efflux_tet_abc | 2e-9 |
| staphylococcus_aureus | target_site_pbp2a_meca | 1e-7 |
| staphylococcus_aureus | target_site_van_a | 1e-8 |
| staphylococcus_aureus | target_site_van_b | 5e-11 |
| staphylococcus_aureus | mutation_gyra_primary | 1.5e-5 |
| staphylococcus_aureus | mutation_gyra_parc_secondary | 1.5e-5 |
| staphylococcus_aureus | enzyme_16s_rrmt | 1e-5 |
| staphylococcus_aureus | target_site_erm_b | 4e-6 |
| staphylococcus_aureus | target_site_cfr | 4e-8 |
| staphylococcus_aureus | enzyme_cat | 1e-8 |
| staphylococcus_aureus | global_efflux_pump | 4e-6 |
| staphylococcus_aureus | mutation_folate_pathway | 1e-5 |
| staphylococcus_aureus | mutation_nitroreductase | 1e-5 |
| staphylococcus_aureus | enzyme_fos_a | 1e-5 |
| staphylococcus_aureus | mutation_mpr_f | 1e-5 |
| staphylococcus_aureus | mutation_rpo_b | 1e-5 |
| staphylococcus_aureus | protection_fus_b | 1e-5 |
| staphylococcus_aureus | protection_tet_m | 4e-6 |
| staphylococcus_aureus | enzyme_aac_aph | 1e-5 |
| staphylococcus_aureus | enzyme_bla_z | 1e-8 |
| staphylococcus_epidermidis | target_site_pbp2a_meca | 1e-6 |
| staphylococcus_epidermidis | target_site_van_a | 1e-8 |
| staphylococcus_epidermidis | target_site_van_b | 1e-8 |
| staphylococcus_epidermidis | mutation_gyra_primary | 1e-5 |
| staphylococcus_epidermidis | mutation_gyra_parc_secondary | 1e-6 |
| staphylococcus_epidermidis | target_site_erm_b | 1e-5 |
| staphylococcus_epidermidis | target_site_cfr | 1e-7 |
| staphylococcus_epidermidis | enzyme_cat | 1e-6 |
| staphylococcus_epidermidis | global_efflux_pump | 1e-6 |
| staphylococcus_epidermidis | mutation_folate_pathway | 1e-5 |
| staphylococcus_epidermidis | mutation_mpr_f | 1e-6 |
| staphylococcus_epidermidis | mutation_rpo_b | 0.001 |
| staphylococcus_epidermidis | protection_fus_b | 3e-6 |
| staphylococcus_epidermidis | protection_tet_m | 3e-5 |
| staphylococcus_epidermidis | enzyme_aac_aph | 1e-9 |
| staphylococcus_epidermidis | enzyme_bla_z | 1e-7 |
| streptococcus_pneumoniae | target_site_pbp2a_meca | 3e-8 |
| streptococcus_pneumoniae | mutation_gyra_primary | 3e-5 |
| streptococcus_pneumoniae | mutation_gyra_parc_secondary | 3e-5 |
| streptococcus_pneumoniae | target_site_erm_b | 1e-6 |
| streptococcus_pneumoniae | enzyme_cat | 3e-4 |
| streptococcus_pneumoniae | global_efflux_pump | 1e-5 |
| streptococcus_pneumoniae | global_porin_loss | 3e-8 |
| streptococcus_pneumoniae | mutation_folate_pathway | 0.003 |
| streptococcus_pneumoniae | mutation_rpo_b | 0.001 |
| streptococcus_pneumoniae | protection_tet_m | 0.001 |
| streptococcus_pneumoniae | enzyme_bla_z | 3e-8 |
| streptococcus_pneumoniae | mutation_23s_rrna | 5e-7 |
| streptococcus_pneumoniae | mutation_pbp_mosaic | 3e-8 |
| salmonella_enterica_serovar_typhi | enzyme_esbl_ctx_m | 0.005 |
| salmonella_enterica_serovar_typhi | enzyme_esbl_tem | 0.15 |
| salmonella_enterica_serovar_typhi | enzyme_esbl_shv | 0.0015 |
| salmonella_enterica_serovar_typhi | enzyme_kpc | 3e-5 |
| salmonella_enterica_serovar_typhi | enzyme_ndm_vim | 3e-5 |
| salmonella_enterica_serovar_typhi | enzyme_oxa_48 | 3e-5 |
| salmonella_enterica_serovar_typhi | enzyme_ampc_cmy | 5e-4 |
| salmonella_enterica_serovar_typhi | enzyme_ampc_dha | 5e-4 |
| salmonella_enterica_serovar_typhi | mutation_gyra_primary | 0.8 |
| salmonella_enterica_serovar_typhi | mutation_gyra_parc_secondary | 0.05 |
| salmonella_enterica_serovar_typhi | protection_qnr | 0.03 |
| salmonella_enterica_serovar_typhi | enzyme_16s_rrmt | 1e-4 |
| salmonella_enterica_serovar_typhi | enzyme_cat | 0.015 |
| salmonella_enterica_serovar_typhi | efflux_acrab_tolc | 0.05 |
| salmonella_enterica_serovar_typhi | modification_mcr_1 | 0.0015 |
| salmonella_enterica_serovar_typhi | global_efflux_pump | 0.005 |
| salmonella_enterica_serovar_typhi | global_porin_loss | 0.001 |
| salmonella_enterica_serovar_typhi | mutation_folate_pathway | 0.15 |
| salmonella_enterica_serovar_typhi | mutation_nitroreductase | 0.0015 |
| salmonella_enterica_serovar_typhi | enzyme_fos_a | 0.0015 |
| salmonella_enterica_serovar_typhi | mutation_rpo_b | 1.5e-4 |
| salmonella_enterica_serovar_typhi | protection_tet_m | 0.015 |
| salmonella_enterica_serovar_typhi | enzyme_aac_aph | 1e-4 |
| salmonella_enterica_serovar_typhi | efflux_tet_abc | 2e-6 |
| salmonella_enterica_serovar_typhi | efflux_mtr_cde | 0.01 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_esbl_ctx_m | 0.002 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_esbl_tem | 0.002 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_esbl_shv | 0.002 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_kpc | 5e-5 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_ndm_vim | 7e-5 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_oxa_48 | 5e-5 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_ampc_cmy | 0.0015 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_ampc_dha | 0.0015 |
| salmonella_enterica_serovar_paratyphi_a | mutation_gyra_primary | 0.05 |
| salmonella_enterica_serovar_paratyphi_a | mutation_gyra_parc_secondary | 0.04 |
| salmonella_enterica_serovar_paratyphi_a | protection_qnr | 0.005 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_16s_rrmt | 0.15 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_cat | 0.02 |
| salmonella_enterica_serovar_paratyphi_a | efflux_acrab_tolc | 0.03 |
| salmonella_enterica_serovar_paratyphi_a | efflux_mexxy_oprm | 8e-4 |
| salmonella_enterica_serovar_paratyphi_a | modification_mcr_1 | 0.003 |
| salmonella_enterica_serovar_paratyphi_a | global_efflux_pump | 0.0045 |
| salmonella_enterica_serovar_paratyphi_a | global_porin_loss | 1e-4 |
| salmonella_enterica_serovar_paratyphi_a | mutation_folate_pathway | 0.03 |
| salmonella_enterica_serovar_paratyphi_a | mutation_nitroreductase | 1.5e-4 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_fos_a | 1.5e-4 |
| salmonella_enterica_serovar_paratyphi_a | mutation_rpo_b | 0.01 |
| salmonella_enterica_serovar_paratyphi_a | protection_tet_m | 0.02 |
| salmonella_enterica_serovar_paratyphi_a | enzyme_aac_aph | 0.15 |
| salmonella_enterica_serovar_paratyphi_a | efflux_tet_abc | 0.007 |
| salmonella_enterica_serovar_paratyphi_a | mutation_pbp_mosaic | 0.001 |
| salmonella_enterica_serovar_paratyphi_a | efflux_mtr_cde | 5e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_esbl_ctx_m | 2.2e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_esbl_tem | 2.2e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_esbl_shv | 2.2e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_kpc | 2.2e-5 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_ndm_vim | 2.2e-5 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_oxa_48 | 2.2e-5 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_ampc_cmy | 2.2e-4 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_ampc_dha | 2.2e-4 |
| invasive_non-typhoidal_salmonella_spp. | mutation_gyra_primary | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | mutation_gyra_parc_secondary | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | protection_qnr | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_16s_rrmt | 0.8 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_cat | 0.015 |
| invasive_non-typhoidal_salmonella_spp. | efflux_acrab_tolc | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | modification_mcr_1 | 0.03 |
| invasive_non-typhoidal_salmonella_spp. | global_efflux_pump | 0.1 |
| invasive_non-typhoidal_salmonella_spp. | global_porin_loss | 1e-5 |
| invasive_non-typhoidal_salmonella_spp. | mutation_folate_pathway | 0.003 |
| invasive_non-typhoidal_salmonella_spp. | mutation_nitroreductase | 0.003 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_fos_a | 0.001 |
| invasive_non-typhoidal_salmonella_spp. | mutation_rpo_b | 0.02 |
| invasive_non-typhoidal_salmonella_spp. | protection_tet_m | 0.004 |
| invasive_non-typhoidal_salmonella_spp. | enzyme_aac_aph | 0.8 |
| invasive_non-typhoidal_salmonella_spp. | efflux_tet_abc | 0.08 |
| invasive_non-typhoidal_salmonella_spp. | efflux_mtr_cde | 2.5e-6 |
| shigella_spp. | enzyme_esbl_ctx_m | 0.001 |
| shigella_spp. | enzyme_esbl_tem | 0.001 |
| shigella_spp. | enzyme_esbl_shv | 0.001 |
| shigella_spp. | enzyme_kpc | 0.001 |
| shigella_spp. | enzyme_ndm_vim | 0.001 |
| shigella_spp. | enzyme_oxa_48 | 0.001 |
| shigella_spp. | enzyme_ampc_cmy | 0.001 |
| shigella_spp. | enzyme_ampc_dha | 0.001 |
| shigella_spp. | mutation_gyra_primary | 5e-4 |
| shigella_spp. | mutation_gyra_parc_secondary | 5e-4 |
| shigella_spp. | protection_qnr | 5e-4 |
| shigella_spp. | enzyme_16s_rrmt | 0.9 |
| shigella_spp. | target_site_erm_b | 0.8 |
| shigella_spp. | enzyme_cat | 2.5e-4 |
| shigella_spp. | efflux_acrab_tolc | 5e-4 |
| shigella_spp. | modification_mcr_1 | 3e-4 |
| shigella_spp. | global_efflux_pump | 0.9 |
| shigella_spp. | global_porin_loss | 3e-5 |
| shigella_spp. | mutation_folate_pathway | 0.003 |
| shigella_spp. | mutation_rpo_b | 0.04 |
| shigella_spp. | protection_tet_m | 0.03 |
| shigella_spp. | enzyme_aac_aph | 0.9 |
| shigella_spp. | mutation_23s_rrna | 0.9 |
| shigella_spp. | efflux_tet_abc | 0.03 |
| shigella_spp. | mutation_pbp_mosaic | 0.001 |
| shigella_spp. | efflux_mtr_cde | 0.001 |
| neisseria_gonorrhoeae | mutation_gyra_primary | 1 |
| neisseria_gonorrhoeae | mutation_gyra_parc_secondary | 1 |
| neisseria_gonorrhoeae | protection_qnr | 1 |
| neisseria_gonorrhoeae | enzyme_16s_rrmt | 1 |
| neisseria_gonorrhoeae | target_site_erm_b | 0.01 |
| neisseria_gonorrhoeae | target_site_cfr | 0.01 |
| neisseria_gonorrhoeae | enzyme_cat | 0.001 |
| neisseria_gonorrhoeae | modification_mcr_1 | 0.005 |
| neisseria_gonorrhoeae | global_efflux_pump | 1 |
| neisseria_gonorrhoeae | global_porin_loss | 0.001 |
| neisseria_gonorrhoeae | mutation_folate_pathway | 0.1 |
| neisseria_gonorrhoeae | mutation_nitroreductase | 0.1 |
| neisseria_gonorrhoeae | enzyme_fos_a | 3e-4 |
| neisseria_gonorrhoeae | mutation_rpo_b | 1 |
| neisseria_gonorrhoeae | protection_tet_m | 1 |
| neisseria_gonorrhoeae | enzyme_aac_aph | 1 |
| neisseria_gonorrhoeae | enzyme_bla_z | 0.02 |
| neisseria_gonorrhoeae | mutation_23s_rrna | 0.03 |
| neisseria_gonorrhoeae | efflux_tet_abc | 0.003 |
| neisseria_gonorrhoeae | mutation_pbp_mosaic | 0.02 |
| neisseria_gonorrhoeae | efflux_mtr_cde | 0.02 |
| streptococcus_pyogenes | mutation_gyra_primary | 3.8e-7 |
| streptococcus_pyogenes | mutation_gyra_parc_secondary | 7.5e-8 |
| streptococcus_pyogenes | target_site_erm_b | 7.5e-6 |
| streptococcus_pyogenes | target_site_cfr | 7.5e-9 |
| streptococcus_pyogenes | enzyme_cat | 3.8e-7 |
| streptococcus_pyogenes | global_efflux_pump | 7.5e-7 |
| streptococcus_pyogenes | mutation_folate_pathway | 3.8e-6 |
| streptococcus_pyogenes | mutation_mpr_f | 7.5e-9 |
| streptococcus_pyogenes | mutation_rpo_b | 7.5e-9 |
| streptococcus_pyogenes | protection_fus_b | 7.5e-9 |
| streptococcus_pyogenes | protection_tet_m | 7.5e-6 |
| streptococcus_pyogenes | mutation_23s_rrna | 5e-11 |
| streptococcus_agalactiae | target_site_pbp2a_meca | 1e-8 |
| streptococcus_agalactiae | target_site_van_a | 1e-8 |
| streptococcus_agalactiae | target_site_van_b | 1e-8 |
| streptococcus_agalactiae | mutation_gyra_primary | 1e-6 |
| streptococcus_agalactiae | mutation_gyra_parc_secondary | 1e-6 |
| streptococcus_agalactiae | target_site_erm_b | 5e-5 |
| streptococcus_agalactiae | target_site_cfr | 1e-8 |
| streptococcus_agalactiae | enzyme_cat | 3e-7 |
| streptococcus_agalactiae | global_efflux_pump | 1e-6 |
| streptococcus_agalactiae | mutation_folate_pathway | 1e-5 |
| streptococcus_agalactiae | mutation_mpr_f | 1e-7 |
| streptococcus_agalactiae | mutation_rpo_b | 1e-7 |
| streptococcus_agalactiae | protection_fus_b | 1e-7 |
| streptococcus_agalactiae | protection_tet_m | 0.001 |
| streptococcus_agalactiae | mutation_23s_rrna | 5e-11 |
| streptococcus_agalactiae | mutation_pbp_mosaic | 1e-6 |
| haemophilus_influenzae | enzyme_esbl_ctx_m | 5e-7 |
| haemophilus_influenzae | enzyme_esbl_tem | 5e-7 |
| haemophilus_influenzae | enzyme_esbl_shv | 5e-7 |
| haemophilus_influenzae | enzyme_kpc | 5e-8 |
| haemophilus_influenzae | enzyme_ndm_vim | 5e-8 |
| haemophilus_influenzae | enzyme_oxa_48 | 5e-8 |
| haemophilus_influenzae | enzyme_ampc_cmy | 5e-7 |
| haemophilus_influenzae | enzyme_ampc_dha | 5e-7 |
| haemophilus_influenzae | mutation_gyra_primary | 3e-4 |
| haemophilus_influenzae | mutation_gyra_parc_secondary | 8e-5 |
| haemophilus_influenzae | protection_qnr | 1e-4 |
| haemophilus_influenzae | enzyme_16s_rrmt | 0.05 |
| haemophilus_influenzae | target_site_erm_b | 3e-5 |
| haemophilus_influenzae | target_site_cfr | 3e-5 |
| haemophilus_influenzae | enzyme_cat | 1e-4 |
| haemophilus_influenzae | modification_mcr_1 | 2e-6 |
| haemophilus_influenzae | global_efflux_pump | 3e-4 |
| haemophilus_influenzae | global_porin_loss | 3e-7 |
| haemophilus_influenzae | mutation_folate_pathway | 8e-4 |
| haemophilus_influenzae | mutation_nitroreductase | 2e-5 |
| haemophilus_influenzae | mutation_rpo_b | 0.015 |
| haemophilus_influenzae | protection_tet_m | 8e-4 |
| haemophilus_influenzae | enzyme_aac_aph | 0.05 |
| haemophilus_influenzae | enzyme_bla_z | 5e-7 |
| haemophilus_influenzae | mutation_23s_rrna | 3e-6 |
| haemophilus_influenzae | mutation_pbp_mosaic | 5e-7 |
| haemophilus_influenzae | efflux_mtr_cde | 5e-7 |
| chlamydia_trachomatis | mutation_gyra_primary | 2e-7 |
| chlamydia_trachomatis | mutation_gyra_parc_secondary | 2e-7 |
| chlamydia_trachomatis | target_site_erm_b | 1e-7 |
| chlamydia_trachomatis | target_site_cfr | 2e-9 |
| chlamydia_trachomatis | enzyme_cat | 1.3e-9 |
| chlamydia_trachomatis | global_efflux_pump | 2e-8 |
| chlamydia_trachomatis | mutation_folate_pathway | 2e-9 |
| chlamydia_trachomatis | mutation_nitroreductase | 2e-9 |
| chlamydia_trachomatis | mutation_rpo_b | 2e-8 |
| chlamydia_trachomatis | protection_tet_m | 2e-7 |
| chlamydia_trachomatis | mutation_23s_rrna | 1e-10 |
| mycoplasma_genitalium | mutation_gyra_primary | 0.3 |
| mycoplasma_genitalium | mutation_gyra_parc_secondary | 0.3 |
| mycoplasma_genitalium | target_site_erm_b | 0.025 |
| mycoplasma_genitalium | target_site_cfr | 0.025 |
| mycoplasma_genitalium | enzyme_cat | 0.01 |
| mycoplasma_genitalium | global_efflux_pump | 0.01 |
| mycoplasma_genitalium | mutation_folate_pathway | 0.003 |
| mycoplasma_genitalium | mutation_nitroreductase | 0.003 |
| mycoplasma_genitalium | mutation_rpo_b | 0.03 |
| mycoplasma_genitalium | protection_tet_m | 0.06 |
| mycoplasma_genitalium | mutation_23s_rrna | 0.03 |
| vibrio_cholerae | enzyme_esbl_ctx_m | 3e-6 |
| vibrio_cholerae | enzyme_esbl_tem | 3e-6 |
| vibrio_cholerae | enzyme_esbl_shv | 1e-6 |
| vibrio_cholerae | enzyme_kpc | 1e-7 |
| vibrio_cholerae | enzyme_ndm_vim | 1e-6 |
| vibrio_cholerae | enzyme_oxa_48 | 1e-7 |
| vibrio_cholerae | enzyme_ampc_cmy | 1e-6 |
| vibrio_cholerae | enzyme_ampc_dha | 1e-6 |
| vibrio_cholerae | mutation_gyra_primary | 3e-5 |
| vibrio_cholerae | mutation_gyra_parc_secondary | 1.5e-5 |
| vibrio_cholerae | protection_qnr | 3e-6 |
| vibrio_cholerae | enzyme_16s_rrmt | 3e-8 |
| vibrio_cholerae | enzyme_cat | 1e-4 |
| vibrio_cholerae | efflux_acrab_tolc | 3e-6 |
| vibrio_cholerae | modification_mcr_1 | 3e-7 |
| vibrio_cholerae | global_efflux_pump | 3e-6 |
| vibrio_cholerae | global_porin_loss | 1.5e-6 |
| vibrio_cholerae | mutation_folate_pathway | 1.5e-4 |
| vibrio_cholerae | mutation_nitroreductase | 3e-7 |
| vibrio_cholerae | enzyme_fos_a | 3e-7 |
| vibrio_cholerae | mutation_rpo_b | 3e-7 |
| vibrio_cholerae | protection_tet_m | 5e-5 |
| vibrio_cholerae | enzyme_aac_aph | 5e-11 |
| vibrio_cholerae | mutation_23s_rrna | 1.5e-4 |
| vibrio_cholerae | efflux_tet_abc | 5e-11 |
| vibrio_cholerae | efflux_mtr_cde | 1e-9 |
| neisseria_meningitidis | enzyme_esbl_ctx_m | 1e-8 |
| neisseria_meningitidis | enzyme_esbl_tem | 1e-8 |
| neisseria_meningitidis | enzyme_esbl_shv | 1e-8 |
| neisseria_meningitidis | enzyme_ampc_cmy | 1e-8 |
| neisseria_meningitidis | enzyme_ampc_dha | 1e-8 |
| neisseria_meningitidis | mutation_gyra_primary | 3e-5 |
| neisseria_meningitidis | mutation_gyra_parc_secondary | 1e-5 |
| neisseria_meningitidis | protection_qnr | 3e-7 |
| neisseria_meningitidis | enzyme_16s_rrmt | 1e-7 |
| neisseria_meningitidis | target_site_erm_b | 1e-7 |
| neisseria_meningitidis | target_site_cfr | 2e-7 |
| neisseria_meningitidis | enzyme_cat | 2e-6 |
| neisseria_meningitidis | efflux_acrab_tolc | 3e-7 |
| neisseria_meningitidis | modification_mcr_1 | 3e-7 |
| neisseria_meningitidis | global_efflux_pump | 2e-6 |
| neisseria_meningitidis | global_porin_loss | 5e-7 |
| neisseria_meningitidis | mutation_folate_pathway | 5e-6 |
| neisseria_meningitidis | mutation_nitroreductase | 3e-6 |
| neisseria_meningitidis | mutation_rpo_b | 3e-6 |
| neisseria_meningitidis | protection_tet_m | 5e-6 |
| neisseria_meningitidis | mutation_23s_rrna | 1e-7 |
| neisseria_meningitidis | efflux_tet_abc | 5e-6 |
| neisseria_meningitidis | mutation_pbp_mosaic | 1e-6 |
| neisseria_meningitidis | efflux_mtr_cde | 1e-9 |
| listeria_monocytogenes | target_site_van_a | 3.8e-6 |
| listeria_monocytogenes | target_site_van_b | 3.8e-6 |
| listeria_monocytogenes | mutation_gyra_primary | 1.9e-4 |
| listeria_monocytogenes | mutation_gyra_parc_secondary | 3.8e-5 |
| listeria_monocytogenes | target_site_erm_b | 3.8e-4 |
| listeria_monocytogenes | target_site_cfr | 3.8e-6 |
| listeria_monocytogenes | enzyme_cat | 1.9e-4 |
| listeria_monocytogenes | global_efflux_pump | 1.9e-4 |
| listeria_monocytogenes | mutation_folate_pathway | 3.8e-4 |
| listeria_monocytogenes | mutation_mpr_f | 3.8e-5 |
| listeria_monocytogenes | mutation_rpo_b | 3.8e-5 |
| listeria_monocytogenes | protection_fus_b | 3.8e-5 |
| listeria_monocytogenes | protection_tet_m | 0.0019 |
| clostridioides_difficile | mutation_gyra_primary | 6e-5 |
| clostridioides_difficile | mutation_gyra_parc_secondary | 1e-5 |
| clostridioides_difficile | enzyme_16s_rrmt | 1e-7 |
| clostridioides_difficile | target_site_erm_b | 6e-5 |
| clostridioides_difficile | target_site_cfr | 1e-6 |
| clostridioides_difficile | enzyme_cat | 6e-6 |
| clostridioides_difficile | global_efflux_pump | 1e-5 |
| clostridioides_difficile | mutation_folate_pathway | 6e-6 |
| clostridioides_difficile | mutation_nitroreductase | 1.2e-4 |
| clostridioides_difficile | mutation_rpo_b | 6e-5 |
| clostridioides_difficile | protection_tet_m | 6e-5 |
| bacteroides_fragilis | enzyme_esbl_ctx_m | 0.003 |
| bacteroides_fragilis | enzyme_esbl_tem | 0.003 |
| bacteroides_fragilis | enzyme_esbl_shv | 0.003 |
| bacteroides_fragilis | enzyme_kpc | 3e-5 |
| bacteroides_fragilis | enzyme_ndm_vim | 3e-5 |
| bacteroides_fragilis | enzyme_oxa_48 | 3e-5 |
| bacteroides_fragilis | enzyme_ampc_cmy | 0.001 |
| bacteroides_fragilis | enzyme_ampc_dha | 0.001 |
| bacteroides_fragilis | mutation_gyra_primary | 5e-4 |
| bacteroides_fragilis | mutation_gyra_parc_secondary | 0.01 |
| bacteroides_fragilis | protection_qnr | 2e-4 |
| bacteroides_fragilis | enzyme_16s_rrmt | 1 |
| bacteroides_fragilis | target_site_erm_b | 0.02 |
| bacteroides_fragilis | target_site_cfr | 2e-4 |
| bacteroides_fragilis | enzyme_cat | 2e-5 |
| bacteroides_fragilis | efflux_acrab_tolc | 1e-4 |
| bacteroides_fragilis | modification_mcr_1 | 0.002 |
| bacteroides_fragilis | global_efflux_pump | 1e-4 |
| bacteroides_fragilis | global_porin_loss | 1e-5 |
| bacteroides_fragilis | mutation_folate_pathway | 0.005 |
| bacteroides_fragilis | mutation_nitroreductase | 5e-4 |
| bacteroides_fragilis | mutation_rpo_b | 0.002 |
| bacteroides_fragilis | protection_tet_m | 0.01 |
| bacteroides_fragilis | enzyme_aac_aph | 1 |
| bacteroides_fragilis | mutation_pbp_mosaic | 0.001 |
| campylobacter_jejuni | mutation_gyra_primary | 0.03 |
| campylobacter_jejuni | mutation_gyra_parc_secondary | 0.03 |
| campylobacter_jejuni | target_site_erm_b | 1e-5 |
| campylobacter_jejuni | target_site_cfr | 1e-5 |
| campylobacter_jejuni | enzyme_cat | 6e-4 |
| campylobacter_jejuni | global_efflux_pump | 0.03 |
| campylobacter_jejuni | global_porin_loss | 2e-5 |
| campylobacter_jejuni | mutation_folate_pathway | 0.6 |
| campylobacter_jejuni | mutation_rpo_b | 0.6 |
| campylobacter_jejuni | protection_tet_m | 0.01 |
| campylobacter_jejuni | enzyme_aac_aph | 0.1 |
| campylobacter_jejuni | mutation_23s_rrna | 3e-4 |
| campylobacter_jejuni | efflux_tet_abc | 0.001 |
| campylobacter_jejuni | efflux_mtr_cde | 2e-4 |
| enterobacter_cloacae | enzyme_esbl_ctx_m | 5e-6 |
| enterobacter_cloacae | enzyme_esbl_tem | 5e-6 |
| enterobacter_cloacae | enzyme_esbl_shv | 5e-6 |
| enterobacter_cloacae | enzyme_kpc | 1e-5 |
| enterobacter_cloacae | enzyme_ndm_vim | 1e-5 |
| enterobacter_cloacae | enzyme_oxa_48 | 1e-5 |
| enterobacter_cloacae | enzyme_ampc_cmy | 5e-6 |
| enterobacter_cloacae | enzyme_ampc_dha | 5e-6 |
| enterobacter_cloacae | mutation_gyra_primary | 2e-4 |
| enterobacter_cloacae | mutation_gyra_parc_secondary | 2e-4 |
| enterobacter_cloacae | protection_qnr | 2e-4 |
| enterobacter_cloacae | enzyme_16s_rrmt | 0.03 |
| enterobacter_cloacae | enzyme_cat | 2e-4 |
| enterobacter_cloacae | efflux_acrab_tolc | 5e-5 |
| enterobacter_cloacae | modification_mcr_1 | 0.01 |
| enterobacter_cloacae | global_porin_loss | 1e-6 |
| enterobacter_cloacae | mutation_folate_pathway | 5e-5 |
| enterobacter_cloacae | mutation_nitroreductase | 0.001 |
| enterobacter_cloacae | enzyme_fos_a | 0.001 |
| enterobacter_cloacae | mutation_rpo_b | 0.01 |
| enterobacter_cloacae | protection_tet_m | 0.002 |
| enterobacter_cloacae | enzyme_aac_aph | 0.1 |
| enterobacter_cloacae | efflux_tet_abc | 1e-4 |
| enterobacter_cloacae | mutation_pbp_mosaic | 5e-6 |
| yersinia_enterocolitica | enzyme_esbl_ctx_m | 3e-10 |
| yersinia_enterocolitica | enzyme_esbl_tem | 3e-10 |
| yersinia_enterocolitica | enzyme_esbl_shv | 1e-10 |
| yersinia_enterocolitica | enzyme_kpc | 3e-11 |
| yersinia_enterocolitica | enzyme_ndm_vim | 3e-11 |
| yersinia_enterocolitica | enzyme_oxa_48 | 3e-11 |
| yersinia_enterocolitica | enzyme_ampc_cmy | 3e-10 |
| yersinia_enterocolitica | enzyme_ampc_dha | 3e-10 |
| yersinia_enterocolitica | mutation_gyra_primary | 3e-10 |
| yersinia_enterocolitica | mutation_gyra_parc_secondary | 3e-10 |
| yersinia_enterocolitica | protection_qnr | 3e-10 |
| yersinia_enterocolitica | enzyme_16s_rrmt | 3e-11 |
| yersinia_enterocolitica | enzyme_cat | 3e-10 |
| yersinia_enterocolitica | efflux_acrab_tolc | 3e-10 |
| yersinia_enterocolitica | modification_mcr_1 | 3e-10 |
| yersinia_enterocolitica | global_efflux_pump | 3e-10 |
| yersinia_enterocolitica | global_porin_loss | 1e-10 |
| yersinia_enterocolitica | mutation_folate_pathway | 3e-10 |
| yersinia_enterocolitica | mutation_nitroreductase | 1e-9 |
| yersinia_enterocolitica | enzyme_fos_a | 3e-10 |
| yersinia_enterocolitica | mutation_rpo_b | 3e-11 |
| yersinia_enterocolitica | protection_tet_m | 3e-10 |
| yersinia_enterocolitica | enzyme_aac_aph | 1.5e-10 |
| yersinia_enterocolitica | efflux_tet_abc | 1.5e-10 |
| yersinia_enterocolitica | efflux_mtr_cde | 3e-10 |
| moraxella_catarrhalis | enzyme_esbl_ctx_m | 2e-7 |
| moraxella_catarrhalis | enzyme_esbl_tem | 1e-6 |
| moraxella_catarrhalis | enzyme_esbl_shv | 2e-8 |
| moraxella_catarrhalis | enzyme_ampc_cmy | 5e-7 |
| moraxella_catarrhalis | enzyme_ampc_dha | 2e-7 |
| moraxella_catarrhalis | mutation_gyra_primary | 2e-6 |
| moraxella_catarrhalis | mutation_gyra_parc_secondary | 5e-7 |
| moraxella_catarrhalis | protection_qnr | 2e-7 |
| moraxella_catarrhalis | enzyme_16s_rrmt | 2e-8 |
| moraxella_catarrhalis | target_site_erm_b | 5e-6 |
| moraxella_catarrhalis | target_site_cfr | 5e-8 |
| moraxella_catarrhalis | enzyme_cat | 5e-7 |
| moraxella_catarrhalis | efflux_acrab_tolc | 5e-7 |
| moraxella_catarrhalis | modification_mcr_1 | 2e-8 |
| moraxella_catarrhalis | global_efflux_pump | 5e-6 |
| moraxella_catarrhalis | global_porin_loss | 2e-7 |
| moraxella_catarrhalis | mutation_folate_pathway | 1e-5 |
| moraxella_catarrhalis | mutation_nitroreductase | 2e-7 |
| moraxella_catarrhalis | mutation_rpo_b | 2e-7 |
| moraxella_catarrhalis | protection_tet_m | 1e-5 |
| moraxella_catarrhalis | mutation_pbp_mosaic | 5e-6 |
| moraxella_catarrhalis | efflux_mtr_cde | 5e-7 |
| treponema_pallidum | mutation_gyra_primary | 1.5e-6 |
| treponema_pallidum | mutation_gyra_parc_secondary | 7.5e-7 |
| treponema_pallidum | enzyme_cat | 1.5e-7 |
| treponema_pallidum | global_efflux_pump | 1.5e-7 |
| treponema_pallidum | mutation_folate_pathway | 1.5e-7 |
| treponema_pallidum | mutation_rpo_b | 1.5e-7 |
| treponema_pallidum | protection_tet_m | 1.5e-6 |
| bordetella_pertussis | enzyme_esbl_ctx_m | 1e-10 |
| bordetella_pertussis | enzyme_esbl_tem | 1e-10 |
| bordetella_pertussis | enzyme_esbl_shv | 1e-10 |
| bordetella_pertussis | enzyme_ampc_cmy | 1e-10 |
| bordetella_pertussis | enzyme_ampc_dha | 1e-10 |
| bordetella_pertussis | mutation_gyra_primary | 2e-12 |
| bordetella_pertussis | mutation_gyra_parc_secondary | 1e-12 |
| bordetella_pertussis | enzyme_16s_rrmt | 1e-10 |
| bordetella_pertussis | target_site_cfr | 5e-11 |
| bordetella_pertussis | enzyme_cat | 3.8e-12 |
| bordetella_pertussis | efflux_acrab_tolc | 1e-12 |
| bordetella_pertussis | global_efflux_pump | 5e-12 |
| bordetella_pertussis | global_porin_loss | 1e-12 |
| bordetella_pertussis | mutation_folate_pathway | 1e-11 |
| bordetella_pertussis | mutation_nitroreductase | 1e-10 |
| bordetella_pertussis | mutation_rpo_b | 1e-12 |
| bordetella_pertussis | protection_tet_m | 1e-11 |
| bordetella_pertussis | efflux_mtr_cde | 2e-12 |
| helicobacter_pylori | mutation_gyra_primary | 1 |
| helicobacter_pylori | mutation_gyra_parc_secondary | 1 |
| helicobacter_pylori | target_site_erm_b | 1 |
| helicobacter_pylori | target_site_cfr | 1 |
| helicobacter_pylori | global_efflux_pump | 1 |
| helicobacter_pylori | global_porin_loss | 1 |
| helicobacter_pylori | mutation_folate_pathway | 1 |
| helicobacter_pylori | mutation_nitroreductase | 1 |
| helicobacter_pylori | mutation_rpo_b | 1 |
| helicobacter_pylori | protection_tet_m | 1 |
| helicobacter_pylori | enzyme_bla_z | 1 |
| helicobacter_pylori | mutation_23s_rrna | 1 |
| helicobacter_pylori | mutation_pbp_mosaic | 1 |
| mycoplasma_pneumoniae | mutation_gyra_primary | 3e-7 |
| mycoplasma_pneumoniae | mutation_gyra_parc_secondary | 1.5e-7 |
| mycoplasma_pneumoniae | target_site_erm_b | 1.5e-5 |
| mycoplasma_pneumoniae | target_site_cfr | 3e-9 |
| mycoplasma_pneumoniae | enzyme_cat | 3e-9 |
| mycoplasma_pneumoniae | global_efflux_pump | 1.5e-7 |
| mycoplasma_pneumoniae | mutation_folate_pathway | 3e-9 |
| mycoplasma_pneumoniae | mutation_nitroreductase | 3e-9 |
| mycoplasma_pneumoniae | mutation_rpo_b | 3e-8 |
| mycoplasma_pneumoniae | protection_tet_m | 3e-7 |
| mycoplasma_pneumoniae | mutation_23s_rrna | 5e-9 |
| legionella_pneumophila | enzyme_esbl_ctx_m | 3e-7 |
| legionella_pneumophila | enzyme_esbl_tem | 3e-7 |
| legionella_pneumophila | enzyme_esbl_shv | 3e-7 |
| legionella_pneumophila | enzyme_ampc_cmy | 3e-7 |
| legionella_pneumophila | enzyme_ampc_dha | 3e-7 |
| legionella_pneumophila | mutation_gyra_primary | 3e-5 |
| legionella_pneumophila | mutation_gyra_parc_secondary | 3e-5 |
| legionella_pneumophila | protection_qnr | 1e-7 |
| legionella_pneumophila | enzyme_16s_rrmt | 3e-7 |
| legionella_pneumophila | target_site_erm_b | 3e-5 |
| legionella_pneumophila | target_site_cfr | 3e-7 |
| legionella_pneumophila | enzyme_cat | 3e-6 |
| legionella_pneumophila | efflux_acrab_tolc | 3e-7 |
| legionella_pneumophila | modification_mcr_1 | 3e-7 |
| legionella_pneumophila | global_efflux_pump | 3e-5 |
| legionella_pneumophila | global_porin_loss | 1e-6 |
| legionella_pneumophila | mutation_folate_pathway | 3e-6 |
| legionella_pneumophila | mutation_nitroreductase | 3e-7 |
| legionella_pneumophila | mutation_rpo_b | 3e-6 |
| legionella_pneumophila | protection_tet_m | 3e-5 |
| legionella_pneumophila | mutation_23s_rrna | 5e-11 |
| burkholderia_cepacia_complex | enzyme_esbl_ctx_m | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_esbl_tem | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_esbl_shv | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_kpc | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_ndm_vim | 3.7e-5 |
| burkholderia_cepacia_complex | enzyme_oxa_48 | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_ampc_cmy | 3.8e-4 |
| burkholderia_cepacia_complex | enzyme_ampc_dha | 7.5e-5 |
| burkholderia_cepacia_complex | mutation_gyra_primary | 3.8e-4 |
| burkholderia_cepacia_complex | mutation_gyra_parc_secondary | 7.5e-5 |
| burkholderia_cepacia_complex | protection_qnr | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_16s_rrmt | 7.5e-6 |
| burkholderia_cepacia_complex | enzyme_cat | 3.7e-5 |
| burkholderia_cepacia_complex | modification_mcr_1 | 7.5e-6 |
| burkholderia_cepacia_complex | global_efflux_pump | 7.5e-4 |
| burkholderia_cepacia_complex | global_porin_loss | 3.8e-4 |
| burkholderia_cepacia_complex | mutation_folate_pathway | 7.5e-5 |
| burkholderia_cepacia_complex | enzyme_fos_a | 3.7e-5 |
| burkholderia_cepacia_complex | mutation_rpo_b | 7.5e-6 |
| burkholderia_cepacia_complex | protection_tet_m | 3.7e-5 |
| burkholderia_cepacia_complex | enzyme_aac_aph | 5e-11 |
| burkholderia_cepacia_complex | efflux_tet_abc | 5e-11 |

### B.11 Horizontal Gene Transfer Matrix

Per-day probability of horizontal gene transfer of resistance between co-colonising bacterial species. Only non-zero entries shown.

See: [§9.1 Transfer compatibility](#91-transfer-compatibility), [§9.2 The HGT process](#92-the-hgt-process).

| Donor | Recipient | Probability/day |
| --- | ---: | ---: |
| acinetobacter_baumannii | citrobacter_spp. | 1e-10 |
| acinetobacter_baumannii | enterobacter_spp. | 1e-10 |
| acinetobacter_baumannii | escherichia_coli | 1e-10 |
| acinetobacter_baumannii | klebsiella_pneumoniae | 1e-10 |
| acinetobacter_baumannii | morganella_spp. | 1e-10 |
| acinetobacter_baumannii | proteus_spp. | 1e-10 |
| acinetobacter_baumannii | serratia_spp. | 1e-10 |
| acinetobacter_baumannii | p_stuartii | 1e-10 |
| acinetobacter_baumannii | pseudomonas_aeruginosa | 1e-9 |
| acinetobacter_baumannii | stenotrophomonas_maltophilia | 1e-9 |
| acinetobacter_baumannii | salmonella_enterica_serovar_typhi | 1e-10 |
| acinetobacter_baumannii | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| acinetobacter_baumannii | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| acinetobacter_baumannii | shigella_spp. | 1e-10 |
| acinetobacter_baumannii | neisseria_gonorrhoeae | 3e-11 |
| acinetobacter_baumannii | haemophilus_influenzae | 3e-11 |
| acinetobacter_baumannii | chlamydia_trachomatis | 3e-11 |
| acinetobacter_baumannii | mycoplasma_genitalium | 3e-11 |
| acinetobacter_baumannii | vibrio_cholerae | 1e-10 |
| acinetobacter_baumannii | neisseria_meningitidis | 3e-11 |
| acinetobacter_baumannii | clostridioides_difficile | 3e-11 |
| acinetobacter_baumannii | bacteroides_fragilis | 3e-11 |
| acinetobacter_baumannii | enterobacter_cloacae | 1e-10 |
| acinetobacter_baumannii | yersinia_enterocolitica | 1e-10 |
| acinetobacter_baumannii | moraxella_catarrhalis | 3e-11 |
| acinetobacter_baumannii | bordetella_pertussis | 3e-11 |
| acinetobacter_baumannii | mycoplasma_pneumoniae | 3e-11 |
| acinetobacter_baumannii | legionella_pneumophila | 3e-11 |
| acinetobacter_baumannii | burkholderia_cepacia_complex | 1e-9 |
| citrobacter_spp. | acinetobacter_baumannii | 1e-10 |
| citrobacter_spp. | enterobacter_spp. | 1e-9 |
| citrobacter_spp. | escherichia_coli | 1e-9 |
| citrobacter_spp. | klebsiella_pneumoniae | 1e-9 |
| citrobacter_spp. | morganella_spp. | 1e-9 |
| citrobacter_spp. | proteus_spp. | 1e-9 |
| citrobacter_spp. | serratia_spp. | 1e-9 |
| citrobacter_spp. | p_stuartii | 1e-9 |
| citrobacter_spp. | pseudomonas_aeruginosa | 1e-10 |
| citrobacter_spp. | stenotrophomonas_maltophilia | 1e-10 |
| citrobacter_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| citrobacter_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| citrobacter_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| citrobacter_spp. | shigella_spp. | 1e-9 |
| citrobacter_spp. | neisseria_gonorrhoeae | 3e-11 |
| citrobacter_spp. | haemophilus_influenzae | 3e-11 |
| citrobacter_spp. | chlamydia_trachomatis | 3e-11 |
| citrobacter_spp. | mycoplasma_genitalium | 3e-11 |
| citrobacter_spp. | vibrio_cholerae | 1e-10 |
| citrobacter_spp. | neisseria_meningitidis | 3e-11 |
| citrobacter_spp. | clostridioides_difficile | 3e-11 |
| citrobacter_spp. | bacteroides_fragilis | 3e-11 |
| citrobacter_spp. | enterobacter_cloacae | 1e-9 |
| citrobacter_spp. | yersinia_enterocolitica | 1e-9 |
| citrobacter_spp. | moraxella_catarrhalis | 3e-11 |
| citrobacter_spp. | bordetella_pertussis | 3e-11 |
| citrobacter_spp. | mycoplasma_pneumoniae | 3e-11 |
| citrobacter_spp. | legionella_pneumophila | 3e-11 |
| citrobacter_spp. | burkholderia_cepacia_complex | 1e-10 |
| enterobacter_spp. | acinetobacter_baumannii | 1e-10 |
| enterobacter_spp. | citrobacter_spp. | 1e-9 |
| enterobacter_spp. | escherichia_coli | 1e-9 |
| enterobacter_spp. | klebsiella_pneumoniae | 1e-9 |
| enterobacter_spp. | morganella_spp. | 1e-9 |
| enterobacter_spp. | proteus_spp. | 1e-9 |
| enterobacter_spp. | serratia_spp. | 1e-9 |
| enterobacter_spp. | p_stuartii | 1e-9 |
| enterobacter_spp. | pseudomonas_aeruginosa | 1e-10 |
| enterobacter_spp. | stenotrophomonas_maltophilia | 1e-10 |
| enterobacter_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| enterobacter_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| enterobacter_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| enterobacter_spp. | shigella_spp. | 1e-9 |
| enterobacter_spp. | neisseria_gonorrhoeae | 3e-11 |
| enterobacter_spp. | haemophilus_influenzae | 3e-11 |
| enterobacter_spp. | chlamydia_trachomatis | 3e-11 |
| enterobacter_spp. | mycoplasma_genitalium | 3e-11 |
| enterobacter_spp. | vibrio_cholerae | 1e-10 |
| enterobacter_spp. | neisseria_meningitidis | 3e-11 |
| enterobacter_spp. | clostridioides_difficile | 3e-11 |
| enterobacter_spp. | bacteroides_fragilis | 3e-11 |
| enterobacter_spp. | enterobacter_cloacae | 1e-9 |
| enterobacter_spp. | yersinia_enterocolitica | 1e-9 |
| enterobacter_spp. | moraxella_catarrhalis | 3e-11 |
| enterobacter_spp. | bordetella_pertussis | 3e-11 |
| enterobacter_spp. | mycoplasma_pneumoniae | 3e-11 |
| enterobacter_spp. | legionella_pneumophila | 3e-11 |
| enterobacter_spp. | burkholderia_cepacia_complex | 1e-10 |
| enterococcus_faecalis | enterococcus_faecium | 1e-9 |
| enterococcus_faecalis | staphylococcus_aureus | 1e-9 |
| enterococcus_faecalis | staphylococcus_epidermidis | 1e-9 |
| enterococcus_faecalis | streptococcus_pneumoniae | 1e-9 |
| enterococcus_faecalis | streptococcus_pyogenes | 1e-9 |
| enterococcus_faecalis | streptococcus_agalactiae | 1e-9 |
| enterococcus_faecalis | listeria_monocytogenes | 1e-9 |
| enterococcus_faecium | enterococcus_faecalis | 1e-9 |
| enterococcus_faecium | staphylococcus_aureus | 1e-9 |
| enterococcus_faecium | staphylococcus_epidermidis | 1e-9 |
| enterococcus_faecium | streptococcus_pneumoniae | 1e-9 |
| enterococcus_faecium | streptococcus_pyogenes | 1e-9 |
| enterococcus_faecium | streptococcus_agalactiae | 1e-9 |
| enterococcus_faecium | listeria_monocytogenes | 1e-9 |
| escherichia_coli | acinetobacter_baumannii | 1e-10 |
| escherichia_coli | citrobacter_spp. | 1e-9 |
| escherichia_coli | enterobacter_spp. | 1e-9 |
| escherichia_coli | klebsiella_pneumoniae | 1e-9 |
| escherichia_coli | morganella_spp. | 1e-9 |
| escherichia_coli | proteus_spp. | 1e-9 |
| escherichia_coli | serratia_spp. | 1e-9 |
| escherichia_coli | p_stuartii | 1e-9 |
| escherichia_coli | pseudomonas_aeruginosa | 1e-10 |
| escherichia_coli | stenotrophomonas_maltophilia | 1e-10 |
| escherichia_coli | salmonella_enterica_serovar_typhi | 1e-9 |
| escherichia_coli | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| escherichia_coli | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| escherichia_coli | shigella_spp. | 1e-9 |
| escherichia_coli | neisseria_gonorrhoeae | 3e-11 |
| escherichia_coli | haemophilus_influenzae | 3e-11 |
| escherichia_coli | chlamydia_trachomatis | 3e-11 |
| escherichia_coli | mycoplasma_genitalium | 3e-11 |
| escherichia_coli | vibrio_cholerae | 1e-10 |
| escherichia_coli | neisseria_meningitidis | 3e-11 |
| escherichia_coli | clostridioides_difficile | 3e-11 |
| escherichia_coli | bacteroides_fragilis | 3e-11 |
| escherichia_coli | enterobacter_cloacae | 1e-9 |
| escherichia_coli | yersinia_enterocolitica | 1e-9 |
| escherichia_coli | moraxella_catarrhalis | 3e-11 |
| escherichia_coli | bordetella_pertussis | 3e-11 |
| escherichia_coli | mycoplasma_pneumoniae | 3e-11 |
| escherichia_coli | legionella_pneumophila | 3e-11 |
| escherichia_coli | burkholderia_cepacia_complex | 1e-10 |
| klebsiella_pneumoniae | acinetobacter_baumannii | 1e-10 |
| klebsiella_pneumoniae | citrobacter_spp. | 1e-9 |
| klebsiella_pneumoniae | enterobacter_spp. | 1e-9 |
| klebsiella_pneumoniae | escherichia_coli | 1e-9 |
| klebsiella_pneumoniae | morganella_spp. | 1e-9 |
| klebsiella_pneumoniae | proteus_spp. | 1e-9 |
| klebsiella_pneumoniae | serratia_spp. | 1e-9 |
| klebsiella_pneumoniae | p_stuartii | 1e-9 |
| klebsiella_pneumoniae | pseudomonas_aeruginosa | 1e-10 |
| klebsiella_pneumoniae | stenotrophomonas_maltophilia | 1e-10 |
| klebsiella_pneumoniae | salmonella_enterica_serovar_typhi | 1e-9 |
| klebsiella_pneumoniae | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| klebsiella_pneumoniae | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| klebsiella_pneumoniae | shigella_spp. | 1e-9 |
| klebsiella_pneumoniae | neisseria_gonorrhoeae | 3e-11 |
| klebsiella_pneumoniae | haemophilus_influenzae | 3e-11 |
| klebsiella_pneumoniae | chlamydia_trachomatis | 3e-11 |
| klebsiella_pneumoniae | mycoplasma_genitalium | 3e-11 |
| klebsiella_pneumoniae | vibrio_cholerae | 1e-10 |
| klebsiella_pneumoniae | neisseria_meningitidis | 3e-11 |
| klebsiella_pneumoniae | clostridioides_difficile | 3e-11 |
| klebsiella_pneumoniae | bacteroides_fragilis | 3e-11 |
| klebsiella_pneumoniae | enterobacter_cloacae | 1e-9 |
| klebsiella_pneumoniae | yersinia_enterocolitica | 1e-9 |
| klebsiella_pneumoniae | moraxella_catarrhalis | 3e-11 |
| klebsiella_pneumoniae | bordetella_pertussis | 3e-11 |
| klebsiella_pneumoniae | mycoplasma_pneumoniae | 3e-11 |
| klebsiella_pneumoniae | legionella_pneumophila | 3e-11 |
| klebsiella_pneumoniae | burkholderia_cepacia_complex | 1e-10 |
| morganella_spp. | acinetobacter_baumannii | 1e-10 |
| morganella_spp. | citrobacter_spp. | 1e-9 |
| morganella_spp. | enterobacter_spp. | 1e-9 |
| morganella_spp. | escherichia_coli | 1e-9 |
| morganella_spp. | klebsiella_pneumoniae | 1e-9 |
| morganella_spp. | proteus_spp. | 1e-9 |
| morganella_spp. | serratia_spp. | 1e-9 |
| morganella_spp. | p_stuartii | 1e-9 |
| morganella_spp. | pseudomonas_aeruginosa | 1e-10 |
| morganella_spp. | stenotrophomonas_maltophilia | 1e-10 |
| morganella_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| morganella_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| morganella_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| morganella_spp. | shigella_spp. | 1e-9 |
| morganella_spp. | neisseria_gonorrhoeae | 3e-11 |
| morganella_spp. | haemophilus_influenzae | 3e-11 |
| morganella_spp. | chlamydia_trachomatis | 3e-11 |
| morganella_spp. | mycoplasma_genitalium | 3e-11 |
| morganella_spp. | vibrio_cholerae | 1e-10 |
| morganella_spp. | neisseria_meningitidis | 3e-11 |
| morganella_spp. | clostridioides_difficile | 3e-11 |
| morganella_spp. | bacteroides_fragilis | 3e-11 |
| morganella_spp. | enterobacter_cloacae | 1e-9 |
| morganella_spp. | yersinia_enterocolitica | 1e-9 |
| morganella_spp. | moraxella_catarrhalis | 3e-11 |
| morganella_spp. | bordetella_pertussis | 3e-11 |
| morganella_spp. | mycoplasma_pneumoniae | 3e-11 |
| morganella_spp. | legionella_pneumophila | 3e-11 |
| morganella_spp. | burkholderia_cepacia_complex | 1e-10 |
| proteus_spp. | acinetobacter_baumannii | 1e-10 |
| proteus_spp. | citrobacter_spp. | 1e-9 |
| proteus_spp. | enterobacter_spp. | 1e-9 |
| proteus_spp. | escherichia_coli | 1e-9 |
| proteus_spp. | klebsiella_pneumoniae | 1e-9 |
| proteus_spp. | morganella_spp. | 1e-9 |
| proteus_spp. | serratia_spp. | 1e-9 |
| proteus_spp. | p_stuartii | 1e-9 |
| proteus_spp. | pseudomonas_aeruginosa | 1e-10 |
| proteus_spp. | stenotrophomonas_maltophilia | 1e-10 |
| proteus_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| proteus_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| proteus_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| proteus_spp. | shigella_spp. | 1e-9 |
| proteus_spp. | neisseria_gonorrhoeae | 3e-11 |
| proteus_spp. | haemophilus_influenzae | 3e-11 |
| proteus_spp. | chlamydia_trachomatis | 3e-11 |
| proteus_spp. | mycoplasma_genitalium | 3e-11 |
| proteus_spp. | vibrio_cholerae | 1e-10 |
| proteus_spp. | neisseria_meningitidis | 3e-11 |
| proteus_spp. | clostridioides_difficile | 3e-11 |
| proteus_spp. | bacteroides_fragilis | 3e-11 |
| proteus_spp. | enterobacter_cloacae | 1e-9 |
| proteus_spp. | yersinia_enterocolitica | 1e-9 |
| proteus_spp. | moraxella_catarrhalis | 3e-11 |
| proteus_spp. | bordetella_pertussis | 3e-11 |
| proteus_spp. | mycoplasma_pneumoniae | 3e-11 |
| proteus_spp. | legionella_pneumophila | 3e-11 |
| proteus_spp. | burkholderia_cepacia_complex | 1e-10 |
| serratia_spp. | acinetobacter_baumannii | 1e-10 |
| serratia_spp. | citrobacter_spp. | 1e-9 |
| serratia_spp. | enterobacter_spp. | 1e-9 |
| serratia_spp. | escherichia_coli | 1e-9 |
| serratia_spp. | klebsiella_pneumoniae | 1e-9 |
| serratia_spp. | morganella_spp. | 1e-9 |
| serratia_spp. | proteus_spp. | 1e-9 |
| serratia_spp. | p_stuartii | 1e-9 |
| serratia_spp. | pseudomonas_aeruginosa | 1e-10 |
| serratia_spp. | stenotrophomonas_maltophilia | 1e-10 |
| serratia_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| serratia_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| serratia_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| serratia_spp. | shigella_spp. | 1e-9 |
| serratia_spp. | neisseria_gonorrhoeae | 3e-11 |
| serratia_spp. | haemophilus_influenzae | 3e-11 |
| serratia_spp. | chlamydia_trachomatis | 3e-11 |
| serratia_spp. | mycoplasma_genitalium | 3e-11 |
| serratia_spp. | vibrio_cholerae | 1e-10 |
| serratia_spp. | neisseria_meningitidis | 3e-11 |
| serratia_spp. | clostridioides_difficile | 3e-11 |
| serratia_spp. | bacteroides_fragilis | 3e-11 |
| serratia_spp. | enterobacter_cloacae | 1e-9 |
| serratia_spp. | yersinia_enterocolitica | 1e-9 |
| serratia_spp. | moraxella_catarrhalis | 3e-11 |
| serratia_spp. | bordetella_pertussis | 3e-11 |
| serratia_spp. | mycoplasma_pneumoniae | 3e-11 |
| serratia_spp. | legionella_pneumophila | 3e-11 |
| serratia_spp. | burkholderia_cepacia_complex | 1e-10 |
| p_stuartii | acinetobacter_baumannii | 1e-10 |
| p_stuartii | citrobacter_spp. | 1e-9 |
| p_stuartii | enterobacter_spp. | 1e-9 |
| p_stuartii | escherichia_coli | 1e-9 |
| p_stuartii | klebsiella_pneumoniae | 1e-9 |
| p_stuartii | morganella_spp. | 1e-9 |
| p_stuartii | proteus_spp. | 1e-9 |
| p_stuartii | serratia_spp. | 1e-9 |
| p_stuartii | pseudomonas_aeruginosa | 1e-10 |
| p_stuartii | stenotrophomonas_maltophilia | 1e-10 |
| p_stuartii | salmonella_enterica_serovar_typhi | 1e-9 |
| p_stuartii | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| p_stuartii | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| p_stuartii | shigella_spp. | 1e-9 |
| p_stuartii | neisseria_gonorrhoeae | 3e-11 |
| p_stuartii | haemophilus_influenzae | 3e-11 |
| p_stuartii | chlamydia_trachomatis | 3e-11 |
| p_stuartii | mycoplasma_genitalium | 3e-11 |
| p_stuartii | vibrio_cholerae | 1e-10 |
| p_stuartii | neisseria_meningitidis | 3e-11 |
| p_stuartii | clostridioides_difficile | 3e-11 |
| p_stuartii | bacteroides_fragilis | 3e-11 |
| p_stuartii | enterobacter_cloacae | 1e-9 |
| p_stuartii | yersinia_enterocolitica | 1e-9 |
| p_stuartii | moraxella_catarrhalis | 3e-11 |
| p_stuartii | bordetella_pertussis | 3e-11 |
| p_stuartii | mycoplasma_pneumoniae | 3e-11 |
| p_stuartii | legionella_pneumophila | 3e-11 |
| p_stuartii | burkholderia_cepacia_complex | 1e-10 |
| pseudomonas_aeruginosa | acinetobacter_baumannii | 1e-9 |
| pseudomonas_aeruginosa | citrobacter_spp. | 1e-10 |
| pseudomonas_aeruginosa | enterobacter_spp. | 1e-10 |
| pseudomonas_aeruginosa | escherichia_coli | 1e-10 |
| pseudomonas_aeruginosa | klebsiella_pneumoniae | 1e-10 |
| pseudomonas_aeruginosa | morganella_spp. | 1e-10 |
| pseudomonas_aeruginosa | proteus_spp. | 1e-10 |
| pseudomonas_aeruginosa | serratia_spp. | 1e-10 |
| pseudomonas_aeruginosa | p_stuartii | 1e-10 |
| pseudomonas_aeruginosa | stenotrophomonas_maltophilia | 1e-9 |
| pseudomonas_aeruginosa | salmonella_enterica_serovar_typhi | 1e-10 |
| pseudomonas_aeruginosa | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| pseudomonas_aeruginosa | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| pseudomonas_aeruginosa | shigella_spp. | 1e-10 |
| pseudomonas_aeruginosa | neisseria_gonorrhoeae | 3e-11 |
| pseudomonas_aeruginosa | haemophilus_influenzae | 3e-11 |
| pseudomonas_aeruginosa | chlamydia_trachomatis | 3e-11 |
| pseudomonas_aeruginosa | mycoplasma_genitalium | 3e-11 |
| pseudomonas_aeruginosa | vibrio_cholerae | 1e-10 |
| pseudomonas_aeruginosa | neisseria_meningitidis | 3e-11 |
| pseudomonas_aeruginosa | clostridioides_difficile | 3e-11 |
| pseudomonas_aeruginosa | bacteroides_fragilis | 3e-11 |
| pseudomonas_aeruginosa | enterobacter_cloacae | 1e-10 |
| pseudomonas_aeruginosa | yersinia_enterocolitica | 1e-10 |
| pseudomonas_aeruginosa | moraxella_catarrhalis | 3e-11 |
| pseudomonas_aeruginosa | bordetella_pertussis | 3e-11 |
| pseudomonas_aeruginosa | mycoplasma_pneumoniae | 3e-11 |
| pseudomonas_aeruginosa | legionella_pneumophila | 3e-11 |
| pseudomonas_aeruginosa | burkholderia_cepacia_complex | 1e-9 |
| stenotrophomonas_maltophilia | acinetobacter_baumannii | 1e-9 |
| stenotrophomonas_maltophilia | citrobacter_spp. | 1e-10 |
| stenotrophomonas_maltophilia | enterobacter_spp. | 1e-10 |
| stenotrophomonas_maltophilia | escherichia_coli | 1e-10 |
| stenotrophomonas_maltophilia | klebsiella_pneumoniae | 1e-10 |
| stenotrophomonas_maltophilia | morganella_spp. | 1e-10 |
| stenotrophomonas_maltophilia | proteus_spp. | 1e-10 |
| stenotrophomonas_maltophilia | serratia_spp. | 1e-10 |
| stenotrophomonas_maltophilia | p_stuartii | 1e-10 |
| stenotrophomonas_maltophilia | pseudomonas_aeruginosa | 1e-9 |
| stenotrophomonas_maltophilia | salmonella_enterica_serovar_typhi | 1e-10 |
| stenotrophomonas_maltophilia | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| stenotrophomonas_maltophilia | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| stenotrophomonas_maltophilia | shigella_spp. | 1e-10 |
| stenotrophomonas_maltophilia | neisseria_gonorrhoeae | 3e-11 |
| stenotrophomonas_maltophilia | haemophilus_influenzae | 3e-11 |
| stenotrophomonas_maltophilia | chlamydia_trachomatis | 3e-11 |
| stenotrophomonas_maltophilia | mycoplasma_genitalium | 3e-11 |
| stenotrophomonas_maltophilia | vibrio_cholerae | 1e-10 |
| stenotrophomonas_maltophilia | neisseria_meningitidis | 3e-11 |
| stenotrophomonas_maltophilia | clostridioides_difficile | 3e-11 |
| stenotrophomonas_maltophilia | bacteroides_fragilis | 3e-11 |
| stenotrophomonas_maltophilia | enterobacter_cloacae | 1e-10 |
| stenotrophomonas_maltophilia | yersinia_enterocolitica | 1e-10 |
| stenotrophomonas_maltophilia | moraxella_catarrhalis | 3e-11 |
| stenotrophomonas_maltophilia | bordetella_pertussis | 3e-11 |
| stenotrophomonas_maltophilia | mycoplasma_pneumoniae | 3e-11 |
| stenotrophomonas_maltophilia | legionella_pneumophila | 3e-11 |
| stenotrophomonas_maltophilia | burkholderia_cepacia_complex | 1e-9 |
| staphylococcus_aureus | enterococcus_faecalis | 1e-9 |
| staphylococcus_aureus | enterococcus_faecium | 1e-9 |
| staphylococcus_aureus | staphylococcus_epidermidis | 1e-9 |
| staphylococcus_aureus | streptococcus_pneumoniae | 1e-9 |
| staphylococcus_aureus | streptococcus_pyogenes | 1e-9 |
| staphylococcus_aureus | streptococcus_agalactiae | 1e-9 |
| staphylococcus_aureus | listeria_monocytogenes | 1e-9 |
| staphylococcus_epidermidis | enterococcus_faecalis | 1e-9 |
| staphylococcus_epidermidis | enterococcus_faecium | 1e-9 |
| staphylococcus_epidermidis | staphylococcus_aureus | 1e-9 |
| staphylococcus_epidermidis | streptococcus_pneumoniae | 1e-9 |
| staphylococcus_epidermidis | streptococcus_pyogenes | 1e-9 |
| staphylococcus_epidermidis | streptococcus_agalactiae | 1e-9 |
| staphylococcus_epidermidis | listeria_monocytogenes | 1e-9 |
| streptococcus_pneumoniae | enterococcus_faecalis | 1e-9 |
| streptococcus_pneumoniae | enterococcus_faecium | 1e-9 |
| streptococcus_pneumoniae | staphylococcus_aureus | 1e-9 |
| streptococcus_pneumoniae | staphylococcus_epidermidis | 1e-9 |
| streptococcus_pneumoniae | streptococcus_pyogenes | 1e-9 |
| streptococcus_pneumoniae | streptococcus_agalactiae | 1e-9 |
| streptococcus_pneumoniae | listeria_monocytogenes | 1e-9 |
| salmonella_enterica_serovar_typhi | acinetobacter_baumannii | 1e-10 |
| salmonella_enterica_serovar_typhi | citrobacter_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | enterobacter_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | escherichia_coli | 1e-9 |
| salmonella_enterica_serovar_typhi | klebsiella_pneumoniae | 1e-9 |
| salmonella_enterica_serovar_typhi | morganella_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | proteus_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | serratia_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | p_stuartii | 1e-9 |
| salmonella_enterica_serovar_typhi | pseudomonas_aeruginosa | 1e-10 |
| salmonella_enterica_serovar_typhi | stenotrophomonas_maltophilia | 1e-10 |
| salmonella_enterica_serovar_typhi | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| salmonella_enterica_serovar_typhi | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | shigella_spp. | 1e-9 |
| salmonella_enterica_serovar_typhi | neisseria_gonorrhoeae | 3e-11 |
| salmonella_enterica_serovar_typhi | haemophilus_influenzae | 3e-11 |
| salmonella_enterica_serovar_typhi | chlamydia_trachomatis | 3e-11 |
| salmonella_enterica_serovar_typhi | mycoplasma_genitalium | 3e-11 |
| salmonella_enterica_serovar_typhi | vibrio_cholerae | 1e-10 |
| salmonella_enterica_serovar_typhi | neisseria_meningitidis | 3e-11 |
| salmonella_enterica_serovar_typhi | clostridioides_difficile | 3e-11 |
| salmonella_enterica_serovar_typhi | bacteroides_fragilis | 3e-11 |
| salmonella_enterica_serovar_typhi | enterobacter_cloacae | 1e-9 |
| salmonella_enterica_serovar_typhi | yersinia_enterocolitica | 1e-9 |
| salmonella_enterica_serovar_typhi | moraxella_catarrhalis | 3e-11 |
| salmonella_enterica_serovar_typhi | bordetella_pertussis | 3e-11 |
| salmonella_enterica_serovar_typhi | mycoplasma_pneumoniae | 3e-11 |
| salmonella_enterica_serovar_typhi | legionella_pneumophila | 3e-11 |
| salmonella_enterica_serovar_typhi | burkholderia_cepacia_complex | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | acinetobacter_baumannii | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | citrobacter_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | enterobacter_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | escherichia_coli | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | klebsiella_pneumoniae | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | morganella_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | proteus_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | serratia_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | p_stuartii | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | pseudomonas_aeruginosa | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | stenotrophomonas_maltophilia | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | salmonella_enterica_serovar_typhi | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | shigella_spp. | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | neisseria_gonorrhoeae | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | haemophilus_influenzae | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | chlamydia_trachomatis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | mycoplasma_genitalium | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | vibrio_cholerae | 1e-10 |
| salmonella_enterica_serovar_paratyphi_a | neisseria_meningitidis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | clostridioides_difficile | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | bacteroides_fragilis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | enterobacter_cloacae | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | yersinia_enterocolitica | 1e-9 |
| salmonella_enterica_serovar_paratyphi_a | moraxella_catarrhalis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | bordetella_pertussis | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | mycoplasma_pneumoniae | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | legionella_pneumophila | 3e-11 |
| salmonella_enterica_serovar_paratyphi_a | burkholderia_cepacia_complex | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | acinetobacter_baumannii | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | citrobacter_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | enterobacter_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | escherichia_coli | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | klebsiella_pneumoniae | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | morganella_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | proteus_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | serratia_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | p_stuartii | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | pseudomonas_aeruginosa | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | stenotrophomonas_maltophilia | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | shigella_spp. | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | neisseria_gonorrhoeae | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | haemophilus_influenzae | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | chlamydia_trachomatis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | mycoplasma_genitalium | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | vibrio_cholerae | 1e-10 |
| invasive_non-typhoidal_salmonella_spp. | neisseria_meningitidis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | clostridioides_difficile | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | bacteroides_fragilis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | enterobacter_cloacae | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | yersinia_enterocolitica | 1e-9 |
| invasive_non-typhoidal_salmonella_spp. | moraxella_catarrhalis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | bordetella_pertussis | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | mycoplasma_pneumoniae | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | legionella_pneumophila | 3e-11 |
| invasive_non-typhoidal_salmonella_spp. | burkholderia_cepacia_complex | 1e-10 |
| shigella_spp. | acinetobacter_baumannii | 1e-10 |
| shigella_spp. | citrobacter_spp. | 1e-9 |
| shigella_spp. | enterobacter_spp. | 1e-9 |
| shigella_spp. | escherichia_coli | 1e-9 |
| shigella_spp. | klebsiella_pneumoniae | 1e-9 |
| shigella_spp. | morganella_spp. | 1e-9 |
| shigella_spp. | proteus_spp. | 1e-9 |
| shigella_spp. | serratia_spp. | 1e-9 |
| shigella_spp. | p_stuartii | 1e-9 |
| shigella_spp. | pseudomonas_aeruginosa | 1e-10 |
| shigella_spp. | stenotrophomonas_maltophilia | 1e-10 |
| shigella_spp. | salmonella_enterica_serovar_typhi | 1e-9 |
| shigella_spp. | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| shigella_spp. | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| shigella_spp. | neisseria_gonorrhoeae | 3e-11 |
| shigella_spp. | haemophilus_influenzae | 3e-11 |
| shigella_spp. | chlamydia_trachomatis | 3e-11 |
| shigella_spp. | mycoplasma_genitalium | 3e-11 |
| shigella_spp. | vibrio_cholerae | 1e-10 |
| shigella_spp. | neisseria_meningitidis | 3e-11 |
| shigella_spp. | clostridioides_difficile | 3e-11 |
| shigella_spp. | bacteroides_fragilis | 3e-11 |
| shigella_spp. | enterobacter_cloacae | 1e-9 |
| shigella_spp. | yersinia_enterocolitica | 1e-9 |
| shigella_spp. | moraxella_catarrhalis | 3e-11 |
| shigella_spp. | bordetella_pertussis | 3e-11 |
| shigella_spp. | mycoplasma_pneumoniae | 3e-11 |
| shigella_spp. | legionella_pneumophila | 3e-11 |
| shigella_spp. | burkholderia_cepacia_complex | 1e-10 |
| neisseria_gonorrhoeae | acinetobacter_baumannii | 3e-11 |
| neisseria_gonorrhoeae | citrobacter_spp. | 3e-11 |
| neisseria_gonorrhoeae | enterobacter_spp. | 3e-11 |
| neisseria_gonorrhoeae | escherichia_coli | 3e-11 |
| neisseria_gonorrhoeae | klebsiella_pneumoniae | 3e-11 |
| neisseria_gonorrhoeae | morganella_spp. | 3e-11 |
| neisseria_gonorrhoeae | proteus_spp. | 3e-11 |
| neisseria_gonorrhoeae | serratia_spp. | 3e-11 |
| neisseria_gonorrhoeae | p_stuartii | 3e-11 |
| neisseria_gonorrhoeae | pseudomonas_aeruginosa | 3e-11 |
| neisseria_gonorrhoeae | stenotrophomonas_maltophilia | 3e-11 |
| neisseria_gonorrhoeae | salmonella_enterica_serovar_typhi | 3e-11 |
| neisseria_gonorrhoeae | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| neisseria_gonorrhoeae | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| neisseria_gonorrhoeae | shigella_spp. | 3e-11 |
| neisseria_gonorrhoeae | haemophilus_influenzae | 1e-9 |
| neisseria_gonorrhoeae | chlamydia_trachomatis | 1e-9 |
| neisseria_gonorrhoeae | mycoplasma_genitalium | 1e-9 |
| neisseria_gonorrhoeae | vibrio_cholerae | 3e-11 |
| neisseria_gonorrhoeae | neisseria_meningitidis | 1e-9 |
| neisseria_gonorrhoeae | enterobacter_cloacae | 3e-11 |
| neisseria_gonorrhoeae | yersinia_enterocolitica | 3e-11 |
| neisseria_gonorrhoeae | moraxella_catarrhalis | 1e-9 |
| neisseria_gonorrhoeae | bordetella_pertussis | 1e-9 |
| neisseria_gonorrhoeae | mycoplasma_pneumoniae | 1e-9 |
| neisseria_gonorrhoeae | legionella_pneumophila | 1e-9 |
| neisseria_gonorrhoeae | burkholderia_cepacia_complex | 3e-11 |
| streptococcus_pyogenes | enterococcus_faecalis | 1e-9 |
| streptococcus_pyogenes | enterococcus_faecium | 1e-9 |
| streptococcus_pyogenes | staphylococcus_aureus | 1e-9 |
| streptococcus_pyogenes | staphylococcus_epidermidis | 1e-9 |
| streptococcus_pyogenes | streptococcus_pneumoniae | 1e-9 |
| streptococcus_pyogenes | streptococcus_agalactiae | 1e-9 |
| streptococcus_pyogenes | listeria_monocytogenes | 1e-9 |
| streptococcus_agalactiae | enterococcus_faecalis | 1e-9 |
| streptococcus_agalactiae | enterococcus_faecium | 1e-9 |
| streptococcus_agalactiae | staphylococcus_aureus | 1e-9 |
| streptococcus_agalactiae | staphylococcus_epidermidis | 1e-9 |
| streptococcus_agalactiae | streptococcus_pneumoniae | 1e-9 |
| streptococcus_agalactiae | streptococcus_pyogenes | 1e-9 |
| streptococcus_agalactiae | listeria_monocytogenes | 1e-9 |
| haemophilus_influenzae | acinetobacter_baumannii | 3e-11 |
| haemophilus_influenzae | citrobacter_spp. | 3e-11 |
| haemophilus_influenzae | enterobacter_spp. | 3e-11 |
| haemophilus_influenzae | escherichia_coli | 3e-11 |
| haemophilus_influenzae | klebsiella_pneumoniae | 3e-11 |
| haemophilus_influenzae | morganella_spp. | 3e-11 |
| haemophilus_influenzae | proteus_spp. | 3e-11 |
| haemophilus_influenzae | serratia_spp. | 3e-11 |
| haemophilus_influenzae | p_stuartii | 3e-11 |
| haemophilus_influenzae | pseudomonas_aeruginosa | 3e-11 |
| haemophilus_influenzae | stenotrophomonas_maltophilia | 3e-11 |
| haemophilus_influenzae | salmonella_enterica_serovar_typhi | 3e-11 |
| haemophilus_influenzae | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| haemophilus_influenzae | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| haemophilus_influenzae | shigella_spp. | 3e-11 |
| haemophilus_influenzae | neisseria_gonorrhoeae | 1e-9 |
| haemophilus_influenzae | chlamydia_trachomatis | 1e-9 |
| haemophilus_influenzae | mycoplasma_genitalium | 1e-9 |
| haemophilus_influenzae | vibrio_cholerae | 3e-11 |
| haemophilus_influenzae | neisseria_meningitidis | 1e-9 |
| haemophilus_influenzae | enterobacter_cloacae | 3e-11 |
| haemophilus_influenzae | yersinia_enterocolitica | 3e-11 |
| haemophilus_influenzae | moraxella_catarrhalis | 1e-9 |
| haemophilus_influenzae | bordetella_pertussis | 1e-9 |
| haemophilus_influenzae | mycoplasma_pneumoniae | 1e-9 |
| haemophilus_influenzae | legionella_pneumophila | 1e-9 |
| haemophilus_influenzae | burkholderia_cepacia_complex | 3e-11 |
| chlamydia_trachomatis | acinetobacter_baumannii | 3e-11 |
| chlamydia_trachomatis | citrobacter_spp. | 3e-11 |
| chlamydia_trachomatis | enterobacter_spp. | 3e-11 |
| chlamydia_trachomatis | escherichia_coli | 3e-11 |
| chlamydia_trachomatis | klebsiella_pneumoniae | 3e-11 |
| chlamydia_trachomatis | morganella_spp. | 3e-11 |
| chlamydia_trachomatis | proteus_spp. | 3e-11 |
| chlamydia_trachomatis | serratia_spp. | 3e-11 |
| chlamydia_trachomatis | p_stuartii | 3e-11 |
| chlamydia_trachomatis | pseudomonas_aeruginosa | 3e-11 |
| chlamydia_trachomatis | stenotrophomonas_maltophilia | 3e-11 |
| chlamydia_trachomatis | salmonella_enterica_serovar_typhi | 3e-11 |
| chlamydia_trachomatis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| chlamydia_trachomatis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| chlamydia_trachomatis | shigella_spp. | 3e-11 |
| chlamydia_trachomatis | neisseria_gonorrhoeae | 1e-9 |
| chlamydia_trachomatis | haemophilus_influenzae | 1e-9 |
| chlamydia_trachomatis | mycoplasma_genitalium | 1e-9 |
| chlamydia_trachomatis | vibrio_cholerae | 3e-11 |
| chlamydia_trachomatis | neisseria_meningitidis | 1e-9 |
| chlamydia_trachomatis | enterobacter_cloacae | 3e-11 |
| chlamydia_trachomatis | yersinia_enterocolitica | 3e-11 |
| chlamydia_trachomatis | moraxella_catarrhalis | 1e-9 |
| chlamydia_trachomatis | bordetella_pertussis | 1e-9 |
| chlamydia_trachomatis | mycoplasma_pneumoniae | 1e-9 |
| chlamydia_trachomatis | legionella_pneumophila | 1e-9 |
| chlamydia_trachomatis | burkholderia_cepacia_complex | 3e-11 |
| mycoplasma_genitalium | acinetobacter_baumannii | 3e-11 |
| mycoplasma_genitalium | citrobacter_spp. | 3e-11 |
| mycoplasma_genitalium | enterobacter_spp. | 3e-11 |
| mycoplasma_genitalium | escherichia_coli | 3e-11 |
| mycoplasma_genitalium | klebsiella_pneumoniae | 3e-11 |
| mycoplasma_genitalium | morganella_spp. | 3e-11 |
| mycoplasma_genitalium | proteus_spp. | 3e-11 |
| mycoplasma_genitalium | serratia_spp. | 3e-11 |
| mycoplasma_genitalium | p_stuartii | 3e-11 |
| mycoplasma_genitalium | pseudomonas_aeruginosa | 3e-11 |
| mycoplasma_genitalium | stenotrophomonas_maltophilia | 3e-11 |
| mycoplasma_genitalium | salmonella_enterica_serovar_typhi | 3e-11 |
| mycoplasma_genitalium | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| mycoplasma_genitalium | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| mycoplasma_genitalium | shigella_spp. | 3e-11 |
| mycoplasma_genitalium | neisseria_gonorrhoeae | 1e-9 |
| mycoplasma_genitalium | haemophilus_influenzae | 1e-9 |
| mycoplasma_genitalium | chlamydia_trachomatis | 1e-9 |
| mycoplasma_genitalium | vibrio_cholerae | 3e-11 |
| mycoplasma_genitalium | neisseria_meningitidis | 1e-9 |
| mycoplasma_genitalium | enterobacter_cloacae | 3e-11 |
| mycoplasma_genitalium | yersinia_enterocolitica | 3e-11 |
| mycoplasma_genitalium | moraxella_catarrhalis | 1e-9 |
| mycoplasma_genitalium | bordetella_pertussis | 1e-9 |
| mycoplasma_genitalium | mycoplasma_pneumoniae | 1e-9 |
| mycoplasma_genitalium | legionella_pneumophila | 1e-9 |
| mycoplasma_genitalium | burkholderia_cepacia_complex | 3e-11 |
| vibrio_cholerae | acinetobacter_baumannii | 1e-10 |
| vibrio_cholerae | citrobacter_spp. | 1e-10 |
| vibrio_cholerae | enterobacter_spp. | 1e-10 |
| vibrio_cholerae | escherichia_coli | 1e-10 |
| vibrio_cholerae | klebsiella_pneumoniae | 1e-10 |
| vibrio_cholerae | morganella_spp. | 1e-10 |
| vibrio_cholerae | proteus_spp. | 1e-10 |
| vibrio_cholerae | serratia_spp. | 1e-10 |
| vibrio_cholerae | p_stuartii | 1e-10 |
| vibrio_cholerae | pseudomonas_aeruginosa | 1e-10 |
| vibrio_cholerae | stenotrophomonas_maltophilia | 1e-10 |
| vibrio_cholerae | salmonella_enterica_serovar_typhi | 1e-10 |
| vibrio_cholerae | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| vibrio_cholerae | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| vibrio_cholerae | shigella_spp. | 1e-10 |
| vibrio_cholerae | neisseria_gonorrhoeae | 3e-11 |
| vibrio_cholerae | haemophilus_influenzae | 3e-11 |
| vibrio_cholerae | chlamydia_trachomatis | 3e-11 |
| vibrio_cholerae | mycoplasma_genitalium | 3e-11 |
| vibrio_cholerae | neisseria_meningitidis | 3e-11 |
| vibrio_cholerae | clostridioides_difficile | 3e-11 |
| vibrio_cholerae | bacteroides_fragilis | 3e-11 |
| vibrio_cholerae | enterobacter_cloacae | 1e-10 |
| vibrio_cholerae | yersinia_enterocolitica | 1e-10 |
| vibrio_cholerae | moraxella_catarrhalis | 3e-11 |
| vibrio_cholerae | bordetella_pertussis | 3e-11 |
| vibrio_cholerae | mycoplasma_pneumoniae | 3e-11 |
| vibrio_cholerae | legionella_pneumophila | 3e-11 |
| vibrio_cholerae | burkholderia_cepacia_complex | 1e-10 |
| neisseria_meningitidis | acinetobacter_baumannii | 3e-11 |
| neisseria_meningitidis | citrobacter_spp. | 3e-11 |
| neisseria_meningitidis | enterobacter_spp. | 3e-11 |
| neisseria_meningitidis | escherichia_coli | 3e-11 |
| neisseria_meningitidis | klebsiella_pneumoniae | 3e-11 |
| neisseria_meningitidis | morganella_spp. | 3e-11 |
| neisseria_meningitidis | proteus_spp. | 3e-11 |
| neisseria_meningitidis | serratia_spp. | 3e-11 |
| neisseria_meningitidis | p_stuartii | 3e-11 |
| neisseria_meningitidis | pseudomonas_aeruginosa | 3e-11 |
| neisseria_meningitidis | stenotrophomonas_maltophilia | 3e-11 |
| neisseria_meningitidis | salmonella_enterica_serovar_typhi | 3e-11 |
| neisseria_meningitidis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| neisseria_meningitidis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| neisseria_meningitidis | shigella_spp. | 3e-11 |
| neisseria_meningitidis | neisseria_gonorrhoeae | 1e-9 |
| neisseria_meningitidis | haemophilus_influenzae | 1e-9 |
| neisseria_meningitidis | chlamydia_trachomatis | 1e-9 |
| neisseria_meningitidis | mycoplasma_genitalium | 1e-9 |
| neisseria_meningitidis | vibrio_cholerae | 3e-11 |
| neisseria_meningitidis | enterobacter_cloacae | 3e-11 |
| neisseria_meningitidis | yersinia_enterocolitica | 3e-11 |
| neisseria_meningitidis | moraxella_catarrhalis | 1e-9 |
| neisseria_meningitidis | bordetella_pertussis | 1e-9 |
| neisseria_meningitidis | mycoplasma_pneumoniae | 1e-9 |
| neisseria_meningitidis | legionella_pneumophila | 1e-9 |
| neisseria_meningitidis | burkholderia_cepacia_complex | 3e-11 |
| listeria_monocytogenes | enterococcus_faecalis | 1e-9 |
| listeria_monocytogenes | enterococcus_faecium | 1e-9 |
| listeria_monocytogenes | staphylococcus_aureus | 1e-9 |
| listeria_monocytogenes | staphylococcus_epidermidis | 1e-9 |
| listeria_monocytogenes | streptococcus_pneumoniae | 1e-9 |
| listeria_monocytogenes | streptococcus_pyogenes | 1e-9 |
| listeria_monocytogenes | streptococcus_agalactiae | 1e-9 |
| clostridioides_difficile | acinetobacter_baumannii | 3e-11 |
| clostridioides_difficile | citrobacter_spp. | 3e-11 |
| clostridioides_difficile | enterobacter_spp. | 3e-11 |
| clostridioides_difficile | escherichia_coli | 3e-11 |
| clostridioides_difficile | klebsiella_pneumoniae | 3e-11 |
| clostridioides_difficile | morganella_spp. | 3e-11 |
| clostridioides_difficile | proteus_spp. | 3e-11 |
| clostridioides_difficile | serratia_spp. | 3e-11 |
| clostridioides_difficile | p_stuartii | 3e-11 |
| clostridioides_difficile | pseudomonas_aeruginosa | 3e-11 |
| clostridioides_difficile | stenotrophomonas_maltophilia | 3e-11 |
| clostridioides_difficile | salmonella_enterica_serovar_typhi | 3e-11 |
| clostridioides_difficile | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| clostridioides_difficile | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| clostridioides_difficile | shigella_spp. | 3e-11 |
| clostridioides_difficile | vibrio_cholerae | 3e-11 |
| clostridioides_difficile | bacteroides_fragilis | 1e-9 |
| clostridioides_difficile | enterobacter_cloacae | 3e-11 |
| clostridioides_difficile | yersinia_enterocolitica | 3e-11 |
| clostridioides_difficile | burkholderia_cepacia_complex | 3e-11 |
| bacteroides_fragilis | acinetobacter_baumannii | 3e-11 |
| bacteroides_fragilis | citrobacter_spp. | 3e-11 |
| bacteroides_fragilis | enterobacter_spp. | 3e-11 |
| bacteroides_fragilis | escherichia_coli | 3e-11 |
| bacteroides_fragilis | klebsiella_pneumoniae | 3e-11 |
| bacteroides_fragilis | morganella_spp. | 3e-11 |
| bacteroides_fragilis | proteus_spp. | 3e-11 |
| bacteroides_fragilis | serratia_spp. | 3e-11 |
| bacteroides_fragilis | p_stuartii | 3e-11 |
| bacteroides_fragilis | pseudomonas_aeruginosa | 3e-11 |
| bacteroides_fragilis | stenotrophomonas_maltophilia | 3e-11 |
| bacteroides_fragilis | salmonella_enterica_serovar_typhi | 3e-11 |
| bacteroides_fragilis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| bacteroides_fragilis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| bacteroides_fragilis | shigella_spp. | 3e-11 |
| bacteroides_fragilis | vibrio_cholerae | 3e-11 |
| bacteroides_fragilis | clostridioides_difficile | 1e-9 |
| bacteroides_fragilis | enterobacter_cloacae | 3e-11 |
| bacteroides_fragilis | yersinia_enterocolitica | 3e-11 |
| bacteroides_fragilis | burkholderia_cepacia_complex | 3e-11 |
| enterobacter_cloacae | acinetobacter_baumannii | 1e-10 |
| enterobacter_cloacae | citrobacter_spp. | 1e-9 |
| enterobacter_cloacae | enterobacter_spp. | 1e-9 |
| enterobacter_cloacae | escherichia_coli | 1e-9 |
| enterobacter_cloacae | klebsiella_pneumoniae | 1e-9 |
| enterobacter_cloacae | morganella_spp. | 1e-9 |
| enterobacter_cloacae | proteus_spp. | 1e-9 |
| enterobacter_cloacae | serratia_spp. | 1e-9 |
| enterobacter_cloacae | p_stuartii | 1e-9 |
| enterobacter_cloacae | pseudomonas_aeruginosa | 1e-10 |
| enterobacter_cloacae | stenotrophomonas_maltophilia | 1e-10 |
| enterobacter_cloacae | salmonella_enterica_serovar_typhi | 1e-9 |
| enterobacter_cloacae | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| enterobacter_cloacae | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| enterobacter_cloacae | shigella_spp. | 1e-9 |
| enterobacter_cloacae | neisseria_gonorrhoeae | 3e-11 |
| enterobacter_cloacae | haemophilus_influenzae | 3e-11 |
| enterobacter_cloacae | chlamydia_trachomatis | 3e-11 |
| enterobacter_cloacae | mycoplasma_genitalium | 3e-11 |
| enterobacter_cloacae | vibrio_cholerae | 1e-10 |
| enterobacter_cloacae | neisseria_meningitidis | 3e-11 |
| enterobacter_cloacae | clostridioides_difficile | 3e-11 |
| enterobacter_cloacae | bacteroides_fragilis | 3e-11 |
| enterobacter_cloacae | yersinia_enterocolitica | 1e-9 |
| enterobacter_cloacae | moraxella_catarrhalis | 3e-11 |
| enterobacter_cloacae | bordetella_pertussis | 3e-11 |
| enterobacter_cloacae | mycoplasma_pneumoniae | 3e-11 |
| enterobacter_cloacae | legionella_pneumophila | 3e-11 |
| enterobacter_cloacae | burkholderia_cepacia_complex | 1e-10 |
| yersinia_enterocolitica | acinetobacter_baumannii | 1e-10 |
| yersinia_enterocolitica | citrobacter_spp. | 1e-9 |
| yersinia_enterocolitica | enterobacter_spp. | 1e-9 |
| yersinia_enterocolitica | escherichia_coli | 1e-9 |
| yersinia_enterocolitica | klebsiella_pneumoniae | 1e-9 |
| yersinia_enterocolitica | morganella_spp. | 1e-9 |
| yersinia_enterocolitica | proteus_spp. | 1e-9 |
| yersinia_enterocolitica | serratia_spp. | 1e-9 |
| yersinia_enterocolitica | p_stuartii | 1e-9 |
| yersinia_enterocolitica | pseudomonas_aeruginosa | 1e-10 |
| yersinia_enterocolitica | stenotrophomonas_maltophilia | 1e-10 |
| yersinia_enterocolitica | salmonella_enterica_serovar_typhi | 1e-9 |
| yersinia_enterocolitica | salmonella_enterica_serovar_paratyphi_a | 1e-9 |
| yersinia_enterocolitica | invasive_non-typhoidal_salmonella_spp. | 1e-9 |
| yersinia_enterocolitica | shigella_spp. | 1e-9 |
| yersinia_enterocolitica | neisseria_gonorrhoeae | 3e-11 |
| yersinia_enterocolitica | haemophilus_influenzae | 3e-11 |
| yersinia_enterocolitica | chlamydia_trachomatis | 3e-11 |
| yersinia_enterocolitica | mycoplasma_genitalium | 3e-11 |
| yersinia_enterocolitica | vibrio_cholerae | 1e-10 |
| yersinia_enterocolitica | neisseria_meningitidis | 3e-11 |
| yersinia_enterocolitica | clostridioides_difficile | 3e-11 |
| yersinia_enterocolitica | bacteroides_fragilis | 3e-11 |
| yersinia_enterocolitica | enterobacter_cloacae | 1e-9 |
| yersinia_enterocolitica | moraxella_catarrhalis | 3e-11 |
| yersinia_enterocolitica | bordetella_pertussis | 3e-11 |
| yersinia_enterocolitica | mycoplasma_pneumoniae | 3e-11 |
| yersinia_enterocolitica | legionella_pneumophila | 3e-11 |
| yersinia_enterocolitica | burkholderia_cepacia_complex | 1e-10 |
| moraxella_catarrhalis | acinetobacter_baumannii | 3e-11 |
| moraxella_catarrhalis | citrobacter_spp. | 3e-11 |
| moraxella_catarrhalis | enterobacter_spp. | 3e-11 |
| moraxella_catarrhalis | escherichia_coli | 3e-11 |
| moraxella_catarrhalis | klebsiella_pneumoniae | 3e-11 |
| moraxella_catarrhalis | morganella_spp. | 3e-11 |
| moraxella_catarrhalis | proteus_spp. | 3e-11 |
| moraxella_catarrhalis | serratia_spp. | 3e-11 |
| moraxella_catarrhalis | p_stuartii | 3e-11 |
| moraxella_catarrhalis | pseudomonas_aeruginosa | 3e-11 |
| moraxella_catarrhalis | stenotrophomonas_maltophilia | 3e-11 |
| moraxella_catarrhalis | salmonella_enterica_serovar_typhi | 3e-11 |
| moraxella_catarrhalis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| moraxella_catarrhalis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| moraxella_catarrhalis | shigella_spp. | 3e-11 |
| moraxella_catarrhalis | neisseria_gonorrhoeae | 1e-9 |
| moraxella_catarrhalis | haemophilus_influenzae | 1e-9 |
| moraxella_catarrhalis | chlamydia_trachomatis | 1e-9 |
| moraxella_catarrhalis | mycoplasma_genitalium | 1e-9 |
| moraxella_catarrhalis | vibrio_cholerae | 3e-11 |
| moraxella_catarrhalis | neisseria_meningitidis | 1e-9 |
| moraxella_catarrhalis | enterobacter_cloacae | 3e-11 |
| moraxella_catarrhalis | yersinia_enterocolitica | 3e-11 |
| moraxella_catarrhalis | bordetella_pertussis | 1e-9 |
| moraxella_catarrhalis | mycoplasma_pneumoniae | 1e-9 |
| moraxella_catarrhalis | legionella_pneumophila | 1e-9 |
| moraxella_catarrhalis | burkholderia_cepacia_complex | 3e-11 |
| bordetella_pertussis | acinetobacter_baumannii | 3e-11 |
| bordetella_pertussis | citrobacter_spp. | 3e-11 |
| bordetella_pertussis | enterobacter_spp. | 3e-11 |
| bordetella_pertussis | escherichia_coli | 3e-11 |
| bordetella_pertussis | klebsiella_pneumoniae | 3e-11 |
| bordetella_pertussis | morganella_spp. | 3e-11 |
| bordetella_pertussis | proteus_spp. | 3e-11 |
| bordetella_pertussis | serratia_spp. | 3e-11 |
| bordetella_pertussis | p_stuartii | 3e-11 |
| bordetella_pertussis | pseudomonas_aeruginosa | 3e-11 |
| bordetella_pertussis | stenotrophomonas_maltophilia | 3e-11 |
| bordetella_pertussis | salmonella_enterica_serovar_typhi | 3e-11 |
| bordetella_pertussis | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| bordetella_pertussis | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| bordetella_pertussis | shigella_spp. | 3e-11 |
| bordetella_pertussis | neisseria_gonorrhoeae | 1e-9 |
| bordetella_pertussis | haemophilus_influenzae | 1e-9 |
| bordetella_pertussis | chlamydia_trachomatis | 1e-9 |
| bordetella_pertussis | mycoplasma_genitalium | 1e-9 |
| bordetella_pertussis | vibrio_cholerae | 3e-11 |
| bordetella_pertussis | neisseria_meningitidis | 1e-9 |
| bordetella_pertussis | enterobacter_cloacae | 3e-11 |
| bordetella_pertussis | yersinia_enterocolitica | 3e-11 |
| bordetella_pertussis | moraxella_catarrhalis | 1e-9 |
| bordetella_pertussis | mycoplasma_pneumoniae | 1e-9 |
| bordetella_pertussis | legionella_pneumophila | 1e-9 |
| bordetella_pertussis | burkholderia_cepacia_complex | 3e-11 |
| mycoplasma_pneumoniae | acinetobacter_baumannii | 3e-11 |
| mycoplasma_pneumoniae | citrobacter_spp. | 3e-11 |
| mycoplasma_pneumoniae | enterobacter_spp. | 3e-11 |
| mycoplasma_pneumoniae | escherichia_coli | 3e-11 |
| mycoplasma_pneumoniae | klebsiella_pneumoniae | 3e-11 |
| mycoplasma_pneumoniae | morganella_spp. | 3e-11 |
| mycoplasma_pneumoniae | proteus_spp. | 3e-11 |
| mycoplasma_pneumoniae | serratia_spp. | 3e-11 |
| mycoplasma_pneumoniae | p_stuartii | 3e-11 |
| mycoplasma_pneumoniae | pseudomonas_aeruginosa | 3e-11 |
| mycoplasma_pneumoniae | stenotrophomonas_maltophilia | 3e-11 |
| mycoplasma_pneumoniae | salmonella_enterica_serovar_typhi | 3e-11 |
| mycoplasma_pneumoniae | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| mycoplasma_pneumoniae | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| mycoplasma_pneumoniae | shigella_spp. | 3e-11 |
| mycoplasma_pneumoniae | neisseria_gonorrhoeae | 1e-9 |
| mycoplasma_pneumoniae | haemophilus_influenzae | 1e-9 |
| mycoplasma_pneumoniae | chlamydia_trachomatis | 1e-9 |
| mycoplasma_pneumoniae | mycoplasma_genitalium | 1e-9 |
| mycoplasma_pneumoniae | vibrio_cholerae | 3e-11 |
| mycoplasma_pneumoniae | neisseria_meningitidis | 1e-9 |
| mycoplasma_pneumoniae | enterobacter_cloacae | 3e-11 |
| mycoplasma_pneumoniae | yersinia_enterocolitica | 3e-11 |
| mycoplasma_pneumoniae | moraxella_catarrhalis | 1e-9 |
| mycoplasma_pneumoniae | bordetella_pertussis | 1e-9 |
| mycoplasma_pneumoniae | legionella_pneumophila | 1e-9 |
| mycoplasma_pneumoniae | burkholderia_cepacia_complex | 3e-11 |
| legionella_pneumophila | acinetobacter_baumannii | 3e-11 |
| legionella_pneumophila | citrobacter_spp. | 3e-11 |
| legionella_pneumophila | enterobacter_spp. | 3e-11 |
| legionella_pneumophila | escherichia_coli | 3e-11 |
| legionella_pneumophila | klebsiella_pneumoniae | 3e-11 |
| legionella_pneumophila | morganella_spp. | 3e-11 |
| legionella_pneumophila | proteus_spp. | 3e-11 |
| legionella_pneumophila | serratia_spp. | 3e-11 |
| legionella_pneumophila | p_stuartii | 3e-11 |
| legionella_pneumophila | pseudomonas_aeruginosa | 3e-11 |
| legionella_pneumophila | stenotrophomonas_maltophilia | 3e-11 |
| legionella_pneumophila | salmonella_enterica_serovar_typhi | 3e-11 |
| legionella_pneumophila | salmonella_enterica_serovar_paratyphi_a | 3e-11 |
| legionella_pneumophila | invasive_non-typhoidal_salmonella_spp. | 3e-11 |
| legionella_pneumophila | shigella_spp. | 3e-11 |
| legionella_pneumophila | neisseria_gonorrhoeae | 1e-9 |
| legionella_pneumophila | haemophilus_influenzae | 1e-9 |
| legionella_pneumophila | chlamydia_trachomatis | 1e-9 |
| legionella_pneumophila | mycoplasma_genitalium | 1e-9 |
| legionella_pneumophila | vibrio_cholerae | 3e-11 |
| legionella_pneumophila | neisseria_meningitidis | 1e-9 |
| legionella_pneumophila | enterobacter_cloacae | 3e-11 |
| legionella_pneumophila | yersinia_enterocolitica | 3e-11 |
| legionella_pneumophila | moraxella_catarrhalis | 1e-9 |
| legionella_pneumophila | bordetella_pertussis | 1e-9 |
| legionella_pneumophila | mycoplasma_pneumoniae | 1e-9 |
| legionella_pneumophila | burkholderia_cepacia_complex | 3e-11 |
| burkholderia_cepacia_complex | acinetobacter_baumannii | 1e-9 |
| burkholderia_cepacia_complex | citrobacter_spp. | 1e-10 |
| burkholderia_cepacia_complex | enterobacter_spp. | 1e-10 |
| burkholderia_cepacia_complex | escherichia_coli | 1e-10 |
| burkholderia_cepacia_complex | klebsiella_pneumoniae | 1e-10 |
| burkholderia_cepacia_complex | morganella_spp. | 1e-10 |
| burkholderia_cepacia_complex | proteus_spp. | 1e-10 |
| burkholderia_cepacia_complex | serratia_spp. | 1e-10 |
| burkholderia_cepacia_complex | p_stuartii | 1e-10 |
| burkholderia_cepacia_complex | pseudomonas_aeruginosa | 1e-9 |
| burkholderia_cepacia_complex | stenotrophomonas_maltophilia | 1e-9 |
| burkholderia_cepacia_complex | salmonella_enterica_serovar_typhi | 1e-10 |
| burkholderia_cepacia_complex | salmonella_enterica_serovar_paratyphi_a | 1e-10 |
| burkholderia_cepacia_complex | invasive_non-typhoidal_salmonella_spp. | 1e-10 |
| burkholderia_cepacia_complex | shigella_spp. | 1e-10 |
| burkholderia_cepacia_complex | neisseria_gonorrhoeae | 3e-11 |
| burkholderia_cepacia_complex | haemophilus_influenzae | 3e-11 |
| burkholderia_cepacia_complex | chlamydia_trachomatis | 3e-11 |
| burkholderia_cepacia_complex | mycoplasma_genitalium | 3e-11 |
| burkholderia_cepacia_complex | vibrio_cholerae | 1e-10 |
| burkholderia_cepacia_complex | neisseria_meningitidis | 3e-11 |
| burkholderia_cepacia_complex | clostridioides_difficile | 3e-11 |
| burkholderia_cepacia_complex | bacteroides_fragilis | 3e-11 |
| burkholderia_cepacia_complex | enterobacter_cloacae | 1e-10 |
| burkholderia_cepacia_complex | yersinia_enterocolitica | 1e-10 |
| burkholderia_cepacia_complex | moraxella_catarrhalis | 3e-11 |
| burkholderia_cepacia_complex | bordetella_pertussis | 3e-11 |
| burkholderia_cepacia_complex | mycoplasma_pneumoniae | 3e-11 |
| burkholderia_cepacia_complex | legionella_pneumophila | 3e-11 |

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



#### Summary numeric columns (per timestep)

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



#### Per-drug columns (~62 each)

| Pattern | Description |
|---------|-------------|
| `{drug}_prescribed` | Courses initiated this day |
| `{drug}_active_treatments` | Currently on this drug |



#### Per-bacteria × per-drug columns

| Pattern | Description |
|---------|-------------|
| `{bacteria}_{drug}_activity_r` and related summaries | Mean resistance/activity summaries derived from mechanism flags |
| `{bacteria}_{drug}_infected_and_standardized_mic_lt2` and related summaries | Treatment-relevant susceptible-infection summaries |
| `{bacteria}_{drug}_currently_on_drug` and related summaries | Drug exposure by organism-drug pair |
| `{bacteria}_{drug}_any_r`, `microbiome_r_positive`, and split hospital/community variants | Mechanism-derived resistance summaries. The summary CSV no longer exports `{bacteria}_{drug}_majority_r`; mechanism-majority state is used internally and exposed through mechanism-derived summaries rather than the old majority-r column family. |



#### Per-region columns (~6 each)

| Pattern | Description |
|---------|-------------|
| `{region}_infected` | Regional infection count |
| `{region}_hospitalized` | Regional hospital count |
| `{region}_deaths` | Regional death count |



### C.4 Total Column Count

With 42 bacteria, 62 drugs, 46 resistance mechanisms, and 6 regions, the CSV contains thousands of columns. The exact schema is generated by `Simulation::export_summary_to_csv` and guarded by `tests/csv_invariants.rs`; it can change when output families are added or retired.

Approximate families:

- summary numeric columns for population state, mortality, treatment, policy branch, and calendar time
- per-bacteria columns for infection, carriage, acquisition route, deaths, and sepsis
- per-drug columns for starts and active treatment
- per-bacteria-per-drug columns for activity, MIC, treatment exposure, susceptible-infection counts, and mechanism-derived resistance summaries
- per-region columns for infection, hospitalization, and death summaries



### C.5 Infection Journey Logs

When enabled, individual infection journeys are logged to the `infection_journeys/` directory as CSV files, capturing:

- Infection acquisition details
- Resistance profile at acquisition and over time
- Treatment episodes
- Clinical outcome (clearance, death, ongoing)
- Mechanism gains and losses

---

*This document describes the model as implemented in the Rust codebase. All variable names correspond to parameter keys used in `src/config.rs`.*

---

## References

- Ali M, Nelson AR, Lopez AL, Sack DA. Updated global burden of cholera in endemic countries. *PLoS Negl Trop Dis.* 2015;9(6):e0003832. doi:10.1371/journal.pntd.0003832

- Andersson DI, Hughes D. Antibiotic resistance and its cost: is it possible to reverse resistance? *Nat Rev Microbiol.* 2010;8(4):260–271. doi:10.1038/nrmicro2319

- Arcilla MS, van Hattem JM, Haverkate MR, et al. Import and spread of extended-spectrum β-lactamase-producing Enterobacteriaceae by international travellers (COMBAT study): a prospective, multicentre cohort study. *Lancet Infect Dis.* 2017;17(1):78–85. doi:10.1016/S1473-3099(16)30319-X

- Barlam TF, Cosgrove SE, Abbo LM, et al. Implementing an antibiotic stewardship program: guidelines by the Infectious Diseases Society of America and the Society for Healthcare Epidemiology of America. *Clin Infect Dis.* 2016;62(10):e51–e77. doi:10.1093/cid/ciw118

- Bassetti M, Vena A, Croxatto A, Righi E, Guery B. How to manage *Pseudomonas aeruginosa* infections. *Drugs Context.* 2018;7:212527. doi:10.7573/dic.212527

- Bauer AW, Kirby WMM, Sherris JC, Turck M. Antibiotic susceptibility testing by a standardized single disk method. *Am J Clin Pathol.* 1966;45(4):493–496. doi:10.1093/ajcp/45.4_ts.493

- Beaber JW, Hochhut B, Waldor MK. SOS response promotes horizontal dissemination of antibiotic resistance genes. *Nature.* 2004;427(6969):72–74. doi:10.1038/nature02241

- Blaser MJ. Epidemiologic and clinical features of *Campylobacter jejuni* infections. *J Infect Dis.* 1997;176(Suppl 2):S103–S105.

- Bogaert D, De Groot R, Hermans PWM. *Streptococcus pneumoniae* colonisation: the key to pneumococcal disease. *Lancet Infect Dis.* 2004;4(3):144–154. doi:10.1016/S1473-3099(04)00938-7

- Borger AL, Abarca AA, Dötsch A, et al. Mobile resistance genes in *Mycobacterium tuberculosis*: current evidence and future perspectives. *Lancet Infect Dis.* 2023;23(7):e268–e278. doi:10.1016/S1473-3099(22)00785-0

- Bratzler DW, Dellinger EP, Olsen KM, et al. Clinical practice guidelines for antimicrobial prophylaxis in surgery. *Am J Health-Syst Pharm.* 2013;70(3):195–283. doi:10.2146/ajhp120568

- Brooke JS. *Stenotrophomonas maltophilia*: an emerging global opportunistic pathogen. *Clin Microbiol Rev.* 2012;25(1):2–41. doi:10.1128/CMR.00019-11

- Browne AJ, Chipeta MG, Haines-Woodhouse G, et al. Global antibiotic consumption and usage in humans, 2000–18: a spatial modelling study. *Lancet Planet Health.* 2021;5(12):e893–e904. doi:10.1016/S2542-5196(21)00280-1

- Brunton LL, Hilal-Dandan R, Knollmann BC, eds. *Goodman & Gilman's: The Pharmacological Basis of Therapeutics.* 13th ed. New York: McGraw-Hill; 2018.

- Buelow E, Gonzalez TB, Versluis D, et al. Effects of selective digestive decontamination on the human gut microbiome and resistome as revealed by a large-scale longitudinal metagenomic study. *Microbiome.* 2017;5(1):154. doi:10.1186/s40168-017-0369-0

- Carapetis JR, Steer AC, Mulholland EK, Weber M. The global burden of group A streptococcal diseases. *Lancet Infect Dis.* 2005;5(11):685–694. doi:10.1016/S1473-3099(05)70267-X

- Carattoli A. Resistance plasmid families in Enterobacteriaceae. *Antimicrob Agents Chemother.* 2009;53(6):2227–2238. doi:10.1128/AAC.01707-08

- Centers for Disease Control and Prevention. *Antibiotic Resistance Threats in the United States, 2019.* Atlanta, GA: U.S. Department of Health and Human Services, CDC; 2019. https://www.cdc.gov/antimicrobial-resistance/data-research/threats/index.html

- Centers for Disease Control and Prevention. *COVID-19: U.S. Impact on Antimicrobial Resistance, Special Report 2022.* Atlanta, GA: U.S. Department of Health and Human Services, CDC; 2022. https://www.cdc.gov/antimicrobial-resistance/data-research/threats/covid-19.html

- Crossman LC, Gould VC, Dow JM, et al. The complete genome, comparative and functional analysis of *Stenotrophomonas maltophilia* reveals an organism heavily shielded by drug resistance determinants. *Genome Biol.* 2008;9(4):R74. doi:10.1186/gb-2008-9-4-r74

- Davey P, Marwick CA, Scott CL, et al. Interventions to improve antibiotic prescribing practices for hospital inpatients. *Cochrane Database Syst Rev.* 2017;(2):CD003543. doi:10.1002/14651858.CD003543.pub4

- Drlica K, Zhao X. Mutant selection window hypothesis updated. *Clin Infect Dis.* 2007;44(5):681–688. doi:10.1086/511025

- Dunne MW, Puttagunta S, Giordano P, Krievins D, Zelasky M, Baldassarre J. A randomized clinical trial of single-dose versus weekly dalbavancin for treatment of acute bacterial skin and skin structure infection. *Clin Infect Dis.* 2016;62(5):545–551. doi:10.1093/cid/ciw005

- European Centre for Disease Prevention and Control. *Antimicrobial consumption in the EU/EEA (ESAC-Net) - Annual Epidemiological Report 2023.* Stockholm: ECDC; 2024. https://www.ecdc.europa.eu/en/publications-data/antimicrobial-consumption-eueea-esac-net-annual-epidemiological-report-2023

- European Centre for Disease Prevention and Control. *Antimicrobial resistance in the EU/EEA (EARS-Net) - Annual Epidemiological Report 2023.* Stockholm: ECDC; 2024. https://www.ecdc.europa.eu/en/publications-data/antimicrobial-resistance-eueea-ears-net-annual-epidemiological-report-2023

- European Committee on Antimicrobial Susceptibility Testing (EUCAST). *Breakpoint tables for interpretation of MICs and zone diameters.* Version 13.0; 2023. https://www.eucast.org/clinical_breakpoints

- Evans L, Rhodes A, Alhazzani W, et al. Surviving sepsis campaign: international guidelines for management of sepsis and septic shock 2021. *Intensive Care Med.* 2021;47(11):1181–1247. doi:10.1007/s00134-021-06506-y

- Fishman JA. Infection in solid-organ transplant recipients. *N Engl J Med.* 2007;357(25):2601–2614. doi:10.1056/NEJMra064928

- Fleming-Dutra KE, Hersh AL, Shapiro DJ, et al. Prevalence of inappropriate antibiotic prescriptions among US ambulatory care visits, 2010–2011. *JAMA.* 2016;315(17):1864–1873. doi:10.1001/jama.2016.4151

- Forslund K, Sunagawa S, Kultima JR, et al. Country-specific antibiotic use practices impact the human gut resistome. *Genome Res.* 2013;23(7):1163–1169. doi:10.1101/gr.155465.113

- GBD 2019 Lower Respiratory Infections Collaborators. Age-sex differences in the global burden of lower respiratory infections and risk factors, 1990–2019: results from the Global Burden of Disease Study 2019. *Lancet Infect Dis.* 2022;22(11):1626–1647. doi:10.1016/S1473-3099(22)00510-2

- GBD 2021 Antimicrobial Resistance Collaborators. Global burden of bacterial antimicrobial resistance 1990–2021: a systematic analysis with forecasts to 2050. *Lancet.* 2024;404(10459):1199–1226. doi:10.1016/S0140-6736(24)01867-1

- Gorrie CL, Mirčeta M, Wick RR, et al. Gastrointestinal carriage is a major reservoir of *Klebsiella pneumoniae* infection in intensive care patients. *Clin Infect Dis.* 2017;65(2):208–215. doi:10.1093/cid/cix270

- Guh AY, Mu Y, Winston LG, et al. Trends in U.S. burden of *Clostridioides difficile* infection and outcomes. *N Engl J Med.* 2020;382(14):1320–1330. doi:10.1056/NEJMoa1910215

- Gupta K, Hooton TM, Naber KG, et al. International clinical practice guidelines for the treatment of acute uncomplicated cystitis and pyelonephritis in women: a 2010 update by the Infectious Diseases Society of America and the European Society for Microbiology and Infectious Diseases. *Clin Infect Dis.* 2011;52(5):e103–e120. doi:10.1093/cid/ciq257

- Havelaar AH, Kirk MD, Torgerson PR, et al. World Health Organization global estimates and regional comparisons of the burden of foodborne disease in 2010. *PLoS Med.* 2015;12(12):e1001923. doi:10.1371/journal.pmed.1001923

- Hooi JKY, Lai WY, Ng WK, et al. Global prevalence of *Helicobacter pylori* infection: systematic review and meta-analysis. *Gastroenterology.* 2017;153(2):420–429. doi:10.1053/j.gastro.2017.04.022

- Human Microbiome Project Consortium. Structure, function and diversity of the healthy human microbiome. *Nature.* 2012;486(7402):207–214. doi:10.1038/nature11234

- Ikuta KS, Swetschinski LR, Robles Aguilar G, et al. Global mortality associated with 33 bacterial pathogens in 2019: a systematic analysis for the Global Burden of Disease Study 2019. *Lancet.* 2022;400(10369):2221–2248. doi:10.1016/S0140-6736(22)02185-7

- International Organization for Standardization. *ISO 20776-2:2021 Clinical laboratory testing and in vitro diagnostic test systems - Susceptibility testing of infectious agents and evaluation of performance of antimicrobial susceptibility test devices - Part 2: Evaluation of performance of antimicrobial susceptibility test devices against reference broth micro-dilution.* Geneva: ISO; 2021.

- Jacobs J, Hardy L, Semret M, et al. Diagnostic bacteriology in district hospitals in sub-Saharan Africa: at the forefront of the containment of antimicrobial resistance. *Front Med (Lausanne).* 2019;6:205. doi:10.3389/fmed.2019.00205

- Klein EY, Van Boeckel TP, Martinez EM, et al. Global increase and geographic convergence in antibiotic consumption between 2000 and 2015. *Proc Natl Acad Sci USA.* 2018;115(15):E3463–E3470. doi:10.1073/pnas.1717295115

- Koning S, van der Sande R, Verhagen AP, et al. Interventions for impetigo. *Cochrane Database Syst Rev.* 2012;(1):CD003261. doi:10.1002/14651858.CD003261.pub3

- Korenromp EL, Rowley J, Alonso M, et al. Global burden of maternal and congenital syphilis and associated adverse birth outcomes — Estimates for 2016 and progress since 2012. *PLoS One.* 2019;14(2):e0211720. doi:10.1371/journal.pone.0211720

- Lee CF, Cowling BJ, Feng S, et al. Impact of antibiotic stewardship programmes in Asia: a systematic review and meta-analysis. *J Antimicrob Chemother.* 2018;73(4):844–851. doi:10.1093/jac/dkx492

- Levy MM, Dellinger RP, Townsend SR, et al. The Surviving Sepsis Campaign: results of an international guideline-based performance improvement program targeting severe sepsis. *Intensive Care Med.* 2010;36(2):222–231. doi:10.1007/s00134-009-1738-3

- Li J, Nation RL, Turnidge JD, et al. Colistin: the re-emerging antibiotic for multidrug-resistant Gram-negative bacterial infections. *Lancet Infect Dis.* 2006;6(9):589–601. doi:10.1016/S1473-3099(06)70580-1

- Llewelyn MJ, Fitzpatrick JM, Darwin E, et al. The antibiotic course has had its day. *BMJ.* 2017;358:j3418. doi:10.1136/bmj.j3418

- Magill SS, O'Leary E, Janelle SJ, et al. Changes in prevalence of health care–associated infections in U.S. hospitals. *N Engl J Med.* 2018;379(18):1732–1744. doi:10.1056/NEJMoa1801550

- Martinson ML, Lapham J. Prevalence of immunosuppression among US adults. *JAMA.* 2024;331(10):880–882. doi:10.1001/jama.2023.28019

- McInnes RS, McCallum GE, Lamberte LE, van Schaik W. Horizontal transfer of antibiotic resistance genes in the human gut microbiome. *Curr Opin Microbiol.* 2020;53:35–43. doi:10.1016/j.mib.2020.02.002

- Metlay JP, Waterer GW, Long AC, et al. Diagnosis and treatment of adults with community-acquired pneumonia: an official clinical practice guideline of the American Thoracic Society and Infectious Diseases Society of America. *Am J Respir Crit Care Med.* 2019;200(7):e45–e67. doi:10.1164/rccm.201908-1581ST

- Murray CJL, Ikuta KS, Sharara F, et al. Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *Lancet.* 2022;399(10325):629–655. doi:10.1016/S0140-6736(21)02724-0

- Nielsen EI, Friberg LE. Pharmacokinetic-pharmacodynamic modeling of antibacterial drugs. *Pharmacol Rev.* 2013;65(3):1053–1090. doi:10.1124/pr.111.005769

- Partridge SR, Kwong SM, Firth N, Jensen SO. Mobile genetic elements associated with antimicrobial resistance. *Clin Microbiol Rev.* 2018;31(4):e00088-17. doi:10.1128/CMR.00088-17

- Pitt TL, Batchelor BI. Antimicrobial susceptibility testing. In: Greenwood D, Barer M, Slack R, Irving W, eds. *Medical Microbiology.* 19th ed. Edinburgh: Churchill Livingstone; 2019.

- Plummer M, Franceschi S, Vignat J, Forman D, de Martel C. Global burden of gastric cancer attributable to *Helicobacter pylori*. *Int J Cancer.* 2015;136(2):487–490. doi:10.1002/ijc.28999

- Poolman JT, Wacker M. Extraintestinal pathogenic *Escherichia coli*, a common human pathogen: challenges for vaccine development and progress in the field. *J Infect Dis.* 2016;213(1):6–13. doi:10.1093/infdis/jiv429

- Rhodes A, Evans LE, Alhazzani W, et al. Surviving Sepsis Campaign: international guidelines for management of sepsis and septic shock: 2016. *Intensive Care Med.* 2017;43(3):304–377. doi:10.1007/s00134-017-4683-6

- Rowley J, Vander Hoorn S, Korenromp EL, et al. Chlamydia, gonorrhoea, trichomoniasis and syphilis: global prevalence and incidence estimates, 2016. *Bull World Health Organ.* 2019;97(8):548–562P. doi:10.2471/BLT.18.228486

- Rudd KE, Johnson SC, Agesa KM, et al. Global, regional, and national sepsis incidence and mortality, 1990–2017: analysis for the Global Burden of Disease Study. *Lancet.* 2020;395(10219):200–211. doi:10.1016/S0140-6736(19)32989-7

- Rybak MJ, Le J, Lodise TP, et al. Therapeutic monitoring of vancomycin for serious methicillin-resistant *Staphylococcus aureus* infections: A revised consensus guideline and review by the American Society of Health-System Pharmacists, the Infectious Diseases Society of America, and the Society of Infectious Diseases Pharmacists. *Am J Health-Syst Pharm.* 2020;77(11):835–864. doi:10.1093/ajhp/zxaa036

- San Millán A, MacLean RC. Fitness costs of plasmids: a limit to plasmid transmission. *Microbiol Spectr.* 2017;5(5):MTBP-0016-2017. doi:10.1128/microbiolspec.MTBP-0016-2017

- Savoldi A, Carrara E, Graham DY, Conti M, Tacconelli E. Prevalence of antibiotic resistance in *Helicobacter pylori*: a systematic review and meta-analysis in World Health Organization regions. *Gastroenterology.* 2018;155(5):1372–1382.e17. doi:10.1053/j.gastro.2018.07.022

- Schuts EC, Hulscher MEJL, Mouton JW, et al. Current evidence on hospital antimicrobial stewardship objectives: a systematic review and meta-analysis. *Lancet Infect Dis.* 2016;16(7):847–856. doi:10.1016/S1473-3099(16)00065-7

- Seale AC, Blencowe H, Zaidi A, et al. Neonatal severe bacterial infection impairment estimates in South Asia, sub-Saharan Africa, and Latin America for 2010. *Pediatr Res.* 2013;74(S1):73–85. doi:10.1038/pr.2013.207

- Singer M, Deutschman CS, Seymour CW, et al. The Third International Consensus Definitions for Sepsis and Septic Shock (Sepsis-3). *JAMA.* 2016;315(8):801–810. doi:10.1001/jama.2016.0287

- Slimings C, Riley TV. Antibiotics and healthcare facility-associated *Clostridioides difficile* infection: updated systematic review and meta-analysis. *J Antimicrob Chemother.* 2021;76(7):1676–1688. doi:10.1093/jac/dkab091

- Solomkin JS, Mazuski JE, Bradley JS, et al. Diagnosis and management of complicated intra-abdominal infection in adults and children: guidelines by the Surgical Infection Society and the Infectious Diseases Society of America. *Clin Infect Dis.* 2010;50(2):133–164. doi:10.1086/649554

- Stanaway JD, Parisi A, Sarber K, et al. The global burden of non-typhoidal salmonella invasive disease: a systematic analysis for the Global Burden of Disease Study 2017. *Lancet Infect Dis.* 2019;19(12):1312–1324. doi:10.1016/S1473-3099(19)30418-9

- Stevens DL, Bisno AL, Chambers HF, et al. Practice guidelines for the diagnosis and management of skin and soft tissue infections: 2014 update by the Infectious Diseases Society of America. *Clin Infect Dis.* 2014;59(2):e10–e52. doi:10.1093/cid/ciu296

- Taplitz RA, Kennedy EB, Bow EJ, et al. Antimicrobial prophylaxis for adult patients with cancer-related immunosuppression: ASCO and IDSA clinical practice guideline update. *J Clin Oncol.* 2018;36(30):3043–3054. doi:10.1200/JCO.18.00374

- Tong SYC, Davis JS, Eichenberger E, Holland TL, Fowler VG Jr. *Staphylococcus aureus* infections: epidemiology, pathophysiology, clinical manifestations, and management. *Clin Microbiol Rev.* 2015;28(3):603–661. doi:10.1128/CMR.00134-14

- Trampuz A, Zimmerli W. Prosthetic joint infections: update in diagnosis and treatment. *Swiss Med Wkly.* 2005;135(17-18):243–251. doi:10.4414/smw.2005.10934

- Troeger C, Blacker BF, Khalil IA, et al. Estimates of the global, regional, and national morbidity, mortality, and aetiologies of diarrhoea in 195 countries: a systematic analysis for the Global Burden of Disease Study 2016. *Lancet Infect Dis.* 2018;18(11):1211–1228. doi:10.1016/S1473-3099(18)30362-1

- Tunkel AR, Hartman BJ, Kaplan SL, et al. Practice guidelines for the management of bacterial meningitis. *Clin Infect Dis.* 2004;39(9):1267–1284. doi:10.1086/425368

- UN Tourism. *UN Tourism Data Dashboard: Global and Regional Tourism Performance.* 2025. Accessed March 24, 2026. https://www.untourism.int/tourism-data/global-and-regional-tourism-performance

- Unemo M, Lahra MM, Escher M, et al. WHO global antimicrobial resistance surveillance (GASP/GLASS) for *Neisseria gonorrhoeae* 2017–2018: a retrospective observational study. *Lancet Microbe.* 2021;2(11):e627–e636. doi:10.1016/S2666-5247(21)00171-3

- Unemo M, Shafer WM. Antimicrobial resistance in *Neisseria gonorrhoeae* in the 21st century: past, evolution, and future. *Clin Microbiol Rev.* 2014;27(3):587–613. doi:10.1128/CMR.00010-14

- United Nations Department of Economic and Social Affairs, Population Division. *World Population Prospects 2024.* Accessed March 24, 2026. https://population.un.org/wpp/

- van Boeckel TP, Gandra S, Ashok A, et al. Global antibiotic consumption 2000 to 2010: an analysis of national pharmaceutical sales data. *Lancet Infect Dis.* 2014;14(8):742–750. doi:10.1016/S1473-3099(14)70780-7

- van de Beek D, Brouwer MC, Thwaites GE, Tunkel AR. Advances in treatment of bacterial meningitis. *Lancet.* 2012;380(9854):1693–1702. doi:10.1016/S0140-6736(12)61186-6

- van Schaik W. The human gut resistome. *Philos Trans R Soc Lond B Biol Sci.* 2015;370(1670):20140087. doi:10.1098/rstb.2014.0087

- Verani JR, McGee L, Schrag SJ; Division of Bacterial Diseases, National Center for Immunization and Respiratory Diseases, Centers for Disease Control and Prevention. Prevention of perinatal group B streptococcal disease - revised guidelines from CDC, 2010. *MMWR Recomm Rep.* 2010;59(RR-10):1–36.

- Wain J, Kilmarx PH, eds. *Practical Laboratory Manual for National Tuberculosis Programmes.* Geneva: WHO; 2006.

- Watkins DA, Johnson CO, Colquhoun SM, et al. Global, regional, and national burden of rheumatic heart disease, 1990–2015. *N Engl J Med.* 2017;377(8):713–722. doi:10.1056/NEJMoa1603693

- Werner G, Coque TM, Hammerum AM, et al. Emergence and spread of vancomycin resistance among enterococci in Europe. *Euro Surveill.* 2008;13(47):19046.

- Wertheim HFL, Melles DC, Vos MC, et al. The role of nasal carriage in *Staphylococcus aureus* infections. *Lancet Infect Dis.* 2005;5(12):751–762. doi:10.1016/S1473-3099(05)70295-4

- WHO Collaborating Centre for Drug Statistics Methodology. *DDD: definition and general considerations.* Accessed June 22, 2026. https://atcddd.fhi.no/ddd/definition_and_general_considera/

- Workowski KA, Bachmann LH, Chan PA, et al. Sexually transmitted infections treatment guidelines, 2021. *MMWR Recomm Rep.* 2021;70(4):1–187. doi:10.15585/mmwr.rr7004a1

- World Bank. *Air transport, passengers carried (IS.AIR.PSGR).* World Development Indicators; source: International Civil Aviation Organization (ICAO). Accessed March 24, 2026. https://data.worldbank.org/indicator/IS.AIR.PSGR

- World Bank. *Hospital beds (per 1,000 people) (SH.MED.BEDS.ZS).* World Development Indicators; source: World Health Organization. Accessed March 24, 2026. https://data.worldbank.org/indicator/SH.MED.BEDS.ZS

- World Bank. *International tourism, number of departures (ST.INT.DPRT).* World Development Indicators; source: UN Tourism. Accessed March 24, 2026. https://data.worldbank.org/indicator/ST.INT.DPRT

- World Health Organization. *Global antibiotic resistance surveillance report 2025.* Geneva: WHO; 2025. ISBN 9789240116337. https://www.who.int/publications/i/item/9789240116337

- World Health Organization. *Global Antimicrobial Resistance and Use Surveillance System (GLASS) report: antibiotic use data for 2022.* Geneva: WHO; 2025. ISBN 9789240108127. https://www.who.int/publications/i/item/9789240108127

- World Health Organization. *Global Antimicrobial Resistance and Use Surveillance System (GLASS).* Accessed June 22, 2026. https://www.who.int/initiatives/glass

- World Health Organization. *Global tuberculosis report 2025.* Geneva: WHO; 2025. https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2025

- World Health Organization. *Sexually transmitted infections (STIs).* Fact sheet. Updated September 10, 2025. Accessed June 22, 2026. https://www.who.int/news-room/fact-sheets/detail/sexually-transmitted-infections-(stis)

- World Health Organization. *Universal health coverage (UHC).* Fact sheet, 2025. Accessed March 24, 2026. https://www.who.int/news-room/fact-sheets/detail/universal-health-coverage-(uhc)

- World Health Organization. *WHO consolidated guidelines on drug-resistant tuberculosis treatment.* Geneva: WHO; 2020. ISBN 978-92-4-155056-7. Available at: https://www.who.int/publications/i/item/9789241550567

- Wunderink RG, Matsunaga Y, Ariyasu M, et al. Cefiderocol versus high-dose, extended-infusion meropenem for the treatment of Gram-negative nosocomial pneumonia (APEKS-NP): a randomised, double-blind, phase 3, non-inferiority trial. *Lancet Infect Dis.* 2021;21(2):213–225. doi:10.1016/S1473-3099(20)30731-3

- Xu L, Sun X, Ma X. Systematic review and meta-analysis of mortality of patients infected with carbapenem-resistant *Klebsiella pneumoniae*. *Ann Clin Microbiol Antimicrob.* 2017;16(1):18. doi:10.1186/s12941-017-0191-3

- Yeung KHT, Duclos P, Nelson EAS, Hutubessy RCW. An update of the global burden of pertussis in children younger than 5 years: a modelling study. *Lancet Infect Dis.* 2017;17(9):974–980. doi:10.1016/S1473-3099(17)30390-0
