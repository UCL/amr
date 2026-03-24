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

The model tracks **42 bacterial species**, **61 antibiotics** (grouped into **39 internal drug classes**), and **40 resistance mechanisms**. The population is distributed across **6 world regions** (North America, Europe, Asia, Oceania, South America, Africa), each with distinct epidemiological, travel, hospitalisation, and healthcare profiles.


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

**Calibration.** The model's parameters (the numbers that control how frequently infections occur, how often drugs are prescribed, how quickly resistance emerges, etc.) are adjusted — *calibrated* — so that the model's outputs match real-world data. For example, the model is calibrated against observed antibiotic consumption rates, resistance prevalence reported by surveillance networks (such as ECDC and CDC), and infection incidence data. Some parameters therefore behave as **effective model parameters** rather than direct one-to-one measurements from surveillance datasets: they are chosen to reproduce the joint behaviour of complex clinical systems that are only partially observed. This is especially true for access modifiers, composite vulnerability states, and behaviourally driven prescribing terms. Sections 2–10 describe what the model does; Appendix B lists all the parameter values.

Throughout this document, region- and age-specific parameter tables should generally be interpreted as **qualitatively constrained model structures**: the ordering and rough scale are informed by global demographic, travel, and health-system datasets, but the exact numeric values remain modelling choices that are subsequently checked against calibration targets rather than copied directly from any single source.


### 1.3 Scope and purpose

The model is specifically designed for reconstructing the historical emergence and growth of AMR over time by mechanistically linking antibiotic consumption, biological mutability, and transmission. It evaluates the potential impact of antibiotic stewardship policies by recreating empirical observations of resistance incidence and separating resistance acquisition across different care settings (e.g., community-acquired versus hospital-acquired).


### 1.4 How to read the rest of this document

The document is structured to follow the journey a person takes through the model:

| Section | What it covers |
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
| **11. Policy Evaluation** | Comparing stewardship interventions and counterfactual scenarios |
| **12. Limitations** | What the model does not capture / caveats for interpretation |
| **Appendices** | Reference tables of all bacteria, drugs, parameters, and outputs |



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



These shares are intended as a coarse world-population partition for simulation purposes rather than a literal census reconstruction of any single year. Their ordering and approximate magnitudes are consistent with the United Nations *World Population Prospects 2024*, which provides official demographic estimates and projections across global regions and countries (UN DESA Population Division, 2024).

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



These age bands are structural groupings rather than claims about sharply separated biological states. They are meant to preserve widely observed global gradients in infection burden and mortality risk, especially the concentration of severe infectious outcomes at the extremes of age and the relative protection of school-age and younger adult groups in many syndromes (GBD 2019 Lower Respiratory Infections Collaborators, 2022).

**Sepsis/mortality age categories** (a separate, finer grouping):

Neonates (0–28 days) have dramatically different infection risks and case-fatality rates compared to older infants — for example, Group B *Streptococcus* sepsis in a 5-day-old neonate is a very different clinical entity from a respiratory infection in a 10-month-old. To avoid averaging over these differences, the model uses a separate age classification for sepsis onset and infection-related mortality:

| Category | Age Range |
|----------|-----------|
| Neonatal | 0–28 days |
| Paediatric | 28 days–18 years |
| Young Adult | 18–50 years |
| Elderly | 50+ years |



### 2.3 Immunodeficiency

In real life, patients with weakened immune systems — whether from HIV, chemotherapy, organ transplantation, or simply old age — are at much higher risk of infection, harder to treat, and more likely to die (Fishman JA, 2007; Taplitz RA et al., 2018). The model captures this through two types of immunosuppression.

At simulation start, a configurable fraction of the population is seeded into this broader higher-risk host state (`immunosuppression_startup_seed_fraction`, baseline 5%). This startup seeding is a calibration device to avoid an unrealistically long burn-in before immunocompromised-host effects become visible in the simulated population. Published US NHIS analyses place self-reported immunosuppression among adults in the low-single-digit to mid-single-digit range over the last decade (2.7% in 2013 and 6.6% in 2021), so a 5% startup seed sits within the right order of magnitude for a broadened composite vulnerability construct while still remaining a model initial-condition choice rather than a direct epidemiologic estimate (Martinson ML et al., 2024).

**Temporary immunosuppression** represents acute episodes such as a short course of steroids, a viral illness that transiently suppresses immunity, or post-surgical immunosuppression. People enter this state at a rate of `0.00005` per day and recover at `0.01` per day (average duration ~100 days).

**Chronic immunosuppression** represents long-term conditions like HIV/AIDS, solid organ transplant, or autoimmune disease requiring ongoing immunosuppressive therapy. It develops at `0.00006` per day and recovers much more slowly at `0.0012` per day.

When a new immunodeficiency episode occurs in the model, the following age-band probabilities determine whether that episode is typed as **chronic** rather than **temporary**. They are therefore best read as a structural mapping from age to chronic-vs-temporary assignment, not as literal age-specific prevalence estimates of diagnosed immunodeficiency in the underlying population:

| Age group | Probability of chronic typing | Interpretation |
|-----------|------------|---------|
| 0–1 year | 30% | Allows some early-life episodes to map to persistent congenital or neonatal high-risk states |
| 1–18 years | 20% | Keeps most childhood episodes temporary while permitting a smaller chronic subgroup |
| 18–65 years | 40% | Shifts more episodes into persistent high-risk states compatible with HIV, transplantation, or long-term immunosuppression |
| 65+ years | 60% | Makes late-life episodes more likely to persist as a composite frailty/immunosenescence-type vulnerability state |



These probabilities should be read as part of a **composite infection-vulnerability state**, not as literal prevalence estimates of formal immunodeficiency diagnoses. In other words, the model deliberately aggregates classic immunodeficiency, transplant medicine, chemotherapy-related neutropenia, advanced HIV, frailty, and other clinically important causes of impaired host defence into one tractable state variable (Fishman JA, 2007; Taplitz RA et al., 2018). In the current implementation, the seeded starting population uses the configurable startup fraction described above, and the chronic-versus-temporary split follows the same age-stratified mapping shown here. The same age-stratified probabilities also govern the typing of newly arising immunodeficiency episodes during the simulation.

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

Hospital admission matters for AMR because hospitals are where the most resistant organisms are found, where the broadest-spectrum antibiotics are used, and where vulnerable patients are concentrated (Magill SS et al., 2018). The model captures this by simulating daily admission decisions, length of stay, and the elevated risks of hospital-acquired (nosocomial) infection.

**Who gets admitted?** Each day, the model calculates a probability of hospital admission for every person using a logistic model (see Section 1.2 for an explanation of log-odds). The key factors are:

| Factor | Log-odds contribution | What it means |
|--------|----------------------|---------------|
| Baseline (healthy person) | −10.4 | Very low daily risk (~0.003%) — most people are not admitted on any given day |
| Age | +0.02 per year | Older patients are progressively more likely to be admitted |
| Sepsis | +4.4 | Sepsis is a strong driver of admission (~80× multiplier) |
| Symptomatic infection (severity > 3.0) | +2.5 | Moderate-to-severe infections prompt admission (~12× multiplier) |
| Regional healthcare access | varies (see below) | Reflects real-world differences in hospital capacity |



**Length of stay:** Once admitted, patients are discharged at a rate of `0.28` per day (average stay ~3.6 days), with a hard maximum of 30 days. Patients with active sepsis or those currently receiving intravenous (IV) antibiotics cannot be discharged — they remain in hospital until resolving or completing their IV course. This should be interpreted as an **effective all-cause discharge hazard** in the model, not as a claim that every real-world admission has the same geometric length-of-stay distribution.

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



These modifiers should be read as a qualitative ordering of effective hospital access rather than literal estimates of admission probabilities. The ranking is consistent with broad cross-country differences in service coverage and infrastructure documented by WHO's universal health coverage monitoring framework and the World Bank's hospital-bed indicator, which show persistent between-country variation in effective access to care and inpatient capacity even as global service coverage has improved over time (WHO, 2025; World Bank, `SH.MED.BEDS.ZS`).

Negative values mean patients are *less* likely to be admitted — not because they are less sick, but because hospital beds are less available. This matters for AMR because patients who cannot access hospital care may not receive appropriate antibiotics or diagnostics, whereas international sepsis-care programmes have associated better structured in-hospital and ICU bundle delivery with lower hospital mortality (Evans L et al., 2021; Levy MM et al., 2010).

**Nosocomial (hospital-acquired) risks:**

Being in hospital dramatically changes a patient's infection risk profile. Hospital patients are exposed to multi-drug-resistant organisms on surfaces, devices, and other patients. The model captures this with pathogen-specific hospital acquisition modifiers:

| Pathogen | Hospital modifier | Approximate risk multiplier | Clinical context |
|----------|------------------|---------------------------|-----------------|
| *A. baumannii* | +3.4 | ~30× | Ventilator-associated pneumonia, ICU pathogen |
| *E. faecium* | +3.3 | ~27× | Line infections, post-surgical |
| *P. aeruginosa* | +3.0 | ~20× | Burns, wounds, ventilators |
| *S. aureus* | +2.3 | ~10× | Surgical site, line-related infections |
| *K. pneumoniae* | +2.0 | ~7× | Carbapenem-resistant strains in ICU |
| Community pathogens (*C. trachomatis*, *T. pallidum*, *Campylobacter*) | −0.6 to −1.5 | Lower in hospital | Sexually transmitted or food-borne — acquired in the community, not in hospital |



Hospital patients also face higher baseline mortality (+0.262 log-odds, ~1.3×) and higher sepsis onset risk (+0.5 log-odds, ~1.6×), but they also have a higher probability of *recovering* from sepsis (+0.8 log-odds) because of access to intensive care.


### 2.5 Travel

International travel is a well-documented driver of AMR spread. Travellers who visit regions with high resistance prevalence can acquire resistant bacteria and bring them home — this is how, for example, ESBL-producing *E. coli* from South and South-East Asia has spread to European populations (Arcilla MS et al., 2017).

The model simulates this by giving each person a small daily probability of travelling to another region (`0.00005` per day, roughly one trip every 55 years per person). This is intentionally a low **effective cross-region mixing rate**, because the model only needs enough travel to reproduce long-run AMR importation and reseeding; it is not intended to represent literal passenger-trip counts. Travel frequency varies by region of origin, reflecting real-world patterns:

| Region | Travel multiplier | Rationale |
|--------|------------------|-----------|
| Europe | ×3.5 | High international travel rates |
| North America | ×3.0 | High travel, large business travel |
| Oceania | ×2.5 | Geographic distance drives air travel |
| Asia | ×1.5 | Rapidly growing travel volumes |
| South America | ×0.8 | Moderate travel rates |
| Africa | ×0.3 | Lowest international travel rates |



These multipliers are intended as a **qualitative ranking of cross-region mixing intensity**, not as literal estimates of per-capita trip counts. The ordering is supported by broad regional patterns in UN Tourism's global and regional tourism dashboard and World Bank indicators for international tourism departures and air passenger volumes, which collectively show very high international mobility in Europe and North America, strong air-travel dependence in Oceania, rapid growth but substantial heterogeneity across Asia, intermediate volumes in South America, and lower outbound tourism and aviation intensity across much of Africa (UN Tourism, 2025; World Bank, `ST.INT.DPRT`; World Bank, `IS.AIR.PSGR`).

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
- **Circulating resistance landscape** — the EWMA-smoothed prevalence of each resistance mechanism across currently infected individuals shapes the probability that a newly acquired bacterium already carries one or more resistance mechanisms (see Section 3.4)

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

These community-acquisition templates should be read as structured relative-risk shapes rather than literal incidence-rate estimates for each age-region-organism cell. They are intended to preserve broad, globally observed patterns such as the concentration of enteric disease in children, respiratory vulnerability at the extremes of age, and young-adult concentration of sexually transmitted infections, while leaving exact organism-level burden to calibration of the bacterium-specific baseline and interaction terms against the model's target outputs.


### 3.2 Hospital acquisition

Hospitalised patients are exposed to a different set of pathogens and at different rates than people in the community. Instead of using community acquisition rates, the model uses separate hospital-specific acquisition parameters (`{bacteria}_log_odds_hospital_acquired`) for each species.

This reflects the clinical reality that hospitals concentrate drug-resistant organisms: patients on ventilators are exposed to *Acinetobacter* and *Pseudomonas*, patients with central lines to *Staphylococcus* and *Enterococcus*, and patients on broad-spectrum antibiotics to *C. difficile*. The specific hospital modifiers for each pathogen are listed in Section 2.4.

These hospital-acquisition terms are best interpreted as qualitative rankings of nosocomial exposure pressure rather than direct ward-level attack-rate measurements. That is consistent with the way global AMR surveillance systems aggregate routine clinical microbiology data: they show that healthcare-associated pathogen mixes differ systematically from community mixes, but with large between-country differences in sampling intensity, bed capacity, case mix, and laboratory coverage (WHO GLASS, 2026).


### 3.3 Carrier-derived infection

People can carry bacteria in their gut, skin, or respiratory tract without being ill — this is called **asymptomatic carriage** (see Section 8 for details). Occasionally, these carried bacteria can cause an active infection in the same person. This is called **endogenous infection** and is extremely important for AMR because:

- The carried bacteria may already be resistant (having been selected by previous antibiotic courses)
- The person's resistance profile passes directly from carriage to infection via **mechanism-bit copying** — each resistance mechanism present in the microbiome compartment (`mechanism_microbiome`) is independently considered for transfer to the infection compartment (`mechanism_any`)

This pathway is governed by two parameters:

| Parameter | Value | What it means |
|-----------|-------|---------------|
| `carrier_resistance_inheritance_probability` | 0.50 | 50% chance that the carrier-derived infection pathway fires at all — when it does, individual mechanisms are copied from the microbiome to the infection compartment |
| `infection_from_microbiome_dampening` | 0.70 | Per-mechanism transfer probability: each mechanism in the microbiome has a 70% chance of being copied to the infection site, reflecting that not all colonising lineages successfully transition to the infection site |



### 3.4 Resistance at acquisition

When a new infection is acquired from the community, the model needs to decide: is this bacterium resistant to any drugs, and if so, which ones?

Rather than rolling a separate dice for each drug independently (which would produce unrealistic resistance patterns), the model uses a three-step process that reflects how resistance actually exists in bacterial populations:

1. **Mechanism-level prevalence (EWMA tracking)**: After every simulated day, the model updates an **exponential moving average** (EWMA) of the fraction of infected individuals carrying each of the 40 resistance mechanisms, tracked separately for every combination of region × care setting (community / hospital) × bacteria × mechanism. The EWMA smoothing factor (`mechanism_cache_ewma_decay` = 0.9) means today's prevalence estimate is 90% the previous estimate and 10% the newly observed fraction — giving the cache a memory that damps day-to-day noise while still following genuine trends. Mathematically:

   $$\text{EWMA}_{t+1} = \alpha \cdot \text{EWMA}_t + (1 - \alpha) \cdot \frac{\text{infected with mechanism}_t}{\text{total infected}_t}$$

   where $\alpha = 0.9$. Hospital and community populations are tracked separately, so a newly hospitalised patient draws from the hospital strain pool and a community infection draws from the community pool.

2. **Community dilution**: Clinical samples tend to over-represent resistant strains. To account for the fact that community bacteria are less resistant than those seen in clinics, the model applies a dilution factor (`community_resistance_dilution_factor` = 0.50). A draw from random determines whether the infection originates from the human (circulating) reservoir at all; if not, the bacterium is treated as wild-type.

3. **Correlated mechanism profiles — profile-cache sampling**: Rather than independently sampling each mechanism from its marginal EWMA prevalence (which would miss the real-world phenomenon of multiple resistance genes co-travelling on the same plasmid), the model maintains a **profile reservoir** (`MechanismProfileCache`) of up to 200 complete resistance genotypes sampled — via reservoir sampling — from currently infected individuals. Each genotype is stored as a compact 64-bit bitmask (one bit per mechanism). When a new infection is acquired from the human reservoir, the model samples a **complete genotype profile** from this reservoir and assigns all mechanisms set in that profile to the newly infected individual simultaneously. The result is that newly acquired *E. coli*, for example, arrives with a resistance profile that mirrors an actual circulating strain — e.g., ESBL CTX-M together with fluoroquinolone resistance, as these co-occur on real plasmids.

   If the profile cache is empty (early in the simulation, before enough infections have accumulated), the model falls back to sampling a single mechanism from the marginal EWMA cache.

   The `any_r` resistance level reported for each bacterium–drug combination is then **derived** from the mechanisms present in the individual, using the multiplicative susceptibility formula:

   $$\text{any\_r} = 1 - \prod_{m : \text{mechanism}_m \text{ present}} (1 - e_m)$$

   where $e_m$ is the enhancement multiplier for mechanism $m$ against the drug in question (see Section 7.2).

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

**Sepsis** is a life-threatening organ dysfunction caused by a dysregulated host response to infection. In clinical practice, it is treated as a medical emergency with high mortality (Singer M et al., 2016; Evans L et al., 2021). In the model, sepsis is a distinct state that dramatically increases both the urgency of treatment and the risk of death.

Each day, the model calculates the probability of a person's infection progressing to sepsis using a logistic model that combines:

| Risk factor | Parameter | Value | What it means |
|-------------|-----------|-------|---------------|
| Bacterial load | `log_odds_sepsis_infection_level` | +0.9 per unit | The sicker the patient (higher bacterial burden), the more likely sepsis becomes — the single strongest driver |
| Duration of infection | `log_odds_sepsis_infection_duration` | +0.005 per day | Untreated infections gradually become more dangerous |
| Neonatal age | `sepsis_age_log_odds_neonatal` | +1.10 | Neonates are ~3× more likely to develop sepsis |
| Elderly age | `sepsis_age_log_odds_elderly` | +0.69 | Over-70s are ~2× more likely |



**Not all bacteria or body sites carry equal sepsis risk.** Clinically, a bloodstream infection with *N. meningitidis* is incomparably more dangerous than a *C. trachomatis* genital infection (van de Beek D et al., 2012). The model captures this through:

- **Per-bacterium baseline**: Ranges from very low (*E. coli* UTI: −21.0, making sepsis extremely rare for routine UTIs) to high (*N. meningitidis*: −1.2, reflecting its aggressive clinical course)
- **Per-syndrome modifier**: Bloodstream (+1.5) and CNS (+1.2) infections are far more likely to cause sepsis; genitourinary (−2.0) and skin (−1.0) infections far less so

**Regional factors** also affect sepsis risk, reflecting differences in healthcare access and sanitation:
- Europe: −0.6 (best mitigation)
- North America, Oceania: −0.5
- Asia: −0.1 (reference)
- Africa: +0.1 (least mitigation — patients present later and with fewer resources)


### 4.4 Natural clearance

This part of the model actually contains two related but distinct processes, and the previous wording blurred them together:

- **Microbiome or carriage clearance**: `default_microbiome_clearance_probability_per_day` = 0.01 is the default daily chance of losing asymptomatic carriage from the microbiome reservoir, with bacteria-specific overrides for organisms that are known to persist much longer or clear more quickly.
- **Duration penalty on carriage clearance**: `carriage_duration_log_odds_coefficient` = −0.01 per day, capped by `carriage_duration_max_log_odds_effect` = −2.0, applies to microbiome carriage rather than directly to symptomatic infection. The idea is that long-established colonization becomes harder to dislodge because organisms have had time to occupy a stable niche, form biofilms, and adapt to the host environment (Trampuz A et al., 2005).
- **Drug-assisted microbiome clearance**: `microbiome_clearance_probability_on_drug_treatment` = 0.80 is the probability that effective treatment also clears carriage once a drug-treated infection resolves.

**Infection resolution itself is modeled separately.** Infection level changes each day according to bacterial growth, host-driven suppression, and any active antibiotic effect. An infection resolves when the simulated bacterial level is driven down to a near-zero threshold in the rules engine, or when an immune-clearance event is triggered; this is not controlled by `default_microbiome_clearance_probability_per_day`.

This distinction matters for AMR because there can be a delay between infection acquisition and symptom-driven treatment. During that untreated interval, bacteria continue replicating, and resistant subclones can emerge or expand within the infecting population before antibiotics are started. In other words, the main process before treatment is untreated growth and diversification rather than antibiotic selection.

---



## 5. Diagnostic Testing

Diagnostic testing is the bridge between empiric prescribing (guessing which antibiotic to use) and targeted prescribing (knowing which antibiotic will work). In real clinical practice, a doctor sends a sample (blood, urine, sputum) to the microbiology laboratory, which first identifies which bacterium is causing the infection (culture), and then tests which antibiotics it is susceptible or resistant to (antimicrobial susceptibility testing, or AST). This process takes days — and during that waiting period, the patient is treated empirically based on clinical judgement.

The model simulates this entire workflow: the decision to send a test, the delay in getting results, the possibility of laboratory errors, and the historical availability of testing technology.


### 5.1 Historical introduction

Modern diagnostic microbiology did not exist in 1930. The model introduces testing capabilities at historically appropriate time points:

| Technology | Available from | ~ Calendar year | Clinical context |
|------------|---------------|-----------------|-----------------|
| **Bacterial culture** | Day 5,478 | ~1945 | Basic culture techniques became routine in the mid-20th century |
| **Antimicrobial susceptibility testing (AST)** | Day 9,131 | ~1955 | Standardised AST methods (e.g., disc diffusion) followed about a decade later (Bauer AW et al., 1966) |



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



**Regional differences:** Laboratory capacity varies dramatically around the world. Many hospitals in sub-Saharan Africa lack the microbiological infrastructure that is routine in European hospitals (Jacobs J et al., 2019). The model captures this with regional testing multipliers:

| Region | Testing multiplier | Context |
|--------|-------------------|---------|
| Europe | ×1.2 | Highest testing density |
| North America | ×1.1 | High infrastructure |
| Oceania | ×0.8 | Good but geographically dispersed |
| Asia | ×0.7 | Highly variable by country |
| South America | ×0.6 | Variable access |
| Africa | ×0.3 | Very limited lab infrastructure in many settings |



These regional differences have direct consequences for AMR: in settings where testing is rare, patients are more likely to continue on ineffective empiric therapy, creating selection pressure for resistance without the feedback loop of culture results to guide narrower prescribing.

As with the admission and travel modifiers above, these testing multipliers should be read as qualitative effective-capacity terms rather than literal claims about national culture rates. They combine laboratory availability, specimen transport, clinician ordering behaviour, turnaround reliability, and AST reporting infrastructure, which is the same bundle of constraints emphasized by WHO's GLASS laboratory-strengthening programme and reviews of district-level bacteriology capacity in resource-limited settings (Jacobs J et al., 2019; WHO GLASS, 2026).

---



## 6. Antibiotic Treatment

This section covers the entire antibiotic prescribing process as the model simulates it — from the decision to start an antibiotic, through drug selection and dosing, to stopping the course. This is the heart of the AMR model: antibiotic use drives the selection pressure that causes resistance to emerge and spread.

The model aims to reproduce how antibiotics are actually prescribed in clinical practice — including imperfect decisions, regional variation in drug access, and the difference between empiric therapy (best-guess prescribing before lab results are available) and targeted therapy (prescribing guided by culture and susceptibility results).


### 6.1 Treatment initiation — deciding to start antibiotics

Each day, the model decides whether to start a new antibiotic course for each person, using a logistic model (see Section 1.2). The probability of starting antibiotics depends on the person's clinical state:

| Factor | Log-odds | Approximate effect | Clinical rationale |
|--------|----------|-------------------|-------------------|
| Baseline (no symptoms) | −5.5 | ~0.4% daily chance | Occasionally antibiotics are prescribed without clear indication — this captures "just in case" prescribing and background inappropriate antibiotic use seen in ambulatory care (Fleming-Dutra KE et al., 2011) |
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



These access barriers have a paradoxical effect on AMR: in settings where antibiotics are hard to obtain, people die of treatable infections, but the selection pressure for resistance is also lower. The model captures both sides of this equation. The regional prescribing modifiers should therefore be read as **effective access-and-behaviour terms** combining healthcare access, affordability, dispensing practice, and care-seeking, rather than as pure pharmacy-supply measurements.


### 6.2 Drug selection — choosing which antibiotic to use

Once the model decides to start an antibiotic, it must choose *which* antibiotic. This is one of the most clinically complex parts of the model, because the choice depends on what the doctor knows at the time.

**Two modes of prescribing:**

1. **Empiric therapy** — the doctor has no lab results yet and must choose a drug based on clinical judgement. "The patient has a UTI; guidelines say to use nitrofurantoin or trimethoprim." The model uses syndrome-specific scoring templates (see below) to replicate this guideline-based prescribing. If there is no meaningful syndrome-specific signal, the candidate drug is treated as ineffective and is heavily penalised rather than being allowed to compete on generic properties alone.

2. **Targeted therapy** — lab results have identified the bacterium and its susceptibility profile. The doctor can now choose a drug known to work. The model strongly rewards narrow-spectrum choices at this stage (×5.0 bonus for narrow-spectrum drugs) and penalises unnecessary broad-spectrum use (×0.1 penalty), reflecting the principle of antibiotic de-escalation that sits at the core of antimicrobial stewardship guidance and is supported by hospital stewardship evidence from Europe and Asia (Barlam TF et al., 2016; Schuts EC et al., 2016; Lee CF et al., 2018).

**How drug scoring works:**

For each candidate drug, the model calculates a score based on several factors. The final candidate scores are placed into a weighted index (probabilistic selection) using a temperature-scaled power function: `Weight = Score^(1.0 / Temperature)`. The baseline `drug_selection_temperature` is **0.55**. A lower temperature makes prescribing more deterministic (strongly favouring the highest score), while a higher temperature reflects stochastic variance (idiosyncratic prescribing habits) in clinical settings. 

| Scoring factor | Empiric phase | Targeted phase | What it captures |
|---------------|---------------|----------------|-----------------|
| Syndrome-specific template score | Primary driver | Secondary | How well this drug matches guidelines for the infection site |
| Spectrum width | Slight bonus (×0.85) for broad-spectrum | Strong penalty (×0.1) for broad-spectrum | Empiric: cast a wide net. Targeted: use the narrowest effective drug |
| Known ineffectiveness | Near-zero score (×0.001) | Near-zero score (×0.001) | Never select a drug that is known to not work |
| Narrow-spectrum bonus | — | ×5.0 | Reward de-escalation to targeted therapy |



**Restricted niche agents:** Some drugs are hard-blocked outside their clinically plausible niche. In particular, **retapamulin** and **fusidic acid** are now restricted to **skin/soft-tissue prescribing contexts** and are excluded from undifferentiated prophylaxis, no-syndrome empiric starts, sepsis, and non-skin systemic infections. In targeted therapy they are only allowed when the identified pathogen set is consistent with the narrow skin-focused niche (currently *Staphylococcus aureus* or *Streptococcus pyogenes*). This is intended to reflect their main clinical role as topical or narrowly targeted anti-staphylococcal/anti-impetigo agents rather than general systemic therapy (Stevens DL et al., 2014; Koning S et al., 2012).

The same site-restriction logic is applied to **nitrofurantoin** and **furazolidone**, which are limited to uncomplicated lower-UTI contexts and excluded from sepsis, bloodstream infection, and non-urinary syndromes because they are not reliable systemic treatment options (Gupta K et al., 2011).

**Regional resistance surveillance:** If population-level resistance data shows that a drug class is failing frequently in the region, the model penalises empiric use of that drug — mimicking real-world guideline updates when local resistance rates exceed thresholds:

| Local resistance rate | Empiric score penalty | Clinical parallel |
|----------------------|----------------------|------------------|
| >60% resistant | ×0.3 | Drug dropped from guidelines (e.g., ciprofloxacin for *E. coli* UTI in South-East Asia) |
| >45% resistant | ×0.5 | Drug used cautiously, alternatives preferred |
| >10% resistant | ×0.8 | Drug still used but with awareness of resistance risk |



The syndrome scoring tables below are therefore stylised prescribing-preference weights, not literal market-share estimates for each antibiotic. They are designed to preserve broad world-recognisable clinical tendencies such as narrower outpatient UTI therapy, broader empiric treatment for sepsis and intra-abdominal infection, and de-escalation after microbiology results, while allowing the realised prescribing mix to emerge from access constraints, testing availability, and resistance feedback.


#### Treatment cessation — stopping antibiotics

Patients stop their antibiotic course based on several factors.

These values are best interpreted as **daily probabilities of prematurely stopping treatment**, not as the inverse of total course length. In other words, they are dropout hazards calibrated so that most patients still remain on therapy through a guideline-like treatment window.

| Scenario | Daily stop probability | Approximate implication | Real-world parallel |
|----------|----------------------|-------------------------|-------------------|
| Default course | 0.45% per day | About 94% of patients are still on treatment by day 14 | Standard course for many infections |
| No active infection found | 15% per day | Rapid discontinuation over the next few days once infection seems absent | Antibiotics stopped when investigation shows no infection |
| Cholera / *E. coli* GI | 2.5% per day | Supports short-course therapy, with most patients still on treatment through about 3-5 days | Short courses per guidelines |
| *S. aureus* / *S. pneumoniae* | 1.5% per day | About 90% of patients are still on treatment by day 7 | Standard courses |
| MDR-TB | 0.06% per day | About 90% of patients are still on treatment by 6 months before regional adherence modifiers | Prolonged anti-TB regimens |

So the previous wording that mapped these directly to "typical course length" was too literal. A constant daily stop probability of 0.45% does not mean an average course of 14 days; it means treatment is only rarely interrupted on any given day, so most courses survive to around two weeks.



#### Syndrome-specific empiric scoring templates

The tables below show which drugs score highest for empiric prescribing in each syndrome. Higher scores mean the drug is more likely to be selected. These templates are calibrated to match real-world prescribing guidelines — for example, nitrofurantoin and trimethoprim-sulfamethoxazole score highest for UTI, while piperacillin-tazobactam and meropenem score highest for bloodstream infections.

**Syndrome 1 — UTI** *(most common bacterial infection; oral `pen`, `bli`, and `c1_2g` agents plus `sulf` are preferred, with `fq`, `nitrofurans`, and `phosphonic_acids` alternatives)*

| Drug | Score |
|------|-------|
| trim_sulf | 14.0 |
| amoxicillin_clavulanate | 14.0 |
| amoxicillin | 12.0 |
| ciprofloxacin | 12.0 |
| ampicillin | 10.0 |
| levofloxacin | 10.0 |
| nitrofurantoin | 8.0 |
| cephalexin | 8.0 |
| ceftriaxone | 8.0 |
| cefazolin | 7.0 |
| cefuroxime | 7.0 |
| piperacillin_tazobactam | 5.0 |
| cefepime | 4.0 |
| ceftazidime | 4.0 |
| meropenem | 4.0 |
| imipenem_c | 4.0 |
| ertapenem | 4.0 |
| meropenem_vaborbactam | 3.0 |
| ceftazidime_avibactam | 3.0 |



**Syndrome 2 — Skin/Soft Tissue** *(anti-staphylococcal and streptococcal coverage; `pen`, `c1_2g`, `glyc`, `lipoglycopeptides`, and `oxa` dominate, while `pleuromutilins` remain niche topical agents)* (Stevens DL et al., 2014)

| Drug | Score |
|------|-------|
| penicillin_g | 16.0 |
| amoxicillin | 14.0 |
| amoxicillin_clavulanate | 14.0 |
| ampicillin | 13.0 |
| cephalexin | 13.0 |
| cefazolin | 12.0 |
| clindamycin | 12.0 |
| vancomycin | 11.0 |
| linezolid | 10.0 |
| trim_sulf | 9.0 |
| doxycycline | 9.0 |
| minocycline | 9.0 |
| tedizolid | 9.0 |
| dalbavancin | 9.0 |
| quinu_dalfo | 8.0 |
| ciprofloxacin | 4.0 |
| piperacillin_tazobactam | 3.0 |



**Syndrome 3 — Respiratory** *(community-acquired pneumonia pattern: a `pen`/`bli` backbone plus `mls` atypical cover; amoxicillin-clavulanate and penicillins score highest)* (Metlay JP et al., 2019)

| Drug | Score |
|------|-------|
| amoxicillin_clavulanate | 20.0 |
| amoxicillin | 17.0 |
| penicillin_g | 16.0 |
| ampicillin | 15.0 |
| azithromycin | 12.0 |
| clarithromycin | 11.0 |
| ceftriaxone | 9.5 |
| erythromycin | 9.0 |
| cefuroxime | 8.5 |
| piperacillin_tazobactam | 8.0 |
| levofloxacin | 8.0 |
| moxifloxacin | 8.0 |
| cefepime | 7.5 |
| cephalexin | 7.0 |
| linezolid | 7.0 |
| doxycycline | 6.5 |
| vancomycin | 6.5 |
| meropenem | 6.0 |
| imipenem_c | 6.0 |
| ofloxacin | 6.0 |
| minocycline | 5.5 |



**Syndrome 4 — Bloodstream** *(medical emergency; `bli`, `bli_anti_pseudomonal`, `bli_sulbactam`, `carb_group2`, and other broad IV agents with strong bactericidal activity dominate)* (Rhodes A et al., 2016)

| Drug | Score |
|------|-------|
| piperacillin_tazobactam | 18.0 |
| ampicillin_sulbactam | 16.0 |
| amoxicillin_clavulanate | 16.0 |
| meropenem | 13.0 |
| imipenem_c | 13.0 |
| meropenem_vaborbactam | 13.0 |
| ceftazidime_avibactam | 12.5 |
| cefepime | 12.0 |
| ceftazidime | 11.0 |
| vancomycin | 11.0 |
| ceftriaxone | 10.0 |
| ampicillin | 10.0 |
| linezolid | 10.0 |
| amoxicillin | 9.5 |
| tedizolid | 9.0 |
| quinu_dalfo | 8.5 |
| dalbavancin | 8.0 |
| penicillin_g | 6.5 |
| cefazolin | 6.0 |
| ciprofloxacin | 6.0 |
| levofloxacin | 5.5 |
| cephalexin | 4.0 |
| gentamicin | 1.0 |
| tobramycin | 1.0 |
| amikacin | 1.0 |



**Syndrome 5 — Intra-abdominal** *(must cover Gram-negatives and anaerobes; `bli`, `bli_anti_pseudomonal`, `bli_sulbactam`, `carb_group1`, and `carb_group2` are preferred)* (Solomkin JS et al., 2010)

| Drug | Score |
|------|-------|
| piperacillin_tazobactam | 13.0 |
| meropenem | 13.0 |
| ampicillin_sulbactam | 12.5 |
| imipenem_c | 12.5 |
| amoxicillin_clavulanate | 11.5 |
| ertapenem | 11.0 |
| ceftazidime_avibactam | 10.0 |
| meropenem_vaborbactam | 10.0 |
| ceftazidime | 9.0 |
| cefepime | 9.0 |
| ceftriaxone | 9.0 |
| ampicillin | 8.0 |
| ciprofloxacin | 7.0 |
| amoxicillin | 7.0 |
| levofloxacin | 6.5 |
| trim_sulf | 4.0 |
| metronidazole | 2.5 |



**Syndrome 6 — CNS** *(meningitis; only drugs that cross the blood-brain barrier are useful — see Section 6.4)* (Tunkel AR et al., 2004)

| Drug | Score |
|------|-------|
| ceftriaxone | 15.0 |
| ampicillin | 13.0 |
| vancomycin | 13.0 |
| ceftazidime | 12.0 |
| cefepime | 12.0 |
| penicillin_g | 11.0 |
| meropenem | 11.0 |
| linezolid | 10.0 |
| imipenem_c | 10.0 |
| piperacillin_tazobactam | 6.0 |
| chloramphenicol | 2.0 |



**Syndrome 7 — Gastrointestinal** *(`fq` and `mls` dominate uncomplicated bacterial gastroenteritis scoring, with oral `pen`, `bli`, and selected `nitrofurans` also present)*

| Drug | Score |
|------|-------|
| ciprofloxacin | 12.0 |
| azithromycin | 12.0 |
| amoxicillin_clavulanate | 11.0 |
| amoxicillin | 10.0 |
| ampicillin | 10.0 |
| levofloxacin | 10.0 |
| ampicillin_sulbactam | 9.0 |
| trim_sulf | 8.5 |
| doxycycline | 8.5 |
| minocycline | 6.5 |
| penicillin_g | 5.0 |
| cephalexin | 5.0 |
| cefuroxime | 5.0 |
| furazolidone | 3.0 |
| metronidazole | 1.0 |



**Syndrome 8 — Genital/Pelvic** *(STI guidelines: ceftriaxone + azithromycin for gonorrhoea; doxycycline for chlamydia; penicillin G for syphilis)* (Workowski KA et al., 2021)

| Drug | Score |
|------|-------|
| penicillin_g | 14.0 |
| azithromycin | 13.0 |
| ceftriaxone | 13.0 |
| doxycycline | 12.0 |
| amoxicillin_clavulanate | 12.0 |
| amoxicillin | 11.0 |
| cefuroxime | 10.0 |
| clindamycin | 9.0 |
| ampicillin | 9.0 |
| ampicillin_sulbactam | 8.0 |
| ciprofloxacin | 7.0 |
| levofloxacin | 6.5 |
| cephalexin | 6.0 |
| trim_sulf | 5.0 |
| metronidazole | 1.0 |



**Syndrome 9 — Bone/Joint** *(prolonged courses required; good bone penetration essential — `pen` including `flucloxacillin`, `c1_2g`/`c3g`, `fq`, `rifamycins`, and `oxa` feature prominently)*

| Drug | Score |
|------|-------|
| penicillin_g | 14.0 |
| cefazolin | 13.0 |
| ampicillin | 12.0 |
| vancomycin | 12.0 |
| cephalexin | 11.0 |
| ceftriaxone | 11.0 |
| linezolid | 11.0 |
| tedizolid | 10.0 |
| dalbavancin | 10.0 |
| clindamycin | 10.0 |
| ciprofloxacin | 9.0 |
| levofloxacin | 9.0 |
| trim_sulf | 8.0 |
| meropenem | 7.0 |
| piperacillin_tazobactam | 6.5 |
| rifampicin | 2.0 |

This table is illustrative rather than exhaustive. `flucloxacillin` is part of the model's `pen` class and participates in the same prescribing logic even though it is not shown as a separate row in this abbreviated top-score summary.



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

Broad-spectrum antibiotics kill not only the target pathogen but also many commensal ("friendly") bacteria in the gut, skin, and respiratory tract. This collateral damage creates ecological niches for resistant organisms to fill.

The current implementation represents this in two related but distinct ways:

1. `spectrum_breadth` is a stewardship-facing drug property used when scoring treatment choices. In empiric therapy it favors broader agents when coverage is uncertain, while in targeted therapy it rewards de-escalation toward narrower agents once the pathogen is identified.
2. The longer ecological consequence is handled through each drug's `microbiome_disruption_log_odds`, which accumulates into a persistent `microbiome_disruption_level` reservoir. That reservoir decays over time rather than disappearing immediately when treatment stops, and it directly raises the log-odds of later microbiome acquisition events.

So the model consequence is not just that "broad drugs are broad". Broader therapy influences prescribing behavior up front, and microbiome disruption leaves an ecological hangover that can increase later carriage risk even after the course has finished.

Illustrative `spectrum_breadth` values:

| Drug | Breadth | Meaning |
|------|---------|---------|
| penicillin_g | 2.0 (Narrow) | Minimal disruption to the microbiome |
| linezolid | 2.0 (Narrow) | Targets Gram-positives only |
| vancomycin | 2.5 (Narrow-medium) | Mainly Gram-positive spectrum |
| trim_sulf | 3.5 (Medium-broad) | Moderate disruption |
| azithromycin | 4.0 (Broad) | Significant microbiome disruption |
| ceftriaxone | 4.0 (Broad) | Major disruption; linked to *C. difficile* risk (Slimings C et al., 2021) |
| ciprofloxacin | 4.5 (Very broad) | Extensive gut microbiome disruption |
| meropenem | 5.0 (Very broad) | Maximum disruption — the "sledgehammer" antibiotic |

Operationally, this means broad-spectrum therapy can affect the simulation in two downstream places: first by making a drug more attractive for empirical cover but less attractive for narrow targeted de-escalation, and second by increasing later colonization pressure through the microbiome-disruption reservoir that feeds carriage acquisition.



### 6.4 Drug penetration by syndrome

A drug can only work if it reaches the infection site at adequate concentrations. This is particularly important for certain anatomical sites:

- **CNS (meningitis):** The blood-brain barrier blocks most antibiotics. Only a few drugs (ceftriaxone, metronidazole, chloramphenicol, linezolid) achieve therapeutic levels in cerebrospinal fluid.
- **Bone/joint:** Drugs must penetrate dense, poorly vascularised tissue. `rifamycins` and `fq` agents penetrate well; `ag_group1` and `ag_group2` do not.
- **Bloodstream:** By definition, any IV drug achieves full levels here (penetration = 1.0 for all drugs).

Penetration values range from 0.0 (no drug reaches the site) to 1.0 (full systemic concentration available):

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

Not all antibiotics work against all bacteria. Penicillin G is highly effective against *Streptococcus pneumoniae* (potency 0.90) but has zero activity against *Pseudomonas aeruginosa* (intrinsically resistant). The model encodes this in a **potency matrix** — a 42×61 table (42 bacteria × 61 named drugs) where each cell represents the intrinsic activity of that drug against that bacterium when no acquired resistance is present. Resistance mechanisms are then applied on top of that baseline through the separate 39-class enhancement system described in Section 7.2.

Values range from 0.0 (no activity — the drug simply does not work against this organism) to 1.0 (maximum activity). These potency values are based on published MIC (minimum inhibitory concentration) data and clinical breakpoints. If an organism is intrinsically resistant to a drug (defined as having a baseline potency $\le 0.1$), the model strictly prevents any *acquired* resistance mechanisms from being erroneously assigned to or tracked for that organism-drug pair (e.g., *Mycoplasma*, which lacks a cell wall, cannot acquire PBP mutations against penicillins).

Key examples:
- Meropenem vs *E. coli*: 0.95 (very high potency — a carbapenem is one of the most effective drugs against Gram-negatives)
- Vancomycin vs *E. coli*: 0.0 (vancomycin does not work against Gram-negative bacteria)
- Ceftriaxone vs *S. pneumoniae*: 0.90 (standard treatment for pneumococcal meningitis)



### 6.6 Drug availability by region and era

The model simulates antibiotics becoming available at their historical introduction dates. Before sulfanilamide was introduced in 1937, there were no antibiotics in the model. Before penicillin G was introduced in 1942, the `pen` era had not yet begun. Before ciprofloxacin was introduced in 1987, the model had no `fq` agents. This historical layering is essential for reproducing the sequential emergence of resistance over the 20th century.

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

These availability tiers should be interpreted as qualitative access strata rather than audited procurement shares. They summarize broad world patterns in which older essential antibiotics are much more widely available than newer reserve agents, and in which stewardship, financing, regulatory approval, supply-chain reliability, and laboratory support jointly determine whether a drug is realistically usable in practice (WHO GLASS, 2026; WHO, 2025).



#### Drug introduction dates

The 61 antibiotics in the model span 88 years of pharmaceutical development:

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

Additional drugs in the current canonical list that were omitted from the earlier summary are `flucloxacillin` (~1970), `cefixime` (~1989), and `aztreonam_avibactam` (~2025).



**Special case — Colistin:** Colistin was introduced in 1952 but withdrawn from routine use between ~1970 and ~1995 due to severe nephrotoxicity. It was then reintroduced as a last-resort agent for multi-drug-resistant Gram-negative infections (Li J et al., 2006). The model reflects this by dropping colistin availability to 5% during the withdrawal window.



### 6.7 Drug toxicity

Antibiotics are not without harm. Some drugs — particularly `ag_group1`/`ag_group2` agents (nephrotoxicity, ototoxicity) and `poly` (`colistin`, nephrotoxicity) — carry significant toxicity risks. The model simulates drug toxicity as a **reservoir** that accumulates with continued use and decays when the drug is stopped.

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

Patients who are already receiving an effective antibiotic are partially protected against acquiring new infections — the drug in their system kills susceptible bacteria before they can establish. This mirrors the real-world concept of antibiotic prophylaxis (e.g., surgical prophylaxis with cefazolin prevents wound infections) (Bratzler DW et al., 2013).

The model applies a 70% reduction in new infection risk for susceptible organisms when the patient is already on an active antibiotic (`antibiotic_infection_prevention_efficacy` = 0.7). This does *not* protect against resistant organisms — a crucial point, because it means patients on antibiotics are selectively more likely to acquire resistant infections relative to susceptible ones, creating further selection pressure for resistance.
---


## 7. Resistance Dynamics

This section describes the heart of the model — how bacteria become resistant to antibiotics. For a clinician, this section explains the mechanisms behind the resistance patterns you see in microbiology reports. For example, when your lab reports "ESBL-producing *E. coli*", the model tracks the specific enzyme (CTX-M, TEM, or SHV) that produces that phenotype, which drugs it affects, and how it spreads.

The model tracks resistance at the level of individual **mechanisms** — the specific biological tools bacteria use to evade antibiotics. This matters because the same phenotype (e.g., "carbapenem-resistant *K. pneumoniae*") can arise from very different mechanisms (KPC, NDM, OXA-48), each with different implications for treatment, spread, and even which novel drugs might still work.

**Mechanism-centric architecture.** All resistance state is stored as a set of boolean flags — one per mechanism — for each individual's active infection (`mechanism_any`), majority strain (`mechanism_majority`), and microbiome carriage (`mechanism_microbiome`). The scalar resistance metrics (`any_r`, `activity_r`) reported in outputs are **derived** from these mechanism flags via the multiplicative susceptibility formula (Section 7.2) rather than being tracked independently. A single unified `MechanismCache` maintains the population-level picture: it holds both the EWMA-smoothed per-mechanism prevalence (used as fallback sampling) and a reservoir of up to 200 complete clinical resistance genotypes (used for profile-based acquisition — see Section 3.4).


### 7.1 Resistance mechanisms

The model explicitly tracks **40** distinct resistance mechanisms. Each mechanism represents a specific biological pathway: an enzyme that destroys the drug, a mutation that changes the drug's target, a pump that ejects the drug from the cell, or a barrier that prevents the drug entering.

The table below lists every mechanism, the drugs it affects, and which bacterial groups can acquire it. You do not need to memorise this table — it is a reference. The key insight is that each mechanism has a defined scope: ESBL enzymes (rows 1–3) hit `pen`, `c1_2g`, `c3g`, `c4g`, and related monobactam-active entries but not `carb_group1`/`carb_group2`, while KPC and NDM/VIM (rows 6–7) compromise the carbapenem classes as well.


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
  | PBP mosaic | `mutation_pbp_mosaic` | Penicillin-binding protein mosaic mutations (PBP2x/2b/1a in pneumococcus, penA in gonococci, PBP3 in *H. influenzae*) — reduced β-lactam affinity | `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `flucloxacillin`, `amoxicillin_clavulanate`, `ampicillin_sulbactam`, `piperacillin_tazobactam`, `ticarcillin_clavulanate`, `cephalexin`, `cefazolin`, `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefixime`, `cefepime`, `ceftaroline`, `ceftolozane_tazobactam`, `ceftazidime_avibactam`, `aztreonam` | All |
  | mtrCDE efflux | `efflux_mtr_cde` | mtrCDE-type broad efflux pump (Neisseria, Haemophilus, Campylobacter CmeABC) | `erythromycin`, `azithromycin`, `clarithromycin`, `penicillin_g`, `ampicillin`, `amoxicillin`, `piperacillin`, `ticarcillin`, `tetracycline`, `doxycycline`, `minocycline`, `chloramphenicol` | Fastidious, Enteric Pathogens |
   | Unknown | `as_yet_unknown` | Placeholder mechanism (dormant) | Evaluates `true` dynamically for all applied overrides | All (Calibration Placeholder) |



### 7.2 Mechanism–drug-class enhancement multipliers

When a bacterium possesses a resistance mechanism, it does not simply become immune to every drug. Instead, each mechanism **reduces** drug efficacy by a specific amount. The "enhancement multiplier" (0.0–1.0) represents **how much** of a drug's effectiveness is knocked out:

- **0.0** = the mechanism has no effect on this drug (e.g., a tetracycline efflux pump does nothing against meropenem)
- **0.95** = the mechanism eliminates 95% of the drug's activity (e.g., NDM metallo-β-lactamase virtually destroys carbapenem efficacy)
- **1.0** = complete resistance (the drug is useless)

In clinical terms, an enhancement multiplier of 0.85 for ESBL CTX-M against cephalosporins means: if a patient has an ESBL-producing *E. coli* UTI turned treated with ceftriaxone, the drug retains only 15% of its normal killing power — enough to provide some marginal activity but not enough to reliably cure the infection.

There are 40 mechanisms × 39 drug classes = 1,560 individual values. The table below shows the **global default** multiplier for each mechanism (used when a specific per-class value has not been configured):

These enhancement multipliers should be interpreted as qualitative within-model effect sizes rather than literal MIC shifts or breakpoint translations. Their role is to preserve the clinically familiar ordering in which carbapenemases, van genes, and key target-site alterations have very large effects, whereas efflux and permeability mechanisms are usually weaker on their own, while final realised resistance still depends on baseline potency, site penetration, and combination with other mechanisms.

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

The probability of resistance emerging is not simply proportional to drug concentration. Instead, it follows a bell-shaped curve that peaks when drug levels are at **roughly half the therapeutic concentration** (Drlica K et al., 2007):

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

The incidence bands and per-mechanism emergence rates are therefore effective simulation terms, not literal mutation-rate measurements. They absorb unmodeled ecological constraints such as transmission bottlenecks, within-host competition, and sampling sparsity so that common organisms do not accumulate implausible resistance too quickly and rare organisms still generate enough events to matter over the model horizon.




### 7.4 Resistance reversion and fitness costs

Resistance is not free. Maintaining resistance mechanisms costs the bacterium energy and resources — like carrying a heavy suitcase through an airport (Andersson DI et al., 2010). In the absence of antibiotics, resistant bacteria grow more slowly than their susceptible competitors and are gradually outcompeted. This is why resistance can decline after antibiotic use is reduced — a key insight for stewardship policy.

The model assigns each mechanism a daily **reversion rate** — the probability of losing resistance per day when no antibiotic pressure is present. Higher rates mean the mechanism is "expensive" and lost quickly; lower rates mean it is nearly cost-free and persists indefinitely. All per-mechanism reversion rates are scaled by a global calibration multiplier (`mechanism_reversion_rate_global_multiplier`, default 1.0) so that the overall speed of resistance decay can be tuned without changing individual mechanism rates.

Reversion operates in **both** compartments: the active infection (`mechanism_any`, `mechanism_majority`) and the microbiome carriage (`mechanism_microbiome`). In each compartment, a mechanism can only revert on a given day if no antibiotic with selective pressure for that mechanism is currently present — i.e., the mechanism only decays when the fitness cost is uncompensated. When a mechanism reverts, the scalar resistance metrics (`any_r` or `microbiome_r`) are re-derived from the updated mechanism flags.

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



*Note: The system reserves one remaining placeholder variable (`as_yet_unknown`, baseline rate `0.001`) designated for future empirical calibration. `mutation_pbp_mosaic` has been activated as **PBP mosaic mutations** (chromosomal target modification affecting penicillins, cephalosporins, and aztreonam — NOT carbapenems), and `efflux_mtr_cde` as **mtrCDE-type broad efflux** (chromosomal efflux affecting macrolides, penicillins, tetracyclines, and chloramphenicol). Neither is HGT-transferable.*

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

These floors are structural guardrails, not claims about immutable global prevalence minima. They are used only where the model would otherwise erase well-established intrinsic or near-intrinsic non-susceptibility because of finite population noise.



### 7.6 Cross-resistance groups

When a bacterium becomes resistant to one antibiotic, it often becomes resistant to related antibiotics at the same time. For instance, if *E. coli* acquires an ESBL enzyme and becomes resistant to ceftriaxone, you would expect it to also be resistant to other cephalosporins (cefazolin, cefuroxime) and penicillins — because the enzyme destroys the same β-lactam ring in all of them.

The model captures this by defining **cross-resistance groups**: sets of drugs for which resistance is always acquired and lost together. When an individual's *E. coli* becomes resistant to any drug in Group 1 (β-lactams), it simultaneously becomes resistant to all drugs in that group.

These groups are deliberately stylised phenotype bundles rather than exhaustive mechanistic truth tables. They are meant to preserve the broad empirical regularity that related agents often rise and fall together once a dominant mechanism is present, while the mechanism-level layer above still carries the main biological detail.

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

Most bacteria in and on the human body are harmless commensals — they live on skin, in the gut, and in the respiratory tract without causing disease. This is **carriage** (or colonisation), and it is the normal state. Only a small fraction of carried bacteria ever cause active infection, but for several clinically important resistant organisms colonisation is the stage from which later infection emerges (Werner G et al., 2008).

However, carriage is critically important for AMR because the microbiome is where resistance is **stored and exchanged** (van Schaik W, 2015; McInnes RS et al., 2020). A patient who was treated with ciprofloxacin last month may still carry ciprofloxacin-resistant *E. coli* in their gut microbiome. If they develop a UTI from that same resistant strain, the empiric therapy may fail.


### 8.1 Carriage compartments

Each bacterium in the model has a designated ecological niche — where it naturally lives in (or on) the body:

| Compartment | Example bacteria | Clinical relevance |
|-------------|-----------------|-------------------|
| Gut | *E. coli*, *K. pneumoniae*, *Enterococcus spp.*, *Shigella*, *Salmonella*, *C. difficile* | Largest reservoir; disrupted by broad-spectrum antibiotics |
| Respiratory | *S. pneumoniae*, *H. influenzae*, *P. aeruginosa*, *A. baumannii*, *M. catarrhalis*, *M. tuberculosis* | Carriage often precedes pneumonia |
| Skin/Soft tissue | *S. aureus*, *S. epidermidis* | Nasal/skin MRSA carriage drives surgical wound infections |
| Genitourinary | *N. gonorrhoeae*, *C. trachomatis*, *M. genitalium*, *T. pallidum*, *S. agalactiae* | Asymptomatic STI carriage enables transmission |



These compartment assignments are simplified ecological defaults rather than a full atlas of colonisation niches. They mainly provide the model with the right qualitative reservoirs for bystander selection, endogenous infection, and HGT opportunity.


### 8.2 Resistance in the microbiome

The microbiome serves as a hidden reservoir of resistance. Each individual carries a per-mechanism boolean resistance array (`mechanism_microbiome`) for every organism, mirroring the structure of the infection compartment (`mechanism_any`). The scalar `microbiome_r` metric is **derived** from these mechanism flags via the same multiplicative susceptibility formula used for infection resistance (Section 7.2). This unified mechanism-centric architecture ensures that resistance in carriage and infection compartments is always coherent — there is no separate "float-based" tracking for the microbiome.

Key dynamics:

| Process | Parameter | Value | What it means |
|---------|-----------|-------|---------------|
| Resistance seeding on acquisition | `microbiome_resistance_multiplier_on_acquisition` | 0.50 | When a person acquires a new carriage episode, there is a 50% probability that the colonising strain inherits the circulating resistance profile (sampled from the mechanism profile cache, just as for infection acquisition — see Section 3.4). If the draw fails, the strain arrives susceptible. |
| Established colonies harder to clear | `carriage_duration_log_odds_coefficient` | −0.01/day (caps at −2.0) | The longer a resistant strain has been carried, the harder it is to eradicate — mature colonies are ~7× harder to clear than newly acquired ones |
| Mechanism-level reversion | `mechanism_reversion_rate_global_multiplier` | 1.0 | Per-mechanism reversion operates in the microbiome compartment using the same rates and global multiplier as in the infection compartment (Section 7.4). Each mechanism can only revert when no selecting antibiotic is present. |
| De-novo emergence under treatment | `microbiome_de_novo_multiplier` | 1.0 | When antibiotics exert selective pressure on carried bacteria, resistance mechanisms can emerge in the microbiome via the same emergence formula used for infections (Section 7.3), scaled by this multiplier. Emergence writes to both `mechanism_microbiome` and `mechanism_any`. |
| Carrier → infection bridge | `carrier_resistance_inheritance_probability` | 0.50 | When a carrier develops an endogenous infection, each mechanism in `mechanism_microbiome` is independently considered for transfer to `mechanism_any` (see Section 3.3) |
| Infection → microbiome transfer | (automatic) | — | When an infected individual also carries the same bacterium, resistance mechanisms present in the infection but absent in the microbiome are copied to `mechanism_microbiome`, reflecting spillover from the active infection back into the commensal reservoir |
| HGT into the microbiome | (see Section 9) | — | When a horizontal gene transfer event fires and the recipient carries the donor's bacterium in the microbiome, the transferred mechanism is written to `mechanism_microbiome` as well as `mechanism_any` |



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



These compatibility probabilities should be read as qualitative transfer-likelihood weights rather than experimentally measured conjugation frequencies. They encode the main biological asymmetries the policy model needs: exchange is easiest within closely related Gram-negative groups, harder across distant taxa, and effectively absent across some structural boundaries.


### 9.2 The HGT process

Each day, for every individual carrying resistant bacteria in their microbiome, the model evaluates potential gene transfer events. The model evaluates HGT dynamically per distinct resistance mechanism, allowing independent plasmids (e.g., KPC and *mcr-1*) to transmit independently rather than as a single all-or-nothing block. Furthermore, bacteria do not restrict plasmid donation to only the dominant strain; minority resistance populations can donate, but face a transfer penalty.

When an HGT event fires, the transferred mechanism is written to the recipient's `mechanism_any` (infection compartment). If the recipient also carries the donor's target bacterium in the microbiome, the mechanism is simultaneously written to `mechanism_microbiome`, ensuring the carriage reservoir stays consistent with the infection compartment. All HGT rates are scaled by a global calibration multiplier (`hgt_multiplier`, default 1.0).

| Step | Parameter | Value | Clinical parallel |
|------|-----------|-------|-------------------|
| Global HGT scaling | hgt_multiplier | 1.0 | Calibration knob — scales all HGT rates up or down uniformly |
| Base transfer rate | microbiome_resistance_transfer_probability_per_day | 0.0001 | Background rate — equivalent to a conjugation event occurring every ~27 years per carrier, reflecting how rare HGT is without antibiotic pressure |
| Amplification during antibiotic therapy | hgt_antibiotic_pressure_multiplier | 1.50 (×1.5) | Antibiotic stress triggers the bacterial SOS response, which activates mobile genetic elements and increases conjugation rates by 50% (Beaber JW et al., 2004) — one of the reasons antibiotic use drives resistance even beyond the target pathogen |
| Hospitalization boost | hgt_hospital_multiplier | 3.0 (×3.0) | Captures increased transmission risks in clinical environments where close physical proximity and shared infrastructure elevate exchange. |
| Co-infection baseline | hgt_coinfection_multiplier | 1.25 (×1.25) | Active multi-pathogen infections slightly increase the probability of genetic collision. |
| Microbiome-only penalty | hgt_microbiome_only_penalty | 0.65 (×0.65) | Asymptomatic carriage interactions are less frequent than active infection environments. |
| Gut compartment boost | hgt_gut_compartment_multiplier | 2.0 (×2.0) | The gut has higher bacterial density and provides more conjugation opportunities compared to skin or respiratory tracts. |
| Minority donor penalty | hgt_minority_donor_multiplier | 0.20 (×0.20) | If a donor bacterium carries resistance as a minority strain (sub-dominant), its probability of successful conjugation is penalized by 80%. |



The absolute HGT probabilities are intentionally low and should be interpreted as effective daily hazards at the model scale, not bedside-measurable event rates. Their main purpose is to preserve plausible relative ordering between low-contact community settings, antibiotic-stressed microbiomes, and high-contact hospital environments.

## 10. Mortality

The model tracks mortality from three sources: background (non-infection) causes, **sepsis** (organ dysfunction from uncontrolled infection), and **non-sepsis infection death** (direct tissue damage, toxin production, or chronic complications of infection that do not involve the sepsis cascade). This dual-pathway architecture reflects the clinical reality that different pathogens kill through fundamentally different mechanisms (Rudd KE et al., 2017).


### 10.1 Background mortality

Everyone faces a baseline mortality risk that increases with age:

| Factor | Parameter | Value | What it means |
|--------|-----------|-------|---------------|
| Aging penalty | `log_odds_mortality_per_year_of_age` | 0.04 | Each year of age adds ~4% relative increase in daily death risk (exp(0.04) = 1.04) |
| Elderly frailty acceleration | `log_odds_mortality_per_year_of_age_squared` | 0.05 | A quadratic term that makes mortality rise faster above ~70 — capturing how an 85-year-old is much frailer than a 65-year-old |



These parameters operate on a log-odds scale, so they compound multiplicatively over time.

They should be read as effective demographic mortality-shape terms rather than direct life-table fits for any single country or year. Their role is to preserve the globally familiar pattern of sharply rising all-cause mortality with age and frailty while allowing the simulation's infection-specific pathways to add the AMR-relevant excess risk on top.


### 10.2 Sepsis mortality

Sepsis is the primary death pathway for classic invasive bacterial pathogens. When a person's infection progresses to sepsis (see Section 4.3), the model applies an aggressively escalated daily death risk using a logistic model. The probability of dying from sepsis each day depends on the patient's age, immune status, bacterial burden, and access to hospital care. Without effective antibiotics, sepsis is rapidly fatal — and resistant organisms that are untreatable with empiric therapy are exactly the scenario where this matters most (Murray CJL et al., 2019).

**Per-bacterium sepsis baseline log-odds** control how likely each species is to cause sepsis. These range from very low for organisms that rarely invade the bloodstream to high for classically invasive pathogens:

| Bacterium | Sepsis baseline | Clinical rationale |
|-----------|----------------|-------------------|
| *S. aureus* | −7.3 | Aggressive bloodstream pathogen; 20–30% mortality in bacteraemia (Tong SYC et al., 2015) |
| *P. aeruginosa* | −6.5 | High mortality in ICU infections; often in immunocompromised hosts (Bassetti M et al., 2018) |
| *S. agalactiae* | −7.0 | Neonatal and pregnancy-associated sepsis (Seale AC et al., 2010) |
| *S. pyogenes* | −7.0 | Invasive GAS disease including necrotising fasciitis and toxic shock (Carapetis JR et al., 2005) |
| *E. faecium* | −7.0 | Hospital-acquired bloodstream infections, especially VRE |
| *K. pneumoniae* | −7.5 | Gram-negative sepsis; carbapenem-resistant strains carry >40% mortality (Xu L et al., 2017) |
| *E. faecalis* | −7.5 | Endocarditis and line-related bacteraemia |
| *E. coli* | −9.5 | Most common Gram-negative bloodstream isolate; UTI-source sepsis usually less severe (Poolman JT et al., 2016) |
| *S. Paratyphi A* | −9.2 | Enteric fever with occasional septic complications |
| *iNTS* | −9.0 | Invasive non-typhoidal salmonellosis; high mortality in sub-Saharan Africa (Stanaway JD et al., 2017) |
| *Y. enterocolitica* | −10.0 | Rare sepsis, mainly in iron-overload or immunosuppressed patients |
| *C. trachomatis* | −19.0 | STI — essentially never causes sepsis |
| *N. gonorrhoeae* | −21.0 | Disseminated gonococcal infection is exceedingly rare |
| *H. pylori* | −250.0 | Gastric pathogen — deaths from cancer, not sepsis |
| Fallback default | −14.0 | Used for any organism without an explicit override |



These per-bacterium sepsis baselines are qualitative severity orderings anchored to widely observed differences between invasive and non-invasive pathogens, not claims of portable case-fatality estimates across all settings. Real-world sepsis mortality depends heavily on time-to-treatment, ICU access, comorbidity structure, and health-system capacity, so the model uses these terms mainly to maintain defensible ranking and then lets care access, treatment effectiveness, and syndrome site shape realised mortality in each branch (Rudd KE et al., 2017; Murray CJL et al., 2019).


### 10.3 Non-sepsis infection death

Not all infection-related deaths involve sepsis. Many pathogens kill through tissue-specific mechanisms: *V. cholerae* through fatal dehydration (Ali M et al., 2015), *B. pertussis* through infantile respiratory failure (Yeung KHT et al., 2017), *H. pylori* through gastric adenocarcinoma (Plummer M et al., 2015), *T. pallidum* through tertiary and congenital syphilis (Korenromp EL et al., 2012), and *C. difficile* through toxic megacolon (Guh AY et al., 2020). These deaths would not be captured by the sepsis pathway alone.

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
| *S. epidermidis* | −6.0 | Very low direct mortality; primarily a device-associated pathogen |
| *S. maltophilia* | −4.0 | Some mortality via pneumonia progression, but limited |
| *B. pertussis* | +4.0 | Deaths from respiratory failure in infants, not sepsis (Yeung KHT et al., 2017) |
| *T. pallidum* | +3.5 | Tertiary/congenital syphilis deaths (Korenromp EL et al., 2012) |
| *V. cholerae* | +2.5 | Death from dehydration, not bacteraemia (Ali M et al., 2015) |
| *C. difficile* | +2.0 | Colitis and toxic megacolon deaths (Guh AY et al., 2020) |
| *S. pyogenes* | +1.5 | Rheumatic heart disease and post-streptococcal complications (Watkins DA et al., 2015) |
| *B. fragilis* | +1.5 | Intra-abdominal abscess mortality |
| *H. pylori* | +1.0 | Gastric cancer deaths; essentially zero sepsis risk (Plummer M et al., 2015) |
| *Shigella* spp. | +1.0 | Dysentery deaths in children; sepsis pathway contributes minimally (Troeger C et al., 2016) |



This dual-pathway design ensures that the model can reproduce both the typical sepsis mortality pattern (where broad-spectrum antibiotics and ICU care determine survival) and the non-sepsis mortality pattern (where the primary driver may be dehydration, organ-specific damage, or chronic sequelae).

The non-sepsis adjustments are therefore best viewed as compensating structural terms for important death pathways that a pure sepsis model would miss, rather than as direct organism-specific fatality estimates. That is especially important for globally important syndromes such as cholera, pertussis, diarrhoeal disease, and chronic sequelae-associated infections, where the pathway to death is real but not well represented by bloodstream invasion alone.


### 10.4 Infection mortality — syndrome multipliers

Both death pathways are modulated by the anatomical site of infection. The syndrome multipliers reflect how dangerous each body site is:

| Syndrome | Multiplier | Rationale |
|----------|-----------|-----------|
| Genital | 0.05 | Rarely fatal (localised mucosal infections) |
| Skin / Ear | 0.1 | Low systemic risk unless secondary bacteraemia |
| UTI | 0.5 | Usually self-limiting but can ascend to urosepsis |
| Bone/Joint | 0.8 | Serious but slow-progressing; mortality from surgical complications |
| Intra-abdominal | 1.5 | Peritonitis carries high mortality even with surgery |
| Respiratory | 1.5 | Pneumonia — leading infectious cause of death globally (GBD  Lower Respiratory Infections Collaborators, 2019) |
| CNS | 3.0 | Meningitis/brain abscess — poor penetration of many antibiotics (Tunkel AR et al., 2004) |
| Bloodstream | 4.0 | Bacteraemia/sepsis — the most immediately life-threatening |



These syndrome multipliers are deliberately qualitative. They encode the broad global ordering in which bloodstream and CNS infections are most lethal, respiratory and intra-abdominal infections are high-risk, and genital or superficial infections are usually much less fatal unless they progress, which is the main pattern needed for policy comparisons in the model.



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

That causal attribution holds within the model's own structure: the branches are counterfactual experiments on the simulation, not claims that the listed intervention multipliers are directly transportable policy effect sizes in every real-world setting.


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



These policy parameters are scenario levers rather than empirically fixed intervention coefficients. They are designed to represent the direction and approximate strength of stewardship packages, diagnostic expansion, and resistance-removal counterfactuals so that branch comparisons answer policy questions about mechanism and tradeoff, not provide a literal forecast for any single program design.


### 11.3 Key simulation constants

| Constant | Value | Purpose |
|----------|-------|---------|
| `SIMULATION_START_YEAR` | 1930.0 | Calendar year at day 0 — early enough to capture the pre-antibiotic era |
| `POLICY_BRANCH_YEAR` | 2027.0 | Year when the three policy branches diverge |
| `INFECTION_EPS` | 0.001 | Minimum meaningful infection level (below this, the infection is treated as cleared) |
| `MICROBIOME_MAJORITY_THRESHOLD` | 0.5 | If >50% of a species in the microbiome carries a resistance mechanism, it is classified as "majority resistant" |
| `MAX_MECHANISM_PROFILES` | 200 | Reservoir sample size per bacteria for mechanism profile caching (performance optimisation) |



These constants are internal modeling choices chosen for numerical stability, interpretability, and runtime feasibility. They should not be read as externally validated biological thresholds unless explicitly noted elsewhere.

---



## 12. Known Limitations

Every model is a simplification of reality. These are the main areas where this model knowingly trades accuracy for tractability:

Several of the appendices that follow list exact configuration values and enum definitions. Those tables are included for transparency and reproducibility, but they should still be read in the context established above: many entries are implementation defaults, calibration targets, or structural coding choices rather than direct empirical measurements.

1. **Abstract drug levels**: Antibiotic concentrations are modelled as dimensionless units rather than true pharmacokinetic concentrations (mg/L). This means we capture the *relative* dynamics of drug accumulation and clearance, but cannot directly compare model values to MIC breakpoints from a clinical microbiology report.

2. **No explicit strain competition**: Within the microbiome, resistant and susceptible strains do not explicitly compete for resources. Instead, resistance is promoted by antibiotic pressure and decays through reversion rates. This means the model cannot capture scenarios where a fitness-cost-free resistant strain permanently outcompetes susceptible strains in the absence of antibiotics.

3. **No within-host spatial structure**: Infections are treated as homogeneous within a body compartment. Biofilm formation, abscess walling-off, and planktonic-vs-sessile distinctions are not modelled. In reality, biofilm-embedded bacteria can survive antibiotic concentrations 100–1000× higher than planktonic cells.

4. **Static vaccine model**: Vaccinated individuals have a fixed proportional reduction in infection risk. Vaccine effects do not depend on background prevalence (no herd immunity dynamics), and vaccine-driven serotype replacement is not captured.

5. **Broad regional groupings**: The model uses continental-level regions (e.g., "Europe", "Africa") rather than country-level or hospital-level variation. Antibiotic consumption patterns and resistance rates can vary dramatically between countries within the same region.

---



## Appendix A — Bacteria, Drugs, Mechanisms and Enums

This appendix lists every entity in the model. Use it as a lookup reference when you encounter a specific bacterium, drug, or mechanism code in the main text.

The appendix is intentionally implementation-facing. Names, groupings, and enum labels are the simulation's internal vocabulary for representing major clinical categories; they are not meant to imply that every organism, drug, or ecological niche is exhaustively or uniquely represented by a single real-world classification scheme.



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



### A.2 Antibiotics (61 drugs)

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
| aztreonam_avibactam | Novel BL/BLI |
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

Unless a section above explicitly ties a parameter to a historical date, laboratory quantity, or external benchmark, the values below should be interpreted as the **current default configuration** of the model rather than as a table of direct empirical estimates. Some are mechanistic approximations, some are calibration levers, and some are numerical thresholds chosen to make the simulation stable and interpretable.



### B.1 Global Scalars

These are the ~120 top-level parameters stored in the `GlobalScalars` struct:

For readability, they are grouped by function rather than by epistemic type. In practice this means structurally different things appear side by side: some values are historical switches, some are effective behavioural terms, and some are purely internal numerical controls.



#### Infection acquisition

These are model-scale infection-dynamics defaults, not direct microbiological counts or surveillance incidence rates.

| Parameter | Baseline Value | Description |
|----------|---------|-------------|
| `infection_growth_rate_per_day` | 0.1 | Daily bacterial growth increment |
| `infection_initial_level` | 1.0 | Starting bacterial load |
| `infection_death_threshold` | 50.0 | Level at which death may occur |
| `symptom_onset_threshold` | 3.0 | Level for symptom development |
| `symptom_recheck_interval_days` | 7.0 | Re-evaluation interval |
| `symptom_onset_rate_per_day` | 0.1 | Base symptom development rate |
| `not_under_care_fraction` | 0.05 | Fraction not seeking medical care |

Note: infection resolution in the current rules implementation is operationally handled when infection burden falls to a near-zero level or an immune-clearance event fires, rather than through a standalone configurable `infection_clearance_threshold` parameter.



#### Antibiotic treatment initiation

These terms are best interpreted as **effective prescribing propensities** within the simulation, calibrated to reproduce aggregate treatment rates rather than estimated from a single prescribing dataset.

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

These are simplifying course-control defaults used by the prescribing engine; they are not intended as a complete statistical description of real-world duration distributions.

| Parameter | Baseline Value |
|----------|---------|
| `treatment_stop_improvement_threshold` | 2.0 |
| `treatment_stop_rate_per_day` | 0.03 |
| `treatment_duration_base_days` | 7.0 |



#### Drug efficacy

These values define the model's abstract treatment-response scale and should not be read as pharmacodynamic parameters in standard laboratory units.

| Parameter | Baseline Value |
|----------|---------|
| `drug_effect_on_bacteria_per_day` | 0.5 |
| `drug_minimum_effective_level` | 0.1 |



#### Drug selection

| Parameter | Baseline Value |
|----------|---------|
| `drug_selection_temperature` | 0.55 |
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

This block mixes historically anchored activation dates with effective workflow terms such as testing propensity, turnaround, and reporting reliability.

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



#### Non-sepsis infection death

| Parameter | Baseline Value |
|----------|---------|
| `infection_non_sepsis_base_log_odds` | −9.0 |
| `infection_non_sepsis_log_odds_per_level` | 0.0 |
| `infection_non_sepsis_minimum_bacteria_level` | 0.5 |
| `infection_non_sepsis_log_odds_age_infant` | 0.0 |
| `infection_non_sepsis_log_odds_age_child` | 0.0 |
| `infection_non_sepsis_log_odds_age_adult` | 0.0 |
| `infection_non_sepsis_log_odds_age_elderly` | 0.0 |
| `infection_non_sepsis_log_odds_immunosuppressed` | 0.0 |
| `infection_non_sepsis_log_odds_in_hospital` | 0.0 |



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

These values define the model's admission and discharge process. In particular, `hospitalization_recovery_rate_per_day` is an effective discharge hazard used to match plausible turnover, not a literal summary of all real-world hospital length-of-stay distributions.

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

These parameters mix mechanistic state-transition rates with a deliberately broad vulnerability construct. The age-specific `chronic_immunodeficiency_probability_*` terms are **not literal prevalence estimates of diagnosed immunodeficiency**; they initialise the model's composite higher-risk host state described in Section 2.3.

| Parameter | Baseline Value |
|----------|---------|
| `immunosuppression_startup_seed_fraction` | 0.05 |
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

This section contains implementation-level control parameters for assigning and tracking resistance states.

| Parameter | Baseline Value |
|----------|---------|
| `mechanism_assignment_probability` | 0.8 |



#### Microbiome and carriage

These values act on the model's abstract carriage reservoir and therefore represent effective ecological behaviour within the simulation rather than directly sampled colonisation-study estimates.

| Parameter | Baseline Value |
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
| `mechanism_cache_ewma_decay` | 0.9 |



#### Calibration multipliers

These global scaling parameters control the major axes of resistance dynamics. All default to values that reproduce the baseline calibration; they can be varied independently to explore alternative calibrations or counterfactual scenarios.

| Parameter | Default | What it scales |
|----------|---------|---------------|
| `mechanism_reversion_rate_global_multiplier` | 1.0 | All per-mechanism reversion rates (Section 7.4) |
| `infection_de_novo_multiplier` | 1.0 | De-novo resistance emergence in active infections (Section 7.3) |
| `microbiome_de_novo_multiplier` | 1.0 | De-novo resistance emergence in microbiome carriage (Section 8.2) |
| `hgt_multiplier` | 1.0 | All horizontal gene transfer rates (Section 9.2) |



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

`travel_probability_per_day` is an effective inter-regional mixing parameter. It is intentionally lower-dimensional than real passenger mobility data and is used to reproduce long-run resistant-strain importation pressure rather than literal trip counts.

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
| `microbiome_clearance_probability_per_day` | Bacteria-specific daily carriage-clearance rate from the microbiome reservoir |
| `growth_rate_multiplier` | Bacteria-specific growth modifier |
| `symptom_onset_override` | Override for symptom onset threshold |
| `hgt_donor_probability` | Probability of being an HGT donor |
| `resistance_floor_enabled` | Whether resistance floors apply |
| `resistance_floor_ramp_years` | Ramp-up period for floors |
| `{drug_class}_resistance_floor` | Per-class floor target |



### B.3 Per-Drug Parameters (61 drugs × N parameters)

Generated with key pattern `drug_{name}_{param}`:

| Parameter suffix | Default | Description |
|------------------|---------|-------------|
| `half_life_days` | Drug-specific | PK half-life |
| `initial_level` | 10.0 | Administration level |
| `double_dose_multiplier` | 2.0 | Double-dose level |
| `spectrum_breadth` | 3.0 | Microbiome disruption breadth |
| `microbiome_disruption_log_odds` | 0.3 | Per-drug contribution to the persistent microbiome-disruption reservoir |
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



### B.8 Drug Penetration Table (10 syndromes × 39 drug classes)

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
| `{bacteria}_{drug}_activity_r` | Mean resistance (activity_r, derived from mechanism flags) |
| `{bacteria}_{drug}_majority_r` | Population-level resistance prevalence (derived from mechanism EWMA) |



#### Per-region columns (~6 each)

| Pattern | Description |
|---------|-------------|
| `{region}_infected` | Regional infection count |
| `{region}_hospitalized` | Regional hospital count |
| `{region}_deaths` | Regional death count |



### C.4 Total Column Count

With 42 bacteria, 61 drugs, and 6 regions, the CSV contains approximately:

- ~16 scalar columns
- ~336 per-bacteria columns (42 × 8)
- ~122 per-drug columns (61 × 2)
- ~5,124 per-bacteria-per-drug columns (42 × 61 × 2)
- ~18 per-region columns (6 × 3)
- **Total: ~5,616 columns**



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

- Bratzler DW, Dellinger EP, Olsen KM, et al. Clinical practice guidelines for antimicrobial prophylaxis in surgery. *Am J Health-Syst Pharm.* 2013;70(3):195–283. doi:10.2146/ajhp120568

- Carapetis JR, Steer AC, Mulholland EK, Weber M. The global burden of group A streptococcal diseases. *Lancet Infect Dis.* 2005;5(11):685–694. doi:10.1016/S1473-3099(05)70267-X

- Drlica K, Zhao X. Mutant selection window hypothesis updated. *Clin Infect Dis.* 2007;44(5):681–688. doi:10.1086/511025

- Evans L, Rhodes A, Alhazzani W, et al. Surviving sepsis campaign: international guidelines for management of sepsis and septic shock 2021. *Intensive Care Med.* 2021;47(11):1181–1247. doi:10.1007/s00134-021-06506-y

- Fishman JA. Infection in solid-organ transplant recipients. *N Engl J Med.* 2007;357(25):2601–2614. doi:10.1056/NEJMra064928

- Fleming-Dutra KE, Hersh AL, Shapiro DJ, et al. Prevalence of inappropriate antibiotic prescriptions among US ambulatory care visits, 2010–2011. *JAMA.* 2016;315(17):1864–1873. doi:10.1001/jama.2016.4151

- GBD 2019 Lower Respiratory Infections Collaborators. Age-sex differences in the global burden of lower respiratory infections and risk factors, 1990–2019: results from the Global Burden of Disease Study 2019. *Lancet Infect Dis.* 2022;22(11):1626–1647. doi:10.1016/S1473-3099(22)00510-2

- Guh AY, Mu Y, Winston LG, et al. Trends in U.S. burden of *Clostridioides difficile* infection and outcomes. *N Engl J Med.* 2020;382(14):1320–1330. doi:10.1056/NEJMoa1910215

- Gupta K, Hooton TM, Naber KG, et al. International clinical practice guidelines for the treatment of acute uncomplicated cystitis and pyelonephritis in women: a 2010 update by the Infectious Diseases Society of America and the European Society for Microbiology and Infectious Diseases. *Clin Infect Dis.* 2011;52(5):e103–e120. doi:10.1093/cid/ciq257

- UN Tourism. *UN Tourism Data Dashboard: Global and Regional Tourism Performance.* 2025. Accessed March 24, 2026. https://www.untourism.int/tourism-data/global-and-regional-tourism-performance

- World Bank. *International tourism, number of departures (ST.INT.DPRT).* World Development Indicators; source: UN Tourism. Accessed March 24, 2026. https://data.worldbank.org/indicator/ST.INT.DPRT

- World Bank. *Air transport, passengers carried (IS.AIR.PSGR).* World Development Indicators; source: International Civil Aviation Organization (ICAO). Accessed March 24, 2026. https://data.worldbank.org/indicator/IS.AIR.PSGR

- Jacobs J, Hardy L, Semret M, et al. Diagnostic bacteriology in district hospitals in sub-Saharan Africa: at the forefront of the containment of antimicrobial resistance. *Front Med (Lausanne).* 2019;6:205. doi:10.3389/fmed.2019.00205

- Koning S, van der Sande R, Verhagen AP, et al. Interventions for impetigo. *Cochrane Database Syst Rev.* 2012;(1):CD003261. doi:10.1002/14651858.CD003261.pub3

- Korenromp EL, Rowley J, Alonso M, et al. Global burden of maternal and congenital syphilis and associated adverse birth outcomes — Estimates for 2016 and progress since 2012. *PLoS One.* 2019;14(2):e0211720. doi:10.1371/journal.pone.0211720

- Lee CF, Cowling BJ, Feng S, et al. Impact of antibiotic stewardship programmes in Asia: a systematic review and meta-analysis. *J Antimicrob Chemother.* 2018;73(4):844–851. doi:10.1093/jac/dkx492

- Levy MM, Dellinger RP, Townsend SR, et al. The Surviving Sepsis Campaign: results of an international guideline-based performance improvement program targeting severe sepsis. *Intensive Care Med.* 2010;36(2):222–231. doi:10.1007/s00134-009-1738-3

- Li J, Nation RL, Turnidge JD, et al. Colistin: the re-emerging antibiotic for multidrug-resistant Gram-negative bacterial infections. *Lancet Infect Dis.* 2006;6(9):589–601. doi:10.1016/S1473-3099(06)70580-1

- Magill SS, O'Leary E, Janelle SJ, et al. Changes in prevalence of health care–associated infections in U.S. hospitals. *N Engl J Med.* 2018;379(18):1732–1744. doi:10.1056/NEJMoa1801550

- Martinson ML, Lapham J. Prevalence of immunosuppression among US adults. *JAMA.* 2024;331(10):880–882. doi:10.1001/jama.2023.28019

- United Nations Department of Economic and Social Affairs, Population Division. *World Population Prospects 2024.* Accessed March 24, 2026. https://population.un.org/wpp/

- World Bank. *Hospital beds (per 1,000 people) (SH.MED.BEDS.ZS).* World Development Indicators; source: World Health Organization. Accessed March 24, 2026. https://data.worldbank.org/indicator/SH.MED.BEDS.ZS

- World Health Organization. *Universal health coverage (UHC).* Fact sheet, 2025. Accessed March 24, 2026. https://www.who.int/news-room/fact-sheets/detail/universal-health-coverage-(uhc)

- World Health Organization. *Global Antimicrobial Resistance and Use Surveillance System (GLASS).* Accessed March 24, 2026. https://www.who.int/initiatives/glass

- McInnes RS, McCallum GE, Lamberte LE, van Schaik W. Horizontal transfer of antibiotic resistance genes in the human gut microbiome. *Curr Opin Microbiol.* 2020;53:35–43. doi:10.1016/j.mib.2020.02.002

- Metlay JP, Waterer GW, Long AC, et al. Diagnosis and treatment of adults with community-acquired pneumonia: an official clinical practice guideline of the American Thoracic Society and Infectious Diseases Society of America. *Am J Respir Crit Care Med.* 2019;200(7):e45–e67. doi:10.1164/rccm.201908-1581ST

- Murray CJL, Ikuta KS, Sharara F, et al. Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *Lancet.* 2022;399(10325):629–655. doi:10.1016/S0140-6736(21)02724-0

- Plummer M, Franceschi S, Vignat J, Forman D, de Martel C. Global burden of gastric cancer attributable to *Helicobacter pylori*. *Int J Cancer.* 2015;136(2):487–490. doi:10.1002/ijc.28999

- Poolman JT, Wacker M. Extraintestinal pathogenic *Escherichia coli*, a common human pathogen: challenges for vaccine development and progress in the field. *J Infect Dis.* 2016;213(1):6–13. doi:10.1093/infdis/jiv429

- Rhodes A, Evans LE, Alhazzani W, et al. Surviving Sepsis Campaign: international guidelines for management of sepsis and septic shock: 2016. *Intensive Care Med.* 2017;43(3):304–377. doi:10.1007/s00134-017-4683-6

- Rudd KE, Johnson SC, Agesa KM, et al. Global, regional, and national sepsis incidence and mortality, 1990–2017: analysis for the Global Burden of Disease Study. *Lancet.* 2020;395(10219):200–211. doi:10.1016/S0140-6736(19)32989-7

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

- van de Beek D, Brouwer MC, Thwaites GE, Tunkel AR. Advances in treatment of bacterial meningitis. *Lancet.* 2012;380(9854):1693–1702. doi:10.1016/S0140-6736(12)61186-6

- van Schaik W. The human gut resistome. *Philos Trans R Soc Lond B Biol Sci.* 2015;370(1670):20140087. doi:10.1098/rstb.2014.0087

- Watkins DA, Johnson CO, Colquhoun SM, et al. Global, regional, and national burden of rheumatic heart disease, 1990–2015. *N Engl J Med.* 2017;377(8):713–722. doi:10.1056/NEJMoa1603693

- Werner G, Coque TM, Hammerum AM, et al. Emergence and spread of vancomycin resistance among enterococci in Europe. *Euro Surveill.* 2008;13(47):19046.

- Workowski KA, Bachmann LH, Chan PA, et al. Sexually transmitted infections treatment guidelines, 2021. *MMWR Recomm Rep.* 2021;70(4):1–187. doi:10.15585/mmwr.rr7004a1

- Xu L, Sun X, Ma X. Systematic review and meta-analysis of mortality of patients infected with carbapenem-resistant *Klebsiella pneumoniae*. *Ann Clin Microbiol Antimicrob.* 2017;16(1):18. doi:10.1186/s12941-017-0191-3

- Yeung KHT, Duclos P, Nelson EAS, Hutubessy RCW. An update of the global burden of pertussis in children younger than 5 years: a modelling study. *Lancet Infect Dis.* 2017;17(9):974–980. doi:10.1016/S1473-3099(17)30390-0
