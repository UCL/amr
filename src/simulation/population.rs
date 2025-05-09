use rand::prelude::*;
use std::collections::HashMap;

pub const BACTERIA_LIST: &[&str] = &[
    "acinetobac_bau", "citrobac_spec", "enterobac_spec", "enterococ_faeca", "enterococ_faeci",
    "esch_coli", "kleb_pneu", "morg_spec", "prot_spec", "serrat_spec", "pseud_aerug", "staph_aureus",
    "strep_pneu", "salm_typhi", "salm_parat_a", "inv_nt_salm", "shig_spec", "n_gonorrhoeae",
    "group_a_strep", "group_b_strep", "haem_infl",
];

pub const DRUG_SHORT_NAMES: &[&str] = &[
    "penicilling", "ampicillin", "amoxicillin",
    "piperacillin", "ticarcillin", "cephalexin", "cefazolin",
    "cefuroxime", "ceftriaxone", "ceftazidime", "cefepime", "ceftaroline", "meropenem", "imipenem_c",
    "ertapenem", "aztreonam", "erythromycin", "azithromycin", "clarithromycin", "clindamycin",
    "gentamicin", "tobramycin", "amikacin", "ciprofloxacin", "levofloxacin", "moxifloxacin",
    "ofloxacin", "tetracycline", "doxyclycline", "minocycline", "vancomycin", "teicoplanin",
    "linezolid", "tedizolid", "quinu_dalfo", "trim_sulf", "chlorampheni", "nitrofurantoin",
    "retapamulin", "fusidic_a", "metronidazole",
];

#[derive(Clone, Copy, Debug)]
pub struct DrugResistance {
    pub microbiome_r: f64,
    pub test_r: f64,
    pub activity_r: f64,
    pub e_r: f64,
    pub c_r: f64,
}

#[derive(Clone, Debug)]
pub struct Individual {
    pub id: usize,
    pub age: usize,
    pub sex_at_birth: String,
    pub resistances: Vec<Vec<DrugResistance>>, // [bacteria][drug]

    // Explicitly named variables
    pub date_last_infected: HashMap<&'static str, i32>, // Days since last infection (integer)
    pub infectious_syndrome: HashMap<&'static str, f64>, // Level of infectious syndrome (0.0 - 1.0)
    pub level: HashMap<&'static str, f64>,                 // Level of bacteria (e.g., count or concentration)
    pub immune_resp: HashMap<&'static str, f64>,           // Immune response level (0.0 - 1.0)
    pub sepsis: HashMap<&'static str, bool>,               // Is the individual septic for this bacteria?
    pub level_microbiome: HashMap<&'static str, f64>,       // Level in the microbiome (e.g., proportion)

    pub haem_infl_vaccination_status: bool,             // true if vaccinated, false otherwise
    pub strep_pneu_vaccination_status: bool,             // true if vaccinated, false otherwise
    pub salm_typhi_vaccination_status: bool,             // true if vaccinated, false otherwise
    pub esch_coli_vaccination_status: bool,              // true if vaccinated, false otherwise

    pub cur_use_drug: Vec<bool>,                         // true if currently using the drug
    pub cur_level_drug: Vec<f64>,                        // current level/dosage of the drug

    pub current_infection_related_death_risk: f64,       // Probability (0.0 - 1.0)
    pub background_all_cause_mortality_rate: f64,         // Rate (e.g., per year)
    pub sexual_contact_level: f64,                      // Frequency or intensity
    pub airborne_contact_level_with_adults: f64,        // Frequency or intensity
    pub airborne_contact_level_with_children: f64,      // Frequency or intensity
    pub oral_exposure_level: f64,                       // Frequency or intensity
    pub mosquito_exposure_level: f64,                   // Frequency or intensity
    pub under_care: bool,                               // true if currently under medical care
    pub infection_hospital_acquired: bool,              // true if infection was hospital-acquired
    pub current_toxicity: f64,                          // Level of current toxicity
    pub mortality_risk_current_toxicity: f64,           // Probability (0.0 - 1.0)
}

impl Individual {
    pub fn new(id: usize, age: usize, sex_at_birth: String) -> Self {
        let mut resistances = Vec::new();
        let num_drugs = DRUG_SHORT_NAMES.len();

        for _bacteria in BACTERIA_LIST.iter() {
            let mut drug_resistances = Vec::new();
            for _ in 0..num_drugs {
                drug_resistances.push(DrugResistance {
                    microbiome_r: 0.0, // Initializing with a default value
                    test_r: 0.0,
                    activity_r: 0.0,
                    e_r: 0.0,
                    c_r: 0.0,
                });
            }
            resistances.push(drug_resistances);
        }

        let mut individual = Individual {
            id,
            age,
            sex_at_birth,
            resistances,
            date_last_infected: HashMap::new(),
            infectious_syndrome: HashMap::new(),
            level: HashMap::new(),
            immune_resp: HashMap::new(),
            sepsis: HashMap::new(),
            level_microbiome: HashMap::new(),
            haem_infl_vaccination_status: false,
            strep_pneu_vaccination_status: false,
            salm_typhi_vaccination_status: false,
            esch_coli_vaccination_status: false,
            cur_use_drug: vec![false; num_drugs],
            cur_level_drug: vec![0.0; num_drugs],
            current_infection_related_death_risk: 0.0,
            background_all_cause_mortality_rate: 0.0,
            sexual_contact_level: 0.0,
            airborne_contact_level_with_adults: 0.0,
            airborne_contact_level_with_children: 0.0,
            oral_exposure_level: 0.0,
            mosquito_exposure_level: 0.0,
            under_care: false,
            infection_hospital_acquired: false,
            current_toxicity: 0.0,
            mortality_risk_current_toxicity: 0.0,
        };
        individual.initialize_state();
        individual
    }

    fn initialize_state(&mut self) {
        let mut rng = rand::thread_rng();

        // Initialize HashMap variables
        for &bacteria in BACTERIA_LIST.iter() {
            self.date_last_infected.insert(bacteria, 0); // Initializing as 0 days since last infection (integer)
            self.infectious_syndrome.insert(bacteria, rng.gen_range(0.0..=1.0));
            self.level.insert(bacteria, rng.gen_range(0.0..1000.0)); // Example level
            self.immune_resp.insert(bacteria, rng.gen_range(0.0..=1.0));
            self.sepsis.insert(bacteria, rng.gen_bool(0.1)); // 10% chance of initial sepsis
            self.level_microbiome.insert(bacteria, rng.gen_range(0.0..=1.0));
        }

        // Initialize boolean vaccination statuses
        self.haem_infl_vaccination_status = rng.gen_bool(0.5);
        self.strep_pneu_vaccination_status = rng.gen_bool(0.5);
        self.salm_typhi_vaccination_status = rng.gen_bool(0.5);
        self.esch_coli_vaccination_status = rng.gen_bool(0.5);

        // Initialize current drug use and level (initially mostly false/zero)
        for i in 0..self.cur_use_drug.len() {
            self.cur_use_drug[i] = rng.gen_bool(0.05); // 5% chance of initial drug use
            self.cur_level_drug[i] = if self.cur_use_drug[i] { rng.gen_range(0.1..1.0) } else { 0.0 };
        }

        // Initialize other scalar variables
        self.current_infection_related_death_risk = rng.gen_range(0.0..=0.1); // Low initial risk
        self.background_all_cause_mortality_rate = rng.gen_range(0.001..=0.01); // Example rate
        self.sexual_contact_level = rng.gen_range(0.0..5.0);
        self.airborne_contact_level_with_adults = rng.gen_range(0.0..10.0);
        self.airborne_contact_level_with_children = rng.gen_range(0.0..10.0);
        self.oral_exposure_level = rng.gen_range(0.0..3.0);
        self.mosquito_exposure_level = rng.gen_range(0.0..2.0);
        self.under_care = rng.gen_bool(0.2); // 20% chance of being under care
        self.infection_hospital_acquired = rng.gen_bool(0.05); // 5% chance of hospital-acquired infection
        self.current_toxicity = rng.gen_range(0.0..=0.2); // Low initial toxicity
        self.mortality_risk_current_toxicity = rng.gen_range(0.0..=0.05); // Low initial risk
    }

    pub fn get_resistance(&self, bacteria: &str, drug: &str, resistance_type: &str) -> Option<f64> {
        if let (Some(bi), Some(di)) = (
            BACTERIA_LIST.iter().position(|&x| x == bacteria),
            DRUG_SHORT_NAMES.iter().position(|&x| x == drug),
        ) {
            let r = &self.resistances[bi][di];
            match resistance_type {
                "e_r" => Some(r.e_r),
                "c_r" => Some(r.c_r),
                "microbiome_r" => Some(r.microbiome_r),
                "test_r" => Some(r.test_r),
                "activity_r" => Some(r.activity_r),
                _ => None,
            }
        } else {
            None
        }
    }
}

#[derive(Debug)]
pub struct Population {
    pub individuals: Vec<Individual>,
    pub num_other_variables: usize, // Still keeping this for now, might be removed later
}

impl Population {
    pub fn new(num_individuals: usize) -> Self {
        let mut rng = rand::thread_rng();
        let individuals = (0..num_individuals)
            .map(|i| {
                let age = rng.gen_range(0..100);
                let sex = if rng.gen_bool(0.5) { "male" } else { "female" }.to_string();
                Individual::new(i, age, sex) // Corrected call: removed the extra argument
            })
            .collect();
        Population { individuals, num_other_variables: 0 }
    }

    pub fn update(&mut self, _time_step: usize, _all_rules: Vec<Vec<(usize, f64)>>) {
        // This function might need to be removed or significantly changed
        // as the rules are now applied directly in apply_rules.
        // Keeping it as a placeholder for now in case there's other logic.
        eprintln!("Warning: Population::update is called but might not be needed anymore.");
    }

    pub fn calculate_summary_statistics(&self) -> HashMap<String, f64> {
        let mut summary = HashMap::new();
        let num_individuals = self.individuals.len() as f64;

        if num_individuals == 0.0 {
            return summary;
        }

        let average_age: f64 = self.individuals.iter().map(|ind| ind.age as f64).sum::<f64>() / num_individuals;
        summary.insert("average_age".to_string(), average_age);

        let male_count = self.individuals.iter().filter(|ind| ind.sex_at_birth == "male").count() as f64;
        summary.insert("proportion_male".to_string(), male_count / num_individuals);

        summary // Returns the HashMap containing the summary statistics.
    }

    pub fn get_individual(&self, index: usize) -> Option<&Individual> {
        self.individuals.get(index)
    }

    pub fn get_num_individuals(&self) -> usize {
        self.individuals.len()
    }
}

