// src/simulation/population.rs
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
    pub age: i32,  // age in days (negative = before reference date)
    pub sex_at_birth: String,
    pub resistances: Vec<Vec<DrugResistance>>, // [bacteria][drug]

    // Explicitly named variables
    pub date_last_infected: HashMap<&'static str, i32>,
    pub infectious_syndrome: HashMap<&'static str, f64>,
    pub level: HashMap<&'static str, f64>,
    pub immune_resp: HashMap<&'static str, f64>,
    pub sepsis: HashMap<&'static str, bool>,
    pub level_microbiome: HashMap<&'static str, f64>,

    pub haem_infl_vaccination_status: bool,
    pub strep_pneu_vaccination_status: bool,
    pub salm_typhi_vaccination_status: bool,
    pub esch_coli_vaccination_status: bool,

    pub cur_use_drug: Vec<bool>,
    pub cur_level_drug: Vec<f64>,

    pub current_infection_related_death_risk: f64,
    pub background_all_cause_mortality_rate: f64,
    pub sexual_contact_level: f64,
    pub airborne_contact_level_with_adults: f64,
    pub airborne_contact_level_with_children: f64,
    pub oral_exposure_level: f64,
    pub mosquito_exposure_level: f64,
    pub under_care: bool,
    pub infection_hospital_acquired: bool,
    pub current_toxicity: f64,
    pub mortality_risk_current_toxicity: f64,
}

impl Individual {
    pub fn new(id: usize, age_days: i32, sex_at_birth: String) -> Self {
        let mut rng = rand::thread_rng();
        // Uniform distribution of age in days [-36500,36500]
        let age_days = rng.gen_range(-36500..=36500);

        // Initialize resistances
        let mut resistances = Vec::new();
        for _ in BACTERIA_LIST.iter() {
            let mut drug_resistances = Vec::new();
            for _ in DRUG_SHORT_NAMES.iter() {
                drug_resistances.push(DrugResistance {
                    microbiome_r: rng.gen(),
                    test_r: rng.gen(),
                    activity_r: rng.gen(),
                    e_r: rng.gen(),
                    c_r: rng.gen(),
                });
            }
            resistances.push(drug_resistances);
        }

        // Initialize maps
        let mut date_last_infected = HashMap::new();
        let mut infectious_syndrome = HashMap::new();
        let mut level = HashMap::new();
        let mut immune_resp = HashMap::new();
        let mut sepsis = HashMap::new();
        let mut level_microbiome = HashMap::new();
        for &b in BACTERIA_LIST.iter() {
            date_last_infected.insert(b, 0);
            infectious_syndrome.insert(b, rng.gen());
            level.insert(b, rng.gen());
            immune_resp.insert(b, rng.gen());
            sepsis.insert(b, rng.gen_bool(0.1));
            level_microbiome.insert(b, rng.gen());
        }

        Individual {
            id,
            age: age_days,
            sex_at_birth,
            resistances,
            date_last_infected,
            infectious_syndrome,
            level,
            immune_resp,
            sepsis,
            level_microbiome,
            haem_infl_vaccination_status: rng.gen_bool(0.5),
            strep_pneu_vaccination_status: rng.gen_bool(0.5),
            salm_typhi_vaccination_status: rng.gen_bool(0.5),
            esch_coli_vaccination_status: rng.gen_bool(0.5),
            cur_use_drug: vec![false; DRUG_SHORT_NAMES.len()],
            cur_level_drug: vec![0.0; DRUG_SHORT_NAMES.len()],
            current_infection_related_death_risk: rng.gen(),
            background_all_cause_mortality_rate: rng.gen(),
            sexual_contact_level: rng.gen(),
            airborne_contact_level_with_adults: rng.gen(),
            airborne_contact_level_with_children: rng.gen(),
            oral_exposure_level: rng.gen(),
            mosquito_exposure_level: rng.gen(),
            under_care: rng.gen_bool(0.5),
            infection_hospital_acquired: rng.gen_bool(0.05),
            current_toxicity: rng.gen(),
            mortality_risk_current_toxicity: rng.gen(),
        }
    }
}

#[derive(Debug)]
pub struct Population {
    pub individuals: Vec<Individual>,
}

impl Population {
    pub fn new(num_individuals: usize) -> Self {
        let individuals = (0..num_individuals)
            .map(|i| Individual::new(i, 0, "".to_string()))
            .collect();
        Population { individuals }
    }

    pub fn get_individual(&self, index: usize) -> Option<&Individual> {
        self.individuals.get(index)
    }
}
