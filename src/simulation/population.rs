// src/simulation/population.rs
use rand::Rng;
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

#[derive(Debug, Clone)]
pub struct Resistance {
    pub microbiome_r: f64,
    pub test_r: f64,
    pub activity_r: f64,
    pub e_r: f64,
    pub c_r: f64,
}

#[derive(Debug, Clone)]
pub struct Individual {
    pub id: usize,
    pub age: i32, // age in days (negative = before reference date)
    pub sex_at_birth: String,
    pub date_last_infected: HashMap<&'static str, i32>,
    pub infectious_syndrome: HashMap<&'static str, i32>, // Changed to i32
    pub level: HashMap<&'static str, f64>,
    pub immune_resp: HashMap<&'static str, f64>,
    pub sepsis: HashMap<&'static str, bool>,
    pub level_microbiome: HashMap<&'static str, f64>,
    pub haem_infl_vaccination_status: bool,
    pub strep_pneu_vaccination_status: bool,
    pub salm_typhi_vaccination_status: bool,
    pub esch_coli_vaccination_status: bool,
    // Note: You might have more vaccination statuses in your actual code
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
    pub resistances: Vec<Vec<Resistance>>,
}

impl Individual {
    pub fn new(id: usize, age_days: i32, sex_at_birth: String) -> Self {
        let mut rng = rand::thread_rng();
        let mut date_last_infected = HashMap::new();
        let mut infectious_syndrome: HashMap<&'static str, i32> = HashMap::new();
        let mut level = HashMap::new();
        let mut immune_resp = HashMap::new();
        let mut sepsis = HashMap::new();
        let mut level_microbiome = HashMap::new();

        for &bacteria in BACTERIA_LIST.iter() {
            date_last_infected.insert(bacteria, 0);
            infectious_syndrome.insert(bacteria, 0);
            level.insert(bacteria, 0.0);
            immune_resp.insert(bacteria, 0.0);
            sepsis.insert(bacteria, false);
            level_microbiome.insert(bacteria, 0.0);
        }

        let num_drugs = DRUG_SHORT_NAMES.len();
        let num_bacteria = BACTERIA_LIST.len();
        let mut resistances = Vec::with_capacity(num_bacteria);
        for _ in 0..num_bacteria {
            let mut drug_resistances = Vec::with_capacity(num_drugs);
            for _ in 0..num_drugs {
                drug_resistances.push(Resistance {
                    microbiome_r: rng.gen_range(0.0..1.0),
                    test_r: rng.gen_range(0.0..=3.0),
                    activity_r: rng.gen_range(0.0..=3.0),
                    e_r: rng.gen_range(0.0..=1.0),
                    c_r: rng.gen_range(0.0..=1.0),
                });
            }
            resistances.push(drug_resistances);
        }

        let background_all_cause_mortality_rate = if age_days < 0 {
            0.0
        } else {
            0.000001
        };

        Individual {
            id,
            age: age_days,
            sex_at_birth,
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
            cur_use_drug: vec![false; num_drugs],
            cur_level_drug: (0..num_drugs).map(|_| rng.gen_range(0.0..=3.0)).collect(),
            current_infection_related_death_risk: 0.0,
            background_all_cause_mortality_rate,
            sexual_contact_level: rng.gen_range(0.0..=10.0),
            airborne_contact_level_with_adults: rng.gen_range(0.0..=10.0),
            airborne_contact_level_with_children: rng.gen_range(0.0..=10.0),
            oral_exposure_level: rng.gen_range(0.0..=10.0),
            mosquito_exposure_level: rng.gen_range(0.0..=10.0),
            under_care: rng.gen_bool(0.1),
            infection_hospital_acquired: rng.gen_bool(0.05),
            current_toxicity: rng.gen_range(0.0..=3.0),
            mortality_risk_current_toxicity: 0.0,
            resistances,
        }
    }

}

pub struct Population {
    pub individuals: Vec<Individual>,
}

impl Population {
    pub fn new(size: usize) -> Self {
        let mut individuals = Vec::with_capacity(size);
        let mut rng = rand::thread_rng();
        for i in 0..size {
            let age = rng.gen_range(-36500..=36500); // Age range from -100 to +100 years in days
            let sex = if rng.gen_bool(0.5) { "male".to_string() } else { "female".to_string() };
            individuals.push(Individual::new(i, age, sex));
        }
        Population { individuals }
    }
}