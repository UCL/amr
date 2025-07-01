// src/simulation/population.rs
use rand::Rng;
use rand::distributions::{Distribution, Standard};
use std::fmt; 


pub const BACTERIA_LIST: &[&str] = &[
    "acinetobac_bau", "citrobac_spec", "enterobac_spec", "enterococ_faeca", "enterococ_faeci",
    "esch_coli", "kleb_pneu", "morg_spec", "prot_spec", "serrat_spec", "pseud_aerug", "staph_aureus",
    "strep_pneu", "salm_typhi", "salm_parat_a", "inv_nt_salm", "shig_spec", "n_gonorrhoeae",
    "group_a_strep", "group_b_strep", "haem_infl", "chlam_trach", 
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


#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HospitalStatus {
    InHospital,  // consider in future whether to have a variable for whether in icu
    NotInHospital,
}

impl HospitalStatus {
    pub fn is_hospitalized(&self) -> bool {
        matches!(self, HospitalStatus::InHospital)
    }
}

// Add Display to the derive attribute and implement it
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Region {
    NorthAmerica,
    SouthAmerica,
    Africa,
    Asia,
    Europe,
    Oceania,
    Home, // This represents the individual's home region, which could be any of the above.
}

impl Distribution<Region> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Region {
        match rng.gen_range(0..6) { // 0 to 5 for the 6 geographic regions
            0 => Region::NorthAmerica,
            1 => Region::SouthAmerica,
            2 => Region::Africa,
            3 => Region::Asia,
            4 => Region::Europe,
            _ => Region::Oceania, // Default for 5
        }
    }
}

// Implement the Display trait for Region
impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use debug format and then convert to lowercase and replace spaces
        write!(f, "{:?}", self) // This will give "NorthAmerica", "SouthAmerica", etc.
    }
}


#[derive(Debug, Clone)]
pub struct Resistance {
    pub microbiome_r: f64,
    pub test_r: f64,
    pub activity_r: f64,
    pub any_r: f64, // Effective Resistance in minority or majority (0-1)
    pub majority_r: f64, // Resistance in majority of bacteria infected with (0-1) - when majority_r is non zero 
                         // it will always take the same value as any_r
}

#[derive(Debug, Clone)]
pub struct Individual {
    pub id: usize,
    pub age: i32, // age in days (negative = not yet born date)
    pub sex_at_birth: String,
    pub region_living: Region,
    pub region_cur_in: Region,
    pub days_visiting: u32, 
    pub hospital_status: HospitalStatus,
    pub days_hospitalized: u32, 
    pub date_last_infected: Vec<i32>,              
    pub infectious_syndrome: Vec<i32>,             
    pub level: Vec<f64>,

    pub immune_resp: Vec<f64>,                     
    pub sepsis: Vec<bool>,                         
    pub presence_microbiome: Vec<bool>,            
    pub vaccination_status: Vec<bool>,             
    pub cur_infection_from_environment: Vec<bool>, 
    pub test_identified_infection: Vec<bool>,      
    pub cur_use_drug: Vec<bool>,
    pub cur_level_drug: Vec<f64>,  // standard level is 10 for a day on which a standard dose is taken / administered 
    pub date_drug_initiated: Vec<i32>, // the time_step when each drug was last initiated
    pub current_infection_related_death_risk: f64,
    pub background_all_cause_mortality_rate: f64,  
    pub sexual_contact_level: f64,
    pub airborne_contact_level_with_adults: f64,
    pub airborne_contact_level_with_children: f64,
    pub oral_exposure_level: f64,
    pub mosquito_exposure_level: f64,
    pub infection_hospital_acquired: Vec<bool>,    
    pub current_toxicity: f64,
    pub mortality_risk_current_toxicity: f64, 
    pub resistances: Vec<Vec<Resistance>>,
    pub date_of_death: Option<usize>,
    pub cause_of_death: Option<String>,
    pub is_severely_immunosuppressed: bool, 

}

impl Individual {
    pub fn new(id: usize, age_days: i32, sex_at_birth: String) -> Self {
        let mut rng = rand::thread_rng();
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();

        let date_last_infected = vec![0; num_bacteria];
        let infectious_syndrome = vec![0; num_bacteria];
        let level = vec![0.0; num_bacteria];
        let immune_resp = vec![0.1; num_bacteria];
        let sepsis = vec![false; num_bacteria];
        let presence_microbiome = vec![false; num_bacteria];
        let infection_hospital_acquired = vec![false; num_bacteria];
        let cur_infection_from_environment = vec![false; num_bacteria];
        let test_identified_infection = vec![false; num_bacteria];
        let vaccination_status = (0..num_bacteria).map(|_| rng.gen_bool(0.5)).collect();

        let mut resistances = Vec::with_capacity(num_bacteria);
        for _ in 0..num_bacteria {
            let mut drug_resistances = Vec::with_capacity(num_drugs);
            for _ in 0..num_drugs {
                drug_resistances.push(Resistance {
                    microbiome_r: 0.0,
                    test_r: 0.0,
                    activity_r: 0.0,
                    any_r: 0.0,
                    majority_r: 0.0,
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
            region_living: rng.gen(), 
            region_cur_in: Region::Home, 
            days_visiting: 0, 
            hospital_status: HospitalStatus::NotInHospital, 
            days_hospitalized: 0, 
            sex_at_birth,
            date_last_infected,
            infectious_syndrome,
            level, 
            immune_resp,
            sepsis,
            presence_microbiome,
            vaccination_status, 
            cur_use_drug: vec![false; num_drugs],
            cur_level_drug: vec![0.0; num_drugs],
            date_drug_initiated: vec![i32::MIN; num_drugs], 
            current_infection_related_death_risk: 0.0,
            background_all_cause_mortality_rate,  
            sexual_contact_level: rng.gen_range(0.0..=10.0),
            airborne_contact_level_with_adults: rng.gen_range(0.0..=10.0),
            airborne_contact_level_with_children: rng.gen_range(0.0..=10.0),
            oral_exposure_level: rng.gen_range(0.0..=10.0),
            mosquito_exposure_level: rng.gen_range(0.0..=10.0),
            infection_hospital_acquired,
            cur_infection_from_environment,
            test_identified_infection,
            current_toxicity: rng.gen_range(0.0..=3.0),
            mortality_risk_current_toxicity: 0.0, // todo: probably should be removed as this death risk is implemented with separate logic
            resistances,
            date_of_death: None,
            cause_of_death: None,
            is_severely_immunosuppressed: false, 
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
            let age = rng.gen_range(0..=36500); // Age range from 0 to 100 years in days - will need to change 0 to -36500
            let sex = if rng.gen_bool(0.5) { "male".to_string() } else { "female".to_string() };
            individuals.push(Individual::new(i, age, sex));
        }
        Population { individuals }
    }
}