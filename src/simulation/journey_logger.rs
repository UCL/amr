//
// Journey Logger for AMR Simulation
//
// This module provides focused logging of individual infection journeys,
        // Journey logging enabled silently capturing detailed daily snapshots only during active infections.
// Much more efficient than full individual logging while providing
// rich data for debugging and patient journey analysis.
//

use crate::simulation::population::{Individual, BACTERIA_LIST, DRUG_SHORT_NAMES};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Write, BufWriter};
use rand::Rng;

#[derive(Clone, Debug)]
pub struct InfectionJourneySnapshot {
    // Journey identification
    pub journey_id: u32,
    pub individual_id: usize,
    pub time_step: usize,
    pub day_of_journey: u32,
    
    // Demographics (set at journey start)
    pub age_at_onset: i32,
    pub sex: String,
    pub region_living: String,
    pub region_current: String,
    pub immunodeficiency: String,
    
    // Primary infection details
    pub primary_bacteria: String,
    pub primary_bacteria_level: f64,
    pub syndrome: i32,
    pub sepsis: bool,
    pub hospital_acquired: bool,
    
    // All active infections
    pub all_bacteria_levels: Vec<(String, f64)>, // (bacteria_name, level)
    
    // Treatment details
    pub current_drugs: Vec<(String, f64)>, // (drug_name, level)
    pub days_on_current_treatment: i32,
    pub treatment_failures_count: u32,
    
    // Resistance status for primary bacteria
    pub resistance_any_r: Vec<(String, f64)>, // (drug_name, any_r_level)
    pub resistance_mechanisms: Vec<String>, // active mechanisms
    
    // Clinical status
    pub hospital_status: String,
    pub immunity_level: f64,
    pub toxicity_level: f64,
    pub background_mortality_risk: f64,
    
    // Testing status
    pub infection_identified: bool,
    pub resistance_testing_done: bool,
    
    // Journey outcome (only on final snapshot)
    pub resolution_type: Option<String>,
}

#[derive(Clone, Debug)]
pub struct ActiveJourney {
    pub journey_id: u32,
    pub individual_id: usize,
    pub onset_time_step: usize,
    pub primary_bacteria_idx: usize,
    pub day_count: u32,
    pub snapshots: Vec<InfectionJourneySnapshot>,
}

pub struct JourneyLogger {
    // Configuration
    pub enabled: bool,
    pub sample_rate: f64,
    
    // Active tracking
    active_journeys: HashMap<usize, ActiveJourney>, // individual_id -> journey
    next_journey_id: u32,
    
    // Output
    csv_writer: Option<BufWriter<File>>,
    output_filename: String,
    
    // Statistics
    pub total_journeys_started: u32,
    pub total_journeys_completed: u32,
    pub total_snapshots_logged: u32,
}

impl JourneyLogger {
    pub fn new() -> Self {
        Self {
            enabled: false,
            sample_rate: 0.0,
            active_journeys: HashMap::new(),
            next_journey_id: 1,
            csv_writer: None,
            output_filename: String::new(),
            total_journeys_started: 0,
            total_journeys_completed: 0,
            total_snapshots_logged: 0,
        }
    }
    
    pub fn enable(&mut self, sample_rate: f64) -> Result<(), Box<dyn std::error::Error>> {
        self.enabled = true;
        self.sample_rate = sample_rate.clamp(0.0, 1.0);
        
        // Create output file
        self.output_filename = "infection_journeys.csv".to_string();
        let file = File::create(&self.output_filename)?;
        let mut writer = BufWriter::new(file);
        
        // Write CSV header
        writeln!(writer, "{}", JourneyLogger::get_csv_header())?;
        
        self.csv_writer = Some(writer);
        
        // Journey logging enabled silently
        
        Ok(())
    }
    
    fn get_csv_header() -> &'static str {
        "journey_id,individual_id,time_step,day_of_journey,age_at_onset,sex,region_living,region_current,immunodeficiency,primary_bacteria,primary_bacteria_level,syndrome,sepsis,hospital_acquired,all_bacteria_levels,current_drugs,days_on_current_treatment,treatment_failures,resistance_any_r,resistance_mechanisms,hospital_status,immunity_level,toxicity_level,background_mortality_risk,infection_identified,resistance_testing_done,resolution_type"
    }
    
    pub fn check_individual(&mut self, individual: &Individual, time_step: usize) {
        if !self.enabled {
            return;
        }
        
        let individual_id = individual.id;
        let has_active_infection = individual.level.iter().any(|&level| level > 0.001);
        let is_currently_tracked = self.active_journeys.contains_key(&individual_id);
        
        match (is_currently_tracked, has_active_infection) {
            (false, true) => {
                // New infection - potentially start tracking
                let mut rng = rand::thread_rng();
                if rng.gen::<f64>() < self.sample_rate {
                    self.start_journey(individual, time_step);
                }
            },
            (true, true) => {
                // Continue existing journey
                self.update_journey(individual, time_step);
            },
            (true, false) => {
                // Journey completed (infection resolved)
                self.complete_journey(individual, time_step);
            },
            (false, false) => {
                // Not tracking, no infection - do nothing
            }
        }
    }
    
    fn start_journey(&mut self, individual: &Individual, time_step: usize) {
        // Find primary bacteria (highest level)
        let mut primary_bacteria_idx = 0;
        let mut highest_level = 0.0;
        
        for (b_idx, &level) in individual.level.iter().enumerate() {
            if level > highest_level {
                highest_level = level;
                primary_bacteria_idx = b_idx;
            }
        }
        
        if highest_level <= 0.001 {
            return; // No significant infection
        }
        
        let journey_id = self.next_journey_id;
        self.next_journey_id += 1;
        
        // Create initial snapshot
        let snapshot = JourneyLogger::create_snapshot(individual, journey_id, time_step, 1, primary_bacteria_idx, None);
        
        let journey = ActiveJourney {
            journey_id,
            individual_id: individual.id,
            onset_time_step: time_step,
            primary_bacteria_idx,
            day_count: 1,
            snapshots: vec![snapshot],
        };
        
        self.active_journeys.insert(individual.id, journey);
        self.total_journeys_started += 1;
        
        // Write initial snapshot to CSV
        if let Some(ref mut writer) = self.csv_writer {
            if let Some(snapshot) = self.active_journeys.get(&individual.id).unwrap().snapshots.last() {
                if let Err(e) = writeln!(writer, "{}", JourneyLogger::snapshot_to_csv(snapshot)) {
                    eprintln!("Error writing snapshot: {}", e);
                } else {
                    self.total_snapshots_logged += 1;
                }
            }
        } else {
            eprintln!("CSV writer is None when trying to write initial snapshot!");
        }
    }
    
    fn update_journey(&mut self, individual: &Individual, time_step: usize) {
        if let Some(journey) = self.active_journeys.get_mut(&individual.id) {
            journey.day_count += 1;
            
            let snapshot = JourneyLogger::create_snapshot(
                individual, 
                journey.journey_id, 
                time_step, 
                journey.day_count, 
                journey.primary_bacteria_idx,
                None
            );
            
            journey.snapshots.push(snapshot);
            
            // Write snapshot to CSV
            if let Some(ref mut writer) = self.csv_writer {
                if let Some(snapshot) = journey.snapshots.last() {
                    if let Err(e) = writeln!(writer, "{}", JourneyLogger::snapshot_to_csv(snapshot)) {
                        eprintln!("Error writing update snapshot: {}", e);
                    } else {
                        self.total_snapshots_logged += 1;
                    }
                }
            } else {
                eprintln!("CSV writer is None when trying to write update snapshot!");
            }
        }
    }
    
    fn complete_journey(&mut self, individual: &Individual, time_step: usize) {
        if let Some(mut journey) = self.active_journeys.remove(&individual.id) {
            journey.day_count += 1;
            
            // Determine resolution type
            let resolution_type = self.determine_resolution_type(individual);
            
            // Create final snapshot with resolution
            let final_snapshot = JourneyLogger::create_snapshot(
                individual,
                journey.journey_id,
                time_step,
                journey.day_count,
                journey.primary_bacteria_idx,
                Some(resolution_type.clone())
            );
            
            journey.snapshots.push(final_snapshot);
            
            // Write final snapshot
            if let Some(ref mut writer) = self.csv_writer {
                if let Some(snapshot) = journey.snapshots.last() {
                    let _ = writeln!(writer, "{}", JourneyLogger::snapshot_to_csv(snapshot));
                    self.total_snapshots_logged += 1;
                }
            }
            
            self.total_journeys_completed += 1;
        }
    }
    
    fn create_snapshot(
        individual: &Individual,
        journey_id: u32,
        time_step: usize,
        day_of_journey: u32,
        primary_bacteria_idx: usize,
        resolution_type: Option<String>
    ) -> InfectionJourneySnapshot {
        
        // Collect all active bacteria
        let all_bacteria_levels: Vec<(String, f64)> = individual.level.iter()
            .enumerate()
            .filter(|(_, &level)| level > 0.001)
            .map(|(idx, &level)| (BACTERIA_LIST[idx].to_string(), level))
            .collect();
        
        // Collect current drugs
        let current_drugs: Vec<(String, f64)> = individual.cur_use_drug.iter()
            .enumerate()
            .filter(|(_, &is_taking)| is_taking)
            .map(|(idx, _)| (DRUG_SHORT_NAMES[idx].to_string(), individual.cur_level_drug[idx]))
            .collect();
        
        // Collect resistance for primary bacteria
        let resistance_any_r: Vec<(String, f64)> = DRUG_SHORT_NAMES.iter()
            .enumerate()
            .filter(|(idx, _)| individual.resistances[primary_bacteria_idx][*idx].any_r > 0.0)
            .map(|(idx, &drug_name)| (drug_name.to_string(), individual.resistances[primary_bacteria_idx][idx].any_r))
            .collect();
        
        // Collect active resistance mechanisms
        let resistance_mechanisms: Vec<String> = individual.resistance_mechanisms[primary_bacteria_idx].iter()
            .enumerate()
            .filter(|(_, &is_active)| is_active)
            .map(|(idx, _)| format!("mechanism_{}", idx)) // Will need to map to actual mechanism names
            .collect();
        
        InfectionJourneySnapshot {
            journey_id,
            individual_id: individual.id,
            time_step,
            day_of_journey,
            age_at_onset: individual.age,
            sex: individual.sex_at_birth.clone(),
            region_living: format!("{:?}", individual.region_living),
            region_current: format!("{:?}", individual.region_cur_in),
            immunodeficiency: format!("{:?}", individual.immunodeficiency_type),
            primary_bacteria: BACTERIA_LIST[primary_bacteria_idx].to_string(),
            primary_bacteria_level: individual.level[primary_bacteria_idx],
            syndrome: individual.infectious_syndrome[primary_bacteria_idx],
            sepsis: individual.sepsis[primary_bacteria_idx],
            hospital_acquired: individual.infection_hospital_acquired[primary_bacteria_idx],
            all_bacteria_levels,
            current_drugs,
            days_on_current_treatment: individual.days_on_current_treatment[primary_bacteria_idx],
            treatment_failures_count: 0, // Will need to track this
            resistance_any_r,
            resistance_mechanisms,
            hospital_status: format!("{:?}", individual.hospital_status),
            immunity_level: individual.immune_resp[primary_bacteria_idx],
            toxicity_level: individual.current_toxicity,
            background_mortality_risk: individual.background_all_cause_mortality_rate,
            infection_identified: individual.test_identified_infection[primary_bacteria_idx],
            resistance_testing_done: individual.test_for_resistance[primary_bacteria_idx],
            resolution_type,
        }
    }
    
    fn determine_resolution_type(&self, individual: &Individual) -> String {
        if individual.date_of_death.is_some() {
            if let Some(ref cause) = individual.cause_of_death {
                match cause.as_str() {
                    "sepsis_related" => "DeathFromSepsis".to_string(),
                    "drug_toxicity_related" => "DeathFromToxicity".to_string(),
                    _ => "DeathFromBackground".to_string(),
                }
            } else {
                "DeathFromBackground".to_string()
            }
        } else {
            // Check if any drugs are active
            let has_active_drugs = individual.cur_use_drug.iter().any(|&taking| taking);
            if has_active_drugs {
                "DrugAssistedClearance".to_string()
            } else {
                "ImmuneClearance".to_string()
            }
        }
    }
    
    fn snapshot_to_csv(snapshot: &InfectionJourneySnapshot) -> String {
        // Format complex fields as semicolon-separated strings
        let all_bacteria_str = snapshot.all_bacteria_levels.iter()
            .map(|(name, level)| format!("{}:{:.6}", name, level))
            .collect::<Vec<_>>()
            .join(";");
        
        let current_drugs_str = snapshot.current_drugs.iter()
            .map(|(name, level)| format!("{}:{:.6}", name, level))
            .collect::<Vec<_>>()
            .join(";");
        
        let resistance_str = snapshot.resistance_any_r.iter()
            .map(|(drug, level)| format!("{}:{:.6}", drug, level))
            .collect::<Vec<_>>()
            .join(";");
        
        let mechanisms_str = snapshot.resistance_mechanisms.join(";");
        
        let resolution_str = snapshot.resolution_type.as_ref().unwrap_or(&String::new()).clone();
        
        format!(
            "{},{},{},{},{},{},{},{},{},{},{:.6},{},{},{},\"{}\",\"{}\",{},{},\"{}\",\"{}\",{},{:.6},{:.6},{:.6},{},{},{}",
            snapshot.journey_id,
            snapshot.individual_id,
            snapshot.time_step,
            snapshot.day_of_journey,
            snapshot.age_at_onset,
            snapshot.sex,
            snapshot.region_living,
            snapshot.region_current,
            snapshot.immunodeficiency,
            snapshot.primary_bacteria,
            snapshot.primary_bacteria_level,
            snapshot.syndrome,
            snapshot.sepsis,
            snapshot.hospital_acquired,
            all_bacteria_str,
            current_drugs_str,
            snapshot.days_on_current_treatment,
            snapshot.treatment_failures_count,
            resistance_str,
            mechanisms_str,
            snapshot.hospital_status,
            snapshot.immunity_level,
            snapshot.toxicity_level,
            snapshot.background_mortality_risk,
            snapshot.infection_identified,
            snapshot.resistance_testing_done,
            resolution_str
        )
    }
    
    pub fn get_stats(&self) -> (usize, u32, u32) {
        (self.active_journeys.len(), self.total_journeys_started, self.total_snapshots_logged)
    }

    pub fn finalize(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        if let Some(ref mut writer) = self.csv_writer {
            writer.flush()?;
        }
        
        // Journey logging summary disabled for cleaner output
        
        // Journey statistics summary disabled for cleaner output
        
        Ok(())
    }

    pub fn close(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        if let Some(writer) = self.csv_writer.take() {
            writer.into_inner()?.flush()?;
        }
        Ok(())
    }
}