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
    pub resistance_majority_r: Vec<(String, f64)>, // (drug_name, majority_r_level)
    pub resistance_activity_r: Vec<(String, f64)>, // (drug_name, activity_r_level)
    pub resistance_mechanisms: Vec<String>, // active mechanisms
    
    // Drug selection information (captured when treatment starts)
    pub drug_selection_bacteria: Option<String>, // bacteria that triggered drug selection
    pub drug_selection_scores: Vec<(String, f64)>, // all drug scores at selection time
    pub selected_drug: Option<String>, // which drug was actually selected
    
    // Clinical status
    pub hospital_status: String,
    pub immunity_level: f64,
    pub toxicity_level: f64,
    pub background_mortality_risk: f64,
    
    // Testing status
    pub infection_identified: bool,
    pub infection_has_caused_symptoms: bool,
    pub resistance_testing_done: bool,
    
    // Journey outcome (only on final snapshot)
    pub resolution_type: Option<String>,
    
    // De novo resistance detection
    pub has_de_novo_resistance: bool,
    
    // Resistance source tracking
    pub resistance_sources: Vec<(String, String)>, // (drug_name, acquisition_source)
}

#[derive(Clone, Debug)]
pub struct ActiveJourney {
    pub journey_id: u32,
    pub individual_id: usize,
    pub onset_time_step: usize,
    pub primary_bacteria_idx: usize,
    pub day_count: u32,
    pub snapshots: Vec<InfectionJourneySnapshot>,
    pub primary_bacteria_cleared_day: Option<u32>, // Day when primary bacteria cleared
    pub has_de_novo_resistance: bool, // Track if de novo resistance emerged during this journey
}

pub struct JourneyLogger {
    // Configuration
    pub enabled: bool,
    pub sample_rate: f64,
    pub bacteria_filter: Option<String>, // Filter to log only specific bacteria
    
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
            bacteria_filter: None,
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
        self.bacteria_filter = None;
        
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
    
    /// Enable journey logging with sample rate and optional bacteria filter
    pub fn enable_with_filter(&mut self, sample_rate: f64, bacteria_filter: Option<String>) -> Result<(), Box<dyn std::error::Error>> {
        self.enabled = true;
        self.sample_rate = sample_rate.clamp(0.0, 1.0);
        self.bacteria_filter = bacteria_filter;
        
        // Create output file with bacteria-specific name if filtering
        self.output_filename = if let Some(ref filter) = self.bacteria_filter {
            format!("infection_journeys_{}.csv", filter)
        } else {
            "infection_journeys.csv".to_string()
        };
        
        let file = File::create(&self.output_filename)?;
        let mut writer = BufWriter::new(file);
        
        // Write CSV header
        writeln!(writer, "{}", JourneyLogger::get_csv_header())?;
        
        self.csv_writer = Some(writer);
        
        // Journey logging enabled silently
        
        Ok(())
    }
    
    fn get_csv_header() -> &'static str {
        "journey_id,individual_id,time_step,day_of_journey,age_at_onset,sex,region_living,region_current,immunodeficiency,primary_bacteria,primary_bacteria_level,syndrome,sepsis,hospital_acquired,all_bacteria_levels,current_drugs,days_on_current_treatment,treatment_failures,resistance_any_r,resistance_majority_r,resistance_activity_r,resistance_mechanisms,drug_selection_bacteria,drug_selection_scores,selected_drug,hospital_status,immunity_level,toxicity_level,background_mortality_risk,infection_identified,infection_has_caused_symptoms,resistance_testing_done,resolution_type,has_de_novo_resistance,resistance_sources"
    }
    
    pub fn check_individual(&mut self, individual: &Individual, time_step: usize) {
        if !self.enabled {
            return;
        }
        
        let individual_id = individual.id;
        let has_active_infection = individual.level.iter().any(|&level| level > 0.001);
        let is_currently_tracked = self.active_journeys.contains_key(&individual_id);
        
        // Check if individual is dead
        let is_dead = individual.date_of_death.is_some();
        
        match (is_currently_tracked, has_active_infection, is_dead) {
            (false, true, false) => {
                // New infection - potentially start tracking
                let has_de_novo_resistance = self.detect_de_novo_resistance_emergence(individual);
                let mut rng = rand::thread_rng();
                
                // Always sample if de novo resistance is detected, otherwise use normal sampling
                if has_de_novo_resistance || rng.gen::<f64>() < self.sample_rate {
                    self.start_journey(individual, time_step);
                }
            },
            (true, _, true) => {
                // Individual died - complete journey immediately
                self.complete_journey(individual, time_step);
            },
            (true, _, false) => {
                // Continue tracking - check if we should terminate based on clearance + time
                if self.should_terminate_journey(individual, time_step) {
                    self.complete_journey(individual, time_step);
                } else {
                    self.update_journey(individual, time_step);
                }
            },
            (false, false, false) | (false, _, true) => {
                // Not tracking - do nothing (covers dead/alive, no infection cases)
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
        
        // Check bacteria filter if specified
        if let Some(ref filter_bacteria) = self.bacteria_filter {
            let bacteria_name = BACTERIA_LIST[primary_bacteria_idx].to_lowercase().replace(" ", "_");
            if bacteria_name != *filter_bacteria {
                return; // Skip this infection - doesn't match filter
            }
        }
        
        let journey_id = self.next_journey_id;
        self.next_journey_id += 1;
        
        // Create initial snapshot
        let snapshot = JourneyLogger::create_snapshot(individual, journey_id, time_step, 1, primary_bacteria_idx, None, false);
        
        let journey = ActiveJourney {
            journey_id,
            individual_id: individual.id,
            onset_time_step: time_step,
            primary_bacteria_idx,
            day_count: 1,
            snapshots: vec![snapshot],
            primary_bacteria_cleared_day: None,
            has_de_novo_resistance: false,
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
        // Check for de novo resistance emergence before getting mutable reference
        let has_new_resistance = self.detect_de_novo_resistance_emergence(individual);
        
        if let Some(journey) = self.active_journeys.get_mut(&individual.id) {
            journey.day_count += 1;
            
            // Update resistance flag if detected
            if !journey.has_de_novo_resistance && has_new_resistance {
                journey.has_de_novo_resistance = true;
            }
            
            let snapshot = JourneyLogger::create_snapshot(
                individual, 
                journey.journey_id, 
                time_step, 
                journey.day_count, 
                journey.primary_bacteria_idx,
                None,
                journey.has_de_novo_resistance
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
    
    fn should_terminate_journey(&mut self, individual: &Individual, _time_step: usize) -> bool {
        if let Some(journey) = self.active_journeys.get_mut(&individual.id) {
            let primary_bacteria_level = individual.level[journey.primary_bacteria_idx];
            
            // Check if primary bacteria has cleared
            if primary_bacteria_level <= 0.001 {
                if journey.primary_bacteria_cleared_day.is_none() {
                    // First time we detected clearance - record the day
                    journey.primary_bacteria_cleared_day = Some(journey.day_count);
                }
                
                // Check if 7 days have passed since clearance
                if let Some(cleared_day) = journey.primary_bacteria_cleared_day {
                    return journey.day_count >= cleared_day + 7;
                }
            } else {
                // Bacteria level increased again - reset clearance tracking
                journey.primary_bacteria_cleared_day = None;
            }
        }
        false
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
                Some(resolution_type.clone()),
                journey.has_de_novo_resistance
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
        resolution_type: Option<String>,
        has_de_novo_resistance: bool
    ) -> InfectionJourneySnapshot {
        
        // Collect all active bacteria
        let all_bacteria_levels: Vec<(String, f64)> = individual.level.iter()
            .enumerate()
            .filter(|(_, &level)| level > 0.001)
            .map(|(idx, &level)| (BACTERIA_LIST[idx].to_string(), level))
            .collect();
        
        // Collect ALL drugs with detectable levels (active AND decaying after cessation)
        let mut current_drugs: Vec<(String, f64)> = Vec::new();
        
        // Include all drugs with levels above detection threshold (0.001)
        // This captures both actively taken drugs and drugs in decay phase after cessation
        for (idx, &level) in individual.cur_level_drug.iter().enumerate() {
            if level > 0.001 {
                let drug_name = DRUG_SHORT_NAMES[idx].to_string();
                current_drugs.push((drug_name, level));
            }
        }
        
        // Also show drugs that were just initiated this time step
        for (idx, &initiation_day) in individual.date_drug_initiated.iter().enumerate() {
            if initiation_day == time_step as i32 {
                let drug_name = DRUG_SHORT_NAMES[idx].to_string();
                // Check if not already added from cur_use_drug
                let already_added = current_drugs.iter().any(|(name, _)| name == &drug_name);
                if !already_added {
                    // Use initial level for newly initiated drugs
                    let initial_level = if individual.cur_level_drug[idx] > 0.0 {
                        individual.cur_level_drug[idx]
                    } else {
                        10.0 // Standard initial level
                    };
                    current_drugs.push((drug_name, initial_level));
                }
            }
        }
        
        // Also include drugs that were initiated in the last few days even if no longer active
        // This provides comprehensive drug history for better timing accuracy
        for (idx, &initiation_day) in individual.date_drug_initiated.iter().enumerate() {
            if initiation_day >= 0 && initiation_day < time_step as i32 - 1 {  // Skip current and previous day (already handled above)
                let drug_name = DRUG_SHORT_NAMES[idx].to_string();
                // Only add if not already in current_drugs
                if !current_drugs.iter().any(|(name, _)| name == &drug_name) {
                    let days_since_initiation = time_step as i32 - initiation_day;
                    if days_since_initiation <= 2 {  // Show drugs used within last 2 days
                        // Show recently used drugs with minimal level to indicate they were used
                        current_drugs.push((drug_name, 0.1));
                    }
                }
            }
        }
        
        // Collect resistance for primary bacteria
        let resistance_any_r: Vec<(String, f64)> = DRUG_SHORT_NAMES.iter()
            .enumerate()
            .filter(|(idx, _)| individual.resistances[primary_bacteria_idx][*idx].any_r > 0.0)
            .map(|(idx, &drug_name)| (drug_name.to_string(), individual.resistances[primary_bacteria_idx][idx].any_r))
            .collect();
        
        let resistance_majority_r: Vec<(String, f64)> = DRUG_SHORT_NAMES.iter()
            .enumerate()
            .filter(|(idx, _)| individual.resistances[primary_bacteria_idx][*idx].majority_r > 0.0)
            .map(|(idx, &drug_name)| (drug_name.to_string(), individual.resistances[primary_bacteria_idx][idx].majority_r))
            .collect();
        
        let resistance_activity_r: Vec<(String, f64)> = DRUG_SHORT_NAMES.iter()
            .enumerate()
            .filter(|(idx, _)| {
                // Only log activity_r when both drug and bacteria are present
                // This covers cases where activity_r = 0.0 (complete resistance) or > 0.0 (partial effectiveness)
                individual.cur_use_drug[*idx] && // Drug is being used
                individual.level[primary_bacteria_idx] > 0.001 // AND bacteria are present
            })
            .map(|(idx, &drug_name)| (drug_name.to_string(), individual.resistances[primary_bacteria_idx][idx].activity_r))
            .collect();
        
        // Collect active resistance mechanisms
        let resistance_mechanisms: Vec<String> = individual.resistance_mechanisms[primary_bacteria_idx].iter()
            .enumerate()
            .filter(|(_, &is_active)| is_active)
            .map(|(idx, _)| format!("mechanism_{}", idx)) // Will need to map to actual mechanism names
            .collect();
        
        // Collect resistance sources for primary bacteria
        let resistance_sources: Vec<(String, String)> = DRUG_SHORT_NAMES.iter()
            .enumerate()
            .filter_map(|(idx, &drug_name)| {
                if let Some(acquisition_type) = &individual.how_resistance_acquired[primary_bacteria_idx][idx] {
                    Some((drug_name.to_string(), acquisition_type.as_str().to_string()))
                } else {
                    None
                }
            })
            .collect();
        
        // Collect drug selection information if available
        let (drug_selection_bacteria, drug_selection_scores, selected_drug) = if individual.bacteria_on_selection_day >= 0 {
            let selection_bacteria_idx = individual.bacteria_on_selection_day as usize;
            let selection_bacteria_name = if selection_bacteria_idx < BACTERIA_LIST.len() {
                Some(BACTERIA_LIST[selection_bacteria_idx].to_string())
            } else {
                None
            };
            
            let scores: Vec<(String, f64)> = DRUG_SHORT_NAMES.iter()
                .enumerate()
                .filter(|(idx, _)| individual.drug_score_on_selection_day[*idx] >= 0.0)
                .map(|(idx, &drug_name)| (drug_name.to_string(), individual.drug_score_on_selection_day[idx]))
                .collect();
            
            // Find which drug was selected (the one currently being taken that was started most recently)
            let selected = current_drugs.iter()
                .max_by_key(|(drug_name, _)| {
                    if let Some(drug_idx) = DRUG_SHORT_NAMES.iter().position(|&d| d == drug_name) {
                        individual.date_drug_initiated[drug_idx]
                    } else {
                        -1
                    }
                })
                .map(|(drug_name, _)| drug_name.clone());
            
            (selection_bacteria_name, scores, selected)
        } else {
            (None, Vec::new(), None)
        };
        
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
            resistance_majority_r,
            resistance_activity_r,
            resistance_mechanisms,
            drug_selection_bacteria,
            drug_selection_scores,
            selected_drug,
            hospital_status: format!("{:?}", individual.hospital_status),
            immunity_level: individual.immune_resp[primary_bacteria_idx],
            toxicity_level: individual.current_toxicity,
            background_mortality_risk: individual.background_all_cause_mortality_rate,
            infection_identified: individual.test_identified_infection[primary_bacteria_idx],
            infection_has_caused_symptoms: individual.infection_has_caused_symptoms[primary_bacteria_idx],
            resistance_testing_done: individual.test_for_resistance[primary_bacteria_idx],
            resolution_type,
            has_de_novo_resistance,
            resistance_sources,
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
        
        let resistance_any_r_str = snapshot.resistance_any_r.iter()
            .map(|(drug, level)| format!("{}:{:.6}", drug, level))
            .collect::<Vec<_>>()
            .join(";");
        
        let resistance_majority_r_str = snapshot.resistance_majority_r.iter()
            .map(|(drug, level)| format!("{}:{:.6}", drug, level))
            .collect::<Vec<_>>()
            .join(";");
        
        let resistance_activity_r_str = snapshot.resistance_activity_r.iter()
            .map(|(drug, level)| format!("{}:{:.6}", drug, level))
            .collect::<Vec<_>>()
            .join(";");
        
        let mechanisms_str = snapshot.resistance_mechanisms.join(";");
        
        let drug_selection_bacteria_str = snapshot.drug_selection_bacteria.as_ref().unwrap_or(&String::new()).clone();
        
        let drug_selection_scores_str = snapshot.drug_selection_scores.iter()
            .map(|(drug, score)| format!("{}:{:.6}", drug, score))
            .collect::<Vec<_>>()
            .join(";");
        
        let selected_drug_str = snapshot.selected_drug.as_ref().unwrap_or(&String::new()).clone();
        
        let resolution_str = snapshot.resolution_type.as_ref().unwrap_or(&String::new()).clone();
        
        let resistance_sources_str = snapshot.resistance_sources.iter()
            .map(|(drug, source)| format!("{}:{}", drug, source))
            .collect::<Vec<_>>()
            .join(";");
        
        format!(
            "{},{},{},{},{},{},{},{},{},{},{:.6},{},{},{},\"{}\",\"{}\",{},{},\"{}\",\"{}\",\"{}\",\"{}\",{},\"{}\",{},{},{:.6},{:.6},{:.6},{},{},{},{},{},\"{}\"",
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
            resistance_any_r_str,
            resistance_majority_r_str,
            resistance_activity_r_str,
            mechanisms_str,
            drug_selection_bacteria_str,
            drug_selection_scores_str,
            selected_drug_str,
            snapshot.hospital_status,
            snapshot.immunity_level,
            snapshot.toxicity_level,
            snapshot.background_mortality_risk,
            snapshot.infection_identified,
            snapshot.infection_has_caused_symptoms,
            snapshot.resistance_testing_done,
            resolution_str,
            snapshot.has_de_novo_resistance,
            resistance_sources_str
        )
    }
    
    pub fn get_stats(&self) -> (usize, u32, u32) {
        (self.active_journeys.len(), self.total_journeys_started, self.total_snapshots_logged)
    }

    // Detect if de novo resistance has emerged during active treatment
    fn detect_de_novo_resistance_emergence(&self, individual: &Individual) -> bool {
        // Check if any resistance acquisition events occurred while on active treatment
        // Look for resistance acquisition types that indicate de novo emergence:
        // - AtInfectionEnv: acquired during infection in environment
        // - AtInfectionTB: acquired during infection in tissue/blood  
        // - FromMicrobiomeR: transferred from resistant microbiome
        // - Hgt: horizontal gene transfer
        
        // Check if individual is currently on drugs
        let on_active_treatment = individual.cur_use_drug.iter().any(|&taking| taking);
        
        if !on_active_treatment {
            return false;
        }
        
        // Check how_resistance_acquired field for de novo patterns
        // This is a 2D Vec indexed by [bacteria][drug]
        for (_bacteria_idx, bacteria_resistances) in individual.how_resistance_acquired.iter().enumerate() {
            for (_drug_idx, acquisition_type_opt) in bacteria_resistances.iter().enumerate() {
                if let Some(acquisition_type) = acquisition_type_opt {
                    // Check if this resistance was acquired during treatment
                    match acquisition_type {
                        crate::simulation::population::ResistanceAcquisitionType::AtInfectionEnv |
                        crate::simulation::population::ResistanceAcquisitionType::AtInfectionTB |
                        crate::simulation::population::ResistanceAcquisitionType::FromMicrobiomeR |
                        crate::simulation::population::ResistanceAcquisitionType::Hgt => {
                            return true;
                        },
                        _ => {
                            // Other types like AtInfectionCommunity are less concerning
                            // as they represent pre-existing resistance
                        }
                    }
                }
            }
        }
        
        false
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