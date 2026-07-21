use amr_project::config::PARAMETER_STORE;
use amr_project::rules::ParameterKeyCache;
use amr_project::simulation::population::{
    ResistanceMechanism, BACTERIA_LIST, DRUG_CLASS_LOOKUP, DRUG_SHORT_NAMES,
};
use csv::WriterBuilder;
use std::env;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

fn write_reachability_matrix<W: Write>(writer: W) -> csv::Result<()> {
    let cache = ParameterKeyCache::new();
    let mut csv = WriterBuilder::new()
        .terminator(csv::Terminator::Any(b'\n'))
        .from_writer(writer);
    csv.write_record([
        "bacteria",
        "drug",
        "resistance_representable",
        "maximum_any_r",
        "positive_effect_mechanisms",
    ])?;

    for (bacteria_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
        for (drug_idx, drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            let mechanisms = ResistanceMechanism::all()
                .iter()
                .enumerate()
                .filter_map(|(mechanism_idx, mechanism)| {
                    let effect = PARAMETER_STORE
                        .resistance_mechanism
                        .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_idx]);
                    (cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                        && effect > 0.0)
                        .then_some(mechanism.as_str())
                })
                .collect::<Vec<_>>();
            let maximum_any_r = ResistanceMechanism::all()
                .iter()
                .enumerate()
                .filter_map(|(mechanism_idx, _)| {
                    let effect = PARAMETER_STORE
                        .resistance_mechanism
                        .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_idx]);
                    (cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                        && effect > 0.0)
                        .then_some(effect)
                })
                .fold(0.0, |combined, effect| {
                    1.0 - (1.0 - combined) * (1.0 - effect)
                });
            csv.write_record([
                *bacteria,
                *drug,
                if mechanisms.is_empty() {
                    "false"
                } else {
                    "true"
                },
                &maximum_any_r.to_string(),
                &mechanisms.join(";"),
            ])?;
        }
    }
    csv.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if let Some(path) = env::args_os().nth(1) {
        let file = File::create(Path::new(&path))?;
        write_reachability_matrix(file)?;
    } else {
        write_reachability_matrix(io::stdout().lock())?;
    }
    Ok(())
}
