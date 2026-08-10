use amr_project::config::PARAMETER_STORE;
use amr_project::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES};
use csv::WriterBuilder;
use std::env;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

fn write_potency_matrix<W: Write>(writer: W) -> csv::Result<()> {
    let mut csv = WriterBuilder::new()
        .terminator(csv::Terminator::Any(b'\n'))
        .from_writer(writer);
    csv.write_record(["bacteria", "drug", "potency_when_no_r"])?;

    for (bacteria_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
        for (drug_idx, drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            let potency = PARAMETER_STORE
                .drug_bacteria
                .potency(bacteria_idx, drug_idx)
                .to_string();
            csv.write_record([*bacteria, *drug, potency.as_str()])?;
        }
    }
    csv.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if let Some(path) = env::args_os().nth(1) {
        let file = File::create(Path::new(&path))?;
        write_potency_matrix(file)?;
    } else {
        write_potency_matrix(io::stdout().lock())?;
    }
    Ok(())
}
