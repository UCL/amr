use rand::SeedableRng;
use rand_chacha::ChaCha12Rng;

pub type ModelRng = ChaCha12Rng;

#[derive(Clone, Copy)]
pub enum RngStream {
    Initialization,
    JourneyLogger,
    RunId,
    TimestepChunk,
    TimestepMerge,
    ProfileRetention,
}

impl RngStream {
    const fn domain(self) -> u64 {
        match self {
            Self::Initialization => 0x18f4_d8a9_9d12_8b7d,
            Self::JourneyLogger => 0x7f39_2d4c_a91e_0865,
            Self::RunId => 0x4f6d_3c51_b572_22c7,
            Self::TimestepChunk => 0x9bc2_7475_3a8f_19e3,
            Self::TimestepMerge => 0xd6e8_58f9_2f0b_c44a,
            Self::ProfileRetention => 0x5b51_7a2e_c93d_f806,
        }
    }
}

pub fn model_rng(base_seed: u64, stream: RngStream, stream_id: u64) -> ModelRng {
    let mut state =
        base_seed ^ stream.domain().rotate_left(17) ^ stream_id.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let mut seed = [0u8; 32];
    for chunk in seed.chunks_exact_mut(8) {
        chunk.copy_from_slice(&splitmix64(&mut state).to_le_bytes());
    }
    ModelRng::from_seed(seed)
}

pub fn model_rng_from_entropy() -> ModelRng {
    ModelRng::from_entropy()
}

pub fn model_stream_seed(base_seed: u64, stream: RngStream, stream_id: u64) -> u64 {
    let mut rng = model_rng(base_seed, stream, stream_id);
    rand::Rng::gen(&mut rng)
}

pub fn timestep_stream_id(timestep: usize, chunk_index: usize) -> u64 {
    ((timestep as u64) << 32) ^ chunk_index as u64
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}
