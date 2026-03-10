import sys

with open('src/rules/mod.rs', 'r', encoding='utf-8') as f:
    content = f.read()

trait_code = """
use std::f32::consts::LN_2;

pub trait FastMath {
    fn fast_exp(self) -> Self;
    fn fast_ln(self) -> Self;
}

impl FastMath for f64 {
    #[inline(always)]
    fn fast_exp(self) -> Self {
        fast_math::exp(self as f32) as f64
    }

    #[inline(always)]
    fn fast_ln(self) -> Self {
        (fast_math::log2(self as f32) * LN_2) as f64
    }
}
"""

# Insert trait near the top after the first few imports/structs
insert_idx = content.find("pub struct ApplyRulesConfig")
content = content[:insert_idx] + trait_code + "\n" + content[insert_idx:]

content = content.replace(".exp()", ".fast_exp()")
content = content.replace(".ln()", ".fast_ln()")

with open('src/rules/mod.rs', 'w', encoding='utf-8') as f:
    f.write(content)
