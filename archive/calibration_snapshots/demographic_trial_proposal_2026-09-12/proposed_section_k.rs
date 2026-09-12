// PROPOSAL ONLY: not applied to src/config.rs.
        // === [K] Demographic distribution defaults ===
        // ---------------- 11) Demographic distribution defaults ----------------
        // There are 108 non-negative sampling weights: 6 regions x 18 age bands.
        // Population initialization normalizes them by their total at sampling time,
        // so they are relative weights and need not sum to 1.0.
        // Age-band bounds are days relative to birth; negative values are future cohorts.
        // Regional weighting pass, 2026-09-09: Asia x1.06, Africa x1.05, Europe x0.74,
        // North America x1.09, South America x1.07, Oceania x0.22.
        // Age-profile pass, 2026-09-09: eight birth cohorts per region fitted to
        // the under-80 age mix. Regional weight totals, the restored zero-ending
        // cohort, initially living cohorts and the unobserved future band are retained.

        // Demographic trial, 2026-09-12: reduce the 1930-1941 birth cohort
        // by 10-25% and redistribute within each region across the eight
        // younger birth cohorts. Regional totals are preserved. These are
        // provisional trial values; achieved age profiles require validation.

        // Asia demographic weights; raw weights emphasize future cohorts.
        map.insert("demo_asia_age_neg40000_neg36000".to_string(), 0.11448); // Heavy weighting on future births
        map.insert("demo_asia_age_neg36000_neg32000".to_string(), 0.075114366935);
        map.insert("demo_asia_age_neg32000_neg28000".to_string(), 0.086791298955);
        map.insert("demo_asia_age_neg28000_neg24000".to_string(), 0.08999365113);
        map.insert("demo_asia_age_neg24000_neg20000".to_string(), 0.081720298624);
        map.insert("demo_asia_age_neg20000_neg16000".to_string(), 0.073316415656);
        map.insert("demo_asia_age_neg16000_neg12000".to_string(), 0.069670286379);
        map.insert("demo_asia_age_neg12000_neg8000".to_string(), 0.05752127039);
        map.insert("demo_asia_age_neg8000_neg4000".to_string(), 0.043996411931);
        map.insert("demo_asia_age_neg4000_0".to_string(), 0.022896);          // Future births tapering
        // Positive-age bands continue the decline from the final future-cohort band.
        map.insert("demo_asia_age_0_4000".to_string(), 0.02332);             // Smooth continuation from neg4000_0
        map.insert("demo_asia_age_4000_8000".to_string(), 0.02014);
        map.insert("demo_asia_age_8000_12000".to_string(), 0.01696);
        map.insert("demo_asia_age_12000_16000".to_string(), 0.01378);
        map.insert("demo_asia_age_16000_20000".to_string(), 0.00954);
        map.insert("demo_asia_age_20000_24000".to_string(), 0.00636);
        map.insert("demo_asia_age_24000_28000".to_string(), 0.00318);
        map.insert("demo_asia_age_28000_32000".to_string(), 0.00106);

        // Africa demographic weights.
        map.insert("demo_africa_age_neg40000_neg36000".to_string(), 0.03465); // Heavy weighting on future births
        map.insert("demo_africa_age_neg36000_neg32000".to_string(), 0.046096263371);
        map.insert("demo_africa_age_neg32000_neg28000".to_string(), 0.040180545816);
        map.insert("demo_africa_age_neg28000_neg24000".to_string(), 0.027589382424);
        map.insert("demo_africa_age_neg24000_neg20000".to_string(), 0.023758334901);
        map.insert("demo_africa_age_neg20000_neg16000".to_string(), 0.018643663382);
        map.insert("demo_africa_age_neg16000_neg12000".to_string(), 0.009656608849);
        map.insert("demo_africa_age_neg12000_neg8000".to_string(), 0.008266245476);
        map.insert("demo_africa_age_neg8000_neg4000".to_string(), 0.009033955781);
        map.insert("demo_africa_age_neg4000_0".to_string(), 0.007875);          // Future births tapering
        // Positive-age bands continue the decline from the final future-cohort band.
        map.insert("demo_africa_age_0_4000".to_string(), 0.00945);             // Smooth continuation from neg4000_0
        map.insert("demo_africa_age_4000_8000".to_string(), 0.00735);
        map.insert("demo_africa_age_8000_12000".to_string(), 0.00525);
        map.insert("demo_africa_age_12000_16000".to_string(), 0.0042);
        map.insert("demo_africa_age_16000_20000".to_string(), 0.0021);
        map.insert("demo_africa_age_20000_24000".to_string(), 0.00105);
        map.insert("demo_africa_age_24000_28000".to_string(), 0.00105);
        map.insert("demo_africa_age_28000_32000".to_string(), 0.00105);

        // Europe demographic weights.
        map.insert("demo_europe_age_neg40000_neg36000".to_string(), 0.0148); // Moderate weighting on future births
        map.insert("demo_europe_age_neg36000_neg32000".to_string(), 0.007456960815);
        map.insert("demo_europe_age_neg32000_neg28000".to_string(), 0.008936564014);
        map.insert("demo_europe_age_neg28000_neg24000".to_string(), 0.010554267374);
        map.insert("demo_europe_age_neg24000_neg20000".to_string(), 0.010334669994);
        map.insert("demo_europe_age_neg20000_neg16000".to_string(), 0.010355560345);
        map.insert("demo_europe_age_neg16000_neg12000".to_string(), 0.012306071377);
        map.insert("demo_europe_age_neg12000_neg8000".to_string(), 0.011488649298);
        map.insert("demo_europe_age_neg8000_neg4000".to_string(), 0.009893256783);
        map.insert("demo_europe_age_neg4000_0".to_string(), 0.005994);          // Future births tapering
        map.insert("demo_europe_age_0_4000".to_string(), 0.00592);             // Larger portion alive in 1930
        map.insert("demo_europe_age_4000_8000".to_string(), 0.00518);
        map.insert("demo_europe_age_8000_12000".to_string(), 0.00444);
        map.insert("demo_europe_age_12000_16000".to_string(), 0.0037);
        map.insert("demo_europe_age_16000_20000".to_string(), 0.00296);
        map.insert("demo_europe_age_20000_24000".to_string(), 0.00222);
        map.insert("demo_europe_age_24000_28000".to_string(), 0.00148);
        map.insert("demo_europe_age_28000_32000".to_string(), 0.00148);

        // North America demographic weights.
        map.insert("demo_north_america_age_neg40000_neg36000".to_string(), 0.01308);
        map.insert("demo_north_america_age_neg36000_neg32000".to_string(), 0.008150724629);
        map.insert("demo_north_america_age_neg32000_neg28000".to_string(), 0.009156941472);
        map.insert("demo_north_america_age_neg28000_neg24000".to_string(), 0.00996263476);
        map.insert("demo_north_america_age_neg24000_neg20000".to_string(), 0.009158956125);
        map.insert("demo_north_america_age_neg20000_neg16000".to_string(), 0.008414060573);
        map.insert("demo_north_america_age_neg16000_neg12000".to_string(), 0.008599559899);
        map.insert("demo_north_america_age_neg12000_neg8000".to_string(), 0.007107454719);
        map.insert("demo_north_america_age_neg8000_neg4000".to_string(), 0.005340167823);
        map.insert("demo_north_america_age_neg4000_0".to_string(), 0.0027795);
        map.insert("demo_north_america_age_0_4000".to_string(), 0.00327);
        map.insert("demo_north_america_age_4000_8000".to_string(), 0.00327);
        map.insert("demo_north_america_age_8000_12000".to_string(), 0.00218);
        map.insert("demo_north_america_age_12000_16000".to_string(), 0.00218);
        map.insert("demo_north_america_age_16000_20000".to_string(), 0.00218);
        map.insert("demo_north_america_age_20000_24000".to_string(), 0.00218);
        map.insert("demo_north_america_age_24000_28000".to_string(), 0.00109);
        map.insert("demo_north_america_age_28000_32000".to_string(), 0.00109);

        // South America demographic weights.
        map.insert("demo_south_america_age_neg40000_neg36000".to_string(), 0.0107);
        map.insert("demo_south_america_age_neg36000_neg32000".to_string(), 0.006205468715);
        map.insert("demo_south_america_age_neg32000_neg28000".to_string(), 0.006971879198);
        map.insert("demo_south_america_age_neg28000_neg24000".to_string(), 0.008287331573);
        map.insert("demo_south_america_age_neg24000_neg20000".to_string(), 0.007330618268);
        map.insert("demo_south_america_age_neg20000_neg16000".to_string(), 0.006184320396);
        map.insert("demo_south_america_age_neg16000_neg12000".to_string(), 0.005354350968);
        map.insert("demo_south_america_age_neg12000_neg8000".to_string(), 0.005391215151);
        map.insert("demo_south_america_age_neg8000_neg4000".to_string(), 0.003815815731);
        map.insert("demo_south_america_age_neg4000_0".to_string(), 0.001819);
        map.insert("demo_south_america_age_0_4000".to_string(), 0.00214);
        map.insert("demo_south_america_age_4000_8000".to_string(), 0.00214);
        map.insert("demo_south_america_age_8000_12000".to_string(), 0.00107);
        map.insert("demo_south_america_age_12000_16000".to_string(), 0.00107);
        map.insert("demo_south_america_age_16000_20000".to_string(), 0.00107);
        map.insert("demo_south_america_age_20000_24000".to_string(), 0.00107);
        map.insert("demo_south_america_age_24000_28000".to_string(), 0.00107);
        map.insert("demo_south_america_age_28000_32000".to_string(), 0.00107);

        // Oceania demographic weights.
        map.insert("demo_oceania_age_neg40000_neg36000".to_string(), 0.00088);
        map.insert("demo_oceania_age_neg36000_neg32000".to_string(), 0.000696533346);
        map.insert("demo_oceania_age_neg32000_neg28000".to_string(), 0.000706705739);
        map.insert("demo_oceania_age_neg28000_neg24000".to_string(), 0.000659715012);
        map.insert("demo_oceania_age_neg24000_neg20000".to_string(), 0.000654535732);
        map.insert("demo_oceania_age_neg20000_neg16000".to_string(), 0.000670450779);
        map.insert("demo_oceania_age_neg16000_neg12000".to_string(), 0.000492262483);
        map.insert("demo_oceania_age_neg12000_neg8000".to_string(), 0.000500262533);
        map.insert("demo_oceania_age_neg8000_neg4000".to_string(), 0.000492534376);
        map.insert("demo_oceania_age_neg4000_0".to_string(), 0.000187);
        map.insert("demo_oceania_age_0_4000".to_string(), 0.00044);
        map.insert("demo_oceania_age_4000_8000".to_string(), 0.00044);
        map.insert("demo_oceania_age_8000_12000".to_string(), 0.00044);
        map.insert("demo_oceania_age_12000_16000".to_string(), 0.00022);
        map.insert("demo_oceania_age_16000_20000".to_string(), 0.00022);
        map.insert("demo_oceania_age_20000_24000".to_string(), 0.00022);
        map.insert("demo_oceania_age_24000_28000".to_string(), 0.00022);
        map.insert("demo_oceania_age_28000_32000".to_string(), 0.00022);

