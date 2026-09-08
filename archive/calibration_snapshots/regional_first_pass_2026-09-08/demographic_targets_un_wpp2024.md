# UN population age targets

`demographic_targets_un_wpp2024.csv` contains 30 regional age shares derived from
the UN Population Division's **World Population Prospects 2024**, both sexes,
**1 July 2023**, `Estimates` worksheet. The year is an estimate year, not a 2025
projection. Retrieval date: 8 September 2026.

Official sources:

- [WPP2024 population by single age, both sexes workbook](https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20%28Standard%29/EXCEL_FILES/2_Population/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx)
- [UN download catalogue](https://population.un.org/wpp/assets/downloads.json)
- [WPP2024 overview](https://population.un.org/wpp/)
- [UN definition of regions](https://population.un.org/wpp/definition-of-regions)

Citation: United Nations, Department of Economic and Social Affairs, Population
Division (2024). *World Population Prospects 2024, Online Edition.* The workbook
identifies its licence as CC BY 3.0 IGO. These extracts and aggregations are derived
data; the UN has not endorsed this project's application of them.

The workbook has one column per completed year of age 0 through 99 and a final
100+ column. Exact aggregation is:

| Model age_group | Source columns summed | Age interval |
| --- | --- | --- |
| 0_5 | 0 through 5 | birth to below age 6 |
| 6_14 | 6 through 14 | age 6 to below age 15 |
| 15_49 | 15 through 49 | age 15 to below age 50 |
| 50_79 | 50 through 79 | age 50 to below age 80 |
| 80plus | 80 through 99 and 100+ | age 80 and above |

For each region, `age_share` is the sum for that band divided by the sum across
all 101 age columns. No within-band uniformity or interpolation is assumed.
The source unit is thousands of persons. `age_population` and
`regional_total_population` multiply those counts by 1,000. Fractional persons
reflect the source's statistical estimates; they are not a census enumeration.
The total is derived from the complete single-age distribution. Shares are
written to 15 decimal places and sum to one within 3e-15.

| Model region | UN source region | UN location code | Workbook row | 2023 population, persons |
| --- | --- | --- | --- | ---: |
| north_america | Northern America | 905 | 19558 | 382,902,741.5 |
| south_america | South America | 931 | 18448 | 433,024,175.0 |
| europe | Europe | 908 | 11491 | 745,602,874.5 |
| asia | Asia | 935 | 7273 | 4,778,004,485.5 |
| africa | Africa | 903 | 2537 | 1,480,770,525.0 |
| oceania | Oceania | 909 | 20002 | 45,562,787.0 |

**Geographic caveat:** UN Northern America excludes Central America and the
Caribbean; South America also excludes those areas. These six source aggregates
therefore do not partition the world population. The model's `north_america`
label is explicitly mapped to UN Northern America for this extract. If the
model intends a different continental boundary, its target should be recomputed
from constituent countries. No region sampling shares are specified in this
file; the separate calibration instructions determine those shares.

`demographic_source_single_ages_un_wpp2024.csv` preserves the 606 source values
in the original thousands-of-persons units, with workbook row identifiers.
`demographic_source_un_wpp2024.manifest.json` records source URL, byte size and
SHA-256. `extract_un_demographic_targets.py` reproduces both CSVs from the
downloaded workbook using Python's standard library:

```text
python extract_un_demographic_targets.py /path/to/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx
```

The 220 MB workbook is omitted from the repository. The extractor verifies all
six official region names/codes, the estimate year and variant, complete age
coverage, nonnegative counts, and normalized age shares.
