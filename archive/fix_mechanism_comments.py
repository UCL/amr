#!/usr/bin/env python3
"""
Fix mechanism emergence rate drug-class comments in src/config.rs.

Each change is annotated with ***changed: ... inline.

Rules derived by cross-referencing:
  - mechanism_applies_to_drug() in src/rules/mod.rs
  - mechanism_allowed_group_mask() in src/simulation/population.rs
"""

import re
import sys
from pathlib import Path

CONFIG_PATH = Path("src/config.rs")

with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    lines = f.readlines()

changed_count = 0

def fix_line(i, line):
    global changed_count
    orig = line

    # ---------------------------------------------------------------
    # 1. ESBL CTX-M / TEM / SHV
    #    Code: pen + flu + ceph(all) + aztreonam + aztreonam_avibactam
    #    Bug: comment has "bli" (wrong — BLIs inhibit ESBLs so BLI combos are
    #    NOT substrates); "flu" missing.
    # ---------------------------------------------------------------
    if re.search(r'_esbl_(ctx_m|tem|shv)_emergence_rate', line):
        line = re.sub(
            r'// classes: pen, bli, ceph, mono',
            '// classes: pen, flu, ceph, mono'
            ' ***changed: bli→removed (BLI combos inhibit ESBLs — not substrates);'
            ' flu→added (ESBLs hydrolyze flucloxacillin)',
            line
        )

    # ---------------------------------------------------------------
    # 2. AmpC CMY / DHA
    #    Code: pen + flu + bli + ceph(all) + ceftolozane_taz + aztreonam + az_avi
    #    Bug: "flu" missing.  "bli" is correct (AmpC is NOT inhibited by BLIs).
    # ---------------------------------------------------------------
    elif re.search(r'_ampc_(cmy|dha)_emergence_rate', line):
        line = re.sub(
            r'// classes: pen, bli, ceph, mono',
            '// classes: pen, flu, bli, ceph, mono'
            ' ***changed: flu→added (AmpC β-lactamases hydrolyze flucloxacillin'
            ' and are NOT inhibited by BLI combos)',
            line
        )

    # ---------------------------------------------------------------
    # 3. KPC
    #    Code: pen + flu + bli + ceph + carb + mono (all aztreonam variants)
    #    Bug: "flu" missing.
    # ---------------------------------------------------------------
    elif re.search(r'_enzyme_kpc_emergence_rate', line):
        line = re.sub(
            r'// classes: pen, bli, ceph, carb, mono',
            '// classes: pen, flu, bli, ceph, carb, mono'
            ' ***changed: flu→added (KPC hydrolyzes all penicillins incl. flucloxacillin)',
            line
        )

    # ---------------------------------------------------------------
    # 4. NDM / VIM (metallo-β-lactamases)
    #    Code: pen + flu + bli + ceph + carb + aztreonam_AVIBACTAM only
    #          (plain aztreonam NOT a substrate — MBLs cannot hydrolyze it)
    #    Bug: "flu" missing; "mono" misleading (plain aztreonam is NOT covered).
    # ---------------------------------------------------------------
    elif re.search(r'_ndm_vim_emergence_rate', line):
        line = re.sub(
            r'// classes: pen, bli, ceph, carb, mono',
            '// classes: pen, flu, bli, ceph, carb (not aztreonam; aztreonam_avibactam covered)'
            ' ***changed: flu→added; mono→clarified'
            ' (MBLs do NOT hydrolyze plain aztreonam; aztreonam_avibactam IS a substrate)',
            line
        )

    # ---------------------------------------------------------------
    # 5. OXA-48
    #    Code: pen + flu + bli + ceph + carb + ceftazidime_avibactam + aztreonam_avibactam
    #    Bug: "flu" missing.
    # ---------------------------------------------------------------
    elif re.search(r'_oxa_48_emergence_rate', line):
        line = re.sub(
            r'// classes: pen, bli, ceph, carb',
            '// classes: pen, flu, bli, ceph, carb'
            ' ***changed: flu→added (OXA-48 hydrolyzes all penicillins incl. flucloxacillin)',
            line
        )

    # ---------------------------------------------------------------
    # 6. PBP2a / mecA  (Gram-positive/Helicobacter group only)
    #
    #    6a. N. gonorrhoeae — Fastidious group, OUTSIDE PBP2a group mask →
    #        mechanism cannot fire; comment "// carb" is also wrong class.
    # ---------------------------------------------------------------
    if re.search(r'bacteria_neisseria_gonorrhoeae_mechanism_target_site_pbp2a_meca_emergence_rate', line):
        line = re.sub(
            r'// carb\s*$',
            '// tier 0 by group; N. gonorrhoeae is Fastidious, not GramPositive/Helicobacter;'
            ' PBP2a/mecA group mask excludes Fastidious'
            ' ***changed: carb→tier 0 by group (mechanism cannot fire for this organism;'
            ' and "carb" alone was incomplete — correct coverage is pen+flu+bli+ceph+carb+mono)',
            line
        )
    #    6b. All other bacteria with PBP2a — "pen, bli, ceph, carb" → add flu + mono
    elif re.search(r'_pbp2a_meca_emergence_rate', line):
        # Match the class list, preserve any trailing inline note after ';'
        line = re.sub(
            r'// classes: pen, bli, ceph, carb(.*?)$',
            lambda m: (
                '// classes: pen, flu, bli, ceph, carb, mono'
                + m.group(1)
                + ' ***changed: flu→added (mecA/PBP2a confers resistance to ALL'
                ' β-lactams incl. flucloxacillin — defining feature of MRSA);'
                ' mono→added (aztreonam also covered per code)'
            ),
            line
        )

    # ---------------------------------------------------------------
    # 7. ErmB
    #
    #    7a. P. aeruginosa — NonFermenter, OUTSIDE ErmB group mask →
    #        mechanism cannot fire.
    # ---------------------------------------------------------------
    if re.search(r'bacteria_pseudomonas_aeruginosa_mechanism_target_site_erm_b_emergence_rate', line):
        line = re.sub(
            r'// classes: mls',
            '// tier 0 by group; P. aeruginosa is NonFermenter;'
            ' ErmB group mask excludes NonFermenter'
            ' ***changed: mechanism cannot fire for this organism',
            line
        )
    #    7b. S. maltophilia — NonFermenter, OUTSIDE ErmB group mask → dead code
    elif re.search(r'bacteria_stenotrophomonas_maltophilia_mechanism_target_site_erm_b_emergence_rate', line):
        line = re.sub(
            r'// classes: mls;.*$',
            '// tier 0 by group; S. maltophilia is NonFermenter;'
            ' ErmB group mask excludes NonFermenter'
            ' ***changed: mechanism cannot fire for this organism'
            ' (non-zero rate is dead code)',
            line
        )

    # ---------------------------------------------------------------
    # 8. Cfr
    #    Code: linezolid + tedizolid (oxa) + chloramphenicol + clindamycin (lin)
    #          + retapamulin (pleuro) — does NOT cover macrolides or streptogramins
    #    Bug: "mls" implies full MLS (macrolides+lincosamides+streptogramins)
    #         but Cfr covers only lincosamides (not macrolides/streptogramins).
    #
    #    8a. P. aeruginosa — NonFermenter, OUTSIDE Cfr group mask → dead code
    # ---------------------------------------------------------------
    if re.search(r'bacteria_pseudomonas_aeruginosa_mechanism_target_site_cfr_emergence_rate', line):
        line = re.sub(
            r'// classes: oxa, mls, chl',
            '// tier 0 by group; P. aeruginosa is NonFermenter;'
            ' Cfr group mask excludes NonFermenter'
            ' ***changed: mechanism cannot fire for this organism'
            ' (and "mls" was also wrong: Cfr covers oxa+lin+chl+pleuro, not macrolides/streptogramins)',
            line
        )
    #    8b. All other Cfr entries with "oxa, mls, chl"
    elif re.search(r'_cfr_emergence_rate', line):
        # Preserve any trailing note after ';'
        line = re.sub(
            r'// classes: oxa, mls, chl(.*?)$',
            lambda m: (
                '// classes: oxa, lin, chl, pleuro'
                + m.group(1)
                + ' ***changed: mls→lin+pleuro'
                ' (Cfr covers oxazolidinones [oxa], lincosamides [lin = clindamycin],'
                ' chloramphenicol [chl], and pleuromutilins [pleuro = retapamulin];'
                ' does NOT cover macrolides or streptogramins)'
            ),
            line
        )

    # ---------------------------------------------------------------
    # 9. 16S rRNA methyltransferase (16sRrmt)
    #    Group mask: Enterobacterales, NonFermenter, EntericPathogen,
    #                Fastidious, Anaerobe — NOT GramPositive
    #
    #    9a. S. aureus — GramPositive → dead code (rate 0.3 is non-trivial!)
    # ---------------------------------------------------------------
    if re.search(r'bacteria_staphylococcus_aureus_mechanism_enzyme_16s_rrmt_emergence_rate', line):
        line = re.sub(
            r'// classes: ag',
            '// tier 0 by group; S. aureus is GramPositive;'
            ' 16sRrmt group mask excludes GramPositive'
            ' ***changed: aminoglycoside resistance in Staphylococci is modelled'
            ' via EnzymeAacAph, not 16sRrmt',
            line
        )
    #    9b. E. faecium — GramPositive → dead code
    elif re.search(r'bacteria_enterococcus_faecium_mechanism_enzyme_16s_rrmt_emergence_rate', line):
        line = re.sub(
            r'// classes: ag',
            '// tier 0 by group; E. faecium is GramPositive;'
            ' 16sRrmt group mask excludes GramPositive'
            ' ***changed: mechanism cannot fire for this organism',
            line
        )

    # ---------------------------------------------------------------
    # 10. GlobalEffluxPump
    #     Code: cipro + oflox + levo + moxi + tet + doxy + mino + tigecycline
    #            + chloramphenicol — NO macrolides/lincosamides/streptogramins
    #     Bug: "mls" in comment is wrong.
    # ---------------------------------------------------------------
    if re.search(r'_global_efflux_pump_emergence_rate', line):
        # Preserve any trailing note after the class list
        line = re.sub(
            r'// classes: fq, mls, tet, chl(.*?)$',
            lambda m: (
                '// classes: fq, tet, chl'
                + m.group(1)
                + ' ***changed: mls→removed'
                ' (GlobalEffluxPump does not cover macrolides/lincosamides/streptogramins;'
                ' covers FQ + all tetracyclines + chloramphenicol)'
            ),
            line
        )

    # ---------------------------------------------------------------
    # 11. EffluxMexxyOprm
    #     Code: tet + doxy + mino + ag + chl + cipro (and oflox, levo, moxi)
    #     Bug: several bacteria have "ceph, carb, fq, ag" which is completely wrong;
    #          ceph and carb are NOT covered; tet and chl ARE covered.
    # ---------------------------------------------------------------
    if re.search(r'_efflux_mexxy_oprm_emergence_rate', line):
        line = re.sub(
            r'// classes: ceph, carb, fq, ag',
            '// classes: fq, ag, tet, chl'
            ' ***changed: ceph+carb→removed; tet+chl→added'
            ' (MexXY-OprM covers FQ, aminoglycosides, tetracyclines, chloramphenicol;'
            ' does NOT cover cephalosporins or carbapenems)',
            line
        )

    # ---------------------------------------------------------------
    # 12. EffluxMtrCde
    #     Code: pen (NOT flu) + mac (erythro+azithro+clarithro, NOT clindamycin)
    #            + tet + doxy + mino + chl
    #     Group mask: Fastidious + EntericPathogen ONLY
    #
    #     12a. Enterobacterales bacteria → dead code by group mask
    #     (S. Typhi, S. Paratyphi A, iNTS, Shigella, Y. enterocolitica)
    # ---------------------------------------------------------------
    enterobacterales_mtr = [
        'salmonella_enterica_serovar_typhi',
        'salmonella_enterica_serovar_paratyphi_a',
        'invasive_non-typhoidal_salmonella_spp.',
        'shigella_spp.',
        'yersinia_enterocolitica',
    ]
    if any(re.search(r'bacteria_' + b + r'_mechanism_efflux_mtr_cde', line)
           for b in enterobacterales_mtr):
        line = re.sub(
            r'// classes: pen, mls, tet, chl.*$',
            '// tier 0 by group; this organism is Enterobacterales;'
            ' EffluxMtrCde group mask is Fastidious+EntericPathogen only'
            ' ***changed: mechanism cannot fire for this organism'
            ' (non-zero rate is dead code; and "mls" was also wrong: MtrCDE covers'
            ' pen+mac+tet+chl, not full MLS)',
            line
        )
    #     12b. C. jejuni (Helicobacter group) → dead code by group mask
    elif re.search(r'bacteria_campylobacter_jejuni_mechanism_efflux_mtr_cde_emergence_rate', line):
        line = re.sub(
            r'// classes: pen, mls, tet, chl.*$',
            '// tier 0 by group; C. jejuni is Helicobacter;'
            ' EffluxMtrCde group mask excludes Helicobacter'
            ' ***changed: CmeABC efflux has no model equivalent here; mechanism cannot fire'
            ' (and "mls" was also wrong: MtrCDE covers pen+mac+tet+chl, not full MLS)',
            line
        )
    #     12c. Active bacteria (Fastidious / EntericPathogen):
    #          "mls" → "mac" (MtrCDE covers macrolides only, not clindamycin/streptogramins)
    elif re.search(r'_efflux_mtr_cde_emergence_rate', line):
        line = re.sub(
            r'// classes: pen, mls, tet, chl(.*?)$',
            lambda m: (
                '// classes: pen, mac, tet, chl'
                + m.group(1)
                + ' ***changed: mls→mac'
                ' (MtrCDE efflux covers macrolides [erythro/azithro/clarithro],'
                ' NOT clindamycin or streptogramins)'
            ),
            line
        )

    # ---------------------------------------------------------------
    # 13. Mutation23sRrna
    #     Code: erythromycin + azithromycin + clarithromycin ONLY
    #           (NOT clindamycin, NOT quinupristin-dalfopristin)
    #     Bug: comment says "mls" which implies full MLS.
    #     Note: some lines also have "disabled" — still correct the class label.
    # ---------------------------------------------------------------
    if re.search(r'_mutation_23s_rrna_emergence_rate', line):
        # Careful: only fix the "// classes: mls" pattern (not "// classes: mls; disabled...")
        # Replace "classes: mls" (with optional trailing note) → "classes: mac"
        line = re.sub(
            r'// classes: mls(.*?)$',
            lambda m: (
                '// classes: mac (erythro, azithro, clarithro only; not clindamycin)'
                + m.group(1)
                + ' ***changed: mls→mac'
                ' (23S rRNA point mutations affect macrolides only,'
                ' NOT lincosamides or streptogramins)'
            ),
            line
        )

    # ---------------------------------------------------------------
    # 14. MutationPbpMosaic
    #     Code: pen + flu + bli + ceph(all incl. ceftaroline+ceftolozane_taz+ceftazidime_avi)
    #            + aztreonam + aztreonam_avibactam — NOT carbapenems
    #     Bug: "flu" missing in "pen, bli, ceph, mono".
    # ---------------------------------------------------------------
    if re.search(r'_mutation_pbp_mosaic_emergence_rate', line):
        line = re.sub(
            r'// classes: pen, bli, ceph, mono(.*?)$',
            lambda m: (
                '// classes: pen, flu, bli, ceph, mono'
                + m.group(1)
                + ' ***changed: flu→added'
                ' (PBP mosaic mutations affect all penicillins incl. flucloxacillin,'
                ' but NOT carbapenems)'
            ),
            line
        )

    # ---------------------------------------------------------------
    # 15. EnzymeAacAph dead-code cases
    #     Group mask: Enterobacterales + NonFermenter + EntericPathogen + Fastidious
    #                 + GramPositive — NOT Helicobacter, NOT Anaerobe
    #
    #     15a. C. jejuni (Helicobacter) → dead code
    # ---------------------------------------------------------------
    if re.search(r'bacteria_campylobacter_jejuni_mechanism_enzyme_aac_aph_emergence_rate', line):
        line = re.sub(
            r'// classes: ag.*$',
            '// tier 0 by group; C. jejuni is Helicobacter;'
            ' EnzymeAacAph group mask excludes Helicobacter'
            ' ***changed: mechanism cannot fire for this organism'
            ' (aminoglycoside resistance in Campylobacter occurs via other routes)',
            line
        )
    #     15b. B. fragilis (Anaerobe) → dead code
    elif re.search(r'bacteria_bacteroides_fragilis_mechanism_enzyme_aac_aph_emergence_rate', line):
        line = re.sub(
            r'// classes: ag.*$',
            '// tier 0 by group; B. fragilis is Anaerobe;'
            ' EnzymeAacAph group mask excludes Anaerobe'
            ' ***changed: mechanism cannot fire for this organism'
            ' (aminoglycosides are intrinsically inactive against anaerobes)',
            line
        )

    if line != orig:
        changed_count += 1

    return line


result_lines = [fix_line(i, line) for i, line in enumerate(lines)]

with open(CONFIG_PATH, "w", encoding="utf-8") as f:
    f.writelines(result_lines)

print(f"Done. {changed_count} lines modified.")
