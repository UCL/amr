import re
with open('src/rules/mod.rs', 'r') as f:
    c = f.read()

c = c.replace(
'''        EffluxAcrabTolc => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // All classical tetracyclines affected by RND efflux
           | "tigecycline"          // AcrAB-TolC overexpression (via ramA/marA) is the primary documented tigecycline resistance in Enterobacterales
           | "chloramphenicol" | "ciprofloxacin"
           | "erythromycin" | "azithromycin" | "clarithromycin" // Macrolides are classic RND efflux substrates
        ),''',
'''        EffluxAcrabTolc => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // All classical tetracyclines affected by RND efflux
           | "tigecycline"          // AcrAB-TolC overexpression (via ramA/marA) is the primary documented tigecycline resistance in Enterobacterales
           | "chloramphenicol" | "ciprofloxacin"
        ),'''
)

c = c.replace(
'''        EffluxMexxyOprm => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // Classical tetracyclines
           | "gentamicin" | "tobramycin" | "amikacin"            // Primary aminoglycoside efflux
           | "chloramphenicol" | "ciprofloxacin"
           | "erythromycin" | "azithromycin" | "clarithromycin" // Extrudes bulky macrolides
        ),''',
'''        EffluxMexxyOprm => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // Classical tetracyclines
           | "gentamicin" | "tobramycin" | "amikacin"            // Primary aminoglycoside efflux
           | "chloramphenicol" | "ciprofloxacin"
        ),'''
)

c = c.replace(
'''        GlobalEffluxPump => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // Classical tetracyclines
           | "tigecycline"          // Tigecycline evades tet-specific efflux but susceptible to broad RND pumps
           | "chloramphenicol" | "ciprofloxacin"
           | "erythromycin" | "azithromycin" | "clarithromycin" // Confers intrinsic macrolide resistance (e.g. in S. maltophilia)
        ),''',
'''        GlobalEffluxPump => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // Classical tetracyclines
           | "tigecycline"          // Tigecycline evades tet-specific efflux but susceptible to broad RND pumps
           | "chloramphenicol" | "ciprofloxacin"
        ),'''
)

with open('src/rules/mod.rs', 'w') as f:
    f.write(c)
