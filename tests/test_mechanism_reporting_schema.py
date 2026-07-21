import ast
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPER_TABLES_SOURCE = ROOT / "amr_simulation_output_analysis" / "make_paper_tables.py"
POPULATION_SOURCE = ROOT / "src" / "simulation" / "population.rs"


def _literal_assignment(name: str):
    tree = ast.parse(PAPER_TABLES_SOURCE.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.AnnAssign) and getattr(node.target, "id", None) == name:
            return ast.literal_eval(node.value)
    raise AssertionError(f"missing reporting assignment: {name}")


def _rust_mechanism_variants() -> set[str]:
    source = POPULATION_SOURCE.read_text(encoding="utf-8")
    enum_body = source.split("pub enum ResistanceMechanism {", 1)[1].split("}", 1)[0]
    return set(re.findall(r"^\s+([A-Z][A-Za-z0-9_]+),", enum_body, re.MULTILINE))


class MechanismReportingSchemaTests(unittest.TestCase):
    def test_exact_reporting_mechanisms_match_the_rust_enum(self) -> None:
        definitions = _literal_assignment("_SF4_EXACT_MECHANISMS")
        variants = [row["variant"].split("::")[-1] for row in definitions]
        slugs = [row["slug"] for row in definitions]

        self.assertEqual(set(variants), _rust_mechanism_variants())
        self.assertEqual(len(variants), len(set(variants)))
        self.assertEqual(len(slugs), len(set(slugs)))

    def test_reporting_families_partition_the_exact_mechanisms(self) -> None:
        definitions = _literal_assignment("_SF4_EXACT_MECHANISMS")
        families = _literal_assignment("_SF4_MECHANISM_FAMILIES")
        exact = {row["variant"] for row in definitions}
        family_variants = [variant for family in families for variant in family["variants"]]

        self.assertEqual(set(family_variants), exact)
        self.assertEqual(len(family_variants), len(set(family_variants)))


if __name__ == "__main__":
    unittest.main()
