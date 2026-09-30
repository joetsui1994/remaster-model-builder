#!/usr/bin/env python3
"""Tests for spatial CSV validation, rate conventions, generation, and audits."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from generate_spatial_reactions import audit_mobility, generate_bdm, generate_coalescent, parse_demes, parse_mobility


class CSVTestCase(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def csv(self, name: str, content: str) -> Path:
        path = self.path / name
        path.write_text(content.strip() + "\n", encoding="utf-8")
        return path


class TestCSVValidation(CSVTestCase):
    def test_valid_bdm_and_mobility(self):
        demes = parse_demes(self.csv("d.csv", "deme,S,I,R,beta,gamma\nA,1000,1,0,0.001,0.1\nB,2000,0,0,0.001,0.1"))
        mobility = parse_mobility(self.csv("m.csv", "origin,destination,rate\nA,B,0.01\nB,A,0.01"), {d["deme"] for d in demes})
        self.assertEqual((len(demes), len(mobility)), (2, 2))

    def test_bdm_requires_all_model_columns(self):
        with self.assertRaisesRegex(ValueError, "missing required column"):
            parse_demes(self.csv("d.csv", "deme,S,I,R\nA,1000,1,0"))

    def test_coalescent_requires_population_size(self):
        with self.assertRaisesRegex(ValueError, "popSize.*N"):
            parse_demes(self.csv("d.csv", "deme\nA"), mode="coalescent")

    def test_duplicate_and_unsafe_demes_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate deme"):
            parse_demes(self.csv("d.csv", "deme,S,I,R,beta,gamma\nA,1,0,0,1,1\nA,1,0,0,1,1"))
        with self.assertRaisesRegex(ValueError, "safe identifier"):
            parse_demes(self.csv("d.csv", "deme,S,I,R,beta,gamma\nA&amp;B,1,0,0,1,1"))

    def test_nonfinite_and_fractional_values_rejected(self):
        with self.assertRaisesRegex(ValueError, "finite"):
            parse_demes(self.csv("d.csv", "deme,S,I,R,beta,gamma\nA,1,0,0,nan,1"))
        with self.assertRaisesRegex(ValueError, "integer count"):
            parse_demes(self.csv("d.csv", "deme,S,I,R,beta,gamma\nA,1.5,0,0,1,1"))
        with self.assertRaisesRegex(ValueError, "finite"):
            parse_mobility(self.csv("m.csv", "origin,destination,rate\nA,B,inf"), {"A", "B"})

    def test_mobility_referential_integrity_and_duplicates(self):
        with self.assertRaisesRegex(ValueError, "unknown origin"):
            parse_mobility(self.csv("m.csv", "origin,destination,rate\nAtlantis,B,0.1"), {"A", "B"})
        with self.assertRaisesRegex(ValueError, "unknown destination"):
            parse_mobility(self.csv("m.csv", "origin,destination,rate\nA,Atlantis,0.1"), {"A", "B"})
        with self.assertRaisesRegex(ValueError, "duplicate directed"):
            parse_mobility(self.csv("m.csv", "origin,destination,rate\nA,B,0.1\nA,B,0.2"), {"A", "B"})


class TestMobilityAudits(unittest.TestCase):
    def test_clean_balanced_bdm_network(self):
        demes = [{"deme": "A", "S": "1000", "I": "0", "R": "0"}, {"deme": "B", "S": "1000", "I": "0", "R": "0"}]
        mobility = [{"origin": "A", "destination": "B", "rate": "0.02"}, {"origin": "B", "destination": "A", "rate": "0.02"}]
        self.assertEqual(audit_mobility(demes, mobility), [])

    def test_self_loop_unidirectional_and_asymmetry(self):
        demes = [{"deme": "A", "S": "1000", "I": "0", "R": "0"}, {"deme": "B", "S": "1000", "I": "0", "R": "0"}]
        self.assertTrue(any("self-loop" in w for w in audit_mobility(demes, [{"origin": "A", "destination": "A", "rate": "0.05"}])))
        self.assertTrue(any("Unidirectional mobility" in w for w in audit_mobility(demes, [{"origin": "A", "destination": "B", "rate": "0.05"}])))
        asymmetric = [{"origin": "A", "destination": "B", "rate": "0.05"}, {"origin": "B", "destination": "A", "rate": "0.01"}]
        self.assertTrue(any("Asymmetric mobility" in w and "5.0x" in w for w in audit_mobility(demes, asymmetric)))

    def test_bdm_flux_rate_and_propensity_warnings(self):
        demes = [{"deme": "A", "S": "50000", "I": "0", "R": "0"}, {"deme": "B", "S": "50000", "I": "0", "R": "0"}]
        mobility = [{"origin": "A", "destination": "B", "rate": "2.5"}, {"origin": "B", "destination": "A", "rate": "0.2"}]
        warnings = audit_mobility(demes, mobility)
        self.assertTrue(any("Forward demographic flux" in w for w in warnings))
        self.assertTrue(any("Large per-capita mobility rate" in w for w in warnings))
        self.assertTrue(any("High aggregate forward migration propensity" in w for w in warnings))

    def test_coalescent_audit_does_not_infer_demographic_flux(self):
        demes = [{"deme": "A", "popSize": "100"}, {"deme": "B", "popSize": "500"}]
        mobility = [{"origin": "A", "destination": "B", "rate": "0.5"}, {"origin": "B", "destination": "A", "rate": "0.01"}]
        warnings = audit_mobility(demes, mobility, mode="coalescent")
        self.assertFalse(any("demographic flux" in w.lower() or "lineage concentration" in w.lower() for w in warnings))


class TestGeneration(unittest.TestCase):
    def setUp(self):
        self.bdm_demes = [{"deme": "A", "S": "1000", "I": "1", "R": "0", "beta": "0.001", "gamma": "0.1"}]
        self.coal_demes = [{"deme": "A", "popSize": "100"}, {"deme": "B", "popSize": "200"}]
        self.mobility = [{"origin": "A", "destination": "B", "rate": "0.02"}]

    def test_bdm_uses_valid_default_ancestry(self):
        fragment = generate_bdm(self.bdm_demes, [], sample_rate=0.05)
        self.assertIn("S_A + I_A -> 2I_A", fragment)
        self.assertIn("I_A -> sample", fragment)

    def test_backward_lineage_rates_are_used_directly(self):
        fragment = generate_coalescent(self.coal_demes, self.mobility, "backward-lineage")
        self.assertIn('rate="0.02"> A -> B', fragment)

    def test_forward_demographic_rates_are_reversed_and_scaled(self):
        fragment = generate_coalescent(self.coal_demes, self.mobility, "forward-demographic")
        self.assertIn('rate="0.01"> B -> A', fragment)


class TestCLIExecution(unittest.TestCase):
    def setUp(self):
        self.script = str(SCRIPTS_DIR / "generate_spatial_reactions.py")
        self.fixtures = REPO_ROOT / "tests" / "advanced" / "fixtures"

    def run_cli(self, *args: str):
        return subprocess.run([sys.executable, self.script, *args], capture_output=True, text=True)

    def test_clean_audit_and_strict_warning(self):
        common = ("--demes", str(self.fixtures / "demes.csv"), "--mobility", str(self.fixtures / "mobility.csv"))
        clean = self.run_cli(*common, "--audit-only")
        self.assertEqual(clean.returncode, 0)
        self.assertIn("PASS: mobility network audit", clean.stdout)
        strict = self.run_cli("--demes", str(self.fixtures / "demes.csv"), "--mobility", str(self.fixtures / "asymmetric_mobility.csv"), "--strict")
        self.assertEqual(strict.returncode, 1)
        self.assertIn("strict mode rejected", strict.stderr)

    def test_bdm_cli_generation(self):
        result = self.run_cli("--demes", str(self.fixtures / "demes.csv"), "--mobility", str(self.fixtures / "mobility.csv"), "--sample-rate", "0.05")
        self.assertEqual(result.returncode, 0)
        self.assertIn("S_North + I_North -> 2I_North", result.stdout)

    def test_coalescent_cli_requires_convention_and_sizes(self):
        missing_convention = self.run_cli("--demes", str(self.fixtures / "demes_coalescent.csv"), "--mobility", str(self.fixtures / "mobility.csv"), "--mode", "coalescent")
        self.assertEqual(missing_convention.returncode, 2)
        result = self.run_cli("--demes", str(self.fixtures / "demes_coalescent.csv"), "--mobility", str(self.fixtures / "mobility.csv"), "--mode", "coalescent", "--coalescent-rate-convention", "forward-demographic")
        self.assertEqual(result.returncode, 0)
        self.assertIn('rate="0.01"> Central -> North', result.stdout)


if __name__ == "__main__":
    unittest.main()
