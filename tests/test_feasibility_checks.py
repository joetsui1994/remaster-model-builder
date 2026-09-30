#!/usr/bin/env python3
"""Automated tests for scientific feasibility, R0 calculations, and model contradiction audits."""

import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from check_model_feasibility import (
    check_xml_feasibility,
    compute_r0,
    compute_takeoff_prob,
)


class TestEpidemiologicalMath(unittest.TestCase):
    """Test R0 and extinction/takeoff probability formulas."""

    def test_r0_calculation(self):
        # beta=0.0005, S0=1000, gamma=0.1, psi=0.0 -> R0 = 0.5 / 0.1 = 5.0
        r0 = compute_r0(beta=0.0005, s0=1000, gamma=0.1)
        self.assertAlmostEqual(r0, 5.0)

        # with sampling psi=0.05 -> R0 = 0.5 / 0.15 = 3.3333
        r0_sampled = compute_r0(beta=0.0005, s0=1000, gamma=0.1, psi=0.05)
        self.assertAlmostEqual(r0_sampled, 3.3333333333)

    def test_takeoff_probability(self):
        # Subcritical R0 <= 1.0 has 0 takeoff probability
        self.assertEqual(compute_takeoff_prob(r0=0.8, i0=1), 0.0)
        self.assertEqual(compute_takeoff_prob(r0=1.0, i0=5), 0.0)

        # Supercritical R0=2.0, I0=1 -> P(takeoff) = 1 - 1/2 = 0.50
        self.assertAlmostEqual(compute_takeoff_prob(r0=2.0, i0=1), 0.50)

        # Supercritical R0=2.0, I0=3 -> P(takeoff) = 1 - (1/2)^3 = 0.875
        self.assertAlmostEqual(compute_takeoff_prob(r0=2.0, i0=3), 0.875)


class TestModelContradictionsAndHazards(unittest.TestCase):
    """Test XML audits for contradictory specifications and simulation compute traps."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def write_xml(self, content: str) -> Path:
        p = self.dir_path / "model.xml"
        p.write_text(content.strip() + "\n", encoding="utf-8")
        return p

    def test_deterministic_with_tree_is_supported(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">
          <run spec="Simulator" nSims="1">
            <simulate id="tree" spec="SimulatedTree">
              <trajectory id="traj" spec="DeterministicTrajectory" maxTime="10">
                <population id="S" spec="RealParameter" value="1000"/>
                <population id="I" spec="RealParameter" value="1"/>
                <reaction spec="Reaction" rate="0.001"> S + I -> 2I </reaction>
              </trajectory>
            </simulate>
            <logger spec="Logger" mode="tree" fileName="tree.trees">
              <log spec="TypedTreeLogger" tree="@tree"/>
            </logger>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertEqual(report.errors, [])

    def test_contradiction_coalescent_with_sample_population(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster:beast.base.evolution.tree.coalescent">
          <run spec="Simulator" nSims="1">
            <simulate id="tree" spec="SimulatedTree">
              <trajectory id="traj" spec="CoalescentTrajectory">
                <population id="pop" spec="ConstantPopulation" popSize="100.0"/>
                <samplePopulation id="sample" spec="RealParameter" value="0"/>
              </trajectory>
            </simulate>
            <logger spec="Logger" mode="tree" fileName="tree.trees">
              <log spec="TypedTreeLogger" tree="@tree"/>
            </logger>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertTrue(any("CoalescentTrajectory" in e and "samplePopulation" in e for e in report.errors))

    def test_gillespie_explosion_hazard(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">
          <run spec="Simulator" nSims="1">
            <simulate id="traj" spec="StochasticTrajectory">
              <population id="X" spec="RealParameter" value="10"/>
              <!-- Supercritical branching with no removal and no maxTime -->
              <reaction spec="Reaction" rate="2.0"> X -> 2X </reaction>
            </simulate>
            <logger spec="Logger" fileName="out.traj"><log idref="traj"/></logger>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertTrue(any("Unbounded stochastic growth risk" in w for w in report.warnings))

    def test_default_ancestry_handles_unlabeled_sir_transmission(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">
          <run spec="Simulator" nSims="1">
            <simulate id="tree" spec="SimulatedTree">
              <trajectory id="traj" spec="StochasticTrajectory" maxTime="10">
                <population id="S" spec="RealParameter" value="999"/>
                <population id="I" spec="RealParameter" value="1"/>
                <samplePopulation id="sample" spec="RealParameter" value="0"/>
                <reaction spec="Reaction" rate="0.001"> S + I -> 2I </reaction>
                <reaction spec="Reaction" rate="0.05"> I -> sample </reaction>
              </trajectory>
            </simulate>
            <logger spec="Logger" mode="tree" fileName="out.trees"><log spec="TypedTreeLogger" tree="@tree"/></logger>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertFalse(any("Lineage ancestry" in w for w in report.warnings))

    def test_commuter_transmission_ancestry_fallback_warning(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">
          <run spec="Simulator" nSims="1">
            <simulate id="tree" spec="SimulatedTree">
              <trajectory id="traj" spec="StochasticTrajectory" maxTime="10">
                <population id="S_B" spec="RealParameter" value="999"/>
                <population id="I_A" spec="RealParameter" value="1"/>
                <population id="I_B" spec="RealParameter" value="0"/>
                <reaction spec="Reaction" rate="0.001"> S_B + I_A -> I_B + I_A </reaction>
              </trajectory>
            </simulate>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertTrue(any("Lineage ancestry fallback" in w and "I_B" in w for w in report.warnings))

    def test_unreachable_condition_is_error(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">
          <run spec="Simulator" nSims="1">
            <simulate id="traj" spec="StochasticTrajectory" endsWhen="sample==50">
              <population id="X" spec="RealParameter" value="10"/>
              <samplePopulation id="sample" spec="RealParameter" value="0"/>
              <reaction spec="Reaction" rate="2.0"> X -> 2X </reaction>
            </simulate>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertTrue(any("Unreachable endsWhen" in error and "sample" in error for error in report.errors))

    def test_exact_sample_tree_requires_acceptance_condition(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">
          <run spec="Simulator" nSims="1">
            <simulate id="tree" spec="SimulatedTree">
              <trajectory id="traj" spec="StochasticTrajectory" endsWhen="sample==10">
                <population id="X" spec="RealParameter" value="1"/>
                <samplePopulation id="sample" spec="RealParameter" value="0"/>
                <reaction spec="Reaction" rate="1.0"> X -> 2X </reaction>
                <reaction spec="Reaction" rate="0.1"> X -> sample </reaction>
              </trajectory>
            </simulate>
            <logger spec="Logger" mode="tree" fileName="out.trees">
              <log spec="TypedTreeLogger" tree="@tree"/>
            </logger>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertTrue(any("Unconditioned sample target" in warning for warning in report.warnings))

    def test_exact_sample_tree_with_acceptance_condition_is_clean(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">
          <run spec="Simulator" nSims="1">
            <simulate id="tree" spec="SimulatedTree">
              <trajectory id="traj" spec="StochasticTrajectory" endsWhen="sample==10" mustHave="sample==10">
                <population id="X" spec="RealParameter" value="1"/>
                <samplePopulation id="sample" spec="RealParameter" value="0"/>
                <reaction spec="Reaction" rate="1.0"> X -> 2X </reaction>
                <reaction spec="Reaction" rate="0.1"> X -> sample </reaction>
              </trajectory>
            </simulate>
            <logger spec="Logger" mode="tree" fileName="out.trees">
              <log spec="TypedTreeLogger" tree="@tree"/>
            </logger>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertFalse(any("Unconditioned sample target" in warning for warning in report.warnings))

    def test_large_population_scale_warning(self):
        xml = """
        <beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">
          <run spec="Simulator" nSims="1">
            <simulate id="traj" spec="StochasticTrajectory" maxTime="10">
              <population id="City" spec="RealParameter" value="1000000"/>
              <reaction spec="Reaction" rate="0.01"> City -> 0 </reaction>
            </simulate>
            <logger spec="Logger" fileName="out.traj"><log idref="traj"/></logger>
          </run>
        </beast>
        """
        report = check_xml_feasibility(self.write_xml(xml))
        self.assertTrue(any("Large population scale" in w for w in report.warnings))

    def test_gold_fixtures_pass_feasibility(self):
        gold_t18 = REPO_ROOT / "tests" / "gold" / "T18.xml"
        report = check_xml_feasibility(gold_t18)
        self.assertEqual(report.errors, [])
        self.assertEqual(report.warnings, [])


if __name__ == "__main__":
    unittest.main()
