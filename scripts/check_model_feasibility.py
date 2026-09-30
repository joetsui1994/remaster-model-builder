#!/usr/bin/env python3
"""Conservative static checks for common ReMASTER model feasibility errors."""

from __future__ import annotations

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


TERM_RE = re.compile(r"^\s*(?:(\d+)\s*)?([A-Za-z_][A-Za-z0-9_]*)(?::([A-Za-z_][A-Za-z0-9_]*))?\s*$")
IDENTIFIER_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
SAMPLE_TARGET_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*==\s*(\d+)\s*$")


class FeasibilityReport:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.metrics: dict[str, float | str] = {}

    @property
    def has_issues(self) -> bool:
        return bool(self.errors or self.warnings)


@dataclass(frozen=True)
class ReactionTerm:
    species: str
    count: int = 1
    label: str | None = None


def compute_r0(beta: float, s0: float, gamma: float, psi: float = 0.0) -> float:
    """Early-phase R0 for homogeneous mass-action SIR with removal rates gamma + psi."""
    values = (beta, s0, gamma, psi)
    if not all(math.isfinite(value) and value >= 0 for value in values):
        raise ValueError("beta, s0, gamma, and psi must be finite and non-negative")
    denominator = gamma + psi
    return float("inf") if denominator == 0 else beta * s0 / denominator


def compute_takeoff_prob(r0: float, i0: int = 1) -> float:
    """Birth-death branching approximation 1-(1/R0)^i0, not a general epidemic formula."""
    if math.isnan(r0) or r0 < 0 or i0 < 1 or not isinstance(i0, int):
        raise ValueError("r0 must be non-negative and i0 must be a positive integer")
    if r0 <= 1.0:
        return 0.0
    return 1.0 - (1.0 / r0) ** i0


def _parse_side(side: str) -> list[ReactionTerm]:
    if side.strip() == "0":
        return []
    terms: list[ReactionTerm] = []
    for raw_term in side.split("+"):
        match = TERM_RE.fullmatch(raw_term)
        if not match:
            return []
        count, species, label = match.groups()
        terms.append(ReactionTerm(species, int(count or 1), label))
    return terms


def _parse_reaction(text: str) -> tuple[list[ReactionTerm], list[ReactionTerm]] | None:
    if text.count("->") != 1:
        return None
    left, right = text.split("->")
    reactants, products = _parse_side(left), _parse_side(right)
    if (left.strip() != "0" and not reactants) or (right.strip() != "0" and not products):
        return None
    return reactants, products


def _counts(terms: list[ReactionTerm]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for term in terms:
        counts[term.species] += term.count
    return counts


def _ancestry_warnings(text: str, reactants: list[ReactionTerm], products: list[ReactionTerm]) -> list[str]:
    """Flag fallback-to-first-reactant products only when a tree reaction has multiple reactants."""
    if len(reactants) < 2:
        return []
    warnings: list[str] = []
    lhs_labels = {term.label for term in reactants if term.label}
    lhs_species = {term.species for term in reactants}
    first_species = reactants[0].species
    for product in products:
        if product.label and product.label in lhs_labels:
            continue
        if product.species in lhs_species:
            continue
        warnings.append(
            f"Lineage ancestry fallback in reaction '{text}': product '{product.species}' has no same-type "
            f"reactant or matching label, so ReMASTER assigns it to the first reactant '{first_species}'. "
            f"Add matching labels if a different reactant is the intended parent."
        )
    return warnings


def _local_name(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def check_xml_feasibility(path: Path) -> FeasibilityReport:
    """Analyze one ReMASTER XML without claiming to prove biological correctness."""
    report = FeasibilityReport()
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as exc:
        report.errors.append(f"Failed to parse XML: {exc}")
        return report

    elements = list(root.iter())
    trajectories = [element for element in elements if element.get("spec", "").endswith("Trajectory")]
    stochastic = [element for element in trajectories if element.get("spec") == "StochasticTrajectory"]
    coalescent = [element for element in trajectories if element.get("spec") == "CoalescentTrajectory"]
    tree_requested = any(element.get("spec") == "SimulatedTree" for element in elements)
    sample_populations = [element for element in elements if _local_name(element) == "samplePopulation"]
    sample_ids = {element.get("id") for element in sample_populations if element.get("id")}

    # DeterministicTrajectory + SimulatedTree is supported: it samples a diffusion-approximation tree.
    if coalescent and sample_populations:
        report.errors.append(
            "Contradictory specification: CoalescentTrajectory is backward-time and should introduce sampled "
            "lineages with punctual reactions, not samplePopulation."
        )

    reaction_elements = [element for element in elements if element.get("spec") in {"Reaction", "PunctualReaction"}]
    parsed_reactions: list[tuple[str, list[ReactionTerm], list[ReactionTerm]]] = []
    produced: set[str] = set()
    for element in reaction_elements:
        text = " ".join((element.text or "").split())
        parsed = _parse_reaction(text)
        if parsed is None:
            continue
        reactants, products = parsed
        produced.update(term.species for term in products)
        parsed_reactions.append((text, reactants, products))
        if tree_requested and element.get("spec") == "Reaction":
            report.warnings.extend(_ancestry_warnings(text, reactants, products))

    initial_values: dict[str, float] = {}
    for element in elements:
        if _local_name(element) not in {"population", "samplePopulation"} or not element.get("id"):
            continue
        value = element.get("value")
        if value is None:
            continue
        try:
            initial_values[element.get("id", "")] = float(value)
        except ValueError:
            pass

    for trajectory in trajectories:
        for attribute in ("endsWhen", "mustHave"):
            condition = trajectory.get(attribute, "")
            for identifier in IDENTIFIER_RE.findall(condition):
                if identifier in initial_values and initial_values[identifier] <= 0 and identifier not in produced:
                    report.errors.append(
                        f"Unreachable {attribute} condition '{condition}': population '{identifier}' starts at "
                        f"{initial_values[identifier]:g} and no reaction produces it."
                    )

        if trajectory.get("spec") == "StochasticTrajectory" and tree_requested:
            target = SAMPLE_TARGET_RE.fullmatch(trajectory.get("endsWhen", ""))
            if target and target.group(1) in sample_ids:
                sample_id, count = target.groups()
                acceptance = trajectory.get("mustHave", "")
                accepted_target = re.search(
                    rf"\b{re.escape(sample_id)}\s*(?:==|>=)\s*{re.escape(count)}\b",
                    acceptance,
                )
                if not accepted_target:
                    report.warnings.append(
                        f"Unconditioned sample target: endsWhen='{sample_id}=={count}' stops at the requested "
                        f"sample count, but without a matching mustHave condition extinction can yield an "
                        f"accepted tree with fewer than {count} tips."
                    )

    # Warn only for a species that can increase but is never consumed and has no explicit stopping bound.
    for trajectory in stochastic:
        if trajectory.get("maxTime") or trajectory.get("endsWhen"):
            continue
        local_reactions: list[tuple[list[ReactionTerm], list[ReactionTerm]]] = []
        for element in trajectory.iter():
            if element.get("spec") != "Reaction":
                continue
            parsed = _parse_reaction(" ".join((element.text or "").split()))
            if parsed:
                local_reactions.append(parsed)
        deltas: dict[str, list[int]] = {}
        for reactants, products in local_reactions:
            left, right = _counts(reactants), _counts(products)
            for species in left.keys() | right.keys():
                deltas.setdefault(species, []).append(right[species] - left[species])
        unbounded = sorted(species for species, values in deltas.items() if any(v > 0 for v in values) and not any(v < 0 for v in values))
        if unbounded:
            report.warnings.append(
                "Unbounded stochastic growth risk: population(s) " + ", ".join(unbounded) +
                " can increase but are never consumed, and the trajectory has neither maxTime nor endsWhen."
            )

    for population in elements:
        if _local_name(population) != "population" or population.get("spec") != "RealParameter":
            continue
        try:
            value = float(population.get("value", "0"))
        except ValueError:
            continue
        if value >= 50000 and stochastic:
            report.warnings.append(
                f"Large population scale: '{population.get('id')}' starts at {value:.0f}. Exact stochastic "
                f"simulation may be expensive; benchmark it or justify rescaling/approximation."
            )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Conservative static checks for ReMASTER feasibility hazards.")
    parser.add_argument("xml", type=Path)
    parser.add_argument("--strict", action="store_true", help="Return nonzero for warnings as well as errors")
    args = parser.parse_args()

    report = check_xml_feasibility(args.xml)
    if report.errors:
        print(f"ERROR: {len(report.errors)} critical issue(s) found in {args.xml}:", file=sys.stderr)
        for error in report.errors:
            print(f"  [!] {error}", file=sys.stderr)
        return 1
    if report.warnings:
        print(f"WARNING: {len(report.warnings)} feasibility hazard(s) found in {args.xml}:", file=sys.stderr)
        for warning in report.warnings:
            print(f"  [*] {warning}", file=sys.stderr)
        return 1 if args.strict else 0
    print(f"PASS: {args.xml} passed the configured static feasibility checks.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
