#!/usr/bin/env python3
"""Generic static checks for a self-contained ReMASTER simulation XML."""

from __future__ import annotations

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


class ValidationError(Exception):
    pass


def tag_name(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def elements(root: ET.Element, name: str) -> list[ET.Element]:
    return [element for element in root.iter() if tag_name(element) == name]


def normalise_reaction(text: str | None) -> str:
    return "".join((text or "").replace("→", "->").split())


def numeric_values(value: str | None) -> list[float]:
    if value is None:
        return []
    try:
        return [float(item) for item in value.split()]
    except ValueError as exc:
        raise ValidationError(f"expected numeric vector, got {value!r}") from exc


def population_names(side: str) -> set[str]:
    names: set[str] = set()
    for term in side.split("+"):
        name = re.sub(r"^\d+", "", term)
        name = name.split(":", 1)[0].split("[", 1)[0]
        if name and name != "0":
            names.add(name)
    return names


def validate(path: Path) -> list[str]:
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError) as exc:
        raise ValidationError(f"could not parse {path}: {exc}") from exc

    if tag_name(root) != "beast":
        raise ValidationError("root element must be <beast>")
    version = root.get("version")
    try:
        if version is None or float(version) < 2.0:
            raise ValidationError("root <beast> requires a version of at least 2.0")
    except ValueError as exc:
        raise ValidationError(f"root <beast> has invalid version {version!r}") from exc
    namespace = root.get("namespace", "")
    if "remaster" not in namespace:
        raise ValidationError("root namespace must include remaster")
    namespace_parts = set(namespace.split(":"))
    if any(element.get("spec") == "RealParameter" for element in root.iter()):
        if "beast.base.inference.parameter" not in namespace_parts:
            raise ValidationError(
                "unqualified RealParameter requires beast.base.inference.parameter "
                "in the root namespace"
            )
    if any(element.get("spec") == "Logger" for element in root.iter()):
        if "beast.base.inference" not in namespace_parts:
            raise ValidationError(
                "unqualified Logger requires beast.base.inference in the root namespace"
            )

    runs = elements(root, "run")
    if len(runs) != 1 or runs[0].get("spec") != "Simulator":
        raise ValidationError("expected exactly one <run spec=\"Simulator\">")
    try:
        if int(runs[0].get("nSims", "0")) < 1:
            raise ValidationError("nSims must be a positive integer")
    except ValueError as exc:
        raise ValidationError("nSims must be a positive integer") from exc

    simulations = elements(root, "simulate")
    if len(simulations) != 1:
        raise ValidationError(f"expected exactly one <simulate>; found {len(simulations)}")

    trajectory_specs = {
        "StochasticTrajectory",
        "DeterministicTrajectory",
        "CoalescentTrajectory",
    }
    trajectories = [
        element for element in root.iter() if element.get("spec") in trajectory_specs
    ]
    if len(trajectories) != 1:
        raise ValidationError(f"expected exactly one trajectory; found {len(trajectories)}")
    trajectory = trajectories[0]
    trajectory_spec = trajectory.get("spec")

    if trajectory.get("initialTime") is not None:
        raise ValidationError(
            f"{trajectory_spec} does not accept initialTime; ReMASTER time starts at zero"
        )
    for condition_name in ("endsWhen", "mustHave"):
        if elements(trajectory, condition_name):
            raise ValidationError(
                f"{condition_name} must be an attribute on the trajectory, not a child element"
            )

    if trajectory_spec == "DeterministicTrajectory":
        if trajectory.get("maxTime") is None:
            raise ValidationError("DeterministicTrajectory requires maxTime")
        if "==" in trajectory.get("endsWhen", ""):
            raise ValidationError("deterministic endsWhen must use an inequality, not equality")

    population_ids: list[str] = []
    sample_ids: set[str] = set()
    for name in ("population", "samplePopulation"):
        for population in elements(trajectory, name):
            population_id = population.get("id")
            if not population_id:
                raise ValidationError(f"{name} is missing id")
            population_ids.append(population_id)
            if name == "samplePopulation":
                sample_ids.add(population_id)
            if trajectory_spec != "CoalescentTrajectory":
                if population.get("spec") is None:
                    raise ValidationError(
                        f"birth-death {name} {population_id!r} requires an explicit concrete spec"
                    )
                if population.get("spec") == "RealParameter" and population.get("value") is None:
                    raise ValidationError(
                        f"RealParameter {name} {population_id!r} requires an initial value"
                    )
    if len(population_ids) != len(set(population_ids)):
        raise ValidationError("population ids must be unique")

    if trajectory_spec == "CoalescentTrajectory":
        if sample_ids:
            raise ValidationError("coalescent models must not use samplePopulation")
        if "beast.base.evolution.tree.coalescent" not in namespace_parts:
            raise ValidationError("coalescent population functions require the coalescent namespace")
        for population in elements(trajectory, "population"):
            spec = population.get("spec")
            required_inputs: tuple[str, ...] = ()
            if spec == "ConstantPopulation":
                required_inputs = ("popSize",)
            elif spec == "ExponentialGrowth":
                required_inputs = ("popSize", "growthRate")
            for input_name in required_inputs:
                as_attribute = population.get(input_name) is not None
                as_child = any(tag_name(child) == input_name for child in population)
                if not as_attribute and not as_child:
                    raise ValidationError(
                        f"{spec} population {population.get('id')!r} requires {input_name}"
                    )

    for reaction in elements(trajectory, "reaction"):
        reaction_text = normalise_reaction(reaction.text)
        if "->" not in reaction_text:
            raise ValidationError(f"reaction lacks an arrow: {reaction_text!r}")
        reactants = population_names(reaction_text.split("->", 1)[0])
        for sample_id in sample_ids:
            if sample_id in reactants:
                raise ValidationError(
                    f"sample population {sample_id!r} appears on the reactant side"
                )

        spec = reaction.get("spec")
        if spec == "Reaction":
            rates = numeric_values(reaction.get("rate"))
            if not rates:
                raise ValidationError(f"continuous reaction {reaction_text!r} lacks rate")
            change_times = numeric_values(reaction.get("changeTimes"))
            if change_times and len(rates) != len(change_times) + 1:
                raise ValidationError(
                    f"reaction {reaction_text!r} needs one more rate than changeTimes"
                )
        elif spec == "PunctualReaction":
            has_n = reaction.get("n") is not None
            has_p = reaction.get("p") is not None
            if has_n == has_p:
                raise ValidationError(
                    f"punctual reaction {reaction_text!r} needs exactly one of n or p"
                )
            if reaction.get("times") is None and not elements(reaction, "times"):
                raise ValidationError(f"punctual reaction {reaction_text!r} lacks times")
        else:
            raise ValidationError(
                f"reaction {reaction_text!r} has unsupported spec {spec!r}"
            )

    loggers = elements(root, "logger")
    if not loggers:
        raise ValidationError("simulation must include at least one output logger")
    id_targets = {
        element.get("id"): element for element in root.iter() if element.get("id")
    }
    for logger in loggers:
        if logger.get("spec") not in {None, "Logger"}:
            raise ValidationError(
                "outer <logger> must be a BEAST Logger (spec=\"Logger\" or omitted)"
            )
        if not any(tag_name(child) == "log" for child in logger):
            raise ValidationError("outer <logger> must contain a nested <log>")
        for child in logger:
            if tag_name(child) != "log":
                continue
            idref = child.get("idref")
            if idref is not None and idref not in id_targets:
                raise ValidationError(f"logger references unknown id {idref!r}")
            tree_ref = child.get("tree")
            if child.get("spec") in {"TypedTreeLogger", "TreeStatLogger"}:
                if child.get("trajectory") is not None:
                    raise ValidationError(
                        f"{child.get('spec')} uses tree=\"@tree-id\", not trajectory=..."
                    )
                if tree_ref is None:
                    raise ValidationError(f"{child.get('spec')} requires a tree reference")
                target_id = tree_ref.removeprefix("@")
                target = id_targets.get(target_id)
                if target is None:
                    raise ValidationError(f"tree logger references unknown id {tree_ref!r}")
                if target.get("spec") != "SimulatedTree":
                    raise ValidationError(
                        f"tree logger must reference a SimulatedTree; {tree_ref!r} references "
                        f"{target.get('spec')!r}"
                    )
        nested_tree_logs = [
            element
            for element in logger.iter()
            if element is not logger
            and element.get("spec") in {"TypedTreeLogger", "TreeStatLogger"}
        ]
        if nested_tree_logs and (
            logger.get("mode") != "tree"
        ):
            raise ValidationError(
                "outer logger containing a tree log must use mode=\"tree\""
            )

    simulated_trees = [element for element in root.iter() if element.get("spec") == "SimulatedTree"]
    if any(
        element.get("spec") in {"TypedTreeLogger", "TreeStatLogger"}
        for element in root.iter()
    ):
        if len(simulated_trees) != 1:
            raise ValidationError(
                f"tree output requires exactly one SimulatedTree; found {len(simulated_trees)}"
            )
        if trajectory not in list(simulated_trees[0].iter()):
            raise ValidationError("SimulatedTree must contain the trajectory it reconstructs")

    return [
        f"trajectory={trajectory_spec}",
        f"populations={len(population_ids)}",
        f"reactions={len(elements(trajectory, 'reaction'))}",
        f"loggers={len(elements(root, 'logger'))}",
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("xml", type=Path)
    args = parser.parse_args()
    try:
        summary = validate(args.xml)
    except ValidationError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print(f"PASS: {args.xml}")
    print(" ".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
