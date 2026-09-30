#!/usr/bin/env python3
"""Static semantic checks for the ReMASTER Model Builder test cases."""

from __future__ import annotations

import argparse
import math
import re
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from validate_remaster_xml import ValidationError as GenericValidationError
from validate_remaster_xml import validate as validate_generic


class VerificationError(Exception):
    pass


def fail(message: str) -> None:
    raise VerificationError(message)


def norm_text(value: str | None) -> str:
    return "".join((value or "").replace("→", "->").split())


def norm_condition(value: str | None) -> str:
    return "".join((value or "").split())


def strip_lineage_labels(value: str) -> str:
    return re.sub(r":[A-Za-z0-9_]+", "", value)


def stoichiometric_norm_reaction(text: str | None) -> str:
    """Normalize labels and equivalent coefficient/expanded term spellings."""
    text = strip_lineage_labels(norm_text(text))
    if "->" not in text:
        return text

    def normalize_side(side: str) -> str:
        order: list[str] = []
        counts: dict[str, int] = {}
        for term in side.split("+"):
            match = re.fullmatch(r"(\d*)([A-Za-z_][A-Za-z0-9_]*)", term)
            if not match:
                return side
            coefficient, population_name = match.groups()
            if population_name not in counts:
                order.append(population_name)
                counts[population_name] = 0
            counts[population_name] += int(coefficient or "1")
        return "+".join(
            f"{counts[name] if counts[name] != 1 else ''}{name}" for name in order
        )

    lhs, rhs = text.split("->", 1)
    return f"{normalize_side(lhs)}->{normalize_side(rhs)}"


def alpha_norm_reaction(text: str | None) -> str:
    text = norm_text(text)
    if "->" not in text:
        return text
    lhs, rhs = text.split("->", 1)
    term_re = re.compile(r"^(\d*)([A-Za-z_][A-Za-z0-9_]*)(?::([A-Za-z0-9_]+))?$")

    def parse_terms(side: str) -> list[tuple[str, str, str | None]]:
        terms = []
        for t in side.split("+"):
            if not t:
                continue
            m = term_re.match(t)
            if m:
                count, pop, label = m.groups()
                terms.append((count or "", pop, label))
            else:
                terms.append(("", t, None))
        return terms

    r_terms = parse_terms(lhs)
    p_terms = parse_terms(rhs)

    r_labels = {lbl for _, _, lbl in r_terms if lbl}
    p_labels = {lbl for _, _, lbl in p_terms if lbl}
    active_labels = r_labels & p_labels

    label_map: dict[str, str] = {}
    idx = 0
    for _, _, lbl in r_terms:
        if lbl in active_labels and lbl not in label_map:
            label_map[lbl] = f"__L{idx}__"
            idx += 1

    def format_terms(terms: list[tuple[str, str, str | None]], is_product: bool = False) -> str:
        order: list[tuple[str, str | None]] = []
        counts: dict[tuple[str, str | None], int] = {}
        for count, pop, lbl in terms:
            if lbl and lbl in label_map:
                normalized_label = label_map[lbl]
            else:
                if is_product and lbl:
                    normalized_label = lbl
                else:
                    normalized_label = None
            key = (pop, normalized_label)
            if key not in counts:
                order.append(key)
                counts[key] = 0
            counts[key] += int(count or "1")
        parts = []
        for pop, label in order:
            count = counts[(pop, label)]
            rendered = f"{count if count != 1 else ''}{pop}"
            if label:
                rendered += f":{label}"
            parts.append(rendered)
        return "+".join(parts)

    return f"{format_terms(r_terms)}->{format_terms(p_terms, is_product=True)}"


def numbers(value: str | None) -> list[float]:
    if value is None:
        return []
    try:
        return [float(item) for item in value.split()]
    except ValueError:
        fail(f"expected numeric vector, got {value!r}")


def same_numbers(actual: str | None, expected: list[float]) -> bool:
    values = numbers(actual)
    return len(values) == len(expected) and all(
        math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)
        for a, b in zip(values, expected)
    )


def elements(root: ET.Element, tag: str) -> list[ET.Element]:
    return [elem for elem in root.iter() if elem.tag.rsplit("}", 1)[-1] == tag]


def one(items: list[ET.Element], description: str) -> ET.Element:
    if len(items) != 1:
        fail(f"expected exactly one {description}; found {len(items)}")
    return items[0]


def by_spec(root: ET.Element, spec: str) -> list[ET.Element]:
    return [elem for elem in root.iter() if elem.get("spec") == spec]


def trajectory(root: ET.Element) -> ET.Element:
    matches = [
        elem
        for elem in root.iter()
        if elem.get("spec")
        in {"StochasticTrajectory", "DeterministicTrajectory", "CoalescentTrajectory"}
    ]
    return one(matches, "trajectory")


def reaction(
    root: ET.Element,
    text: str,
    spec: str = "Reaction",
    ignore_lineage_labels: bool = False,
    **attrs: str,
) -> ET.Element:
    if ignore_lineage_labels:
        wanted = stoichiometric_norm_reaction(text)
    else:
        wanted = alpha_norm_reaction(text)
    candidates = [
        elem
        for elem in elements(root, "reaction")
        if elem.get("spec") == spec
        and (
            stoichiometric_norm_reaction(elem.text)
            if ignore_lineage_labels
            else alpha_norm_reaction(elem.text)
        )
        == wanted
    ]
    numeric_attrs = {"rate", "changeTimes", "n", "p", "times"}
    for elem in candidates:
        if all(
            same_numbers(elem.get(key), numbers(value))
            if key in numeric_attrs
            else elem.get(key) == value
            for key, value in attrs.items()
        ):
            return elem
    rendered = ", ".join(f"{key}={value!r}" for key, value in attrs.items())
    fail(f"missing {spec} {text!r} with {rendered}")


def population(
    root: ET.Element,
    pop_id: str,
    *,
    tag: str = "population",
    spec: str | None = None,
    value: str | None = None,
    pop_size: str | None = None,
    growth_rate: str | None = None,
) -> ET.Element:
    candidates = [elem for elem in elements(root, tag) if elem.get("id") == pop_id]
    elem = one(candidates, f"{tag} {pop_id!r}")
    expected = {
        "spec": spec,
        "value": value,
        "popSize": pop_size,
        "growthRate": growth_rate,
    }
    numeric_attrs = {"value", "popSize", "growthRate"}
    for key, wanted in expected.items():
        actual = elem.get(key)
        if actual is None and key in {"popSize", "growthRate"}:
            children = [child for child in elem if child.tag.rsplit("}", 1)[-1] == key]
            if len(children) == 1:
                actual = children[0].get("value") or children[0].text
        matches = (
            same_numbers(actual, numbers(wanted))
            if wanted is not None and key in numeric_attrs
            else actual == wanted
        )
        if wanted is not None and not matches:
            fail(f"{tag} {pop_id!r} must have {key}={wanted!r}; got {actual!r}")
    return elem


def logger_with(
    root: ET.Element,
    *,
    file_name: str,
    mode: str | None = None,
    log_idref: str | None = None,
) -> ET.Element:
    candidates = [elem for elem in elements(root, "logger") if elem.get("fileName") == file_name]
    elem = one(candidates, f"logger for {file_name!r}")
    if elem.get("spec") not in {None, "Logger"}:
        fail(
            f"logger {file_name!r} must be an outer BEAST Logger; "
            f"got spec={elem.get('spec')!r}"
        )
    if mode is not None and elem.get("mode") != mode:
        fail(f"logger {file_name!r} must have mode={mode!r}")
    if log_idref is not None:
        matches = [
            child
            for child in elem
            if child.tag.rsplit("}", 1)[-1] == "log"
            and child.get("idref") == log_idref
        ]
        one(matches, f"nested <log idref={log_idref!r}> in logger {file_name!r}")
    return elem


def require_attr(elem: ET.Element, key: str, value: str) -> None:
    if elem.get(key) != value:
        fail(f"{elem.tag} must have {key}={value!r}; got {elem.get(key)!r}")


def require_condition(elem: ET.Element, key: str, *accepted: str) -> None:
    actual = norm_condition(elem.get(key))
    wanted = {norm_condition(value) for value in accepted}
    if actual not in wanted:
        fail(f"{elem.tag} must have {key} in {sorted(wanted)!r}; got {actual!r}")


def require_number(elem: ET.Element, key: str, value: float) -> None:
    if not same_numbers(elem.get(key), [value]):
        fail(f"{elem.tag} must have numeric {key}={value}; got {elem.get(key)!r}")


def require_n_sims(root: ET.Element, value: int) -> None:
    require_attr(one(elements(root, "run"), "run"), "nSims", str(value))


def require_counts(
    root: ET.Element,
    *,
    populations: int,
    sample_populations: int,
    reactions: int,
) -> None:
    traj = trajectory(root)
    actual = (
        len(elements(traj, "population")),
        len(elements(traj, "samplePopulation")),
        len(elements(traj, "reaction")),
    )
    expected = (populations, sample_populations, reactions)
    if actual != expected:
        fail(
            "expected "
            f"{populations} populations, {sample_populations} sample populations, "
            f"and {reactions} reactions; got {actual}"
        )


def require_simulated_tree(root: ET.Element) -> ET.Element:
    return one(by_spec(root, "SimulatedTree"), "SimulatedTree")


def require_trajectory_logger(root: ET.Element, file_name: str) -> ET.Element:
    traj_id = trajectory(root).get("id")
    if not traj_id:
        fail("trajectory must have an id when it is logged")
    return logger_with(root, file_name=file_name, log_idref=traj_id)


def require_typed_tree_logger(root: ET.Element, file_name: str) -> ET.Element:
    tree = require_simulated_tree(root)
    tree_id = tree.get("id")
    if not tree_id:
        fail("SimulatedTree must have an id when it is logged")
    outer = logger_with(root, file_name=file_name, mode="tree")
    typed = one(by_spec(outer, "TypedTreeLogger"), "TypedTreeLogger")
    require_attr(typed, "tree", f"@{tree_id}")
    return typed


def population_names(side: str) -> set[str]:
    names: set[str] = set()
    for term in side.split("+"):
        name = re.sub(r"^\d+", "", term)
        name = name.split(":", 1)[0].split("[", 1)[0]
        if name and name != "0":
            names.add(name)
    return names


def general_checks(root: ET.Element) -> None:
    if root.tag.rsplit("}", 1)[-1] != "beast":
        fail("root element must be <beast>")
    if "remaster" not in (root.get("namespace") or ""):
        fail("root namespace must include remaster")
    namespace_parts = set((root.get("namespace") or "").split(":"))
    if any(elem.get("spec") == "RealParameter" for elem in root.iter()):
        if "beast.base.inference.parameter" not in namespace_parts:
            fail(
                "unqualified RealParameter requires beast.base.inference.parameter "
                "in the root namespace"
            )

    run = one(elements(root, "run"), "run element")
    require_attr(run, "spec", "Simulator")
    try:
        if int(run.get("nSims", "0")) < 1:
            fail("Simulator nSims must be a positive integer")
    except ValueError as exc:
        raise VerificationError("Simulator nSims must be an integer") from exc

    one(elements(root, "simulate"), "simulate element")
    traj = trajectory(root)
    if traj.get("spec") == "DeterministicTrajectory" and traj.get("maxTime") is None:
        fail("DeterministicTrajectory requires maxTime")
    if traj.get("spec") == "DeterministicTrajectory" and "==" in traj.get("endsWhen", ""):
        fail("DeterministicTrajectory endsWhen must not use equality")

    ids: list[str] = []
    sample_ids: set[str] = set()
    for tag in ("population", "samplePopulation"):
        for elem in elements(traj, tag):
            pop_id = elem.get("id")
            if not pop_id:
                fail(f"{tag} is missing an id")
            ids.append(pop_id)
            if tag == "samplePopulation":
                sample_ids.add(pop_id)
    if len(ids) != len(set(ids)):
        fail("population ids must be unique")

    for elem in elements(traj, "reaction"):
        text = norm_text(elem.text)
        if "->" not in text:
            fail(f"reaction is missing an arrow: {text!r}")
        lhs_names = population_names(text.split("->", 1)[0])
        for sample_id in sample_ids:
            if sample_id in lhs_names:
                fail(f"sample population {sample_id!r} appears on a reactant side")
        spec = elem.get("spec")
        if spec == "Reaction":
            rates = numbers(elem.get("rate"))
            if not rates:
                fail(f"continuous reaction {text!r} is missing rate")
            changes = numbers(elem.get("changeTimes"))
            if changes and len(rates) != len(changes) + 1:
                fail(f"reaction {text!r} needs one more rate than changeTimes")
        elif spec == "PunctualReaction":
            has_n = elem.get("n") is not None
            has_p = elem.get("p") is not None
            if has_n == has_p:
                fail(f"punctual reaction {text!r} must have exactly one of n or p")
            if elem.get("times") is None:
                fail(f"punctual reaction {text!r} is missing times")
        else:
            fail(f"reaction {text!r} has unsupported spec {spec!r}")

    if traj.get("spec") == "CoalescentTrajectory":
        if elements(traj, "samplePopulation"):
            fail("coalescent models must not contain samplePopulation")
        if "beast.base.evolution.tree.coalescent" not in namespace_parts:
            fail("coalescent namespace is missing")

    for logger in elements(root, "logger"):
        if logger.get("spec") not in {None, "Logger"}:
            fail("outer <logger> must use spec='Logger' or omit spec")
        if not any(child.tag.rsplit("}", 1)[-1] == "log" for child in logger):
            fail("outer <logger> must contain a nested <log>")


def verify_T01(root: ET.Element) -> None:
    require_n_sims(root, 25)
    require_counts(root, populations=1, sample_populations=1, reactions=3)
    require_simulated_tree(root)
    traj = trajectory(root)
    require_attr(traj, "spec", "StochasticTrajectory")
    require_condition(traj, "endsWhen", "sample==12", "sample>=12")
    require_condition(traj, "mustHave", "sample==12", "sample>=12")
    population(root, "X", spec="RealParameter", value="2")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    reaction(root, "X -> 2X", rate="1.6", ignore_lineage_labels=True)
    reaction(root, "X -> 0", rate="0.5", ignore_lineage_labels=True)
    reaction(root, "X -> sample", rate="0.2", ignore_lineage_labels=True)
    require_trajectory_logger(root, "serial-bds.traj")
    require_typed_tree_logger(root, "serial-bds.trees")


def verify_T02(root: ET.Element) -> None:
    require_n_sims(root, 20)
    require_counts(root, populations=3, sample_populations=0, reactions=2)
    traj = trajectory(root)
    require_attr(traj, "spec", "StochasticTrajectory")
    require_number(traj, "maxTime", 60)
    for pop_id, value in (("S", "999"), ("I", "1"), ("R", "0")):
        population(root, pop_id, spec="RealParameter", value=value)
    inf = reaction(root, "S + I -> 2I")
    if not same_numbers(inf.get("rate"), [0.0004, 0.0001]):
        fail("infection rates must be 0.0004 then 0.0001")
    if not same_numbers(inf.get("changeTimes"), [20]):
        fail("infection changeTimes must be 20")
    reaction(root, "I -> R", rate="0.2")
    require_trajectory_logger(root, "sir.traj")


def verify_T03(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=1, sample_populations=0, reactions=2)
    traj = trajectory(root)
    require_attr(traj, "spec", "DeterministicTrajectory")
    require_number(traj, "maxTime", 10)
    population(root, "X", spec="RealParameter", value="1000")
    reaction(root, "X -> 2X", rate="1.2")
    reaction(root, "X -> 0", rate="1")
    require_trajectory_logger(root, "deterministic.traj")


def verify_T04(root: ET.Element) -> None:
    require_n_sims(root, 1)
    traj = trajectory(root)
    if len(elements(traj, "population")) != 1 or len(elements(traj, "samplePopulation")) != 1:
        fail("T04 requires one population and one sample population")
    events = [
        elem
        for elem in elements(traj, "reaction")
        if elem.get("spec") == "PunctualReaction"
        and norm_text(elem.text) == norm_text("X -> sample")
    ]
    if len(events) not in {1, 3} or len(elements(traj, "reaction")) != len(events):
        fail("T04 requires either one vectorized or three scalar punctual sampling reactions")
    require_attr(traj, "spec", "StochasticTrajectory")
    require_number(traj, "maxTime", 5)
    population(root, "X", spec="RealParameter", value="100")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    if len(events) == 1:
        if not same_numbers(events[0].get("n"), [5, 7, 10]):
            fail("fixed sampling counts must be 5 7 10")
        if not same_numbers(events[0].get("times"), [1, 2, 5]):
            fail("fixed sampling times must be 1 2 5")
    else:
        observed = sorted(
            (numbers(event.get("times"))[0], numbers(event.get("n"))[0])
            for event in events
            if len(numbers(event.get("times"))) == 1 and len(numbers(event.get("n"))) == 1
        )
        if observed != [(1.0, 5.0), (2.0, 7.0), (5.0, 10.0)]:
            fail("scalar punctual reactions must pair times 1,2,5 with counts 5,7,10")
    require_trajectory_logger(root, "fixed-samples.traj")


def verify_T05(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=1, sample_populations=1, reactions=1)
    traj = trajectory(root)
    require_attr(traj, "spec", "StochasticTrajectory")
    require_number(traj, "maxTime", 5)
    population(root, "X", spec="RealParameter", value="100")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    event = reaction(root, "X -> X + sample", spec="PunctualReaction")
    require_attr(event, "p", "0.5")
    if not same_numbers(event.get("times"), [1, 2, 5]):
        fail("probabilistic sampling times must be 1 2 5")
    if any(norm_text(elem.text) == norm_text("X -> sample") for elem in elements(root, "reaction")):
        fail("non-removing sampling must not use X -> sample")
    require_trajectory_logger(root, "probability-samples.traj")


def verify_T06(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=2, sample_populations=1, reactions=2)
    traj = trajectory(root)
    require_attr(traj, "spec", "StochasticTrajectory")
    require_condition(traj, "endsWhen", "sample==10", "sample>=10")
    population(root, "S", spec="RealParameter", value="999")
    population(root, "I", spec="RealParameter", value="1")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    reaction(root, "S:a + I:b -> 2I:b", rate="0.1")
    reaction(root, "I -> sample", rate="0.05")
    require_simulated_tree(root)
    require_typed_tree_logger(root, "labelled.trees")


def verify_T07(root: ET.Element) -> None:
    require_n_sims(root, 10)
    require_counts(root, populations=1, sample_populations=1, reactions=3)
    traj = trajectory(root)
    require_attr(traj, "spec", "StochasticTrajectory")
    require_condition(traj, "endsWhen", "sample==20", "sample>=20")
    require_condition(traj, "mustHave", "X>=1")
    require_number(traj, "maxTime", 50)
    population(root, "X", spec="RealParameter", value="1")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    reaction(root, "X -> 2X", rate="2")
    reaction(root, "X -> 0", rate="0.5")
    reaction(root, "X -> X + sample", rate="0.2")
    require_trajectory_logger(root, "conditioned.traj")


def verify_T08(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=2, sample_populations=1, reactions=4)
    require_simulated_tree(root)
    require_number(trajectory(root), "maxTime", 10)
    for pop_id, value in (("X", "1"), ("Y", "0")):
        population(root, pop_id, spec="RealParameter", value=value)
    population(root, "samp", tag="samplePopulation", spec="RealParameter", value="0")
    reaction(root, "X -> 2X", rate="1.4", ignore_lineage_labels=True)
    reaction(root, "X -> Y", rate="0.1", ignore_lineage_labels=True)
    reaction(root, "X -> 0", rate="0.5", ignore_lineage_labels=True)
    reaction(root, "Y -> samp", rate="0.5", ignore_lineage_labels=True)
    require_typed_tree_logger(root, "typed.trees")


def verify_T09(root: ET.Element) -> None:
    require_n_sims(root, 10)
    require_counts(root, populations=1, sample_populations=1, reactions=3)
    require_simulated_tree(root)
    traj = trajectory(root)
    require_number(traj, "maxTime", 10)
    require_condition(traj, "mustHave", "sample>=2")
    population(root, "X", spec="RealParameter", value="1")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    reaction(root, "X -> 2X", rate="1.4")
    reaction(root, "X -> 0", rate="0.5")
    typed = require_typed_tree_logger(root, "variable.trees")
    require_attr(typed, "noLabels", "true")
    reaction(root, "X -> X + sample", rate="0.5")


def verify_T10(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=1, sample_populations=0, reactions=1)
    traj = trajectory(root)
    require_attr(traj, "spec", "CoalescentTrajectory")
    population(root, "pop", spec="ConstantPopulation", pop_size="1.0")
    tips = reaction(root, "0 -> pop", spec="PunctualReaction")
    require_attr(tips, "n", "10")
    require_attr(tips, "times", "0")
    require_simulated_tree(root)
    require_typed_tree_logger(root, "constant-coalescent.trees")


def verify_T11(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=2, sample_populations=0, reactions=3)
    require_attr(trajectory(root), "spec", "CoalescentTrajectory")
    require_simulated_tree(root)
    population(root, "C", spec="ConstantPopulation", pop_size="10.0")
    population(root, "E", spec="ExponentialGrowth", pop_size="100.0", growth_rate="1")
    reaction(root, "C -> E", rate="0.5")
    c_tips = reaction(root, "0 -> C", spec="PunctualReaction")
    e_tips = reaction(root, "0 -> E", spec="PunctualReaction")
    for event in (c_tips, e_tips):
        require_attr(event, "n", "50")
        require_attr(event, "times", "0")
    require_typed_tree_logger(root, "structured-coalescent.trees")


def verify_T12(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=1, sample_populations=0, reactions=1)
    population(root, "X", spec="RealParameter", value="100")
    event = reaction(root, "2X -> 3X")
    if not same_numbers(event.get("rate"), [0.4]):
        fail("MASTER rate 0.8 with 2 identical reactants must become ReMASTER rate 0.4")
    traj = trajectory(root)
    require_attr(traj, "spec", "StochasticTrajectory")
    require_number(traj, "maxTime", 1)
    require_attr(traj, "endsWhen", "X==110")
    require_trajectory_logger(root, "converted.traj")


def verify_T13(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=9, sample_populations=1, reactions=17)
    require_simulated_tree(root)
    traj = trajectory(root)
    require_number(traj, "maxTime", 30)
    # Check within-deme populations
    for d, s_val, i_val in (("North", "1000", "1"), ("Central", "2000", "0"), ("South", "1500", "0")):
        population(root, f"S_{d}", spec="RealParameter", value=s_val)
        population(root, f"I_{d}", spec="RealParameter", value=i_val)
        population(root, f"R_{d}", spec="RealParameter", value="0")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    # Same-type ancestry defaults already assign both I products to the I reactant.
    reaction(root, "S_North + I_North -> 2I_North", rate="0.0005", ignore_lineage_labels=True)
    reaction(root, "S_Central + I_Central -> 2I_Central", rate="0.0003", ignore_lineage_labels=True)
    reaction(root, "S_South + I_South -> 2I_South", rate="0.0004", ignore_lineage_labels=True)
    # Recoveries & sampling
    for d in ("North", "Central", "South"):
        reaction(root, f"I_{d} -> R_{d}", rate="0.1")
        reaction(root, f"I_{d} -> sample", rate="0.05")
    # Between-deme migrations
    for orig, dest, rate in (
        ("North", "Central", "0.02"),
        ("Central", "North", "0.02"),
        ("Central", "South", "0.015"),
        ("South", "Central", "0.015"),
    ):
        reaction(root, f"S_{orig} -> S_{dest}", rate=rate, ignore_lineage_labels=True)
        reaction(root, f"I_{orig} -> I_{dest}", rate=rate, ignore_lineage_labels=True)
    require_trajectory_logger(root, "spatial.traj")
    require_typed_tree_logger(root, "spatial.trees")


def verify_T14(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=3, sample_populations=0, reactions=7)
    traj = trajectory(root)
    require_attr(traj, "spec", "CoalescentTrajectory")
    require_simulated_tree(root)
    population(root, "North", spec="ConstantPopulation", pop_size="100.0")
    population(root, "Central", spec="ConstantPopulation", pop_size="200.0")
    population(root, "South", spec="ConstantPopulation", pop_size="150.0")
    # Convert forward demographic i->j rates to backward lineage j->i rates Ni*mij/Nj.
    reaction(root, "Central -> North", rate="0.01")
    reaction(root, "North -> Central", rate="0.04")
    reaction(root, "South -> Central", rate="0.02")
    reaction(root, "Central -> South", rate="0.01125")
    # Sampling 10 lineages in each deme at time 0
    for d in ("North", "Central", "South"):
        tips = reaction(root, f"0 -> {d}", spec="PunctualReaction")
        require_attr(tips, "n", "10")
        require_attr(tips, "times", "0")
    require_typed_tree_logger(root, "coalescent_spatial.trees")


def verify_T15(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=3, sample_populations=1, reactions=2)
    require_simulated_tree(root)
    traj = trajectory(root)
    require_attr(traj, "endsWhen", "sample==10")
    population(root, "S", spec="RealParameter", value="999")
    population(root, "I", spec="RealParameter", value="1")
    population(root, "E", spec="RealParameter", value="0")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    # E has no same-type reactant, so an explicit label is required to make I its parent.
    reaction(root, "S:a + I:b -> E:b + I:b", rate="0.1")
    reaction(root, "E -> sample", rate="0.05")
    require_typed_tree_logger(root, "labelled.trees")


def verify_T16(root: ET.Element) -> None:
    require_n_sims(root, 1)
    traj = trajectory(root)
    require_attr(traj, "endsWhen", "sample==50")
    # Must have a positive maxTime safeguard to prevent infinite hang
    val = traj.get("maxTime")
    if val is None or float(val) <= 0:
        fail("repaired model must include a positive maxTime safeguard on trajectory")
    population(root, "X", spec="RealParameter", value="10")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    reaction(root, "X -> 2X", rate="2.0")
    # Must have a sampling reaction so sample==50 is reachable
    sample_rxn = [
        e for e in elements(root, "reaction")
        if norm_text(e.text) in {"X->sample", "X->X+sample"}
    ]
    if not sample_rxn:
        fail("repaired model must include a reaction generating sample (e.g. X -> sample)")
    require_trajectory_logger(root, "gillespie.traj")


def verify_T17(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=1, sample_populations=0, reactions=1)
    population(root, "X", spec="RealParameter", value="100")
    event = reaction(root, "3X -> 4X")
    # Rate 0.6 per ordered sequence divided by 3! = 6 equals 0.1
    if not same_numbers(event.get("rate"), [0.1]):
        fail("MASTER rate 0.6 with 3 identical reactants must become ReMASTER rate 0.1 (0.6 / 3!)")
    traj = trajectory(root)
    require_attr(traj, "spec", "StochasticTrajectory")
    require_number(traj, "maxTime", 1)
    require_attr(traj, "endsWhen", "X==120")
    require_trajectory_logger(root, "trimolecular.traj")


def verify_T18(root: ET.Element) -> None:
    require_n_sims(root, 1)
    require_counts(root, populations=4, sample_populations=1, reactions=4)
    require_simulated_tree(root)
    traj = trajectory(root)
    require_attr(traj, "endsWhen", "sample==10")
    require_attr(traj, "mustHave", "sample==10")
    population(root, "S_A", spec="RealParameter", value="500")
    population(root, "I_A", spec="RealParameter", value="1")
    population(root, "S_B", spec="RealParameter", value="500")
    population(root, "I_B", spec="RealParameter", value="0")
    population(root, "sample", tag="samplePopulation", spec="RealParameter", value="0")
    # Commuter transmission must have explicit ancestry
    reaction(root, "S_B:a + I_A:b -> I_B:b + I_A:b", rate="0.0002")
    reaction(root, "S_A + I_A -> 2I_A", rate="0.0003", ignore_lineage_labels=True)
    reaction(root, "I_A -> sample", rate="0.05", ignore_lineage_labels=True)
    reaction(root, "I_B -> sample", rate="0.05", ignore_lineage_labels=True)
    require_typed_tree_logger(root, "commuter.trees")


CASE_VERIFIERS = {
    case_id: globals()[f"verify_{case_id}"] for case_id in [f"T{i:02d}" for i in range(1, 19)]
}


NEGATIVE_MUTATIONS = {
    "T01": ("rate=\"0.2\"", "rate=\"0.21\""),
    "T02": ("changeTimes=\"20\"", "changeTimes=\"21\""),
    "T03": ("DeterministicTrajectory", "StochasticTrajectory"),
    "T04": ("n=\"5 7 10\"", "n=\"5 7 9\""),
    "T05": ("X -> X + sample", "X -> sample"),
    "T06": ("S:a + I:b -> 2I:b", "S:a + I:b -> 2I:a"),
    "T07": ("mustHave=\"X&gt;=1\"", "mustHave=\"X&gt;=2\""),
    "T08": ("spec=\"TypedTreeLogger\"", "spec=\"TreeStatLogger\""),
    "T09": ("noLabels=\"true\"", "noLabels=\"false\""),
    "T10": ("popSize=\"1.0\"", "popSize=\"2.0\""),
    "T11": ("C -> E", "E -> C"),
    "T12": ("rate=\"0.4\"", "rate=\"0.8\""),
    "T13": ("rate=\"0.02\"", "rate=\"0.03\""),
    "T14": ("rate=\"0.01\"> Central -> North <", "rate=\"0.02\"> Central -> North <"),
    "T15": ("> S:a + I:b -> E:b + I:b <", "> S + I -> E + I <"),
    "T16": ("maxTime=\"20\"", "maxTime=\"-1\""),
    "T17": ("rate=\"0.1\"", "rate=\"0.6\""),
    "T18": ("S_B:a + I_A:b -> I_B:b + I_A:b", "S_B + I_A -> I_B + I_A"),
}


def parse_xml(path: Path) -> ET.Element:
    try:
        return ET.parse(path).getroot()
    except (ET.ParseError, OSError) as exc:
        raise VerificationError(f"could not parse {path}: {exc}") from exc


def verify(path: Path, case_id: str | None = None) -> None:
    try:
        validate_generic(path)
    except GenericValidationError as exc:
        fail(str(exc))
    root = parse_xml(path)
    general_checks(root)
    if case_id:
        CASE_VERIFIERS[case_id](root)


def self_test(skill_dir: Path) -> None:
    gold_dir = skill_dir / "tests" / "gold"
    passed_gold = 0
    rejected_mutations = 0
    for case_id in CASE_VERIFIERS:
        path = gold_dir / f"{case_id}.xml"
        verify(path, case_id)
        passed_gold += 1

        source = path.read_text(encoding="utf-8")
        old, new = NEGATIVE_MUTATIONS[case_id]
        if old not in source:
            fail(f"self-test mutation target missing for {case_id}: {old!r}")
        mutated = source.replace(old, new, 1)
        with tempfile.TemporaryDirectory(prefix=f"remaster-{case_id}-") as temp_dir:
            bad_path = Path(temp_dir) / f"{case_id}-bad.xml"
            bad_path.write_text(mutated, encoding="utf-8")
            try:
                verify(bad_path, case_id)
            except VerificationError:
                rejected_mutations += 1
            else:
                fail(f"{case_id} verifier accepted its known-bad mutation")

    print(f"PASS: {passed_gold} gold fixtures accepted")
    print(f"PASS: {rejected_mutations} known-bad mutations rejected")
    alternative_dir = skill_dir / "tests" / "alternatives"
    alternatives = sorted(alternative_dir.glob("T*.xml"))
    for path in alternatives:
        case_id = path.stem[:3]
        verify(path, case_id)
    print(f"PASS: {len(alternatives)} executable alternative fixtures accepted")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml", type=Path, help="ReMASTER XML to verify")
    parser.add_argument("--case", choices=sorted(CASE_VERIFIERS), help="test case")
    parser.add_argument("--self-test", action="store_true", help="test all case oracles")
    args = parser.parse_args()

    try:
        if args.self_test:
            self_test(Path(__file__).resolve().parents[1])
        elif args.xml:
            verify(args.xml, args.case)
            suffix = f" satisfies {args.case}" if args.case else ""
            print(f"PASS: {args.xml}{suffix}")
        else:
            parser.error("--xml is required unless --self-test is used")
    except VerificationError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
