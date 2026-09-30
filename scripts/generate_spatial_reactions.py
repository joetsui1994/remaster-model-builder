#!/usr/bin/env python3
"""Generate audited ReMASTER spatial population and reaction fragments."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path


IDENTIFIER_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")


def _finite_nonnegative(value: str, *, field: str, context: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context}: {field} must be a non-negative finite number, got {value!r}") from exc
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{context}: {field} must be a non-negative finite number, got {value!r}")
    return number


def _positive(value: str, *, field: str, context: str) -> float:
    number = _finite_nonnegative(value, field=field, context=context)
    if number == 0:
        raise ValueError(f"{context}: {field} must be greater than zero")
    return number


def parse_demes(path: Path, mode: str = "bdm") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        required = {"deme", "S", "I", "R", "beta", "gamma"} if mode == "bdm" else {"deme"}
        if not required.issubset(fields):
            missing = ", ".join(sorted(required - fields))
            raise ValueError(f"demes file {path} is missing required column(s): {missing}")
        if mode == "coalescent" and not ({"popSize", "N"} & fields):
            raise ValueError(f"demes file {path} must contain a 'popSize' or 'N' column in coalescent mode")

        rows = list(reader)
        if not rows:
            raise ValueError(f"demes file {path} contains no data rows")

        seen: set[str] = set()
        for row_number, row in enumerate(rows, 2):
            context = f"demes file {path} row {row_number}"
            name = (row.get("deme") or "").strip()
            if not IDENTIFIER_RE.fullmatch(name):
                raise ValueError(
                    f"{context}: deme must be a ReMASTER-safe identifier matching "
                    f"[A-Za-z_][A-Za-z0-9_]*, got {name!r}"
                )
            if name in seen:
                raise ValueError(f"{context}: duplicate deme {name!r}")
            seen.add(name)
            row["deme"] = name

            if mode == "bdm":
                for field in ("S", "I", "R"):
                    value = _finite_nonnegative(row[field], field=field, context=context)
                    if not value.is_integer():
                        raise ValueError(f"{context}: {field} must be an integer count, got {row[field]!r}")
                for field in ("beta", "gamma"):
                    _finite_nonnegative(row[field], field=field, context=context)
            else:
                field = "popSize" if "popSize" in fields else "N"
                _positive(row[field], field=field, context=context)
        return rows


def parse_mobility(path: Path, valid_demes: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"origin", "destination", "rate"}
        if not required.issubset(set(reader.fieldnames or [])):
            missing = ", ".join(sorted(required - set(reader.fieldnames or [])))
            raise ValueError(f"mobility file {path} is missing required column(s): {missing}")

        rows = list(reader)
        seen_edges: set[tuple[str, str]] = set()
        for row_number, row in enumerate(rows, 2):
            context = f"mobility file {path} row {row_number}"
            origin = (row.get("origin") or "").strip()
            destination = (row.get("destination") or "").strip()
            if origin not in valid_demes:
                raise ValueError(f"{context}: unknown origin deme {origin!r}")
            if destination not in valid_demes:
                raise ValueError(f"{context}: unknown destination deme {destination!r}")
            edge = (origin, destination)
            if edge in seen_edges:
                raise ValueError(f"{context}: duplicate directed mobility edge {origin!r} -> {destination!r}")
            seen_edges.add(edge)
            _finite_nonnegative(row["rate"], field="rate", context=context)
            row["origin"] = origin
            row["destination"] = destination
        return rows


def _population_sizes(demes: list[dict[str, str]], mode: str) -> dict[str, float]:
    if mode == "bdm":
        return {d["deme"]: float(d["S"]) + float(d["I"]) + float(d["R"]) for d in demes}
    return {d["deme"]: float(d.get("popSize") or d.get("N") or "") for d in demes}


def audit_mobility(
    demes: list[dict[str, str]],
    mobility: list[dict[str, str]],
    mode: str = "bdm",
    asymmetry_ratio_threshold: float = 3.0,
    flux_imbalance_threshold: float = 0.10,
    large_rate_threshold: float = 1.0,
    large_propensity_threshold: float = 10000.0,
) -> list[str]:
    """Return diagnostic warnings; only BDM rows receive demographic-flux interpretations."""
    warnings: list[str] = []
    pop_sizes = _population_sizes(demes, mode)
    rates: dict[tuple[str, str], float] = {}

    for row_number, row in enumerate(mobility, 2):
        origin = row["origin"]
        destination = row["destination"]
        rate = float(row["rate"])
        if origin == destination:
            warnings.append(f"Row {row_number}: self-loop for deme '{origin}' is redundant.")
            continue
        rates[(origin, destination)] = rate

    names = sorted(pop_sizes)
    for index, left in enumerate(names):
        for right in names[index + 1 :]:
            forward = rates.get((left, right), 0.0)
            reverse = rates.get((right, left), 0.0)
            if forward > 0 and reverse == 0:
                warnings.append(f"Unidirectional mobility: '{left}' -> '{right}' has no reverse edge.")
            elif reverse > 0 and forward == 0:
                warnings.append(f"Unidirectional mobility: '{right}' -> '{left}' has no reverse edge.")
            elif forward > 0 and reverse > 0:
                ratio = max(forward, reverse) / min(forward, reverse)
                if ratio >= asymmetry_ratio_threshold:
                    warnings.append(
                        f"Asymmetric mobility: '{left}' <-> '{right}' rates are {forward:g} vs {reverse:g} "
                        f"({ratio:.1f}x). Confirm this is intentional."
                    )

    total_propensity = 0.0
    if mode == "bdm":
        for name in names:
            size = pop_sizes[name]
            outflow = size * sum(rates.get((name, other), 0.0) for other in names if other != name)
            inflow = sum(pop_sizes[other] * rates.get((other, name), 0.0) for other in names if other != name)
            net = inflow - outflow
            total_propensity += outflow
            if size > 0 and abs(net / size) >= flux_imbalance_threshold:
                direction = "decline" if net < 0 else "growth"
                warnings.append(
                    f"Forward demographic flux in deme '{name}' is {net:+.2f}/time "
                    f"({net / size * 100:+.1f}% of N={size:g} per time unit), implying rapid migration-driven {direction}."
                )

    for (origin, destination), rate in rates.items():
        if rate >= large_rate_threshold:
            warnings.append(
                f"Large per-capita mobility rate: '{origin}' -> '{destination}' is {rate:g}/time. "
                f"Interpret this relative to the stated time unit and other process rates."
            )

    if mode == "bdm" and total_propensity >= large_propensity_threshold:
        warnings.append(
            f"High aggregate forward migration propensity ({total_propensity:.0f} events/time). "
            f"This may make an exact Gillespie simulation expensive; compare it with other reaction propensities."
        )
    return warnings


def generate_bdm(
    demes: list[dict[str, str]], mobility: list[dict[str, str]], sample_rate: float | None = None
) -> str:
    lines = ["<!-- Populations per Deme -->"]
    for deme in demes:
        name = deme["deme"]
        lines.extend(
            [
                f'<population id="S_{name}" spec="RealParameter" value="{deme["S"]}"/>',
                f'<population id="I_{name}" spec="RealParameter" value="{deme["I"]}"/>',
                f'<population id="R_{name}" spec="RealParameter" value="{deme["R"]}"/>',
            ]
        )
    lines.append('<samplePopulation id="sample" spec="RealParameter" value="0"/>')
    lines.append("")
    lines.append("<!-- Within-Deme Transmission and Recovery Reactions -->")
    for deme in demes:
        name = deme["deme"]
        lines.append(f'<reaction spec="Reaction" rate="{deme["beta"]}"> S_{name} + I_{name} -> 2I_{name} </reaction>')
        lines.append(f'<reaction spec="Reaction" rate="{deme["gamma"]}"> I_{name} -> R_{name} </reaction>')
        if sample_rate is not None and sample_rate > 0:
            lines.append(f'<reaction spec="Reaction" rate="{sample_rate:g}"> I_{name} -> sample </reaction>')
    lines.append("")
    lines.append("<!-- Between-Deme Movement / Migration Reactions -->")
    for row in mobility:
        origin, destination, rate = row["origin"], row["destination"], row["rate"]
        lines.append(f'<reaction spec="Reaction" rate="{rate}"> S_{origin} -> S_{destination} </reaction>')
        lines.append(f'<reaction spec="Reaction" rate="{rate}"> I_{origin} -> I_{destination} </reaction>')
    return "\n".join(lines)


def generate_coalescent(
    demes: list[dict[str, str]],
    mobility: list[dict[str, str]],
    rate_convention: str,
) -> str:
    sizes = _population_sizes(demes, "coalescent")
    lines = ["<!-- Coalescent Deme Populations -->"]
    for deme in demes:
        name = deme["deme"]
        lines.append(f'<population id="{name}" spec="ConstantPopulation" popSize="{sizes[name]:g}"/>')
    lines.append("")
    lines.append("<!-- Backward-in-Time Coalescent Lineage Migration Reactions -->")
    for row in mobility:
        origin, destination = row["origin"], row["destination"]
        input_rate = float(row["rate"])
        if rate_convention == "backward-lineage":
            backward_origin, backward_destination, backward_rate = origin, destination, input_rate
            lines.append(f"<!-- Input is backward lineage mobility {origin} -> {destination}. -->")
        else:
            backward_origin, backward_destination = destination, origin
            backward_rate = sizes[origin] * input_rate / sizes[destination]
            lines.append(
                f"<!-- Forward demographic {origin} -> {destination} at {input_rate:g}; "
                f"backward rate = N_{origin} * m / N_{destination}. -->"
            )
        lines.append(
            f'<reaction spec="Reaction" rate="{backward_rate:.12g}"> '
            f'{backward_origin} -> {backward_destination} </reaction>'
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate audited ReMASTER spatial populations and reactions.")
    parser.add_argument("--demes", type=Path, required=True, help="CSV containing deme configurations")
    parser.add_argument("--mobility", type=Path, required=True, help="CSV containing origin,destination,rate")
    parser.add_argument("--mode", choices=["bdm", "coalescent"], default="bdm")
    parser.add_argument("--sample-rate", type=float, default=None, help="Optional per-infected sampling rate (BDM)")
    parser.add_argument(
        "--coalescent-rate-convention",
        choices=["backward-lineage", "forward-demographic"],
        help="Required in coalescent mode; states what the mobility CSV rates mean",
    )
    parser.add_argument("--output", type=Path, help="Optional output XML fragment file")
    parser.add_argument("--audit-only", action="store_true", help="Print diagnostics without XML")
    parser.add_argument("--strict", action="store_true", help="Fail if any audit warning is detected")
    parser.add_argument("--quiet", action="store_true", help="Suppress audit warnings")
    args = parser.parse_args()

    if args.mode == "coalescent" and not args.coalescent_rate_convention:
        parser.error("--coalescent-rate-convention is required when --mode coalescent")
    if args.mode == "bdm" and args.coalescent_rate_convention:
        parser.error("--coalescent-rate-convention is only valid with --mode coalescent")
    if args.sample_rate is not None:
        try:
            _finite_nonnegative(str(args.sample_rate), field="sample rate", context="command line")
        except ValueError as exc:
            parser.error(str(exc))

    try:
        demes = parse_demes(args.demes, args.mode)
        mobility = parse_mobility(args.mobility, {d["deme"] for d in demes})
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    warnings = audit_mobility(demes, mobility, mode=args.mode)
    if warnings and (not args.quiet or args.audit_only):
        print(f"WARNING: Detected {len(warnings)} mobility network hazard(s):", file=sys.stderr)
        for warning in warnings:
            print(f"  * {warning}", file=sys.stderr)
    if warnings and args.strict:
        print(f"FAIL: strict mode rejected {len(warnings)} warning(s).", file=sys.stderr)
        return 1
    if args.audit_only:
        if not warnings:
            print("PASS: mobility network audit found no configured hazards.")
        return 0

    if args.mode == "bdm":
        fragment = generate_bdm(demes, mobility, sample_rate=args.sample_rate)
    else:
        fragment = generate_coalescent(demes, mobility, args.coalescent_rate_convention)
    if args.output:
        args.output.write_text(fragment + "\n", encoding="utf-8")
        print(f"Wrote XML fragment to {args.output}")
    else:
        print(fragment)
    return 0


if __name__ == "__main__":
    sys.exit(main())
