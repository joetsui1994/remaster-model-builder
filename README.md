# ReMASTER Model Builder

ReMASTER Model Builder is an agent skill for translating scientific phylodynamic model
specifications into executable [ReMASTER](https://github.com/tgvaughan/remaster) XML. It
also reviews existing models for scientific, structural, and computational problems before
they are run.

## What it does

- turns a prose model description into an explicit simulation contract;
- builds stochastic birth-death, deterministic, spatial metapopulation, and coalescent
  simulations;
- checks reaction-rate conventions, sampling with or without removal, stopping and
  acceptance conditions, lineage ancestry, migration direction, and logger structure;
- audits spatial CSV inputs and generates ReMASTER population and reaction fragments;
- detects common feasibility hazards such as unreachable endpoints, unconditioned sample
  targets, explosive event workloads, and contradictory coalescent sampling constructs;
- validates XML statically and, when BEAST 2 with ReMASTER is available, guides runtime
  execution and output inspection.

The skill is designed to make the mapping from scientific assumptions to executable XML
reviewable. A file that merely parses is not considered sufficient.

## Repository layout

```text
.
├── SKILL.md                     Skill entry point and workflow
├── references/                  ReMASTER semantics, XML patterns, debugging, and spatial guidance
├── scripts/
│   ├── validate_remaster_xml.py Structural XML validator
│   ├── check_model_feasibility.py Scientific and computational preflight checks
│   ├── generate_spatial_reactions.py CSV validation, mobility audit, and XML generation
│   └── verify_remaster_xml.py   Executable oracle for the model-neutral test cases
├── tests/                       Unit tests, fixtures, reference XML, and valid alternatives
└── docs/                        Core and advanced test-case specifications
```

## Installation

Clone or copy this repository into your agent's skills directory under the name
`remaster-model-builder`. The directory containing `SKILL.md` is the skill root.

Invoke it as `$remaster-model-builder`, or describe a ReMASTER model-building or review
task and allow normal skill discovery to select it.

## Example requests

- “Build a ReMASTER stochastic SIR tree simulation that stops after 50 samples.”
- “Review this ReMASTER XML for ancestry and sampling errors.”
- “Generate spatial reactions from these deme and mobility CSV files.”
- “Convert these forward demographic migration rates into backward structured-coalescent
  lineage rates.”
- “Diagnose why this Gillespie simulation does not terminate.”

## Command-line helpers

All bundled scripts use the Python standard library and require Python 3.10 or newer.

Validate XML structure:

```bash
python3 scripts/validate_remaster_xml.py path/to/model.xml
```

Run the feasibility audit:

```bash
python3 scripts/check_model_feasibility.py path/to/model.xml
```

Audit spatial inputs without generating XML:

```bash
python3 scripts/generate_spatial_reactions.py \
  --demes path/to/demes.csv \
  --mobility path/to/mobility.csv \
  --audit-only
```

Generate birth-death spatial fragments:

```bash
python3 scripts/generate_spatial_reactions.py \
  --demes path/to/demes.csv \
  --mobility path/to/mobility.csv \
  --mode bdm
```

For structured-coalescent fragments, choose the rate convention explicitly with
`--coalescent-rate-convention backward-lineage` or
`--coalescent-rate-convention forward-demographic`.

## Testing

Run the unit tests:

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
```

Self-test the semantic oracle against the reference, alternative, and deliberately invalid
cases:

```bash
python3 scripts/verify_remaster_xml.py --self-test
```

The repository's continuous-integration workflow runs these static checks on Python 3.10
and a current Python release. It does not install BEAST or ReMASTER. The bundled reference
and alternative XML cases have also been executed separately during development with
BEAST 2.7.7 and ReMASTER 2.7.4; runtime validation of a newly generated model still
requires a local BEAST/ReMASTER installation.

The complete model-neutral case specifications are documented in
[core test cases](docs/core-test-cases.md) and
[advanced test cases](docs/advanced-test-cases.md).

## Runtime dependencies

The static tools have no third-party Python dependencies. Executing generated simulations
requires [BEAST 2](https://www.beast2.org/) with the
[ReMASTER package](https://github.com/tgvaughan/remaster) installed.

## Scope and limitations

- Static validation cannot prove that a model answers the intended scientific question.
- Feasibility warnings are conservative diagnostics, not universal biological thresholds.
- Structured-coalescent migration conversion depends on the stated demographic convention
  and population-size assumptions.
- Stochastic output should be assessed through appropriate invariants or distributions,
  not by expecting one fixed trajectory.

For ReMASTER syntax and semantics, consult the
[official manual](https://tgvaughan.github.io/remaster/) and the references bundled with
this skill.

## License

This skill and its bundled tooling are available under the [MIT License](LICENSE).
ReMASTER and BEAST are separate dependencies, are not bundled here, and remain subject
to their own licenses.
