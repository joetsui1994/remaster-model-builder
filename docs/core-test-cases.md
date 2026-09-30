# Core model-building test cases

These tests evaluate whether the skill produces a scientifically faithful ReMASTER model, not whether its prose resembles a reference answer. Run each prompt with the skill loaded in a clean directory and ask the agent to save its complete answer as `model.xml`. From this skill directory, run the listed command against that file.

The executable oracle is `scripts/verify_remaster_xml.py`. It parses XML, applies general ReMASTER invariants, and checks the case-specific model semantics. `tests/gold/` contains one reviewable reference answer per case.

To prove that the test oracles themselves work:

```bash
python3 scripts/verify_remaster_xml.py --self-test
```

This must accept all 18 gold fixtures, reject 18 deliberately mutated answers, and accept
10 structurally different but semantically valid alternatives.

## T01 — Conditioned serial birth-death-sampling simulation

**Prompt**

> Create a complete ReMASTER XML file named `model.xml` that simulates 25 exact stochastic reconstructed trees from a linear birth-death-sampling process. Start with `X=2` and `sample=0`; use `X -> 2X` at rate 1.6, `X -> 0` at rate 0.5, and sampling with removal `X -> sample` at rate 0.2. Stop each trajectory when `sample==12`, and reject trajectories that terminate before reaching 12 samples. Write the trajectory to `serial-bds.traj` and a typed tree to `serial-bds.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T01 --xml model.xml
```

The check requires `nSims=25`, a `SimulatedTree` around a `StochasticTrajectory`, the exact initial values and three reactions, both `endsWhen="sample==12"` and `mustHave="sample==12"`, a trajectory logger, and a `TypedTreeLogger`. The acceptance condition prevents extinction-shortened simulations from being accepted as twelve-tip trees. Ground truth: ReMASTER manual §§4.2.3 and 4.3.

## T02 — Piecewise-constant transmission rate

**Prompt**

> Build `model.xml` for 20 stochastic ReMASTER trajectories of an SIR epidemic over 60 days. Initial values are `S=999`, `I=1`, and `R=0`. Infection is `S + I -> 2I`, with a per-pair rate of 0.0004 before day 20 and 0.0001 from day 20 onward. Recovery is `I -> R` at rate 0.2. Save trajectories as `sir.traj`. No tree is required. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T02 --xml model.xml
```

The oracle requires `rate="0.0004 0.0001" changeTimes="20"` on the infection reaction, not two competing reactions or a population-wide rate. Ground truth: ReMASTER manual §4.1.2.1, which defines continuous rates per reactant combination and piecewise rates using `rate` plus `changeTimes`.

## T03 — Deterministic linear trajectory

**Prompt**

> Create `model.xml` for one deterministic ReMASTER trajectory from time 0 to 10. Start with `X=1000`; use birth `X -> 2X` at rate 1.2 and death `X -> 0` at rate 1. Save the trajectory to `deterministic.traj`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T03 --xml model.xml
```

The oracle requires `DeterministicTrajectory`, `maxTime="10"`, `nSims="1"`, and the exact linear reactions. Ground truth: ReMASTER manual §4.2.2.

## T04 — Fixed-count punctual sampling with removal

**Prompt**

> Build `model.xml` for one stochastic ReMASTER trajectory lasting 5 time units. Start with `X=100` and a sample population `sample=0`. At times 1, 2, and 5, sample and remove exactly 5, 7, and 10 individuals respectively. Save the trajectory to `fixed-samples.traj`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T04 --xml model.xml
```

The oracle accepts either one vectorized `PunctualReaction` with `n="5 7 10"` and `times="1 2 5"`, or three scalar punctual reactions pairing those counts and times. It rejects `p`, continuous sampling, and non-removing sampling. Ground truth: ReMASTER manual §4.1.2.2; both accepted encodings execute successfully in ReMASTER 2.7.4.

## T05 — Probabilistic punctual sampling without removal

**Prompt**

> Build `model.xml` for one stochastic ReMASTER trajectory lasting 5 time units. Start with `X=100` and sample population `sample=0`. At each of times 1, 2, and 5, independently sample every current `X` individual with probability 0.5 without removing it from `X`. Save the trajectory to `probability-samples.traj`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T05 --xml model.xml
```

The oracle requires `PunctualReaction p="0.5" times="1 2 5"` and `X -> X + sample`. It explicitly rejects `X -> sample`. Ground truth: ReMASTER manual §§4.1.1 and 4.1.2.2.

## T06 — Explicit tree ancestry labels

**Prompt**

> Create `model.xml` for a stochastic ReMASTER tree simulation with `S=999`, `I=1`, and sample population `sample=0`. Transmission occurs at rate 0.1 and must explicitly make the original `I` reactant the parent of both resulting `I` children, while the `S` reactant has no child. Sampling removes `I` into `sample` at rate 0.05. Stop at 10 samples and write a typed tree to `labelled.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T06 --xml model.xml
```

The oracle requires the explicit labelled transmission reaction `S:a + I:b -> 2I:b`, making both infectious products descendants of the infecting lineage. It rejects the alternative `2I:a`, which can leave samples descended from several initially susceptible individuals and prevent a single-root reconstructed tree. Ground truth: the label semantics in ReMASTER manual §4.1.2.3, confirmed by executing the model in ReMASTER 2.7.4.

## T07 — Termination versus acceptance conditions

**Prompt**

> Build `model.xml` for 10 stochastic ReMASTER trajectories with `X=1` and sample population `sample=0`. Use birth `X -> 2X` at rate 2, death `X -> 0` at rate 0.5, and non-removing sampling `X -> X + sample` at rate 0.2. Stop when `sample==20`, impose a maximum time of 50, and accept a finished trajectory only when `X>=1`. Save trajectories to `conditioned.traj`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T07 --xml model.xml
```

The oracle requires `endsWhen="sample==20"` and `mustHave="X>=1"` as distinct trajectory attributes, plus `maxTime="50"`. Ground truth: ReMASTER manual §4.2.3.

## T08 — Typed tree from a two-type model

**Prompt**

> Create `model.xml` for one stochastic ReMASTER tree simulation lasting 10 time units. Start with `X=1`, `Y=0`, and sample population `samp=0`. Use `X -> 2X` at rate 1.4, `X -> Y` at rate 0.1, `X -> 0` at rate 0.5, and `Y -> samp` at rate 0.5. Write a population-type-annotated tree to `typed.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T08 --xml model.xml
```

The oracle requires the four exact reactions and a `TypedTreeLogger` referencing `@tree`; an ordinary tree logger is not accepted. Ground truth: ReMASTER manual §4.3.2.

## T09 — Variable numbers of leaves across replicates

**Prompt**

> Build `model.xml` for 10 stochastic ReMASTER tree simulations lasting 10 time units. Start with `X=1` and sample population `sample=0`. Use `X -> 2X` at rate 1.4, `X -> 0` at rate 0.5, and non-removing sampling `X -> X + sample` at rate 0.5. Accept only trajectories with at least two samples. Because replicate trees can have different leaf counts, write them safely to one file named `variable.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T09 --xml model.xml
```

The oracle requires `TypedTreeLogger noLabels="true"`, non-removing sampling, and a tree logger. Ground truth: ReMASTER manual §4.3.4, which prescribes `noLabels="true"` for tree files containing varying leaf counts.

## T10 — Constant-size coalescent tree

**Prompt**

> Create `model.xml` for one ReMASTER coalescent tree. Use one population `pop` with constant effective population size 1.0. Generate exactly 10 sampled lineages at the present, time 0, and write the tree to `constant-coalescent.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T10 --xml model.xml
```

The oracle requires a `CoalescentTrajectory`, `ConstantPopulation popSize="1.0"`, no `samplePopulation`, the punctual backward-time reaction `0 -> pop` with `n="10" times="0"`, and the BEAST coalescent namespace. Ground truth: ReMASTER manual §§5.1 and 5.3.

## T11 — Structured coalescent with migration

**Prompt**

> Build `model.xml` for one typed ReMASTER coalescent tree. Population `C` has constant effective size 10.0. Population `E` has present effective size 100.0 and `ExponentialGrowth` growth rate 1. Backward in time, lineages transition `C -> E` at rate 0.5. Generate 50 lineages in each population at time 0. Write a typed tree to `structured-coalescent.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T11 --xml model.xml
```

The oracle requires the two correct BEAST population functions, backward-time `C -> E`, two `0 -> population` punctual reactions, the coalescent namespace, and `TypedTreeLogger`. Ground truth: ReMASTER manual §§5.2–5.3.

## T12 — Convert an old MASTER repeated-reactant rate

**Prompt**

> An old MASTER model contains the reaction `2X -> 3X` with rate 0.8 per ordered reactant sequence. Create `model.xml` for the equivalent ReMASTER stochastic trajectory, preserving the same dynamics under ReMASTER's per-combination convention. Start with `X=100`, stop when `X==110` or at time 1, whichever happens first, run one replicate, and save `converted.traj`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T12 --xml model.xml
```

The oracle requires a ReMASTER rate of `0.4`, because the old rate must be multiplied by `1/2!` for the two identical `X` reactants. It rejects retaining `0.8`. Ground truth: the MASTER compatibility note in ReMASTER manual §4.1.2.1.

## Running the suite

Run the complete oracle self-test from the repository root:

```bash
python3 scripts/verify_remaster_xml.py --self-test
```

The oracle checks model semantics and requested invariants rather than requiring a fixed
stochastic trajectory.
