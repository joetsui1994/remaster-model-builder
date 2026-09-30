# Advanced model-building test cases

These cases introduce computational phylodynamics and scientific modeling challenges beyond basic XML boilerplate. They evaluate:
1. Multi-patch spatial metapopulation modeling from tabular data (`demes.csv`, `mobility.csv`).
2. Conversion between forward demographic mobility and backward lineage rates.
3. Model debugging and repair of broken or hanging simulations.
4. Higher-order mass-action kinetics scaling ($1/k!$).
5. Commuter cross-patch transmission ancestry tracking.

The executable oracle is `scripts/verify_remaster_xml.py`. Run self-test via:
```bash
python3 scripts/verify_remaster_xml.py --self-test
```

---

## T13 — Spatial Metapopulation from Tabular Files

**Prompt**

> Build `model.xml` for a stochastic spatial metapopulation epidemic lasting 30 time units using the configuration in `tests/advanced/fixtures/demes.csv` and `tests/advanced/fixtures/mobility.csv`. Initial susceptible, infected, recovery, and transmission parameters for each deme (`North`, `Central`, `South`) must match `demes.csv`. Individuals in `I` sample with removal into `sample` at per-capita rate 0.05. Between-deme movements for both $S$ and $I$ must follow `mobility.csv`. Write the trajectory to `spatial.traj` and a population-type-annotated tree to `spatial.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T13 --xml model.xml
```

The oracle requires 9 populations ($S, I, R$ per deme), 1 sample population, 3 within-deme transmission reactions, 3 recoveries, 3 samplings, 8 between-deme migrations, and loggers for trajectory and typed tree. It accepts either ReMASTER's correct same-type ancestry default or equivalent explicit labels for within-deme transmission.

---

## T14 — Structured Coalescent Tree with Mobility Matrix

**Prompt**

> Build `model.xml` for one typed ReMASTER coalescent tree using the constant effective population sizes in `tests/advanced/fixtures/demes_coalescent.csv`. The rates in `tests/advanced/fixtures/mobility.csv` are forward per-capita demographic movement rates. Under the fixed-size migration convention, convert each forward rate $m_{ij}$ to the backward lineage transition $j \to i$ at rate $N_i m_{ij}/N_j$. Sample 10 lineages in each deme at time 0. Write a typed tree to `coalescent_spatial.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T14 --xml model.xml
```

The oracle requires `CoalescentTrajectory`, 3 `ConstantPopulation` elements, and the four converted backward rates: `Central -> North` at 0.01, `North -> Central` at 0.04, `South -> Central` at 0.02, and `Central -> South` at 0.01125. It also requires 3 punctual sample reactions generating 10 lineages at time 0 and a `TypedTreeLogger`.

---

## T15 — Model Repair: Tree Ancestry Hijacking

**Prompt**

> The model file `tests/advanced/fixtures/t15_broken_ancestry.xml` attempts to simulate a birth-death-sampling tree, but tree reconstruction produces trees with multiple roots because of how transmission is declared. Diagnose and repair `model.xml` so that the tree simulation produces a single-root monophyletic tree descended from the infecting lineage. Preserve all other rates and stopping conditions. Return only the repaired XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T15 --xml model.xml
```

The broken reaction is `S + I -> E + I`. Because `E` has no same-type reactant, ReMASTER otherwise assigns it to the first reactant, `S`. The oracle requires the repair `S:a + I:b -> E:b + I:b`, so sampled `E` lineages descend from the infectious donor. This deliberately tests a case where labels change the default; ordinary `S + I -> 2I` does not require labels.

---

## T16 — Model Repair: Gillespie Simulation Hang

**Prompt**

> The model file `tests/advanced/fixtures/t16_broken_gillespie.xml` consumes 100% CPU and hangs indefinitely during execution without completing. Diagnose the root cause and provide a repaired `model.xml` that terminates safely and saves `gillespie.traj`. Add a protective `maxTime="20"` safeguard and ensure `sample==50` is reachable by adding case sampling at rate 0.5. Return only the repaired XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T16 --xml model.xml
```

The oracle requires `endsWhen="sample==50"`, a positive `maxTime` safeguard, and a reaction generating `sample` (e.g. `X -> sample` at rate 0.5).

---

## T17 — Higher-Order Mass-Action Rate Conversion

**Prompt**

> An older chemical kinetics simulation specifies the trimolecular reaction `3X -> 4X` with rate constant 0.6 per ordered triplet. Create `model.xml` for the equivalent ReMASTER stochastic trajectory under ReMASTER's per-combination continuous convention. Start with `X=100`, stop when `X==120` or at time 1, whichever happens first, run 1 replicate, and save `trimolecular.traj`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T17 --xml model.xml
```

The oracle requires ReMASTER rate `0.1` ($0.6 \times \frac{1}{3!} = 0.6 / 6 = 0.1$), rejecting `0.6`.

---

## T18 — Commuter Metapopulation Cross-Patch Transmission

**Prompt**

> Build `model.xml` for a two-deme stochastic tree simulation lasting until 10 samples are collected. Deme A has `S_A=500, I_A=1`. Deme B has `S_B=500, I_B=0`. An infectious individual in A commutes to B and infects a susceptible individual in B at rate 0.0002 without relocating (`S_B + I_A -> I_B + I_A`). Within Deme A, transmission is `S_A + I_A -> 2I_A` at rate 0.0003. Sampling removes infected individuals at rate 0.05 into `sample`. The commuter reaction must explicitly preserve lineage ancestry so that the newly infected individual in B is a child of the visiting $I_A$ donor. Save the tree to `commuter.trees`. Return only the XML.

**Verify**

```bash
python3 scripts/verify_remaster_xml.py --case T18 --xml model.xml
```

The oracle requires `mustHave="sample==10"` so extinction cannot yield a tree with fewer than the requested 10 samples, `S_B:a + I_A:b -> I_B:b + I_A:b`, both sampling reactions, and `TypedTreeLogger`. For within-deme `S_A + I_A -> 2I_A`, it accepts the correct same-type default and equivalent labelled or expanded-product spellings.
