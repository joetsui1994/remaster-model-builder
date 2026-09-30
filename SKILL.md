---
name: remaster-model-builder
description: Build, explain, review, repair, and validate executable ReMASTER XML for birth-death, spatial metapopulation, or coalescent phylodynamic simulations. Use when translating a scientific model description or tabular mobility files into ReMASTER, diagnosing or fixing flawed ReMASTER models, or verifying reactions, spatial migration, lineage semantics, stopping conditions, and loggers.
license: MIT
---

# ReMASTER Model Builder

Translate the scientific model into an explicit simulation contract before generating XML. Treat an XML file that merely parses as insufficient: the reactions, event rates, sampling process, tree ancestry, spatial migration, conditions, and loggers must also match the stated model.

Read [references/remaster-semantics.md](references/remaster-semantics.md) before authoring or reviewing ReMASTER XML. It contains the non-obvious semantics and links to the authoritative manual and examples.

When authoring a new file, also read [references/xml-patterns.md](references/xml-patterns.md) and begin from the matching executable pattern: trajectory only, birth-death tree, or coalescent tree. Preserve its object nesting, required types, references, and logger structure while replacing the scientific content.

When working with spatial metapopulations, mobility matrices, or multi-deme CSV files, read [references/spatial-metapopulations.md](references/spatial-metapopulations.md). Use `scripts/generate_spatial_reactions.py` for validated CSV-to-XML fragments when its input conventions match the user's data.

When debugging or repairing a flawed, hanging, or misspecified ReMASTER XML, consult [references/model-debugging.md](references/model-debugging.md).

## Establish the model contract

Record the following information in a compact table or structured block:

- simulation family: stochastic birth-death, deterministic birth-death approximation, spatial metapopulation, or coalescent;
- time direction and units;
- populations, subpopulations, demes, sample populations, and initial values;
- every reaction, its rate interpretation, and any rate-change times;
- sampling times or rates and whether sampling removes the individual;
- parent-child mapping when a simulated tree is requested;
- between-patch migration rates and forward vs. backward time direction;
- maximum time, stopping conditions, acceptance conditions, and replicate count;
- required trajectory, statistics, typed-tree, untyped-tree, or pruned-tree outputs.

Do not silently decide a scientifically consequential ambiguity. In particular, clarify or explicitly label assumptions about sampling with removal, continuous versus punctual sampling, time direction, rate units, lineage ancestry, and whether mobility values are forward demographic rates or backward lineage rates. When a draft remains useful, choose a conservative assumption, label it, and isolate it so it is easy to change.

## Translate the contract

- Use `StochasticTrajectory` for exact continuous-time discrete-state birth-death simulations and for exact birth-death tree simulation.
- Use `DeterministicTrajectory` for the ODE limit. It requires `maxTime`; equality-based `endsWhen` conditions are inappropriate because its state is floating point. A deterministic trajectory may be wrapped in `SimulatedTree`; the resulting tree is a diffusion-approximation sample, not an exact birth-death realization. Check numerical convergence by varying the forward and backward relative step sizes when dynamics are fast.
- Use `CoalescentTrajectory` when effective population-size functions are fixed and lineages evolve backward in time. Do not create a `samplePopulation` inside a coalescent model; generate tips using punctual reactions such as `0 -> pop`.
- Interpret a continuous `Reaction` rate per available combination of reactants under mass action. Never substitute a population-wide propensity.
- Divide by factorials only when converting a MASTER or explicitly per-ordered-sequence rate. For each reactant type repeated $n_i$ times, multiply the source rate by $1/n_i!$. Do not rescale a rate already specified per combination.
- Represent a piecewise-constant continuous rate with one more `rate` value than `changeTimes` values.
- Use `PunctualReaction n="..."` for fixed event counts and `p="..."` for independent per-combination probabilities. Keep `n` and `p` mutually exclusive.
- Declare every birth-death population and sample population completely as a concrete `Function`, normally a `RealParameter` with an initial `value`. A sample population must not appear on a reactant side. Distinguish removal (`I -> sample`) from non-removal (`I -> I + sample`).
- Apply ReMASTER's ancestry defaults before adding labels: each product is assigned to the first same-type reactant, otherwise the first reactant, otherwise no parent. Thus `S + I -> 2I` already assigns both products to `I`; labels are optional clarification. In `S_B + I_A -> I_B + I_A`, however, `I_B` would fall back to `S_B`, so use `S_B:a + I_A:b -> I_B:b + I_A:b` when the commuter is the intended donor.
- For forward birth-death migration, represent an individual's relocation $i -> j$ directly with per-capita reactions such as `I_i -> I_j` at rate $m_{ij}$.
- Do not blindly reuse forward demographic migration rates in a structured coalescent. Under a fixed-size migration model, a forward per-capita rate $m_{ij}$ implies backward lineage rate $q_{j i}=N_i m_{ij}/N_j$. If the input already contains backward lineage rates, use them directly. Record the convention and population-size basis.
- With spatial CSVs, run the generator's audit first. In coalescent mode, pass an explicit `--coalescent-rate-convention`; do not infer it from column names.
- For unqualified `RealParameter`, use `beast.base.inference.parameter:beast.base.inference:remaster`; append `beast.base.evolution.tree.coalescent` for BEAST coalescent population functions.
- Use `<beast version="2.0" ...>`. Put `maxTime`, `endsWhen`, and `mustHave` on the trajectory as attributes. Do not emit `initialTime`.
- When a stochastic simulation must return exactly $N$ samples, normally use both `endsWhen="sample==N"` and `mustHave="sample==N"`. `endsWhen` stops at the target; `mustHave` rejects trajectories that terminate earlier through extinction. Omit the acceptance condition only when fewer samples are scientifically acceptable, and state that choice.
- Use a BEAST `Logger` as every outer logger. Nest `<log idref="trajectory-id"/>` for trajectories. Do not use `TrajectoryLogger` as the outer logger.
- Wrap every requested tree in `<simulate id="tree" spec="SimulatedTree">`. Put tree output in an outer logger with `mode="tree"` and a nested `TypedTreeLogger` or requested tree logger referencing `tree="@tree"`. Use `noLabels="true"` when replicate leaf counts vary.

Keep runnable XML self-contained unless the user asks for a fragment: a `beast` root, one `Simulator` run with positive `nSims`, one simulated object, and loggers for every requested artifact.

## Validate before reporting completion

Save the XML and run the bundled static validator:

```bash
python3 <skill-directory>/scripts/validate_remaster_xml.py path/to/model.xml
```

For repair, spatial, or conditioning tasks, also run the feasibility audit and review warnings in scientific context:

```bash
python3 <skill-directory>/scripts/check_model_feasibility.py path/to/model.xml
```

If BEAST 2 and ReMASTER are available, execute the XML in an isolated directory and confirm the process exits successfully and every requested output exists and is non-empty:
```bash
beast -seed 17 -overwrite path/to/model.xml
```

Do not claim runtime validation after parsing or static checks alone. For stochastic models, assess invariants or distributions across replicates rather than expecting one fixed trajectory.

## Deliver the result

Provide:

1. the model contract and labelled assumptions;
2. the complete ReMASTER XML;
3. a mapping from scientific events to XML reactions;
4. the exact validation commands and their outcomes;
5. any remaining scientific or runtime limitations.

Follow the user's requested format. If they request XML only, write or return only the XML while still performing the applicable validation.
