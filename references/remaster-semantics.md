# ReMASTER semantics used by this skill

This is a focused implementation reference, not a substitute for the ReMASTER manual.

## Authoritative sources

- [ReMASTER manual](https://tgvaughan.github.io/remaster/)
- [Official ReMASTER repository](https://github.com/tgvaughan/remaster)
- [Official `BDserial.xml` example](https://raw.githubusercontent.com/tgvaughan/remaster/master/examples/BDserial.xml)
- [Vaughan 2024 application note](https://doi.org/10.1093/bioinformatics/btae015)

## Birth-death models

- A normal population is a `RealParameter` (or another `Function`) supplied as a `<population>` with a unique `id`; its value is the initial population size.
- Models that use the short spec name `RealParameter` must include `beast.base.inference.parameter` in the root namespace. The standard birth-death namespace is `beast.base.inference.parameter:beast.base.inference:remaster`.
- An explicit short `spec="Logger"` requires `beast.base.inference`; retain it when extending the namespace for a coalescent model.
- A `<samplePopulation>` records samples and is forbidden on the reactant side of a reaction. When it is a simple counter, declare it as a fully initialized `RealParameter`, for example `<samplePopulation id="sample" spec="RealParameter" value="0"/>`; omitting its type causes BEAST to attempt to instantiate the abstract input type.
- A continuous reaction with rate `lambda` has total propensity `lambda * product(combinations(N_i, n_i))`. The `rate` is therefore per available reactant combination, not the total population-wide event rate.
- For a time-varying reaction, `rate="r0 r1 ..." changeTimes="t1 ..."` is piecewise constant: `r0` applies before `t1`, `r1` after `t1`, and so forth.
- A punctual reaction uses either `n` for a fixed number of firings or `p` for an independent probability for each available reactant at the listed times.
- Sampling with removal is represented by `I -> sample`; sampling without removal is `I -> I + sample`.
- `StochasticTrajectory` uses a Gillespie-style algorithm and is an exact realization of the specified continuous-time discrete-state Markov process.
- `DeterministicTrajectory` integrates the corresponding ODE system. It requires `maxTime`, and its `endsWhen` condition should use inequalities rather than equality. It gives the exact expected population dynamics only for linear models. It can be wrapped by `SimulatedTree`; such trees are diffusion-approximation samples rather than exact realizations of the discrete branching process. Fast dynamics require convergence checks using smaller `forwardRelativeStepSize` and `backwardRelativeStepSize` values.
- `endsWhen` stops the trajectory when its predicate becomes true. `mustHave` is an acceptance condition on the finished stochastic trajectory. Both are attributes on the trajectory element, not nested elements. ReMASTER trajectory classes do not accept an `initialTime` input.
- If the requested output must contain exactly $N$ samples, pair `endsWhen="sample==N"` with `mustHave="sample==N"`. Without the acceptance condition, extinction can finish and accept a trajectory with fewer samples than requested.

## Tree ancestry and logging

- Every output uses an outer BEAST `<logger>` (`spec="Logger"` may be omitted). A trajectory file references the trajectory using a nested `<log idref="trajectory-id"/>`; `TrajectoryLogger` is not the outer logger type.
- Reaction notation alone may not fully specify parent-child relationships. By default, each product is assigned to the first reactant of the same population type, otherwise the first reactant, otherwise no parent. Multiplicity does not consume that match: in `S + I -> 2I`, both `I` products default to the `I` reactant.
- Matching `:label` suffixes explicitly connect a reactant parent to product children. In an epidemic transmission tree, `S:a + I:b -> 2I:b` makes both post-transmission infectious lineages descend from the infecting `I` lineage. The manual's contrasting `2I:a` example instead makes them descend from `S` and is useful for explaining the syntax, but is usually not the intended transmission ancestry.
- Wrap any trajectory in `SimulatedTree` when a reconstructed tree is requested, including `CoalescentTrajectory`; the trajectory itself is not accepted as the `tree` input of a tree logger.
- Tree files use an outer `<logger spec="Logger" mode="tree" ...>` containing the requested tree log. `TypedTreeLogger` belongs inside that outer logger and records population types. `removeSingletonNodes="true"` removes degree-two type-change nodes when an untyped downstream tree is needed.
- Use `noLabels="true"` when a single tree log contains replicate trees with different numbers of leaves.

## Coalescent models

- Coalescent time increases into the past. Populations are BEAST `PopulationFunction` objects, such as `ConstantPopulation` and `ExponentialGrowth`.
- Effective population size has units of time and is the inverse pairwise coalescent rate.
- There is no coalescent `samplePopulation`. Generate sampled lineages backward in time with reactions such as `<reaction spec="PunctualReaction" n="10" times="0"> 0 -> pop </reaction>`.
- Reactions act backward in time. In `S -> L`, `L` is the parent and `S` is the child.
- If a table contains forward per-capita demographic migration rates, direction reversal alone is generally insufficient. Under a fixed-size migration model, forward $i -> j$ at rate $m_{ij}$ contributes a backward lineage transition $j -> i$ at rate $N_i m_{ij}/N_j$. A table of backward lineage rates should instead be used directly. State which convention is being used.
- XML using `ConstantPopulation` or `ExponentialGrowth` must make `beast.base.evolution.tree.coalescent` discoverable through the root namespace.

## MASTER conversion

MASTER rates were expressed per ordered sequence of reactants, whereas ReMASTER rates are per combination. To reproduce a MASTER model in ReMASTER, multiply the old rate by `1 / n_i!` for each reactant type `i` that appears `n_i` times. Thus an old MASTER rate of `0.8` for a reaction with `2X` on the left becomes `0.8 / 2! = 0.4`.
