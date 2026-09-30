# ReMASTER Model Debugging & Troubleshooting Guide

This guide maps common runtime failures, simulation hangs, and scientific misspecifications to their exact XML fixes.

## 1. Unexpected ancestry or multi-root trees

* **Symptom**: Runtime throws `Tree contains multiple roots` or downstream phylogenetic analysis fails because trees have multiple root lineages.
* **Diagnosis**: Apply ReMASTER's rule to each product: first same-type reactant, otherwise first reactant, otherwise orphan. `S + I -> 2I` already maps both products to `I` and is not, by itself, an ancestry bug. By contrast, in `S_B + I_A -> I_B + I_A`, `I_B` has no same-type reactant and falls back to `S_B`.
* **Fix**: Add labels where the default differs from the intended donor, for example `S_B:a + I_A:b -> I_B:b + I_A:b`. Also check multiple initial lineages, incomplete coalescence before time zero, and acceptance conditions; these can independently produce multiple roots.

## 2. Gillespie Propensity Explosion (Simulation Hang)

* **Symptom**: Simulator consumes 100% CPU and never finishes, without advancing simulation time.
* **Root Cause**: Supercritical branching ($X \to 2X$) combined with an unreachable stopping event (for example, `endsWhen="sample==50"` when no reaction produces `sample`) can make population and event counts explode while the condition remains false.
* **Fix**:
  1. Add the missing sampling or removal reaction (`X -> sample`).
  2. Add a scientifically justified `maxTime` safeguard. This bounds simulated time but is not a substitute for checking event-count feasibility.

## 3. Subcritical rejection loop

* **Symptom**: `nSims` replicates are requested with `mustHave="sample>=50"`, but the simulation spins indefinitely.
* **Root Cause**: A stringent `mustHave` condition can repeatedly reject trajectories when the target event is extremely unlikely. A plain `endsWhen` condition does not itself request rejection; when no reactions can fire, ReMASTER stops the trajectory.
* **Fix**: Estimate feasibility under the actual model, lower or remove the acceptance condition if scientifically justified, or revise parameters/initial seeds. The shortcut $R_0=\beta S_0/(\gamma+\psi)$ and $1-(1/R_0)^{I_0}$ takeoff probability apply only to a simple early-phase homogeneous birth-death approximation, not arbitrary structured models.

## 4. Deterministic Trajectory Floating-Point Overshoot

* **Symptom**: An ODE simulation (`DeterministicTrajectory`) with `endsWhen="X==100"` never terminates at 100 and runs until `maxTime`.
* **Root Cause**: ODE numerical integrators use floating-point states. Continuous state variable $X(t)$ will almost never hit exact integer `100.0000000000`.
* **Fix**: Use inequality conditions: `endsWhen="X>=100"`.

## 5. Corrupted Nexus Tree Logs Across Replicates

* **Symptom**: Multi-replicate tree file crashes tree parsers or fails to parse in FigTree / Tracer.
* **Root Cause**: Replicates have varying numbers of sampled leaves, but `TypedTreeLogger` wrote an initial fixed translation table.
* **Fix**: Add `noLabels="true"` to `<log spec="TypedTreeLogger" tree="@tree" noLabels="true"/>`.

## 6. Higher-Order mass-action rate conversion

* **Symptom**: Simulated rates for reactions with multiple identical reactants are too fast or too slow.
* **Rule**: In continuous ReMASTER, rates are defined **per combination** of available reactants under mass action. Apply factorial scaling only when the source rate is explicitly per ordered sequence, as in MASTER; do not rescale an already per-combination ReMASTER rate.
  - For $k$ identical reactants ($kX \to \dots$), the combinatorial count of available sets is $\binom{X}{k} = \frac{X!}{k!(X-k)!}$.
  - When converting a rate specified per ordered sequence (or macroscopic kinetic constant $k_{\text{ord}}$), the ReMASTER rate must be multiplied by $\frac{1}{k!}$.
  - Examples:
    - $2X \to 3X$: factor is $1/2! = 1/2 = 0.5$.
    - $3X \to 4X$: factor is $1/3! = 1/6 \approx 0.166667$.
    - $2A + B \to C$: factor for $A$ is $1/2! = 0.5$.
