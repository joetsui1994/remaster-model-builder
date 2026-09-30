# Spatial Metapopulation and Multi-Deme Modeling in ReMASTER

This guide specifies how to model multi-patch epidemics and spatial phylodynamics in ReMASTER.

## 1. Population Structure

For an epidemic across $K$ demes (e.g. `North`, `Central`, `South`):
- Declare discrete population sets per deme: `S_deme`, `I_deme`, `R_deme`.
- Keep population types as `RealParameter` for birth-death models.
- Sampling can be directed to a shared `sample` population or per-deme sample populations (`sample_deme`).

```xml
<population id="S_North" spec="RealParameter" value="1000"/>
<population id="I_North" spec="RealParameter" value="5"/>
<population id="R_North" spec="RealParameter" value="0"/>
<samplePopulation id="sample" spec="RealParameter" value="0"/>
```

## 2. Within-Deme vs. Between-Deme Reactions

### A. Within-Deme Dynamics
- **Transmission**: In `S_North + I_North -> 2I_North`, ReMASTER's same-type default already assigns both products to `I_North`. Explicit labels are optional but can make the intent reviewable:
  `<reaction spec="Reaction" rate="0.0005"> S_North:a + I_North:b -> 2I_North:b </reaction>`
- **Recovery**:
  `<reaction spec="Reaction" rate="0.1"> I_North -> R_North </reaction>`

### B. Between-Deme Relocation / Migration (Forward BDM)
- For individuals physically relocating from origin $i$ to destination $j$ at per-capita rate $m_{ij}$:
  `<reaction spec="Reaction" rate="0.02"> S_North -> S_Central </reaction>`
  `<reaction spec="Reaction" rate="0.02"> I_North -> I_Central </reaction>`
- In the simulated tree, single-reactant single-product reactions of different types create degree-2 **singleton nodes** marking location transitions.
- Use `<log spec="TypedTreeLogger" tree="@tree"/>` to preserve geographic state annotations on branches.

### C. Cross-Deme Commuter Transmission
- When an infectious host in Deme $A$ visits or interacts with Deme $B$ and transmits to a susceptible host without relocating:
  `S_B + I_A -> I_B + I_A`
- **Crucial Ancestry Rule**: By default, ReMASTER assigns $I_B$ to the first reactant ($S_B$), corrupting the transmission tree. You **must** annotate ancestry:
  `<reaction spec="Reaction" rate="..."> S_B:a + I_A:b -> I_B:b + I_A:b </reaction>`

## 3. Backward Structured Coalescent Migration

When modeling structured populations under `CoalescentTrajectory`:
- Time runs **backward from present to past** ($0 \to \infty$).
- State variables represent **ancestral lineages** residing in each deme.
- If the inputs are already backward lineage transition rates, use their directions and values directly.
- If the inputs are forward per-capita demographic rates, reverse the direction and account for population sizes. Under a fixed-size migration model:
  $$q_{j i}=\frac{N_i m_{ij}}{N_j}$$
  where $m_{ij}$ is forward movement $i -> j$ and $q_{ji}$ is backward lineage movement $j -> i$. Reversing the arrow while retaining the same number is only valid under additional balance/equal-size assumptions.
- Lineages are initialized at time 0 using punctual reactions:
  `<reaction spec="PunctualReaction" n="20" times="0"> 0 -> North </reaction>`

## 4. Helper Tool: `generate_spatial_reactions.py`

When provided with `demes.csv` and `mobility.csv`, generate XML blocks using:
```bash
python3 <skill-dir>/scripts/generate_spatial_reactions.py \
    --demes path/to/demes.csv \
    --mobility path/to/mobility.csv \
    --mode bdm
```

For coalescent fragments, also supply `--coalescent-rate-convention backward-lineage` when rows already describe backward lineage transitions, or `--coalescent-rate-convention forward-demographic` to convert forward demographic rates using deme sizes.

## 5. Mobility Network Auditing & Feasibility Safeguards

Before generating XML from user-supplied mobility tables, audit the network for data integrity, flow asymmetry, and compute hazards using `--audit-only`:

```bash
python3 <skill-dir>/scripts/generate_spatial_reactions.py \
    --demes path/to/demes.csv \
    --mobility path/to/mobility.csv \
    --audit-only
```

### A. Referential Integrity & Schema Checks
* **Unknown Demes**: Every `origin` and `destination` in `mobility.csv` must be declared in `demes.csv`. Unrecognized demes indicate typos and will produce undeclared populations.
* **Non-Negativity**: Migration rates must be non-negative numbers ($m_{ij} \ge 0$).
* **Self-Loops**: Migration from a deme to itself ($i \to i$) is redundant and flagged as a warning.

### B. Asymmetrical Movement Flows & Demographic Flux Imbalance
* **Pairwise Asymmetry**:
  - Unidirectional flow ($m_{ij} > 0$ while $m_{ji} = 0$) or severe rate asymmetry ($\max(m_{ij}, m_{ji}) / \min(m_{ij}, m_{ji}) \ge 3$) causes directional population draining.
* **Net Demographic Flux (forward BDM only)**:
  $$\text{NetFlux}_i = \sum_{j \neq i} N_j m_{ji} - N_i \sum_{j \neq i} m_{ij}$$
  - If $|\text{NetFlux}_i| / N_i \ge 10\%$ per time unit, the population will rapidly drain or accumulate in forward BDM simulations, distorting transmission contact rates ($\beta S_i I_i$).
  - Do not interpret this demographic-flux diagnostic as a structured-coalescent result; coalescent behavior depends on the explicitly supplied backward lineage rates and effective population sizes.

### C. Unusually Large Movement Flows & Compute Hazards
* **High mobility relative to the chosen time unit**:
  - A per-capita rate or total outflow near 1 per time unit deserves review, but its meaning depends on the time unit and competing epidemiological rates. Treat the threshold as a heuristic warning, not proof of panmixia.
* **Gillespie Migration Propensity Trap**:
  - Aggregate migration propensity is $A_{\text{mob}} = \sum_{i,j} N_i m_{ij}$.
  - If $A_{\text{mob}} \ge 10^4$ events/time unit, the absolute event load may be expensive. The fraction of work spent on migration cannot be inferred without comparing other propensities.
  - **Mitigation**: Review units and model intent; consider justified rescaling or a deterministic approximation. `maxTime` bounds simulation time but does not guarantee a cheap event simulation when propensities are extreme.
