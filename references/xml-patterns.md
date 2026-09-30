# Executable ReMASTER XML patterns

Use the smallest pattern matching the requested artifact. These are structural patterns, not scientific defaults: replace populations, reactions, conditions, replicate counts, and filenames with the user's model. Keep the required object types, nesting, IDs, references, and logger shape.

## Trajectory only

```xml
<beast version="2.0"
       namespace="beast.base.inference.parameter:beast.base.inference:remaster">
  <run spec="Simulator" nSims="1">
    <simulate id="traj" spec="StochasticTrajectory" maxTime="10">
      <population id="X" spec="RealParameter" value="1"/>
      <samplePopulation id="sample" spec="RealParameter" value="0"/>
      <reaction spec="Reaction" rate="1"> X -> 2X </reaction>
      <reaction spec="Reaction" rate="0.1"> X -> sample </reaction>
    </simulate>
    <logger spec="Logger" fileName="output.traj">
      <log idref="traj"/>
    </logger>
  </run>
</beast>
```

Omit `samplePopulation` when the model has no sampling. Put `endsWhen="..."` and `mustHave="..."` beside `maxTime` on `<simulate>` when required. Do not create child `<endsWhen>` or `<mustHave>` elements, and do not add `initialTime`.

## Birth-death tree

```xml
<beast version="2.0"
       namespace="beast.base.inference.parameter:beast.base.inference:remaster">
  <run spec="Simulator" nSims="1">
    <simulate id="tree" spec="SimulatedTree">
      <trajectory id="traj" spec="StochasticTrajectory"
                  endsWhen="sample==10" mustHave="sample==10">
        <population id="X" spec="RealParameter" value="1"/>
        <samplePopulation id="sample" spec="RealParameter" value="0"/>
        <reaction spec="Reaction" rate="1"> X -> 2X </reaction>
        <reaction spec="Reaction" rate="0.1"> X -> sample </reaction>
      </trajectory>
    </simulate>
    <logger spec="Logger" mode="tree" fileName="output.trees">
      <log spec="TypedTreeLogger" tree="@tree"/>
    </logger>
  </run>
</beast>
```

The `tree` reference targets the `SimulatedTree`, never the nested trajectory. Add a separate trajectory logger only when trajectory output is requested. Here `mustHave` ensures that extinction cannot yield an accepted tree with fewer than ten tips; omit it only when fewer samples are acceptable.

## Coalescent tree

```xml
<beast version="2.0"
       namespace="beast.base.inference.parameter:beast.base.inference:remaster:beast.base.evolution.tree.coalescent">
  <run spec="Simulator" nSims="1">
    <simulate id="tree" spec="SimulatedTree">
      <trajectory id="traj" spec="CoalescentTrajectory">
        <population id="pop" spec="ConstantPopulation" popSize="1.0"/>
        <reaction spec="PunctualReaction" n="10" times="0"> 0 -> pop </reaction>
      </trajectory>
    </simulate>
    <logger spec="Logger" mode="tree" fileName="output.trees">
      <log spec="TypedTreeLogger" tree="@tree"/>
    </logger>
  </run>
</beast>
```

Coalescent models have no `samplePopulation`. Generate tips with punctual `0 -> population` reactions. A BEAST population-size input may be written as an attribute, as above, or as a correctly typed nested parameter when needed.

## Preflight invariants

Before running the validator, confirm:

- the root is `<beast version="2.0">` and its namespace resolves every unqualified `spec`;
- the run is exactly `<run spec="Simulator" nSims="...">`;
- every birth-death `population` and `samplePopulation` has an explicit concrete `spec` and initial `value`;
- conditions are trajectory attributes and only supported inputs are present;
- tree output has a `SimulatedTree` wrapper and the nested tree logger points to its ID;
- each requested output has an outer BEAST `Logger` and a nested `<log>` object or reference.
