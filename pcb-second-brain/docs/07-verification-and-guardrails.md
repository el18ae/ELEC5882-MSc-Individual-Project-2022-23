# 07 — Verification and guard rails

## The five tiers

| Tier | Nature | Examples | Who decides | Cost |
|---|---|---|---|---|
| **T0 Resolution** | Existence | Every BOM line → part row with pins and verified footprint; every pin in netlist exists on that part; footprint pad count = pin count | Code | ms |
| **T1 Hard checks** | Deterministic EDA | ERC; DRC against the chosen fab's DRU; courtyard overlaps; unconnected pads; silkscreen on pads; board outline; mask slivers; drill-to-copper | kicad-cli | seconds |
| **T2 Rule checks** | Declarative rules × design graph | Part rules (decoupling, bootstrap, feedback divider range, thermal pad vias), interface rules (USB series R, ESD near connector, I2C pull-ups, CAN termination), fab non-DRU rules (fiducials, panel rails, part orientation), EMC/SI hygiene (return path continuity, stitching vias, loop area), style rules | Rule engine | seconds |
| **T3 Precedent** | Statistical | Motif present in ≥ 80% of corpus designs using this part but absent here; track width for this net class below corpus 10th percentile; decoupler distance above corpus 90th percentile | Code + corpus | seconds |
| **T4 Simulation** | Physics | ngspice operating point/transient for regulators and filters; impedance of controlled nets vs stackup; thermal estimate from copper area and dissipation; openEMS on flagged structures | Simulators | seconds–minutes |
| **T5 Model review** | Judgement | Intent-level errors, cross-block issues, things no rule encodes; *plus* interpretation of T0–T4 output | Opus 5.5 then Fable 5.1 critic | tokens |

T5 may only *add* findings with cited evidence or *dismiss* lower-tier
findings with cited justification. It may not lower the severity of a T1
finding (the fab's rules are the fab's rules).

## Finding contract

See `schemas/finding.schema.json`. Essentials: `id, tier, severity
(blocker/error/warning/info), title, location {refdes[], nets[], sheet,
xy, layer}, evidence (what was observed), provenance[] (doc+page, rule id,
corpus stats, sim run id), fix (structured patch or null), confidence
(high/medium/low/unverified), status (open/accepted/dismissed/fixed)`.

## Guard rails against the model itself

1. Tool allow-lists per workflow and per sub-agent.
2. The model cannot write EDA files; only compilers and `eda.edit.*` can,
   and those only on working copies.
3. Schema-validated IR; a part reference that fails T0 aborts generation.
4. Repair loops are bounded (≤ 5) and must converge on fewer findings, else
   the workflow stops and reports.
5. Citation enforcement: the report generator drops any finding whose
   `provenance` is empty to an "Unverified observations" appendix.
6. Run ledger makes every review reproducible; a nondeterministic
   difference between two runs on the same input is itself a bug report.

## Trust ladder (autonomy earned per capability)

| Rung | Capability unlocked | Precondition (measured by docs/08) |
|---|---|---|
| 0 | Review reports only; nothing applied | — |
| 1 | Auto-apply T1/T2 fixes with unambiguous patches (add missing decoupler, fix silkscreen, swap to in-stock equivalent of the same part) | Precision ≥ 0.95 on those finding classes over ≥ 10 owner boards; zero regressions on re-check |
| 2 | Auto-apply all accepted-class fixes; BOM alternates | Alternates proven footprint- and rule-compatible in 100% of a 100-sample eval |
| 3 | Generated blocks enter the library without review | Block-generation eval: 0 T0–T2 findings and ≤ 1 T3 deviation per block over 50 blocks; critic agreement |
| 4 | Whole designs to schematic without review | Design eval on 10 re-created owner boards: functional-equivalence score (docs/08) ≥ 0.9 |
| 5 | Layout to fab-ready without review | Layout eval: 0 DRC, impedance within ±10%, return-path check clean, your blind review finds no blocker on 5 boards |

Rungs are per capability class, not global; e.g. "passives around an MCU"
can be at rung 3 while "switching power stage" is at rung 1.
