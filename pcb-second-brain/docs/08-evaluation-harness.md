# 08 — Evaluation harness

"As good as me or better" becomes measurable here. Everything is in
`pcbbrain-eval/` and runs in CI on every change to rules, checkers, prompts
or models.

## 8.1 Golden set (sourced and vetted from the web; no owner boards yet)

Selection criteria (all must hold):
1. Licence permits use as test data (CC-BY/CC-BY-SA/CERN-OHL/MIT/TAPR/OCP);
   recorded per board.
2. Complete CAD: schematic + board + libraries (or vendor libraries
   resolvable), not just PDFs. KiCad, Altium or Eagle.
3. Proven: shipped as a product or an official dev board, or built in
   quantity by a community, so the design is known to work.
4. Maintained: has a revision history (git history or revision notes)
   and ideally an issue tracker or errata list.
5. Covers one of the domains in docs/11 Q4; the set is balanced across
   all of them (target 4–6 boards per domain, 40–60 boards total).

Vetting procedure (automated, with a short human look):
- Parse through `kernel.io`; record component/net counts; must parse
  with zero unresolved structural errors.
- Run `oracle.kicad.erc/drc` with the project's own rules; record the
  baseline (a vetted board may have known violations; they are recorded,
  not hidden).
- Resolve BOM to the parts store; record resolve rate.
- Mine **revision-diff ground truth**: for each pair of consecutive
  revisions, diff the DesignIR; changes that touch decouplers, pull-ups,
  terminations, footprints, part swaps, clearances or silkscreen are
  candidate "real findings a human made on this board". Cross-reference
  commit messages, issue trackers and errata. Each confirmed item becomes
  a ground-truth finding with provenance (commit, issue URL). This is the
  substitute for owner-supplied known issues and is often richer.
- Claude review of the board summary to tag domain, interfaces present,
  and any obviously unusual practices (recorded, not judged).

Candidate pools: Adafruit, SparkFun, Olimex, Raspberry Pi (Pico/Pico 2,
CM IO boards), BeagleBoard, Arduino, Seeed, Pine64, Open Compute Project,
TI/ADI/NXP/ST/Espressif/Nordic eval boards with CAD, OSHWA directory,
KiCad demos, Altium's installed Examples (licence per EULA: usable
locally as test data).

- **Synthetic boards**: small projects generated for specific rules.

## 8.2 Fault injection

A library of mutations applied to golden boards, each tagged with the tier
and rule expected to catch it:

remove decoupler · move decoupler 15 mm away · drop bootstrap cap · wrong
feedback divider · I2C without pull-ups · CAN without termination · USB
diff pair split across plane gap · crystal load caps wrong · thermal pad
unconnected · silkscreen over pad · trace below fab minimum · via in pad
where fab disallows · EOL part substituted · part with no stock in EU ·
footprint with pad count mismatch · power pin left floating · ESD diode
far from connector · polarity swap on electrolytic · LDO with ceramic cap
outside stability range · missing fiducials for assembly.

Each mutation is a function `mutate(DesignIR) -> (DesignIR, expected_finding)`.

## 8.3 Metrics

| Metric | Definition | Target by Phase 5 |
|---|---|---|
| Recall per finding class | expected findings caught / injected | ≥ 0.95 T1–T2, ≥ 0.8 T3, ≥ 0.7 T5 |
| Precision per finding class | confirmed findings / reported | ≥ 0.9 overall; ≥ 0.95 for auto-apply classes |
| Hallucination rate | findings/claims with no supporting provenance or whose citation does not support them (checked by critic + sampled human review) | ≤ 1% |
| Ground-truth parity | fraction of revision-diff/errata ground-truth findings caught on golden boards, plus new valid findings you accept on sampled review | ≥ 0.9, rising to 1.0 |
| Oracle agreement | kernel ERC/DRC/netlist vs `kicad-cli` on golden boards (identical violation sets) | 100% |
| Repair convergence | fraction of auto-fix loops that reach zero T1/T2 findings within 5 iterations | ≥ 0.9 |
| Determinism | identical reports on identical inputs (modulo timestamps) | 100% for T0–T4 |
| Cost and time per review | tokens and wall-clock per board, per 100 components | tracked, budgeted |

Generation and layout evals (Phases 6–7) add: functional-equivalence score
(netlist graph similarity to the original golden board after re-creation from
its requirements), DRC-clean rate, impedance error, and blind human review
(you score 5 boards without knowing which were machine-made).

## 8.4 Process

- Weekly 20-minute owner review of a random sample of findings and rules
  (the only recurring human input the plan asks for).
- A finding you dismiss becomes a negative example; a finding you accept
  that the checkers missed becomes a new rule candidate.
- Prompt and rule changes are gated on no metric regression.
