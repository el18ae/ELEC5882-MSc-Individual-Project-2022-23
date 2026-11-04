# 08 — Evaluation harness

"As good as me or better" becomes measurable here. Everything is in
`pcbbrain-eval/` and runs in CI on every change to rules, checkers, prompts
or models.

## 8.1 Golden set

- **Owner boards**: all your past designs (Altium → KiCad imported), with the
  issues you know they had (or found at bring-up) recorded as ground-truth
  findings. This is the most valuable asset in the project; start collecting
  in Phase 0.
- **Open boards**: 30–50 well-regarded open designs across categories
  (MCU breakout, SBC, buck/boost power, USB device, Ethernet, motor driver,
  sensor analog front end, RF module carrier).
- **Synthetic boards**: small KiCad projects generated for specific rules.

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
| Owner-parity score | on owner boards: fraction of your ground-truth findings caught, plus new valid findings you accept | ≥ 1.0 (i.e. finds everything you did) |
| Repair convergence | fraction of auto-fix loops that reach zero T1/T2 findings within 5 iterations | ≥ 0.9 |
| Determinism | identical reports on identical inputs (modulo timestamps) | 100% for T0–T4 |
| Cost and time per review | tokens and wall-clock per board, per 100 components | tracked, budgeted |

Generation and layout evals (Phases 6–7) add: functional-equivalence score
(netlist graph similarity to the original owner board after re-creation from
its requirements), DRC-clean rate, impedance error, and blind human review
(you score 5 boards without knowing which were machine-made).

## 8.4 Process

- Weekly 20-minute owner review of a random sample of findings and rules
  (the only recurring human input the plan asks for).
- A finding you dismiss becomes a negative example; a finding you accept
  that the checkers missed becomes a new rule candidate.
- Prompt and rule changes are gated on no metric regression.
