# ADR-0003 — Deterministic checkers are the guard rails; the model interprets

Status: accepted

## Context
A model asked "is this design OK?" will say yes with confidence. ERC, DRC,
rule engines, API lookups and simulators do not have that failure mode.

## Decision
Verification tiers T0–T4 are code. The model (T5) may add findings only with
cited evidence and may dismiss lower-tier findings only with cited
justification; it may never lower the severity of a T1 finding. Trust
(auto-apply, unreviewed generation) is granted per capability only after
the evaluation harness shows precision/recall above thresholds.

## Consequences
- Most engineering effort goes into checkers, rules and the eval harness,
  not prompts. This is intentional.
- New knowledge enters the system preferentially as a rule or a corpus
  statistic, not as prose for the model to remember.
