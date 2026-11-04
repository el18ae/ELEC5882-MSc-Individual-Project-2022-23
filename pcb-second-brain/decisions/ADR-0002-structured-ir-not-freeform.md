# ADR-0002 — The model emits structured IR; compilers emit EDA files

Status: accepted

## Context
The dominant failure of LLM-generated hardware is invented parts, pins and
footprints, and files that do not parse. Free-form generation has no point
at which reality can be enforced.

## Decision
The model only ever produces schema-validated intermediate representations:
atopile modules or JSON netlist IR for circuits, JSON intents for
placement, structured patches for fixes. Deterministic compilers and the
KiCad IPC API turn these into files. Every part reference in IR must
resolve (T0) before compilation proceeds.

## Consequences
- Generation is slower to build (compiler path first) but every output is
  guaranteed parseable and every part real.
- Capabilities of the IR bound the capabilities of generation; extend the
  IR deliberately, with tests.
