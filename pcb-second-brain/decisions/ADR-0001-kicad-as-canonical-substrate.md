# ADR-0001 — KiCad is the canonical substrate; Altium is an import/export edge

Status: accepted (proposed by plan; confirm in Phase 0)

## Context
The owner works in both Altium and KiCad. Automation needs one file format
to parse, one CLI to run checks, and one API to edit boards. Altium's
scripting is proprietary, Windows-bound and not headless; KiCad 10 has a
headless CLI, an open s-expression format, a Protobuf IPC API with official
Python bindings, importers for Altium and Eagle, and the whole open
schematic-as-code ecosystem (atopile, SKiDL) targets it.

## Decision
All analysis, checking, generation and layout operate on KiCad projects.
Altium inputs are imported to KiCad once (recorded as `imported_from`).
Altium output, when required, is produced by export at the end of a
workflow and verified by opening in Altium.

## Consequences
- One parser, one DRC engine, one edit API.
- Some Altium-specific data (rooms, certain rule types, variants) may be
  lost on import; the importer's coverage is tested in Phase 0 on owner
  boards and gaps are recorded.
- Clients who require Altium deliverables get them via export; the owner
  should expect a short manual tidy-up in Altium for those until the export
  path is proven.
