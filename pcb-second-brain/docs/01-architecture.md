# 01 — Architecture

## Three layers

```
┌─────────────────────────────────────────────────────────────────────┐
│ AGENT LAYER                                                         │
│   Claude Code (local) / Claude Agent SDK / Managed Agents           │
│   Workflows: review · source · generate-block · generate-design ·   │
│              layout                                                 │
│   Sub-agents: critic, sourcing, layout-intent, explainer            │
└───────────────▲──────────────────────────────────────▲──────────────┘
                │ MCP tool calls                       │ findings, diffs
┌───────────────┴──────────────────────────────────────┴──────────────┐
│ TOOL LAYER  (one MCP server: `pcbbrain-mcp`)                        │
│   Kernel (ours): data model, geometry, KiCad/Altium/Eagle readers,  │
│     KiCad writer, Gerber X3 / IPC-2581 / STEP exporters, ERC, DRC & │
│     rule engine, netlist, renderer, placer, router (Phase 7)        │
│   Analysis: netlist graph, decoupling, diff-pair, return path,      │
│     thermal, stackup/impedance, SPICE (ngspice), EM (openEMS opt.)  │
│   Sourcing: Digi-Key, Mouser, Nexar, LCSC/JLCPCB, element14, TME    │
│   KB: part lookup, doc search, rule lookup, motif lookup, capability│
│   Oracle (tests only): kicad-cli ERC/DRC/netlist for cross-checking │
└───────────────▲──────────────────────────────────────────────────────┘
                │ SQL / vector / graph queries
┌───────────────┴──────────────────────────────────────────────────────┐
│ KNOWLEDGE LAYER  (PostgreSQL + pgvector; files on NVMe)              │
│   parts · documents · design corpus · rules · fab capabilities ·    │
│   EDA knowledge (KiCad/Altium, docs/13) · project memory            │
│   fed by: ingestion pipeline (crawl → parse → extract → index)       │
└──────────────────────────────────────────────────────────────────────┘
```

## The kernel (`pcbbrain-kernel`) — what "our own software" means here

The agent needs an engine, not an editor. The kernel is a Python package
(hot paths in Rust via `pyo3` later if profiling demands) with these
modules, each specified from the design notes produced in docs/13:

| Module | Responsibility | Primary knowledge source |
|---|---|---|
| `model` | Schematic + board data model (DesignIR, schemas/design-graph.md): sheets, symbols, components, pins, nets, netclasses, footprints, pads, tracks, vias, zones, stackup, rules, variants | KiCad board/schematic classes; Altium object model (scripting docs) |
| `geom` | 2D geometry: polygons, arcs, offsetting, boolean ops, spatial index (`shapely` + R-tree) | KiCad geometry/polygon code notes |
| `io.kicad` | Read/write `.kicad_pro/.kicad_sch/.kicad_pcb/.kicad_sym/.kicad_mod` (s-expression) | KiCad file-format docs + source |
| `io.altium` | Read `.SchDoc/.PcbDoc/.SchLib/.PcbLib/.PrjPcb` (OLE compound documents, record streams) | KiCad's Altium importer notes + Altium file structure as observed on installed libraries/examples |
| `io.eagle`, `io.spectra` | Read Eagle XML; DSN/SES for external routers | KiCad importer notes; Freerouting |
| `export` | Gerber X3, Excellon, IPC-2581, drill/pos/BOM, STEP (via `cadquery`/OCC) | Ucamco + IPC-2581 specs; KiCad exporter notes |
| `erc` | Pin-type matrix, unconnected/driven checks, power-flag logic | KiCad ERC + Altium ERC matrix |
| `drc` | Geometric checks (clearance, width, annular, drill, courtyard, mask, silk) over the spatial index; rule language that is a superset of KiCad DRU and Altium rule types | KiCad DRC engine + Altium rule taxonomy (docs/13 § 13.2) |
| `rules` | Declarative rule loader (schemas/rule.schema.json) + check bindings; T2 | — |
| `netlist` | Graph build, net classes, power tree, diff-pair detection, motif extraction | — |
| `render` | SVG/PNG of boards, sheets, regions; layer colours and highlights for findings | KiCad plotter notes |
| `place` | Constraint-driven placer from JSON intents (Phase 6–7) | KiCad/Altium placement tools notes; corpus statistics |
| `route` | Own router (Phase 7): rip-up-and-reroute + push-and-shove, diff pairs, length tuning; DSN export to Freerouting as interim | Freerouting + KiCad PNS notes |
| `sim` | ngspice driver, impedance calculators, thermal estimate | KiCad PCB calculator formulas (documented) |

Design constraints on the kernel:
- Deterministic and headless; every operation is a pure function of the
  design + rules, so runs are reproducible (ADR-0004).
- The rule taxonomy is designed first (`schemas/rule-taxonomy.md`), as the
  union of KiCad DRU and Altium rule types, so existing projects import
  without losing constraints.
- Round-trip tests: read KiCad project → write → `kicad-cli` parses it and
  reports identical ERC/DRC (oracle test). Same for Altium read → KiCad
  write → open in Altium via its KiCad importer (manual, per release).

## Key data flows

**Review flow (Workflow A)**

```
Altium/KiCad project
  → kernel.io reads → DesignIR
  → resolve every BOM line to parts store (unresolved = finding)
  → kernel.erc / kernel.drc (fab rules) → T1
  → rules engine → T2 ; corpus precedent → T3 ; sim → T4
  → context pack: findings + cited KB snippets + kernel.render regions
  → Claude review pass (Opus 5.5) → structured findings with fixes
  → Claude critic pass (Fable 5.1) challenges each finding against sources
  → report ; optional auto-apply via kernel edits → re-run → diff
  → export to KiCad (and Altium via importer) for the owner to open
```

**Generation flow (Workflow C/D)**

```
requirements / "block around U3"
  → retrieve: datasheet sections, app notes, rules for part, motifs
  → Claude emits IR (JSON netlist IR, or atopile module compiled to IR)
  → kernel builds schematic + netlist (deterministic)
  → T1/T2 → repair loop (≤ 5) → T3 precedent diff → human gate (trust ladder)
```

## Repository layout (implementation repo `pcb-brain`)

```
pcb-brain/
  pyproject.toml            uv workspace
  packages/
    pcbbrain-core/          pydantic models, provenance types, config, logging
    pcbbrain-kernel/        the engine (modules above)
    pcbbrain-kb/            DB schema, migrations, query API
    pcbbrain-ingest/        crawlers, PDF parsing, EDA-knowledge extraction (docs/13), extractors
    pcbbrain-checks/        T2–T4 checker suites and rule bindings
    pcbbrain-sourcing/      distributor & lifecycle clients
    pcbbrain-mcp/           MCP server
    pcbbrain-agent/         workflows, prompts, sub-agent configs
    pcbbrain-eval/          golden boards, fault injection, oracle cross-checks, scoring
  design-notes/             subsystem notes written from KiCad/Altium/Freerouting knowledge (docs/13 § 13.4)
  rules/                    YAML rules (guard rails)
  fabs/                     fab capability YAML + generated rule sets
  golden/                   vetted open boards + ground truth (docs/08)
  data/                     NOT in git: raw/, parsed/, index/
  docs/plan/                this folder
```

## Interfaces between layers

- `Finding` (schemas/finding.schema.json), `Rule` (schemas/rule.schema.json),
  `DesignIR` (schemas/design-graph.md), `PartRecord` (schemas/parts.sql)
  are the contracts; everything downstream of `kernel.io` uses DesignIR,
  never raw files.

## Why one MCP server

One server with namespaced tools (`kb.*`, `kernel.*`, `check.*`,
`source.*`, `gen.*`, `render.*`, `oracle.*`) keeps auth, logging and
caching in one place and serves local Claude Code now and Managed Agents
later without change.
