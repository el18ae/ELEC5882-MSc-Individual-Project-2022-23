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
│   EDA: kicad-cli (ERC/DRC/export/import), kipy IPC (live board),    │
│        atopile/SKiDL compile, Freerouting, ngspice, openEMS (opt.)  │
│   Analysis: netlist graph, decoupling, diff-pair, thermal, stackup  │
│   Sourcing: Digi-Key, Mouser, Nexar, LCSC/JLCPCB, element14, TME    │
│   KB: part lookup, doc search, rule lookup, motif lookup, capability│
│   Render: board/schematic PNG+SVG for vision review                 │
└───────────────▲──────────────────────────────────────────────────────┘
                │ SQL / vector / graph queries
┌───────────────┴──────────────────────────────────────────────────────┐
│ KNOWLEDGE LAYER  (PostgreSQL + pgvector; files on NVMe)              │
│   parts · documents · design corpus · rules · fab capabilities ·    │
│   project memory                                                     │
│   fed by: ingestion pipeline (crawl → parse → extract → index)       │
└──────────────────────────────────────────────────────────────────────┘
```

## Key data flows

**Review flow (Workflow A)**

```
Altium/KiCad project
  → normalise (kicad-cli import if Altium)
  → extract: netlist, BOM, placement, stackup, DRU in use
  → resolve every BOM line to parts store (unresolved = finding)
  → run checker suites (T1 hard, T2 rules, T3 precedent, T4 sim)
  → assemble context pack: findings + cited KB snippets + board renders
  → Claude review pass (Opus 5.5) → structured findings with fixes
  → Claude critic pass (Fable 5.1) challenges each finding against sources
  → report (Markdown + JSON) ; optional auto-apply → re-run checks → diff
```

**Generation flow (Workflow C/D)**

```
requirements / "block around U3"
  → retrieve: datasheet sections, app notes, rules for part, motifs from corpus
  → Claude emits IR (atopile module or JSON netlist IR)
  → compiler: IR → KiCad schematic + netlist (deterministic)
  → T1/T2 checks → repair loop (max N iterations) → T3 precedent diff
  → human gate (until trust ladder says otherwise)
```

## Repository layout (for the implementation repo)

```
pcb-brain/
  pyproject.toml            uv-managed; packages below are workspace members
  packages/
    pcbbrain-core/          shared models (pydantic), provenance types, config
    pcbbrain-kb/            DB schema, migrations (alembic), query API
    pcbbrain-ingest/        crawlers, parsers (PDF, KiCad, Altium-via-KiCad, Eagle), extractors
    pcbbrain-eda/           kicad-cli & kipy wrappers, atopile/SKiDL compile, Freerouting, ngspice
    pcbbrain-checks/        checker suites T1–T4, rule engine, finding emitters
    pcbbrain-sourcing/      distributor & lifecycle API clients, normalisation, caching
    pcbbrain-mcp/           the MCP server exposing everything above
    pcbbrain-agent/         workflow definitions, prompts, sub-agent configs
    pcbbrain-eval/          golden boards, fault injection, scoring, reports
  data/                     NOT in git: raw/, parsed/, index/ (see docs/04)
  rules/                    YAML rule files (in git; these are the guard rails)
  fabs/                     manufacturer capability YAML + generated .kicad_dru
  golden/                   owner's boards + selected open designs used for evals (private submodule or ignored)
  docs/                     this plan, kept current
```

## Interfaces between layers (the contracts that make this testable)

- `Finding` (schemas/finding.schema.json): every checker, rule and LLM pass
  emits the same object. Severity, location (ref des / net / coordinate),
  evidence, provenance, suggested fix as a structured patch.
- `Rule` (schemas/rule.schema.json): declarative rule with a `check` binding to
  a Python callable in `pcbbrain-checks`. Rules are data; checkers are code.
- `DesignIR` (schemas/design-graph.md): the in-memory graph of a design
  (components, pins, nets, sheets, placement, stackup). Everything downstream
  of the parsers reads this, never raw files.
- `PartRecord` (schemas/parts.sql): the canonical part row; symbols, footprints
  and pin tables hang off it.

## Why one MCP server rather than many

Claude Code and the Agent SDK both speak MCP. One server with namespaced tools
(`kb.*`, `eda.*`, `check.*`, `source.*`, `gen.*`, `render.*`) keeps auth,
logging and caching in one place and lets the same tool surface be used from
local Claude Code today and from Managed Agents later without change.
