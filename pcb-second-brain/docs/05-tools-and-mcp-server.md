# 05 — Tools and the MCP server

One Python MCP server (`pcbbrain-mcp`, built on the official MCP Python SDK)
exposes namespaced tools. All tools return JSON with a `provenance` field
where applicable. All tools are side-effect-free except those in `eda.edit.*`
and `gen.*`, which operate on a working copy and return a diff.

## 5.1 `kb.*` — knowledge base

See docs/02 § Query patterns. Plus `kb.cite(ids[])` which expands citation
ids into human-readable references for reports.

## 5.2 `eda.*` — EDA operations (KiCad 10)

| Tool | Backing | Notes |
|---|---|---|
| `eda.import(path)` | kicad-cli / KiCad importers | Altium, Eagle → KiCad project in working dir |
| `eda.netlist(project)` | `kicad-cli sch export netlist` | → DesignIR |
| `eda.erc(project)`, `eda.drc(project, dru?)` | `kicad-cli sch erc`, `kicad-cli pcb drc` with chosen fab DRU | JSON reports → Findings |
| `eda.export.gerbers/ipc2581/step/bom/pos` | kicad-cli | For fab submission and 3D checks |
| `eda.stackup(project)` | board file parse | Layers, thicknesses, impedance table if present |
| `eda.edit.*` (place, move, add_track, add_via, set_netclass, set_property) | kipy (KiCad IPC API) against a running KiCad instance, or offline s-expression edits via kiutils when no GUI is available | Board editing is supported by the IPC API in KiCad 10; schematic editing through IPC is expected in KiCad 11, so schematic edits go through the generation path (5.5) or kicad-sch-api style file edits |
| `eda.route(project, nets[], engine)` | Freerouting (DSN/SES round trip); optional commercial engines (Quilter, DeepPCB) behind the same interface if you license them | Always followed by `eda.drc` |
| `eda.sim.spice(netlist, analysis)` | ngspice (bundled with KiCad) | Operating point, transient, AC for analog/power blocks with models from the KB |
| `eda.sim.impedance(stackup, geometry)` | KiCad PCB calculator formulas / `scikit-rf` | Microstrip/stripline/diff-pair |
| `eda.sim.em(...)` | openEMS (optional, Phase 7) | Only for flagged high-speed structures |

## 5.3 `check.*` — checker suites (see docs/07 for tiers)

`check.run(project, suites=[T1,T2,T3,T4], fab, profile)` → Finding[].
Individual checkers callable too (`check.decoupling`, `check.diffpairs`,
`check.return_path`, `check.thermal_pad`, `check.unresolved_parts`,
`check.availability`, `check.lifecycle`, `check.assembler_library`,
`check.silkscreen`, `check.fiducials`, `check.courtyard`, `check.testpoints`,
`check.esd_placement`, `check.connector_keepouts`, ...).

## 5.4 `source.*` — sourcing

| Tool | Behaviour |
|---|---|
| `source.offers(part_id or mpn, regions=[us,eu,uk,apac])` | Fan out to Digi-Key, Mouser, Nexar, element14, TME, LCSC; normalise; cache 24 h; return per-region stock, MOQ, price breaks, lead time |
| `source.lifecycle(part_id)` | Union of distributor flags + manufacturer PCN crawl + Nexar lifecycle (if tier allows); returns status + evidence |
| `source.alternates(part_id, constraints)` | KB alternates + Nexar "similar parts" + parametric search in Digi-Key, filtered by footprint compatibility and the part's rules |
| `source.bom_health(bom)` | Per line: resolved?, lifecycle, worldwide-available?, assembler class, single-source risk, price @ qty |

"Worldwide available" definition (configurable): in stock ≥ `qty × 3` at ≥ 2
distributors covering ≥ 2 of {Americas, EMEA, APAC}, lifecycle active,
not single-sourced to one factory where known.

## 5.5 `gen.*` — generation (structured IR → KiCad)

| Tool | Behaviour |
|---|---|
| `gen.block(part_id, intent, constraints)` | Retrieves datasheet sections, rules and motifs; Claude emits an **atopile module** (or JSON netlist IR); compiler produces a KiCad hierarchical sheet + netlist; T1/T2 checks run; repair loop ≤ 5 iterations; returns sheet + Finding[] + citations |
| `gen.design(requirements_doc)` | Decomposes into blocks (power tree, MCU, interfaces, connectors, protection); calls `gen.block` per block; stitches; runs full checks |
| `gen.placement_intent(project)` | Claude emits placement groups, keep-outs, proximity constraints (decouplers ≤ 2 mm from pins, crystal ≤ 5 mm, etc.) as JSON; a deterministic placer applies them via `eda.edit.*` |
| `gen.fix(finding_id)` | Produces a structured patch for a finding (add component, change value, move, add via) and applies it on a working copy |

Why atopile as the IR: it is a typed, modular hardware description language
with a compiler that emits KiCad projects and keeps layout in sync; MIT
licensed; modules are reusable, so verified blocks accumulate into your own
library. SKiDL (Python) is the fallback for anything atopile cannot express.
The compiler step is what guarantees the model never writes KiCad files
directly.

## 5.6 `render.*`

`render.board(project, layers, region)` and `render.schematic(sheet)` →
PNG/SVG via kicad-cli, for the vision review pass and for reports.

## 5.7 Cross-cutting

- Every tool logs `{tool, args_hash, duration, result_hash}` to a run ledger
  so a review is reproducible.
- Tool results that include KB content include citation ids; the agent
  prompt requires citations to be carried into findings.
- Working copies: every mutating workflow runs in a git-tracked copy of the
  project so the output is a diff you can read and revert.
