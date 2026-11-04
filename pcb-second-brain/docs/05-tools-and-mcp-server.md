# 05 — Tools and the MCP server

One Python MCP server (`pcbbrain-mcp`, official MCP Python SDK) exposes
namespaced tools over the kernel and the knowledge base. All tools return
JSON with `provenance` where applicable. Mutating tools operate on a
git-tracked working copy and return a diff.

## 5.1 `kb.*` — knowledge base
See docs/02 § Query patterns, plus `kb.cite(ids[])` and
`kb.eda.search(query, tool: kicad|altium)` for EDA-tool knowledge (docs/13).

## 5.2 `kernel.*` — our engine

| Tool | Backing |
|---|---|
| `kernel.open(path)` | `io.kicad` / `io.altium` / `io.eagle` → DesignIR handle |
| `kernel.save(handle, format: kicad)` | `io.kicad` writer (Altium delivery = KiCad file opened by Altium's importer, until an Altium writer is justified) |
| `kernel.netlist(handle)` | `netlist` |
| `kernel.erc(handle)`, `kernel.drc(handle, ruleset)` | `erc`, `drc` → Finding[] |
| `kernel.stackup(handle)` | model + impedance table |
| `kernel.edit.*` (add/remove component, set value/part, move, rotate, add track/via/zone, set netclass, set property, set rule) | `model` + `geom`; every edit validated against T0 |
| `kernel.export.gerbers/ipc2581/step/bom/pos` | `export` |
| `kernel.route(handle, nets[], engine: own|freerouting)` | `route` or DSN/SES round trip; always followed by `kernel.drc` |
| `kernel.sim.spice(...)`, `kernel.sim.impedance(...)`, `kernel.sim.em(...)` | `sim` |

## 5.3 `check.*` — suites (docs/07)
`check.run(handle, suites=[T0..T4], fab, profile)` → Finding[]; individual
checkers also exposed (`check.decoupling`, `check.diffpairs`,
`check.return_path`, `check.thermal_pad`, `check.unresolved_parts`,
`check.availability`, `check.lifecycle`, `check.assembler_library`,
`check.silkscreen`, `check.fiducials`, `check.courtyard`,
`check.testpoints`, `check.esd_placement`, ...).

## 5.4 `source.*` — sourcing
`source.offers`, `source.lifecycle`, `source.alternates`,
`source.bom_health` as before. "Worldwide available" (configurable): in
stock ≥ `qty × 3` at ≥ 2 distributors covering ≥ 2 of {Americas, EMEA,
APAC}, lifecycle active, not single-factory where known.

## 5.5 `gen.*` — generation (structured IR → kernel)
`gen.block`, `gen.design`, `gen.placement_intent`, `gen.fix` as in the
original plan, except the compile target is the kernel's DesignIR rather
than KiCad files. atopile remains an optional front-end IR (MIT; its
compiler output is imported through `io.kicad`), JSON netlist IR is the
native path.

## 5.6 `render.*`
`render.board(handle, layers, region, highlight: Finding[])`,
`render.schematic(handle, sheet)` → SVG/PNG from `kernel.render`.

## 5.7 `oracle.*` — cross-checks (evaluation only)
`oracle.kicad.erc/drc/netlist(path)` wrap `kicad-cli` from the local
KiCad install. `oracle.altium.erc/drc/netlist/export(path)` drive the
installed Altium Designer through its scripting interface (docs/13
§ 13.2b). Both are used by `pcbbrain-eval` to diff the kernel's results
against each tool on every golden board. Not available to workflows.

## 5.8 Cross-cutting
Run ledger per tool call; citation ids carried into findings; working
copies for every mutating workflow.
