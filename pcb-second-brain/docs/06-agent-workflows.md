# 06 — Agent workflows

All workflows run as Claude Code sessions locally (Phase 1–6) using the MCP
server, and can later be lifted to Managed Agents with the same tools.
Each workflow is a named prompt + tool allow-list + gates, stored in
`pcbbrain-agent/workflows/<name>.md`.

## Workflow A — Design review (first deliverable, highest value)

Input: a KiCad or Altium project, target fab/process, target markets
(regions, regulatory regime), production quantity, optional focus list.

1. `eda.import` if needed → `eda.netlist` → DesignIR.
2. `check.unresolved_parts` → every BOM line resolved or flagged.
3. `check.run` T1 (ERC/DRC with fab DRU) → T2 (rules for every part and
   interface present) → T3 (precedent: missing motifs, outlier layout stats)
   → T4 (targeted simulation where rules request it: LDO stability with the
   chosen output cap, buck ripple, RC filter corner).
4. `source.bom_health` for availability/lifecycle/assembler class.
5. Context pack assembled: findings + cited KB chunks + renders of regions
   with findings + stackup + your preferences.
6. **Review pass** (Opus 5.5, effort high/xhigh): for each finding, confirm or
   dismiss with reasoning that cites evidence; add findings the checkers
   cannot see (intent-level: "this LDO is fed from the same rail it
   regulates"); propose fixes as structured patches; rank by impact.
7. **Critic pass** (Fable 5.1): adversarially challenge each confirmed
   finding against the cited sources; a finding whose citation does not
   support it is downgraded to "unverified". Also checks for missing
   coverage ("no finding about USB ESD, but there is a USB connector").
8. Report: Markdown for you, JSON for tooling, each finding with severity,
   location, evidence, provenance, fix, and confidence.
9. Optional `gen.fix` on accepted findings → re-run 3–4 → diff.

Exit: a report you would be comfortable sending to a colleague as your own.

## Workflow B — Sourcing and BOM hardening

Input: BOM (from A or standalone), quantity, regions, assembler.
Produces: per-line health, alternates with footprint compatibility proof,
cost at quantity, risk summary (single-source, EOL, long lead), proposed
BOM changes as a patch. Scheduled re-runs on owner boards (nightly) raise
alerts on lifecycle changes.

## Workflow C — Circuit block generation

Input: anchor part (or function: "3.3 V 2 A buck from 12 V"), constraints.
`gen.block` as in docs/05 § 5.5, followed by a mini Workflow A on the output
and a precedent diff ("your block vs the 12 corpus designs using this part:
differences are X, Y"). Output lands in your atopile/KiCad block library
only after the human gate (until trust ladder rung 3).

## Workflow D — Whole design from requirements

Input: requirements document (interfaces, power budget, envelope, environment,
quantity, fab). Steps: requirements → block list → part selection (Workflow
B logic) → `gen.block` per block → stitch → hierarchical schematic → full
Workflow A → BOM → report. Human gates after block list and after schematic.

## Workflow E — Layout

Input: schematic-complete project, fab, stackup, envelope, connector
positions, mechanical keep-outs.
1. `gen.placement_intent` → placer → `check.run` (placement-stage rules:
   decoupler distance, crystal distance, thermal spacing, connector edges).
2. Net classes and constraints from rules (widths, clearances, diff-pair
   gap/impedance from stackup).
3. Power nets and critical nets routed first by intent (pours, wide tracks,
   via stitching) via `eda.edit.*`; remaining nets by `eda.route`.
4. DRC, return-path check, impedance check, length matching, EM sim on
   flagged structures.
5. Repair loop; then Workflow A on the result; renders for your eyes.
Human gate before fab export (until rung 5).

## Sub-agents

- **critic**: reads sources, never tools that mutate; tasked to disprove.
- **sourcing**: only `source.*` + `kb.part.*`.
- **layout-intent**: only `kb.*`, `render.*`, emits JSON intents.
- **explainer**: turns findings into the Markdown report in your voice
  (learned from your past reports in project memory).

## Prompting notes (Opus 5.5 / Fable 5.1 era)

- Give the whole task up front with the gates and the definition of done;
  these models run long autonomous turns well.
- Require tool-backed evidence per claim; the system prompt states that an
  uncited claim is reported as opinion.
- Do not over-prescribe steps; prescribe outputs and constraints.
- Structured outputs for anything another program reads.
