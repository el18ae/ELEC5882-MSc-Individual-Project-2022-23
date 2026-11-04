# Session log

One entry per working session. Newest at the bottom. The next session
starts by reading this file, then `README.md`, then whatever the work
queue points at. Keep entries short: decisions, artefacts, open items.

## How sessions work (cloud phase, before the local build)

- Work happens on branch `pcb-second-brain-plan` of
  `el18ae/ELEC5882-MSc-Individual-Project-2022-23`, inside
  `pcb-second-brain/`. Every session ends with a commit, a push, and a new
  entry here plus an updated `prompts/next-session-prompt.md`.
- Cloud sessions cannot see the owner's Windows machine, KiCad or Altium
  installs. Anything that needs them is queued for the local phase.
- Cloud sessions can read the public web, clone public repositories
  (network policy permitting) and write files here. That covers most of
  Phase 0 preparation and a large part of K1 (design notes from KiCad
  source) and Phase 3 preparation (rules from public app notes).

## Standing decisions (do not re-litigate; change only by owner instruction)

| # | Decision | Where recorded |
|---|---|---|
| D1 | Own engine (kernel); KiCad and Altium are knowledge sources, interchange formats and test oracles; no GUI editor | ADR-0001 |
| D2 | Model emits structured IR only; the kernel writes files | ADR-0002 |
| D3 | Deterministic checkers are the guard rails; model interprets; trust ladder per capability | ADR-0003, docs/07 |
| D4 | Provenance on every record, rule and finding | ADR-0004 |
| D5 | Personal project: never sold, distributed or open-sourced; may be shown to Altium one day subject to legality | ADR-0001, owner-answers.yaml |
| D6 | Never decompile, disassemble or instrument Altium binaries; use scripting-based black-box observation instead | docs/13 § 13.2b, § 13.5 |
| D7 | KiCad code may be ported only on a "port" verdict from a per-module best-known-method review, with provenance tags | docs/13 § 13.6, CLAUDE.md |
| D8 | All domains equal; quotas balanced across the nine-domain list | docs/11 Q4, owner-answers.yaml |
| D9 | No owner boards; golden boards are web-sourced and vetted; revision-diff ground truth | docs/08 § 8.1 |
| D10 | Platform is Windows 11 with KiCad and Altium installed; their directories may be mined | docs/10, docs/13 |
| D11 | Models: Opus 5.5 agent, Fable 5.1 critic/hard generation, Sonnet/Haiku 5.5 extraction via Batches | CLAUDE.md |

## Open items (owner)

Unanswered questions still on defaults: Q1 (existing second-brain stack),
Q5 (fabs), Q6 (cloud embeddings), Q7 (regions/regulatory), Q8 (API
accounts), Q9 (owned standards), Q10 (budget), Q20 (kicad-cli as interim
T1 engine), Q21 (repo name). See docs/11 and config/owner-answers.yaml.

## Work queue for cloud sessions (in order; each is one or two sessions)

| # | Task | Needs local? | Output |
|---|---|---|---|
| W1 | Golden board selection: build the candidate list across all nine domains with licence, CAD format, revision history and errata availability; attempt downloads where network policy allows; produce the vetting records that do not need the kernel (licence, completeness, revision count) | No | `golden/CANDIDATES.md`, `golden/<board>/PROVENANCE.yaml` |
| W2 | Rule taxonomy: union of KiCad DRU rule language and Altium rule types from public documentation; bidirectional mapping table | No | `schemas/rule-taxonomy.md` |
| W3 | K1 design notes, part 1: clone KiCad source (GitLab) and write notes with verdicts for `io.kicad` (file formats), ERC matrix, DRC engine and rule evaluator, zone filler | No | `design-notes/*.md` |
| W4 | K1 design notes, part 2: Altium importer internals (from KiCad source) as the spec for `io.altium`; exporters (Gerber X3, IPC-2581); footprint generators; calculator formulas | No | `design-notes/*.md` |
| W5 | K1 design notes, part 3 (fresh-design modules): survey routing, placement, constraint models, impedance solvers from Freerouting, literature and other open tools; verdicts and chosen approaches | No | `design-notes/*.md` |
| W6 | Kernel data model: pydantic models for DesignIR per `schemas/design-graph.md`, plus the JSON netlist IR spec for generation | No | `schemas/design-ir.py`, `schemas/netlist-ir.schema.json` |
| W7 | Fab capability YAML for JLCPCB, Eurocircuits, PCBWay from public capability pages, with `verified_at` and URLs | No | `fabs/*.yaml` |
| W8 | First rule set: interface rules (USB 2, I2C, SPI, CAN, RS-485, Ethernet 10/100, crystal, LDO/buck basics) from public app notes, with provenance | No | `rules/interfaces/*.yaml` |
| W9 | MCP tool specification: every tool in docs/05 with input/output JSON schemas | No | `schemas/mcp-tools.md` |
| W10 | Local kickoff package: final `prompts/local-kickoff-prompt.md`, install checklist for Windows (uv, Docker Desktop, Java, ngspice, Docling), and the Altium-oracle script skeletons | No | `prompts/`, `pcbbrain-eval/altium-oracle/README.md` |
| L1 | Everything in docs/09 Phase 0 that touches the machine: inventory KiCad and Altium trees, measure hardware, Docker, oracle baselines | **Yes** | local repo |

## Model and budget guidance for cloud sessions

- Owner has about $58 of cloud credits, then moves to Claude Ultra
  subscription limits. Sessions must be economical: shallow clones, no
  re-reading of large files already summarised in the log, batched web
  lookups, one work item per session.
- Model per work item: **Opus 5.5** for W1, W2, W7, W8, W9, W10 (research,
  cataloguing, specification). **Fable 5.1** for W3, W4, W5 (design notes
  with architectural verdicts) and W6 if budget allows. Fable costs about
  2.5× Opus per token and consumes subscription allowance faster.
- Continuity does not depend on the model; it depends on this log and
  the committed files.

## Entries

### Session 1 — 2026-11-04 (cloud)

Decisions D1–D11 taken. Plan written from scratch: README, CLAUDE.md,
docs 00–13, four ADRs, schemas (parts SQL, rule and finding JSON Schema,
design graph), kickoff prompt, env example, owner answers.

Commits on `pcb-second-brain-plan`:
- 7418375 Add PCB Second Brain project plan
- da4f8c3 Revise plan: own engine, all domains equal, web-sourced golden boards, Windows 11
- 4374d95 Record personal-use status; decide against decompiling Altium
- 7439ffd Add port/adapt/design-fresh policy for kernel modules

Owner clarifications received this session: build our own software rather
than wrap KiCad; no domain priority; no golden boards available; Windows
11 with both tools installed; personal never-distributed project;
decompiling Altium judged not worth it and ruled out; porting KiCad code
only where it is the best-known method.

Next: W1 (golden board candidates) unless the owner reorders the queue.
