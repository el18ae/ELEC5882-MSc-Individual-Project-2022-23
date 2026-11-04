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
| D12 | Third-party CAD is never committed to this public repository; `golden/<board>/src/` is git-ignored and re-fetched at pinned commits (follows from D5; session 2, owner may override via Q25) | `golden/.gitignore`, docs/11 Q25 |

## Open items (owner)

Unanswered questions still on defaults: Q1 (existing second-brain stack),
Q5 (fabs), Q6 (cloud embeddings), Q7 (regions/regulatory), Q8 (API
accounts), Q9 (owned standards), Q10 (budget), Q20 (kicad-cli as interim
T1 engine), Q21 (repo name), and from session 2: Q23 (ATEX golden boards),
Q24 (extra golden-board licences), Q25 (public repo vs private for CAD).
See docs/11 and config/owner-answers.yaml.

## Work queue for cloud sessions (in order; each is one or two sessions)

| # | Task | Needs local? | Output |
|---|---|---|---|
| W1 ✅ | Golden board selection (done, session 2): build the candidate list across all nine domains with licence, CAD format, revision history and errata availability; attempt downloads where network policy allows; produce the vetting records that do not need the kernel (licence, completeness, revision count) | No | `golden/CANDIDATES.md`, `golden/<board>/PROVENANCE.yaml` |
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

### Session 2 — 2026-11-04 (cloud)

W1 (golden board selection) is complete. Commit 3757364 holds the artefacts.

Decisions:
- D12 (new, follows from D5). The repo is public, so third-party CAD is
  not committed. `golden/<slug>/src/` is git-ignored.
  `golden/fetch-phase0.sh` re-fetches each board at the commit pinned in
  its `PROVENANCE.yaml`.
- Golden-board licence list widened to Apache-2.0, 0BSD, BSD, GPL and
  Solderpad. Non-commercial licences are accepted with a flag
  (docs/11 Q24, default applied).
- ATEX/IS: no open-licensed design with complete CAD exists. The quota
  will be met by synthetic boards from IEC 60079-11 clause-referenced
  rules. Vendor IS reference designs are local reference only, not
  golden (docs/11 Q23, default applied).
- Phase 0 ten: hackrf-one, glasgow, moteus, bms-c1,
  jetson-nano-baseboard, sh-esp32, speeduino-v04, esp32-evb, odrive-v3
  (Altium), eez-psu-h24005 (Eagle). They cover all eight domains that
  have candidates and all three planned readers.

Artefacts:
- `golden/CANDIDATES.md` lists 55 candidates. Each has a quoted
  licence, CAD format, completeness, history, tracker status and
  proven status, ranked within its domain. The file also covers the
  Phase 0 rationale, the ATEX shortfall, leads at blocked hosts and
  rejected repositories.
- `golden/phase0.list`, `golden/fetch-phase0.sh`, `golden/.gitignore`.
- `golden/<slug>/PROVENANCE.yaml` for all ten. Each holds the URL,
  pinned commit, licence quote, fetched_at and fetch commands. Its
  vetting fields stay null until the kernel exists.
- `golden/<slug>/GROUND-TRUTH-CANDIDATES.md` for all ten, all
  unjudged. In total: 26 revision-directory pairs, 600 design-file
  commits whose message matches a fix keyword, and 137 issue-tracker
  leads. Verbatim revision notes come from ESP32-EVB, ODrive v3,
  HackRF and Glasgow, and the BMS C1 test report.
- docs/08 § 8.1 licence list widened; docs/11 has Q23–Q25; docs/12 has
  corrections; `config/owner-answers.yaml` has defaults for Q23–Q25.

Environment facts (cloud):
- Reachable: GitHub (git and `gh api search/*`). Refused:
  `source.mnt.re` and `datasheets.raspberrypi.com`, so MNT Reform and
  the Raspberry Pi IO boards are local-phase leads. Not tried:
  gitlab.com, docs.kicad.org, altium.com.
- About 500 MB of working trees plus 190 MB of history were fetched for
  the ten boards. They were not committed (D12).

Open items:
- Q23, Q24 and Q25 are on their defaults.
- Industrial comms has no 4–20 mA or HART board. RS-485 rests on one
  weakly-proven dongle. Revisit when W8 rules exist, or use synthetic
  boards.
- The ground-truth keyword categories are noisy (Glasgow has 216
  leads). Triage needs the DesignIR diff in Phase 0/K2, not more cloud
  work.
- Vetting steps that need the kernel or the oracle (parse, ERC/DRC
  baseline, BOM resolve) remain for L1/K2.

Next: W2, rule taxonomy (`schemas/rule-taxonomy.md`).

