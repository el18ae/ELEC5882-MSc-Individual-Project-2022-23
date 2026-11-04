# 09 — Roadmap

Two parallel tracks: **K** (kernel: our own engine) and **N** (knowledge,
checks, agent). N depends on K for anything that reads a design. Effort
assumes Claude Code implements and you give ~1–2 hours a week of decisions
and sampling. Durations are calendar and conservative.

| Phase | Track | Goal | Duration | Exit criteria |
|---|---|---|---|---|
| **0 Foundations** | both | Answer docs/11; repo, `uv`, Postgres+pgvector (Docker Desktop); inventory and licence-tag the KiCad and Altium install trees (docs/13 § 13.4 step 1); clone KiCad and Freerouting source; pick and vet the first 10 golden boards from the web (docs/08 § 8.1) | 1–2 weeks | Inventory complete; 10 golden boards downloaded with licence records; CI runs; ADRs accepted |
| **K1 Design notes** | K | Read KiCad source, Altium docs/examples, Freerouting and the literature; write `design-notes/` per subsystem, each ending in a port/adapt/design-fresh verdict (docs/13 § 13.6); draft `schemas/rule-taxonomy.md` (union of DRU + Altium rule types) | 3–4 weeks | Every kernel module has a note with a verdict; taxonomy covers 100% of rule types in golden boards |
| **K2 Model + I/O** | K | `model`, `geom`, `io.kicad` (read/write), `io.altium` (read), `io.eagle`; round-trip tests against `oracle.*` | 5–7 weeks | All golden boards round-trip with identical `kicad-cli` ERC/DRC and netlist; Altium examples parse with 100% component/net recovery |
| **K3 ERC/DRC/netlist/render** | K | `erc`, `drc` with the rule language, `netlist`, `render`, `export` (Gerber X3, IPC-2581, BOM/pos) | 5–7 weeks | Kernel DRC matches oracle DRC on golden boards (same violations ±0); Gerbers accepted by a fab's online checker |
| **1 Deterministic review (T0–T1)** | N | Finding schema, fab capability YAML + rule sets, report generator, MCP server with `kernel.*`, `check.*` T0–T1 | 2–3 weeks (overlaps K3) | Workflow A T0–T1 on every golden board; fault-injection recall 100% for T1 classes |
| **2 Parts and sourcing** | N | Parts store seeded from KiCad + Altium libraries; distributor clients; `source.*`; nightly lifecycle job | 3–4 weeks (overlaps K2) | Golden BOM lines resolved or stubbed; availability verified on 50-part sample |
| **3 Documents and rules** | N | Datasheet/app-note crawler; Docling; extraction via Batches; rules store; T2 engine; first 300 rules spanning all domains equally | 5–7 weeks | T2 recall ≥ 0.9; pin tables ≥ 90% match; weekly rule review loop |
| **4 Design corpus and precedent** | N | Corpus crawlers (licence-gated), all through `kernel.io`; motif mining; layout stats; T3; **revision-diff ground truth** (docs/08) | 5–7 weeks | ≥ 2,000 designs; motifs for top-200 parts across domains; T3 judged useful ≥ 70% |
| **5 Model review and evals** | N | Workflow A end-to-end with Opus review + Fable critic; eval harness; trust rungs 1–2 | 3–4 weeks | docs/08 metrics met |
| **6 Generation** | N+K | JSON IR, `gen.block`, block library, Workflow C/D, T4 SPICE, `place` v1 | 7–9 weeks | 50-block eval; one golden board re-created from its requirements, equivalence ≥ 0.9 |
| **7 Layout** | K+N | `place` v2, own `route` (Freerouting via DSN as interim), return-path and impedance checks, Workflow E | 10–14 weeks | Layout eval; blind review on 5 boards |
| **8 Hardening and hand-off** | both | Scheduled runs, nightly BOM watch, backups, docs | 2–3 weeks | Every workflow from one prompt |

Total: roughly 14–18 months of calendar time at the stated human
involvement. The review assistant (through Phase 5) lands at roughly 7–8
months because the kernel sits in front of it. If that is too slow, the
one lever that shortens it without changing the end state is temporarily
using `oracle.*` (kicad-cli) as the T1 engine inside Workflow A until K3
replaces it; that is an owner decision (docs/11 Q20).

## Ordering rationale
- Kernel I/O and DRC come first because every later tier reads designs
  through them and because the oracle cross-check gives a hard,
  independent correctness target for free.
- Review before generation: review produces rules and ground truth that
  generation is later checked against.
- All domains equal: ingestion is ordered by data availability and
  licence clarity, not by domain; rule extraction quotas are balanced
  across the domain list in docs/11 Q4.
- Layout last: it has the weakest open knowledge and the highest
  dependence on everything above.

## Phase 0 task list (executed by `prompts/local-kickoff-prompt.md`)
1. Create `pcb-brain` with the layout in docs/01; copy this plan to `docs/plan/`.
2. `uv`, packages, `ruff`, `mypy --strict`, `pytest`, pre-commit; Docker Desktop with Postgres 16 + pgvector; Alembic from `schemas/parts.sql`.
3. Record the machine (CPU, RAM, GPU, disk, Windows build, KiCad and Altium versions and install paths) in docs/10 § Measured.
4. Inventory the KiCad and Altium install trees (docs/13); write `data/inventory/*.jsonl` with licence tags; parse KiCad symbol and footprint libraries into the parts/footprint stores as a first smoke test.
5. Clone KiCad and Freerouting source into `data/raw/source/` (GPL-3 recorded); begin `design-notes/` with the file-format and DRC subsystems.
6. Golden boards: run the selection and vetting procedure in docs/08 § 8.1 for the first 10 boards across the domain list; store under `golden/` with licence and provenance records.
7. Run `oracle.kicad.drc` on each golden board to capture the baseline the kernel must match.
8. Write `config/owner-answers.yaml` from docs/11 with defaults filled.
