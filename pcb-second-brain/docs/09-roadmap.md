# 09 — Roadmap

Effort assumes Claude Code does the implementation and you give ~1–2 hours
per week of decisions and review. Durations are calendar, sequential, and
deliberately conservative; several phases can overlap once Phase 1 is in.

| Phase | Goal | Duration | Exit criteria |
|---|---|---|---|
| **0 Foundations** | Answer docs/11; stand up repo, `uv`, Postgres+pgvector, KiCad 10 + kicad-cli + kipy, atopile, Freerouting, ngspice; collect owner boards into `golden/`; import Altium → KiCad for all of them | 1 week | All owner boards parse to DesignIR; CI runs; ADRs accepted |
| **1 Deterministic review (T0–T1)** | Parsers, DesignIR, fab capability YAML + DRU for your 2–3 fabs, kicad-cli ERC/DRC wrapper, Finding schema, report generator, MCP server skeleton with `eda.*` and `check.*` | 2–3 weeks | Workflow A runs T0–T1 on every owner board; fault-injection recall 100% for T1 classes |
| **2 Parts and sourcing** | Parts store; Digi-Key, Mouser, Nexar, LCSC clients; MPN normalisation; `source.*` tools; BOM health; nightly lifecycle job | 3–4 weeks | Every owner BOM line resolved or explicitly stubbed; worldwide-availability check correct on a hand-verified 50-part sample |
| **3 Documents and rules** | Datasheet/app-note crawler for all resolved parts; Docling parsing; pin-table and rule extraction via Batches; rules store; T2 engine; first 200 rules (your top parts + interfaces + fab non-DRU) | 4–6 weeks | T2 recall ≥ 0.9 on injected faults; pin tables match package counts ≥ 90%; weekly rule review loop running |
| **4 Design corpus and precedent** | Crawlers for Adafruit/SparkFun/Olimex/RPi/Beagle/OCP/vendor refs (licence-gated); import to KiCad; motif mining; layout stats; T3 checks | 4–6 weeks | ≥ 2,000 designs ingested; motifs for your top-100 parts; T3 findings on owner boards judged useful by you ≥ 70% |
| **5 Model review and evals** | Workflow A end-to-end with Opus review + Fable critic; eval harness complete; trust ladder rungs 1–2 earned on owner boards | 3–4 weeks | Metrics in docs/08 met; you accept auto-apply for rung-1 classes |
| **6 Generation** | atopile IR, compiler path, `gen.block`, block library, Workflow C; T4 SPICE for regulators/filters; Workflow D | 6–8 weeks | 50-block eval passes; one owner board re-created from requirements with functional-equivalence ≥ 0.9 |
| **7 Layout** | Placement intents + placer, net classes from rules, Freerouting integration, return-path and impedance checks, Workflow E; optional commercial router | 6–10 weeks | Layout eval; blind review on 5 boards |
| **8 Hardening and hand-off** | Managed Agents or scheduled local runs; nightly BOM watch; Altium export path; docs; backups | 2–3 weeks | You can run every workflow from a single command/prompt |

Total: roughly 8–11 months of calendar time at the stated human involvement;
Phases 1–5 (the review assistant, which already removes most of the second-
eye burden) in roughly 4 months.

## Ordering rationale

- Review before generation: review produces the ground truth and the rules
  that generation is later checked against. It is also useful from week 3.
- Your boards first: they define "as good as me".
- Rules before corpus: rules are high-precision; precedent is high-recall.
  Precision first keeps early trust.
- Layout last: it has the weakest open tooling and the highest dependence on
  everything above.

## Phase 0 task list (what the local kickoff prompt will execute)

1. Create `pcb-brain` repo with the layout in docs/01; copy this folder into `docs/`.
2. `uv init`; add packages; `ruff`, `mypy`, `pytest`; pre-commit.
3. `docker compose` for Postgres 16 + pgvector; alembic migrations from `schemas/parts.sql`.
4. Install KiCad 10 (verify `kicad-cli --version`, enable IPC API, test `kipy` connection), atopile, Freerouting jar, ngspice; record versions in `docs/10`.
5. Ingest owner boards: import Altium → KiCad, parse to DesignIR, snapshot in `golden/`.
6. Write the first fab capability YAML (your main fab) and generate its DRU; run `kicad-cli pcb drc` on an owner board; emit Findings; render report.
7. Open docs/11 questions as a checklist file and fill defaults where you did not answer.
