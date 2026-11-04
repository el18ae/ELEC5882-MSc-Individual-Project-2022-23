# Kickoff prompt for local Claude Code (Phase 0)

Paste the following into a Claude Code session opened in an empty directory
on your local machine (the plan folder should be copied there first, or
reachable by path).

---

You are starting the PCB Second Brain project. The complete plan is in
`./pcb-second-brain/` (copied from my GitHub repo). Read `README.md`,
`CLAUDE.md`, every file in `docs/`, `decisions/` and `schemas/` before doing
anything else.

Then execute Phase 0 exactly as listed in `docs/09-roadmap.md` § "Phase 0
task list":

1. Create the implementation repo `pcb-brain` with the layout in
   `docs/01-architecture.md` § Repository layout, move the plan into
   `pcb-brain/docs/plan/`, and keep it current as you learn things.
2. Set up `uv`, packages, `ruff`, `mypy --strict`, `pytest`, pre-commit, and
   a `docker compose` with PostgreSQL 16 + pgvector. Apply
   `schemas/parts.sql` through an Alembic migration.
3. Detect and record my machine (CPU, RAM, GPU, free disk, OS) in
   `docs/10-stack-and-hardware.md` under "Measured".
4. Install and verify: KiCad 10 with `kicad-cli`, enable the IPC API and
   confirm a `kipy` connection; atopile; Freerouting; ngspice; Docling.
   Record exact versions. If anything differs from the plan's assumptions,
   fix the plan in the same commit and log it in `docs/12-sources.md`
   § Corrections.
5. Ask me once for the paths to my past boards (Altium or KiCad) and import
   them into `golden/` as KiCad projects and DesignIR snapshots. If I have
   not provided them yet, proceed with three open-source boards (one
   Adafruit, one SparkFun, one Raspberry Pi Pico) as placeholders.
6. Write the fab capability YAML for my main fab (ask once; default
   JLCPCB 4-layer), generate the `.kicad_dru`, run `kicad-cli pcb drc` on a
   golden board, emit Findings per `schemas/finding.schema.json`, and render
   a Markdown report.
7. Copy the questions in `docs/11-open-questions-for-owner.md` into
   `config/owner-answers.yaml` with the defaults filled in; I will edit it.

Work autonomously. Do not stop to ask questions that have a default in
docs/11; apply the default, note it, continue. Stop only for: credentials,
paths to my private boards, or anything destructive. When Phase 0's exit
criteria in docs/09 are met and tests pass, write `docs/status.md` with
what was built, what deviated from the plan and why, and the exact first
command for Phase 1.
