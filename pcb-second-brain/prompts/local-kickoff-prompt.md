# Kickoff prompt for local Claude Code (Phase 0, Windows 11)

Paste into a Claude Code session opened in an empty directory on your
machine, with this plan folder copied there.

---

You are starting the PCB Second Brain project on this Windows 11 machine.
The complete plan is in `./pcb-second-brain/`. Read `README.md`,
`CLAUDE.md`, every file in `docs/`, `decisions/` and `schemas/` before
doing anything else. Key decisions: we build our own engine (ADR-0001);
KiCad and Altium are installed here and are to be mined for knowledge
(docs/13), used as interchange formats and as test oracles; all domains
are equal; there are no owner boards, so golden boards come from the web
and are vetted per docs/08 § 8.1.

Execute Phase 0 as listed in `docs/09-roadmap.md` § "Phase 0 task list":

1. Create the implementation repo `pcb-brain` with the layout in
   `docs/01-architecture.md`; move the plan to `pcb-brain/docs/plan/` and
   keep it current as you learn things.
2. Set up `uv`, packages, `ruff`, `mypy --strict`, `pytest`, pre-commit;
   Docker Desktop compose with PostgreSQL 16 + pgvector; Alembic migration
   from `schemas/parts.sql`.
3. Measure and record this machine and the installed tools (CPU, RAM, GPU,
   free disk, Windows build, KiCad version and path, Altium version and
   path, Docker/WSL2, Java) in `docs/10` § Measured.
4. Inventory the KiCad and Altium install trees per docs/13 § 13.4 step 1:
   every file with type, size, sha256 and licence tag, to
   `data/inventory/`. Touch no executables beyond recording their names.
   Then parse the KiCad symbol and footprint libraries into the parts,
   symbol and footprint stores as the first smoke test, and parse one
   Altium `.SchLib` and one `.PcbLib` with `olefile` to confirm the record
   structure matches the design notes you will write.
5. Clone the KiCad and Freerouting source into `data/raw/source/` with
   their GPL-3 licences recorded. Start `design-notes/` with two notes:
   KiCad file formats (s-expression, schematic and board) and the DRC
   engine and rule language. Cite file paths and line ranges.
6. Golden boards: apply the selection and vetting procedure in docs/08 §
   8.1 to pick the first 10 boards balanced across the domain list in
   docs/11 Q4; download them with licence records into `golden/`; run
   `kicad-cli` ERC/DRC on each as the oracle baseline and store the
   reports.
7. Write `config/owner-answers.yaml` from docs/11 with the defaults filled;
   the owner will edit it.

Work autonomously. Apply docs/11 defaults instead of asking. Stop only
for credentials or destructive actions. When Phase 0's exit criteria are
met and tests pass, write `docs/status.md` with what was built, what
deviated from the plan and why, and the exact first command for K1.
