# CLAUDE.md — standing instructions for building the PCB Second Brain

You are building the system specified in `docs/`. At the start of every
session read `SESSION-LOG.md` first (standing decisions D1–D12, open items,
work queue), then `README.md`, then the docs the current work item points
at. At the end of every session: commit, push, append a session entry to
`SESSION-LOG.md`, and rewrite `prompts/next-session-prompt.md`.
The owner is an electronic and electrical engineer who uses Altium and KiCad
professionally and wants to touch as little as possible. Assume competence in
electronics; do not assume time.

## Non-negotiable rules

1. **Deterministic tools are the guard rails, not the model.** Never implement a
   check as "ask the LLM whether this looks right" when a parser, a rule engine,
   `kicad-cli`, SPICE or an API lookup can answer it. See ADR-0003.
2. **No part, pin, footprint or value may be emitted that does not resolve to a
   record in the parts store.** Hallucinated parts are a release blocker. See ADR-0002.
3. **Every finding, rule and generated element carries provenance**: source
   document, page or section, retrieval timestamp. A finding without provenance
   is reported as "unverified opinion", never as a fact. See ADR-0004.
4. **We build our own engine.** KiCad and Altium are knowledge sources (read
   their installed libraries, docs, templates, examples; read KiCad's source
   from GitLab), interchange formats, and test oracles. This is a personal,
   never-distributed project (ADR-0001): porting KiCad algorithms is
   allowed only on a "port" verdict from the module's best-known-method
   review (docs/13 § 13.6), and the module header must carry
   `provenance: ported_from_kicad`. Modules judged on engineering quality
   (routing, placement, constraints, high-speed, impedance) are designed
   fresh from the best sources, not ported.
   **Never decompile, disassemble or instrument Altium binaries**; use
   Altium's scripting interface for black-box observation instead
   (docs/13 § 13.2b and § 13.5).
5. **Respect licences.** Do not ingest IPC, IEC, CISPR or other paid standards
   text. Encode rules derived from them with clause references only. Check each
   source's terms in `docs/03` before a crawler is written.
6. **Every phase ends with the evaluation harness passing its exit criteria**
   (`docs/09`). Do not start the next phase on a "mostly working" one.
7. **Golden boards are sourced and vetted from the web** (docs/08 § 8.1);
   the owner has none to give yet. All domains are equal; balance every
   quota across the domain list in docs/11 Q4.
8. **Platform is Windows 11.** Native Python via `uv`; Docker Desktop for
   Postgres; `kicad-cli.exe` from the local KiCad install as oracle only.

## Working conventions

- Language: Python 3.12+, managed with `uv`. Type hints everywhere, `ruff` and
  `mypy --strict` clean before any commit.
- Repository layout is defined in `docs/01-architecture.md` § Repository layout.
  Do not invent a different one.
- Tests: `pytest`. Every checker gets (a) a unit test on a synthetic netlist and
  (b) a fault-injection test on a golden board. No checker ships without both.
- Secrets live in `.env` (see `config/env.example`), never in code or docs.
- Commit small and often, one concern per commit, imperative mood subject line.
- When a tool, API or library behaves differently from what `docs/` assumes,
  update the doc in the same commit and note the discrepancy in
  `docs/12-sources.md` under "Corrections".
- When you need an owner decision, append it to
  `docs/11-open-questions-for-owner.md` with your recommended default, apply
  the default, and continue. Do not block.

## Model usage

- Agent orchestration and design review: `claude-opus-5-5`, adaptive thinking,
  effort `high` (raise to `xhigh` for review passes on boards with more than
  ~300 components).
- Hard second-opinion review passes and generation of new circuit blocks:
  `claude-fable-5-1`.
- Bulk extraction during ingestion (pin tables, parametrics, rule candidates):
  `claude-haiku-5-5` via the Message Batches API, with `claude-sonnet-5-5` for
  pages Haiku flags as low-confidence.
- Always use structured outputs (`output_config.format`) for extraction; never
  parse free text for data.
- Enable server-side refusal fallbacks (`fallbacks: "default"`) on every
  Opus/Fable call.
