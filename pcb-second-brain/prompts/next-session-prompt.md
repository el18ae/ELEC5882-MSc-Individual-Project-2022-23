# Next-session prompt (cloud, session 3)

Paste this as the first message of the next Claude Code cloud session on
the repository `el18ae/ELEC5882-MSc-Individual-Project-2022-23`.
Recommended model: Opus 5.5 (see SESSION-LOG.md, model and budget guidance).

---

Continue the PCB Second Brain project. First run:

```
git fetch origin pcb-second-brain-plan && git checkout pcb-second-brain-plan
```

Then read, in order: `pcb-second-brain/SESSION-LOG.md` (standing
decisions, open items, work queue), `pcb-second-brain/README.md`,
`pcb-second-brain/CLAUDE.md`, and the docs the work item points at. Do not
re-open settled decisions D1–D12.

This session's work item is **W2 — rule taxonomy** (see the work queue in
SESSION-LOG.md; K1 exit criterion in `docs/09-roadmap.md`: "taxonomy covers
100% of rule types in golden boards"; rule record format in
`schemas/rule.schema.json`; Altium and KiCad knowledge sources in
`docs/13-eda-knowledge-ingestion.md`):

1. Enumerate the KiCad custom-rule language (`.kicad_dru`): every
   constraint type, every condition function and property usable in
   `(condition ...)`, layer and severity clauses, and the built-in board
   setup constraints that are not DRU rules (`.kicad_pro`
   `design_settings`, net classes). Source of truth is KiCad's source:
   the GitHub mirror `KiCad/kicad-source-mirror` (if gitlab.com is
   refused), files under `pcbnew/drc/` (rule parser, `drc_rule.h`
   constraint enum, `drc_test_provider_*`) and the DRU syntax help text.
   Use a blobless, sparse clone of those paths only.
2. Enumerate Altium Designer design-rule types (electrical, routing, SMT,
   mask, plane, testpoint, manufacturing, high speed, placement, signal
   integrity) with their parameters and query-language scoping. Use
   public Altium documentation if reachable. Otherwise use KiCad's Altium
   importer source (`pcbnew/pcb_io/altium/` rule parsing) as the record of
   rule kinds and fields. Record which source each entry came from.
3. Write `pcb-second-brain/schemas/rule-taxonomy.md` with:
   - one canonical taxonomy of rule kinds for the kernel (id, meaning,
     parameters with units, scope model);
   - a bidirectional mapping table: canonical ↔ KiCad DRU / board setup ↔
     Altium rule type, with "no equivalent" and "lossy" marked and the
     loss explained;
   - how scoping (KiCad conditions vs Altium queries) maps onto one
     canonical scope expression;
   - how Eagle DRU settings (`.dru` / `<designrules>` in `.brd`) map, as
     a short third column, because two Phase 0 boards are Eagle.
   Check that the taxonomy fits `schemas/rule.schema.json`. Update the
   schema in the same commit if it does not.
4. Coverage check against the golden boards. For each Phase 0 board in
   `golden/phase0.list`, fetch only the rule-bearing files (`*.kicad_dru`,
   `*.kicad_pro`; for Eagle, the `<designrules>` block of the `.brd`)
   with `git show <commit>:<path>` from a blobless, no-checkout clone at
   the commit pinned in its `PROVENANCE.yaml`. Do not check out the full
   CAD. List every rule or constraint used and confirm each maps to a
   canonical kind. The Altium boards (odrive-v3) need the local phase,
   so record them as pending.
5. Note in the taxonomy where fab capability YAML (W7) and interface rules
   (W8) will plug in, without writing them.

Budget: the owner is on limited credits. Be economical: sparse, blobless
clones of the needed paths only; fetch single files with `git show`, not
full checkouts. Do not re-read files the session log already summarises.
Batch web lookups. Stop at the end of W2 rather than starting W3.

Work autonomously; apply docs/11 defaults rather than asking. At the end:
commit on `pcb-second-brain-plan`, push, append a Session 3 entry to
`SESSION-LOG.md` (decisions, artefacts, open items, next work item), and
rewrite this file for session 4 (W3 — K1 design notes part 1, Fable 5.1
recommended) in the same format. Finish with a short recap of what was
found, what was done and what is next.
