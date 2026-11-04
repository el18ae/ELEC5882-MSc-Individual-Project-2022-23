# Next-session prompt (cloud, session 4)

Paste this as the first message of the next Claude Code cloud session on
the repository `el18ae/ELEC5882-MSc-Individual-Project-2022-23`.
Recommended model: Fable 5.1 (design notes with architectural verdicts; see
SESSION-LOG.md, model and budget guidance). Opus 5.5 is acceptable if
budget is tight.

---

Continue the PCB Second Brain project. First run:

```
git fetch origin pcb-second-brain-plan && git checkout pcb-second-brain-plan
```

Then read, in order: `pcb-second-brain/SESSION-LOG.md` (standing
decisions, open items, work queue), `pcb-second-brain/README.md`,
`pcb-second-brain/CLAUDE.md`, and the docs the work item points at. Do not
re-open settled decisions D1–D12.

This session's work item is **W3 — K1 design notes, part 1** (work queue
in SESSION-LOG.md; the port / adapt / design-fresh verdict policy is in
`docs/13-eda-knowledge-ingestion.md` § 13.6; K1 exit criteria are in
`docs/09-roadmap.md`; kernel module names and layout are in
`docs/01-architecture.md`):

1. Clone KiCad from GitLab as a blobless, sparse clone of only the paths
   each note needs (`gitlab.com` is reachable; docs.kicad.org is not).
   Pin to commit `4c23271608ae7fa9e584c6f63b57e54163832c48` so the notes
   line up with `schemas/rule-taxonomy.md` and
   `design-notes/w2-extracts/`. Reuse those extracts rather than re-reading
   what they already cite (DRU language, board-setup keys, DRC and ERC
   codes, pin matrix, test providers).
2. Write one design note per subsystem under `pcb-second-brain/design-notes/`:
   - `io-kicad.md`: S-expression `.kicad_pcb` / `.kicad_sch` /
     `.kicad_sym` / `.kicad_mod` / `.kicad_pro` readers and writers,
     version handling (KiCad 5 legacy `.sch`/`.lib` and KiCad 5 `(setup)`
     included, since six Phase 0 boards use them), round-trip fidelity
     requirements.
   - `erc.md`: pin-type matrix, net-driver resolution, ERC checks and
     severities.
   - `drc-engine.md`: rule compilation, condition evaluator (libeval),
     constraint resolution (implicit rules, netclass, item overrides,
     last-match-per-value precedence), the test-provider architecture,
     spatial indexing, performance.
   - `zone-filler.md`: fill algorithm, thermal relief, island removal,
     clearance inflation, polygon library dependencies, determinism.
   Each note: what KiCad does (file:line citations, paraphrased; no large
   GPL copies), its known weaknesses (from source comments, issue tracker
   where reachable), alternatives from the literature or other open tools,
   and a **verdict** (port / adapt / design fresh) with reasons, a
   proposed kernel interface (Python signatures, DesignIR types from
   `schemas/design-graph.md`, canonical rule kinds from
   `schemas/rule-taxonomy.md`), and test strategy (oracle = `kicad-cli`
   on golden boards in L1). Any "port" verdict must name the files and
   the `provenance: ported_from_kicad` header requirement (CLAUDE.md
   rule 4).
3. Record cross-cutting decisions (geometry/polygon library choice,
   integer nm coordinates, determinism) in the notes. If a decision is the
   owner's, add it to docs/11 with a default.

Budget: the owner is on limited credits, and session 3 overspent (about
3.9 M subagent tokens in one workflow). Work single-agent unless the owner
explicitly opts in to a multi-agent workflow this session. Use sparse,
blobless clones of needed paths only. Do not re-read files the session log
or the W2 extracts already summarise. Batch web lookups. Stop at the end
of W3 rather than starting W4.

Work autonomously; apply docs/11 defaults rather than asking. At the end:
commit on `pcb-second-brain-plan`, push, append a Session 4 entry to
`SESSION-LOG.md` (decisions, artefacts, open items, next work item), and
rewrite this file for session 5 (W4 — K1 design notes part 2: Altium
importer internals as the spec for `io.altium`, exporters, footprint
generators, calculator formulas) in the same format. Finish with a short
recap of what was found, what was done and what is next.
