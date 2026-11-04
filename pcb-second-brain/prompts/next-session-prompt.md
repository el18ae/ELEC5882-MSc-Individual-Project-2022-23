# Next-session prompt (cloud, session 2)

Paste this as the first message of the next Claude Code cloud session on
the repository `el18ae/ELEC5882-MSc-Individual-Project-2022-23`.

---

Continue the PCB Second Brain project. First run:

```
git fetch origin pcb-second-brain-plan && git checkout pcb-second-brain-plan
```

Then read, in order: `pcb-second-brain/SESSION-LOG.md` (standing
decisions, open items, work queue), `pcb-second-brain/README.md`,
`pcb-second-brain/CLAUDE.md`, and the docs the work item points at. Do not
re-open settled decisions D1–D11.

This session's work item is **W1 — golden board selection** (see the work
queue in SESSION-LOG.md and the selection and vetting criteria in
`docs/08-evaluation-harness.md` § 8.1):

1. Build a candidate list of 40–60 open-hardware boards balanced across the
   nine domains in `config/owner-answers.yaml`, drawing from the pools in
   docs/08 and docs/03 § 3.3. For each: name, source URL, licence (quoted
   from the repo), CAD format, whether schematic + board + libraries are
   complete, number of revisions or tagged releases, presence of an issue
   tracker or errata, and the domain it covers. Verify each entry against
   its live page; do not list a board you could not confirm.
2. Rank within each domain and mark the first 10 overall for Phase 0.
3. Where the network policy allows, clone or download the first 10 into
   `pcb-second-brain/golden/<board>/` with a `PROVENANCE.yaml` (URL,
   commit or release, licence, fetched_at). If downloads are blocked,
   record the exact fetch commands for the local phase instead.
4. For boards with git history, list the revision pairs and commit
   messages that look like design fixes (decouplers, pull-ups,
   terminations, footprints, part swaps, clearances, silkscreen) as
   candidate ground-truth findings in `golden/<board>/GROUND-TRUTH-CANDIDATES.md`.
   Do not judge them yet; that needs the kernel.
5. Write `golden/CANDIDATES.md` with the full table and the Phase 0 pick.

Budget: the owner is on limited credits. Be economical: shallow clones
only (`--depth 1`, or `--filter=blob:none` when history is needed for
step 4), do not re-read files the session log already summarises, batch
web lookups, and stop at the end of W1 rather than starting W2.

Work autonomously; apply docs/11 defaults rather than asking. At the end:
commit on `pcb-second-brain-plan`, push, append a Session 2 entry to
`SESSION-LOG.md` (decisions, artefacts, open items, next work item), and
rewrite this file for session 3 (W2 — rule taxonomy) in the same format.
Finish with a short recap of what was found, what was done and what is
next.
