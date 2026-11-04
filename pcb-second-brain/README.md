# PCB Second Brain — Project Plan

An AI-operated PCB design and review system built on a curated, provenance-tracked
knowledge base, with deterministic engineering checkers as the guard rails and
Claude (Opus 5.5 / Fable 5.1) as the operator.

This folder is the complete plan. It is written so that a local Claude Code
session can execute it phase by phase with minimal human input. Nothing here is
code yet; it is the specification, the decisions, the schemas and the roadmap.

## Read in this order

| # | File | What it answers |
|---|------|-----------------|
| 0 | `docs/00-vision-and-principles.md` | What we are building, why current AI PCB tools fail, the five principles that fix that |
| 1 | `docs/01-architecture.md` | The three layers (Knowledge, Tools, Agent) and how data flows between them |
| 2 | `docs/02-knowledge-base.md` | The six stores, what goes in each, how they are queried |
| 3 | `docs/03-data-sources-and-licensing.md` | Where the data comes from and what you are allowed to do with it |
| 4 | `docs/04-ingestion-pipeline.md` | Acquire → parse → normalise → extract → index, with provenance |
| 5 | `docs/05-tools-and-mcp-server.md` | The tool surface the agent gets (EDA, sourcing, analysis, generation) |
| 6 | `docs/06-agent-workflows.md` | The five workflows, from design review to full autonomous design |
| 7 | `docs/07-verification-and-guardrails.md` | The five verification tiers and the trust ladder |
| 8 | `docs/08-evaluation-harness.md` | How we measure the system so "better than me" is a number, not a feeling |
| 9 | `docs/09-roadmap.md` | Phases, milestones, exit criteria, order of work |
| 10 | `docs/10-stack-and-hardware.md` | Concrete software stack, local hardware, model choice, cost model |
| 11 | `docs/11-open-questions-for-owner.md` | **The instructions you have not given yet.** Answer these before the local kickoff |
| 12 | `docs/12-sources.md` | Links consulted while writing this plan |
| 13 | `docs/13-eda-knowledge-ingestion.md` | What to extract from the installed KiCad and Altium, and from KiCad's source, to build our own engine |
| — | `decisions/ADR-*.md` | The four architecture decisions that everything else depends on |
| — | `schemas/` | First-cut data schemas (parts, rules, findings, design graph) |
| — | `prompts/local-kickoff-prompt.md` | The prompt to paste into local Claude Code to start Phase 0 |
| — | `CLAUDE.md` | Standing instructions for the local Claude Code session that builds this |

## One-paragraph summary

The system never lets a language model free-hand a circuit. The model proposes
structured intermediate representations (parts, nets, placement intents) that
our own engine (the kernel) builds into designs, then checks with its own
ERC/DRC, manufacturer rule sets, datasheet-derived rules, sourcing APIs,
simulation and precedent lookups against a corpus of vetted open reference
designs. KiCad and Altium are ingested as knowledge to build that engine, read
and written as interchange formats, and used as independent test oracles. Every finding
carries a citation back to a datasheet page, a rule, a manufacturer document or
a reference design. The model's job is orchestration, judgement and explanation;
the tools' job is truth. Trust is extended one rung at a time as the evaluation
harness proves each capability on vetted real boards.

## Where this lives

This folder sits inside your MSc project repository for now. It is
self-contained; move it to a dedicated repository (recommended name:
`pcb-brain`) whenever convenient. Nothing in it references the MATLAB files.
