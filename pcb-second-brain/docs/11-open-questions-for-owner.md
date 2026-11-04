# 11 — Open questions for the owner

Answered items are recorded in `config/owner-answers.yaml`. Unanswered
items proceed on the default and are logged.

## Answered (from the owner, November 2026)
| # | Question | Answer |
|---|---|---|
| Q-A | Wrap KiCad or build our own? | **Build our own engine.** KiCad and Altium are ingested as knowledge, used as interchange formats and as test oracles (ADR-0001). |
| Q4 | Domain priority | **None; all domains equal.** Ingestion and rule quotas are balanced across: power electronics, mixed-signal/instrumentation, MCU/IoT, high-speed digital, RF, motor drive, industrial comms (RS-485/CAN/4–20 mA/HART), intrinsically safe/ATEX, automotive. |
| Q3 | Golden boards | **None available.** The system selects and vets its own from the web (docs/08 § 8.1). |
| Q2 | Machine | **Windows 11** with KiCad and Altium installed; their directories may be used to extract all relevant knowledge (docs/13). Specs to be measured in Phase 0. |
| Q19 | Distribution | **Never sold, distributed or open-sourced; personal project; may be shown to Altium one day subject to legality.** GPL clean-room relaxed with per-module provenance tags (ADR-0001). |
| Q22 | Decompile Altium binaries? | **No.** Gain is near zero versus documentation, installed data, KiCad's open importers and black-box observation through Altium scripting; cost includes losing the Altium option (docs/13 § 13.5). |

## Must answer before Phase 0 (defaults apply otherwise)
| # | Question | Why | Default |
|---|---|---|---|
| Q1 | Stack of your two existing second brains (database, vector store, embedding model, orchestration) | The plumbing (ingestion, chunking, embeddings, retrieval, provenance, agent harness) is domain-agnostic and reusable even though the knowledge is not | Stand-alone Postgres + pgvector |
| Q5 | Target fabs/assemblers and processes | Fab rule sets are per fab | JLCPCB 2/4-layer + Eurocircuits + PCBWay |
| Q6 | Privacy: cloud embeddings for public documents OK? | Chooses embedding path | Yes |
| Q7 | Regions for "worldwide availability"; regulatory targets | Availability rule; EMC rule sets | US + EU + UK + APAC; CE/UKCA/FCC Class B; industrial EMC where relevant |
| Q8 | API accounts (Digi-Key, Mouser, Nexar tier, element14, TME, LCSC) | Sourcing coverage | Digi-Key + Mouser + Nexar evaluation tier |
| Q9 | Standards you own (IPC-2221/2222/7351/A-610, IEC/CISPR) for clause-referenced rule encoding | Legal rule path | Rules from public app notes only |
| Q10 | Monthly API budget ceiling | Batch sizes, effort | $300/month during ingestion phases, $150/month use |
| Q20 | May `kicad-cli` be used as the T1 engine *inside* Workflow A temporarily, until the kernel's DRC (K3) replaces it? Shortens time-to-first-review by months; end state unchanged | Schedule | No (pure own-engine path) |
| Q21 | Name of the software/repo | Naming | `pcb-brain` |
| Q23 | Intrinsically safe / ATEX golden boards: no open-licensed IS/ATEX design with complete CAD exists (W1, golden/CANDIDATES.md § 4). Supply or approve real ones? | ATEX quota of the golden set | Synthetic boards built from IEC 60079-11 clause-referenced rules, plus vendor IS reference designs (TI TIDA-01504, ADI CN0270) held locally under their own licences as reference only, not golden |
| Q24 | Licences accepted for golden boards beyond the docs/08 list (CC-BY/CC-BY-SA/CERN-OHL/MIT/TAPR/OCP) | Size of the candidate pool | Also accept Apache-2.0, 0BSD, BSD, GPL-2.0/3.0, Solderpad; accept non-commercial (CC-BY-NC-SA) with a flag, since the project is never distributed (D5) |
| Q25 | This repository is public, so third-party CAD is kept out of git (D12) and re-fetched by `golden/fetch-phase0.sh`. Move the project to a private `pcb-brain` repository? | Where golden CAD can live | Keep CAD out of git; move to a private repo when Q21 is settled |

## Should answer before Phase 2–3
| # | Question | Default |
|---|---|---|
| Q11 | Preferred passive series, connector families, MCU and regulator families | Learn from corpus statistics |
| Q12 | Schematic conventions (sheet structure, refdes, net naming, DNP, test points) | Adopt KiCad KLC + Altium template conventions; codify as style rules |
| Q13 | Deliverable format: KiCad only, Altium via importer, or an Altium writer in the kernel | KiCad files; Altium opens them via its importer |
| Q14 | Books you own for critic background | None |
| Q15 | Units | Metric, mil accepted |

## Should answer before Phase 6–7
| # | Question | Default |
|---|---|---|
| Q16 | Licence a commercial router if own router + Freerouting underperform | Decide on eval evidence |
| Q17 | Mechanical inputs (STEP/DXF envelopes, connector positions) | Assume yes |
| Q18 | Definition of "tidy" (grid, orientation, silkscreen, test points, panels) | Derive from vetted golden boards; codify |

## Recommended explicit instructions
- Definition of done for "takes my seat": workflows, trust rungs, dates.
- Permission scope for autonomous runs (suggest: anything on a working branch; never main).
- Projects/data that must never be used as corpus or sent to any API.
- No part orders or fab submissions without per-order approval.
