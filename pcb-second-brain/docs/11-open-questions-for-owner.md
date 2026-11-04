# 11 — Open questions for the owner (the instructions you have not given yet)

Each question has a default. If you do not answer, the local build proceeds
with the default and records that it did. Answer in this file, in the chat
with local Claude Code, or in `config/owner-answers.yaml`.

## Must answer before Phase 0

| # | Question | Why it matters | Default if unanswered |
|---|---|---|---|
| Q1 | What stack did your two existing second brains use (database, vector store, embedding model, orchestration, OS)? Should this one reuse it or stand alone? | Reusing proven plumbing saves weeks; mixing stacks costs maintenance | Stand-alone Postgres + pgvector as in docs/10 |
| Q2 | Local machine spec: CPU, RAM, GPU, free NVMe, OS. Is Altium on the same machine? | Sizes the corpus and decides local vs cloud embeddings; Windows vs Linux affects headless KiCad | Linux, 64 GB, 2 TB, no GPU |
| Q3 | Which of your own past boards can be used as golden boards, in what format (Altium/KiCad), and which may be sent to the Anthropic API? List known issues per board if you remember them | They are the benchmark; confidentiality decides T5 usage | None; evals use open boards only until provided |
| Q4 | Your domains and their priority: power electronics (your MSc work), mixed-signal/instrumentation (your refinery job), MCU/IoT, high-speed digital, RF, motor drive, industrial comms (RS-485/CAN/4–20 mA/HART), intrinsically safe/ATEX? | Decides which parts, rules, standards and corpus to ingest first | Power + mixed-signal instrumentation + MCU first |
| Q5 | Target fabs/assemblers (e.g. JLCPCB, PCBWay, Eurocircuits, a local house) and typical processes (layer count, finish, impedance control) | Fab capability YAML and DRU are per fab | JLCPCB 2/4-layer + Eurocircuits |
| Q6 | Privacy stance: may datasheets/corpus be embedded via a cloud embedding API, or must embeddings be local? May owner boards be sent to the Anthropic API? | Chooses embedding path and whether T5 runs on your boards | Cloud embeddings for public docs; owner boards to Anthropic API: yes |
| Q7 | Regions for "worldwide availability" and regulatory targets (CE/UKCA/FCC, industrial EMC class, ATEX/IECEx, automotive AEC-Q) | Defines the availability rule and which EMC rule sets apply | US + EU + UK + APAC; CE/UKCA/FCC Class B |
| Q8 | API accounts you already have or will create: Digi-Key, Mouser, Nexar tier, element14, TME, LCSC; any SiliconExpert/Z2Data access | Sourcing coverage; Nexar tier determines lifecycle field availability | Digi-Key + Mouser + Nexar evaluation tier |
| Q9 | Standards you own copies of (IPC-2221/2222/7351/A-610, IEC/CISPR) and whether you want their rules encoded (clause-referenced, no text ingestion) | Legal path for standards-derived rules | Rules from public app notes only |
| Q10 | Monthly API budget ceiling for ingestion and for daily use | Sets batch sizes and effort levels | $150/month ingestion during Phases 3–4, $100/month use |

## Should answer before Phase 2–3

| # | Question | Default |
|---|---|---|
| Q11 | Preferred passive series and vendors (e.g. Yageo RC / Samsung CL / Murata GRM), preferred connector families, preferred MCU families and regulators | Learn from owner boards |
| Q12 | Your schematic conventions: sheet structure, refdes scheme, net naming, title block, DNP marking, test point policy | Learn from owner boards; codify as style rules |
| Q13 | Do you want Altium output at the end of every workflow, or KiCad only with Altium export on request? | KiCad only; export on request |
| Q14 | Books you own that may be used as background for the critic (personal use) | None |
| Q15 | Languages/units: metric only? Imperial for legacy boards? | Metric, mil accepted on input |

## Should answer before Phase 6–7

| # | Question | Default |
|---|---|---|
| Q16 | Willingness to license a commercial router (Quilter/DeepPCB) if Freerouting underperforms on the layout eval | Decide at Phase 7 on eval evidence |
| Q17 | Mechanical inputs: do your projects come with STEP/DXF envelopes and connector positions? | Assume yes for Workflow E |
| Q18 | How "tidy" is defined for you: grid, orientation rules, silkscreen style, test-point placement, panel conventions | Learn from owner boards; codify |

## Instructions the owner may want to give explicitly (recommended)

- A written definition of done for "takes my seat": which workflows, at
  which trust rung, by when.
- Permission scope for autonomous runs: what the system may change in a
  repo without asking (suggest: anything on a working branch; never main).
- A list of clients/projects whose data must never be used as corpus or
  sent to any API.
- Whether the system may order parts or submit fab orders (suggest: no,
  ever, without explicit per-order approval).
