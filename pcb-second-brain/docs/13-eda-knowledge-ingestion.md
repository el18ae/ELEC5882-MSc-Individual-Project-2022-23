# 13 — EDA knowledge ingestion (what to extract from KiCad and Altium)

Purpose: give the system the engineering knowledge embedded in two mature
EDA tools so that `pcbbrain-kernel` can be built correctly and so that the
agent can reason the way an experienced user of either tool does. This is
a parallel ingestion track to datasheets and the design corpus (docs/04).

## 13.1 KiCad (installed at e.g. `C:\Program Files\KiCad\10.0\`, plus source from GitLab)

| Asset | Where | Licence | What we extract | Used for |
|---|---|---|---|---|
| Symbol libraries (`.kicad_sym`) | `share\kicad\symbols\` | CC-BY-SA 4.0 with the KiCad library exception (designs made with them carry no obligation) | Every symbol: pin numbers, names, electrical types, units, alternate pin functions; the *conventions* (naming, pin type vocabulary, unit splitting) | Parts store seed (pin tables with `verified_by: kicad_lib`); symbol conventions for the kernel's own library format |
| Footprint libraries (`.pretty/*.kicad_mod`) | `share\kicad\footprints\` | same | Pad geometry, courtyard, fab/silk layers, 3D model refs; IPC-7351 naming; density levels | Footprint store; footprint verification (pad count vs pins); kernel footprint model |
| 3D models | `share\kicad\3dmodels\` | same | STEP/WRL refs | Height/keep-out checks (Phase 7) |
| Demo projects & templates | `share\kicad\demos\`, `share\kicad\template\` | GPL / CC | Clean example projects, hierarchical sheet idioms, DRU examples, stackup examples | Corpus seed; style reference |
| Documentation | `share\doc\kicad\` and docs.kicad.org | GPL-3 / CC-BY | Full user manuals, DRC rule language reference, file-format specs, PCB calculator formulas | Document store (`kind: eda_doc`); rule-language design input |
| Python bindings + IPC API protobufs | `bin\`, `lib\python3\`, `api\` | GPL-3 | API surface and data model | Kernel data-model design input; interchange tooling |
| **Source code** (clone `gitlab.com/kicad/code/kicad`) | not installed; clone | GPL-3 | Architecture notes extracted by reading: board/schematic data model, s-expression I/O, Altium/Eagle/CADSTAR/Allegro importers (how the binary Altium formats are parsed), DRC engine and rule evaluator, ERC pin-type matrix, PNS push-and-shove router, zone filler, Gerber/IPC-2581/ODB++ exporters, PCB calculator (impedance, track width, via) formulas, footprint wizards (IPC-7351 generators) | Written up as design notes with file:line citations; the kernel is implemented from the notes, clean-room (ADR-0001) |
| Library conventions (KLC) | kicad.org/libraries/klc | CC | The library convention rules themselves | Encoded as style rules in `rules/` |

## 13.2 Altium Designer (installed; version to be recorded)

Touch data files and documentation only, plus the scripting interface
(§ 13.2b). Never decompile, disassemble or instrument executables (§ 13.5).

| Asset | Where (typical) | What we extract | Used for |
|---|---|---|---|
| Installed libraries (`.IntLib`, `.SchLib`, `.PcbLib`) | `C:\Users\Public\Documents\Altium\AD<ver>\Library\` | Symbols, footprints, pin maps, parameters, supplier links. `.IntLib` is a compiled container; `.SchLib/.PcbLib` are OLE compound documents (parseable with `olefile`; KiCad's importer source documents the record formats) | Parts store seed; Altium file readers in the kernel |
| Templates | `...\Templates\` | Sheet templates, project templates, output job files (`.OutJob`), DRC rule sets inside `.PcbDoc` templates | Export presets; rule-set comparison against our rule store |
| Examples | `...\Examples\` | Complete reference projects with rules, stackups, variants, output jobs, draftsman docs | Corpus seed (Altium idioms); rule examples |
| Design rule definitions | inside `.PcbDoc` (rules section) and help | The full taxonomy of Altium rule types (clearance, width, via, diff pair, length matching, matched lengths, height, component clearance, daisy chain, SMT-to-corner, solder mask sliver, silk-to-mask, etc.) with their parameters and query language | The kernel's rule taxonomy is designed to be a superset of KiCad's DRU and Altium's rule types so both can be imported losslessly |
| Help / documentation | `C:\Program Files\Altium\AD<ver>\Help\` and techdocs | Reference for every rule, every editor command, PCB/schematic file concepts, scripting API | Document store (`kind: eda_doc`) |
| Scripting examples (`.PrjScr`, DelphiScript/JS/VB) | `...\Examples\Scripts\` | The object model (`PCB_Board`, `Sch_Document`, iterators) | Interchange and data-model design input |
| Project/panel defaults, preferences | `%APPDATA%\Altium\` | Default DRC, grids, layer stacks | Style defaults |

## 13.2b Black-box observation of Altium (legal, high value)

Altium Designer has a scripting system (DelphiScript, JavaScript, VBScript,
and Python via community bridges) that can open projects, run DRC, export
reports, iterate every primitive and rule, and drive polygon pours and
output jobs. Using it is ordinary licensed use. We use it to build an
**Altium oracle** and to learn Altium's *behaviour* without touching its
binaries:

| Experiment | What we learn | Feeds |
|---|---|---|
| Run DRC on every golden board and on synthetic rule-probe boards (one rule type per board, parameters swept) | Exact semantics of each rule type: how clearance treats arcs and thermal reliefs, how mask slivers and silk-to-mask are measured, what "matched length" tolerance means, how diff-pair gap is enforced on corners | Kernel `drc` semantics; `rule-taxonomy.md` mapping; `oracle.altium.drc` |
| Pour polygons on probe boards; export the resulting copper | Pour behaviour: thermal relief geometry, island removal, neck width handling | Kernel zone filler |
| Export Gerber X2/ODB++/IPC-2581 for golden boards | Reference outputs the kernel's exporters must match | Export tests |
| Iterate primitives on every installed example and library through the scripting object model | The complete object model and attribute set, as data | `io.altium` reader validation (what the reader parses must equal what the scripting API reports) |
| Run ERC on probe schematics | Altium's ERC matrix and defaults | Kernel `erc` |
| Export netlists and BOMs for golden boards | Reference netlists | `netlist` tests |

The scripts live in `pcbbrain-eval/altium-oracle/` and run on the owner's
machine. Their outputs are data about *our* designs produced by software
we license; nothing of Altium's code is copied.

## 13.5 Decompiling Altium: gain versus cost (decision: no)

| What decompiling could reveal | Is it available another way? | Quality of what decompiling would give |
|---|---|---|
| File format record layouts | Yes: KiCad's Altium importers (open source, documented), installed libraries/examples as samples, scripting object model | Redundant |
| Rule types, parameters, query language | Yes: help, templates, scripting docs | Redundant |
| DRC semantics and tolerances | Yes: black-box observation (§ 13.2b) gives the observable behaviour exactly | Compiled code without symbols is slower and less reliable than measuring the behaviour |
| Polygon pour, teardrop, mask algorithms | Mostly: observation plus KiCad's open implementations of the same algorithms | Low |
| Interactive/auto router internals | Partly: no open equivalent of Altium's router, but decompiled routing code is close to unreadable; Freerouting, KiCad PNS and the published literature teach the algorithms far better | Very low |
| Footprint wizard math (IPC-7351) | Yes: IPC-7351 itself, KiCad generators | Redundant |
| Library/vault internals, licensing, UI | Not needed | None |

Costs: breach of the EULA (a contract; private breach is still a breach,
enforcement risk aside), weeks of effort on unreadable output, and the
permanent loss of the option to show the project to Altium. Net: the gain
is near zero and the cost is the one thing the owner said they might want
to keep. Not done, and recorded as a standing rule in CLAUDE.md.

Note: this is an engineering assessment, not legal advice. Some
jurisdictions (EU, UK) have statutory interoperability exceptions for
decompilation that a licence cannot override, but they are narrow, require
that the information is not otherwise available, and the table above shows
it is otherwise available. If the owner ever wants to rely on such an
exception or to approach Altium, a lawyer should look first.

## 13.6 Port, adapt or design fresh (per-module policy)

Porting KiCad code is permitted (ADR-0001) but is never the default.
Each kernel module gets a **best-known-method review** in K1 before any
implementation: the design note compares at least three sources (KiCad,
Altium observed behaviour, Freerouting, other open tools such as Horizon
EDA, LibrePCB and pcb-rnd, and academic or vendor literature) and records
one verdict with reasons. The oracle cross-check applies to every module
regardless of verdict.

| Verdict | Meaning | Expected modules |
|---|---|---|
| **port** | KiCad's implementation is the convention or the best known; translate it, keep its tests, tag `provenance: ported_from_kicad` | `io.kicad`, ERC pin matrix, DRC check semantics (first version), zone filler (Clipper2), Gerber/IPC-2581/drill exporters, footprint generators, calculator formulas, Altium/Eagle importers, integer-nanometre coordinates |
| **adapt** | KiCad's approach is sound but the data model or interface must change for our needs | `model` (adds provenance, findings, intents, variants), `render`, netlist extraction |
| **design fresh** | KiCad is weak or absent; better sources exist | Rule engine and constraint model (Altium-style scoped rules with queries), `route` (Freerouting, negotiated congestion, rip-up-and-reroute, published RL approaches), `place` (analytic and constraint-driven placement), length/skew matching and diff-pair tuning, return-path and plane analysis, impedance (2D field solver rather than closed form), SI/PI and thermal estimates, everything AI-specific |

Rule of thumb: port what fabs and files already agree on; design fresh
where the result is judged by engineering quality rather than by
conformance.

## 13.3 Other engines ingested the same way

Freerouting (GPL-3): routing algorithms, DSN/SES formats. ngspice (BSD):
model and netlist semantics. atopile (MIT), SKiDL (MIT): IR design input.
`kiutils`, `kicad-sch-api`: s-expression handling patterns. Docling (MIT).

## 13.4 Method

1. Inventory: a script walks both install trees and records every file
   with type, size, sha256 and a licence tag from the table above.
2. Libraries → parsed into the parts/footprint/symbol stores with
   `verified_by: kicad_lib | altium_lib`.
3. Documentation → document store, chunked, embedded, cited like datasheets.
4. Source code → `design-notes/` written by Claude Code reading the code:
   one note per subsystem, with purpose, data structures, algorithms,
   edge cases, and file:line references. Notes are the input to kernel
   implementation; the kernel's tests include the edge cases found here.
5. Rule taxonomies → `schemas/rule-taxonomy.md`, the union of KiCad DRU
   and Altium rule types, with a mapping table in each direction.
6. Everything above is re-run when either tool is updated (version
   recorded per record).
