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

Touch data files and documentation only. Never decompile, disassemble or
instrument executables (EULA).

| Asset | Where (typical) | What we extract | Used for |
|---|---|---|---|
| Installed libraries (`.IntLib`, `.SchLib`, `.PcbLib`) | `C:\Users\Public\Documents\Altium\AD<ver>\Library\` | Symbols, footprints, pin maps, parameters, supplier links. `.IntLib` is a compiled container; `.SchLib/.PcbLib` are OLE compound documents (parseable with `olefile`; KiCad's importer source documents the record formats) | Parts store seed; Altium file readers in the kernel |
| Templates | `...\Templates\` | Sheet templates, project templates, output job files (`.OutJob`), DRC rule sets inside `.PcbDoc` templates | Export presets; rule-set comparison against our rule store |
| Examples | `...\Examples\` | Complete reference projects with rules, stackups, variants, output jobs, draftsman docs | Corpus seed (Altium idioms); rule examples |
| Design rule definitions | inside `.PcbDoc` (rules section) and help | The full taxonomy of Altium rule types (clearance, width, via, diff pair, length matching, matched lengths, height, component clearance, daisy chain, SMT-to-corner, solder mask sliver, silk-to-mask, etc.) with their parameters and query language | The kernel's rule taxonomy is designed to be a superset of KiCad's DRU and Altium's rule types so both can be imported losslessly |
| Help / documentation | `C:\Program Files\Altium\AD<ver>\Help\` and techdocs | Reference for every rule, every editor command, PCB/schematic file concepts, scripting API | Document store (`kind: eda_doc`) |
| Scripting examples (`.PrjScr`, DelphiScript/JS/VB) | `...\Examples\Scripts\` | The object model (`PCB_Board`, `Sch_Document`, iterators) | Interchange and data-model design input |
| Project/panel defaults, preferences | `%APPDATA%\Altium\` | Default DRC, grids, layer stacks | Style defaults |

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
