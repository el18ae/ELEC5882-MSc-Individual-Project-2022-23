# ADR-0001 — Own EDA engine; KiCad and Altium are knowledge sources, interchange formats and test oracles

Status: accepted (owner decision, replaces the earlier "KiCad as canonical substrate")

## Context
The owner wants the software built their way, not as a wrapper around
KiCad. The agent does not need an interactive editor; it needs a data
model, parsers and writers, a rule/DRC engine, a renderer, a router and
exporters that it can drive programmatically and that can be shaped to the
verification architecture in docs/07. KiCad (GPL-3, ~30 years of C++) and
Altium (proprietary, installed locally) contain most of the engineering
knowledge needed to build those pieces well.

## Decision
1. The system has its own kernel (`pcbbrain-kernel`): geometry, design
   data model, file I/O, DRC/rule engine, netlist, renderer, exporters,
   and later a router. It is the only thing the agent edits through.
2. KiCad and Altium are **ingested as knowledge** (docs/13): KiCad source,
   documentation and libraries; Altium's installed libraries, templates,
   examples, rule definitions, help and scripting documentation. The
   kernel is written clean-room from that knowledge, not by copying code.
3. KiCad and Altium are **interchange formats**: the kernel reads and
   writes KiCad project files and reads Altium files, so existing designs
   and the open corpus flow in, and the owner can open results in either
   tool.
4. `kicad-cli` (and Altium, manually) are **test oracles only**: in the
   evaluation harness, the kernel's DRC/ERC/netlist output is cross-checked
   against theirs on every golden board. They are never a runtime
   dependency of a workflow.
5. No interactive GUI editor is built. A static/web renderer for review is
   enough; the owner opens exports in Altium or KiCad when they want to
   look around.

## Project status (owner statement, November 2026)
This is a personal project. It will never be sold, distributed or
open-sourced. It stays on the owner's machine and produces designs for
the owner's own hobby use only. The owner may one day show it to Altium,
subject to legality and permissions, so the code must stay in a state
that could be shown to Altium without embarrassment or infringement.

## Decision on source licences and decompiling (given the status above)
- **KiCad (GPL-3):** because nothing is distributed, GPL obligations do
  not trigger, so porting KiCad algorithms directly into the kernel is
  permitted and faster than strict clean-room. Every module that ports
  KiCad code is tagged `provenance: ported_from_kicad` in its header, so
  those modules can be identified and rewritten clean-room if the owner
  ever decides to show or hand the code to Altium (at that moment it
  becomes distribution and GPL terms apply; Altium would not accept
  GPL-derived code).
- **Altium binaries: not decompiled.** See docs/13 § 13.5 for the
  gain-versus-cost analysis. Summary: the EULA prohibits it regardless of
  distribution; the engineering knowledge it would yield is almost
  entirely available from documentation, installed data files, KiCad's
  open importers and black-box observation through Altium's own
  scripting; compiled code without symbols is a poor teacher for
  algorithms; and decompiling is the single act that would close the
  "take it to Altium" door permanently.
- **Altium black-box observation is used instead** (docs/13 § 13.2b):
  running Altium on test designs through its scripting interface and
  recording what its DRC, polygon pour, rule evaluator and exporters
  produce. That is ordinary licensed use and makes Altium a second
  automated oracle.

## Consequences
- A large kernel workstream is added (docs/09 Phase K). It runs in
  parallel with the knowledge work and is on the critical path for the
  first deterministic review (T0–T1).
- Licensing is handled by the section above: relaxed for KiCad with
  per-module provenance tags; no decompiling of Altium; black-box
  observation of Altium through scripting.
- Early capability lags a wrapper approach by a few months; in exchange
  every check, format and rule is the owner's and can be made to fit the
  verification tiers exactly.
