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

## Consequences
- A large kernel workstream is added (docs/09 Phase K). It runs in
  parallel with the knowledge work and is on the critical path for the
  first deterministic review (T0–T1).
- Licensing: learning from GPL-3 source is fine; copying KiCad code into
  the kernel would make the kernel GPL-3 on distribution. If the owner
  ever intends to sell or distribute the software, the clean-room rule is
  mandatory; if it stays private, GPL obligations do not trigger, but the
  rule is kept anyway so the option stays open (docs/11 Q19).
- Altium's EULA prohibits decompiling or reverse-engineering its binaries.
  Reading installed data files, documentation, templates, libraries and
  scripting examples for personal use is not that. The extraction tooling
  in docs/13 touches no executables.
- Early capability lags a wrapper approach by a few months; in exchange
  every check, format and rule is the owner's and can be made to fit the
  verification tiers exactly.
