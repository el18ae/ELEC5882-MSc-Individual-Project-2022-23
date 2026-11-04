# 00 — Vision and principles

## What we are building

A system that can take your seat for PCB work, in this order of ambition:

1. **Review** an existing design (KiCad or Altium) and produce a prioritised,
   cited list of realistic changes needed for worldwide component availability,
   fab/assembly manufacturability, EMI/EMC and signal-integrity hygiene, and
   datasheet-mandated circuitry around every major component.
2. **Source**: pick parts and alternates that are in stock across regions,
   not end-of-life, and compatible with the chosen assembler's library.
3. **Generate circuit blocks** around a given component from its datasheet,
   app notes and the corpus of open reference designs that use it.
4. **Generate whole designs** from a requirements document, as a tidy
   hierarchical schematic with a verified BOM.
5. **Lay out** the board (placement and routing) to rules, using routing
   engines, and prove it with DRC, impedance and simulation checks.

Every step is executed by Claude with a tool set, grounded in a knowledge base
you own, and gated by checks that do not depend on the model's opinion.

## Why existing AI PCB attempts produce unusable output

| Failure mode | Root cause | Countermeasure in this plan |
|---|---|---|
| Hallucinated parts, pin names, footprints | Model emits free text; nothing validates it against reality | Structured IR only; every part/pin must resolve to the parts store (ADR-0002) |
| Design "looks right" but fails ERC/DRC or won't build | No closed loop; no compile-check-repair cycle | Generate → build in our kernel → check → repair loop with hard gates (docs/06) |
| Missing mandatory circuitry (bootstrap caps, pull-ups, crystal load caps, RC snubbers) | Datasheet "Typical application" and "Layout guidelines" never read | Rule extraction from datasheets into machine-checkable rules with citations (docs/04, docs/07) |
| Parts obsolete or unavailable | Static knowledge | Live distributor/lifecycle APIs at review time; availability is a check, not a memory (docs/05) |
| Layout ignores physics | LLMs do not reason about current loops, return paths, impedance | Placement/routing delegated to engines; physics checks run as tools (openEMS/SPICE/impedance calculators); LLM only sets intents and reviews results |
| Unfamiliar with how experienced designers actually do it | No precedent | Design corpus: thousands of parsed open reference designs, queryable by part and motif (docs/02) |
| Nobody knows how good it is | No measurement | Evaluation harness with fault injection, precision/recall per checker, hallucination rate (docs/08) |
| Human still has to re-check everything | Trust extended all at once | Trust ladder: autonomy is granted per capability after it passes evals on your boards (docs/07) |

## The five principles

1. **The model proposes; the tools dispose.** Any claim that can be checked
   deterministically is checked deterministically.
2. **Structured, never free-form.** The model writes IR (JSON, atopile, SKiDL),
   compilers write KiCad files. The model never hand-writes an s-expression.
3. **Provenance everywhere.** A finding cites a page. A rule cites a clause.
   A generated block cites the datasheet figure and the reference designs it
   was patterned on.
4. **Vetted real boards are the benchmark.** The system is "as good as you"
   when it finds every change the original designers made between
   revisions of proven open boards, plus valid things they missed, with no
   false claims. That is a measurable statement (docs/08).
5. **Autonomy is earned per capability.** Review is auto-applied before
   generation is trusted; passive-component generation is trusted before
   power-stage generation; and so on.

## Explicit non-goals (for now)

- An interactive GUI editor. The system is an engine plus an agent; you
  open exports in Altium or KiCad when you want to look around (ADR-0001).
- RF/microwave layout above a few GHz, and HDI/microvia stackups. The checks
  will flag these as "out of calibrated scope" rather than guess.
- Training or fine-tuning a model. All intelligence comes from the knowledge
  base, the tools and prompting. (Revisit after Phase 6 if evals plateau.)
