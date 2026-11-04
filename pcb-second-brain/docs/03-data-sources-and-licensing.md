# 03 — Data sources and licensing

Principle: a source is ingested only after its licence/terms have been read
and recorded in the `document.licence` or `project.licence` column. "Found on
the internet" is not a licence. Where a source permits personal/internal use
but not redistribution, it may be ingested and must never be published,
shared or pushed to a public repo. The implementation stores a `licence_note`
per source so this is auditable.

## 3.1 Component and sourcing data (APIs)

| Source | What it gives | Access / terms | Use in plan |
|---|---|---|---|
| Digi-Key API | Parametrics, stock, pricing, datasheet links, lifecycle flags, regional stock | Free developer registration; rate limits; ToS forbids bulk re-publication | Primary parametric + availability source (Americas, EU, APAC endpoints) |
| Mouser API | Stock, pricing, lifecycle, datasheet links | Free API key; limited calls/day | Second availability source |
| Nexar (Octopart) API | Aggregated stock across many distributors, lifecycle, specs, CAD, alternates | GraphQL; free evaluation tier is small; paid tiers vary in which fields (lifecycle, lead time) are included. Verify current plan page before relying on it | Aggregation + alternates. Treat tier-dependent fields as optional |
| LCSC / JLCPCB parts | Stock, pricing, basic/extended assembly class | Public web + community CSV dumps; no formal API ToS for scraping, so throttle and cache | Assembler library check for JLCPCB builds |
| element14/Farnell, TME, Arrow, RS | Regional stock (EU/UK), pricing | APIs with registration | "Worldwide" means at least 2 regions × 2 distributors in stock |
| SiliconExpert / Z2Data / IHS | Authoritative lifecycle, PCN/EOL notices, risk | Paid, enterprise | Not in scope unless you already have access. Lifecycle otherwise = union of distributor flags + manufacturer PCN pages |
| Manufacturer PCN/EOL pages | Change notices | Public | Crawl for parts in your BOMs only |

## 3.2 Datasheets, app notes, design guides

| Source | Terms | Notes |
|---|---|---|
| Manufacturer datasheets (TI, ADI, NXP, ST, Microchip, Infineon, Renesas, onsemi, Nexperia, Diodes, Vishay, Murata, TDK, Würth, Samtec, Molex, JST, Hirose, Amphenol, Espressif, Nordic, Silicon Labs, Lattice, Microchip FPGA, AMD/Xilinx, Intel/Altera, NVIDIA Jetson docs) | Generally downloadable for use; most ToS prohibit redistribution | Ingest for parts you resolve; keep local; cite by document number + page |
| Application notes and design guides (TI SLxA/SNVA series, ADI AN-xxx and CN-0xxx, NXP AN, ST AN, Würth design guides, Murata, TDK, Samtec SI guides, Intel/AMD platform design guides where public) | As above | Richest source of *rules*. Prioritise "layout guidelines", "PCB design", "EMC" titles |
| IPC standards (IPC-2221/2222 generic & rigid design, IPC-7351 land patterns, IPC-A-610 acceptability, IPC-6012) | **Paid, copyrighted. Do not ingest text.** | Buy the ones you use; encode the rules you apply as YAML with clause references; the system cites the clause, never quotes the text |
| IEC/CISPR/EN EMC standards (CISPR 32, EN 55032, IEC 61000-4-x), FCC Part 15 (public), UKCA guidance (public) | FCC text is public domain; IEC/CISPR paid | Same approach: limits and test setups as structured data with clause refs. Open explanations (e.g. Würth/TI EMC app notes) carry the prose |
| JEDEC (packages, JESD standards) | Free with registration; redistribution restricted | Package outlines for footprint verification |
| Ucamco Gerber X2/X3 spec, IPC-2581 schema | Free, open | Export formats; ingest fully |
| Books (Ott, Johnson & Graham, Bogatin, Archambeault) | Copyrighted; you may own copies | Personal-use only; strong background for the *critic* pass; never quote at length. Decide in docs/11 whether to include |
| Forums (KiCad forum, EEVblog, EDAboard, Reddit r/PrintedCircuitBoard) | Mixed | Low priority, high noise; only the KiCad forum and vendor E2E/EngineerZone forums (TI, ADI) for errata-style knowledge, tagged `confidence: low` |

## 3.3 Open design corpus (CAD files)

| Source | Format | Licence | Why |
|---|---|---|---|
| Adafruit (github.com/adafruit, "Adafruit-*-PCB" repos) | Eagle, KiCad | CC-BY-SA / MIT (varies) | Hundreds of clean breakout designs; excellent motif source |
| SparkFun (github.com/sparkfun) | Eagle, KiCad | CC-BY-SA 4.0 | Same |
| Olimex | Eagle, KiCad | GPL / CC-BY-SA | SBCs, more complex boards |
| Raspberry Pi (Pico, Pico 2, Compute Module IO boards) | KiCad, Altium PDFs | CC-BY / design-file licence per repo | Reference for RP2040/RP2350, DDR/PCIe-adjacent routing on CM IO boards |
| BeagleBoard.org (BeagleBone, BeaglePlay, BeagleV) | Altium / KiCad / Allegro | CC-BY-SA | Complex SoC + DDR designs |
| Arduino (official boards) | Eagle / Altium | CC-BY-SA | Common patterns |
| Seeed Studio, Pine64, Pimoroni, Waveshare (selected) | Various | Varies; check each | Volume |
| Open Compute Project (DC-SCM, OAM, server boards from Meta/Microsoft) | Allegro/Altium/PDF | OCP licences (OCPHL) | High-speed, power-delivery precedent at the big-guys level |
| KiCad demo projects, KiCad "Made with KiCad" | KiCad | Various | Clean KiCad idioms |
| OSHWA certified projects directory | Various | Open hardware licences | Discovery index |
| Vendor eval boards and reference designs: TI (TIDA, EVM design files), ADI (CN-0xxx, EVAL boards), NXP, ST (Nucleo/Discovery CAD), Microchip, Espressif (DevKit schematics, some KiCad), Nordic (nRF DK), NVIDIA Jetson Orin Nano/NX carrier reference (design files via developer account), AMD/Xilinx eval board schematics (PDF), Intel platform schematics where public | Altium / Allegro / KiCad / PDF | **Vendor-specific licences; TI's design-file licence and NVIDIA's developer terms restrict redistribution and sometimes derivative use. Record and respect each.** | The "big guys" precedent. Schematic PDFs without CAD still yield netlists via PDF schematic parsing (Phase 4 stretch) |
| Your own past boards | Altium / KiCad | Yours | The most important corpus of all |

Eagle and Altium files are imported through KiCad's importers (`kicad-cli`
and the GUI importers; see docs/04) so that one parser (KiCad s-expression)
covers the whole corpus. Allegro files are not importable; use their PDFs.

## 3.4 Fab and assembly capability

JLCPCB, PCBWay, Eurocircuits, Aisler, OSH Park, Würth, Multi-CB, NextPCB,
Seeed Fusion: capability pages, stackup tables, impedance calculators,
assembly part libraries, design-for-assembly guides. All public; most change
quarterly, so each record carries `verified_at` and the crawler re-checks
monthly. Community-maintained KiCad DRU files for JLCPCB exist and are a
starting point but must be re-validated against the current capability page.

## 3.5 Standards-derived rule encoding (how to use paid standards legally)

For each standard you own: the implementer (Claude Code, locally, from your
licensed copy) produces YAML rules of the form
`{id, scope, params, provenance: {std: "IPC-2221B", clause: "6.3"}}`.
No quoted text beyond a few words for the title. The rule file is yours
(it is your application of the standard), the standard text stays on your
shelf. This is the same thing an engineer does when they write a company
design-rule document.

## 3.6 Sizing

Start narrow and widen: ingest datasheets and corpus for (a) every part in
your own BOMs, (b) the top ~5,000 parts by Octopart/Digi-Key popularity in
your categories, (c) every part used in the open corpus you pull. That is
on the order of 20–50 k documents and 2–5 k designs, roughly 100–300 GB raw.
Full-catalogue ingestion (millions of datasheets) is a later, optional step
and is only worth it for parametric search, which the distributor APIs
already provide live.
