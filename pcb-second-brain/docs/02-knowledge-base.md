# 02 — Knowledge base

Six stores, one database engine (PostgreSQL 16+ with `pgvector`), raw files on
disk. Graph questions are answered with recursive SQL and NetworkX in memory;
a dedicated graph database is deliberately deferred (see docs/10).

## 2.1 Parts store (relational)

Purpose: the single source of truth for "does this part exist, what is it,
what are its pins, can I buy it".

| Table | Key columns |
|---|---|
| `manufacturer` | id, name, aliases[], website |
| `part` | id, mpn, manufacturer_id, category (taxonomy path), description, package_id, lifecycle (active/nrnd/eol/obsolete/unknown), lifecycle_source, rohs, reach, aec_q, datasheet_doc_id, parametrics (jsonb), updated_at |
| `package` | id, name (e.g. QFN-48 7x7 0.5mm), jedec_ref, ipc7351_name, pitch, body dims, thermal pad |
| `pin` | part_id, number, name, type (power_in, power_out, io, passive, nc, ...), bank/domain, notes, source_page |
| `symbol`, `footprint` | part_id, library (kicad lib ref / generated), checksum, verified_by (ipc-generator / vendor / human), verification_date |
| `alternate` | part_id, alt_part_id, kind (drop-in / functional / pin-compatible), notes, source |
| `offer` | part_id, distributor, sku, region, stock, moq, price_breaks (jsonb), lead_days, fetched_at |
| `assembler_part` | part_id, assembler (jlcpcb, pcbway, ...), assembler_pn, class (basic/extended), fetched_at |

Rules:
- A part may only be referenced by a design if `pin` rows exist and the
  footprint is `verified_by != none`. Otherwise it is "stub" and triggers a
  finding.
- Offers are snapshots with timestamps; "available worldwide" is computed at
  query time from the latest snapshot across regions (docs/05 § sourcing).

## 2.2 Document store + vector index

Purpose: everything written in prose that the system may need to cite.

| Table | Key columns |
|---|---|
| `document` | id, kind (datasheet / app_note / design_guide / fab_capability / standard_summary / errata / forum_thread / book_chapter), title, publisher, url, local_path, sha256, licence, fetched_at, part_ids[] |
| `section` | doc_id, path (e.g. "8.2.3 Layout Guidelines"), page_start, page_end, text, tables (jsonb), figures[] |
| `chunk` | section_id, text, embedding (vector), token_count |

Chunks are section-aware (never split a table) and carry page numbers so a
citation is `doc_id + page`. Hybrid retrieval: BM25 (Postgres full text) +
vector, reranked.

## 2.3 Design corpus (graph-in-relational)

Purpose: precedent. "How do competent designers wire and lay out this part?"

| Table | Key columns |
|---|---|
| `project` | id, name, source_repo, licence, eda_tool, version, fetched_at, quality_score, tags[] |
| `d_component` | project_id, refdes, part_id (nullable if unresolved), value, footprint, sheet, x, y, rot, layer |
| `d_net` | project_id, name, class (power / gnd / signal / diffpair / clock / analog), width_hint |
| `d_pin_net` | project_id, refdes, pin, net_id |
| `d_stackup` | project_id, layers (jsonb), thickness, cu weights, impedance table |
| `d_track_stats` | project_id, net_id, total_len, min_width, via_count, layers_used |
| `motif` | id, anchor_part_id, kind (decoupling / crystal / usb / ldo / buck / level_shift / ...), subgraph (jsonb), occurrences, canonical_example_project_id |

Motifs are mined after ingestion: for each anchor part, take the 1–2-hop
neighbourhood in the pin-net graph, canonicalise it, and count occurrences.
A motif seen in 40 of 45 designs using a part is strong precedent; its absence
in a reviewed design is a T3 finding.

## 2.4 Rules store (the guard rails)

Rules are YAML files in git (`rules/`), loaded into the DB for querying.

```yaml
id: TPS5430.bootstrap_cap
scope: {part: "TPS5430"}                    # or {package: ...}, {interface: usb2}, {fab: jlcpcb}, {std: "IPC-2221B"}
severity: error                              # error | warning | info
title: Bootstrap capacitor required between BOOT and PH
check: pcbbrain_checks.rules.cap_between_pins
params: {pin_a: BOOT, pin_b: PH, min_nF: 10, max_nF: 100}
provenance:
  - {doc: "sluS7B9", section: "8.2.2.3", page: 14}
confidence: high                             # high = human-reviewed; medium = extracted, consistent across sources; low = single extraction
```

Rule classes:
- **Part rules** (from datasheets): required external components, pin
  strapping, decoupling, thermal pad connection, exposed pad vias.
- **Interface rules** (USB, Ethernet, DDR, I2C, SPI, CAN, RS-485, LVDS):
  termination, impedance, length matching, ESD placement.
- **Fab/assembly rules** (per manufacturer): compiled to `.kicad_dru` plus
  non-DRU constraints (silkscreen over pads, panelisation, fiducials, part
  orientation for assembly).
- **EMC/SI hygiene rules** (from open design guides and clause references of
  paid standards): return-path continuity under high-speed traces, stitching
  via spacing, loop area for switching nodes, split-plane crossings.
- **Project/style rules** (yours): naming, sheet structure, refdes
  conventions, preferred passives series.

## 2.5 Fab capability store

Per manufacturer × process class (e.g. JLCPCB 2-layer standard, JLCPCB 4-layer
JLC04161H-7628, PCBWay standard, Eurocircuits, Aisler, OSH Park):
min trace/space, min drill, annular ring, via-in-pad allowed, solder mask
dam, silkscreen min width, max board size, available stackups with
impedance per geometry, copper weights, surface finishes, assembly
constraints. Source: capability pages, scraped and human-checked, with a
`verified_at` date and the URL. Exported as `.kicad_dru` + stackup JSON.

## 2.6 Project memory (yours)

Your past boards (parsed into the design corpus with `tags: [owner]`), past
review reports and what you accepted or rejected, preferred parts and series,
preferred fabs, naming conventions. This is what makes the system behave like
*you* rather than like the average of the internet. Stored in the same
tables with an `owner` flag; never leaves the machine except in API calls you
authorise.

## Query patterns the agent needs (these become MCP tools)

- `kb.part.resolve(mpn | description | parametric filter)` → PartRecord(s)
- `kb.part.pins(part_id)`, `kb.part.alternates(part_id, constraints)`
- `kb.doc.search(query, filters: part_id | kind | publisher)` → cited chunks
- `kb.doc.section(doc_id, path)` → full section with tables
- `kb.rules.for(part_ids[], interfaces[], fab)` → applicable rules
- `kb.corpus.designs_using(part_id)`, `kb.corpus.motifs(part_id)`
- `kb.corpus.layout_stats(part_id | net_class)` → distributions (track widths, via counts, distances)
- `kb.fab.capabilities(fab, process)` → constraints + DRU path
- `kb.memory.preferences()`, `kb.memory.past_findings(board)`
