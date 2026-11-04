# 04 — Ingestion pipeline ("ingest and hydrate")

Five stages, each idempotent, each writing provenance. Re-running a stage on
an already-processed item is a no-op unless the source changed (sha256).

```
acquire → parse → normalise → extract → index
  raw/      parsed/   DesignIR /     rules,     pgvector,
  (files)   (JSON)    PartRecord     motifs,    FTS, graph
                                     pins       tables
```

## 4.1 Acquire

- One crawler per source family (`pcbbrain-ingest/crawlers/*`), each with:
  licence check (hard-coded allow list from docs/03), rate limiting, ETag/
  Last-Modified caching, a manifest (`url, sha256, fetched_at, licence_note`).
- Git-based corpora are cloned shallow; only CAD + licence + README files
  are retained.
- API-based data (offers, lifecycle) is not "acquired" in bulk; it is
  fetched on demand and cached with a TTL (docs/05 § sourcing).
- Raw files live in `data/raw/<source>/<sha256-prefix>/...`. Never edited.

## 4.2 Parse

| Input | Parser | Output |
|---|---|---|
| PDF datasheet / app note | Docling (layout-aware, table extraction) → fallback pdfplumber → OCR (Tesseract) for scanned pages | `section[]` with page ranges, tables as row/col JSON, figure crops |
| KiCad project (`.kicad_pro/.kicad_sch/.kicad_pcb`) | s-expression parser (`kiutils` or `sexpdata`); `kicad-cli sch export netlist`; `kicad-cli pcb export` for stats | DesignIR |
| Altium (`.SchDoc/.PcbDoc/.PrjPcb`) | KiCad importer (headless: `kicad-cli` where supported, otherwise a scripted GUI import run once per project in a batch session) → KiCad project → above | DesignIR, flagged `imported_from: altium` |
| Eagle (`.sch/.brd`) | KiCad importer | DesignIR, flagged |
| Schematic-only PDFs (vendor eval boards) | Phase 4 stretch: vision model (Claude) page-by-page → netlist candidates → cross-checked against part pin tables; `confidence: low` | DesignIR (partial) |
| Fab capability pages | HTML → table extraction → YAML (human-verified once) | `fabs/*.yaml` |

## 4.3 Normalise

- MPN normalisation (strip packaging suffixes like `-ND`, `T&R`, `/NOPB`
  into a `packaging` field; keep the orderable code as `sku`).
- Manufacturer alias resolution (e.g. "Texas Instruments" / "TI" / "National
  Semiconductor" legacy parts).
- Footprint name canonicalisation to IPC-7351 style; map vendor names.
- Net class inference in DesignIR: power (name regex + connected to
  regulator outputs), ground, diff pairs (`_P/_N`, `+/-`), clocks, analog.

## 4.4 Extract (where Claude does the reading)

Batch jobs via the Message Batches API (50% cost), structured outputs with
strict schemas, model tiering:

| Extraction | Model | Schema | Validation |
|---|---|---|---|
| Pin table from datasheet | Haiku 5.5 → Sonnet 5.5 on low confidence | `pin[]` | Pin count must match package; power pins must exist; cross-check with vendor CAD if available |
| Parametrics (voltage ranges, currents, temp, package) | Haiku 5.5 | `parametrics` | Range sanity; agreement with distributor parametrics |
| Rule candidates from "Typical application", "Layout guidelines", "Recommended operating conditions", app notes | Sonnet 5.5; Opus 5.5 for complex power/RF parts | `rule_candidate[]` (see schemas/rule.schema.json) | Must cite page; must reference pins that exist; duplicates merged; promoted to `confidence: medium` when two sources agree, `high` only after human or eval confirmation |
| Motif mining from corpus | No LLM: graph canonicalisation code | `motif[]` | Occurrence counts |
| Motif *naming/explanation* | Haiku 5.5 | label + one-line rationale | Cosmetic only |

Every extracted record carries `{doc_id, page, model, prompt_version,
extracted_at}`. Re-extraction is triggered when `prompt_version` changes.

## 4.5 Index

- Chunk sections (≤ 800 tokens, table-atomic), embed (see docs/10 for
  model), store in `chunk.embedding`; FTS index on `section.text`.
- Materialised views: `v_part_summary` (part + pin count + rule count +
  corpus usage count + latest availability), `v_motif_frequency`.
- A nightly job refreshes availability for every part in any owner board and
  in the top-N popular list, and records lifecycle transitions as events
  (an EOL transition on a part in one of your boards generates a standing
  finding).

## 4.6 Quality gates for the pipeline itself

- Parse coverage: % of PDFs with ≥ 1 table extracted, % with a pin table
  that matches package pin count. Target ≥ 90% on the top-5k set.
- Corpus resolve rate: % of BOM lines in ingested designs that resolve to a
  `part` row. Target ≥ 70% (open designs use many generic passives).
- Rule precision sample: a random 50 extracted rules per week are reviewed
  (by you, 20 minutes) → feeds the confidence promotion.
