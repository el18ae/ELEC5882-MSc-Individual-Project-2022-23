# 10 — Stack, hardware, models, cost

## Platform: Windows 11 (owner's machine, KiCad and Altium installed)

| Concern | Choice | Notes |
|---|---|---|
| Language | Python 3.12+ via `uv`; `pydantic` v2; `ruff`; `mypy --strict`; Rust via `pyo3`/`maturin` only for kernel hot paths that profiling proves (geometry booleans, DRC spatial queries, router) | Native Windows Python; no WSL requirement for the kernel |
| Database | PostgreSQL 16 + `pgvector` + `pg_trgm` in Docker Desktop (WSL2 backend) | One engine for relational, vector, text; recursive CTEs for graph queries |
| Analytics | DuckDB over Parquet | Corpus statistics |
| Graph (in-memory) | NetworkX, `rustworkx` if slow | Motifs, return-path analysis |
| Geometry | `shapely` 2 + `rtree`/STRtree; `numpy` | Kernel `geom` |
| Altium file reading | `olefile` for OLE compound documents + own record parsers (from docs/13 notes) | No Altium automation required |
| PDF parsing | Docling (primary), pdfplumber (fallback), Tesseract OCR | Datasheets, app notes, Altium/KiCad help PDFs |
| Oracle | KiCad 10 `kicad-cli.exe` from the local install (tests only) | ADR-0001 |
| External router (interim) | Freerouting jar (Java 21) via DSN/SES | Until own `route` lands |
| Simulation | ngspice (Windows build), `scikit-rf`; openEMS optional | |
| STEP export | `cadquery` / OCP | Phase 7 |
| Embeddings | Voyage AI or local BGE-M3 (`sentence-transformers`); decide in docs/11 Q6; evaluate both on a 200-query set in Phase 3 | |
| LLM | Anthropic API: Opus 5.5 (agent), Fable 5.1 (critic, hard generation, source-reading for design notes), Sonnet 5.5 / Haiku 5.5 (extraction via Batches) | CLAUDE.md § Model usage |
| Agent harness | Claude Code on Windows with the MCP server (Phases 0–7); Agent SDK / Managed Agents for unattended runs (Phase 8) | |
| MCP server | Official MCP Python SDK, stdio transport locally | |
| Observability | Run ledger in Postgres; OpenTelemetry optional | |
| Viewing results | Kernel renders SVG/PNG; exports to KiCad format; Altium opens KiCad projects through its importer; no GUI editor is built | ADR-0001 |

Install paths are recorded in § Measured by Phase 0, typically
`C:\Program Files\KiCad\10.0\` and `C:\Program Files\Altium\AD<ver>\` with
user data under `C:\Users\Public\Documents\Altium\AD<ver>\`.

## Hardware (minimum for the local build)
- CPU 8+ cores; RAM 64 GB; 2 TB NVMe for `data/` (raw ≈ 150–400 GB for an
  all-domain scoped corpus, parsed ≈ 2×, Postgres ≈ 100 GB); GPU optional
  (local embeddings/OCR speed). External backup for `data/raw`, `golden/`
  and nightly `pg_dump`.

## Measured (filled by Phase 0)
- CPU / RAM / GPU / free disk / Windows build:
- KiCad version and path:
- Altium version and path:
- Docker Desktop / WSL2 version:
- Java version (Freerouting):

## Cost model (API)
- Extraction: `documents × pages × ~900 tokens/page` at batch Haiku
  rates; an all-domain scoped corpus of ~80 k documents × 40 pages ≈ 2.9 B
  tokens ≈ $150 input at batch Haiku; Sonnet re-passes on ~15% ≈ $450.
- Design notes from source (K1): KiCad source is large; reading the
  relevant subsystems with Fable 5.1 at list prices is on the order of
  $200–600 total, done once per KiCad major version.
- Review: 300-component board ≈ 150–300 k tokens of context; Opus review +
  Fable critic ≈ $5–15 per full review; less with prompt caching.
- Generation/layout loops: bounded by the repair cap; ≈ 2–4× a review per block.

## Backups and privacy
Nightly local backup, weekly encrypted off-site copy. No owner boards
exist yet; when they do, docs/11 Q3 governs whether they may be sent to the
Anthropic API. Public corpus and tool knowledge may be embedded via a
cloud API unless Q6 says otherwise.

> Owner answer Q10 (2026-11-04): budget is the $200/month Claude Ultra subscription; API spend is $0 unless approved. The cost model above assumes API billing and must be revised in W10: build and run through Claude Code under the subscription, and keep the API path as an option.
