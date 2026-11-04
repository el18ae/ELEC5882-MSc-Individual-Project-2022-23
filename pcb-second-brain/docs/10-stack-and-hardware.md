# 10 — Stack, hardware, models, cost

## Software stack (local machine)

| Concern | Choice | Why |
|---|---|---|
| Language | Python 3.12+, `uv`, `pydantic` v2, `ruff`, `mypy` | Ecosystem of every EDA/PDF/ML tool used here |
| Database | PostgreSQL 16 + `pgvector` (+ `pg_trgm`, FTS) in Docker | One engine for relational, vector and text search; recursive CTEs cover graph queries at this scale. Reuse whatever your two existing second brains use if it is Postgres-compatible (see docs/11 Q1) |
| Analytics | DuckDB over Parquet exports of corpus stats | Fast layout-statistics queries |
| Graph (in-memory) | NetworkX; `rustworkx` if slow | Motif mining, return-path analysis |
| PDF parsing | Docling (primary), pdfplumber (fallback), Tesseract OCR | Layout-aware tables and section structure |
| EDA | KiCad 10 (`kicad-cli`, IPC API + `kipy`), `kiutils`/`sexpdata` for offline parsing, `kicad-sch-api` for schematic file edits | Canonical substrate; headless CLI; live board editing via IPC. Schematic editing via IPC arrives in KiCad 11; plan for the upgrade |
| Schematic-as-code | atopile (primary), SKiDL (fallback) | Typed IR → KiCad |
| Routing | Freerouting (open); Quilter/DeepPCB optional behind the same interface | Open first; commercial if evals justify |
| Simulation | ngspice; `scikit-rf`; openEMS optional | Free, scriptable |
| Embeddings | Voyage AI (Anthropic's recommended embedding partner) for quality, or a local open model (BGE-M3 / nomic-embed) if you want zero cloud for documents. Decide in docs/11 Q6 | Datasheet prose is domain-heavy; evaluate both on a 200-query set in Phase 3 |
| LLM | Anthropic API: Opus 5.5 (default agent), Fable 5.1 (critic, hard generation), Sonnet 5.5 / Haiku 5.5 (extraction via Batches) | See CLAUDE.md § Model usage |
| Agent harness | Claude Code locally with the MCP server (Phases 0–7); Claude Agent SDK or Managed Agents for scheduled/unattended runs (Phase 8) | Same tools either way |
| MCP server | Official MCP Python SDK | Namespaced tools (docs/05) |
| Observability | Run ledger in Postgres; OpenTelemetry optional | Reproducibility |
| Altium | Your licensed Altium Designer for final export/verification; KiCad's Altium importers for ingestion. Altium's own KiCad importer for the reverse path (verify on your version) | Altium stays the editor you may still open for clients |

## Hardware (recommended minimum for the local build)

- CPU: 8+ cores (parsing, simulation, routing are CPU-bound).
- RAM: 64 GB (Docling on large PDFs, NetworkX on big designs, Postgres).
- Storage: 2 TB NVMe for `data/` (raw ≈ 100–300 GB for the scoped corpus,
  parsed JSON ≈ 2×, Postgres ≈ 50–100 GB, headroom for growth), plus
  external backup of `data/raw` and a nightly `pg_dump`.
- GPU: optional. Only needed for local embeddings/OCR speed. A 12–24 GB
  card makes local embedding of 50 k documents an overnight job rather than
  a week.
- OS: Linux preferred (KiCad headless + Docker + Freerouting are smoothest);
  Windows with WSL2 works if Altium must coexist on the same machine.

## Cost model (API)

Formula per phase, so you can plug your own scope:

- Extraction: `documents × pages × ~900 tokens/page` input at Haiku batch
  rates (`$0.10 / MTok` list, halved by Batches) plus structured output
  tokens. 50 k datasheets × 40 pages ≈ 1.8 B tokens ≈ $90 input at batch
  Haiku; Sonnet re-passes on ~15% add ≈ $300. Order of magnitude: low
  hundreds of dollars for the scoped corpus.
- Review: a 300-component board context pack ≈ 150–300 k tokens; Opus 5.5
  review + Fable critic ≈ $5–15 per full review at list prices, less with
  prompt caching on the stable KB preamble.
- Generation/layout loops: bounded by the repair-loop cap; budget ≈ 2–4× a
  review per block.

Enable prompt caching on the system prompt + rule preamble; use the Batches
API for everything offline; keep the per-review cost in the run ledger.

## Backups and privacy

- `data/raw`, `golden/`, Postgres: nightly local backup + weekly off-site
  encrypted copy.
- Owner boards are sent only in API calls to Anthropic; nothing is uploaded
  to distributor APIs beyond MPNs. If any client board must never leave the
  machine, run the review with T0–T4 only (no T5) and local embeddings.
