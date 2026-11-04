# ADR-0004 — Provenance on every record, rule and finding

Status: accepted

## Context
The owner must be able to trust a finding without re-deriving it. A finding
that cannot be traced to a datasheet page, a rule with a clause, a corpus
statistic or a simulation run is an opinion.

## Decision
`provenance[]` is a required field on documents, extracted records, rules,
corpus entries and findings. Reports demote unprovenanced findings to an
appendix. The critic pass verifies that cited sources support the claim.

## Consequences
- Ingestion keeps page numbers and section paths; chunking never discards
  them.
- Licensed standards are cited by clause, never quoted.
- Slightly more storage and plumbing; large gain in auditability and in the
  measured hallucination rate.
