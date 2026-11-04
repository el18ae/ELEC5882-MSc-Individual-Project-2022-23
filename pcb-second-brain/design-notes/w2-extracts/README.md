# W2 source extracts (KiCad)

Raw, cited extracts that `schemas/rule-taxonomy.md` was built from. They are kept
because W3 (DRC engine, rule evaluator, ERC matrix design notes) needs the same
facts. The Altium, Eagle and golden-board extracts are not kept: their content
is in the taxonomy, and the Altium extract cited sources that D6 rules out
(docs/11 Q26).

| File | Content | Source |
|---|---|---|
| `A1-kicad-dru.md` | `.kicad_dru` language: 42 constraint types, clauses, resolution order, 34 expression functions, 343 item properties, severities, 9.0/10.0 differences | KiCad GitLab `4c23271608ae7fa9e584c6f63b57e54163832c48` (2026-11-03), file:line on every row |
| `A2-kicad-setup-erc.md` | Board-setup constraints (30 JSON keys), implicit rules, netclass fields, 31 DRC test providers, 78 DRC codes, 63 ERC codes, the 12 × 12 ERC pin matrix | same commit |

Paths marked † were fetched read-only from GitLab at the same commit, because
the sparse checkout used for W2 did not contain them.
