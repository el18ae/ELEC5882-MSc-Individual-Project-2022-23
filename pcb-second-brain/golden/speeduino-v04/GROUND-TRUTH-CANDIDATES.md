# Ground-truth candidates — speeduino-v04

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/speeduino/Hardware
- Pinned commit: `ee6a9242834891877744ff2137663b6ca8d79ed6` (2026-04-04T17:33:26+11:00)
- Paths analysed: `v0.4`
- Domain: automotive
- Commits touching these paths: 8; touching design files: 6; design-file commits whose message matches a fix keyword: 1

## Revision pairs

Consecutive revision directories present at the pinned commit (diff each pair):

- `v0.4/SMD/Prior Versions/`: `0.4.4b` → `0.4.4c`
  - pair `0.4.4b` → `0.4.4c`
- Also `v0.4/SMD/Prior Versions/0.4.4c` → `v0.4/SMD/Latest` (v0.4.4d), and `v0.4/THT/Prior versions/v0.4.3` → `v0.4/THT/Latest` (v0.4.3d).

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [1](https://github.com/speeduino/Hardware/issues/1) | closed | 2021-11-01 | 0.4.3d 22uf ceramic capacitor in reset protection |
| [5](https://github.com/speeduino/Hardware/issues/5) | open | 2021-12-21 | Issues or not ? On schematic v0.4.3.png |
| [6](https://github.com/speeduino/Hardware/issues/6) | open | 2022-05-12 | Latest na6 pnp board schematics? |
| [8](https://github.com/speeduino/Hardware/issues/8) | closed | 2022-11-14 | Resistors 25 27 31 32 |
| [9](https://github.com/speeduino/Hardware/issues/9) | open | 2022-11-22 | Hardware versions schematics & PCB images |
| [12](https://github.com/speeduino/Hardware/issues/12) | open | 2023-09-20 | v0.4.4.d BoM update needed |
| [15](https://github.com/speeduino/Hardware/issues/15) | open | 2024-02-01 | Gerber outline (gml/gko/gbr) file mising |
| [16](https://github.com/speeduino/Hardware/issues/16) | closed | 2024-04-05 | 0.3 Schematic PDF |
| [20](https://github.com/speeduino/Hardware/issues/20) | open | 2024-11-16 | Missing v0.4 SMD symbols after importing in KiCAD 8 |
| [21](https://github.com/speeduino/Hardware/issues/21) | closed | 2024-11-30 | Unable to view A17 and A18 in TunerStudio for DropBear |
| [23](https://github.com/speeduino/Hardware/issues/23) | open | 2025-08-23 | Best drop-in alternative for STP75NS04Z (clamped MOSFET) |
| [24](https://github.com/speeduino/Hardware/issues/24) | closed | 2026-02-23 | Dropbear v2 PCB ground point |
| [25](https://github.com/speeduino/Hardware/issues/25) | open | 2026-09-23 | v0.4.3d U8: Value and MPN name different regulators (LM2940 LDO vs NCV7805) |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`21e9fde79`](https://github.com/speeduino/Hardware/commit/21e9fde796021f3dbbc06fd6b149187b6d01ed6d) | 2021-06-23 | generic fix | Minor fix for v0.4.4c | v0.4.4c.kicad_pcb |
