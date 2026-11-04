# Ground-truth candidates — moteus

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/mjbots/moteus
- Pinned commit: `aec03338883359fae7b32df6d56974958448c67e` (2026-11-03T14:31:27-05:00)
- Paths analysed: `hw/c1`, `hw/controller`
- Domain: motor_drive
- Commits touching these paths: 127; touching design files: 102; design-file commits whose message matches a fix keyword: 48

## Revision pairs

Consecutive revision directories present at the pinned commit (diff each pair):

- `hw/c1/`: `r1.0` → `r1.1` → `r1.2`
  - pair `r1.0` → `r1.1`
  - pair `r1.1` → `r1.2`
- `hw/controller/`: `r4.3` → `r4.5` → `r4.8` → `r4.11`
  - pair `r4.3` → `r4.5`
  - pair `r4.5` → `r4.8`
  - pair `r4.8` → `r4.11`

Tags in repository (77): `0.1-20260428`, `0.1-20260430`, `0.1-20260505`, `firmware/v1.0.0-rc.1`, `firmware/v1.0.0-rc.2`, `python/v1.0.0rc1`, `cpp/v1.0.0-rc.1`, `python/v1.0.0rc2`, `python/v1.0.0rc3`, `python/v1.0.0rc4`, `firmware/v1.0.0`, `python/v1.0.0`, `cpp/v1.0.0`, `rust/v0.5.1`, `rust/v0.5.2`, `rust/v0.5.3`, `firmware/v1.1.0-rc.1`, `python/v1.1.0rc1`, `firmware/v1.1.0`, `python/v1.1.0`, `firmware/v1.1.1`, `python/v1.1.1`, `firmware/v1.1.2`, `rust/v0.5.4`, `rust/v0.5.5` (last 25)

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [3](https://github.com/mjbots/moteus/issues/3) | closed | 2020-10-18 | Schematics are empty |
| [49](https://github.com/mjbots/moteus/issues/49) | closed | 2023-11-04 | Cross-platform Builds |
| [68](https://github.com/mjbots/moteus/issues/68) | closed | 2024-08-15 | Can Moteus work at -45 ° C-100 ° C |
| [79](https://github.com/mjbots/moteus/issues/79) | closed | 2025-01-04 | Question regarding C1 schematic |
| [87](https://github.com/mjbots/moteus/issues/87) | closed | 2025-05-13 | I2C AS5600 Encoder and Digital Inputs |
| [97](https://github.com/mjbots/moteus/issues/97) | closed | 2025-08-15 | GPIO pins won't enable |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`0d4217bb4`](https://github.com/mjbots/moteus/commit/0d4217bb44041fda26bafba41db50f11fc93c064) | 2022-04-29 | footprint, part swap, power, generic fix | moteus r4.11 * Remove snubbers * Switch gate drive resistors to 7.5 ohm * Turn off paste on the SOT23-8 TCAN f | controller.brd, controller.lbr, controller.sch |
| [`f8591bb49`](https://github.com/mjbots/moteus/commit/f8591bb49bd5e607371ac771abcfd8a6caebab83) | 2022-04-29 | decoupling, footprint, part swap, power, routing | moteus r4.9  * Add an alternate SOIC-8 footprint to the CAN transceiver  * Increase gate drive traces to 15 mi | controller.brd, controller.lbr, controller.sch |
| [`c6526177b`](https://github.com/mjbots/moteus/commit/c6526177b7ac0afbbd4d479a2871bea1433dbeec) | 2021-09-09 | decoupling | Add back in attributes for the cap that lost its border | controller.brd, controller.sch |
| [`f68d46e5e`](https://github.com/mjbots/moteus/commit/f68d46e5e1cdb8bf119f3721b0cae6052f37e811) | 2021-09-09 | protection, generic fix | Fix TVS SOT23-3 layout | controller.brd, controller.lbr, controller.sch |
| [`5b19e46e1`](https://github.com/mjbots/moteus/commit/5b19e46e1bb3a47c50d1b0b359ae64b0c3d746af) | 2021-09-09 | clearance, generic fix | Fix DRC violations | controller.brd |
| [`c916c20f9`](https://github.com/mjbots/moteus/commit/c916c20f97f480de0c0dda2a9a5c28cca0a013fb) | 2021-09-09 | routing | More work on ground planes | controller.brd |
| [`238c9b4f8`](https://github.com/mjbots/moteus/commit/238c9b4f868499b53837aaa6a32daa5a405f6b0b) | 2021-09-09 | decoupling | Add a vsense filter capacitor, initially tuned for 1kHz cutoff | controller.brd, controller.sch |
| [`fef1c17d4`](https://github.com/mjbots/moteus/commit/fef1c17d43f64619d37c90b988bd2cca84b98c30) | 2021-09-09 | silkscreen | r4.6: Yikes, more silk tweaks | controller.brd |
| [`7d22e3a54`](https://github.com/mjbots/moteus/commit/7d22e3a54e2ac811e8234b0f1abd6b0744f11b32) | 2021-09-09 | routing | r4.6: More trace and poly tweaks | controller.brd |
| [`e18af2a7e`](https://github.com/mjbots/moteus/commit/e18af2a7e6f153f30f5fa5a910b8fc2f65125349) | 2021-09-09 | silkscreen | Tweak silk positions | controller.brd |
| [`dee5aa80a`](https://github.com/mjbots/moteus/commit/dee5aa80aa0315e10fcbaf90146255ce63e8acc2) | 2021-09-09 | footprint, silkscreen, generic fix | Scoot debug pads around to give more room for silks | controller.brd, controller.sch |
| [`081641a2e`](https://github.com/mjbots/moteus/commit/081641a2eaaeb49c95725035ea1e110b08da38bc) | 2021-09-09 | routing | r4.6: And actually connect the AGND to the logic ground | controller.brd, controller.sch |
| [`85b189cce`](https://github.com/mjbots/moteus/commit/85b189cce6e03a5b01f9e0377935b0c3705567de) | 2021-09-09 | routing | r4.6: Connect the DRV8353F gnd to the power ground | controller.brd, controller.sch |
| [`8183027fb`](https://github.com/mjbots/moteus/commit/8183027fbf46955307ff554b7a4d167b9f896d11) | 2021-09-09 | power | r4.6: More regulator layout tweaking | controller.brd |
| [`aac9bb13a`](https://github.com/mjbots/moteus/commit/aac9bb13aad2ab2d58ddf47fc7effcd7af032762) | 2021-09-09 | footprint, protection, power, routing, generic fix | Beginnings of r4.6  * Fix grounding around MAX15062AATA+T  * Add a TVS diode on the CAN lines  * Add test pads | controller.brd, controller.lbr, controller.sch |
| [`a2aea544f`](https://github.com/mjbots/moteus/commit/a2aea544f9339bb7c7ee33afce8562a0a37d44e4) | 2020-09-15 | clearance | Add test through holes for the CAN bus on a 2mm spacing | controller.brd |
| [`4852c8644`](https://github.com/mjbots/moteus/commit/4852c8644a8b62a4e7d928d0f842b3cb0500b4dc) | 2020-04-11 | silkscreen | Another silk tweak | controller.brd |
| [`7036afa25`](https://github.com/mjbots/moteus/commit/7036afa25eecc7e10fe23d529c06e8d1ed9d7446) | 2020-04-11 | silkscreen | More silk and layout tweaks | controller.brd |
| [`cf6515e51`](https://github.com/mjbots/moteus/commit/cf6515e51dfba4055eb503bc5e14f8e63286e61b) | 2020-04-09 | silkscreen, routing | Tweak some silks and traces | controller.brd |
| [`3895c1504`](https://github.com/mjbots/moteus/commit/3895c15048f0cbd9fba24f5bf4b5ca2bcf074331) | 2020-04-07 | routing | r4.4 autorouted | controller.brd, controller.sch |
| [`1d600dadc`](https://github.com/mjbots/moteus/commit/1d600dadc9a7c8f7febe2c4f487f751aaf03ac06) | 2020-02-13 | routing | Tweak the plane around the B phase | controller.brd |
| [`8f02fe5a1`](https://github.com/mjbots/moteus/commit/8f02fe5a10b91531f209d3b2481a3fc04a0456ed) | 2020-01-22 | silkscreen | Start on r4.3 with updated silk and moved phase outputs | controller.brd, controller.sch |
| [`ef3f0f307`](https://github.com/mjbots/moteus/commit/ef3f0f307fa3836953730d81388f4be16811c141) | 2019-11-04 | footprint, silkscreen | controller: Sigh, I managed to tweak that silk onto a pad | controller.brd |
| [`12236f7fc`](https://github.com/mjbots/moteus/commit/12236f7fcc967da0dcc8a9647065b8b5617073f5) | 2019-11-04 | silkscreen | controller: Tweak another silk position | controller.brd |
| [`268c2ec41`](https://github.com/mjbots/moteus/commit/268c2ec4194c80892a3db39f8c1543a0fdce28d6) | 2019-11-04 | pull-up/down | controller: Actually use a 2k resistor for I2C pullup | controller.brd, controller.sch |
| [`7dd33937a`](https://github.com/mjbots/moteus/commit/7dd33937a08dd1b481f8db55199967e2736434ac) | 2019-11-04 | silkscreen | controller: Tweak a silk position | controller.brd |
| [`05c45948a`](https://github.com/mjbots/moteus/commit/05c45948a41e6cb2e55e13ced557753abdfc47ae) | 2019-09-30 | routing, generic fix | Add some more ground debug pins | controller.brd, controller.sch |
| [`7ec08d918`](https://github.com/mjbots/moteus/commit/7ec08d918b90de573e300ff741c9ccc01f7f0197) | 2019-09-27 | footprint, generic fix | Added some more debugging pads | controller.brd, controller.sch |
| [`c7bcb7c3b`](https://github.com/mjbots/moteus/commit/c7bcb7c3b35f0029a97ee569658d8be8e34a9147) | 2019-09-27 | part swap, silkscreen, generic fix | All silks rough placed and all part numbers fixed up | controller.brd, controller.sch |
| [`3e30973c1`](https://github.com/mjbots/moteus/commit/3e30973c1cbf07429571f47053257646221e1af0) | 2019-09-27 | part swap | Add parts roughly placed and named | controller.brd, controller.sch |
| [`28752a6e0`](https://github.com/mjbots/moteus/commit/28752a6e058e0b86dc23d0f70fa877b549ce9252) | 2019-02-26 | part swap | controller: Transcribe r3.1 part selection back into eagle | controller.brd, controller.sch |
| [`029dc59a3`](https://github.com/mjbots/moteus/commit/029dc59a380474ae10071b1d27510e230129b032) | 2019-02-18 | silkscreen, generic fix | controller: Fix some minor silk issues | controller.brd |
| [`45ae1591f`](https://github.com/mjbots/moteus/commit/45ae1591f8a2c3a706dc5602e22a1c7f975f4d66) | 2019-02-06 | silkscreen | Put stm32f4 placement silk on the right layer | controller.brd, controller.lbr, controller.sch |
| [`d368b30b0`](https://github.com/mjbots/moteus/commit/d368b30b0c40bcb1eb80c08ebf92630d152ec9f2) | 2019-02-04 | silkscreen | Tweak silks | controller.brd |
| [`569f88b34`](https://github.com/mjbots/moteus/commit/569f88b34c9a6cdaeaede30a5fca1077955fe9ad) | 2018-12-19 | generic fix | Fix more errata in controller board | controller.brd, controller.sch |
| [`1f02004f1`](https://github.com/mjbots/moteus/commit/1f02004f119f67b5095b4df6ba144c7ac6f3c758) | 2018-12-19 | generic fix | Fix a number of r1 errata | controller.brd, controller.lbr, controller.sch |
| [`a9fe9689e`](https://github.com/mjbots/moteus/commit/a9fe9689eda5513296aacea724ed4ad48af9854a) | 2018-11-04 | silkscreen | YAST: Yet another silk tweak | controller.brd |
| [`f91c4689b`](https://github.com/mjbots/moteus/commit/f91c4689b4f08baa31bda0dfa8a8928fc16e5685) | 2018-11-04 | silkscreen | More silk tweaks | controller.brd |
| [`271f1c21d`](https://github.com/mjbots/moteus/commit/271f1c21d0f9446dc3f632800b23de67f2f2bdc6) | 2018-11-04 | part swap | Down to 20!  * Swapped in 3x 100k resistors for the 33.2k | controller.brd, controller.sch |
| [`642b90122`](https://github.com/mjbots/moteus/commit/642b90122467a584ce69ec75721332dfe32ab5a7) | 2018-11-04 | part swap | Reduce unique part count to 21 | controller.brd, controller.sch |
| [`cfa1836f9`](https://github.com/mjbots/moteus/commit/cfa1836f9c73a50cf62e304959079e12d6de34f0) | 2018-11-03 | decoupling | MF doesn't have a 0402 10pF cap | controller.brd, controller.sch |
| [`ee1c8b2d5`](https://github.com/mjbots/moteus/commit/ee1c8b2d5cc4d1d3876148a3b9f84b86253119fa) | 2018-11-03 | part swap | More BOM tweaks | controller.brd, controller.sch |
| [`59faa76a8`](https://github.com/mjbots/moteus/commit/59faa76a8e9bc0d50ccbe1fe386948b3e1ff2a73) | 2018-11-03 | footprint | Actually make the wire pads have soldermask stop | controller.brd, controller.lbr, controller.sch |
| [`2407f2326`](https://github.com/mjbots/moteus/commit/2407f232670b0752d7d0dab5efdb4c5ef666d631) | 2018-11-03 | silkscreen | Autorouting and final silk placement | controller.brd, controller.sch |
| [`f64eddcc5`](https://github.com/mjbots/moteus/commit/f64eddcc507d041a0b343aebe0dffebebfc2fec1) | 2018-11-02 | silkscreen | Partway through resizing silks | controller.brd |
| [`723aea90b`](https://github.com/mjbots/moteus/commit/723aea90b984f722e4e661792fc3ee7a94126ac6) | 2018-11-02 | routing | Nominally all power is routed | controller.brd, controller.sch |
| [`4fbeb3038`](https://github.com/mjbots/moteus/commit/4fbeb303861c4528a2b54ce54692d67fb430f509) | 2018-11-02 | routing | Most of the power side is now routed | controller.brd, controller.sch, shorts.lbr |
| [`04b6d65a1`](https://github.com/mjbots/moteus/commit/04b6d65a1e42806cad4729757f9161faa25b45e9) | 2018-11-02 | generic fix | Checkpoint after fixing board consistency errors | controller.brd, controller.sch |
