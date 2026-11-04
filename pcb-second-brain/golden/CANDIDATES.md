# Golden board candidates (W1)

Produced in cloud session 2 (2026-11-04) under docs/08 § 8.1. Every entry
below was checked against its live GitHub repository on that date: a
blobless clone of the default branch gave the file tree, the licence text
was read from the repository itself, tag and commit counts came from the
remote, and issue-tracker status came from the GitHub search API. Boards
that could not be confirmed this way are not in the table (§ 5 and § 6).

## 1. Summary

| Domain | Candidates | Phase 0 picks |
|---|---|---|
| power_electronics | 9 | bms-c1, eez-psu-h24005 |
| mixed_signal_instrumentation | 6 | glasgow |
| mcu_iot | 6 | esp32-evb |
| high_speed_digital | 9 | jetson-nano-baseboard |
| rf | 8 | hackrf-one |
| motor_drive | 6 | moteus, odrive-v3 |
| industrial_comms | 6 | sh-esp32 |
| intrinsically_safe_atex | **0** | — (see § 4) |
| automotive | 5 | speeduino-v04 |
| **Total** | **55** | **10** |

Eight of the nine domains have 5–9 candidates. The ninth, ATEX/IS, has
no open design that meets the criteria (§ 4). Industrial comms is thin
on 4–20 mA and HART: none of its six boards has either.

Format mix across the 55 (by the CAD column): KiCad 6+ only 13, KiCad
legacy (≤5) only 13, Eagle only 11, Altium 5, and 13 that mix formats
across revisions or modules.
That covers all three readers planned for K2 (`io.kicad`, `io.eagle`,
`io.altium`).

## 2. How the criteria were applied

1. **Licence.** docs/08 lists CC-BY, CC-BY-SA, CERN-OHL, MIT, TAPR and OCP.
   This session also accepted other open licences that allow use and
   modification: Apache-2.0, 0BSD, BSD, GPL-2.0/3.0 and Solderpad. This
   is docs/11 Q24. Non-commercial licences
   (Open Ephys, CC-BY-NC-SA) are accepted and flagged, because the project
   is personal and never distributed (D5). A repository with no licence
   at its root qualifies only if its README or a per-folder licence file
   states one. The quote in the table names the file it came from.
2. **Complete CAD.** "yes" means native schematic and board files are both
   present at the pinned commit, not just PDFs or Gerbers. KiCad 6+,
   Eagle XML and Altium files carry their own symbols and footprints.
   KiCad legacy projects are complete only when a `*-cache.lib` or
   project library is present. The two without one are marked "needs
   KiCad 4/5 libs".
3. **Proven.** Each board is a vendor product, a crowdfunded product, a
   vendor's open reference design, or a design with a documented
   community install base. The one exception is the RS-485 dongle,
   marked "proven: weak" and kept only to fill the RS-485 gap.
4. **Maintained.** The history column gives commits on the default branch,
   the year of the first commit, tagged releases, and the number of
   revision directories at HEAD (`rev*`, `r1.x`, `Rev.X`, `0.4.4b`, …).
   Revision directories give the most direct revision-diff ground truth.
5. **Domain balance.** Target 4–6 per domain; extra candidates beyond 6
   are listed so substitutes exist if vetting rejects a board.

Ranking within a domain weights, in order: revision depth (revision
directories, changelogs and tags beat raw commit count), proven status,
issue-tracker activity, licence clarity, and size for Phase 0 (prefer
boards the kernel can parse early).

## 3. Phase 0 pick (first 10)

The picks cover all eight domains that have candidates. A second board
goes to power and to motor drive, so the Eagle (EEZ) and Altium (ODrive)
readers get a test case in Phase 0 as well as KiCad.
All ten are downloaded under `golden/<slug>/src/` (git-ignored, decision
D12) with `PROVENANCE.yaml` and `GROUND-TRUTH-CANDIDATES.md`.
`golden/fetch-phase0.sh` re-fetches them at the pinned commits on the
local machine.

| # | Slug | Domain | Format | Why it is in Phase 0 | Ground-truth source | Leads (commits / issues) |
|---|---|---|---|---|---|---|
| 1 | hackrf-one | rf | KiCad 6+ | Shipped since 2014, 1,144 issues, r1–r10 revision notes | `list_of_hardware_revisions.rst` + commits | 79 / 23 |
| 2 | glasgow | mixed_signal_instrumentation | KiCad legacy → 6+ | Eight revision directories (revA→revD1) in one tree | revision dirs + revC3 notes + commits | 216 / 19 |
| 3 | moteus | motor_drive | Eagle (r4.x) + KiCad 6+ (c1) | Two revision ladders and two formats in one repo | revision dirs + commits | 48 / 6 |
| 4 | bms-c1 | power_electronics | KiCad 6+ | Recent KiCad, high-current BMS, bring-up test report in repo | v0.3 test report + commits | 9 / 21 |
| 5 | jetson-nano-baseboard | high_speed_digital | KiCad 6+ | SoM carrier with high-speed interfaces, clean Apache-2.0 licence | commits + tags | 71 / 13 |
| 6 | sh-esp32 | industrial_comms | KiCad | CAN (NMEA 2000) product, 10 tagged hardware releases | commits + tags | 108 / 1 |
| 7 | speeduino-v04 | automotive | KiCad | ECU with >1000 installs, revision directories 0.4.4b→d | revision dirs + issues | 1 / 13 |
| 8 | esp32-evb | mcu_iot | KiCad legacy → 6+ | Revisions A→L, each listing its changes (part swaps, value changes) | `HARDWARE/README.md` per-revision changes | 0 / 18 |
| 9 | odrive-v3 | motor_drive | Altium | Altium reader test case; CHANGELOG lists v3.x design changes | `v3/CHANGELOG.md` + commits | 14 / 4 |
| 10 | eez-psu-h24005 | power_electronics | Eagle XML | Eagle reader test case; 245 commits on a lab PSU | commits + issues | 54 / 19 |

Downloaded size: about 500 MB of working trees plus 190 MB of blobless
history. That is too large, and the wrong licence position, to commit
to this public repository (D12).

## 4. Intrinsically safe / ATEX: no qualifying open board

GitHub search (`"intrinsically safe" kicad`, `"intrinsically safe" pcb`,
`ATEX hardware pcb`) and two web searches found no open-licensed,
complete-CAD design for an IS or ATEX product. Nobody publishes one,
because certification belongs to a specific build assessed by a notified
body. The vendor reference designs that come closest are not open
licensed and were not verified here:

- TI TIDA-01504 (loop-powered 4–20 mA transmitter with HART). Found by
  web search only; licence and CAD availability not checked.
- ADI CN0270 (4–20 mA transmitter with HART, AD5420 + AD5700). Found by
  web search only; licence and CAD availability not checked.

Default applied (docs/11 Q23): cover ATEX with (a) **synthetic boards**
(docs/08 § 8.1, "Synthetic boards") built from IEC 60079-11 clause-referenced
rules (energy limitation, segregation distances, component derating,
Zener barrier topologies) once W8-style rules exist, and (b) vendor IS
reference designs downloaded locally under their own licences, used as
**non-golden reference material only**. The ATEX quota is met by
synthetic boards until the owner supplies or approves real ones.

## 5. Leads that could not be verified from the cloud

The session's network policy refused these hosts, so the designs are not
in the table. Check them in the local phase:

| Lead | Domain | Blocked host | Local fetch command |
|---|---|---|---|
| MNT Reform 2 motherboard (KiCad; licence to confirm) | high_speed_digital | source.mnt.re | `git clone --filter=blob:none https://source.mnt.re/reform/reform.git` |
| Raspberry Pi CM4 IO board (KiCad; licence to confirm) | high_speed_digital | datasheets.raspberrypi.com | `curl -LO https://datasheets.raspberrypi.com/cm4io/CM4IO-KiCAD.zip` |
| Raspberry Pi CM5 IO board (KiCad) | high_speed_digital | datasheets.raspberrypi.com | `curl -LO https://datasheets.raspberrypi.com/cm5/cm5io-kicad.zip` (path unverified) |
| RP2040 minimal design (KiCad) | mcu_iot | datasheets.raspberrypi.com | `curl -LO https://datasheets.raspberrypi.com/rp2040/Minimal-KiCAD.zip` |

Not attempted, to save budget: BeagleBoard (git.beagleboard.org), OHWR
(ohwr.org), Arduino official boards (docs.arduino.cc), vendor eval
boards behind logins. All are listed in docs/03 § 3.3 for the local
phase.

## 6. Rejected after inspection

| Repository | Reason |
|---|---|
| LowPowerLab/CurrentRanger | No CAD in repository (firmware, enclosure, UF2 only) |
| sqfmi/Watchy | Firmware only; hardware is in `sqfmi/watchy-hardware` (listed) |
| Seeed-Studio/OSHW-XIAO-Series | KiCad symbol/footprint library only, no boards |
| greatscottgadgets/greatfet | No CAD; `greatfet-hardware` is archived and deprecated |
| foxBMS/foxbms-2 | No CAD in repository (software and docs only); foxBMS 1 hardware listed instead |
| macchina/m2-hardware | PDFs and Gerbers only |
| damienmaguire/Tesla-Model-3-Charge-Port-Controller | PDFs only |
| meatpiHQ/WiCAN | No licence |
| rusefi/hellen* boards | No licence |
| VESC (vedderb) | No official open hardware repository found; `craigg96/VESC6_OPEN_HARDWARE` is an unproven re-implementation |
| CANable (protofusion) | Hardware repository not found on GitHub |
| Jana-Marie/OtterIron | Repository not found |
| Fescron/2ch-4-20mA-converter and similar 4–20 mA / RS-485 hobby repos | No licence, or new with no evidence of builds |

## 7. Full candidate table

All commit hashes are the default-branch HEAD on 2026-11-04. "P0" marks
the Phase 0 picks.

### power_electronics (9)

| Rank | Board | Source | Licence (quoted) | CAD format | Sch + PCB + libs | History (commits since / tags / rev dirs) | Issue tracker / errata | Proven | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **1 · P0** | **Libre Solar BMS C1 (16s/100 A, bq769x2)** | [LibreSolar/bms-c1](https://github.com/LibreSolar/bms-c1) @`0ca09706f4` | CERN-OHL-W-2.0; LICENCE: “CERN Open Hardware Licence Version 2 - Weakly Reciprocal” | KiCad 6+ | yes; libs embedded | 57 since 2021 / 0 / — | yes (24 open) | pilot builds; v0.3 test report in repo | test report lists bugs found during bring-up |
| **2 · P0** | **EEZ H24005 bench PSU** | [eez-open/psu-hw](https://github.com/eez-open/psu-hw) @`eef089135b` | TAPR-OHL-1.0; LICENSE.TXT: “The TAPR Open Hardware License Version 1.0 (May 25, 2007)” | Eagle XML + KiCad import | yes; libs embedded (Eagle) | 245 since 2015 / 4 / — | yes (5 open) | crowdfunded product | Eagle modules plus a consolidated KiCad import |
| 3 | Hat Labs SH-RPi (supercap UPS HAT) | [hatlabs/SH-RPi-hardware](https://github.com/hatlabs/SH-RPi-hardware) @`e82e3c9823` | CC-BY-4.0; LICENSE.md: “Creative Commons Attribution 4.0 International” | KiCad 6+ | yes; libs project libs present | 243 since 2018 / 16 / — | yes (2 open) | vendor product | 16 tagged releases |
| 4 | Libre Solar MPPT 2420 LC | [LibreSolar/mppt-2420-lc](https://github.com/LibreSolar/mppt-2420-lc) @`c2211beb36` | CC-BY-SA-4.0; LICENSE: “Attribution-ShareAlike 4.0 International” | KiCad legacy (5) | yes; libs cache lib present | 29 since 2016 / 0 / — | yes (11 open) | community builds |  |
| 5 | Libre Solar MPPT 1210 HUS | [LibreSolar/mppt-1210-hus](https://github.com/LibreSolar/mppt-1210-hus) @`a0c386ace7` | CERN-OHL-W-2.0; LICENCE: “CERN Open Hardware Licence Version 2 - Weakly Reciprocal” | KiCad legacy (5) | yes; libs cache lib present | 14 since 2019 / 2 / — | yes (3 open) | community builds |  |
| 6 | Libre Solar MPPT 2420 HC | [LibreSolar/mppt-2420-hc](https://github.com/LibreSolar/mppt-2420-hc) @`b1c11fee82` | CERN-OHL-W-2.0; LICENCE: “CERN Open Hardware Licence Version 2 - Weakly Reciprocal” | KiCad legacy (5) | yes; libs cache lib present | 6 since 2020 / 0 / — | yes (19 open) | community builds | short history (6 commits) |
| 7 | SparkFun Power Delivery Board USB-C (STUSB4500) | [sparkfun/Power_Delivery_Board-USB-C](https://github.com/sparkfun/Power_Delivery_Board-USB-C) @`4213498524` | CC-BY-SA-4.0; LICENSE.md: “SparkFun hardware is released under Creative Commons Share-alike 4.0 International” | Eagle XML | yes; libs embedded (Eagle) | 13 since 2019 / 0 / — | yes (0 open) | vendor product |  |
| 8 | Adafruit PowerBoost 1000C | [adafruit/Adafruit-PowerBoost-1000C](https://github.com/adafruit/Adafruit-PowerBoost-1000C) @`683ff0cdcd` | CC-BY-SA-3.0; license.txt: “Attribution-ShareAlike 3.0 Unported” | Eagle XML | yes; libs embedded (Eagle) | 4 since 2015 / 0 / — | yes (2 open) | vendor product | 4 commits |
| 9 | Libre Solar BMS 8S50 IC | [LibreSolar/bms-8s50-ic](https://github.com/LibreSolar/bms-8s50-ic) @`3904858673` | CC-BY-SA-4.0; LICENSE: “Attribution-ShareAlike 4.0 International” | KiCad legacy (5) | yes; libs cache lib present | 5 since 2019 / 0 / — | yes (8 open) | community builds | 5 commits |

### mixed_signal_instrumentation (6)

| Rank | Board | Source | Licence (quoted) | CAD format | Sch + PCB + libs | History (commits since / tags / rev dirs) | Issue tracker / errata | Proven | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **1 · P0** | **Glasgow Interface Explorer** | [GlasgowEmbedded/glasgow](https://github.com/GlasgowEmbedded/glasgow) `hardware/boards/glasgow` @`cc9a083c82` | 0BSD OR Apache-2.0; LICENSE-0BSD.txt: “Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted.” (+ LICENSE-Apache-2.0.txt) | KiCad legacy (early revisions) → KiCad 6+ (later) | yes; libs cache libs present | 3128 since 2018 / 0 / 8 | yes (81 open) | crowdfunded product | revA→revD1 directories + revision notes |
| 2 | EEZ BB3 modular bench box | [eez-open/modular-psu](https://github.com/eez-open/modular-psu) @`c3d30a6210` | TAPR-OHL-1.0; LICENSE.TXT: “The TAPR Open Hardware License Version 1.0 (May 25, 2007)” | Eagle XML + KiCad 6+ | yes; libs embedded / project libs | 434 since 2018 / 4 / — | yes (17 open) | crowdfunded product | many modules (PSU, SMU, DMM) |
| 3 | ThunderScope (4-ch 1 GS/s PCIe/TB scope) | [EEVengers/ThunderScope](https://github.com/EEVengers/ThunderScope) @`df2b845927` | MIT; LICENSE: “MIT License Copyright (c) 2021 EEVengers” | KiCad 6+ | yes; libs project libs present | 1320 since 2019 / 0 / — | yes (7 open) | crowdfunded product | also high-speed digital (PCIe) |
| 4 | FOSSASIA PSLab | [fossasia/pslab-hardware](https://github.com/fossasia/pslab-hardware) @`404aae2169` | Apache-2.0; LICENSE: “Apache License Version 2.0, January 2004” | KiCad legacy + 6+ | yes; libs cache libs present | 303 since 2017 / 3 / — | yes (20 open) | vendor product (v5/v6) |  |
| 5 | Open Ephys acquisition board | [open-ephys/acquisition-board](https://github.com/open-ephys/acquisition-board) @`8bf3c0c4da` | CC-BY-NC-SA-3.0-IGO; README: “Creative Commons Attribution-NonCommercial-ShareAlike 3.0 IGO License” | Eagle XML + KiCad legacy | yes; libs embedded (Eagle) | 298 since 2012 / 4 / — | yes (13 open) | vendor product (lab use) | NC licence: fine for personal test use (D5); flagged |
| 6 | Bus Pirate 5 | [DangerousPrototypes/BusPirate5-hardware](https://github.com/DangerousPrototypes/BusPirate5-hardware) @`acb9cd667c` | MIT; LICENSE.TXT: “MIT License Copyright (c) 2023 Ian Lesnet, Where Labs LLC” | KiCad 6+ | yes; libs project libs present | 103 since 2023 / 0 / — | yes (0 open) | vendor product |  |

### mcu_iot (6)

| Rank | Board | Source | Licence (quoted) | CAD format | Sch + PCB + libs | History (commits since / tags / rev dirs) | Issue tracker / errata | Proven | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **1 · P0** | **Olimex ESP32-EVB (Wi-Fi/BLE, Ethernet, CAN, relays)** | [OLIMEX/ESP32-EVB](https://github.com/OLIMEX/ESP32-EVB) `HARDWARE` @`a3ec2f4481` | Apache-2.0; LICENSE: “Apache License Version 2.0, January 2004” | KiCad legacy → KiCad 6+ | yes; libs cache libs present | 76 since 2017 / 0 / 12 | yes (0 open) | vendor product | Rev A→L directories + per-revision change list |
| 2 | Watchy (e-paper smartwatch) | [sqfmi/watchy-hardware](https://github.com/sqfmi/watchy-hardware) @`5f147972aa` | MIT; LICENSE: “MIT License Copyright (c) 2026 SQFMI” | KiCad 6+ | yes; libs project libs present | 27 since 2020 / 3 / — | yes (6 open) | vendor product | 3 tags |
| 3 | Adafruit Feather ESP32-S3 | [adafruit/Adafruit-Feather-ESP32-S3-PCB](https://github.com/adafruit/Adafruit-Feather-ESP32-S3-PCB) @`2c95ed51a5` | CC-BY-SA-3.0; license.txt: “Attribution-ShareAlike 3.0 Unported” | Eagle XML | yes; libs embedded (Eagle) | 9 since 2022 / 0 / — | yes (1 open) | vendor product |  |
| 4 | Adafruit Feather RP2040 | [adafruit/Adafruit-Feather-RP2040-PCB](https://github.com/adafruit/Adafruit-Feather-RP2040-PCB) @`ea88166891` | CC-BY-SA (version not stated); README: “Creative Commons Attribution, Share-Alike license, check license.txt” (license.txt absent at HEAD) | Eagle XML | yes; libs embedded (Eagle) | 11 since 2021 / 0 / — | yes (0 open) | vendor product | licence version unstated: flagged |
| 5 | Adafruit nRF52 Bluefruit Feather | [adafruit/Adafruit-nRF52-Bluefruit-Feather-PCB](https://github.com/adafruit/Adafruit-nRF52-Bluefruit-Feather-PCB) @`20cfdfdcf6` | CC-BY-SA (version not stated); README: “Creative Commons Attribution, Share-Alike license, check license.txt” (license.txt absent at HEAD) | Eagle XML | yes; libs embedded (Eagle) | 10 since 2017 / 0 / — | yes (5 open) | vendor product | licence version unstated: flagged |
| 6 | Adafruit HUZZAH32 ESP32 Feather | [adafruit/Adafruit-HUZZAH32-ESP32-Feather-PCB](https://github.com/adafruit/Adafruit-HUZZAH32-ESP32-Feather-PCB) @`e23514c6f9` | CC-BY-SA-3.0; license.txt: “Attribution-ShareAlike 3.0 Unported” | Eagle XML | yes; libs embedded (Eagle) | 7 since 2017 / 0 / — | yes (0 open) | vendor product |  |

### high_speed_digital (9)

| Rank | Board | Source | Licence (quoted) | CAD format | Sch + PCB + libs | History (commits since / tags / rev dirs) | Issue tracker / errata | Proven | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **1 · P0** | **Antmicro Jetson Nano/TX2 NX baseboard** | [antmicro/jetson-nano-baseboard](https://github.com/antmicro/jetson-nano-baseboard) @`d8d0b2d71a` | Apache-2.0; README: “This project is licensed under the Apache-2.0 license.” | KiCad 6+ | yes; libs embedded | 238 since 2019 / 5 / — | yes (6 open) | vendor open reference design | 5 tags |
| 2 | Olimex A20-OLinuXino-LIME2 (DDR3, GbE) | [OLIMEX/OLINUXINO](https://github.com/OLIMEX/OLINUXINO) `HARDWARE/A20-OLinuXino-LIME2` @`85eee9c17e` | CC-BY-SA-3.0-US; README: “licensed under permissive Creative Commons Attribution-Share Alike 3.0 United States License (CC-BY-SA 3.0 US)” | Eagle → KiCad (varies by revision) | yes; libs embedded / cache libs | 510 since 2012 / 0 / 114 | yes (2 open) | vendor product | revisions B→M (14 directories) + latest |
| 3 | ULX3S (ECP5, SDRAM, USB, HDMI) | [emard/ulx3s](https://github.com/emard/ulx3s) @`6a92cec6b1` | MIT (modified); LICENSE.md: “Permission is hereby granted, free of charge, to any person obtaining a copy of this hardware, software, and associated documentation files (the "Product")…” | KiCad legacy (5) | yes; libs cache lib present | 2732 since 2016 / 16 / — | yes (5 open) | crowdfunded product | 16 tags, 2.7k commits |
| 4 | ButterStick (ECP5, DDR3, GigE, SYZYGY) | [butterstick-fpga/butterstick-hardware](https://github.com/butterstick-fpga/butterstick-hardware) @`0694ddfa98` | CERN-OHL-1.2; README: “Hardware in this repository is licensed under CERN OHL v1.2” | KiCad legacy → 6+ | yes; libs cache libs present | 96 since 2019 / 1 / — | yes (3 open) | crowdfunded product |  |
| 5 | OrangeCrab (ECP5, DDR3) | [orangecrab-fpga/orangecrab-hardware](https://github.com/orangecrab-fpga/orangecrab-hardware) @`c511d569fe` | CERN-OHL-1.2; README: “Hardware in this repository is licensed under CERN OHL v1.2” | KiCad legacy (5) | yes; libs cache libs present | 147 since 2019 / 0 / — | yes (17 open) | crowdfunded product |  |
| 6 | Antmicro Kria K26 SOM devboard | [antmicro/kria-k26-devboard](https://github.com/antmicro/kria-k26-devboard) @`9ead4dfc5b` | Apache-2.0; README: “This project is published under the Apache-2.0 license.” | KiCad 6+ | yes; libs project libs present | 203 since 2022 / 0 / — | yes (1 open) | vendor open reference design |  |
| 7 | Hat Labs HALPI2 (CM5 carrier) | [hatlabs/HALPI2-hardware](https://github.com/hatlabs/HALPI2-hardware) @`43a04926cf` | CERN-OHL-S-2.0; LICENSE.md: “This source describes Open Hardware and is licensed under the CERN-OHL-S v2.” | KiCad 6+ | yes; libs project libs present | 149 since 2025 / 17 / — | yes (1 open) | vendor product | 17 tags |
| 8 | Antmicro LPDDR4 test board | [antmicro/lpddr4-test-board](https://github.com/antmicro/lpddr4-test-board) @`5c8b5d3c27` | Apache-2.0; README: “This project is published under the Apache-2.0 license.” | KiCad 6+ | yes; libs embedded | 264 since 2020 / 0 / — | yes (0 open) | vendor open reference design |  |
| 9 | Precursor / Betrusted | [betrusted-io/betrusted-hardware](https://github.com/betrusted-io/betrusted-hardware) @`fecea1061a` | CERN-OHL-1.2; README: “Licensed under the CERN OHL v1.2” | Altium (.PcbDoc/.SchDoc) | yes; libs embedded (Altium) | 31 since 2019 / 0 / — | yes (3 open) | crowdfunded product | Altium reader test case |

### rf (8)

| Rank | Board | Source | Licence (quoted) | CAD format | Sch + PCB + libs | History (commits since / tags / rev dirs) | Issue tracker / errata | Proven | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **1 · P0** | **HackRF One (1 MHz–6 GHz SDR)** | [greatscottgadgets/hackrf](https://github.com/greatscottgadgets/hackrf) `hardware/hackrf-one` @`468c5d93b6` | GPL-2.0; COPYING: “GNU GENERAL PUBLIC LICENSE Version 2, June 1991” | KiCad 6+ (older boards KiCad legacy) | yes; libs cache libs present | 3482 since 2012 / 18 / — | yes (86 open) | vendor product | r1→r10 documented in docs/source/list_of_hardware_revisions.rst |
| 2 | LibreVNA (100 kHz–6 GHz VNA) | [jankae/LibreVNA](https://github.com/jankae/LibreVNA) @`c477819102` | GPL-3.0; LICENSE: “GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007” | KiCad 6+ (some Eagle files) | yes; libs project libs present | 1172 since 2020 / 21 / — | yes (5 open) | community / third-party builds | 21 tags |
| 3 | Ubertooth One (2.4 GHz) | [greatscottgadgets/ubertooth](https://github.com/greatscottgadgets/ubertooth) @`c9dfdbd76a` | GPL-2.0; COPYING: “GNU GENERAL PUBLIC LICENSE Version 2, June 1991” | KiCad legacy | yes; libs cache libs present | 1355 since 2010 / 23 / — | yes (53 open) | vendor product | 23 tags |
| 4 | LimeSDR Mini v2 | [myriadrf/LimeSDR-Mini-v2](https://github.com/myriadrf/LimeSDR-Mini-v2) @`c90822b5d3` | SHL-2.1; README: “The hardware designs are licensed under the Solderpad Hardware License v2.1.” | Altium | yes; libs embedded + SchLib/PcbLib | 92 since 2022 / 0 / — | yes (6 open) | crowdfunded product | Altium reader test case |
| 5 | YARD Stick One (sub-1 GHz) | [greatscottgadgets/yardstick](https://github.com/greatscottgadgets/yardstick) @`2e6516c330` | GPL-2.0; COPYING: “GNU GENERAL PUBLIC LICENSE Version 2, June 1991” | KiCad 6+ | yes; libs cache lib present | 179 since 2012 / 0 / — | yes (4 open) | vendor product |  |
| 6 | CatSniffer (multi-band IoT sniffer) | [ElectronicCats/CatSniffer](https://github.com/ElectronicCats/CatSniffer) @`f7af615e3b` | CERN-OHL-1.2; LICENSE.md: “CERN Open Hardware Licence v1.2” | KiCad 6+ | yes; libs project libs present | 465 since 2020 / 7 / — | yes (2 open) | vendor product | 7 tags |
| 7 | PortaPack H1 (HackRF add-on) | [sharebrained/portapack-hackrf](https://github.com/sharebrained/portapack-hackrf) @`df26d8a99a` | GPL-2.0-or-later; README: “GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version” | KiCad legacy | yes; libs no cache lib (needs KiCad 4/5 libs) | 1618 since 2015 / 7 / — | yes (69 open); repo archived | vendor product (discontinued) | archived repo |
| 8 | LimeSDR Mini (v1) | [myriadrf/LimeSDR-Mini](https://github.com/myriadrf/LimeSDR-Mini) @`ade9b06639` | CC-BY-3.0; README: “The hardware designs are licensed under a Creative Commons Attribution 3.0 Unported licence.” | Altium | yes; libs embedded + SchLib/PcbLib | 7 since 2017 / 0 / — | yes (2 open) | crowdfunded product | 7 commits |

### motor_drive (6)

| Rank | Board | Source | Licence (quoted) | CAD format | Sch + PCB + libs | History (commits since / tags / rev dirs) | Issue tracker / errata | Proven | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **1 · P0** | **mjbots moteus r4.x / c1 (FOC BLDC)** | [mjbots/moteus](https://github.com/mjbots/moteus) `hw/c1, hw/controller` @`aec0333888` | Apache-2.0; README: “available under an Apache 2.0 License” | Eagle XML (r4.x) + KiCad 6+ (c1) | yes; libs embedded / project libs | 2296 since 2018 / 77 / 12 | yes (2 open) | vendor product | r4.3→r4.11 and c1 r1.0→r1.2 directories |
| **2 · P0** | **ODrive v3.x (dual BLDC)** | [odriverobotics/ODriveHardware](https://github.com/odriverobotics/ODriveHardware) `v3` @`079f30ac58` | MIT; LICENSE: “The MIT License (MIT) Copyright (c) 2016 Oskar Weigl (madcowswe)” | Altium | yes; libs embedded + SchLib/PcbLib | 60 since 2016 / 5 / 2 | yes (9 open) | vendor product | v3/CHANGELOG.md (up to v3.5); Altium reader test case |
| 3 | Watterott SilentStepStick (TMC2xxx) | [watterott/SilentStepStick](https://github.com/watterott/SilentStepStick) @`098b7518b0` | CC-BY-SA-4.0; README: “released under the Creative Commons Attribution Share-Alike License” (link to by-sa/4.0) | Eagle XML | yes; libs embedded (Eagle) | 106 since 2017 / 0 / — | disabled | vendor product | issues disabled; 24 Eagle board files (variants) |
| 4 | SimpleFOC Shield | [simplefoc/Arduino-SimpleFOCShield](https://github.com/simplefoc/Arduino-SimpleFOCShield) @`2a83626b86` | MIT; LICENSE: “MIT License Copyright (c) 2020 SimpleFOCproject” | Altium (+ EasyEDA) | yes; libs embedded (Altium) | 42 since 2020 / 10 / — | yes (6 open) | vendor product | 10 tags |
| 5 | Ben Katz 3-phase integrated controller | [bgkatz/3phase_integrated](https://github.com/bgkatz/3phase_integrated) @`c15b16dc92` | MIT; LICENSE: “MIT License Copyright (c) 2020 Ben Katz” | Eagle XML | yes; libs embedded (Eagle) | 14 since 2018 / 0 / — | yes (12 open) | community builds (Mini Cheetah lineage) |  |
| 6 | SparkFun Qwiic Motor Driver | [sparkfun/Qwiic_Motor_Driver](https://github.com/sparkfun/Qwiic_Motor_Driver) @`4ee98aac48` | CC-BY-SA-4.0; License.md: “SparkFun hardware is released under Creative Commons Share-alike 4.0 International” | Eagle XML | yes; libs embedded (Eagle) | 14 since 2019 / 0 / — | yes (0 open) | vendor product |  |

### industrial_comms (6)

| Rank | Board | Source | Licence (quoted) | CAD format | Sch + PCB + libs | History (commits since / tags / rev dirs) | Issue tracker / errata | Proven | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **1 · P0** | **Hat Labs SH-ESP32 (CAN / NMEA 2000 marine interface)** | [hatlabs/SH-ESP32-hardware](https://github.com/hatlabs/SH-ESP32-hardware) @`167981ac8c` | CC-BY-4.0; LICENSE.md: “Creative Commons Attribution 4.0 International” | KiCad legacy → KiCad 6+ | yes; libs project libs present | 241 since 2020 / 10 / — | yes (1 open) | vendor product | 10 tags |
| 2 | Olimex ESP32-POE-ISO (isolated PoE Ethernet) | [OLIMEX/ESP32-POE-ISO](https://github.com/OLIMEX/ESP32-POE-ISO) `HARDWARE` @`c19307ecd8` | CERN-OHL-1.2; README: “It is distributed with CERN OSH Licensee v1.2” | KiCad legacy → KiCad 6+ | yes; libs cache libs present | 24 since 2019 / 0 / 14 | yes (0 open) | vendor product | Rev B→N1 revision directories |
| 3 | candleLight (USB-CAN) | [HubertD/candleLight](https://github.com/HubertD/candleLight) @`ada932f5b3` | CERN-OHL-1.2; LICENSE.txt: “CERN Open Hardware Licence v1.2” | KiCad legacy | yes; libs cache lib present | 25 since 2016 / 1 / — | yes (1 open) | community builds, widely cloned |  |
| 4 | Adafruit RP2040 CAN Bus Feather | [adafruit/Adafruit-RP2040-CAN-Bus-Feather-PCB](https://github.com/adafruit/Adafruit-RP2040-CAN-Bus-Feather-PCB) @`499de161e9` | CC-BY-SA-3.0; license.txt: “Attribution-ShareAlike 3.0 Unported” | Eagle XML | yes; libs embedded (Eagle) | 9 since 2023 / 0 / — | yes (1 open) | vendor product |  |
| 5 | Olimex A20-CAN-ISO (isolated CAN) | [OLIMEX/A20-CAN-ISO](https://github.com/OLIMEX/A20-CAN-ISO) @`ca00943e63` | CERN-OHL-1.2; README: “Hardware is with CERN OHL v1.2 license” | KiCad legacy | yes; libs cache lib present | 3 since 2020 / 0 / — | yes (0 open) | vendor product | 3 commits |
| 6 | RS-485/TTL USB dongle (CH340E) | [TG9541/rs485-ch340e](https://github.com/TG9541/rs485-ch340e) @`9b5cb4c93c` | MIT; LICENSE: “MIT License Copyright (c) 2019 Thomas” | KiCad legacy | yes; libs no cache lib (needs KiCad 5 libs) | 5 since 2019 / 0 / — | yes (0 open) | hobby build only | proven: weak; RS-485 coverage gap filler |

### intrinsically_safe_atex (0)

No board met the selection criteria; see § 4.

### automotive (5)

| Rank | Board | Source | Licence (quoted) | CAD format | Sch + PCB + libs | History (commits since / tags / rev dirs) | Issue tracker / errata | Proven | Notes |
|---|---|---|---|---|---|---|---|---|---|
| **1 · P0** | **Speeduino v0.4 ECU** | [speeduino/Hardware](https://github.com/speeduino/Hardware) `v0.4` @`ee6a924283` | CERN-OHL-S-2.0; README: “These designs are release under the CERN Open Hardware License v2 Strong (CERN-OHL-S).” | KiCad legacy + 6+ | yes; libs cache libs present | 28 since 2021 / 0 / 8 | yes (19 open) | community product (>1000 installs per README) | 0.4.4b→c→d and 0.4.3→0.4.3d directories |
| 2 | openinverter (Huebner) inverter main board + gate drivers | [jsphuebner/inverter-hardware](https://github.com/jsphuebner/inverter-hardware) @`355c9b6b0e` | CC-BY-SA-4.0; README: “# License Attribution-ShareAlike 4.0 International” | Eagle XML + KiCad | yes; libs embedded / project libs | 37 since 2018 / 0 / — | yes (0 open) | community EV-conversion builds | also motor drive |
| 3 | foxBMS 1 master/slave (Fraunhofer IISB) | [foxBMS/foxBMS-hardware](https://github.com/foxBMS/foxBMS-hardware) @`190b2dd54b` | CC-BY-4.0; LICENSE.md: “This work is licensed under the Creative Commons Attribution 4.0 International License.” | Eagle XML | yes; libs embedded (Eagle) | 17 since 2017 / 2 / — | yes (0 open); repo archived | research reference platform | archived; 2 tags; also power |
| 4 | CANtact (USB-CAN car tool) | [linklayer/cantact-hw](https://github.com/linklayer/cantact-hw) @`470252258a` | CC-BY-SA-4.0; LICENSE.md: “Hardware files released under the Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0)” | KiCad legacy | yes; libs no cache lib (needs KiCad 4 libs) | 9 since 2015 / 3 / — | yes (6 open) | vendor product | 3 tags |
| 5 | rusEFI classic modules (1A injector 6/12-ch, power supply, analog input) | [rusefi/rusefi-hardware](https://github.com/rusefi/rusefi-hardware) `classic-designs/*` @`ea01307410` | MIT (per module); classic-designs/1A_injector_6-channels/MIT_License.txt: “Permission is hereby granted, free of charge…”; no repo-level licence | KiCad legacy | yes; libs project libs present | 1118 since 2021 / 5 / 2 | yes (114 open) | community builds | only the MIT-licensed module folders qualify; in-repo rev00.10→rev00.20 history |

