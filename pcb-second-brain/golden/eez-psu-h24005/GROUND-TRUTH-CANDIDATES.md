# Ground-truth candidates — eez-psu-h24005

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/eez-open/psu-hw
- Pinned commit: `eef089135b58d398419a0d1138970fe71320d6a7` (2021-11-05T23:07:47+01:00)
- Paths analysed: `.`
- Domain: power_electronics
- Commits touching these paths: 241; touching design files: 126; design-file commits whose message matches a fix keyword: 54

## Revision pairs

No revision directories; revisions exist only as git history (pairs = consecutive commits below).

Tags in repository (4): `1.0`, `2.0`, `3.0`, `r5B13a`

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [19](https://github.com/eez-open/psu-hw/issues/19) | closed | 2017-01-27 | Additional Vin capacitor on Power board |
| [31](https://github.com/eez-open/psu-hw/issues/31) | closed | 2017-03-11 | 1M resistor in parallel with Ethernet controller crystal |
| [33](https://github.com/eez-open/psu-hw/issues/33) | closed | 2017-03-20 | Cannot reach max. current in LOW range mode |
| [37](https://github.com/eez-open/psu-hw/issues/37) | open | 2017-04-01 | Migrate Eagle to free EDA such as KiCad |
| [40](https://github.com/eez-open/psu-hw/issues/40) | closed | 2017-06-17 | High range sense resistor |
| [41](https://github.com/eez-open/psu-hw/issues/41) | closed | 2017-06-17 | Additional current range(s) |
| [42](https://github.com/eez-open/psu-hw/issues/42) | closed | 2017-07-25 | Output enters unregulated mode too early (below ~1.8 V)  |
| [43](https://github.com/eez-open/psu-hw/issues/43) | closed | 2017-07-25 | X12, X14 connection improvement on PCB  |
| [44](https://github.com/eez-open/psu-hw/issues/44) | closed | 2017-07-25 | Remote sense resistors |
| [46](https://github.com/eez-open/psu-hw/issues/46) | closed | 2017-07-27 | Overvoltage protection for V+ output |
| [47](https://github.com/eez-open/psu-hw/issues/47) | closed | 2017-08-13 | AC power switch failure and TR1 +12V rail failure |
| [52](https://github.com/eez-open/psu-hw/issues/52) | closed | 2017-08-26 | Difficult or unable to attach power modules to case... |
| [53](https://github.com/eez-open/psu-hw/issues/53) | closed | 2017-08-26 | Arduino Due issue with USB "native" port |
| [56](https://github.com/eez-open/psu-hw/issues/56) | closed | 2017-09-01 | USB isolator speed selection |
| [57](https://github.com/eez-open/psu-hw/issues/57) | closed | 2017-09-14 | X12 X14 Pin confusion |
| [60](https://github.com/eez-open/psu-hw/issues/60) | closed | 2017-11-04 | R93 PCB layout spacing |
| [61](https://github.com/eez-open/psu-hw/issues/61) | closed | 2017-11-15 | Consider presenting at FOSDEM in 2018 |
| [70](https://github.com/eez-open/psu-hw/issues/70) | closed | 2020-01-18 | Expected self test failures |
| [71](https://github.com/eez-open/psu-hw/issues/71) | open | 2026-03-01 | Bare pcb board (03/2017), r5B12 |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`ba9722a91`](https://github.com/eez-open/psu-hw/commit/ba9722a91d3e521bfd66fb86f9155893684e3980) | 2018-08-29 | generic fix | Fixed schematic Layout | AC_soft-start,_+5V_power.sch, Bias_power.sch, Binding_posts_relays,_LED_driver.sch, CC%2fC |
| [`4e7e1c1a4`](https://github.com/eez-open/psu-hw/commit/4e7e1c1a4cda6621e515751d5dd8e978b38d4789) | 2018-08-29 | routing, generic fix | Fixed Ground Planes & Vias | EEZ PSU consolidated r5B13a.kicad_pcb |
| [`5be86e0fc`](https://github.com/eez-open/psu-hw/commit/5be86e0fc6577dc627628d8978013c382751ce21) | 2017-07-09 | footprint, generic fix | IC22 package fix | EEZ PSU consolidated r5B12.brd, EEZ PSU consolidated r5B12.sch |
| [`927faed3d`](https://github.com/eez-open/psu-hw/commit/927faed3dd3ba0b2afc0c17278cc0e03d7c2d79c) | 2017-04-17 | footprint, generic fix | SKL16 diode package fix | EEZ PSU consolidated r5B12.sch |
| [`ad4d7b106`](https://github.com/eez-open/psu-hw/commit/ad4d7b106f6824d6ebbb1ec822d5750545ee72fe) | 2017-04-13 | silkscreen | top silkscreen clean up | EEZ PSU consolidated r5B12.brd, EEZ PSU r5B12 panel.brd |
| [`ead7c2145`](https://github.com/eez-open/psu-hw/commit/ead7c2145f7aeca1e3371a04cb19f7d57b9596f7) | 2017-03-15 | generic fix | wire harness corrections | EEZ PSU consolidated r5B12.sch |
| [`a734a2b20`](https://github.com/eez-open/psu-hw/commit/a734a2b20f5fd6d4a0f2c00e9604f9869de59048) | 2017-03-15 | footprint | IC15 package, wire harness update | EEZ PSU consolidated r5B12.brd, EEZ PSU consolidated r5B12.sch, EEZ PSU r5B12 panel.brd |
| [`741c6a0b2`](https://github.com/eez-open/psu-hw/commit/741c6a0b2ba0c878f50558fb919a4114474fada3) | 2017-03-13 | part swap | X14 pin 1, 2 swap, R177 | EEZ PSU consolidated r5B12.brd, EEZ PSU consolidated r5B12.sch, EEZ PSU r5B12 panel.brd |
| [`34213f167`](https://github.com/eez-open/psu-hw/commit/34213f1675ad50fe9ddbdd3816e7f6977ebe5383) | 2017-03-11 | part swap, power | Q24-> mosfet, swap outputs of Q32, Q33 | EEZ PSU consolidated r5B12.brd, EEZ PSU consolidated r5B12.sch, EEZ PSU r5B12 panel.brd |
| [`5dd149043`](https://github.com/eez-open/psu-hw/commit/5dd1490436dd035bc8d9ec9bb09381222bc4e75a) | 2017-01-11 | silkscreen, generic fix | Corrections on bottom silkscreen | EEZ PSU consolidated r5B11.brd, EEZ PSU r5B11 panel.brd |
| [`33409e514`](https://github.com/eez-open/psu-hw/commit/33409e51480bfd9215bba32036ed8b7fd0065c84) | 2017-01-03 | part swap | optimized values for CC resistors | EEZ PSU consolidated r5B11.brd, EEZ PSU consolidated r5B11.sch |
| [`3157a65ce`](https://github.com/eez-open/psu-hw/commit/3157a65cea9034a196802cbc56a6675d3ea6a236) | 2016-11-04 | footprint | J6 earth PCB pad | EEZ PSU r5B10 panel.brd |
| [`d1bc38c9b`](https://github.com/eez-open/psu-hw/commit/d1bc38c9b04ea741563084b014eefeab8b976261) | 2016-10-25 | silkscreen, generic fix | r5B10 bottom silkscreen another minor correction | EEZ PSU consolidated r5B10.brd, EEZ PSU r5B10 panel.brd |
| [`ebdbc8fa8`](https://github.com/eez-open/psu-hw/commit/ebdbc8fa8ed63521d84e44c3d3a7dbac4a3f18b6) | 2016-10-25 | part swap, silkscreen, generic fix | r5B10 BOM, top silkscreen minor correction | EEZ PSU consolidated r5B10.brd, EEZ PSU consolidated r5B10.sch, EEZ PSU r5B10 panel.brd |
| [`bb7e5d793`](https://github.com/eez-open/psu-hw/commit/bb7e5d793cb2e483802e241715afc77ed6c30246) | 2016-09-30 | part swap | NTC values updated | EEZ PSU consolidated r5B9.brd, EEZ PSU consolidated r5B9.sch |
| [`4e6cab46b`](https://github.com/eez-open/psu-hw/commit/4e6cab46bf51934201ef824b536ea7c29c974fa9) | 2016-09-29 | silkscreen, generic fix | Changed wrong designators in schematics on page 3, 6 and 8 | EEZ PSU consolidated r5B9.sch |
| [`8a169d0e1`](https://github.com/eez-open/psu-hw/commit/8a169d0e17479bb5d632157819804b3f5cb5c54b) | 2016-09-29 | part swap | new BOM, added RS and DigiKeyne | EEZ PSU AUX PS r5B9.brd, EEZ PSU AUX PS r5B9.sch, EEZ PSU Arduino shield+BP r3B4.brd, EEZ  |
| [`18af1a7bf`](https://github.com/eez-open/psu-hw/commit/18af1a7bfac31f6d22f063b26c7c513cb686c1e2) | 2016-09-26 | decoupling, part swap | few caps values changed | EEZ PSU H24005 enclosure r5.brd |
| [`7152e2891`](https://github.com/eez-open/psu-hw/commit/7152e2891854bcfe131882bf3c205a09c304991d) | 2016-09-17 | power | Removed thermal bridge from enclosure drawings | EEZ PSU H24005 enclosure r5.1.brd |
| [`913163542`](https://github.com/eez-open/psu-hw/commit/913163542f71b5ac2ffbcfe4b44181ad44f6d7e4) | 2016-09-04 | part swap | R97 and R102 values swapped | EEZ PSU consolidated r5B9.brd, EEZ PSU consolidated r5B9.sch, EEZ PSU power board r5B9.brd |
| [`380dc518d`](https://github.com/eez-open/psu-hw/commit/380dc518d0d0bc3ed64398654b1b82f62eaaaca3) | 2016-08-15 | generic fix | binding post holes corrected | EEZ PSU H24005 enclosure r4.6.brd |
| [`4646373a9`](https://github.com/eez-open/psu-hw/commit/4646373a92ee0c704dff1f22dbd105bfb2ea993e) | 2016-07-28 | silkscreen, generic fix | Small corrections on top silkscreen | EEZ PSU AUX PS r5B9.brd, EEZ PSU Arduino shield+BP r3B4.brd, EEZ PSU power board r5B9.brd, |
| [`c246063bf`](https://github.com/eez-open/psu-hw/commit/c246063bf48edd5d16177906f2c1437690804a01) | 2016-07-25 | part swap, power, routing, generic fix | BOM corrections, Q1 gate trace improved | EEZ PSU consolidated r5B9.brd, EEZ PSU power board r5B9.brd, EEZ PSU r5B9 panel.brd |
| [`d4d3a63df`](https://github.com/eez-open/psu-hw/commit/d4d3a63dfd0a6d5d846f0649c2f030b203f02501) | 2016-07-25 | silkscreen, generic fix | small correction on Arduino shield PCB silkscreens | EEZ PSU Arduino shield+BP r3B4.brd, EEZ PSU consolidated r5B9.brd, EEZ PSU r5B9 panel.brd |
| [`2d062569d`](https://github.com/eez-open/psu-hw/commit/2d062569d4df020b9ae0a537359f745d5fbd230c) | 2016-07-24 | part swap | AUX PS: R107 value changed | EEZ PSU AUX PS r5B9.brd, EEZ PSU AUX PS r5B9.sch |
| [`da017610d`](https://github.com/eez-open/psu-hw/commit/da017610da2ea852f9f681443fa3bdd02533b3f3) | 2016-06-30 | part swap | Rear panel for alternative IEC inlet, LTspice models | EEZ PSU H24005 enclosure r3.brd, EEZ PSU H24005 rear panel (4304.6090).brd |
| [`4fec135a4`](https://github.com/eez-open/psu-hw/commit/4fec135a42f4b20b5da3e6129db565142328cfd0) | 2016-06-27 | part swap | Myrra AC/DC module for 115VAC, BOM updated | EEZ AUX PS r5B8.brd, EEZ AUX PS r5B8.sch, EEZ PSU consolidated r5B8.brd, EEZ PSU power boa |
| [`bbb01e5b2`](https://github.com/eez-open/psu-hw/commit/bbb01e5b24709745fd7fe42bdb3f1c1248af690a) | 2016-06-26 | part swap | Consolidated BOM for PCBA, AUX PS BOM update, Power module PCB update | EEZ PSU consolidated r5B8.brd, EEZ PSU power board r5B8.brd, EEZ PSU r5B8 panel.brd |
| [`ee4305b76`](https://github.com/eez-open/psu-hw/commit/ee4305b76338ad054c2aad5829af8abf74ea5152) | 2016-06-25 | silkscreen | Missed hole added to Power board, minor bottom silk mask changes | EEZ AUX PS r5B8.brd, EEZ Arduino shield+BP r3B3.brd, EEZ PSU consolidated r5B8.brd, EEZ PS |
| [`3de0c1e9d`](https://github.com/eez-open/psu-hw/commit/3de0c1e9d985ddfe53f28a055fe7fbef33832d74) | 2016-06-25 | part swap | Power board PCB layout and BOM changes | EEZ PSU power board r5B8.brd, EEZ PSU power board r5B8.sch |
| [`85ed365ee`](https://github.com/eez-open/psu-hw/commit/85ed365ee5d2eb57f67e7aa65c7d189c7b9c4831) | 2016-06-25 | part swap, generic fix | Fixed fan control, miniMELF replaced with SOD-323 | EEZ PSU consolidated r5B8.brd, EEZ PSU consolidated r5B8.sch, EEZ PSU r5B8 panel.brd |
| [`91c8e45fc`](https://github.com/eez-open/psu-hw/commit/91c8e45fc40a496b3c426fa9a7e24ce6d29ccefd) | 2016-06-19 | footprint | Power board: increase D3 bottom layer pad | EEZ PSU power board r5B8.brd, EEZ PSU power board r5B8.sch |
| [`949774cdb`](https://github.com/eez-open/psu-hw/commit/949774cdb2bd307de8e86576bdaee1b7d2ee6863) | 2016-06-19 | part swap, generic fix | Rt value corrections for AUX PS and Power board | EEZ AUX PS r5B8.brd, EEZ AUX PS r5B8.sch, EEZ PSU power board r5B8.brd, EEZ PSU power boar |
| [`357ff0eee`](https://github.com/eez-open/psu-hw/commit/357ff0eee383591c08cc58ede94731c850580488) | 2016-06-19 | part swap, generic fix | Arduino shield D5 PCB layout, BOM corrections | EEZ Arduino shield+BP r3B3.brd |
| [`314a510aa`](https://github.com/eez-open/psu-hw/commit/314a510aaeadf6d973a2bc24f8b0536b293da5cc) | 2016-06-19 | generic fix | Enclosure corrections | EEZ PSU H24005 enclosure.brd |
| [`115cb5202`](https://github.com/eez-open/psu-hw/commit/115cb5202aefac7ed3dd34f9ef411851f8a5fcde) | 2016-06-18 | part swap, generic fix | Power board R28, R29 PCB layout, BOM corrections | EEZ PSU power board r5B8.brd |
| [`702ff816a`](https://github.com/eez-open/psu-hw/commit/702ff816a7e4d2e01af7baf5a00ddc55898d472c) | 2016-06-17 | generic fix | Arduino Shield correction | EEZ Arduino shield+BP r3B3.brd |
| [`b9cd1839d`](https://github.com/eez-open/psu-hw/commit/b9cd1839d933e7418278101a0e00806e5558f4bb) | 2016-06-17 | generic fix | Screen masks correction | EEZ AUX PS r5B8.brd, EEZ Arduino shield+BP r3B3.brd |
| [`c0d458659`](https://github.com/eez-open/psu-hw/commit/c0d4586597c500af1dacba31816adafd4e071b4e) | 2016-06-16 | generic fix | Power board tiny text corrections | EEZ PSU power board r5B8.sch |
| [`853b948ab`](https://github.com/eez-open/psu-hw/commit/853b948ab630f0b11dc70561776ca3411285a6c7) | 2016-05-26 | part swap, generic fix | AUX PS r4B3, various BOM corrections | EEZ AUX PS r4B43.brd, EEZ AUX PS r4B43.sch, PCB panel M1 (2016-05-26).brd |
| [`6c03d7db1`](https://github.com/eez-open/psu-hw/commit/6c03d7db18a69b35e70ac8119c23f399a9ae54bb) | 2016-05-01 | decoupling, power | AUX PS ramp capacitor and fan LDO changed | AUX PS board side support (Economica 320x200x115).brd, EEZ AUX PS r4B42.brd, EEZ AUX PS r4 |
| [`16fa5668e`](https://github.com/eez-open/psu-hw/commit/16fa5668e7cf69efab5ba9afe4e7da3107de6962) | 2016-04-17 | power, generic fix | SMPS pre-regulator small corrections | EEZ PSU post-regulator r4B45.sch, EEZ SMPS pre-regulator r1B1.brd, EEZ SMPS pre-regulator  |
| [`26677c5ee`](https://github.com/eez-open/psu-hw/commit/26677c5ee3a519f13c452098eed972783a6c8d92) | 2016-04-09 | part swap | Rubber feet added in Arduino shield BOM | EEZ PSU enclosure (Economica 320x200x115, Arduino shield r1b14).brd |
| [`48cc7aa1b`](https://github.com/eez-open/psu-hw/commit/48cc7aa1b6fe0e756499b076362e567aca6e513c) | 2016-04-08 | footprint | Pads for DC input added into EEZ PSU power board r5B6c | EEZ PSU power board r5B6c.brd, EEZ PSU power board r5B6c.sch |
| [`3791d2e24`](https://github.com/eez-open/psu-hw/commit/3791d2e2406385753acd5744a0d3387b78bd9e61) | 2016-04-07 | part swap | Replaced IEC inlet in preliminary rev5 | EEZ PSU Galaxy 230x80 panels r2.brd |
| [`ee25eb77a`](https://github.com/eez-open/psu-hw/commit/ee25eb77a59a0ba96bd0cd5d086924a9d29b8d00) | 2016-04-07 | part swap | AUX PS r4B41 BOM update | EEZ Arduino shield+BP r1B14.brd |
| [`0c225cd5d`](https://github.com/eez-open/psu-hw/commit/0c225cd5d3bcdd99add535e4644158923433be24) | 2016-03-28 | power, generic fix | Few PCB layout correction for SMPS pre-regulator r1B1 | EEZ SMPS pre-regulator r1B1.brd, EEZ SMPS pre-regulator r1B1.sch |
| [`d43356e6f`](https://github.com/eez-open/psu-hw/commit/d43356e6ffc79f290adb1bbeeb1b0f226ca8f46a) | 2016-03-27 | part swap | Arduino Shield r1B13 BOM updated | EEZ Arduino shield+BP r1B12c BOM.csv, EEZ Arduino shield+BP r1B13.brd, EEZ Arduino shield+ |
| [`76fb76185`](https://github.com/eez-open/psu-hw/commit/76fb7618539bf4398baf95c4953fc719014bb034) | 2016-03-27 | power | Arduino Shield r1B13, Post-regulator r45B45 | AUX PS board side support (Economica 320x200x115).brd, EEZ AUX PS r4B40.brd, EEZ Arduino s |
| [`4e68e211a`](https://github.com/eez-open/psu-hw/commit/4e68e211a716f0cc3f04ab1801dff7f921dd5d0d) | 2016-02-04 | generic fix | minor textual corrections | EEZ AUX PS r4B40.brd, EEZ AUX PS r4B40.sch, EEZ Mosfet pre-regulator r3B3b (SS).sch, EEZ M |
| [`c8500e07a`](https://github.com/eez-open/psu-hw/commit/c8500e07aa6d4b861b2d5351618c8287c3312874) | 2015-12-09 | footprint, generic fix | Small fix on bottom solder mask | EEZ PSU post-regulator r4B44 (SS).brd, EEZ PSU post-regulator r4B44.brd |
| [`ebbc648ee`](https://github.com/eez-open/psu-hw/commit/ebbc648ee23e9c582a6fc29fca3974b423a08b0a) | 2015-12-07 | generic fix | Arduino Shield without BP option folder name fixed | EEZ Arduino shield r1B12b.brd, EEZ Arduino shield r1B12b.sch |
| [`b593e3ca6`](https://github.com/eez-open/psu-hw/commit/b593e3ca646081ef3d9ebf92b88c0ef482eac09e) | 2015-11-10 | power | mosfet pre-regulator r3B3b | EEZ Mosfet pre-regulator r3B3b (SS).brd, EEZ Mosfet pre-regulator r3B3b (SS).sch, EEZ Mosf |
| [`ea96f2b74`](https://github.com/eez-open/psu-hw/commit/ea96f2b749d6fd8a336a2e4a08c80367e75ed963) | 2015-11-10 | power | mosfet post-regulator r4B44 | EEZ PSU post-regulator r4B44 (SS).brd, EEZ PSU post-regulator r4B44 (SS).sch, EEZ PSU post |
