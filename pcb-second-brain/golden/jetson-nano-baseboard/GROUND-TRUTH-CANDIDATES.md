# Ground-truth candidates — jetson-nano-baseboard

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/antmicro/jetson-nano-baseboard
- Pinned commit: `d8d0b2d71ab4cc4bcdb204a09f78be40d4c8cda3` (2025-12-31T13:57:23+01:00)
- Paths analysed: `.`
- Domain: high_speed_digital
- Commits touching these paths: 231; touching design files: 172; design-file commits whose message matches a fix keyword: 71

## Revision pairs

No revision directories; revisions exist only as git history (pairs = consecutive commits below).

Tags in repository (5): `v1.1`, `v1.2`, `v1.3`, `v1.4.2`, `v1.4.6`

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [3](https://github.com/antmicro/jetson-nano-baseboard/issues/3) | closed | 2019-04-03 | PCB files |
| [4](https://github.com/antmicro/jetson-nano-baseboard/issues/4) | closed | 2019-10-15 | Capacitors on the BOM  |
| [5](https://github.com/antmicro/jetson-nano-baseboard/issues/5) | closed | 2020-03-05 | The new PCB doesn't match with the schematic |
| [7](https://github.com/antmicro/jetson-nano-baseboard/issues/7) | closed | 2020-03-09 | HDMI pull down resistors |
| [8](https://github.com/antmicro/jetson-nano-baseboard/issues/8) | closed | 2020-03-18 | the width of the VIN to TPS54561 |
| [10](https://github.com/antmicro/jetson-nano-baseboard/issues/10) | closed | 2020-03-31 | Power sequensing |
| [13](https://github.com/antmicro/jetson-nano-baseboard/issues/13) | closed | 2020-08-13 | The thickness of PCB |
| [14](https://github.com/antmicro/jetson-nano-baseboard/issues/14) | closed | 2020-09-01 | buck capacitors for the Jetson NX SOM |
| [15](https://github.com/antmicro/jetson-nano-baseboard/issues/15) | closed | 2020-09-22 | Any update on Kicad files? |
| [17](https://github.com/antmicro/jetson-nano-baseboard/issues/17) | closed | 2021-01-29 | Routing layers copper pours |
| [19](https://github.com/antmicro/jetson-nano-baseboard/issues/19) | closed | 2021-04-26 | Use of below circuitry |
| [21](https://github.com/antmicro/jetson-nano-baseboard/issues/21) | closed | 2021-07-16 | Undetected CSI camera |
| [38](https://github.com/antmicro/jetson-nano-baseboard/issues/38) | closed | 2025-12-13 | Question: Recommended Manufacturing & PCBA Services for Modified Jetson Baseboard |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`ff7f8198c`](https://github.com/antmicro/jetson-nano-baseboard/commit/ff7f8198cf409fd677fb1bc13e9119190f782bd2) | 2025-09-18 | footprint, generic fix | [#83358] Repair broken footprints after symbols update | jetson-nano-baseboard.kicad_pcb |
| [`18506fe7b`](https://github.com/antmicro/jetson-nano-baseboard/commit/18506fe7b30bef35d32d83020fcacb1b7f11c12d) | 2025-09-18 | generic fix | [#83358] Repair broken symbol after update | display.kicad_sch, id-serial.kicad_sch, interfaces.kicad_sch, som.kicad_sch, supply.kicad_ |
| [`82b2443d6`](https://github.com/antmicro/jetson-nano-baseboard/commit/82b2443d61432b9cf40f896b705359bd875ae9c9) | 2024-04-12 | generic fix | [#57913] fix 2309413-1 | BarrelJack_3167288_Horizontal.kicad_mod, BarrelJack_Pin2mm_MJ-179PH.kicad_mod, Bus_M.2_Soc |
| [`4878ac168`](https://github.com/antmicro/jetson-nano-baseboard/commit/4878ac168f800a6f5a824c893c5a70a7a588a5f2) | 2024-04-11 | generic fix | [#56834] Fix components links | BarrelJack_Pin2mm_MJ-179PH.kicad_mod, C_Pol_EIA-B.kicad_mod, Conn_JST_SH_1x4_SM04B-SRSS-TB |
| [`56339e273`](https://github.com/antmicro/jetson-nano-baseboard/commit/56339e273c61afba2808f72b2d18834f7d096b56) | 2023-05-08 | routing | Route USB-C_ID signal around 3.3V plane | jetson-nano-baseboard.kicad_pcb |
| [`1a7dcee33`](https://github.com/antmicro/jetson-nano-baseboard/commit/1a7dcee33588e3bbe83fc60b44121e2639d46593) | 2023-05-05 | part swap | Replace D5 and D10 with a copy of D12 | display.kicad_sch, jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_sch, usb-d |
| [`e00133872`](https://github.com/antmicro/jetson-nano-baseboard/commit/e00133872bdecc02f1d164dedf37c5fe38d9f6e2) | 2023-05-05 | pull-up/down | Add #RESET pull-down | jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_sch, som.kicad_sch |
| [`0bec1cc8c`](https://github.com/antmicro/jetson-nano-baseboard/commit/0bec1cc8cb88b2fd11b49330bfdff3fb1cd778c6) | 2023-05-05 | generic fix | Fix USB-OTG backfeed issue #44022 | csi.kicad_sch, display.kicad_sch, id-serial.kicad_sch, interfaces.kicad_sch, jetson-nano-b |
| [`407091e2f`](https://github.com/antmicro/jetson-nano-baseboard/commit/407091e2fac5c9885477110b2e59c0798b6a601f) | 2023-05-05 | silkscreen | Align designators | jetson-nano-baseboard.kicad_pcb |
| [`0010bff78`](https://github.com/antmicro/jetson-nano-baseboard/commit/0010bff78d061473f37ccb4cb41755233f74dc10) | 2023-05-05 | footprint | Update local footprints | 0603-cap.kicad_mod, 0603.kicad_mod, DC_SOCKET.kicad_mod, HDMI_Micro-D_Molex_46765-1x01.kic |
| [`385f31916`](https://github.com/antmicro/jetson-nano-baseboard/commit/385f31916804eaabb84158c54c309d526d5769e0) | 2023-05-05 | footprint | Update SC70-3 footprint in the loclib | SC70-3.kicad_mod |
| [`f7523f9c8`](https://github.com/antmicro/jetson-nano-baseboard/commit/f7523f9c81ff3ae4a50e1c425310fcc7907218da) | 2023-05-05 | footprint, generic fix | Fixed D1 and D9 footprints | jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_pro, jetson-nano-baseboard.ki |
| [`122271be9`](https://github.com/antmicro/jetson-nano-baseboard/commit/122271be9b837ba9d89b3bcc120a6a4201bd0967) | 2023-05-05 | part swap | Update parts enumeration | jetson-nano-baseboard.kicad_sch, usb-debug.kicad_sch |
| [`1beef7471`](https://github.com/antmicro/jetson-nano-baseboard/commit/1beef7471a45f76ccfb4568933fd37415453f392) | 2023-05-05 | part swap | Replace USB power switches with version without output discharge | csi.kicad_sch, display.kicad_sch, id-serial.kicad_sch, interfaces.kicad_sch, jetson-nano-b |
| [`d535e8e94`](https://github.com/antmicro/jetson-nano-baseboard/commit/d535e8e94cacabef61e11311777fee11aa487ffa) | 2023-05-05 | footprint, routing | Adjust via and pads to have min 0.125mm annular ring | 2129320-3.kicad_mod, jetson-nano-baseboard.kicad_pcb |
| [`410a9db6c`](https://github.com/antmicro/jetson-nano-baseboard/commit/410a9db6c5a4d5bdf043d9037afaf04951b582a5) | 2023-05-05 | silkscreen | Add labels to the LEDs and buttons | jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_pro |
| [`ede5b1f4e`](https://github.com/antmicro/jetson-nano-baseboard/commit/ede5b1f4edb088fa2dff6e02f7f4716845b6f7d9) | 2023-05-05 | clearance | Clear DRC warnings | IHLP3232DZ.kicad_mod, SI-51005-F.kicad_mod, SM04B-SRSS-TB-LF-SN.kicad_mod, SOD-323.kicad_m |
| [`f3646efc8`](https://github.com/antmicro/jetson-nano-baseboard/commit/f3646efc8f2885494e7c7c5799a02b438c17fb2e) | 2023-05-05 | routing | Auto cleanup tracks and vias | jetson-nano-baseboard.kicad_pcb |
| [`24e0c8042`](https://github.com/antmicro/jetson-nano-baseboard/commit/24e0c804251512667eff9ca695d4db558362e627) | 2023-05-05 | clearance, generic fix | Fix DRC errors. 3 connection TBD | 5051101097.kicad_mod, jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_pro, je |
| [`4776690be`](https://github.com/antmicro/jetson-nano-baseboard/commit/4776690be4e736b4b705189d4c80977171533b2e) | 2023-05-05 | silkscreen | Update VBUS_DET label | display.kicad_sch |
| [`97d178760`](https://github.com/antmicro/jetson-nano-baseboard/commit/97d1787601def91e3eba46f400ad223467d13fe1) | 2023-05-05 | generic fix | Fix SIC477ED note | supply.kicad_sch |
| [`8ff244c3b`](https://github.com/antmicro/jetson-nano-baseboard/commit/8ff244c3b3b3d41f6985b19d2663c7f5bf1b2c22) | 2023-05-05 | routing | Route OTG USB signals. Connection to SOM tbd | jetson-nano-baseboard.kicad_pcb |
| [`92fcb6d86`](https://github.com/antmicro/jetson-nano-baseboard/commit/92fcb6d8614422d2fa459f86799969bee3bd0dfa) | 2023-05-05 | routing | Route I2C and CAN to external connectors | jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_pro |
| [`6efde0b6c`](https://github.com/antmicro/jetson-nano-baseboard/commit/6efde0b6c9ab19b4684fbb56eaea2c58e067d069) | 2023-05-05 | routing | Improve ground plane under USB 3.x signals | jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_pro |
| [`b8e1c0701`](https://github.com/antmicro/jetson-nano-baseboard/commit/b8e1c0701bef08c0d275774f6619ade4d1151dd7) | 2023-05-05 | part swap, silkscreen | Unify reference designator size. Hide values | jetson-nano-baseboard.kicad_pcb |
| [`876aa1cbb`](https://github.com/antmicro/jetson-nano-baseboard/commit/876aa1cbb9bbf236da29c36d85bf158b26c269f6) | 2023-05-05 | routing | Improve USB 3.x routing. Route as DIFF 85Ohm on the top layer | jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_pro |
| [`9b66ac02c`](https://github.com/antmicro/jetson-nano-baseboard/commit/9b66ac02cb959d58b9267c517006aa219f3681b7) | 2023-05-05 | part swap | Sync PCB with schematics, place I2C/CAN connectors and part of the passives | jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.kicad_pro |
| [`72660a7c6`](https://github.com/antmicro/jetson-nano-baseboard/commit/72660a7c684fac8c9ab6fda3443fff5e2e1f6f88) | 2023-05-05 | protection | Improved input protection - added TVSes, reverse polarity protection and input fuse | DO-214AC.kicad_mod, SOD-523.kicad_mod, jetson-nano-baseboard.kicad_sch, jetson-nano-basebo |
| [`6c263a6cd`](https://github.com/antmicro/jetson-nano-baseboard/commit/6c263a6cd4d4150fdfd2d866694c20370c38f88e) | 2023-05-05 | power | Added resistor between RBOOT and CBOOT on 3V3_SYS regulator for reduced EMI | jetson-nano-baseboard.kicad_sch, jetson-nano-baseboard.kicad_sym, supply.kicad_sch |
| [`f42b94f7c`](https://github.com/antmicro/jetson-nano-baseboard/commit/f42b94f7c26484a1770e6d8b33be26a44528ee66) | 2023-05-05 | part swap | Replaced 742792022 ferrites with 74279220181 for higher current capability | jetson-nano-baseboard.kicad_sch, jetson-nano-baseboard.kicad_sym, som.kicad_sch, supply.ki |
| [`1e5fa23b7`](https://github.com/antmicro/jetson-nano-baseboard/commit/1e5fa23b7d73c7fbc95cda526a237c606a8a2651) | 2023-05-05 | generic fix | Fixed schematic connection between USB OTG shield and L7 | jetson-nano-baseboard.kicad_sch, usb-debug.kicad_sch |
| [`02db42c92`](https://github.com/antmicro/jetson-nano-baseboard/commit/02db42c92ef00b2e3029581728c9ba018c2531ac) | 2022-04-06 | part swap | Replace AP22652AW6-7 with AP22653W6-7 | jetson-nano-baseboard-cache.lib, jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.li |
| [`f97a3bffc`](https://github.com/antmicro/jetson-nano-baseboard/commit/f97a3bffcf9c42aed05856647e3aa5257385654d) | 2022-03-22 | pull-up/down | Add i2c-mux reset pull-up | csi.sch, jetson-nano-baseboard.kicad_pcb |
| [`d69e69af0`](https://github.com/antmicro/jetson-nano-baseboard/commit/d69e69af0a7d00a6882b64f7e5fb726bcf5ab92f) | 2022-03-22 | pull-up/down | Add T7 pull-down | jetson-nano-baseboard-cache.lib, jetson-nano-baseboard.kicad_pcb, usb-debug.sch |
| [`0af035427`](https://github.com/antmicro/jetson-nano-baseboard/commit/0af03542778880ea9f4cafef99b627e38ede9202) | 2022-02-25 | generic fix | Fix connections between resistors and DCDCs | jetson-nano-baseboard.kicad_pcb, supply.sch |
| [`82952e21b`](https://github.com/antmicro/jetson-nano-baseboard/commit/82952e21bfa8301becfbb3094c1d0dd578afcb2c) | 2022-02-25 | routing, generic fix | Fix ethernet connector CT and GND connections | interfaces.sch, jetson-nano-baseboard-cache.lib, jetson-nano-baseboard.kicad_pcb |
| [`6d95efe60`](https://github.com/antmicro/jetson-nano-baseboard/commit/6d95efe608778dc5ff80b0fdff504cc6e047da55) | 2022-02-07 | generic fix | Fix annotation on J10 | som.sch |
| [`28c78351c`](https://github.com/antmicro/jetson-nano-baseboard/commit/28c78351c6eed547dbf25f430c137c8f077533f1) | 2022-02-04 | silkscreen | Add pin-1 marking to CSI connectors | jetson-nano-baseboard.kicad_pcb |
| [`f8a23c5b9`](https://github.com/antmicro/jetson-nano-baseboard/commit/f8a23c5b9af8c09e5321912ad32f400a0eeb56c3) | 2022-02-04 | silkscreen | Tweak designators | jetson-nano-baseboard.kicad_pcb |
| [`d7057dc5e`](https://github.com/antmicro/jetson-nano-baseboard/commit/d7057dc5e8d90708654fb41180706a3738d73098) | 2022-02-04 | clearance | Update layout after couryard tunning, DRC-Clean | jetson-nano-baseboard.kicad_pcb |
| [`091261905`](https://github.com/antmicro/jetson-nano-baseboard/commit/09126190540ee773eb62a7ec89a2421f01f1b87b) | 2022-02-03 | footprint | Resolve overlapping courtyards | jetson-nano-baseboard.kicad_pcb |
| [`272c27842`](https://github.com/antmicro/jetson-nano-baseboard/commit/272c2784220087d439c24ab92c1b5fcb8e64abff) | 2022-01-07 | silkscreen | Improving silks | jetson-nano-baseboard.kicad_pcb |
| [`82f1249ad`](https://github.com/antmicro/jetson-nano-baseboard/commit/82f1249ad6add353abd7754528ebd29432d62a65) | 2022-01-07 | silkscreen, power, generic fix | Fix power zone with thermal reliefs, minor silks cleanup | jetson-nano-baseboard.kicad_pcb |
| [`dcac9b3fb`](https://github.com/antmicro/jetson-nano-baseboard/commit/dcac9b3fb5cd91085cd0a6ca0c4e5885c1486379) | 2022-01-05 | clearance, routing, generic fix | Fix spacing issue on diff pair | jetson-nano-baseboard.kicad_pcb |
| [`e9ceebf76`](https://github.com/antmicro/jetson-nano-baseboard/commit/e9ceebf76d1e4d916aa525199e965c093793127f) | 2022-01-04 | generic fix | Fix errors from EC | jetson-nano-baseboard.kicad_pcb |
| [`8b75c51d3`](https://github.com/antmicro/jetson-nano-baseboard/commit/8b75c51d38bec3d10f5cac835083ebb725c41778) | 2022-01-04 | generic fix | Adjust for EC errors | jetson-nano-baseboard.kicad_pcb |
| [`1908e5171`](https://github.com/antmicro/jetson-nano-baseboard/commit/1908e5171636d03018aad551711e5145a585a18f) | 2022-01-03 | silkscreen | Finish routing, cleanup silks | jetson-nano-baseboard-cache.lib, jetson-nano-baseboard.kicad_pcb, jetson-nano-baseboard.pr |
| [`d4f269ed6`](https://github.com/antmicro/jetson-nano-baseboard/commit/d4f269ed6ad0ffeeee4a6578473f239f38a93950) | 2021-12-31 | part swap, routing, generic fix | Replace ethernet connector, fix traces and tuning | SI-51005-F.kicad_mod, interfaces.sch, jetson-nano-baseboard-cache.lib, jetson-nano-baseboa |
| [`9fa7e8f01`](https://github.com/antmicro/jetson-nano-baseboard/commit/9fa7e8f01e4f1e6a76c9d13222ed38a2e010551d) | 2021-12-31 | silkscreen, generic fix | Fix placement of some components, minor fsilk changes | jetson-nano-baseboard.kicad_pcb |
| [`8971d6565`](https://github.com/antmicro/jetson-nano-baseboard/commit/8971d6565f3a12b72e4eac1b86802c86e1ef626c) | 2021-12-31 | generic fix | Minor placement fixes | jetson-nano-baseboard.kicad_pcb |
| [`3ef4cda63`](https://github.com/antmicro/jetson-nano-baseboard/commit/3ef4cda635b65b5a5b8ac3b88eb916aee20ef27b) | 2021-12-31 | routing | Move OSHW logo away from track | jetson-nano-baseboard.kicad_pcb |
| [`04d51dd49`](https://github.com/antmicro/jetson-nano-baseboard/commit/04d51dd4909115b6e1e2c415c63f00aa659bb68d) | 2021-12-31 | silkscreen | Move silks back into PCB | jetson-nano-baseboard.kicad_pcb |
| [`98ef74486`](https://github.com/antmicro/jetson-nano-baseboard/commit/98ef744864b94363a60785f969f235802d9a761e) | 2021-12-31 | routing, generic fix | Layout 5V DCDC, check traces, fix ffab references size | jetson-nano-baseboard.kicad_pcb |
| [`64baa546a`](https://github.com/antmicro/jetson-nano-baseboard/commit/64baa546abf92f6cbf3ab1862a6891ee42c2e7b7) | 2021-12-31 | generic fix | Rework 3V3 and 5V DCDC on schematic, layout 3V3 DCDC | 0603-Fuse-Littelfuse-0467.kicad_mod, Coilcraft-XEL4030.kicad_mod, HVQFN24.kicad_mod, VQFN- |
| [`b452ba3cd`](https://github.com/antmicro/jetson-nano-baseboard/commit/b452ba3cdd67541de097abb871098054fd485526) | 2021-12-31 | routing, generic fix | Fix traces after moving ethernet | jetson-nano-baseboard.kicad_pcb |
| [`edd21b331`](https://github.com/antmicro/jetson-nano-baseboard/commit/edd21b331197cd3f275674c1f5442dfec37ecd71) | 2021-12-31 | part swap | Grouping parts around DCDC | jetson-nano-baseboard.kicad_pcb |
| [`79998184e`](https://github.com/antmicro/jetson-nano-baseboard/commit/79998184ecb2a48c4daad73ad514db3efdaa950a) | 2021-12-31 | decoupling, generic fix | Fix symbols, change DCDCC output capacitors to 0805 | supply.sch |
| [`9352feda6`](https://github.com/antmicro/jetson-nano-baseboard/commit/9352feda6ca04b827bf123bc7261935b719f89b0) | 2021-12-31 | part swap | Replace 5V0 and 3V3 DCDC | supply.sch |
| [`ec352887c`](https://github.com/antmicro/jetson-nano-baseboard/commit/ec352887c7cff0da17949874f2bd5e2e88b45e60) | 2021-12-31 | silkscreen | Update default silk size | jetson-nano-baseboard.pro |
| [`a4d03864d`](https://github.com/antmicro/jetson-nano-baseboard/commit/a4d03864d48c3fac905f2b2e1ca45557ebb8a46c) | 2021-12-31 | clearance, generic fix | Fix DRC errors | jetson-nano-baseboard.kicad_pcb |
| [`23c125eba`](https://github.com/antmicro/jetson-nano-baseboard/commit/23c125ebaa41ad8d0f673b654e72228b72314c38) | 2021-12-31 | silkscreen, protection | Add T7 to reverse logic, cleanup silks | jetson-nano-baseboard.kicad_pcb |
| [`6937e7a7a`](https://github.com/antmicro/jetson-nano-baseboard/commit/6937e7a7ad53ea0fb07c43915abf091ff6b2ed96) | 2021-12-31 | part swap | Replace NCP380HS with AP22652W6-7 | usb-debug.sch |
| [`4fe7b9d47`](https://github.com/antmicro/jetson-nano-baseboard/commit/4fe7b9d476a46cad2901b03314c19ef7cef68d85) | 2021-12-31 | footprint | Changed package of SN74LV1T125 to SOT32 | display.sch |
| [`b2a29cabd`](https://github.com/antmicro/jetson-nano-baseboard/commit/b2a29cabdd6ed25473778189c689616e9b88987d) | 2021-12-31 | clearance | Import DRC rules | jetson-nano-baseboard-cache.lib, jetson-nano-baseboard.dcm, jetson-nano-baseboard.kicad_pc |
| [`c7feb21e1`](https://github.com/antmicro/jetson-nano-baseboard/commit/c7feb21e1acc4f26ef413879dd49f37a975334ec) | 2021-12-31 | part swap | Clean copied parts | display.sch |
| [`998929f92`](https://github.com/antmicro/jetson-nano-baseboard/commit/998929f9248c28cc34c37cc91a2c4529a4897f9c) | 2021-06-09 | silkscreen | Clean up B.SilkS layer | jetson-nano-baseboard.kicad_pcb |
| [`069f053de`](https://github.com/antmicro/jetson-nano-baseboard/commit/069f053de0053c97fd7302e3a51b8b0925c8030b) | 2021-06-09 | clearance, generic fix | Fix DRC errors | jetson-nano-baseboard.kicad_pcb, usb-debug.sch |
| [`31117ba53`](https://github.com/antmicro/jetson-nano-baseboard/commit/31117ba53c7088f2a98d18c0026a90d1bcf027ad) | 2021-06-09 | silkscreen | Update 3d models and clean up designators | jetson-nano-baseboard.kicad_pcb |
| [`46726cd7d`](https://github.com/antmicro/jetson-nano-baseboard/commit/46726cd7d00735b91aa969491fd561ba4248fcad) | 2021-06-09 | generic fix | Fix R113 connection | som.sch |
| [`8016f5a2f`](https://github.com/antmicro/jetson-nano-baseboard/commit/8016f5a2fa41b763edbf86b466c9cb74a734e010) | 2020-03-06 | generic fix | Fix PrjPcb | jetson-nano-baseboard.PrjPcb, usb-debug.SchDoc |
| [`9ec2a8a5e`](https://github.com/antmicro/jetson-nano-baseboard/commit/9ec2a8a5efde07806b7702a3a0cbfdde592e7fa4) | 2019-04-17 | pull-up/down, footprint, generic fix | Updated pads for Displayport connector Corrected I2C connection for camera interface Added Q1 that inverts con | Jetson-Nano-Baseboard-BOM.csv, PCB.PcbDoc, SoM.SchDoc, csi.SchDoc, display.SchDoc, etherne |
