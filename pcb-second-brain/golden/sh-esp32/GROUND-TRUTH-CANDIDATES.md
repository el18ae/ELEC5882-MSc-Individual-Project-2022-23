# Ground-truth candidates — sh-esp32

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/hatlabs/SH-ESP32-hardware
- Pinned commit: `167981ac8c3e1a7a0104c20cf9b487b4d722fc60` (2024-10-16T14:35:59+03:00)
- Paths analysed: `.`
- Domain: industrial_comms
- Commits touching these paths: 183; touching design files: 162; design-file commits whose message matches a fix keyword: 108

## Revision pairs

No revision directories; revisions exist only as git history (pairs = consecutive commits below).

Tags in repository (10): `0.1.0`, `0.2.0`, `0.3.0`, `0.3.0-2`, `0.3.1`, `latest`, `2.0.0`, `2.0.1`, `v2.0.2`, `v2.2.1`

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [66](https://github.com/hatlabs/SH-ESP32-hardware/issues/66) | open | 2026-10-28 | LICENSE.md says CC BY 4.0, design files say CC BY-SA 4.0 |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`e26fdf9e6`](https://github.com/hatlabs/SH-ESP32-hardware/commit/e26fdf9e6a120ce24da2c20aa3f798da7036f65d) | 2024-10-16 | part swap | Adjust parts according to availability | SH-ESP32.kicad_pcb, power-input.kicad_sch |
| [`c03c4e8f8`](https://github.com/hatlabs/SH-ESP32-hardware/commit/c03c4e8f86376c71f6ec3ab31b830e8fc9682b7a) | 2024-10-16 | generic fix | Fix copper fills | SH-ESP32.kicad_pcb |
| [`4d8977890`](https://github.com/hatlabs/SH-ESP32-hardware/commit/4d897789075562bbaf448f3babb7bb949cb47084) | 2024-09-24 | power | Change buck converter to LMR54410 | ESP32-WROOM-32D.kicad_mod, ESP32.kicad_sch, I2C.kicad_sch, PCB.kicad_sch, SH-ESP32.kicad_p |
| [`7f29cec9c`](https://github.com/hatlabs/SH-ESP32-hardware/commit/7f29cec9cd76f3746352344ff86954fd276b1885) | 2024-01-29 | part swap | Change due to availability | optocouplers.kicad_sch |
| [`3a23d65d6`](https://github.com/hatlabs/SH-ESP32-hardware/commit/3a23d65d61c8ec8fac1e6c4db99e0304c548acec) | 2023-08-10 | generic fix | Fix component rotations | ESP32.kicad_sch, I2C.kicad_sch, SH-ESP32.kicad_pcb, USB.kicad_sch, canbus.kicad_sch, onewi |
| [`f8b8e844a`](https://github.com/hatlabs/SH-ESP32-hardware/commit/f8b8e844a4703d2cae772ed7b1b07be556f7ebac) | 2023-08-10 | part swap | Update BAT54AW part # for availability | ESP32.kicad_sch, SH-ESP32.kicad_pcb |
| [`4a3458026`](https://github.com/hatlabs/SH-ESP32-hardware/commit/4a3458026710b241a3103c42b0790843048c40bf) | 2023-08-08 | decoupling, routing | Bulk-resize 0.2 mm drill vias to 0.3 mm | SH-ESP32.kicad_pcb |
| [`db153efd4`](https://github.com/hatlabs/SH-ESP32-hardware/commit/db153efd4cc1f2fe703670c0232853e0d996afca) | 2022-10-26 | part swap, power | Change the filter inductor according to JLC availability | SH-ESP32.kicad_sch, power-input.kicad_sch |
| [`62274d63d`](https://github.com/hatlabs/SH-ESP32-hardware/commit/62274d63dd8aa1bab20625c62fbc19bd006948a2) | 2022-10-26 | part swap | Change parts according to JLCPCB availability | SH-ESP32.kicad_pcb, SH-ESP32.kicad_sch, power-input.kicad_sch |
| [`8213cc201`](https://github.com/hatlabs/SH-ESP32-hardware/commit/8213cc201d2f24646b839d81d5a91c1c4b15861f) | 2022-10-26 | protection | Improve connector ESD protection; bump version to 2.0.2 | ESP32.kicad_sch, I2C.kicad_sch, PCB.kicad_sch, SH-ESP32.kicad_pcb, SH-ESP32.kicad_pro, SH- |
| [`e114cdd14`](https://github.com/hatlabs/SH-ESP32-hardware/commit/e114cdd1405fda9d48a94d89dc65a27a8c965ea0) | 2022-07-07 | protection | Only use the edge signal pins of RCLAMPs if possible | I2C.kicad_sch, SH-ESP32.kicad_pcb, onewire.kicad_sch |
| [`b5b8e4339`](https://github.com/hatlabs/SH-ESP32-hardware/commit/b5b8e4339867efd690a4244c55c891c85d27f382) | 2022-02-07 | pull-up/down | Add a pull-up resistor to CAN TX line | SH-ESP32.kicad_pcb, SH-ESP32.kicad_pro, SH-ESP32.kicad_sch, canbus.kicad_sch |
| [`dd8633d15`](https://github.com/hatlabs/SH-ESP32-hardware/commit/dd8633d1539c680a82aea6de76528451484b29b0) | 2022-01-25 | silkscreen | Hide the ref. designator of J702 | SH-ESP32.kicad_pcb |
| [`6b9b51c99`](https://github.com/hatlabs/SH-ESP32-hardware/commit/6b9b51c99dabfc06103a8c2ff45bd60a15bf1cc9) | 2021-12-31 | part swap | One more component swap due to part availability | ESP32.kicad_sch, SH-ESP32.kicad_sch, power-input.kicad_sch |
| [`b13498f07`](https://github.com/hatlabs/SH-ESP32-hardware/commit/b13498f07b6b07fe05819ab392c5ea2b41c5fcce) | 2021-12-31 | clearance, generic fix | Fix DRC errors reported by KiKit CLI | SH-ESP32.kicad_pcb |
| [`952513401`](https://github.com/hatlabs/SH-ESP32-hardware/commit/952513401772cc5706fdf5260cdb9d568ccc6aca) | 2021-12-31 | clearance, generic fix | Fix PCB DRC errors | SH-ESP32.kicad_pcb |
| [`2597f6fc5`](https://github.com/hatlabs/SH-ESP32-hardware/commit/2597f6fc5dcd08cb3f3f2e48f961fd520d922fef) | 2021-12-31 | part swap | Bump version to 2.0.1, update the components for availability | ESP32.kicad_sch, I2C.kicad_sch, PCB.kicad_sch, SH-ESP32.kicad_pcb, SH-ESP32.kicad_sch, UI. |
| [`d91e83660`](https://github.com/hatlabs/SH-ESP32-hardware/commit/d91e8366020c82dfe4b70e8b0618a4771ec63ecd) | 2021-11-14 | part swap | Change pinheader PN for availability | onewire.sch |
| [`bdcd2a292`](https://github.com/hatlabs/SH-ESP32-hardware/commit/bdcd2a292f589e9a990700136b52d8e73e3afcaf) | 2021-11-14 | decoupling, part swap, power, generic fix | Fix buck FF capacitor value and part number | SH-ESP32.kicad_pcb, power-input.sch |
| [`139e849f9`](https://github.com/hatlabs/SH-ESP32-hardware/commit/139e849f9de441a3eb76ce986d3f12a450c86f7a) | 2021-08-17 | power | Add provisions for adjustable buck converter | SH-ESP32.kicad_pcb, power-input.sch |
| [`d9c9cd862`](https://github.com/hatlabs/SH-ESP32-hardware/commit/d9c9cd862ac411353526990bcdafd1d5219eaad6) | 2021-08-16 | part swap | Change CAN Vcc to 3.3V; change components for availability | ESP32.sch, SH-ESP32.kicad_pcb, canbus.sch, power-input.sch |
| [`9d3474e92`](https://github.com/hatlabs/SH-ESP32-hardware/commit/9d3474e9209f75c653f62dc90081f62e5552ba50) | 2021-08-13 | part swap, power | Replace buck converter and CAN transceiver ICs | ESP32.sch, I2C.sch, PCB.sch, SH-ESP32.kicad_pcb, SH-ESP32.lib, SMMS0650-1R0M.kicad_mod, ca |
| [`947697097`](https://github.com/hatlabs/SH-ESP32-hardware/commit/94769709764022293121e90567c75fda2687941a) | 2021-04-26 | part swap | USB conx is now available at JLCPCB - make non-virtual | SH-ESP32.kicad_pcb |
| [`2fb13fdb5`](https://github.com/hatlabs/SH-ESP32-hardware/commit/2fb13fdb57bc63273e78a103ce0da9694e56b0a6) | 2021-04-26 | part swap | Replace L302 for better availability | SH-ESP32.kicad_pcb, power-input.sch |
| [`10fff1634`](https://github.com/hatlabs/SH-ESP32-hardware/commit/10fff1634d43a41c0cd5cdd5f7644d03b3df92bc) | 2021-04-26 | part swap | Swap BAT54As to SOT-323 to improve availability | ESP32.sch, SH-ESP32.kicad_pcb |
| [`4d157b7ef`](https://github.com/hatlabs/SH-ESP32-hardware/commit/4d157b7ef884c74ce98b4954747a71603bc77c76) | 2021-04-26 | decoupling | Move bypass caps further away from mounting holes | SH-ESP32.kicad_pcb |
| [`e5f51c567`](https://github.com/hatlabs/SH-ESP32-hardware/commit/e5f51c567a1098b79eca542069b5f089e24562cc) | 2021-04-08 | silkscreen | Clarify the silkscreen for the 2x3 pin power conx | SH-ESP32.kicad_pcb, SH-ESP32.pro |
| [`c5f9333de`](https://github.com/hatlabs/SH-ESP32-hardware/commit/c5f9333dea0eb8460584a91185db6868e8933d5c) | 2021-03-04 | part swap, generic fix | Fix 2-pin header part number | power-input.sch |
| [`cef299259`](https://github.com/hatlabs/SH-ESP32-hardware/commit/cef299259f2bc786d7f21f850565b5c47ac46bbb) | 2021-02-21 | footprint, clearance | Remove the keep-out zone from ESP32-WROOM-32 footprint | ESP32-WROOM-32-HandSolder.kicad_mod, SH-ESP32.kicad_pcb |
| [`92ab7da5b`](https://github.com/hatlabs/SH-ESP32-hardware/commit/92ab7da5b28afb5e188ffa3f234dc6a85b478d9b) | 2021-02-18 | decoupling | Update the 1-wire RC cap on the PCB as well | SH-ESP32.kicad_pcb |
| [`4e68ec4b3`](https://github.com/hatlabs/SH-ESP32-hardware/commit/4e68ec4b3ee54056c7e8374fab9af65f8323ac01) | 2021-02-18 | decoupling, part swap | Reduce 1-wire RC filter capacitor value | onewire.sch |
| [`cf5b57911`](https://github.com/hatlabs/SH-ESP32-hardware/commit/cf5b579111f2a25060c52ba12d04fdc8516ac557) | 2021-02-18 | part swap | Swap 1-Wire and QWICC Connectors so 1-wire is on the edge, and can access when LCD on. | SH-ESP32.kicad_pcb |
| [`924ddbb95`](https://github.com/hatlabs/SH-ESP32-hardware/commit/924ddbb954771e35adb8b69e35c1e00ac7dc5f50) | 2021-02-17 | decoupling, clearance, generic fix | Fix a bypass cap location (DRC error) | SH-ESP32.kicad_pcb |
| [`355836d70`](https://github.com/hatlabs/SH-ESP32-hardware/commit/355836d7015001ee28816fada709b5ee38de6724) | 2021-02-17 | routing | Tweak the plane fills around the PSU-CAN gap | SH-ESP32.kicad_pcb |
| [`fc8997cc6`](https://github.com/hatlabs/SH-ESP32-hardware/commit/fc8997cc60b90a88bc4af545c26a01a6719cddf5) | 2021-02-17 | silkscreen | Nudge the bottom silk screen items around | SH-ESP32.kicad_pcb |
| [`6d2cb0beb`](https://github.com/hatlabs/SH-ESP32-hardware/commit/6d2cb0beb3ba96a5138369dc970c7341111622e3) | 2021-02-17 | part swap, routing | Replace the 3.3V plane under the PSU section with a ground fill | SH-ESP32.kicad_pcb |
| [`61e45fb56`](https://github.com/hatlabs/SH-ESP32-hardware/commit/61e45fb568e4799ba051dfb7dcc27be8750e5321) | 2021-02-17 | routing | Move L302 ground fill to IN net | SH-ESP32.kicad_pcb |
| [`fb3fc4838`](https://github.com/hatlabs/SH-ESP32-hardware/commit/fb3fc4838e1940178fddaee7a07272cb999d31ef) | 2021-02-15 | silkscreen | Bump revision on silk screens | SH-ESP32.kicad_pcb |
| [`67b831d3e`](https://github.com/hatlabs/SH-ESP32-hardware/commit/67b831d3eaef1e034bfde76873234bd6713b0817) | 2021-02-15 | silkscreen, routing, generic fix | Fix back silkscreen, remove a redundant via | SH-ESP32.kicad_pcb |
| [`e81e7dedd`](https://github.com/hatlabs/SH-ESP32-hardware/commit/e81e7dedd46680dbb0caac05dbd676b55190553a) | 2021-02-15 | part swap, routing | Swap Power and CAN connector pins to have GND on pin 1 | SH-ESP32.kicad_pcb, canbus.sch, power-input.sch |
| [`de7a42d84`](https://github.com/hatlabs/SH-ESP32-hardware/commit/de7a42d840629e4755c9cc76817dd23a85907586) | 2021-02-15 | power | Change the buck converter inductor back to a shielded one | SH-ESP32.kicad_pcb, power-input.sch |
| [`815bd4c94`](https://github.com/hatlabs/SH-ESP32-hardware/commit/815bd4c94d61257afd9b8e5410badb1bc731533d) | 2021-01-31 | part swap | Update parts according to JLC availability | SH-ESP32.kicad_pcb, canbus.sch, optocouplers.sch, power-input.sch |
| [`769c43391`](https://github.com/hatlabs/SH-ESP32-hardware/commit/769c4339108a48b993f50452e670912b0cb0bd9f) | 2021-01-30 | part swap, generic fix | Fix D1101 and D1102 part numbers | ESP32.sch |
| [`de4c6cc81`](https://github.com/hatlabs/SH-ESP32-hardware/commit/de4c6cc816431bc1127f6a1398ccd38a42a8fd90) | 2021-01-29 | power, generic fix | Fix regulator layout according to comments | SH-ESP32.kicad_pcb, SH-ESP32.pro |
| [`45f96fe22`](https://github.com/hatlabs/SH-ESP32-hardware/commit/45f96fe227e4552d63148ef31dc1df7c88d90f4a) | 2021-01-28 | power | Add an LC pi filter to buck converter input | ESP32.sch, I2C.sch, PCB.sch, SH-ESP32.kicad_pcb, SH-ESP32.sch, UI.sch, USB.sch, canbus.sch |
| [`483b15795`](https://github.com/hatlabs/SH-ESP32-hardware/commit/483b1579585edfc8b15fd4627a725d9a95c8ae81) | 2021-01-18 | silkscreen, generic fix | Fix Revision on back silkscreen | SH-ESP32.kicad_pcb |
| [`22d9f6e53`](https://github.com/hatlabs/SH-ESP32-hardware/commit/22d9f6e53f526955884eeb85b403fd0a85f56099) | 2021-01-18 | silkscreen | Added GPIO Ledged to rear silkscreen | SH-ESP32.kicad_pcb |
| [`e208fe071`](https://github.com/hatlabs/SH-ESP32-hardware/commit/e208fe071b2cb694eefa40aed4b28736530eca1a) | 2021-01-18 | protection | Remove RClamp on CAN bus input | SH-ESP32.kicad_pcb, canbus.sch |
| [`3e986b96d`](https://github.com/hatlabs/SH-ESP32-hardware/commit/3e986b96d30ccd100a16c8ff69a1230562e2f40e) | 2021-01-04 | silkscreen, protection | Moved Fuse etc on CAN-VIN to better show  Connector Silkscreen. | SH-ESP32.kicad_pcb |
| [`ef0f44f82`](https://github.com/hatlabs/SH-ESP32-hardware/commit/ef0f44f82f2df0faabdfd353228020f28c6576b2) | 2021-01-04 | part swap, protection, routing, generic fix | Power system bug fix. Changed fuses & Diode - Swapping incorrect part  D1101 and D1102  to BAT54A   /  LCSC:   | ESP32.sch, I2C.sch, PCB.sch, SH-ESP32.kicad_pcb, SH-ESP32.pro, SH-ESP32.sch, UI.sch, USB.s |
| [`7b9c3091c`](https://github.com/hatlabs/SH-ESP32-hardware/commit/7b9c3091cf589efea21167af46a5573ceba36665) | 2020-12-22 | routing | Move N2K Connector to avoid pins hitting case PCB Standoff (underneath). - Increase mask on hole for screw fla | SH-ESP32.kicad_pcb |
| [`258f82820`](https://github.com/hatlabs/SH-ESP32-hardware/commit/258f82820b138dae2eae7044fa219efcdab5b03f) | 2020-12-19 | clearance, routing | Groundpllane on Hat-Proto & moved KeepOut on Blank to user Cmts | SH-ESP32-Blank.kicad_pcb, SH-ESP32-HAT-Proto.kicad_pcb |
| [`0c830aca7`](https://github.com/hatlabs/SH-ESP32-hardware/commit/0c830aca76d9afeec00b2925e98dd1d0af37bcfd) | 2020-12-19 | footprint, clearance, generic fix | Change 'unconnected' holes to Square pads - Fixed DRC errors | IOHeader.sch, SH-ESP32-HAT-Proto.kicad_pcb |
| [`f141d8800`](https://github.com/hatlabs/SH-ESP32-hardware/commit/f141d88008920113ced420b7eb9abecf4e390617) | 2020-12-19 | generic fix | Update HAT-Proto to V0.3  Pin header locs and fix outline for box corners | IOHeader.sch, PCB.sch, SH-ESP32-HAT-Proto.kicad_pcb |
| [`253df3696`](https://github.com/hatlabs/SH-ESP32-hardware/commit/253df36960fff07a379bdfe7139422da55a194a6) | 2020-12-19 | clearance, generic fix | Fixed spacing on  TXD/RXD  0-Ohm Links | SH-ESP32.kicad_pcb |
| [`0a3da929d`](https://github.com/hatlabs/SH-ESP32-hardware/commit/0a3da929d88b67e19699fd78e748489135c56e93) | 2020-12-19 | footprint, silkscreen, routing | Modify ESP32 footprint for Hand Soldering - Extended all pads 1mm outwards - Rearranged 0-Ohm links and EN/CLK | ESP32-WROOM-32-HandSolder.kicad_mod, ESP32.sch, SH-ESP32.kicad_pcb, THTPad_2.1x2.1mm_Drill |
| [`4dd2f9e40`](https://github.com/hatlabs/SH-ESP32-hardware/commit/4dd2f9e40656e3d44f70d4493a992cf0b70683de) | 2020-12-18 | silkscreen, routing | Change mounting holes to slots to suit multiple box vendors with variable mounting positions. - Move GPIO Conn | MountingSlot_3.2mm_M3_Unplated.kicad_mod, PCB.sch, SH-ESP32.kicad_pcb |
| [`595012df6`](https://github.com/hatlabs/SH-ESP32-hardware/commit/595012df62a085448ebf53ec428145049c89b0d7) | 2020-12-17 | part swap | Change an NPN LCSC part no. according to availability | SH-ESP32.kicad_pcb, optocouplers.sch |
| [`087e718b7`](https://github.com/hatlabs/SH-ESP32-hardware/commit/087e718b754a06968ef21002377eeda642f7538d) | 2020-12-17 | clearance, power | Increase thermal clearances of copper fills to improve solderability | SH-ESP32.kicad_pcb, canbus.sch |
| [`a40410177`](https://github.com/hatlabs/SH-ESP32-hardware/commit/a40410177622eb79ba003a3186775eab45927504) | 2020-12-03 | silkscreen | - Removed Extra 'NEMA' pinout from mid board silkscreen - Cutout board under ESP32 Antenna as per HW Design Sp | SH-ESP32.kicad_pcb |
| [`f378d2264`](https://github.com/hatlabs/SH-ESP32-hardware/commit/f378d226480e3fb48e3393baf5926e6f765781a4) | 2020-12-03 | silkscreen | Switch power and  NEMA to JST-XH Connectors and update silkscreen accordingly | SH-ESP32.kicad_pcb, canbus.sch, power-input.sch |
| [`078200baf`](https://github.com/hatlabs/SH-ESP32-hardware/commit/078200baf6c498f5ca909c18d46768ed51c40c8f) | 2020-11-30 | footprint, part swap, generic fix | Fix C603 part number to match the footprint | SH-ESP32.kicad_pcb, canbus.sch |
| [`bcad36f75`](https://github.com/hatlabs/SH-ESP32-hardware/commit/bcad36f75748519ad1f412bdd872a23657457a1e) | 2020-11-29 | part swap | Update pin header LCSC part numbers | ESP32.sch |
| [`443cbd23b`](https://github.com/hatlabs/SH-ESP32-hardware/commit/443cbd23b7aef30765dfa3d6793d8f215149dcc0) | 2020-11-27 | part swap | Swap D602 to B5819W which is a standard component | SH-ESP32.kicad_pcb, canbus.sch |
| [`44d726282`](https://github.com/hatlabs/SH-ESP32-hardware/commit/44d72628277900b21528ee586b41d3b151b878bc) | 2020-11-26 | routing | Adjust track width of Hi_Freq traces for Ethernet.   Reduce CLK to 2mm trace | SH-ESP32.kicad_pcb |
| [`21192ed93`](https://github.com/hatlabs/SH-ESP32-hardware/commit/21192ed93daacf9a445133ce42e13148f0bdfb83) | 2020-11-26 | part swap, power | Change buck converter inductor to something available at JLCPCB | SH-ESP32.kicad_pcb, power-input.sch |
| [`e3fb75389`](https://github.com/hatlabs/SH-ESP32-hardware/commit/e3fb753891d1fe00b6715e4e099fda8c29ef156e) | 2020-11-26 | part swap, silkscreen, generic fix | - Fixed LCSC part number of USB Connector - Change CPU to ESP32-WROOM32E  part (was D)  and fixed LCSC part nu | ESP32.sch, SH-ESP32-Blank.kicad_pcb, SH-ESP32-Blank.sch, USB.sch |
| [`ebe4d3fea`](https://github.com/hatlabs/SH-ESP32-hardware/commit/ebe4d3fea3002a19fb8cef5817c1a4fe6af6f2a8) | 2020-11-26 | generic fix | Convert USB to Vertical Thru-Hole - Fix Git-ignore for ALL cache files | SH-ESP32.kicad_pcb, USB.sch |
| [`be84eb7ce`](https://github.com/hatlabs/SH-ESP32-hardware/commit/be84eb7ce17bded1b2e91d5538933e729a14cb71) | 2020-11-25 | footprint | Make the USB 5V conx a single-pin test pad for safety | ESP32.sch, SH-ESP32.kicad_pcb |
| [`3ea4bf291`](https://github.com/hatlabs/SH-ESP32-hardware/commit/3ea4bf29191e0891652b6dabed7579a3da132bcc) | 2020-11-25 | routing | ReRoute EN track to be cleaner / away from Antenna | SH-ESP32.kicad_pcb |
| [`490544055`](https://github.com/hatlabs/SH-ESP32-hardware/commit/49054405504655dbe16f882756d27bdeceb8f2d5) | 2020-11-25 | pull-up/down, footprint | Add an unpopulated pad for extra opto output pullup | PCB.sch, SH-ESP32.kicad_pcb, optocouplers.sch |
| [`80ce84620`](https://github.com/hatlabs/SH-ESP32-hardware/commit/80ce846200d22cdd0592d0bc74831b505b258c26) | 2020-11-25 | part swap, silkscreen | Swap power plug pins and reorganize components and silk screen around the area | SH-ESP32.kicad_pcb, power-input.sch |
| [`f2c41f173`](https://github.com/hatlabs/SH-ESP32-hardware/commit/f2c41f173e9115623a92c2e9d4e0904dd9a684b9) | 2020-11-25 | decoupling, silkscreen | Nudged GPIO header labels to not overlap decoupling caps | SH-ESP32.kicad_pcb |
| [`f0192479b`](https://github.com/hatlabs/SH-ESP32-hardware/commit/f0192479b706e2328b889e230733260e3bea4fb9) | 2020-11-25 | routing | Tweak some vias under the ESP32 to improve power plane fill | SH-ESP32.kicad_pcb |
| [`496a8d703`](https://github.com/hatlabs/SH-ESP32-hardware/commit/496a8d703d42451b60d161bf580d7ecdf2ebc0c8) | 2020-11-24 | silkscreen, routing | Added a bunch of vias and tweaked the silkscreens | ESP32.sch, SH-ESP32.kicad_pcb |
| [`f73b71aec`](https://github.com/hatlabs/SH-ESP32-hardware/commit/f73b71aec95cea0def66c51b59875ee714abec71) | 2020-11-24 | routing | Moved Sensor VP to allow Edge plane and added many GND Via's to support top / bottom planes | SH-ESP32.kicad_pcb |
| [`5cbd05c54`](https://github.com/hatlabs/SH-ESP32-hardware/commit/5cbd05c547f5672518e2c5cca774a5d00c456896) | 2020-11-24 | silkscreen | Solder jumper silkscreen | SH-ESP32.kicad_pcb |
| [`d1806b296`](https://github.com/hatlabs/SH-ESP32-hardware/commit/d1806b2964c9465bdc8c77cbac7277e38d4abf6e) | 2020-11-24 | silkscreen, routing | Updated Silkscreen for new GPIO pinout - Moved power header onto 0.1" grid with IO header. - Scattered GND via | ESP32.sch, SH-ESP32.kicad_pcb |
| [`700a372e6`](https://github.com/hatlabs/SH-ESP32-hardware/commit/700a372e6d4ce79b19e603f2fc607c2d2cc132a1) | 2020-11-24 | routing | Initial version with 40Pin header for Ethernet and re-route. | ESP32.sch, SH-ESP32.kicad_pcb, SH-ESP32.sch, UI.sch, USB.sch, ideal_diode.sch, power-input |
| [`57e4bfdfc`](https://github.com/hatlabs/SH-ESP32-hardware/commit/57e4bfdfc1442931ea0fbb57206d0b9827975f09) | 2020-11-23 | footprint, part swap | Changed R905 to 0805 FP and value to 2.2k | SH-ESP32.kicad_pcb, optocouplers.sch |
| [`e4e01f384`](https://github.com/hatlabs/SH-ESP32-hardware/commit/e4e01f38471708d84feebda740ccdb1784a3b661) | 2020-11-20 | routing, generic fix | Brought in all corner arcs by 1.5mm and re-arranged / re-routed accordingly. - Modeled Box outlines on Drawing | SH-ESP32.kicad_pcb |
| [`d62b43d6f`](https://github.com/hatlabs/SH-ESP32-hardware/commit/d62b43d6f446005515337692052ab05e9f9aa1b8) | 2020-11-20 | part swap | - Swapped TXD0 & RXD0 between ESP32 and CH340C (U401) - (at the 0 Ohm Links) - U401 -> Pin 16 (VCC switched to | ESP32.sch, USB.sch |
| [`19a752382`](https://github.com/hatlabs/SH-ESP32-hardware/commit/19a752382c25db95bab6addee54d4101e45aa88b) | 2020-11-19 | part swap | Add LCSC part number to Q1102 | SH-ESP32.kicad_pcb, ideal_diode.sch |
| [`74ac8d2f1`](https://github.com/hatlabs/SH-ESP32-hardware/commit/74ac8d2f11b16daa9dd9ed5b4cbd19bdbe4c0a7b) | 2020-11-19 | power, routing | Tweak the buck converter traces to better conform to datasheet recommendations | SH-ESP32.kicad_pcb |
| [`aded3b7b0`](https://github.com/hatlabs/SH-ESP32-hardware/commit/aded3b7b0376a8c70fa296df6606b126853e87d4) | 2020-11-19 | silkscreen | Rotate silkscreen vertical texts to be oriented down-to-up | SH-ESP32.kicad_pcb |
| [`561cf0116`](https://github.com/hatlabs/SH-ESP32-hardware/commit/561cf01166d0ead3163c23e120d47ab95726bdb0) | 2020-11-18 | silkscreen, routing, generic fix | Final track tidyup and image generation - Ensured all power tracks are consistent 1mm (fixed a couple) - Re-ro | SH-ESP32-cache.lib, SH-ESP32-rescue.dcm, SH-ESP32-rescue.lib, SH-ESP32.kicad_pcb, USB.sch |
| [`caf2925cd`](https://github.com/hatlabs/SH-ESP32-hardware/commit/caf2925cd5c80823a369a17fe0c34e138d4cfec0) | 2020-11-18 | protection | Add a missing 3D shape for the USB fuse | SH-ESP32.kicad_pcb |
| [`b89dea520`](https://github.com/hatlabs/SH-ESP32-hardware/commit/b89dea520b4cac1fedaa6ff15c89e30e5a8e9428) | 2020-11-18 | decoupling | Add a filter cap in front of the USB 5V bead | I2C.sch, SH-ESP32.kicad_pcb, USB.sch, onewire.sch |
| [`0d712fa3c`](https://github.com/hatlabs/SH-ESP32-hardware/commit/0d712fa3cf0cf45b85e8a62b16af4acedd5d3226) | 2020-11-18 | protection, power | Use schottkys instead of MOSFETs for polarity protection | ESP32.sch, I2C.sch, PCB.sch, SH-ESP32-cache.lib, SH-ESP32.dcm, SH-ESP32.kicad_pcb, SH-ESP3 |
| [`b3d9a3fb4`](https://github.com/hatlabs/SH-ESP32-hardware/commit/b3d9a3fb46fde9c7f4691dcd938219ea5cc3d113) | 2020-11-17 | routing | Play with the USB differential pair to make it more straight | SH-ESP32.kicad_pcb, SH-ESP32.pro, USB.sch |
| [`bd462b6dc`](https://github.com/hatlabs/SH-ESP32-hardware/commit/bd462b6dc9cb6b6e85ba8e5c7c0be0b8cafbe65f) | 2020-11-17 | silkscreen, routing, generic fix | Support for OLED/LCD mini displays to direct plug to I2C header - Moved Pin header and Qwiic connectors and re | I2C.sch, SH-ESP32.kicad_pcb |
| [`2340e5af5`](https://github.com/hatlabs/SH-ESP32-hardware/commit/2340e5af53a25d471583aacb18bb7619020c1f50) | 2020-11-16 | silkscreen, routing | Add missing 3d shapes and GND label for J302 | SH-ESP32.kicad_pcb |
| [`68987a1f0`](https://github.com/hatlabs/SH-ESP32-hardware/commit/68987a1f0a26c964241bca93b8a78854258b8dcf) | 2020-11-15 | silkscreen | Added missing JTAG label, Labled JTAG Pins,  changed Header for user I/O Connector to make Pin 1 'round' | PinHeader_1x07_P2.54mm_Vertical-RoundP1.kicad_mod, SH-ESP32.kicad_pcb |
| [`4be154bfe`](https://github.com/hatlabs/SH-ESP32-hardware/commit/4be154bfee494e740d14905f3dc3a97c9457d7de) | 2020-11-07 | generic fix | Fix the library for HT75xx | canbus.sch |
| [`5db33a12d`](https://github.com/hatlabs/SH-ESP32-hardware/commit/5db33a12dec4915e4bca216750c91914ef35adf5) | 2020-11-02 | pull-up/down | Apply pullup change to the PCB design as well | SH-ESP32.kicad_pcb |
| [`5011400e2`](https://github.com/hatlabs/SH-ESP32-hardware/commit/5011400e237b11bbf0d167245ca1bc25f14688c2) | 2020-11-01 | pull-up/down, part swap | Change I2C pullups to 2k due to missing JLCPCB basic parts(!) | I2C.sch |
| [`da9f2016e`](https://github.com/hatlabs/SH-ESP32-hardware/commit/da9f2016ebb7868d0b4294aab8b8f6dd6cf68aff) | 2020-11-01 | part swap | Add a missing LCSC part number for U901 and U902 | optocouplers.sch |
| [`6bba5a80f`](https://github.com/hatlabs/SH-ESP32-hardware/commit/6bba5a80f46adbfe262310c44c962adfb08e3ae6) | 2020-11-01 | footprint | Change D304 footprint | SH-ESP32.kicad_pcb |
| [`2d1e5e34d`](https://github.com/hatlabs/SH-ESP32-hardware/commit/2d1e5e34d2f9384073fec79d842d05fc90df0eea) | 2020-10-30 | routing | Add copious amounts of ground vias | SH-ESP32.kicad_pcb |
| [`372da7c77`](https://github.com/hatlabs/SH-ESP32-hardware/commit/372da7c77eac4ffedd78ba4eca5e5143c3125381) | 2020-10-30 | part swap | Swap optocoupler GPIO pins to connect DAC to output | ESP32.sch, SH-ESP32.kicad_pcb, power-input.sch |
| [`44b17aa5e`](https://github.com/hatlabs/SH-ESP32-hardware/commit/44b17aa5e3c4e70f8bbe7febecc03bc64862e6d0) | 2020-10-30 | routing, generic fix | Change optocoupler types, improve ground planes and a lot of other smaller fixes | ESP32.sch, PCB.sch, SH-ESP32.dcm, SH-ESP32.kicad_pcb, SH-ESP32.lib, USB.sch, canbus.sch, o |
| [`0c3c902c0`](https://github.com/hatlabs/SH-ESP32-hardware/commit/0c3c902c081790fa04674344544e8bf005c3b01e) | 2020-10-28 | silkscreen | Add logo and licensing information + URL on the bottom silkscreen | HL-logo-neg.kicad_mod, MountingHole_3.2mm_M3_Unplated.kicad_mod, SH-ESP32.kicad_pcb, Solde |
| [`a84834fe8`](https://github.com/hatlabs/SH-ESP32-hardware/commit/a84834fe83b5182c5f56c8051f48dfbdb73570ed) | 2020-10-28 | decoupling, routing | Added decoupling caps next to trace vias and copper fill islands | ESP32.sch, I2C.sch, PCB.sch, SH-ESP32.kicad_pcb, SH-ESP32.sch, UI.sch, USB.sch, canbus.sch |
| [`2e7907a14`](https://github.com/hatlabs/SH-ESP32-hardware/commit/2e7907a1454e296a62ea75d4ac811eb3879a50e5) | 2020-10-27 | clearance, routing | Decreased ground fill clearance to fill header gaps | SH-ESP32.kicad_pcb |
| [`26a97296d`](https://github.com/hatlabs/SH-ESP32-hardware/commit/26a97296d349cff8ab7f7262569fe98448c08476) | 2020-10-27 | routing | Rerouted J201 | SH-ESP32.kicad_pcb |
| [`f684bd8ee`](https://github.com/hatlabs/SH-ESP32-hardware/commit/f684bd8eef174af4b3b26ddd844e44f30ca2f2fa) | 2020-10-26 | routing | Route GPIO traces with min. 0.2mm trace width | SH-ESP32.kicad_pcb, SH-ESP32.pro |
| [`d5b3d26b7`](https://github.com/hatlabs/SH-ESP32-hardware/commit/d5b3d26b73fcd3533a325e7bb537b3685fc582f7) | 2020-10-26 | power, routing | Trace USB 5V to buck converter input | power-input.sch |
| [`96d130924`](https://github.com/hatlabs/SH-ESP32-hardware/commit/96d13092433221feb17ef9c2bf7334e91394f0e4) | 2020-10-26 | footprint | Specified JLCPCB footprints for components | ESP32.sch, I2C.sch, SH-ESP32.kicad_pcb, UI.sch, USB.sch, canbus.sch, onewire.sch, optocoup |
