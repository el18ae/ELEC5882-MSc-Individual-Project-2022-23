# Ground-truth candidates — glasgow

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/GlasgowEmbedded/glasgow
- Pinned commit: `cc9a083c827ba729d390c925af2af1e07ee1de97` (2026-11-03T03:17:27+00:00)
- Paths analysed: `hardware/boards/glasgow`
- Domain: mixed_signal_instrumentation
- Commits touching these paths: 356; touching design files: 343; design-file commits whose message matches a fix keyword: 216

## Revision pairs

Consecutive revision directories present at the pinned commit (diff each pair):

- `hardware/boards/glasgow/`: `revA` → `revB` → `revC0` → `revC1` → `revC2` → `revC3` → `revD0` → `revD1`
  - pair `revA` → `revB`
  - pair `revB` → `revC0`
  - pair `revC0` → `revC1`
  - pair `revC1` → `revC2`
  - pair `revC2` → `revC3`
  - pair `revC3` → `revD0`
  - pair `revD0` → `revD1`

## Revision notes and errata in the repository

Source: `docs/manual/src/revisions/revC3.rst` at `cc9a083c82` (lines matching fix/change keywords only)

```text
There are two main independent VIO rails - one for port A, and one for port B. These rails are adjustable and set the IO voltage levels for the ports. Each VIO rail voltage can be independently configured, which allows you to work with different logic levels on ports A and B simultaneously. The VIOA and VIOB rails are exposed as pins on the port A and port B connectors.
The TPS73101DBV regulators feature reverse current blocking which prevents current from being sunk into the regulator instead of sourced from it. They also have a unique foldback current limit characteristic which provides excellent protection against short circuits on the VIO rails - see the "Internal Current Limit" section and Figure 12 in the datasheet.
Additional 330 mΩ resistors (R56, R7) in series with the regulator outputs ensure recovery in the event that a VIO rail is shorted to ground - see `GitHub issue #135 <https://github.com/GlasgowEmbedded/glasgow/issues/135>`__ for details.
The LVDS port will be replaced with different connectors in future hardware revisions of Glasgow, so its use is not preferred for addon boards.
The pins can be used in differential mode (N/P pairs) or in single-ended mode (independent signals on N and P).
No termination resistors are included. You should include termination resistors on your board if you use the LVDS connector. See the Lattice document `FPGA-TN-02213 "Using Differential I/O (LVDS, Sub-LVDS) in iCE40 LP/HX Devices" <https://www.latticesemi.com/view_document?document_id=47960>`__ for details.
Two pads can be found on the board marked "RECOVER", next to the FX2 EEPROM (U2). This footprint is R40 in the schematic. To initiate recovery, short these pads together and press the reset button, then remove the short. This temporarily changes the I2C address of the FX2 EEPROM so that it boots without firmware, placing it into a recovery mode where it enumerates with the default FX2 device descriptor. The ``fx2tool`` utility can then be used to make a backup copy of the FX2 EEPROM, and the ``glasgow factory`` command can be used to re-provision the device configuration block and the firmware
```

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [68](https://github.com/GlasgowEmbedded/glasgow/issues/68) | open | 2018-10-15 | Expose on-board I2C for accessories |
| [151](https://github.com/GlasgowEmbedded/glasgow/issues/151) | open | 2019-08-19 | Selftest is not supported on revC |
| [152](https://github.com/GlasgowEmbedded/glasgow/issues/152) | open | 2019-08-19 | Add selftest for sync and aux pins |
| [169](https://github.com/GlasgowEmbedded/glasgow/issues/169) | open | 2019-10-06 | revD minispec |
| [187](https://github.com/GlasgowEmbedded/glasgow/issues/187) | open | 2020-03-01 | applet.interface.i2c_initiator: timeout & error message when SCL stuck low |
| [234](https://github.com/GlasgowEmbedded/glasgow/issues/234) | open | 2020-11-14 | Long-term maintainability goals |
| [283](https://github.com/GlasgowEmbedded/glasgow/issues/283) | closed | 2021-05-16 | Test 1.2 V capable level shifters 74AXP1T45 |
| [287](https://github.com/GlasgowEmbedded/glasgow/issues/287) | open | 2021-11-13 | GPIB applet |
| [312](https://github.com/GlasgowEmbedded/glasgow/issues/312) | open | 2023-03-30 | Documentation |
| [765](https://github.com/GlasgowEmbedded/glasgow/issues/765) | open | 2025-03-31 | SP3012 ESD protection part is no longer available |
| [965](https://github.com/GlasgowEmbedded/glasgow/issues/965) | open | 2025-07-20 | Some notes on using `uv` for Glasgow firmware tasks |
| [1096](https://github.com/GlasgowEmbedded/glasgow/issues/1096) | open | 2026-01-01 | cannot configure pulls for port B: Vio is off |
| [1154](https://github.com/GlasgowEmbedded/glasgow/issues/1154) | open | 2026-04-19 | Timing sequence in custom applet |
| [1157](https://github.com/GlasgowEmbedded/glasgow/issues/1157) | open | 2026-04-25 | (add --unsafe-hold-vio) Preserve VIO rails across applet load? |
| [1176](https://github.com/GlasgowEmbedded/glasgow/issues/1176) | open | 2026-06-02 | revD0 pre-launch feature checklist |
| [1184](https://github.com/GlasgowEmbedded/glasgow/issues/1184) | closed | 2026-06-23 | revD1 change requests |
| [1190](https://github.com/GlasgowEmbedded/glasgow/issues/1190) | open | 2026-07-03 | Collect and classify Octal SPI RAM properties for revD0 |
| [1229](https://github.com/GlasgowEmbedded/glasgow/issues/1229) | open | 2026-08-08 | selftest appears to have bitrotted |
| [1302](https://github.com/GlasgowEmbedded/glasgow/issues/1302) | open | 2026-09-25 | revD1 errata |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`24155b54a`](https://github.com/GlasgowEmbedded/glasgow/commit/24155b54a57b62ce8ba2424e33d1ad50c9a6e029) | 2023-10-08 | part swap | revC3: Updated bom to drop the excluded parts. | bom.csv |
| [`1fabd6f23`](https://github.com/GlasgowEmbedded/glasgow/commit/1fabd6f23a9f539afdc9fdbe1ad63fbec4fd5653) | 2023-10-08 | part swap, generic fix | revC3: Marked non schematic/bom parts correctly. | glasgow.kicad_pcb, glasgow.kicad_sch, io_banks.kicad_sch, io_buffer.kicad_sch |
| [`d36fdf183`](https://github.com/GlasgowEmbedded/glasgow/commit/d36fdf183500783fa84fee89f9814507dc21b071) | 2023-10-08 | part swap | revC3: Use the new KiCad DNP flag for the DNP parts. | glasgow.kicad_sch, io_banks.kicad_sch, io_buffer.kicad_sch |
| [`160114c26`](https://github.com/GlasgowEmbedded/glasgow/commit/160114c26e8c2619de935fbce4b47d675ea36fae) | 2023-10-08 | part swap | revC3: Added generated BOM file. | bom.csv |
| [`e90ca13b2`](https://github.com/GlasgowEmbedded/glasgow/commit/e90ca13b21b2305dede5b1f5a078b8e6f4894134) | 2023-10-06 | part swap, generic fix | revC3: Corrected U2 (CAT24M01) part number. | glasgow.kicad_sch |
| [`23153b96e`](https://github.com/GlasgowEmbedded/glasgow/commit/23153b96e0ca9fcf5c81107e4f6c8f4e0e3b4ac6) | 2023-10-06 | generic fix | revC3: Corrected the USB-C 3D model position. | glasgow.kicad_pcb |
| [`2f9781df7`](https://github.com/GlasgowEmbedded/glasgow/commit/2f9781df7114127df08034bd2de5a4e57ea4f4af) | 2023-10-06 | footprint | revC3: Added paste to USB-C shield anchors. | glasgow.kicad_pcb |
| [`13a5ba7fb`](https://github.com/GlasgowEmbedded/glasgow/commit/13a5ba7fb5c0409f2151c5a8de7d57b78cd85b88) | 2023-10-06 | silkscreen, generic fix | revC3: Hid U12&U21 labels. Corrected APX 1b2 key. | glasgow.kicad_pcb, glasgow.kicad_pro, glasgow.kicad_sch |
| [`6ccf99cec`](https://github.com/GlasgowEmbedded/glasgow/commit/6ccf99cec56771d3da39bedc018c1891aecb0303) | 2023-10-06 | footprint | revC3: Increased min pad gap to 0.22mm. | glasgow.kicad_pcb, glasgow.kicad_sch, io_buffer.kicad_sch |
| [`ecdd0bd40`](https://github.com/GlasgowEmbedded/glasgow/commit/ecdd0bd40798aa6382a2e7ce54095182501ba038) | 2023-10-06 | footprint | revC3: Adjusted button footprint to reflect selected tact switch. | glasgow.kicad_pcb |
| [`601ac0ef2`](https://github.com/GlasgowEmbedded/glasgow/commit/601ac0ef2d2dcd6e16e17ff613d88709acdd32a9) | 2023-10-06 | footprint | revC3: Updated resistor pack footprint. | glasgow.kicad_pcb |
| [`7956511f3`](https://github.com/GlasgowEmbedded/glasgow/commit/7956511f3e8b6dc4cd1390f2142440b6d0398543) | 2023-10-06 | part swap | revC3: Filled in the blanks for the Reset button part number. | glasgow.sch |
| [`997a00fd0`](https://github.com/GlasgowEmbedded/glasgow/commit/997a00fd055eefcbbe3244ac6acfd0f25561a21c) | 2023-10-06 | decoupling | revC3: add an 4.7µF capacitor to VIO_AUX for improved decoupling | glasgow.kicad_pcb, io_banks.sch |
| [`d9ddf781e`](https://github.com/GlasgowEmbedded/glasgow/commit/d9ddf781ea1c4eb057abd4dc16a389146f232ddb) | 2023-10-06 | silkscreen | revC3: add pin 1 marker to silkscreen for LVDS header | glasgow.kicad_pcb |
| [`06f4129f8`](https://github.com/GlasgowEmbedded/glasgow/commit/06f4129f8be701e6867af578a73677cfb306aecc) | 2023-10-06 | routing | revC3: connect the 3v3 pin of the LVDS connector with a second via | glasgow.kicad_pcb |
| [`0c1282c44`](https://github.com/GlasgowEmbedded/glasgow/commit/0c1282c44c25a80561d57bdce51333e7214be862) | 2023-10-06 | silkscreen | revC3: change URL on back silkscreen to https | glasgow.kicad_pcb |
| [`cca2e9dfe`](https://github.com/GlasgowEmbedded/glasgow/commit/cca2e9dfe6d105db29937c6960578b19c8585023) | 2023-10-06 | silkscreen | revC3: add URL to back silkscreen, resize VDAC labels on back silkscreen | glasgow.kicad_pcb |
| [`f51a08c8a`](https://github.com/GlasgowEmbedded/glasgow/commit/f51a08c8a9c0b90aaa6094950dae01bfdf4ad4fa) | 2023-10-06 | clearance | revC3: increase spacing of the 3V3 flood fill on the internal power layer around the mounting holes | glasgow.kicad_pcb |
| [`09486fcaa`](https://github.com/GlasgowEmbedded/glasgow/commit/09486fcaae91c648ad0314a5d77ef7a0351b1039) | 2023-10-06 | silkscreen, generic fix | revC3: fix internal i2c bus testpoint silkscreen labels on top side | glasgow.kicad_pcb |
| [`114a04da9`](https://github.com/GlasgowEmbedded/glasgow/commit/114a04da9d6b59c3812e6792e3704ecb8d11fb3c) | 2023-10-06 | silkscreen | revC3: remove port pinout marking on top silkscreen | glasgow.kicad_pcb |
| [`c19f3ca68`](https://github.com/GlasgowEmbedded/glasgow/commit/c19f3ca68472c3ca30e08eff7752acf00ce77451) | 2023-10-06 | silkscreen | revC3: change silkscreen to "owGglas revC?" to emphasise that it is a development version for now | glasgow.kicad_pcb |
| [`7de4e867a`](https://github.com/GlasgowEmbedded/glasgow/commit/7de4e867ab683db30e3af5bde87901302aa89b25) | 2023-07-31 | part swap, generic fix | revC2: Corrected the CAT24M01 part number. | bom.csv, glasgow.sch |
| [`6e6b9c1cd`](https://github.com/GlasgowEmbedded/glasgow/commit/6e6b9c1cdd28e8cdef3e9b99b67cc0fd8c3d98b4) | 2022-04-29 | part swap, generic fix | revC2: Fix bom having not enough columns on multiple lines | bom.csv |
| [`c740d43d9`](https://github.com/GlasgowEmbedded/glasgow/commit/c740d43d95e2dc6019235ee96b9658cadc217341) | 2021-01-06 | generic fix | Fix typos. NFC. | io_banks.sch |
| [`678b8972a`](https://github.com/GlasgowEmbedded/glasgow/commit/678b8972a88c4011b532f7b2cd820d00587bb919) | 2020-12-27 | part swap | revC2: Generated gerber, bom and renderings. | bom.csv, glasgow.pro |
| [`832f6ff36`](https://github.com/GlasgowEmbedded/glasgow/commit/832f6ff362261ebe629e480f365e3e761b951603) | 2020-12-24 | generic fix | revC2: fix level shifter MPN | io_banks.sch, io_buffer.sch |
| [`1e002b2c5`](https://github.com/GlasgowEmbedded/glasgow/commit/1e002b2c5588dd7d96d3405bbadc0edaf1c6259a) | 2020-12-24 | footprint | revC2: Removed the solder mask minimum width setting. | glasgow.kicad_pcb, glasgow.pro |
| [`34845b82b`](https://github.com/GlasgowEmbedded/glasgow/commit/34845b82b795ba0f32efffb64090ab48c04f15dd) | 2020-12-24 | footprint, generic fix | revC2: Corrected resistor pack footprints. | glasgow.kicad_pcb |
| [`0958a4d60`](https://github.com/GlasgowEmbedded/glasgow/commit/0958a4d60714f6cf45d35b57945410184b08b09e) | 2020-12-24 | footprint, generic fix | revC2: Corrected the SOT-143 footprint. | glasgow.kicad_pcb |
| [`b92753cd9`](https://github.com/GlasgowEmbedded/glasgow/commit/b92753cd9fad9535bccc0499fbb48f551bd8c409) | 2020-12-24 | footprint, generic fix | revC2: Minor footprint type corrections. | glasgow-cache.lib, glasgow.kicad_pcb |
| [`73e73e425`](https://github.com/GlasgowEmbedded/glasgow/commit/73e73e42507e33da4c0c266bdaa2daefb85424f2) | 2020-12-24 | power, routing | revC2: re-route ~MR trace on the bottom layer to not split the GND plane below voltage regulators | glasgow.kicad_pcb |
| [`db28ff51b`](https://github.com/GlasgowEmbedded/glasgow/commit/db28ff51bbd7029c4d7e32020390b45d3ffbec78) | 2020-12-24 | footprint | revC2: Updated the tact switch footprint to include 3d model. | glasgow.kicad_pcb |
| [`ac3a5b1ec`](https://github.com/GlasgowEmbedded/glasgow/commit/ac3a5b1ec8a4fe90d282d6eb5115f6e097c29ba8) | 2020-12-24 | part swap | revC2: Replaced the small SMD tact switch with a slighly bigger TH version. | glasgow.kicad_pcb, glasgow.sch |
| [`91da78478`](https://github.com/GlasgowEmbedded/glasgow/commit/91da784786cf12692f76d12fa578939d7780f9c5) | 2020-12-24 | footprint, routing | revC2: Connected a few remaining ground pad via solid to ground plane. | glasgow.kicad_pcb |
| [`9d1ae6e29`](https://github.com/GlasgowEmbedded/glasgow/commit/9d1ae6e29e06247ea62450f0e91d9c917d905164) | 2020-12-24 | footprint, clearance, generic fix | revC2: fix paste-in-hole and wave solder keepout | glasgow.kicad_pcb |
| [`ae6874945`](https://github.com/GlasgowEmbedded/glasgow/commit/ae6874945847c51e3806d1ff93cf3cddeb576a5a) | 2020-12-24 | generic fix | revC2: Fix review notes | glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`929334e5b`](https://github.com/GlasgowEmbedded/glasgow/commit/929334e5ba80b6006101ecb73e4e6abd084030a6) | 2020-12-24 | silkscreen | revC2: More work on silkscreen | glasgow.kicad_pcb |
| [`427f1b66f`](https://github.com/GlasgowEmbedded/glasgow/commit/427f1b66f5e6745a6fc3121a0d00fbd2d7314f52) | 2020-12-24 | silkscreen | revC2: Clean up silkscreen & ADC routing | glasgow.kicad_pcb, glasgow.sch |
| [`eef7829a7`](https://github.com/GlasgowEmbedded/glasgow/commit/eef7829a70af002fbb2c62dd22ccdeaba6f69a57) | 2020-12-24 | footprint, protection, generic fix | revC2: Fix ESD protection footprint & clean up routing | glasgow.kicad_pcb, io_buffer.sch |
| [`977f056a1`](https://github.com/GlasgowEmbedded/glasgow/commit/977f056a18e6815bb319f57cf9eb6e0815705f1b) | 2020-12-24 | footprint, clearance, generic fix | revC2: fix drc issues at the lvds connector with chamfered pads | glasgow.kicad_pcb |
| [`afc96df91`](https://github.com/GlasgowEmbedded/glasgow/commit/afc96df91b0f387a7b41f80081b158ca0e0f8bc1) | 2020-12-24 | part swap | revC2: Consolidated sub 1k resistor value notation. | glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`28202daa9`](https://github.com/GlasgowEmbedded/glasgow/commit/28202daa9b2ccf0440d6c31ad2c1f78d3d24ea4d) | 2020-12-24 | clearance, generic fix | revC2: Addressed some DRC issues. | glasgow.kicad_pcb, glasgow.pro |
| [`7a09b9139`](https://github.com/GlasgowEmbedded/glasgow/commit/7a09b9139a25eb1fd45185f0b21ff9def472f963) | 2020-12-24 | part swap, generic fix | revC2: Updated 1b2-bom-keys and corrected levelshifter PN. | glasgow.kicad_pcb, glasgow.pro, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`9ff5f366d`](https://github.com/GlasgowEmbedded/glasgow/commit/9ff5f366da361bde56c11871a0613e404f84838b) | 2020-12-24 | generic fix | revC2: Corrected the reset button 3d model. | glasgow.kicad_pcb |
| [`070715070`](https://github.com/GlasgowEmbedded/glasgow/commit/0707150709d81c0640fe410548b96fb2411be6e2) | 2020-12-24 | footprint | revC2: Increased the solder paste apertures for paste in hole. | glasgow.kicad_pcb |
| [`b460dfdaf`](https://github.com/GlasgowEmbedded/glasgow/commit/b460dfdaf9182d88cb4d04b080c2293a1be50a1f) | 2020-12-24 | part swap | revC2: Replaced top actuated tact switch with a side actuated one. | glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`0274a671f`](https://github.com/GlasgowEmbedded/glasgow/commit/0274a671fe09d1d36664bcd88fd0c15ca8ebccb6) | 2020-12-24 | routing | revC2: add more vias to layer changes of the main +5V rail | glasgow.kicad_pcb |
| [`6c5712989`](https://github.com/GlasgowEmbedded/glasgow/commit/6c571298971e3182acb91a638dbea5e2fefdc12c) | 2020-12-24 | part swap, protection | revC2: add ESD warning to LVDS connector, add note about it not being available on future revs | glasgow-cache.lib, io_banks.sch |
| [`948c60010`](https://github.com/GlasgowEmbedded/glasgow/commit/948c6001055eaad88b96c93e91e6fb2e628f97e4) | 2020-12-24 | footprint, protection | revC2: Added and assigned the UDFN-14 3d model for the IO ESD Diode packages. | glasgow.kicad_pcb |
| [`eb2e790f2`](https://github.com/GlasgowEmbedded/glasgow/commit/eb2e790f2163746419f9be540f5b035137f06c22) | 2020-12-24 | footprint, part swap | revC2: Replaced 0603 HandSolder footprints with regular ones. | glasgow.kicad_pcb, glasgow.pro, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`ee28161d8`](https://github.com/GlasgowEmbedded/glasgow/commit/ee28161d82a9277b2d8a3b37f10d96484f0d74a4) | 2020-12-24 | part swap | revC2: add bom keys and datasheet for button + reset ic | glasgow.sch |
| [`84729c4bf`](https://github.com/GlasgowEmbedded/glasgow/commit/84729c4bfb6c9005f6c572973c7f7fcc1cbf96a0) | 2020-12-24 | silkscreen | revC2: adapt text about recovery in schematics to match current silkscreen + reset button | glasgow.sch |
| [`b48c859eb`](https://github.com/GlasgowEmbedded/glasgow/commit/b48c859ebb90ecb6419259d093100d03d8106f5e) | 2020-12-24 | decoupling, footprint, generic fix | revC2: move polymer capacitor to allow reworking the pads of the FX2 | glasgow.kicad_pcb |
| [`df0005cff`](https://github.com/GlasgowEmbedded/glasgow/commit/df0005cff1d3e1b195891b2b1bd56d2ea58728f1) | 2020-12-24 | silkscreen, protection | revC2: improve silkscreen labes for LEDs, pulls and ESD protection | glasgow.kicad_pcb |
| [`38295629f`](https://github.com/GlasgowEmbedded/glasgow/commit/38295629f16285c346aa603f915dbb9fb17d1b01) | 2020-12-24 | protection, routing | revC2: center the vias to GND for the TVS diode arrays | glasgow.kicad_pcb |
| [`7659ca962`](https://github.com/GlasgowEmbedded/glasgow/commit/7659ca962aeabb339da461515fd14f9defee086e) | 2020-12-24 | generic fix | revC2: fix distance between port headers to multiple of 2.54 mm | glasgow.kicad_pcb |
| [`20d14396f`](https://github.com/GlasgowEmbedded/glasgow/commit/20d14396f9683f74f64d21d3b77a58ea87ae3035) | 2020-12-24 | silkscreen | revC2: Minor silkscreen work. Added tact switch 3d model. | glasgow.kicad_pcb |
| [`3bbef73de`](https://github.com/GlasgowEmbedded/glasgow/commit/3bbef73de7cb864e674ee1eaf43815b5120f0081) | 2020-12-24 | decoupling, footprint | revC2: Bigger footprint bulk cap and reset button. | glasgow.kicad_pcb, glasgow.sch |
| [`7c8027f17`](https://github.com/GlasgowEmbedded/glasgow/commit/7c8027f17be2321ca56d90b0d715c55ce5d94efb) | 2020-12-24 | decoupling, footprint | revC2: Switched bulk cap footprint, Added system reset button circuit. | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`1d9010a1d`](https://github.com/GlasgowEmbedded/glasgow/commit/1d9010a1d5289e2003241fb301244e52d9ee63e6) | 2020-12-24 | footprint | revC2: Switched eeprom packages to TSSOP-8 to gain more space. | glasgow.kicad_pcb, glasgow.sch |
| [`ce87e43fc`](https://github.com/GlasgowEmbedded/glasgow/commit/ce87e43fcc03652ded725270fa0b50ad7498494d) | 2020-12-24 | footprint, routing, generic fix | revC2: fix track too close to pads of U32 | glasgow.kicad_pcb |
| [`8750308b6`](https://github.com/GlasgowEmbedded/glasgow/commit/8750308b61a43ac823c2de422b9cfa8871686232) | 2020-12-24 | part swap | revC2: save one BOM line by replacing 1k resistors with 2k2 | glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`32f1e30b7`](https://github.com/GlasgowEmbedded/glasgow/commit/32f1e30b7c53c5bbde177e338d4adfda326bab99) | 2020-12-24 | power | revC2: separate current shunt for the INA233 | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`d45c87689`](https://github.com/GlasgowEmbedded/glasgow/commit/d45c87689397e0649e16c3c5d4f63b0478087b44) | 2020-12-24 | part swap | revC2: save one BOM line by replacing a 220k resistor against 2x 100k | glasgow.kicad_pcb, glasgow.sch |
| [`5543524de`](https://github.com/GlasgowEmbedded/glasgow/commit/5543524de3044c341ee51a15d87d779c69697b35) | 2020-12-24 | silkscreen, protection | revC2: re-shuffle silkscreen labels ("Glasgow", creators, ports, OSHW, Sync) to work with the new TVS diodes | glasgow.kicad_pcb |
| [`95cd21a12`](https://github.com/GlasgowEmbedded/glasgow/commit/95cd21a1236f43da1ed0a0cdfc1e2edb38260144) | 2020-12-24 | protection | revC2: place new TVS diodes for port B | glasgow.kicad_pcb |
| [`171e74e03`](https://github.com/GlasgowEmbedded/glasgow/commit/171e74e035d5bd2276ed26768e9cd1d1f710d2a3) | 2020-12-24 | protection | revC2: place new TVS diodes for port A | glasgow.kicad_pcb, io_buffer.sch |
| [`b55d7e779`](https://github.com/GlasgowEmbedded/glasgow/commit/b55d7e77989f7bc37231b402038c004f0a0789b8) | 2020-12-24 | footprint, protection | revC2: add symbol & footprint for new TVS diode array SP3012-06UTG | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`26790eeb0`](https://github.com/GlasgowEmbedded/glasgow/commit/26790eeb02410e368a83871e93aa0ddfb37b558a) | 2020-12-24 | silkscreen, generic fix | revC2: fix & improve silkscreen | glasgow.kicad_pcb |
| [`1c7678d14`](https://github.com/GlasgowEmbedded/glasgow/commit/1c7678d1496038bcb228cd600dd4a5ffd7b2ee68) | 2020-12-24 | generic fix | revC2: roundrect and fixed 3d model for USB-C socket | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`b9a70e7ab`](https://github.com/GlasgowEmbedded/glasgow/commit/b9a70e7ab8b815107221cdb02c3d0f59fd3ba8af) | 2020-12-24 | generic fix | revC2: schematics fixup | io_buffer.sch |
| [`1ac12b1fa`](https://github.com/GlasgowEmbedded/glasgow/commit/1ac12b1fad80ebdd0446f1eb2b1c3b24e7a18489) | 2020-12-24 | footprint, routing | revC2: remove left-over traces, place new footprints | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`8cd8e69fd`](https://github.com/GlasgowEmbedded/glasgow/commit/8cd8e69fd643cb31659358407fc19b42d1cafa02) | 2020-12-24 | part swap | revC2: replace ADC with INA233 voltage & current monitor | glasgow-cache.lib, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`43b0b09ed`](https://github.com/GlasgowEmbedded/glasgow/commit/43b0b09edca24af4dece2b572c89f432ceb9bf9e) | 2020-12-24 | silkscreen | revC2: Made back side silkscreen black on white compatible. | glasgow.kicad_pcb |
| [`f9250a576`](https://github.com/GlasgowEmbedded/glasgow/commit/f9250a576db7c5188a5dc18a0a0f255cc6bc3690) | 2020-12-24 | silkscreen | revC2: Added top silk bank pin legend. | glasgow.kicad_pcb |
| [`dcfd20695`](https://github.com/GlasgowEmbedded/glasgow/commit/dcfd206959abb72c29a253bf0a930c1008e016e7) | 2020-12-24 | routing | revC2: Rerouted and cleaned the USB current limit switch section. | glasgow.kicad_pcb |
| [`1f9c114da`](https://github.com/GlasgowEmbedded/glasgow/commit/1f9c114da097cdefef6dac5ec86521d20693e231) | 2020-12-24 | part swap | revC2: Swapped Micro USB for USB-C USB 2.0 connector. | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`87051f949`](https://github.com/GlasgowEmbedded/glasgow/commit/87051f9495475e3e1800351a5807af15ede76f94) | 2020-12-24 | part swap | revC2: Replaced HiZ with P- to better reflect the actual function. | glasgow.kicad_pcb |
| [`5242e82b5`](https://github.com/GlasgowEmbedded/glasgow/commit/5242e82b5b31c0b5295e3c9d15b1de6e13f0359a) | 2020-12-24 | silkscreen | revC2: Some minor silkscreen legend improvements. | glasgow.kicad_pcb |
| [`d3de585dd`](https://github.com/GlasgowEmbedded/glasgow/commit/d3de585dd654a6f2706ec11e061d8a73c86ca370) | 2020-12-24 | part swap, silkscreen | revC2: Part ID removed from silkscreen, added logic block legends. | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`5a4ed5d86`](https://github.com/GlasgowEmbedded/glasgow/commit/5a4ed5d86bf20e40a67c6a8981b8710bac249381) | 2020-12-24 | footprint | revC2: Spread level shifters to respect courtyards. | glasgow.kicad_pcb |
| [`186f37893`](https://github.com/GlasgowEmbedded/glasgow/commit/186f378933a8e45dd347275a8fce5bbf0a975920) | 2020-12-24 | part swap | revC2: Initial routing pass after changing out parts. | glasgow.kicad_pcb |
| [`aaf907753`](https://github.com/GlasgowEmbedded/glasgow/commit/aaf907753975a423f3a48567abf6b136a1b7f572) | 2020-12-24 | footprint, part swap | revC2: All swapped all remaining missing roundrect footprints. | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`41936941f`](https://github.com/GlasgowEmbedded/glasgow/commit/41936941f10fea0434e9fb3b8efd6c741c8702f7) | 2020-12-24 | footprint | revC2: Updated UDFN-18 to roundrect pads. | glasgow.kicad_pcb |
| [`7f0bdffda`](https://github.com/GlasgowEmbedded/glasgow/commit/7f0bdffdad08c788800d73b51122245b02546427) | 2020-12-24 | footprint | revC2: Updated TSSOP footprint to get roundrect pads. | glasgow.kicad_pcb |
| [`7b225a6b2`](https://github.com/GlasgowEmbedded/glasgow/commit/7b225a6b2b753444c5b6725a8c0a6110738118fb) | 2020-12-24 | footprint, part swap | revC2: Replaced rpack8->4 and level shifter package. | glasgow-cache.lib, glasgow.kicad_pcb, io_banks.sch, io_buffer.sch |
| [`52867b64e`](https://github.com/GlasgowEmbedded/glasgow/commit/52867b64e5a58a0be48bde567f168be481ee7033) | 2020-12-24 | footprint | revC2: Updated a bunch of the footprints to the upstream versions. | glasgow.kicad_pcb |
| [`6f80c3221`](https://github.com/GlasgowEmbedded/glasgow/commit/6f80c3221bd30bf38f2ab41c09c42bc061e573ad) | 2020-12-24 | footprint | revC2: Updated the VREG DFN Footprints to match datasheet. | glasgow.kicad_pcb |
| [`efe49ef6e`](https://github.com/GlasgowEmbedded/glasgow/commit/efe49ef6e7969b4b4c3967ba95ffd45f8e2a9b3b) | 2020-12-24 | footprint | revC2: Changed SOT-563 footprint to improve assembly reliability. | glasgow.kicad_pcb |
| [`9dc35df19`](https://github.com/GlasgowEmbedded/glasgow/commit/9dc35df19e8ae895a371b5e02c49ecd2707cd40d) | 2020-12-24 | footprint, protection, routing | revC2: Bump revision. Remove ESD diode ground pads. | glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`4e0e4e338`](https://github.com/GlasgowEmbedded/glasgow/commit/4e0e4e338db2fb5ceb98426cb56f91d7ed1c449f) | 2020-03-03 | generic fix | revC schematics: fix recommended LVDS connector | io_banks.sch |
| [`d6b3a1c0d`](https://github.com/GlasgowEmbedded/glasgow/commit/d6b3a1c0d3c04e6a7011fe29f292c47f6c3bdb53) | 2020-01-19 | part swap, generic fix | revC1: fix formatting issue in BOM. | bom.csv |
| [`b43653018`](https://github.com/GlasgowEmbedded/glasgow/commit/b4365301885e96e25032f63ef11f387053ba1c7f) | 2019-10-31 | part swap, power | revC1: update fabrication outputs with regulator swap | bom.csv |
| [`70646eb11`](https://github.com/GlasgowEmbedded/glasgow/commit/70646eb11b521e80def833d12eeeb0baf59ca1c9) | 2019-10-31 | part swap, power, generic fix | revC1: update PCB from schematic (fix regulator values) | glasgow.kicad_pcb |
| [`5366d8a45`](https://github.com/GlasgowEmbedded/glasgow/commit/5366d8a45aa464cb3b20432e8b576ef3e3288e77) | 2019-08-23 | part swap, power | revC1: swap regulator library references | glasgow.sch |
| [`f1772cacd`](https://github.com/GlasgowEmbedded/glasgow/commit/f1772cacd46bd7a4121ebad178c7da7128d3a528) | 2019-07-26 | part swap, silkscreen, power | revC1: Swap miss labelled 1V2/3V3 regulator MPNs | glasgow.sch |
| [`2fad90676`](https://github.com/GlasgowEmbedded/glasgow/commit/2fad90676c118e7ec1c90a439bec97425d689809) | 2019-07-24 | part swap | RevC1: Added missing component MPN/Mfg values | glasgow.sch, io_banks.sch |
| [`dcfe0d488`](https://github.com/GlasgowEmbedded/glasgow/commit/dcfe0d488e6d51d2ba02a11cbb49805199872378) | 2019-05-02 | silkscreen | revC1: hide 1V2 label from B.Fab layer | glasgow.kicad_pcb |
| [`3693fd56c`](https://github.com/GlasgowEmbedded/glasgow/commit/3693fd56cd554b1f25e74ce4bea7d17ee1ea532c) | 2019-05-02 | clearance | revC1: mirror "Wave soldering keepout" on Fab.B to match rest of text | glasgow.kicad_pcb, glasgow.pro |
| [`5351475c6`](https://github.com/GlasgowEmbedded/glasgow/commit/5351475c63264edf7bdaa38f7191799d09326072) | 2019-05-02 | routing | revC1: move via away from Sync | glasgow.kicad_pcb |
| [`4fc1b6d2c`](https://github.com/GlasgowEmbedded/glasgow/commit/4fc1b6d2c54d81cad81be9979cf2d2cb56215a6a) | 2019-05-02 | silkscreen, routing | revC1: add G/S labels under Sync connector, judge surrounding vias | glasgow.kicad_pcb |
| [`148f8976f`](https://github.com/GlasgowEmbedded/glasgow/commit/148f8976f71ec6262859e56a9903f41a419a7e0a) | 2019-05-02 | silkscreen | revC1: move AUX and J10 silkscreen labels | glasgow.kicad_pcb |
| [`2f5726aa2`](https://github.com/GlasgowEmbedded/glasgow/commit/2f5726aa2a3c39fd45ac2819c71575f221a7bcb2) | 2019-05-02 | footprint, part swap | revC1: replace Sync connector with proper Molex KK footprint | glasgow.kicad_pcb, io_banks.sch |
| [`b7839673e`](https://github.com/GlasgowEmbedded/glasgow/commit/b7839673e9233466ed58dd1efb4c289bc496e541) | 2019-05-02 | protection, routing | revC1: nudge two vias around the ESD diode arrays | glasgow.kicad_pcb |
| [`e1775ab25`](https://github.com/GlasgowEmbedded/glasgow/commit/e1775ab2587629fea46bc4714f4ed66aaf03a3e5) | 2019-05-01 | part swap | revC1: reduce count of different resistor values | glasgow.kicad_pcb, glasgow.sch, io_banks.sch |
| [`9dc13d39c`](https://github.com/GlasgowEmbedded/glasgow/commit/9dc13d39cbe502053ac0f2d0487a7d18f64bb1c0) | 2019-05-01 | decoupling | revC1: change bulk decoupling cap to polymer in 1b2 key | glasgow.sch |
| [`6486e617c`](https://github.com/GlasgowEmbedded/glasgow/commit/6486e617ce66d5bb939c46265da2107792f3a168) | 2019-05-01 | power, routing | revC1: get rid of spurious thermal connections to GND layer | glasgow.kicad_pcb |
| [`e028678fc`](https://github.com/GlasgowEmbedded/glasgow/commit/e028678fc08088d7781a92a37c23b3cfdf3fecf3) | 2019-05-01 | footprint | revC1: disable solder paste for DNP footprints | glasgow.kicad_pcb |
| [`aaed14145`](https://github.com/GlasgowEmbedded/glasgow/commit/aaed14145c939cf2e444a7cb8f2a6f515eb4cb11) | 2019-04-30 | part swap | revC1: Updated BOM management keys in the schematic. | glasgow.sch, io_buffer.sch |
| [`12b81a937`](https://github.com/GlasgowEmbedded/glasgow/commit/12b81a937c6432f11a0ae64abbaf8c36256b0c07) | 2019-05-01 | part swap | revC1: remove DNP part models & add DNP field in sch | glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`f09055329`](https://github.com/GlasgowEmbedded/glasgow/commit/f09055329a8d0915dcf9bf5fe37469b203553b4e) | 2019-05-01 | routing | revC1: adjust "Serial no" text in the rear to clash less with via | glasgow.kicad_pcb |
| [`50b6be708`](https://github.com/GlasgowEmbedded/glasgow/commit/50b6be708912b5fbf6d4f83ba051d9489b1ca287) | 2019-05-01 | routing | revC1: nudge some traces | glasgow.kicad_pcb |
| [`d7392a517`](https://github.com/GlasgowEmbedded/glasgow/commit/d7392a51751691bf60a1be0790d9cb48ed335b19) | 2019-05-01 | silkscreen | revC1: move C87 up slightly, relocate CLKIF label | glasgow.kicad_pcb, glasgow.sch |
| [`b9dd91f6b`](https://github.com/GlasgowEmbedded/glasgow/commit/b9dd91f6b432ef8424e99bb7ad1c89b48a84bad2) | 2019-05-01 | decoupling, footprint | revC1: update 5V capacitor footprint & datasheet | glasgow.kicad_pcb, glasgow.sch |
| [`7303692ff`](https://github.com/GlasgowEmbedded/glasgow/commit/7303692ff907c4d87eadfcd46298245843af7edc) | 2019-05-01 | decoupling | revC1: bump 5V cap to 150µF | glasgow.kicad_pcb, glasgow.sch |
| [`e741c1703`](https://github.com/GlasgowEmbedded/glasgow/commit/e741c1703ac521decfd3224034e12f67c8774510) | 2019-04-30 | part swap | revC1: rotate U15 90deg, swap D1/D2 pins | glasgow.kicad_pcb, glasgow.sch |
| [`3068ca106`](https://github.com/GlasgowEmbedded/glasgow/commit/3068ca1062f88e89fe6de82f4507d00fde2b54b3) | 2019-04-30 | routing | revC1: add additional current-carrying vias to USB GND | glasgow.kicad_pcb |
| [`4500338ba`](https://github.com/GlasgowEmbedded/glasgow/commit/4500338bacc478baa57e5abeb3feeb28ed7768b5) | 2019-04-30 | pull-up/down, power | revC1: change regulator EN pulldowns 100k -> 10k | glasgow.kicad_pcb, glasgow.sch |
| [`f8da530c2`](https://github.com/GlasgowEmbedded/glasgow/commit/f8da530c281869a04680817b3f54ca18d841d2bc) | 2019-04-30 | decoupling | revC1: add 100uF bulk decoupling to 5V rail | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`2edcf4a44`](https://github.com/GlasgowEmbedded/glasgow/commit/2edcf4a442bd66e9ffa43ab955d7ee6a1b48af6c) | 2019-04-30 | generic fix | revC1: fix LM3880MF EN divider for 4.1V brownout | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`7c1833e7d`](https://github.com/GlasgowEmbedded/glasgow/commit/7c1833e7dc28ba5ebc2199749e23f5d521897934) | 2019-04-25 | silkscreen, generic fix | revC1: fix NONO_HANA silk art geometry | glasgow.kicad_pcb |
| [`f70e56d0c`](https://github.com/GlasgowEmbedded/glasgow/commit/f70e56d0c5456d9561270dce22785f00cbf42351) | 2019-04-25 | routing | revC1: run a trace/via cleanup pass, adjust 1V2 plane settings | glasgow.kicad_pcb |
| [`b79ee6250`](https://github.com/GlasgowEmbedded/glasgow/commit/b79ee6250c1b66132832d6978c927ad4e8cbc41e) | 2019-04-25 | silkscreen, routing | revC1: tweak some traces, move vias away from Sync silkscreen | glasgow.kicad_pcb |
| [`7d2f53cb3`](https://github.com/GlasgowEmbedded/glasgow/commit/7d2f53cb3235cf83fd3375fd782683db720684e0) | 2019-04-25 | silkscreen | revC1: make space for and add silkscreen label for VIO LEDs | glasgow.kicad_pcb |
| [`4bf105207`](https://github.com/GlasgowEmbedded/glasgow/commit/4bf10520723955469a375dc783d0f34641efdc46) | 2019-04-25 | routing | revC1: get rid of 3V3 plane stitch across top side | glasgow.kicad_pcb |
| [`1db6e874e`](https://github.com/GlasgowEmbedded/glasgow/commit/1db6e874ef18b6212fceae03511f95792207b3af) | 2019-04-25 | power, routing | revC1: Add extra vias to the 1V2 regulator top side | glasgow.kicad_pcb |
| [`f3fa13ae0`](https://github.com/GlasgowEmbedded/glasgow/commit/f3fa13ae0584a85b1801de3244eb9d514e6e17aa) | 2019-04-25 | part swap, power | revC1: Replace 1V2/3V3 regulator with TLV755/733 series | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.pro, glasgow.sch |
| [`604c3520e`](https://github.com/GlasgowEmbedded/glasgow/commit/604c3520ee9f38a9895e6e4b76a34f7550d80ce7) | 2019-04-24 | generic fix | revC1: fix fab layer too | glasgow.kicad_pcb |
| [`90ee8994d`](https://github.com/GlasgowEmbedded/glasgow/commit/90ee8994d1a96898b4b52ae3e3868f5de70b9034) | 2019-04-24 | footprint, silkscreen, generic fix | revC1: fix silkscreen borkage due to footprint updates | glasgow.kicad_pcb |
| [`104b36495`](https://github.com/GlasgowEmbedded/glasgow/commit/104b364952a0244178da4c36894f7c16fcaee09f) | 2019-04-24 | generic fix | revC1: add Vio active LEDs, misc fixes | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`2dff08f60`](https://github.com/GlasgowEmbedded/glasgow/commit/2dff08f600a7887b2394ecd1ec296a44ad0ef83a) | 2019-04-23 | generic fix | revC1: fix fab layers | glasgow.kicad_pcb |
| [`455da1545`](https://github.com/GlasgowEmbedded/glasgow/commit/455da1545120d19b0f30831076e0037975c0455d) | 2019-04-23 | silkscreen | revC1: relocate VSENSE divider silkscreen away from screws | glasgow.kicad_pcb |
| [`95e1c267d`](https://github.com/GlasgowEmbedded/glasgow/commit/95e1c267da15014ecad12f48c5bfb5235f19dc3b) | 2019-04-22 | footprint | revC1: move R5 and R6 out of the mounting hole courtyard | glasgow.kicad_pcb |
| [`4289cd5f8`](https://github.com/GlasgowEmbedded/glasgow/commit/4289cd5f80f9ad9f76d636a16f2817c53c002094) | 2019-04-22 | footprint | revC1: increase mounting hole size from 3.2mm to 3.5mm, add courtyard space around | glasgow.kicad_pcb, glasgow.sch |
| [`8ed17e1d5`](https://github.com/GlasgowEmbedded/glasgow/commit/8ed17e1d5acd2b7d015bf0e2eadd06f0a0443f89) | 2019-04-21 | part swap | revC1: Updated part information and datasheet links | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`e601fab7c`](https://github.com/GlasgowEmbedded/glasgow/commit/e601fab7c0eef20e8ef4979bf755236ec989821d) | 2019-04-16 | part swap | revC1: Added BOM management keys to the schematic. | glasgow-cache.lib, glasgow.sch, io_banks.sch, io_buffer.sch |
| [`a0ddfddfa`](https://github.com/GlasgowEmbedded/glasgow/commit/a0ddfddfade4ace668e15811a696041fe3ecf777) | 2019-04-13 | protection | Add ESD protection diodes for the AUX header too | glasgow.kicad_pcb, io_banks.sch |
| [`403f118b2`](https://github.com/GlasgowEmbedded/glasgow/commit/403f118b2c415564e4538f74b024ea1d642a9434) | 2019-04-11 | generic fix | revC: correct the explanation for GBIN0/1 balls being shared. | io_banks.sch |
| [`73b1b1872`](https://github.com/GlasgowEmbedded/glasgow/commit/73b1b187263ee5936d388d196a5d11eaa4c18761) | 2019-04-10 | routing | revC1: simplify a trace | glasgow.kicad_pcb |
| [`d571f7bc3`](https://github.com/GlasgowEmbedded/glasgow/commit/d571f7bc39187170f63b6a339ea0cbe9a6e6d233) | 2019-04-10 | silkscreen, generic fix | revC1: fix USB testpoint silkscreen/fab | glasgow.kicad_pcb |
| [`df2ce4a94`](https://github.com/GlasgowEmbedded/glasgow/commit/df2ce4a94905fe88db6ca704a6329d2aa018093e) | 2019-04-09 | generic fix | revC1: fix U32 reference size again. | glasgow.kicad_pcb |
| [`56c17893c`](https://github.com/GlasgowEmbedded/glasgow/commit/56c17893c32087d2598af18d8400e5daad06652c) | 2019-04-09 | pull-up/down | revC1: add 10k pull-up to SYNC pin | glasgow.kicad_pcb, io_banks.sch |
| [`2ee47442f`](https://github.com/GlasgowEmbedded/glasgow/commit/2ee47442f75d9c00b32db7563954a8d8c67ad902) | 2019-04-08 | protection | revC1: make the routing around the ESD protection more uniform | glasgow.kicad_pcb |
| [`ed9ac321d`](https://github.com/GlasgowEmbedded/glasgow/commit/ed9ac321db8dce9b3fcc1d17e6a48c272fbd6d24) | 2019-04-08 | footprint | revC1: add UDFN-18 package model | glasgow.kicad_pcb |
| [`59c960b1d`](https://github.com/GlasgowEmbedded/glasgow/commit/59c960b1dd66305fbc67da76256224aa00a89035) | 2019-04-07 | protection | revC1: add a TVS protection diode to SYNC | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_banks.sch |
| [`e8fcc6ff8`](https://github.com/GlasgowEmbedded/glasgow/commit/e8fcc6ff8ebe65c08325831bb802440fc6b088be) | 2019-04-07 | generic fix | revC1: fix U32 reference size | glasgow.kicad_pcb |
| [`b9794b542`](https://github.com/GlasgowEmbedded/glasgow/commit/b9794b5428ee635f1f0d5e189049a678958e7abb) | 2019-04-07 | protection | revC1: add ESD protection to IO banks | glasgow-cache.lib, glasgow.kicad_pcb, io_buffer.sch |
| [`7bb28430b`](https://github.com/GlasgowEmbedded/glasgow/commit/7bb28430bdc5d54626f53d2601ab6993494df73a) | 2019-04-07 | generic fix | revC1: fix some ERC errors | glasgow-cache.lib, glasgow.sch, io_banks.sch |
| [`f2202ca51`](https://github.com/GlasgowEmbedded/glasgow/commit/f2202ca51ed4a783cc05943f76aaa7c5755a72a8) | 2019-04-05 | silkscreen | revC1: add io port pinout to back silkscreen | glasgow.kicad_pcb |
| [`f5b00708a`](https://github.com/GlasgowEmbedded/glasgow/commit/f5b00708ade789152f3b53cde2254f4aec06d555) | 2019-04-01 | routing | revC1: tweak 1V2 plane and prettify some traces | glasgow.kicad_pcb |
| [`fe2ecfe07`](https://github.com/GlasgowEmbedded/glasgow/commit/fe2ecfe076c370d71e041c139949304e5d6849be) | 2019-04-01 | routing | revC1: route FPGA pin D5 to CLKREF | glasgow.kicad_pcb, glasgow.sch, io_banks.sch |
| [`8cd719d7d`](https://github.com/GlasgowEmbedded/glasgow/commit/8cd719d7dc9547d9999bdd6daf76c2f0735affc3) | 2019-04-01 | silkscreen, generic fix | revC1: improve Sync change layout, fix silk/fab layers | glasgow.kicad_pcb |
| [`14462ee04`](https://github.com/GlasgowEmbedded/glasgow/commit/14462ee04e9100c16037b93c7acb26942656d7fa) | 2019-03-29 | part swap | revC1: note 10/12-bit ADC/DACs as alternates | io_buffer.sch |
| [`97be25e8d`](https://github.com/GlasgowEmbedded/glasgow/commit/97be25e8d3cce3330301717a5d1e3613d5267d38) | 2019-03-29 | pull-up/down | revC1: make Sync pull-up 1K | glasgow.kicad_pcb, io_banks.sch |
| [`d2af97c88`](https://github.com/GlasgowEmbedded/glasgow/commit/d2af97c888a82de2b2505306e4c7ded6df2eb5a7) | 2019-03-29 | part swap, generic fix | revC1: fix old value/datasheet for PCA6408, style fixes | glasgow.kicad_pcb, io_buffer.sch |
| [`88ee1fe1a`](https://github.com/GlasgowEmbedded/glasgow/commit/88ee1fe1ad456c5076a5c2a75dbcbcb96cad0e26) | 2019-03-29 | part swap | revC1: replace pulls with PCA6408A and remove shifters | glasgow-cache.lib, glasgow.kicad_pcb, io_buffer.sch |
| [`b209775ef`](https://github.com/GlasgowEmbedded/glasgow/commit/b209775efd6ce5ae1e39efe2a8a597271f341474) | 2019-03-08 | generic fix | revC1: fix ugliness in user LED schematic section | io_banks.sch |
| [`e767dfa9f`](https://github.com/GlasgowEmbedded/glasgow/commit/e767dfa9f17ef97259fd2850267091acff6bcea4) | 2019-03-07 | part swap | revC1: adjust LED resistor values for eye safety | glasgow.kicad_pcb, glasgow.sch, io_banks.sch |
| [`daa6a2210`](https://github.com/GlasgowEmbedded/glasgow/commit/daa6a221030bb4860372d27c8204c0f2a22c71a0) | 2019-02-05 | silkscreen | revC1: change C78 silkscreen text to 0.8 for consistency | glasgow.kicad_pcb |
| [`fe6073c1c`](https://github.com/GlasgowEmbedded/glasgow/commit/fe6073c1cc4f26b84503b8b6c74b197e558d2f18) | 2018-12-13 | silkscreen, generic fix | revC: fix silkscreen for FB1 | glasgow.kicad_pcb |
| [`8ae10d255`](https://github.com/GlasgowEmbedded/glasgow/commit/8ae10d25574f583828203a85be58edfd92cbc8c7) | 2018-12-05 | routing | revC: line up some traces better | glasgow.kicad_pcb |
| [`ecd4f6a48`](https://github.com/GlasgowEmbedded/glasgow/commit/ecd4f6a485627cba18e2d110b7bc64a3a0bb1756) | 2018-12-05 | silkscreen | revC: another silkscreen tweak near USB connector | glasgow.kicad_pcb |
| [`6ffb5e554`](https://github.com/GlasgowEmbedded/glasgow/commit/6ffb5e5548b7d1e181169cf8ca16e277800e7837) | 2018-12-05 | silkscreen | revC: minor nudging and tweaking to improve silkscreen around reg area | glasgow.kicad_pcb |
| [`6e30c4398`](https://github.com/GlasgowEmbedded/glasgow/commit/6e30c43988e72cdcc507d0f1ef3e92923dbd97f8) | 2018-12-04 | decoupling, routing | revC: simplify VUSB routing, move input cap to top layer, double vias | glasgow.kicad_pcb |
| [`8015b780e`](https://github.com/GlasgowEmbedded/glasgow/commit/8015b780e77950584ef088aed83984b4900279c1) | 2018-12-04 | silkscreen | revC: remove additional I/O port silkscreen markings | glasgow.kicad_pcb |
| [`aefc5cea9`](https://github.com/GlasgowEmbedded/glasgow/commit/aefc5cea965924ed8e5e6a5530f753dd7cb003d3) | 2018-12-04 | silkscreen, generic fix | revC: schematic style fix (add wire for net label) | glasgow.sch |
| [`6423d6e7d`](https://github.com/GlasgowEmbedded/glasgow/commit/6423d6e7d0b43a3802fc518c94b0e2048927405d) | 2018-12-04 | generic fix | revC: update incorrect comment re: Vbus threshold | glasgow.sch |
| [`2acbf944e`](https://github.com/GlasgowEmbedded/glasgow/commit/2acbf944e0c6a1666a38f081e7a8770f98770213) | 2018-12-04 | footprint | revC: remove spurious square in F.Paste | glasgow.kicad_pcb |
| [`d1cc2bf37`](https://github.com/GlasgowEmbedded/glasgow/commit/d1cc2bf374fbcebc31c77aadc46dc17ecad3bf60) | 2018-12-04 | generic fix | revC: Fix a glitch | glasgow.kicad_pcb |
| [`0700da1c7`](https://github.com/GlasgowEmbedded/glasgow/commit/0700da1c7c6f15042ce38679dcbbbf0edae3d795) | 2018-12-04 | footprint, silkscreen | revC: use proper footprint for I/O connectors, tweak silkscreen | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`e36206501`](https://github.com/GlasgowEmbedded/glasgow/commit/e36206501c4d4b69481dee17c070e0a31eeb148a) | 2018-12-04 | silkscreen | revC: restore Glasgow label text | glasgow.kicad_pcb |
| [`636121dd7`](https://github.com/GlasgowEmbedded/glasgow/commit/636121dd70e2f43a48bd20d786677b7733c3cb79) | 2018-12-04 | silkscreen | revC: silkscreen wrapup, add Hana | glasgow.kicad_pcb |
| [`ec5f4da1f`](https://github.com/GlasgowEmbedded/glasgow/commit/ec5f4da1f9cbdf6c3159fa9bfacf8e640b6b6a29) | 2018-11-28 | generic fix | revC: fix some more diffpair issues | glasgow.kicad_pcb |
| [`b50bd6945`](https://github.com/GlasgowEmbedded/glasgow/commit/b50bd6945d1a3eb7eea7ccf7340f02217b71e300) | 2018-11-28 | generic fix | revC: fix up diffpairs | glasgow-cache.lib, glasgow.kicad_pcb |
| [`147314506`](https://github.com/GlasgowEmbedded/glasgow/commit/147314506397112cc78860d89f1ac69570bd96d8) | 2018-11-28 | clearance | revC: update wave soldering keepouts | glasgow.kicad_pcb |
| [`422462cee`](https://github.com/GlasgowEmbedded/glasgow/commit/422462cee7f7d96bd4f199793903930932df7352) | 2018-11-28 | routing | revC: pull +5V a bit further away from GND | glasgow.kicad_pcb |
| [`528864812`](https://github.com/GlasgowEmbedded/glasgow/commit/528864812ba42dd2a63bbc019815d2f307c3e518) | 2018-11-28 | part swap, clearance, generic fix | revC: swap U2/U3, fix DRC issue | glasgow.kicad_pcb |
| [`bad3e59dc`](https://github.com/GlasgowEmbedded/glasgow/commit/bad3e59dc99a10fb2ed5e45482e42b3f8d12f3df) | 2018-11-28 | routing | revC: length match LVDS pairs (#82) | glasgow.kicad_pcb |
| [`9fa367d48`](https://github.com/GlasgowEmbedded/glasgow/commit/9fa367d486b458ead9225a7a1b1a50cb57f40564) | 2018-11-26 | silkscreen, routing | revC: all routed, top silkscreen cleaned up | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`8f7930656`](https://github.com/GlasgowEmbedded/glasgow/commit/8f7930656b0136b1e41f9d3f0c099652d1f977f7) | 2018-11-26 | routing | revC: I/O buffers fully routed | glasgow.kicad_pcb, io_buffer.sch |
| [`2f1bb88c1`](https://github.com/GlasgowEmbedded/glasgow/commit/2f1bb88c1e4aea56ef1167f4191837f893b89065) | 2018-11-21 | routing | revC: route LVDS | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`e2ccf9fe3`](https://github.com/GlasgowEmbedded/glasgow/commit/e2ccf9fe32291127393212ea4ea14cc9f65734b4) | 2018-11-21 | generic fix | revC: fix connector placement, LVDS connector testing | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`c8b7fceca`](https://github.com/GlasgowEmbedded/glasgow/commit/c8b7fceca9efaa5a0601c204410a39f014a22c89) | 2018-11-21 | footprint, clearance, generic fix | revC: Update SOT563 footprints to Glasgow version and more DRC fixes | glasgow.kicad_pcb, io_buffer.sch |
| [`631f2a0ac`](https://github.com/GlasgowEmbedded/glasgow/commit/631f2a0ac16a87431df5460e2abe49be58b6903d) | 2018-11-21 | clearance, generic fix | revC: fix some misc DRC errors | glasgow.kicad_pcb |
| [`895cffb14`](https://github.com/GlasgowEmbedded/glasgow/commit/895cffb1494d1319d486dbba0facf2c6d2b36be8) | 2018-11-21 | routing | revC: route missing D0 and D1 | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`3b602a50f`](https://github.com/GlasgowEmbedded/glasgow/commit/3b602a50f5cf6661c27fc41db35d17fe10175e9f) | 2018-07-21 | pull-up/down | Further reduce SYNC pullup to 130 Ohms. | glasgow.sch |
| [`6dd57c344`](https://github.com/GlasgowEmbedded/glasgow/commit/6dd57c344ce476ec30139231c9b7966d20e7df7d) | 2018-07-20 | footprint, part swap | Use upstream QFN footprint for Cypress part | glasgow.kicad_pcb, glasgow.sch |
| [`0ecd8ce41`](https://github.com/GlasgowEmbedded/glasgow/commit/0ecd8ce41fe9091416454d65247f53e68d8b5fd1) | 2018-07-20 | pull-up/down, part swap | Reduce SYNC pull-up value to 390 Ohms. | glasgow.sch |
| [`6cb4f9972`](https://github.com/GlasgowEmbedded/glasgow/commit/6cb4f9972d681bcb62dff08dbc81986b0b11e50b) | 2018-06-24 | silkscreen | Minor silkscreen cleanups. | glasgow.kicad_pcb |
| [`de22844e9`](https://github.com/GlasgowEmbedded/glasgow/commit/de22844e9447705d42606db0405bc92e1bf6dc50) | 2018-05-24 | decoupling | Update MPNs for out of stock capacitors and resistors. | glasgow.sch, io_buffer.sch |
| [`7603bbc1c`](https://github.com/GlasgowEmbedded/glasgow/commit/7603bbc1cacc52c9a810be09ee2d881bc6b96353) | 2018-05-24 | generic fix | Fix accidentally used 0201 resistor. | glasgow.sch, io_buffer.sch |
| [`adac2ce59`](https://github.com/GlasgowEmbedded/glasgow/commit/adac2ce5972c1814ccddcc799732834a1a02ab05) | 2018-05-18 | silkscreen | Relabel CYP→FX2, FPGA→ICE. | glasgow.kicad_pcb, glasgow.sch |
| [`1f6ff5db5`](https://github.com/GlasgowEmbedded/glasgow/commit/1f6ff5db559d74836fcc8e9da1f6f277e28fc15a) | 2018-05-14 | routing | Move QB traces around a bit to free up space near to level shifter. | glasgow.kicad_pcb |
| [`f159e6a50`](https://github.com/GlasgowEmbedded/glasgow/commit/f159e6a5047a00f50c10a5c71dd680161241abc9) | 2018-05-12 | footprint | Remove all footprints, symbols and 3D packages merged in KiCAD. | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`479013c69`](https://github.com/GlasgowEmbedded/glasgow/commit/479013c69b0963e352099ae9bfb358a6f10520e3) | 2018-05-12 | silkscreen, routing | Make a few cosmetic changes to the PCB, both traces and silk. | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`d8baf5efb`](https://github.com/GlasgowEmbedded/glasgow/commit/d8baf5efb3267cf8ba190d94dabec0ead78541d3) | 2018-05-11 | silkscreen, generic fix | Fix silkscreen issues with resistors moved to top. | glasgow.kicad_pcb |
| [`b122273eb`](https://github.com/GlasgowEmbedded/glasgow/commit/b122273eb699f6b1ea6af1542ddc6efc5b01c79d) | 2018-05-07 | generic fix | Fix LED current comment on schematic. | glasgow.sch |
| [`16aeae915`](https://github.com/GlasgowEmbedded/glasgow/commit/16aeae915bfe401451bb0ef3e10178b6610063ba) | 2018-05-05 | pull-up/down | Change ENVA and ENVB from pullups to pulldowns. | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
| [`f1afd846d`](https://github.com/GlasgowEmbedded/glasgow/commit/f1afd846dfc9302087bff990a64b0866f4a26cf3) | 2018-05-05 | footprint, generic fix | Correct crystal orientation and shrink footprint. | glasgow.kicad_pcb, glasgow.sch |
| [`4b17418f2`](https://github.com/GlasgowEmbedded/glasgow/commit/4b17418f204dfb253f2e7ee1ac3f9810aee94daf) | 2018-05-03 | generic fix | Fix Mouser_PN for R19, R22. | io_buffer.sch |
| [`9958261ed`](https://github.com/GlasgowEmbedded/glasgow/commit/9958261ed31fddcda8fd56d32ff1e327127d0b9a) | 2018-05-03 | protection, generic fix | Fix ENVA/ENVB pin polarity (oops!) | glasgow.kicad_pcb, glasgow.pro, glasgow.sch, io_buffer.sch |
| [`690d0465b`](https://github.com/GlasgowEmbedded/glasgow/commit/690d0465b658b203f40680ad8a1ac3d878b7ba6a) | 2018-05-03 | generic fix | Fix DAC addresses on schematic. | glasgow.pro, glasgow.sch, io_buffer.sch |
| [`ab949d0f5`](https://github.com/GlasgowEmbedded/glasgow/commit/ab949d0f527a91314e2e70c1c3ed5f9271ce6fc4) | 2018-05-02 | part swap, generic fix | Fix R16 value. | glasgow-cache.lib, glasgow.kicad_pcb, glasgow.sch |
| [`a9274fca4`](https://github.com/GlasgowEmbedded/glasgow/commit/a9274fca428e4a93e14e9ef8e9a31f62f3bb6596) | 2018-05-01 | generic fix | Fix layer visibility. | glasgow.kicad_pcb |
| [`9e5878091`](https://github.com/GlasgowEmbedded/glasgow/commit/9e5878091c16ef62f693c4240095aca1aad9cf40) | 2018-04-30 | footprint | Remove hardware/packages3D/LED_0603_1608Metric_Castellated.{step,wrl}. | glasgow.kicad_pcb |
| [`bc3d7f1ea`](https://github.com/GlasgowEmbedded/glasgow/commit/bc3d7f1eaefa702c2304ccd71329c2071e26ee26) | 2018-04-25 | part swap | Specify part number for J4. | glasgow.sch |
| [`297b69b62`](https://github.com/GlasgowEmbedded/glasgow/commit/297b69b62e130c9e36ccb18bb3433bd3ee1061eb) | 2018-04-25 | decoupling, generic fix | Fix u1 capacitors to be 0603, not 0201. | glasgow.sch, io_buffer.sch |
| [`434dbef0a`](https://github.com/GlasgowEmbedded/glasgow/commit/434dbef0a01f4441d2b0f0cbba04a322aa2d1d4c) | 2018-04-18 | power, routing | Enlarge CY7C thermal vias to 0.3mm to satisfy board house capabilities. | glasgow.kicad_pcb |
| [`ae857332c`](https://github.com/GlasgowEmbedded/glasgow/commit/ae857332c332ea19b1edad8fbe3ea5b355b45f48) | 2018-04-18 | footprint | Restore paste layer for back side pads. | glasgow.kicad_pcb |
| [`8190114d1`](https://github.com/GlasgowEmbedded/glasgow/commit/8190114d19272869687e8c1bf6c5f9882abe56fe) | 2018-04-17 | footprint, power, routing | DFM tweaks to MSOP paste and thermal vias. | glasgow.kicad_pcb |
| [`97c7eef2c`](https://github.com/GlasgowEmbedded/glasgow/commit/97c7eef2ca8c2cd7d4bedae72bcd77c7ca1aaac0) | 2018-04-18 | generic fix | Fix LM3880 EN divider to have correct threshold (thanks @awygle!). | glasgow.kicad_pcb, glasgow.sch |
| [`e52160d45`](https://github.com/GlasgowEmbedded/glasgow/commit/e52160d45377ac36d05f619468c4bdc675f85293) | 2018-04-18 | power | Leave only ADC on Vsense, and connect LDO and FXMA both to Vio. | glasgow.kicad_pcb, io_buffer.sch |
| [`e33e830a5`](https://github.com/GlasgowEmbedded/glasgow/commit/e33e830a541b17e24d5afd708d34dc05e82a25ec) | 2018-04-18 | generic fix | Fix some unconnected nets. | glasgow.kicad_pcb |
| [`aa13ed2bd`](https://github.com/GlasgowEmbedded/glasgow/commit/aa13ed2bdd56b9992f2957cf43a43f9c86536c9d) | 2018-04-17 | routing | Add explicit pin1 mark for CY7C as the IPC one is obscured by via. | glasgow.kicad_pcb |
| [`3ea846fa1`](https://github.com/GlasgowEmbedded/glasgow/commit/3ea846fa146803d47164fd846f0a5287d789d3b4) | 2018-04-17 | decoupling | Remove decoupling capacitor from non-power (WAKEUP) pin. | glasgow.kicad_pcb, glasgow.sch |
| [`b172ac20a`](https://github.com/GlasgowEmbedded/glasgow/commit/b172ac20aa2afbaf9910aca37934f35326a46432) | 2018-04-17 | part swap | Swap INT0 and INT1 to give ~ALERT highest natural priority. | glasgow.kicad_pcb, glasgow.sch, io_buffer.sch |
