# Ground-truth candidates — bms-c1

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/LibreSolar/bms-c1
- Pinned commit: `0ca09706f49906bdc73f6ccf99635de349d47474` (2026-03-29T16:03:38+02:00)
- Paths analysed: `.`
- Domain: power_electronics
- Commits touching these paths: 57; touching design files: 24; design-file commits whose message matches a fix keyword: 9

## Revision pairs

No revision directories; revisions exist only as git history (pairs = consecutive commits below).

## Revision notes and errata in the repository

Source: `testing/v0.3/README.md` at `0ca09706f4` (lines matching fix/change keywords only)

```text
During testing some bugs were fixed. The final working firmware commit after testing is [408d783](
Previous issues with disconnected serial port during Bluetooth operation were resolved in the firmware.
There were several issues with the Zephyr I2C driver for the ESP32-C3:
1. [Clock stretching not long enough](https://github.com/zephyrproject-rtos/zephyr/issues/51351)
    The bq76952 requires long clock stretching in the I2C communication, which was working fine with STM32 MCUs, but not with the ESP32-C3 driver. The fix was provided to Zephyr upstream, see above issue for links and further details.
2. [Pin glitches during start-up](https://github.com/zephyrproject-rtos/zephyr/issues/51250)
    This issue is still not fixed and the root cause is unknown. However, the BMS makes a second attempt to communicate with the bq76952 chip, which will succeed if the first one failed because of the glitches.
Output on serial debug console:
- The pre-discharge MOSFET should charge the bus through the pre-discharge resistor before switching on the main MOSFETs.
The BMS switches off the charge MOSFET when the charge overcurrent threshold `Conf/sChgOvercurrent_A` is reached and switches the MOSFET on again after the current changes direction.
The BMS switches off the discharge MOSFET when the discharge overcurrent threshold `Conf/sDisOvercurrent_A` is reached and switches the MOSFET on again after the current changes direction.
The MOSFET and shunt temperatures should stay below 80 °C at 20 °C ambient temperature for the maximum rated current over a period of 60 minutes.
```

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [1](https://github.com/LibreSolar/bms-c1/issues/1) | open | 2021-11-29 | Choosing internal power supply buck regulator |
| [2](https://github.com/LibreSolar/bms-c1/issues/2) | closed | 2021-11-29 | WSMS5515 as shunt resistor |
| [3](https://github.com/LibreSolar/bms-c1/issues/3) | closed | 2021-11-29 | Necessity for two WP-BUCF connectors in parallel? |
| [4](https://github.com/LibreSolar/bms-c1/issues/4) | closed | 2021-12-07 | Plan for board to board connections |
| [5](https://github.com/LibreSolar/bms-c1/issues/5) | closed | 2022-01-05 | Footprints for thermal jumpers? |
| [9](https://github.com/LibreSolar/bms-c1/issues/9) | open | 2022-03-02 | 5V Buck converter |
| [10](https://github.com/LibreSolar/bms-c1/issues/10) | open | 2022-03-02 | Creepage distance Battery voltage |
| [16](https://github.com/LibreSolar/bms-c1/issues/16) | open | 2022-03-04 | Via for all decoupling caps to GND layer |
| [20](https://github.com/LibreSolar/bms-c1/issues/20) | open | 2022-03-10 | What dose this P-MOSFET do? |
| [22](https://github.com/LibreSolar/bms-c1/issues/22) | open | 2022-03-12 | Internal power supply paths |
| [29](https://github.com/LibreSolar/bms-c1/issues/29) | open | 2022-06-28 | Charge pump capacitor value estimation   |
| [40](https://github.com/LibreSolar/bms-c1/issues/40) | open | 2023-04-25 | D16 Schottky diode inverted? |
| [42](https://github.com/LibreSolar/bms-c1/issues/42) | open | 2023-10-09 | Thermal Performance Tests |
| [47](https://github.com/LibreSolar/bms-c1/issues/47) | open | 2023-12-22 | Clarification needed on copper thickness of inner PCB layers  |
| [49](https://github.com/LibreSolar/bms-c1/issues/49) | closed | 2023-12-28 | Minor - correction on schematics |
| [55](https://github.com/LibreSolar/bms-c1/issues/55) | closed | 2024-03-06 | Fix NPTH clearance in v0.4 |
| [65](https://github.com/LibreSolar/bms-c1/issues/65) | open | 2024-11-30 | Power Terminal blocks |
| [66](https://github.com/LibreSolar/bms-c1/issues/66) | open | 2024-12-02 | Use unused GPIOs when RS485 is deactivated |
| [68](https://github.com/LibreSolar/bms-c1/issues/68) | open | 2024-12-09 | Voltage spike damaging LMR during power-on with power supply and electronic load |
| [69](https://github.com/LibreSolar/bms-c1/issues/69) | open | 2025-03-28 | Clarification on the necessity of Diode D9 in the circuit |
| [70](https://github.com/LibreSolar/bms-c1/issues/70) | open | 2025-03-28 | Clarification on C34, C48, C51, and C52 in the BMS-C1 Circuit |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`a0d5743db`](https://github.com/LibreSolar/bms-c1/commit/a0d5743dbbdbb726f75ec7a23e205214f2573f99) | 2026-03-29 | decoupling, protection, generic fix | Add 10mm electrolytic cap to fix overvoltage spikes (#74) | bms-c1.kicad_pcb, bms-c1.kicad_pro, bms-c1.kicad_sch, bms-c1_bom.csv, bq76952.kicad_sch, e |
| [`7545e5262`](https://github.com/LibreSolar/bms-c1/commit/7545e5262a9590aeb754e305fbd84b60d81aa0b5) | 2024-03-17 | clearance | Increase clearance around NPTH (v0.4.1) | bms-c1.kicad_pcb, bms-c1.kicad_pro, bms-c1.kicad_sch, bq76952.kicad_sch, esp32-c3.kicad_sc |
| [`6bcc7ba5e`](https://github.com/LibreSolar/bms-c1/commit/6bcc7ba5e7ce8f8830c14931eb6b16da07300632) | 2024-03-08 | generic fix | fix: Field Remarks of component Q2-Q8 contains extra spaces | power-part.kicad_sch |
| [`21134523d`](https://github.com/LibreSolar/bms-c1/commit/21134523db798ce14475ad67f1700f2478f38876) | 2024-01-17 | power, generic fix | Fix schematic page numbers and MOSFET time constant calc | bms-c1.kicad_sch, power-part.kicad_sch |
| [`cff3aff5c`](https://github.com/LibreSolar/bms-c1/commit/cff3aff5c758026d0cbb1a0389c5dae3b64399e7) | 2022-11-08 | part swap, generic fix | All schematics on A4 and small fixes for BOM | bms-c1.kicad_sch, bms-c1_bom.csv, power-part.kicad_sch |
| [`ca4628db0`](https://github.com/LibreSolar/bms-c1/commit/ca4628db0a9914685ccab9d7d2b0ef831857b408) | 2022-11-03 | protection | Remove FUSE_DETECT circuit | bms-c1.kicad_pcb, bms-c1.kicad_sch, bq76952.kicad_sch |
| [`13485ee33`](https://github.com/LibreSolar/bms-c1/commit/13485ee3315776861d500324442db5cf8893d2f3) | 2022-10-28 | part swap | Swap D and S of Q12 | bms-c1.kicad_pcb, power-part.kicad_sch |
| [`24033fa86`](https://github.com/LibreSolar/bms-c1/commit/24033fa86f239aaabd0f209f86e4758324fe9423) | 2022-10-28 | part swap | Add ESP32-C3 and switch part numbers | bq76952.kicad_sch, esp32-c3.kicad_sch |
| [`02cf5c7f5`](https://github.com/LibreSolar/bms-c1/commit/02cf5c7f53853a723a550bbb7dd0cb178a9694f6) | 2022-09-27 | part swap | BOM updates with a few more alternative parts | bms-c1.kicad_pcb, bms-c1.kicad_pro, bms-c1.kicad_sch, bms-c1_bom.csv, bq76952.kicad_sch, e |
