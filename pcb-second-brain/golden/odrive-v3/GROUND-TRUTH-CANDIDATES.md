# Ground-truth candidates — odrive-v3

Status: **unjudged candidates** produced in W1 (cloud session 2). Each item is a
lead for a revision-diff ground-truth finding (docs/08 § 8.1). Confirmation needs
the kernel's DesignIR diff of the two revisions plus a human look; do not treat
any line below as a confirmed finding.

- Repository: https://github.com/odriverobotics/ODriveHardware
- Pinned commit: `079f30ac58f193226fcc552e4d6ba24553027ee1` (2018-11-15T19:48:32-08:00)
- Paths analysed: `v3`
- Domain: motor_drive
- Commits touching these paths: 46; touching design files: 34; design-file commits whose message matches a fix keyword: 14

## Revision pairs

No revision directories; revisions exist only as git history (pairs = consecutive commits below).

Tags in repository (5): `v3.1`, `v3.2`, `v3.3`, `v3.4`, `v3.4b`

## Revision notes and errata in the repository

Source: `v3/CHANGELOG.md` at `079f30ac58` (verbatim)

```text
## [3.5] - 2018-04-24
### Added
* Ground planes: 4 layer board
  * Logic side planes are GND, AGND
  * Power side planes are PGND, DC_BUS
* GPIO 6,7,8
  * Expose TIM2 CH 1,2
  * One extra ADC capable pin
  * See: https://docs.google.com/spreadsheets/d/1QXDCs1IRtUyG__M_9WruWOheywb-GhOwFtfPcHuN2Fg/edit#gid=949495736
* Two channel DIP switch:
  * Actuate BOOT0 pin to enter bootloader
  * Connect/disconnect CAN termination resistor
* Optional RC filters on step/dir/en and UART lines: GPIO 1,2,6,7,8.
* Strong pull-down resistors on the control signals to the AUX half-H to counteract weak pullups actiavated in DFU mode.
* Two GND pins on J3
* 5V pin on J2
### Changed
* Moved DVDD cap from AGND to GND
* Moved GPIO_5 to PC4
* Change C1 from 2.2uF to 270pF for cycle by cycle Vbus check
* ~~Board voltage determined by stickers (not silkscreen)~~
### Removed
* DC_CAL signals
* AUX current sense amplifier
### Fixed
* Clearance of shunt resistor kelvin connections now enforced with net-tie
* Silkscreen on J2 VCC -> 3.3V
## [3.4] - 2017-11-04
### Changed
* Fixed grounding issue: split digital and power grounds
* Fixed extra apparent shunt resistance, which caused a current measurement scale error
* Change VCC silkscreen to 3.3V
## [3.3] - 2017-06-17
### Added
* Power LED
* Pads on DC rail, for secondary power connection, on the opposite side of main power connector
### Changed
* GPIO Pinout. See: https://docs.google.com/spreadsheets/d/1QXDCs1IRtUyG__M_9WruWOheywb-GhOwFtfPcHuN2Fg/edit#gid=404444347
  * This enables the use of TIM5 on all GPIO pins, which is good for step/dir capture.
  * All GPIO pins are 5v tolerant.
  * 5th GPIO pin now available
* Fliped LDO and C29 to bot side, for better clearance on top side for enclosure style heatsinking.
* Changed all screw terminals to 7.62mm pitch
* Encoder pullup resistance 1k -> 3.3k, to be compatible with 2mA push-pull encoders without desoldering.
### Removed
* AUX_V voltage measurement no longer available.
## [3.2] - 2017-04-19
### Added
* Test-points for most interesting analogue signals
* Break out of the SPI port
* Precision LDO for AVCC (analog VCC)
* Fiducials
### Changed
* Fix current sensor filter capacitor values
* Consolidate all 0.1" headers into two long rows
* Change MOSFETs to a slightly better Rdson version (30% better)
* Adjust heatsink holes
* Improve labling
### Removed
* The confusing EV (Encoder Voltage) rail.
```

## Issue-tracker leads

GitHub issues matching hardware words (search API, 2026-11-04); titles only, unjudged.

| # | State | Opened | Title |
|---|---|---|---|
| [6](https://github.com/odriverobotics/ODriveHardware/issues/6) | open | 2018-10-10 | Top Level not current |
| [9](https://github.com/odriverobotics/ODriveHardware/issues/9) | open | 2020-03-15 | RESC1608X55N footprint not found |
| [10](https://github.com/odriverobotics/ODriveHardware/issues/10) | closed | 2020-03-30 | Schematics and PCB components mismatch |
| [13](https://github.com/odriverobotics/ODriveHardware/issues/13) | open | 2021-10-18 | where is version 3.5  hardware opensource about ODrive？ |

## Commits that look like design fixes

Design-file commits whose subject matches a fix keyword. Categories are keyword
guesses for triage, not judgements. Pair = this commit vs its parent on the same path.

| Commit | Date | Categories | Subject | Design files touched |
|---|---|---|---|---|
| [`3c5231cff`](https://github.com/odriverobotics/ODriveHardware/commit/3c5231cffdbe60e06dfd6f34688c1aca5a59d212) | 2018-01-25 | decoupling, routing | Move GVDD caps from PGND to GND | MotorCell.SchDoc, MotorCell_noPow.SchDoc, PCB.PcbDoc, two_ax_PCB.PrjPcb |
| [`e61b6aa54`](https://github.com/odriverobotics/ODriveHardware/commit/e61b6aa5401e505fc1b7b4019dbe7d4d1fd30499) | 2018-01-25 | power | add second shunt net-ties | MotorCell.SchDoc, MotorCell_noPow.SchDoc, PCB.PcbDoc |
| [`91345da18`](https://github.com/odriverobotics/ODriveHardware/commit/91345da18fae3a5268cd8c8b943fbd2e2e9be276) | 2017-11-07 | silkscreen | add 48V silkscreen | PCB.PcbDoc |
| [`081523bb6`](https://github.com/odriverobotics/ODriveHardware/commit/081523bb6268c7ffd63bfaefb29a0ec39bc27309) | 2017-11-04 | power, routing, generic fix | fix ground issue, shunt via resistance | AuxHalfH.SchDoc, MotorCell.SchDoc, MotorCell_noPow.SchDoc, PCB.PcbDoc, Top.SchDoc |
| [`1eef6e1f6`](https://github.com/odriverobotics/ODriveHardware/commit/1eef6e1f646f29e72c566338595d6b42b40c1d98) | 2017-06-16 | pull-up/down | increase encoder pullups to 3k3 for better compatibility with push-pull encoders | PCB.PcbDoc, Top.SchDoc |
| [`d91cf2c65`](https://github.com/odriverobotics/ODriveHardware/commit/d91cf2c654dca39a4e208ca42587e2ec7a1b117d) | 2017-06-07 | decoupling, footprint | flip cap and vreg, add secondary DC pads | PCB.PcbDoc, two_ax_PCB.PrjPcb |
| [`7e490e93f`](https://github.com/odriverobotics/ODriveHardware/commit/7e490e93fa7b4e56b8f8ae85e8f8da0cdfcb86aa) | 2017-03-25 | part swap, power | change mosfets to more available version (NTMFS4935NT1G) | AuxHalfH.SchDoc, MotorCell.SchDoc, MotorCell_noPow.SchDoc, PCB.PcbDoc |
| [`06d0355c9`](https://github.com/odriverobotics/ODriveHardware/commit/06d0355c91e7fec8dc4314858e72658bdd737771) | 2017-03-24 | power | add precision LDO for AVCC | PCB.PcbDoc, Top.SchDoc |
| [`d6031c8d7`](https://github.com/odriverobotics/ODriveHardware/commit/d6031c8d7883ffc44c3b34d317d116640fc6903b) | 2017-03-22 | part swap, generic fix | Fix current filter values | AuxHalfH.SchDoc, MotorCell.SchDoc, MotorCell_noPow.SchDoc, PCB.PcbDoc |
| [`964f18ce3`](https://github.com/odriverobotics/ODriveHardware/commit/964f18ce3e10ad9284918a2d1e26525948f780f3) | 2017-01-23 | part swap, silkscreen, generic fix | Corrected Silk Screen Mistake ( M0 M1 Designators on Encoder inputs swapped) | PCB.PcbDoc |
| [`2b0d94b72`](https://github.com/odriverobotics/ODriveHardware/commit/2b0d94b72755649d22349dd416736732d2cbcac1) | 2016-12-03 | part swap | schematic pdf and BOM updated to v3.1rc2 (add USB test points) | BOM_v3.0_nonsplit.xlsx, BOM_v3.0_split.xlsx, BOM_v3.1rc.csv, BOM_v3.1rc.xlsx, BOM_v3.1rc_b |
| [`796034994`](https://github.com/odriverobotics/ODriveHardware/commit/79603499418e49d9d148c3f38760d3a611c1784a) | 2016-07-24 | part swap | add schematic and bom | BOM_v3.1rc.xlsx |
| [`6f7674688`](https://github.com/odriverobotics/ODriveHardware/commit/6f7674688501767a2984523492e57630fb41185a) | 2016-06-25 | generic fix | Done with first round of fixes for v3.1 | AuxHalfH.SchDoc, MotorCell.SchDoc, MotorCell_noPow.SchDoc, PCB.PcbDoc, Top.SchDoc |
| [`80b2d348a`](https://github.com/odriverobotics/ODriveHardware/commit/80b2d348ae9ccf957b1a599164ea565a42c27e81) | 2016-06-25 | part swap | rename v3.0 bom saves | v3.0_BOM_nonsplit.xlsx, v3.0_BOM_split.xlsx |
